# rules/cage_ctss.smk

rule build_ctss_and_bigwig:
    """
    Count 5' ends (CTSS) from mapped reads using bedtools coverage,
    then generate raw and CPM-normalized 5' bigWig tracks.
    """
    input:
        bam=lambda wc: final_bam_path(wc.sample),
        bai=lambda wc: final_bai_path(wc.sample),
        chrom_sizes=config["reference"]["chrom_sizes"]
    output:
        ctss=f"{OUTDIR}/ctss/{{sample}}/{{sample}}.ctss.tsv",
        raw_bg=temp(f"{OUTDIR}/ctss/{{sample}}/{{sample}}.5prime.raw.bedgraph"),
        cpm_bg=temp(f"{OUTDIR}/ctss/{{sample}}/{{sample}}.5prime.cpm.bedgraph"),
        raw_bw=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.5prime.raw.bw",
        cpm_bw=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.5prime.cpm.bw"
    params:
        exclude_flags=config.get("samtools_exclude_flags", 1804)
    threads: int(config["threads"]["samtools"])
    log:
        "logs/ctss/{sample}.log"
    conda:
        "envs/cage.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.ctss}) $(dirname {output.raw_bw}) $(dirname {log})

        tmp_dir=$(mktemp -d)
        trap 'rm -rf "$tmp_dir"' EXIT

        five_prime_all="$tmp_dir/5prime_all.bed"
        five_prime_unique="$tmp_dir/5prime_unique.bed"

        # Create 1bp, strand-aware intervals at the most 5' nucleotide of each read.
        samtools view -@ {threads} -F {params.exclude_flags} {input.bam} | \
        awk 'BEGIN{{OFS="\t"}}
            function ref_len(cigar,   len, n, op, chunk) {{
                len=0;
                while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {{
                    chunk=substr(cigar, RSTART, RLENGTH);
                    op=substr(chunk, length(chunk), 1);
                    n=substr(chunk, 1, length(chunk)-1) + 0;
                    if (op ~ /[MDN=X]/) len += n;
                    cigar=substr(cigar, RLENGTH+1);
                }}
                return len;
            }}
        {{
            chrom=$3;
            start=$4-1;
            strand="+";
            if (int($2/16)%2) strand="-";
            aln_len=ref_len($6);
            if (aln_len<1) next;
            if (strand=="+") {{s=start; e=start+1;}}
            else {{s=start+aln_len-1; e=start+aln_len;}}
            if (s>=0) print chrom, s, e, ".", 1, strand;
        }}' | sort -k1,1 -k2,2n -k3,3n -k6,6 > "$five_prime_all"

        cut -f1-3,6 "$five_prime_all" | awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,".",0,$4}}' | \
            sort -k1,1 -k2,2n -k3,3n -k6,6 -u > "$five_prime_unique"

        # Count reads per 5' site using bedtools coverage -counts
        bedtools coverage -a "$five_prime_unique" -b "$five_prime_all" -s -counts > "$tmp_dir/ctss_cov.bed"

        awk 'BEGIN{{OFS="\t"}} {{
            chr=$1; start=$2; end=$3; strand=$6; count=$7;
            pos1=start+1;
            print chr, pos1, strand, count;
            print chr, start, end, count > "{output.raw_bg}";
            total+=count;
        }} END{{print total > "'$tmp_dir'/total.txt"}}' "$tmp_dir/ctss_cov.bed" > {output.ctss}

        total=$(cat "$tmp_dir/total.txt")
        if [ "$total" -eq 0 ]; then
            awk 'BEGIN{{OFS="\t"}} {{print $1,0,1,0}}' {input.chrom_sizes} > {output.cpm_bg}
        else
            awk -v total="$total" 'BEGIN{{OFS="\t"}} {{cpm=($4*1000000)/total; print $1,$2,$3,cpm}}' \
                {output.raw_bg} > {output.cpm_bg}
        fi

        sort -k1,1 -k2,2n {output.raw_bg} -o {output.raw_bg}
        sort -k1,1 -k2,2n {output.cpm_bg} -o {output.cpm_bg}

        bedGraphToBigWig {output.raw_bg} {input.chrom_sizes} {output.raw_bw} >> {log} 2>&1
        bedGraphToBigWig {output.cpm_bg} {input.chrom_sizes} {output.cpm_bw} >> {log} 2>&1
        """


rule cager_cluster:
    """Use CAGEr to build CAGE tag clusters from CTSS."""
    input:
        ctss=f"{OUTDIR}/ctss/{{sample}}/{{sample}}.ctss.tsv",
        chrom_sizes=config["reference"]["chrom_sizes"]
    output:
        tag_clusters=f"{OUTDIR}/cager/{{sample}}/{{sample}}_tagClusters.bed"
    params:
        sample="{sample}",
        min_count=config.get("cager", {}).get("min_count", 1),
        max_dist=config.get("cager", {}).get("max_dist", 20)
    log:
        "logs/cager/{sample}.log"
    conda:
        "envs/cager.yaml"
    script:
        "scripts/run_cager.R"
