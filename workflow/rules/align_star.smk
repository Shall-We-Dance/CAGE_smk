# rules/align_star.smk

# ----------------------------------------------------------------------------
# Alignment Rules (STAR)
# ----------------------------------------------------------------------------

if not is_single_end():
    rule star_pe:
        input:
            r1=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz",
            r2=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R2.fastq.gz"
        params:
            index=config["reference"]["star_index"],
            extra=config.get("star", {}).get("extra", "")
        output:
            bam=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}_pe.bam"),
            bai=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}_pe.bam.bai"),
            flagstat=aligned_flagstat_path("{sample}")
        log:
            "logs/star/{sample}_pe.log"
        threads: int(config["threads"]["star"])
        conda:
            "envs/star.yaml"
        shell:
            r"""
            set -euo pipefail
            sample_dir="{OUTDIR}/star/{wildcards.sample}"
            tmp_prefix="$sample_dir/{wildcards.sample}."
            mkdir -p "$sample_dir" $(dirname {output.flagstat}) $(dirname {log})

            STAR \
                --runThreadN {threads} \
                --genomeDir {params.index} \
                --readFilesIn {input.r1} {input.r2} \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix "$tmp_prefix" \
                {params.extra} \
                > {log} 2>&1

            mv "$tmp_prefix"Aligned.sortedByCoord.out.bam {output.bam}
            samtools index -@ {threads} {output.bam}
            samtools flagstat {output.bam} > {output.flagstat}
            """
else:
    rule star_se:
        input:
            r1=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz"
        params:
            index=config["reference"]["star_index"],
            extra=config.get("star", {}).get("extra", "")
        output:
            bam=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}_se.bam"),
            bai=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}_se.bam.bai"),
            flagstat=aligned_flagstat_path("{sample}")
        log:
            "logs/star/{sample}_se.log"
        threads: int(config["threads"]["star"])
        conda:
            "envs/star.yaml"
        shell:
            r"""
            set -euo pipefail
            sample_dir="{OUTDIR}/star/{wildcards.sample}"
            tmp_prefix="$sample_dir/{wildcards.sample}."
            mkdir -p "$sample_dir" $(dirname {output.flagstat}) $(dirname {log})

            STAR \
                --runThreadN {threads} \
                --genomeDir {params.index} \
                --readFilesIn {input.r1} \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix "$tmp_prefix" \
                {params.extra} \
                > {log} 2>&1

            mv "$tmp_prefix"Aligned.sortedByCoord.out.bam {output.bam}
            samtools index -@ {threads} {output.bam}
            samtools flagstat {output.bam} > {output.flagstat}
            """


# ----------------------------------------------------------------------------
# Filtering Rules
# ----------------------------------------------------------------------------

rule filter_unique_mappers:
    """
    Keep uniquely mapped reads (MAPQ >= threshold; STAR unique default MAPQ=255).
    """
    input:
        bam=lambda wc: aligned_bam_path(wc.sample)
    output:
        bam=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.bam"),
        bai=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.bam.bai"),
        flagstat=unique_flagstat_path("{sample}")
    params:
        mapq=config.get("star", {}).get("mapq_threshold", 255),
        exclude_flags=config.get("samtools_exclude_flags", 1804)
    log:
        f"logs/star/{{sample}}.unique.log"
    threads: int(config["threads"]["samtools"])
    conda:
        "envs/star.yaml"
    shell:
        r"""
        set -euo pipefail
        samtools view -b -h -q {params.mapq} -F {params.exclude_flags} \
            -@ {threads} {input.bam} | \
        samtools sort -@ {threads} -o {output.bam} -

        samtools index -@ {threads} {output.bam}

        echo "Original BAM:" > {log}
        samtools flagstat {input.bam} >> {log}
        echo -e "\nFiltered BAM (unique mappers):" >> {log}
        samtools flagstat {output.bam} >> {log}
        samtools flagstat {output.bam} > {output.flagstat}
        """


if FILTER_BLACKLIST:
    rule filter_blacklist_bam:
        """
        Remove reads overlapping blacklist regions.
        """
        input:
            bam=unique_bam_path("{sample}"),
            blacklist=blacklist_path()
        output:
            bam=maybe_temp(filtered_bam_path("{sample}")),
            bai=maybe_temp(filtered_bam_path("{sample}") + ".bai"),
            flagstat=filtered_flagstat_path("{sample}")
        params:
            nonamecheck=config.get("bedtools_nonamecheck", True)
        threads: int(config["threads"]["samtools"])
        log:
            "logs/star/{sample}.filter_blacklist.log"
        conda:
            "envs/star.yaml"
        shell:
            r"""
            set -euo pipefail
            echo "Input BAM:" > {log}
            samtools flagstat {input.bam} >> {log}

            nonamecheck_opt=""
            if [ "{params.nonamecheck}" = "True" ] && \
                bedtools intersect -h 2>&1 | grep -q -- '--nonamecheck'; then
                nonamecheck_opt="--nonamecheck"
            fi

            if [ ! -s "{input.blacklist}" ]; then
                echo -e "\nBlacklist BED is missing or empty; skipping blacklist filtering." >> {log}
                samtools view -b -@ {threads} {input.bam} -o {output.bam}
            else
                bedtools intersect -v \
                    -abam {input.bam} \
                    -b {input.blacklist} \
                    $nonamecheck_opt 2>> {log} | \
                samtools sort -@ {threads} -o {output.bam} -
            fi

            samtools index -@ {threads} {output.bam} 2>> {log}
            echo -e "\nAfter blacklist filtering:" >> {log}
            samtools flagstat {output.bam} >> {log}
            samtools flagstat {output.bam} > {output.flagstat}
            """


if REMOVE_DUPLICATES:
    rule remove_duplicates:
        """
        Remove PCR duplicates using Picard or samtools.
        """
        input:
            bam=lambda wc: filtered_bam_path(wc.sample) if FILTER_BLACKLIST else unique_bam_path(wc.sample)
        output:
            bam=maybe_temp(dedup_bam_path("{sample}")),
            bai=maybe_temp(dedup_bam_path("{sample}") + ".bai"),
            metrics=f"{OUTDIR}/star/{{sample}}/{{sample}}.dedup_metrics.txt",
            flagstat=dedup_flagstat_path("{sample}")
        params:
            method=config.get("dedup_method", "picard"),
            picard_java_opts=config.get("picard_java_opts", "-Xmx4g")
        threads: int(config["threads"]["samtools"])
        log:
            "logs/star/{sample}.dedup.log"
        conda:
            "envs/star.yaml"
        shell:
            r"""
            set -euo pipefail
             mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat}) $(dirname {log})

            if [ "{params.method}" = "picard" ]; then
                set +e
                picard --java-options "{params.picard_java_opts}" MarkDuplicates \
                    I={input.bam} \
                    O={output.bam} \
                    M={output.metrics} \
                    REMOVE_DUPLICATES=true \
                    VALIDATION_STRINGENCY=LENIENT \
                    > {log} 2>&1
                picard_status=$?
                set -e

                if [ $picard_status -ne 0 ]; then
                    echo "Picard MarkDuplicates failed (exit code: $picard_status). Falling back to samtools markdup." >> {log}
                    samtools markdup -r -@ {threads} {input.bam} {output.bam} >> {log} 2>&1
                    samtools flagstat {output.bam} > {output.metrics}
                fi
            else
                samtools markdup -r -@ {threads} {input.bam} {output.bam} 2> {log}
                samtools flagstat {output.bam} > {output.metrics}
            fi

            samtools index -@ {threads} {output.bam}
            samtools flagstat {output.bam} > {output.flagstat}
            """