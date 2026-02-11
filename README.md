# CAGE Snakemake Pipeline

**CAGE** workflow, with **STAR** as the aligner and **single-end (SE)** mode as the default.

## Pipeline Overview

1. **FASTQ QC and trimming** with `fastp`
2. **Read mapping** with `STAR`
3. **Unique-mapper filtering** using `samtools view -q` (default MAPQ >= 255)
4. *(Optional)* **Blacklist filtering** with `bedtools intersect -v`
5. *(Optional)* **PCR duplicate removal** with `picard MarkDuplicates` or `samtools markdup`
6. **5' end (CTSS) counting**:
   - Extract the most 5' base (1 bp interval) from each mapped read
   - Count per-site signal using **bedtools coverage (2.27.1)**
7. **CAGE cluster construction** with **CAGEr (1.20.0)**
8. **Visualization tracks**: generate 5' end `raw` and `CPM` bigWig files

## Repository Structure

- `workflow/Snakefile`: Main workflow entry point
- `workflow/rules/align_star.smk`: STAR alignment and BAM filtering
- `workflow/rules/cage_ctss.smk`: CTSS counting, 5' bigWig generation, CAGEr clustering
- `workflow/rules/scripts/run_cager.R`: CAGEr execution script
- `workflow/rules/envs/*.yaml`: Conda environments
- `config.yaml`: Parameters and sample definitions

## Requirements

Recommended:

- `snakemake >= 9.3.3`
- `mamba` or `conda`

## Configuration

Edit `config.yaml` before running:

- `read_type`: defaults to `SE`
- `reference.star_index`: STAR genome index directory
- `reference.chrom_sizes`: chromosome sizes file (required by `bedGraphToBigWig`)
- `samples`: sample FASTQ paths
- `star.extra`: extra STAR arguments (default includes `--outFilterMultimapNmax 1 --alignIntronMax 1`)
- `cager.min_count` / `cager.max_dist`: CAGEr clustering parameters

## Run

```bash
snakemake -s workflow/Snakefile --use-conda -j 8
```

## Key Outputs

For sample `sampleA`:

- `results/ctss/sampleA/sampleA.ctss.tsv`
  - Columns: `chr`, `pos (1-based)`, `strand`, `count`
- `results/bigwig/sampleA/sampleA.5prime.raw.bw`
- `results/bigwig/sampleA/sampleA.5prime.cpm.bw`
- `results/cager/sampleA/sampleA_tagClusters.bed`
- `results/qc/multiqc/multiqc_report.html`

## Notes and Caveats

1. **bedGraphToBigWig dependency**: `ucsc-bedgraphtobigwig` is included in `envs/cage.yaml`, but package naming/availability may differ across platforms.
2. **STAR parameter tuning**: `--alignIntronMax 1` is suitable for many TSS-focused CAGE use cases, but should be adjusted if your library/species requires intron-spanning alignments.
3. **CAGEr memory usage**: very large samples may require increased memory and/or splitting strategies.
4. **SE/PE compatibility**: workflow defaults to SE; for PE set `read_type: "PE"` and provide `R2` for each sample.

## Potential Future Improvements

- Add strand-specific bigWig outputs (`plus`/`minus`) for clearer promoter-direction interpretation.
- Add cross-sample consensus cluster construction.
- Provide scheduler profiles (e.g., Slurm/SGE) for cluster deployment.
