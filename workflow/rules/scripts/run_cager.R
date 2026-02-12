suppressPackageStartupMessages({
  library(CAGEr)
})

ctss_file <- snakemake@input[["ctss"]]
out_file <- snakemake@output[["tag_clusters"]]
sample_name <- snakemake@params[["sample"]]
min_count <- as.integer(snakemake@params[["min_count"]])
max_dist <- as.integer(snakemake@params[["max_dist"]])

ctss <- read.table(
  ctss_file,
  sep = "\t",
  header = FALSE,
  col.names = c("chr", "pos", "strand", "counts"),
  stringsAsFactors = FALSE
)

if (nrow(ctss) == 0) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  file.create(out_file)
  quit(save = "no", status = 0)
}

ctss_gr <- GenomicRanges::GRanges(
  seqnames = ctss$chr,
  ranges = IRanges::IRanges(start = ctss$pos, end = ctss$pos),
  strand = ctss$strand,
  score = ctss$counts
)

ce <- if (exists("CAGEset", where = asNamespace("CAGEr"), mode = "function")) {
  CAGEset(genomeName = "unknown")
} else if (exists("CAGEexp", where = asNamespace("CAGEr"), mode = "function")) {
  CAGEexp(genomeName = "unknown")
} else {
  stop("Neither CAGEset() nor CAGEexp() is available in the installed CAGEr version")
}

set_ctss_success <- FALSE

if (methods::is(ce, "CAGEexp") &&
    exists("CTSStagCountSE", where = asNamespace("CAGEr"), mode = "function")) {
  CTSStagCountSE(ce, sample_name) <- ctss_gr
  set_ctss_success <- TRUE
}

if (!set_ctss_success && exists("CTSStagCount", where = asNamespace("CAGEr"), mode = "function")) {
  CTSStagCount(ce, sample_name) <- ctss_gr
  set_ctss_success <- TRUE
}

if (!set_ctss_success) {
  stop("Unable to assign CTSS counts: neither CTSStagCountSE nor CTSStagCount worked")
}

sampleLabels(ce) <- sample_name

clusterCTSS(
  object = ce,
  threshold = min_count,
  nrPassThreshold = 1,
  method = "distclu",
  maxDist = max_dist,
  removeSingletons = FALSE,
  keepSingletonsAbove = min_count
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
exportToBed(
  object = ce,
  what = "tagClusters",
  qLow = 0,
  qUp = 1,
  oneFile = TRUE,
  file = out_file,
  score = "raw"
)
