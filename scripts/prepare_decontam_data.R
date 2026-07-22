#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    "Usage: prepare_decontam_data.R <MUClite.rds> <output-dir> <summary.json>",
    call. = FALSE
  )
}

expected_package_version <- "1.24.0"
observed_package_version <- as.character(
  utils::packageVersion("decontam")
)
stopifnot(identical(
  observed_package_version,
  expected_package_version
))

source_arg <- args[[1]]
if (identical(source_arg, "installed")) {
  source_arg <- system.file(
    "extdata",
    "MUClite.rds",
    package = "decontam"
  )
}
source_rds <- normalizePath(source_arg, mustWork = TRUE)
output_dir <- args[[2]]
summary_path <- args[[3]]

expected_sha256 <- paste0(
  "801456d5780d5f51308d04c23f447710",
  "2cc5fbe9543ffa62fad2240a99c43d62"
)
observed_sha256 <- digest::digest(
  file = source_rds,
  algo = "sha256",
  serialize = FALSE
)
stopifnot(identical(observed_sha256, expected_sha256))

ps <- readRDS(source_rds)
stopifnot(inherits(ps, "phyloseq"))

otutab <- as(phyloseq::otu_table(ps), "matrix")
if (!phyloseq::taxa_are_rows(ps)) {
  otutab <- t(otutab)
}
taxonomy <- as(phyloseq::tax_table(ps), "matrix")
metadata <- as.data.frame(phyloseq::sample_data(ps))

stopifnot(
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata)),
  all(otutab >= 0),
  all(otutab == floor(otutab)),
  all(c("quant_reading", "Sample_or_Control") %in% colnames(metadata))
)

metadata$IsNegative <- metadata$Sample_or_Control == "Control Sample"

write_tsv <- function(x, id_name, path) {
  out <- data.frame(
    setNames(list(rownames(x)), id_name),
    x,
    check.names = FALSE
  )
  utils::write.table(
    out,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

write_tsv(
  otutab,
  "FeatureID",
  file.path(output_dir, "otutab.tsv")
)
write_tsv(
  taxonomy,
  "FeatureID",
  file.path(output_dir, "taxonomy.tsv")
)
write_tsv(
  metadata,
  "SampleID",
  file.path(output_dir, "metadata.tsv")
)

summary <- list(
  source = "decontam::MUClite.rds",
  decontam_version = observed_package_version,
  source_sha256 = observed_sha256,
  features = unname(nrow(otutab)),
  samples = unname(ncol(otutab)),
  biological_samples = unname(sum(!metadata$IsNegative)),
  negative_controls = unname(sum(metadata$IsNegative)),
  plates = unname(length(unique(metadata$PlateNumber))),
  minimum_library_size = unname(min(colSums(otutab))),
  median_library_size = unname(stats::median(colSums(otutab))),
  maximum_library_size = unname(max(colSums(otutab)))
)
jsonlite::write_json(
  summary,
  summary_path,
  auto_unbox = TRUE,
  pretty = TRUE
)

message(
  sprintf(
    "Exported %s features x %s samples (%s controls).",
    format(nrow(otutab), big.mark = ","),
    format(ncol(otutab), big.mark = ","),
    sum(metadata$IsNegative)
  )
)
