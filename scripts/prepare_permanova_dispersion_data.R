#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    paste(
      "Usage: prepare_permanova_dispersion_data.R",
      "<soilrep.RData|installed> <output-dir> <summary-json>"
    ),
    call. = FALSE
  )
}

required_packages <- c("digest", "jsonlite", "phyloseq")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

expected_phyloseq_version <- "1.48.0"
observed_phyloseq_version <- as.character(
  utils::packageVersion("phyloseq")
)
stopifnot(identical(
  observed_phyloseq_version,
  expected_phyloseq_version
))

source_arg <- args[[1L]]
if (identical(source_arg, "installed")) {
  source_arg <- system.file(
    "data",
    "soilrep.RData",
    package = "phyloseq"
  )
}
source_rdata <- normalizePath(source_arg, mustWork = TRUE)
output_dir <- args[[2L]]
summary_path <- args[[3L]]

expected_source_sha256 <- paste0(
  "47b25a0a033b794fca1db30ae6fdb68",
  "c55a44e6c4cd0d1be76788d1d77b51f36"
)
observed_source_sha256 <- digest::digest(
  file = source_rdata,
  algo = "sha256",
  serialize = FALSE
)
stopifnot(identical(
  observed_source_sha256,
  expected_source_sha256
))

data_env <- new.env(parent = emptyenv())
loaded_names <- load(source_rdata, envir = data_env)
stopifnot("soilrep" %in% loaded_names)
ps <- data_env$soilrep
stopifnot(inherits(ps, "phyloseq"))

otutab <- as(phyloseq::otu_table(ps), "matrix")
if (!phyloseq::taxa_are_rows(ps)) {
  otutab <- t(otutab)
}
metadata_source <- as.data.frame(phyloseq::sample_data(ps))
metadata_source <- data.frame(
  lapply(metadata_source, as.character),
  row.names = rownames(metadata_source),
  check.names = FALSE
)

stopifnot(
  identical(colnames(otutab), rownames(metadata_source)),
  nrow(otutab) == 16825L,
  ncol(otutab) == 56L,
  sum(otutab) == 98022,
  all(is.finite(otutab)),
  all(otutab >= 0),
  all(otutab == floor(otutab)),
  all(rowSums(otutab) > 0),
  all(colSums(otutab) > 0),
  all(c("Treatment", "warmed", "clipped", "Sample") %in%
    colnames(metadata_source))
)

biological_id <- metadata_source$Sample
block_id <- sub("^([0-9]+).*$", "\\1", biological_id)
warming <- ifelse(
  metadata_source$warmed == "yes",
  "Warmed",
  "Ambient"
)
clipping <- ifelse(
  metadata_source$clipped == "yes",
  "Clipped",
  "Unclipped"
)
treatment_label <- paste(warming, clipping, sep = " + ")
main_plot_id <- paste0("Block ", block_id, " / ", warming)
library_count <- unname(table(biological_id)[biological_id])
technical_index <- ave(
  seq_along(biological_id),
  biological_id,
  FUN = seq_along
)

metadata <- data.frame(
  SourceTreatment = metadata_source$Treatment,
  SourceWarmed = metadata_source$warmed,
  SourceClipped = metadata_source$clipped,
  BiologicalSampleID = biological_id,
  Block = paste("Block", block_id),
  MainPlot = main_plot_id,
  Warming = warming,
  Clipping = clipping,
  Treatment = treatment_label,
  TechnicalLibraryIndex = as.integer(technical_index),
  TechnicalLibraryCount = as.integer(library_count),
  row.names = rownames(metadata_source),
  check.names = FALSE
)

treatment_levels <- c(
  "Ambient + Unclipped",
  "Ambient + Clipped",
  "Warmed + Unclipped",
  "Warmed + Clipped"
)
stopifnot(
  length(unique(metadata$BiologicalSampleID)) == 24L,
  length(unique(metadata$Block)) == 6L,
  length(unique(metadata$MainPlot)) == 12L,
  setequal(unique(metadata$Treatment), treatment_levels),
  identical(sort(as.integer(table(metadata$SourceTreatment))), rep(14L, 4L)),
  identical(
    unname(as.integer(table(table(metadata$BiologicalSampleID)))),
    c(2L, 12L, 10L)
  )
)

taxonomy <- data.frame(
  Kingdom = rep(NA_character_, nrow(otutab)),
  Phylum = rep(NA_character_, nrow(otutab)),
  Class = rep(NA_character_, nrow(otutab)),
  Order = rep(NA_character_, nrow(otutab)),
  Family = rep(NA_character_, nrow(otutab)),
  Genus = rep(NA_character_, nrow(otutab)),
  Species = rep(NA_character_, nrow(otutab)),
  TaxonomyStatus = rep(
    "Not distributed with phyloseq::soilrep; not used in this analysis",
    nrow(otutab)
  ),
  row.names = rownames(otutab),
  check.names = FALSE
)

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

otutab_path <- file.path(output_dir, "otutab.tsv")
taxonomy_path <- file.path(output_dir, "taxonomy.tsv")
metadata_path <- file.path(output_dir, "metadata.tsv")
write_tsv(otutab, "FeatureID", otutab_path)
write_tsv(taxonomy, "FeatureID", taxonomy_path)
write_tsv(metadata, "LibraryID", metadata_path)

output_sha256 <- vapply(
  c(otutab_path, taxonomy_path, metadata_path),
  function(path) {
    digest::digest(
      file = path,
      algo = "sha256",
      serialize = FALSE
    )
  },
  character(1)
)
names(output_sha256) <- basename(names(output_sha256))

replicate_distribution <- table(table(metadata$BiologicalSampleID))
summary <- list(
  source = "phyloseq::soilrep",
  source_file = "phyloseq/data/soilrep.RData",
  phyloseq_version = observed_phyloseq_version,
  source_sha256 = observed_source_sha256,
  source_bytes = unname(file.info(source_rdata)$size),
  primary_study = list(
    title = paste(
      "Reproducibility and quantitation of amplicon",
      "sequencing-based detection"
    ),
    citation = paste(
      "Zhou J, Wu L, Deng Y, et al.",
      "ISME Journal. 2011;5:1303-1313."
    ),
    doi = "10.1038/ismej.2011.11"
  ),
  license = "phyloseq package AGPL-3",
  features = unname(nrow(otutab)),
  technical_libraries = unname(ncol(otutab)),
  biological_samples = unname(length(unique(metadata$BiologicalSampleID))),
  blocks = unname(length(unique(metadata$Block))),
  main_plots = unname(length(unique(metadata$MainPlot))),
  treatments = unname(length(unique(metadata$Treatment))),
  reads = unname(sum(otutab)),
  minimum_library_size = unname(min(colSums(otutab))),
  median_library_size = unname(stats::median(colSums(otutab))),
  maximum_library_size = unname(max(colSums(otutab))),
  biological_samples_with_1_library = unname(
    as.integer(replicate_distribution[["1"]])
  ),
  biological_samples_with_2_libraries = unname(
    as.integer(replicate_distribution[["2"]])
  ),
  biological_samples_with_3_libraries = unname(
    as.integer(replicate_distribution[["3"]])
  ),
  taxonomy_status = paste(
    "The packaged soilrep object has no tax_table; seven rank columns",
    "are deliberately empty and TaxonomyStatus records this boundary."
  ),
  output_sha256 = as.list(output_sha256)
)
jsonlite::write_json(
  summary,
  summary_path,
  auto_unbox = TRUE,
  pretty = TRUE
)

message(
  sprintf(
    paste(
      "Exported %s features x %s technical libraries (%s biological",
      "samples; %s blocks)."
    ),
    format(nrow(otutab), big.mark = ","),
    ncol(otutab),
    length(unique(metadata$BiologicalSampleID)),
    length(unique(metadata$Block))
  )
)
