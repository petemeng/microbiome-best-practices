#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

required_packages <- c("digest", "jsonlite", "phyloseq", "reconsi")
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

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "data/small/absolute-quantification"
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(reconsi))
data_object <- get("Vandeputte", pos = "package:reconsi")
counts_samples_by_taxa <- as(
  phyloseq::otu_table(data_object),
  "matrix"
)
stopifnot(!phyloseq::taxa_are_rows(data_object))
counts <- t(counts_samples_by_taxa)
storage.mode(counts) <- "integer"

metadata_source <- as(
  phyloseq::sample_data(data_object),
  "data.frame"
)
metadata <- data.frame(
  SubjectID = as.character(metadata_source$ID),
  Cohort = as.character(metadata_source$Cohort),
  HealthStatus = as.character(metadata_source$Health.status),
  stringsAsFactors = FALSE,
  row.names = rownames(metadata_source),
  check.names = FALSE
)
cell_load <- data.frame(
  AbsoluteCellLoad = as.numeric(metadata_source$absCountFrozen),
  Unit = "cells per gram frozen feces",
  stringsAsFactors = FALSE,
  row.names = rownames(metadata_source),
  check.names = FALSE
)

taxonomy <- data.frame(
  Kingdom = NA_character_,
  Phylum = NA_character_,
  Class = NA_character_,
  Order = NA_character_,
  Family = NA_character_,
  Genus = rownames(counts),
  Species = NA_character_,
  stringsAsFactors = FALSE,
  row.names = rownames(counts),
  check.names = FALSE
)

stopifnot(
  identical(colnames(counts), rownames(metadata)),
  identical(colnames(counts), rownames(cell_load)),
  identical(rownames(counts), rownames(taxonomy)),
  nrow(counts) == 234L,
  ncol(counts) == 135L,
  sum(counts) == 4080996L,
  all(counts >= 0L),
  all(is.finite(cell_load$AbsoluteCellLoad)),
  all(cell_load$AbsoluteCellLoad > 0)
)

write_keyed_tsv <- function(x, id_name, path) {
  out <- data.frame(
    setNames(list(rownames(x)), id_name),
    x,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(
    out,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

paths <- c(
  otutab = file.path(output_dir, "otutab.tsv"),
  taxonomy = file.path(output_dir, "taxonomy.tsv"),
  metadata = file.path(output_dir, "metadata.tsv"),
  cell_load = file.path(output_dir, "cell-load.tsv")
)
write_keyed_tsv(counts, "FeatureID", paths[["otutab"]])
write_keyed_tsv(taxonomy, "FeatureID", paths[["taxonomy"]])
write_keyed_tsv(metadata, "SampleID", paths[["metadata"]])
write_keyed_tsv(cell_load, "SampleID", paths[["cell_load"]])

package_data_files <- c(
  data_index = system.file("data", "Rdata.rds", package = "reconsi"),
  data_database = system.file("data", "Rdata.rdb", package = "reconsi")
)
sha256 <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

summary_payload <- list(
  dataset = "Vandeputte",
  description = paste(
    "Crohn's disease and healthy-control genus counts with parallel",
    "flow-cytometry cell loads"
  ),
  primary_publication = list(
    citation = paste(
      "Vandeputte D, Kathagen G, D'hoe K, et al.",
      "Nature. 2017;551:507-511."
    ),
    doi = "10.1038/nature24460",
    url = "https://www.nature.com/articles/nature24460"
  ),
  package = list(
    name = "reconsi",
    version = as.character(utils::packageVersion("reconsi")),
    bioconductor_release = "3.19",
    source_url = paste0(
      "https://bioconductor.org/packages/3.19/bioc/src/contrib/",
      "reconsi_1.16.0.tar.gz"
    ),
    license = utils::packageDescription("reconsi", fields = "License")
  ),
  dimensions = list(
    features = nrow(counts),
    samples = ncol(counts),
    reads = unname(sum(counts)),
    crohns_samples = unname(sum(metadata$HealthStatus == "CD")),
    healthy_samples = unname(sum(metadata$HealthStatus == "Healthy"))
  ),
  package_data_files = lapply(package_data_files, function(path) {
    list(
      file = basename(path),
      bytes = unname(file.info(path)$size),
      sha256 = sha256(path)
    )
  }),
  exported_files = lapply(paths, function(path) {
    list(
      file = basename(path),
      bytes = unname(file.info(path)$size),
      sha256 = sha256(path)
    )
  })
)
jsonlite::write_json(
  summary_payload,
  path = file.path(output_dir, "source-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Prepared Vandeputte absolute-quantification data: ",
  nrow(counts), " genera x ", ncol(counts), " samples; ",
  format(sum(counts), big.mark = ","), " reads."
)
