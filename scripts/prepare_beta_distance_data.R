#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(digest)
  library(jsonlite)
  library(phyloseq)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value after ", flag)
  args[[hit + 1L]]
}

project_root <- normalizePath(
  arg_value("--project-root", "."),
  mustWork = TRUE
)
output_dir <- file.path(
  project_root,
  arg_value("--output-dir", "data/small/beta-distances")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data("GlobalPatterns", package = "phyloseq", envir = environment())
stopifnot(identical(as.character(packageVersion("phyloseq")), "1.48.0"))
ps <- GlobalPatterns
source_features <- ntaxa(ps)
zero_sum_features <- taxa_names(ps)[taxa_sums(ps) == 0]
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

counts <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) counts <- t(counts)
storage.mode(counts) <- "integer"

taxonomy <- as(tax_table(ps), "matrix")
metadata <- data.frame(
  sample_data(ps),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
tree <- phy_tree(ps)

stopifnot(
  identical(rownames(counts), rownames(taxonomy)),
  setequal(colnames(counts), rownames(metadata)),
  setequal(rownames(counts), tree$tip.label),
  ape::is.rooted(tree),
  !is.null(tree$edge.length),
  all(is.finite(tree$edge.length)),
  all(tree$edge.length >= 0),
  all(counts >= 0),
  all(counts == floor(counts))
)

metadata <- metadata[colnames(counts), , drop = FALSE]
habitat_map <- c(
  "Feces" = "Human-associated",
  "Skin" = "Human-associated",
  "Tongue" = "Human-associated",
  "Freshwater" = "Freshwater",
  "Freshwater (creek)" = "Freshwater",
  "Ocean" = "Marine/estuarine",
  "Sediment (estuary)" = "Marine/estuarine",
  "Soil" = "Soil",
  "Mock" = "Mock"
)
metadata$SampleType <- as.character(metadata$SampleType)
metadata$HabitatClass <- unname(habitat_map[metadata$SampleType])
stopifnot(!anyNA(metadata$HabitatClass))

otutab_out <- tibble::as_tibble(counts, rownames = "FeatureID")
taxonomy_out <- tibble::as_tibble(taxonomy, rownames = "FeatureID")
metadata_out <- data.frame(
  SampleID = rownames(metadata),
  SampleType = as.character(metadata[, "SampleType"]),
  HabitatClass = as.character(metadata[, "HabitatClass"]),
  Primer = as.character(metadata[, "Primer"]),
  Final_Barcode = as.character(metadata[, "Final_Barcode"]),
  Barcode_truncated_plus_T = as.character(
    metadata[, "Barcode_truncated_plus_T"]
  ),
  Barcode_full_length = as.character(metadata[, "Barcode_full_length"]),
  Description = as.character(metadata[, "Description"]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(
  nrow(metadata_out) == ncol(counts),
  identical(metadata_out$SampleID, colnames(counts)),
  !anyDuplicated(metadata_out$SampleID)
)

paths <- c(
  otutab = file.path(output_dir, "otutab.tsv"),
  taxonomy = file.path(output_dir, "taxonomy.tsv"),
  metadata = file.path(output_dir, "metadata.tsv"),
  tree = file.path(output_dir, "rooted-tree.nwk.gz")
)

readr::write_tsv(otutab_out, paths[["otutab"]])
readr::write_tsv(taxonomy_out, paths[["taxonomy"]], na = "")
readr::write_tsv(metadata_out, paths[["metadata"]], na = "")

tree_connection <- gzfile(paths[["tree"]], open = "wt", compression = 9)
writeLines(ape::write.tree(tree), tree_connection, useBytes = TRUE)
close(tree_connection)

source_rdata <- system.file(
  "data", "GlobalPatterns.RData", package = "phyloseq"
)
stopifnot(nzchar(source_rdata), file.exists(source_rdata))

file_summary <- lapply(paths, function(path) {
  list(
    path = file.path("data", "small", "beta-distances", basename(path)),
    bytes = unname(file.info(path)$size),
    sha256 = digest::digest(path, algo = "sha256", file = TRUE)
  )
})

summary <- list(
  dataset = "GlobalPatterns",
  source = "phyloseq 1.48.0 packaged data",
  source_file = list(
    path = "phyloseq/data/GlobalPatterns.RData",
    bytes = unname(file.info(source_rdata)$size),
    sha256 = digest::digest(source_rdata, algo = "sha256", file = TRUE)
  ),
  primary_study = list(
    citation = paste(
      "Caporaso JG et al. Global patterns of 16S rRNA diversity at a",
      "depth of millions of sequences per sample. PNAS. 2011;108:4516-4522."
    ),
    doi = "10.1073/pnas.1000080107",
    pmcid = "PMC3063599"
  ),
  package = list(
    name = "phyloseq",
    version = as.character(packageVersion("phyloseq")),
    license = packageDescription("phyloseq")$License
  ),
  dimensions = list(
    source_features = source_features,
    zero_sum_features_removed = length(zero_sum_features),
    features = nrow(counts),
    samples = ncol(counts),
    reads = unname(sum(counts)),
    min_depth = unname(min(colSums(counts))),
    median_depth = unname(stats::median(colSums(counts))),
    max_depth = unname(max(colSums(counts))),
    zero_fraction = unname(mean(counts == 0)),
    sample_types = length(unique(metadata$SampleType)),
    habitat_classes = length(unique(metadata$HabitatClass))
  ),
  phylogeny = list(
    rooted = ape::is.rooted(tree),
    tips = length(tree$tip.label),
    internal_nodes = tree$Nnode,
    branch_lengths_present = !is.null(tree$edge.length),
    negative_branch_lengths = sum(tree$edge.length < 0),
    tip_set_matches_counts = setequal(tree$tip.label, rownames(counts))
  ),
  files = file_summary
)

jsonlite::write_json(
  summary,
  file.path(output_dir, "source-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

message(
  sprintf(
    "Prepared GlobalPatterns beta-distance bundle: %s features x %s samples; %s reads.",
    format(nrow(counts), big.mark = ","),
    format(ncol(counts), big.mark = ","),
    format(sum(counts), big.mark = ",")
  )
)
