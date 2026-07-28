#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: prepare_pilot_data.R <microeco.tar.gz> <output-dir> <summary.json>")
}

source_tar <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
summary_path <- args[[3]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

required <- c(
  "microeco/data/otu_table_16S.RData",
  "microeco/data/taxonomy_table_16S.RData",
  "microeco/data/sample_info_16S.RData",
  "microeco/data/env_data_16S.RData",
  "microeco/data/phylo_tree_16S.RData"
)

extract_dir <- tempfile("microeco-data-")
dir.create(extract_dir)
on.exit(unlink(extract_dir, recursive = TRUE, force = TRUE), add = TRUE)
utils::untar(source_tar, files = required, exdir = extract_dir)

data_env <- new.env(parent = emptyenv())
for (relative_path in required) {
  extracted <- file.path(extract_dir, relative_path)
  if (!file.exists(extracted)) {
    stop(sprintf("Required package data is missing: %s", relative_path))
  }
  load(extracted, envir = data_env)
}

otu <- get("otu_table_16S", envir = data_env)
taxonomy <- get("taxonomy_table_16S", envir = data_env)
metadata <- get("sample_info_16S", envir = data_env)
environment <- get("env_data_16S", envir = data_env)
phylo_tree <- get("phylo_tree_16S", envir = data_env)

feature_ids_match <- identical(rownames(otu), rownames(taxonomy))
sample_ids_match <- identical(colnames(otu), rownames(metadata))
tree_contains_all_features <- all(rownames(otu) %in% phylo_tree$tip.label)
tree_extra_tips <- sum(!phylo_tree$tip.label %in% rownames(otu))
if (!feature_ids_match || !sample_ids_match || !tree_contains_all_features) {
  stop("microeco feature or sample identifiers do not align")
}

write_tsv <- function(x, id_name, path) {
  ids <- rownames(x)
  if (id_name %in% colnames(x) && identical(as.character(x[[id_name]]), ids)) {
    out <- x
  } else {
    out <- data.frame(ids, x, check.names = FALSE)
    colnames(out)[1] <- id_name
  }
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

write_tsv(otu, "FeatureID", file.path(output_dir, "otutab.tsv"))
write_tsv(taxonomy, "FeatureID", file.path(output_dir, "taxonomy.tsv"))
write_tsv(metadata, "SampleID", file.path(output_dir, "metadata.tsv"))
write_tsv(environment, "SampleID", file.path(output_dir, "environment.tsv"))

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required to export phylo_tree_16S")
}
tree_text <- paste0(ape::write.tree(phylo_tree), "\n")
tree_path <- file.path(output_dir, "rooted-tree.nwk.gz")
tree_connection <- gzfile(tree_path, open = "wb", compression = 9)
writeLines(tree_text, tree_connection, useBytes = TRUE)
close(tree_connection)

common_environment_samples <- sum(rownames(metadata) %in% rownames(environment))
json_lines <- c(
  "{",
  sprintf('  "features": %d,', nrow(otu)),
  sprintf('  "samples": %d,', ncol(otu)),
  sprintf('  "taxonomy_features": %d,', nrow(taxonomy)),
  sprintf('  "metadata_samples": %d,', nrow(metadata)),
  sprintf('  "environment_samples": %d,', nrow(environment)),
  sprintf('  "environment_variables": %d,', ncol(environment)),
  sprintf('  "common_environment_samples": %d,', common_environment_samples),
  sprintf('  "tree_tips": %d,', ape::Ntip(phylo_tree)),
  sprintf('  "tree_internal_nodes": %d,', ape::Nnode(phylo_tree)),
  sprintf('  "tree_rooted": %s,', tolower(as.character(ape::is.rooted(phylo_tree)))),
  sprintf('  "tree_contains_all_features": %s,', tolower(as.character(tree_contains_all_features))),
  sprintf('  "tree_extra_tips": %d,', tree_extra_tips),
  sprintf('  "feature_ids_match": %s,', tolower(as.character(feature_ids_match))),
  sprintf('  "sample_ids_match": %s', tolower(as.character(sample_ids_match))),
  "}"
)
writeLines(json_lines, summary_path, useBytes = TRUE)

message(sprintf(
  "Exported %d features across %d samples; %d metadata samples overlap environment data.",
  nrow(otu), ncol(otu), common_environment_samples
))
