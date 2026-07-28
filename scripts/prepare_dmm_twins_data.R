#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value after ", flag, call. = FALSE)
  args[[hit + 1L]]
}

required_packages <- c("digest", "jsonlite", "DirichletMultinomial")
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

project_root <- normalizePath(
  arg_value("--project-root", "."),
  mustWork = TRUE
)
output_dir <- file.path(
  project_root,
  arg_value("--output-dir", "data/small/community-typing-dmm")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expected_version <- "1.46.0"
observed_version <- as.character(
  utils::packageVersion("DirichletMultinomial")
)
stopifnot(identical(observed_version, expected_version))

twins_path <- system.file(
  "extdata", "Twins.csv", package = "DirichletMultinomial"
)
phenotype_path <- system.file(
  "extdata", "TwinStudy.t", package = "DirichletMultinomial"
)
description_path <- system.file(
  "DESCRIPTION", package = "DirichletMultinomial"
)
stopifnot(
  nzchar(twins_path),
  nzchar(phenotype_path),
  nzchar(description_path),
  file.exists(twins_path),
  file.exists(phenotype_path),
  file.exists(description_path)
)

expected_source_sha256 <- c(
  "Twins.csv" = paste0(
    "b6a231acded317b49b7fd800bf8b4ce1",
    "2a75340c5b88499b51e92d72528f93ba"
  ),
  "TwinStudy.t" = paste0(
    "4286c3f9255f7f1ce7f5c5ac5e60d4e",
    "27284fe358ceb553108748c67ab4d8b05"
  )
)
observed_source_sha256 <- c(
  "Twins.csv" = digest::digest(
    twins_path, algo = "sha256", file = TRUE
  ),
  "TwinStudy.t" = digest::digest(
    phenotype_path, algo = "sha256", file = TRUE
  )
)
stopifnot(identical(observed_source_sha256, expected_source_sha256))

source_counts <- as.matrix(utils::read.csv(
  twins_path,
  row.names = 1L,
  check.names = FALSE
))
storage.mode(source_counts) <- "integer"
bmi_code <- scan(phenotype_path, what = integer(), quiet = TRUE)

stopifnot(
  nrow(source_counts) == 130L,
  ncol(source_counts) == 278L,
  length(bmi_code) == ncol(source_counts),
  all(is.finite(source_counts)),
  all(source_counts >= 0L),
  all(source_counts == floor(source_counts)),
  all(rowSums(source_counts) > 0L),
  all(colSums(source_counts) > 0L),
  sum(source_counts) == 570851L,
  identical(sort(unique(bmi_code)), 0:2),
  identical(unname(as.integer(table(bmi_code))), c(61L, 193L, 24L)),
  !anyDuplicated(rownames(source_counts)),
  !anyDuplicated(colnames(source_counts))
)

sample_ids <- colnames(source_counts)
participant_key <- sub("\\.2$", "", sample_ids)
has_dot2_suffix <- grepl("\\.2$", sample_ids)
samples_per_key <- unname(table(participant_key)[participant_key])
participant_tables <- split(has_dot2_suffix, participant_key)

stopifnot(
  length(unique(participant_key)) == 154L,
  sum(table(participant_key) == 2L) == 124L,
  sum(table(participant_key) == 1L) == 30L,
  sum(has_dot2_suffix) == 128L,
  sum(!has_dot2_suffix) == 150L,
  sum(vapply(
    participant_tables,
    function(x) length(x) == 1L && isTRUE(x),
    logical(1)
  )) == 4L,
  sum(vapply(
    participant_tables,
    function(x) length(x) == 1L && !isTRUE(x),
    logical(1)
  )) == 26L
)

bmi_label <- c("Lean", "Obese", "Overweight")[bmi_code + 1L]
repeat_availability <- ifelse(
  samples_per_key == 2L,
  "Two source samples",
  "One source sample"
)
source_suffix <- ifelse(has_dot2_suffix, "Dot2 suffix", "No suffix")

otutab <- source_counts
taxonomy <- data.frame(
  Kingdom = rep(NA_character_, nrow(otutab)),
  Phylum = rep(NA_character_, nrow(otutab)),
  Class = rep(NA_character_, nrow(otutab)),
  Order = rep(NA_character_, nrow(otutab)),
  Family = rep(NA_character_, nrow(otutab)),
  Genus = rownames(otutab),
  Species = rep(NA_character_, nrow(otutab)),
  TaxonomyStatus = rep(
    paste(
      "Source genus category; higher and lower ranks are not",
      "distributed with the packaged table"
    ),
    nrow(otutab)
  ),
  row.names = rownames(otutab),
  check.names = FALSE
)
metadata <- data.frame(
  BMIClass = bmi_label,
  BMIClassCode = bmi_code,
  ParticipantKey = participant_key,
  SourceSuffix = source_suffix,
  SamplesPerParticipantKey = as.integer(samples_per_key),
  RepeatAvailability = repeat_availability,
  row.names = sample_ids,
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

paths <- c(
  otutab = file.path(output_dir, "otutab.tsv"),
  taxonomy = file.path(output_dir, "taxonomy.tsv"),
  metadata = file.path(output_dir, "metadata.tsv")
)
write_tsv(otutab, "FeatureID", paths[["otutab"]])
write_tsv(taxonomy, "FeatureID", paths[["taxonomy"]])
write_tsv(metadata, "SampleID", paths[["metadata"]])

file_summary <- lapply(paths, function(path) {
  list(
    path = file.path(
      "data", "small", "community-typing-dmm", basename(path)
    ),
    bytes = unname(file.info(path)$size),
    sha256 = digest::digest(path, algo = "sha256", file = TRUE)
  )
})

summary <- list(
  dataset = "Holmes Twins genus-count table",
  source = "DirichletMultinomial 1.46.0 packaged extdata",
  source_files = list(
    list(
      path = "DirichletMultinomial/extdata/Twins.csv",
      bytes = unname(file.info(twins_path)$size),
      sha256 = unname(observed_source_sha256[["Twins.csv"]])
    ),
    list(
      path = "DirichletMultinomial/extdata/TwinStudy.t",
      bytes = unname(file.info(phenotype_path)$size),
      sha256 = unname(observed_source_sha256[["TwinStudy.t"]])
    )
  ),
  package = list(
    name = "DirichletMultinomial",
    version = observed_version,
    license = utils::packageDescription("DirichletMultinomial")$License,
    description_sha256 = digest::digest(
      description_path, algo = "sha256", file = TRUE
    ),
    doi = "10.18129/B9.bioc.DirichletMultinomial"
  ),
  primary_study = list(
    citation = paste(
      "Turnbaugh PJ, Hamady M, Yatsunenko T, et al.",
      "Nature. 2009;457:480-484."
    ),
    doi = "10.1038/nature07540"
  ),
  dmm_study = list(
    citation = paste(
      "Holmes I, Harris K, Quince C.",
      "PLOS ONE. 2012;7(2):e30126."
    ),
    doi = "10.1371/journal.pone.0030126"
  ),
  dimensions = list(
    features = nrow(otutab),
    samples = ncol(otutab),
    reads = unname(sum(otutab)),
    minimum_library_size = unname(min(colSums(otutab))),
    median_library_size = unname(stats::median(colSums(otutab))),
    maximum_library_size = unname(max(colSums(otutab))),
    zero_fraction = unname(mean(otutab == 0))
  ),
  metadata = list(
    bmi_class_counts = as.list(table(bmi_label)),
    participant_keys = length(unique(participant_key)),
    keys_with_two_source_samples = sum(table(participant_key) == 2L),
    keys_with_one_source_sample = sum(table(participant_key) == 1L),
    no_suffix_only_keys = sum(vapply(
      participant_tables,
      function(x) length(x) == 1L && !isTRUE(x),
      logical(1)
    )),
    dot2_suffix_only_keys = sum(vapply(
      participant_tables,
      function(x) length(x) == 1L && isTRUE(x),
      logical(1)
    )),
    suffix_boundary = paste(
      "ParticipantKey and SourceSuffix are deterministic string-derived",
      "audit fields. The package does not distribute a clinical visit table,",
      "so suffix order is not independently validated as visit timing."
    )
  ),
  taxonomy_boundary = paste(
    "The source table contains 129 named genus labels plus its original",
    "misspelled 'Uknown' category. Source labels are preserved in taxonomy.tsv;",
    "figures may display that category as 'Unclassified' without altering counts."
  ),
  count_contract = paste(
    "otutab.tsv contains raw non-negative integer counts with features in rows",
    "and samples in columns; no relative-abundance conversion, pseudocount,",
    "filtering, or rarefaction is applied during export."
  ),
  files = file_summary
)

summary_path <- file.path(output_dir, "source-summary.json")
jsonlite::write_json(
  summary,
  summary_path,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

message(sprintf(
  paste(
    "Prepared Holmes Twins DMM bundle: %s genus categories x %s samples;",
    "%s reads across %s suffix-derived participant keys."
  ),
  nrow(otutab),
  ncol(otutab),
  format(sum(otutab), big.mark = ","),
  length(unique(participant_key))
))
