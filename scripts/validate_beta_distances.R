#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "Asia/Shanghai")
set.seed(20260722)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c(
  "project-root", "input-dir", "output-dir", "figure-dir",
  "gemelli-python"
)
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0L) {
  stop(
    "Missing required arguments: ",
    paste(paste0("--", missing_args), collapse = ", "),
    call. = FALSE
  )
}

as_absolute <- function(path) {
  if (grepl("^/", path)) path else file.path(getwd(), path)
}
project_root <- normalizePath(as_absolute(args[["project-root"]]), mustWork = TRUE)
input_dir <- normalizePath(as_absolute(args[["input-dir"]]), mustWork = TRUE)
output_dir <- as_absolute(args[["output-dir"]])
figure_dir <- as_absolute(args[["figure-dir"]])
gemelli_candidate <- as_absolute(args[["gemelli-python"]])
gemelli_python <- file.path(
  normalizePath(dirname(gemelli_candidate), mustWork = TRUE),
  basename(gemelli_candidate)
)
if (!file.exists(gemelli_python)) {
  stop("gemelli Python does not exist: ", gemelli_python, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "ape", "digest", "dplyr", "ggplot2", "jsonlite", "phyloseq",
  "ragg", "readr", "scales", "svglite", "systemfonts", "tibble",
  "tidyr", "vegan"
)
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

source(file.path(project_root, "R", "theme_pub.R"), local = FALSE)

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

sanitize_text <- function(x) {
  x <- gsub(project_root, "<PROJECT_ROOT>", x, fixed = TRUE)
  x <- gsub(input_dir, "<INPUT_DIR>", x, fixed = TRUE)
  x <- gsub(output_dir, "<OUTPUT_DIR>", x, fixed = TRUE)
  x <- gsub(figure_dir, "<FIGURE_DIR>", x, fixed = TRUE)
  x <- gsub(dirname(gemelli_python), "<GEMELLI_BIN>", x, fixed = TRUE)
  x <- gsub(path.expand("~"), "<HOME>", x, fixed = TRUE)
  x
}

write_tsv <- function(x, filename) {
  readr::write_tsv(x, file.path(output_dir, filename), na = "")
}

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    na = character()
  )
  ids <- as.character(x[[1L]])
  if (anyDuplicated(ids)) {
    stop("Duplicated identifiers in ", basename(path), call. = FALSE)
  }
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

check_rows <- list()
add_check <- function(id, category, observed, expected, passed) {
  check_rows[[length(check_rows) + 1L]] <<- data.frame(
    check_id = id,
    category = category,
    observed = paste(observed, collapse = ","),
    expected = paste(expected, collapse = ","),
    status = ifelse(isTRUE(passed), "PASS", "FAIL"),
    check.names = FALSE
  )
}

command_path <- function(command) {
  path <- Sys.which(command)
  if (!nzchar(path)) "" else normalizePath(path, mustWork = TRUE)
}

run_command <- function(command, args = character(), env = character()) {
  status <- 0L
  output <- tryCatch(
    system2(
      command,
      args = args,
      stdout = TRUE,
      stderr = TRUE,
      env = env
    ),
    error = function(e) {
      status <<- 1L
      conditionMessage(e)
    }
  )
  command_status <- attr(output, "status")
  if (!is.null(command_status)) status <- as.integer(command_status)
  list(status = status, output = output)
}

input_paths <- c(
  otutab = file.path(input_dir, "otutab.tsv"),
  taxonomy = file.path(input_dir, "taxonomy.tsv"),
  metadata = file.path(input_dir, "metadata.tsv"),
  tree = file.path(input_dir, "rooted-tree.nwk.gz"),
  source_summary = file.path(input_dir, "source-summary.json"),
  prepare_script = file.path(project_root, "scripts", "prepare_beta_distance_data.R"),
  gemelli_runner = file.path(project_root, "scripts", "run_gemelli_rpca.py"),
  gemelli_env = file.path(project_root, "env", "gemelli.yml"),
  gemelli_requirements = file.path(
    project_root, "env", "gemelli-requirements.txt"
  )
)
expected_sha256 <- c(
  otutab = "727073a3f26aff6730ad54b55c11e1d52bd5d3d34f99dbf75d125d2e84d9722d",
  taxonomy = "1189e6e2a02474a665283f6b9962ec25716b6c7faf66011a6506e87a69ddf0ff",
  metadata = "c708b3b89edbb26a4aab26b5e575cdefef1eae28b9b7f93cea5f8f2edfc94221",
  tree = "535f512b3018da462a416b8f04ef36369bf1c1e383510494a790508361fc8c15",
  source_summary = "fac90c60d56ad18361145b46854cd7090b3fb1436e78a14017211ec13e41c94a",
  prepare_script = "10f92d70866995e1d5f1050b13bd3303c8957b27f150c177b6aef935d4ca02cf",
  gemelli_runner = "451e8e67f97429ddaca63126725517251c4bdeb7be75ddade9699abbd26f98cc",
  gemelli_env = "62b8895c264f8960573641729880b0c9094e72c9a195e785156864fca26a0ee2",
  gemelli_requirements = "1c1b612c8aa8c1f140cc9a176c67561af78869d9597fa9bb34d774f5ba7cfb2c"
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Article 21 input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

observed_sha256 <- vapply(input_paths, sha256_file, character(1))
input_audit <- data.frame(
  asset = names(input_paths),
  path = vapply(
    input_paths,
    function(path) sanitize_text(normalizePath(path, mustWork = TRUE)),
    character(1)
  ),
  bytes = unname(file.info(input_paths)$size),
  observed_sha256 = unname(observed_sha256),
  expected_sha256 = unname(expected_sha256[names(input_paths)]),
  status = ifelse(
    observed_sha256 == expected_sha256[names(input_paths)],
    "PASS",
    "FAIL"
  ),
  check.names = FALSE
)
for (i in seq_len(nrow(input_audit))) {
  add_check(
    paste0("sha256-", input_audit$asset[[i]]),
    "input",
    input_audit$observed_sha256[[i]],
    input_audit$expected_sha256[[i]],
    input_audit$status[[i]] == "PASS"
  )
}

expected_versions <- c(
  ape = "5.8",
  digest = "0.6.36",
  dplyr = "1.1.4",
  ggplot2 = "3.5.2",
  jsonlite = "1.8.8",
  phyloseq = "1.48.0",
  ragg = "1.3.2",
  readr = "2.1.5",
  scales = "1.3.0",
  svglite = "2.1.3",
  systemfonts = "1.1.0",
  tibble = "3.2.1",
  tidyr = "1.3.1",
  vegan = "2.6-6.1"
)
observed_versions <- vapply(
  names(expected_versions),
  function(package) packageDescription(package, fields = "Version"),
  character(1)
)
package_audit <- data.frame(
  package = names(expected_versions),
  expected_version = unname(expected_versions),
  observed_version = unname(observed_versions),
  status = ifelse(observed_versions == expected_versions, "PASS", "FAIL"),
  check.names = FALSE
)
r_version <- paste(R.version$major, R.version$minor, sep = ".")
add_check("r-version", "environment", r_version, "4.4.1", r_version == "4.4.1")
for (i in seq_len(nrow(package_audit))) {
  add_check(
    paste0("package-", package_audit$package[[i]]),
    "environment",
    package_audit$observed_version[[i]],
    package_audit$expected_version[[i]],
    package_audit$status[[i]] == "PASS"
  )
}

command_names <- c("pdfinfo", "pdffonts", "identify")
command_audit <- data.frame(
  command = c(command_names, "gemelli-python"),
  path = c(
    vapply(command_names, command_path, character(1)),
    gemelli_python
  ),
  stringsAsFactors = FALSE
)
command_audit$status <- ifelse(
  nzchar(command_audit$path) & file.exists(command_audit$path),
  "PASS",
  "FAIL"
)
for (i in seq_len(nrow(command_audit))) {
  add_check(
    paste0("command-", command_audit$command[[i]]),
    "environment",
    sanitize_text(command_audit$path[[i]]),
    "available",
    command_audit$status[[i]] == "PASS"
  )
}

font_family <- "DejaVu Sans"
font_matches <- systemfonts::match_fonts(font_family)[1L, , drop = FALSE]
font_audit <- data.frame(
  requested_family = font_family,
  resolved_path = sanitize_text(font_matches$path),
  exists = file.exists(font_matches$path),
  status = ifelse(file.exists(font_matches$path), "PASS", "FAIL"),
  stringsAsFactors = FALSE
)
add_check(
  "font-dejavu-sans",
  "environment",
  font_audit$resolved_path,
  "existing font file",
  font_audit$status == "PASS"
)

otu_frame <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])
source_summary <- jsonlite::read_json(input_paths[["source_summary"]])
tree <- ape::read.tree(gzfile(input_paths[["tree"]]))

counts <- suppressWarnings(data.matrix(otu_frame))
storage.mode(counts) <- "numeric"
feature_ids <- rownames(counts)
sample_ids <- colnames(counts)

add_check("feature-count", "data", nrow(counts), 18988, nrow(counts) == 18988L)
add_check("sample-count", "data", ncol(counts), 26, ncol(counts) == 26L)
add_check("read-count", "data", sum(counts), 28216678, sum(counts) == 28216678)
add_check(
  "taxonomy-feature-count", "data", nrow(taxonomy), 18988,
  nrow(taxonomy) == 18988L
)
add_check(
  "metadata-sample-count", "data", nrow(metadata), 26,
  nrow(metadata) == 26L
)
add_check(
  "taxonomy-ranks", "data", paste(colnames(taxonomy), collapse = ","),
  paste(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), collapse = ","),
  identical(
    colnames(taxonomy),
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  )
)
add_check("feature-ids-unique", "data", anyDuplicated(feature_ids), 0, !anyDuplicated(feature_ids))
add_check("sample-ids-unique", "data", anyDuplicated(sample_ids), 0, !anyDuplicated(sample_ids))
add_check(
  "taxonomy-feature-ids", "data",
  length(intersect(feature_ids, rownames(taxonomy))), length(feature_ids),
  setequal(feature_ids, rownames(taxonomy))
)
add_check(
  "metadata-sample-ids", "data",
  length(intersect(sample_ids, rownames(metadata))), length(sample_ids),
  setequal(sample_ids, rownames(metadata))
)
add_check(
  "counts-finite", "data", sum(is.finite(counts)), length(counts),
  all(is.finite(counts))
)
add_check("counts-nonnegative", "data", min(counts), ">=0", all(counts >= 0))
add_check(
  "counts-integer", "data", max(abs(counts - round(counts))), 0,
  all(abs(counts - round(counts)) < .Machine$double.eps^0.5)
)
add_check(
  "nonempty-features", "data", sum(rowSums(counts) > 0), nrow(counts),
  all(rowSums(counts) > 0)
)
add_check(
  "nonempty-samples", "data", sum(colSums(counts) > 0), ncol(counts),
  all(colSums(counts) > 0)
)
add_check(
  "habitat-classes", "data", length(unique(metadata$HabitatClass)), 5,
  length(unique(metadata$HabitatClass)) == 5L && !anyNA(metadata$HabitatClass)
)
add_check(
  "source-study-doi", "provenance", source_summary$primary_study$doi,
  "10.1073/pnas.1000080107",
  identical(source_summary$primary_study$doi, "10.1073/pnas.1000080107")
)

taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]

tree_tip_match <- setequal(tree$tip.label, feature_ids)
tree_audit <- data.frame(
  metric = c(
    "rooted", "tips", "unique_tips", "internal_nodes", "edges",
    "branch_lengths_present", "finite_branch_lengths",
    "negative_branch_lengths", "zero_branch_lengths", "minimum_branch_length",
    "maximum_branch_length", "tip_set_matches_otutab"
  ),
  observed = c(
    ape::is.rooted(tree),
    length(tree$tip.label),
    length(unique(tree$tip.label)),
    tree$Nnode,
    nrow(tree$edge),
    !is.null(tree$edge.length),
    !is.null(tree$edge.length) && all(is.finite(tree$edge.length)),
    sum(tree$edge.length < 0),
    sum(tree$edge.length == 0),
    min(tree$edge.length),
    max(tree$edge.length),
    tree_tip_match
  ),
  expected = c(
    TRUE, 18988, 18988, 18987, 37974, TRUE, TRUE, 0, 0,
    ">0", ">0", TRUE
  ),
  stringsAsFactors = FALSE
)
tree_audit$status <- "PASS"
add_check("tree-rooted", "tree", ape::is.rooted(tree), TRUE, ape::is.rooted(tree))
add_check("tree-tip-count", "tree", length(tree$tip.label), 18988, length(tree$tip.label) == 18988L)
add_check("tree-tip-unique", "tree", anyDuplicated(tree$tip.label), 0, !anyDuplicated(tree$tip.label))
add_check("tree-tip-match", "tree", tree_tip_match, TRUE, tree_tip_match)
add_check(
  "tree-branch-lengths", "tree", length(tree$edge.length), nrow(tree$edge),
  !is.null(tree$edge.length) && length(tree$edge.length) == nrow(tree$edge)
)
add_check(
  "tree-branch-finite", "tree", sum(is.finite(tree$edge.length)), nrow(tree$edge),
  all(is.finite(tree$edge.length))
)
add_check(
  "tree-branch-nonnegative", "tree", min(tree$edge.length), ">=0",
  all(tree$edge.length >= 0)
)

# Run the maintained gemelli implementation twice and require byte-identical
# canonical outputs. The runner fixes ARPACK initialization, PC signs and the
# precision written to disk.
runner_path <- input_paths[["gemelli_runner"]]
replay_dir <- tempfile("article21-gemelli-replay-")
dir.create(replay_dir, recursive = TRUE, showWarnings = FALSE)
gemelli_env_vars <- c(
  "PYTHONHASHSEED=20260722",
  "OMP_NUM_THREADS=1",
  "OPENBLAS_NUM_THREADS=1",
  "MKL_NUM_THREADS=1"
)
gemelli_args <- function(target_dir) {
  c(
    shQuote(runner_path),
    "--input", shQuote(input_paths[["otutab"]]),
    "--output-dir", shQuote(target_dir),
    "--n-components", "3",
    "--min-sample-count", "500",
    "--min-feature-count", "10",
    "--min-feature-frequency", "10",
    "--max-iterations", "5"
  )
}
gemelli_primary_run <- run_command(
  gemelli_python,
  gemelli_args(output_dir),
  env = gemelli_env_vars
)
add_check(
  "gemelli-primary-run", "execution", gemelli_primary_run$status, 0,
  gemelli_primary_run$status == 0L
)
if (gemelli_primary_run$status != 0L) {
  stop(paste(gemelli_primary_run$output, collapse = "\n"), call. = FALSE)
}
gemelli_replay_run <- run_command(
  gemelli_python,
  gemelli_args(replay_dir),
  env = gemelli_env_vars
)
add_check(
  "gemelli-replay-run", "execution", gemelli_replay_run$status, 0,
  gemelli_replay_run$status == 0L
)
if (gemelli_replay_run$status != 0L) {
  stop(paste(gemelli_replay_run$output, collapse = "\n"), call. = FALSE)
}

gemelli_files <- c(
  "robust-aitchison-distance.tsv",
  "robust-aitchison-sample-scores.tsv",
  "robust-aitchison-feature-loadings.tsv",
  "robust-aitchison-eigenvalues.tsv",
  "robust-aitchison-ordination.txt"
)
primary_hashes <- vapply(
  file.path(output_dir, gemelli_files), sha256_file, character(1)
)
replay_hashes <- vapply(
  file.path(replay_dir, gemelli_files), sha256_file, character(1)
)
gemelli_determinism_audit <- data.frame(
  file = gemelli_files,
  primary_sha256 = unname(primary_hashes),
  replay_sha256 = unname(replay_hashes),
  status = ifelse(primary_hashes == replay_hashes, "PASS", "FAIL"),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(gemelli_determinism_audit))) {
  add_check(
    paste0("gemelli-deterministic-", tools::file_path_sans_ext(gemelli_files[[i]])),
    "determinism",
    gemelli_determinism_audit$primary_sha256[[i]],
    gemelli_determinism_audit$replay_sha256[[i]],
    gemelli_determinism_audit$status[[i]] == "PASS"
  )
}

gemelli_summary <- jsonlite::read_json(
  file.path(output_dir, "gemelli-summary.json")
)
expected_python_versions <- c(
  python = "3.12.13",
  gemelli = "0.0.13",
  numpy = "2.4.2",
  pandas = "2.3.3",
  scipy = "1.17.1",
  scikit_learn = "1.7.1",
  scikit_bio = "0.7.2",
  biom_format = "2.1.17"
)
observed_python_versions <- c(
  python = gemelli_summary$versions$python,
  gemelli = gemelli_summary$versions$gemelli,
  numpy = gemelli_summary$versions$numpy,
  pandas = gemelli_summary$versions$pandas,
  scipy = gemelli_summary$versions$scipy,
  scikit_learn = gemelli_summary$versions$scikit_learn,
  scikit_bio = gemelli_summary$versions$scikit_bio,
  biom_format = gemelli_summary$versions$biom_format
)
python_package_audit <- data.frame(
  package = names(expected_python_versions),
  expected_version = unname(expected_python_versions),
  observed_version = unname(observed_python_versions),
  status = ifelse(
    observed_python_versions == expected_python_versions,
    "PASS",
    "FAIL"
  ),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(python_package_audit))) {
  add_check(
    paste0("python-package-", python_package_audit$package[[i]]),
    "environment",
    python_package_audit$observed_version[[i]],
    python_package_audit$expected_version[[i]],
    python_package_audit$status[[i]] == "PASS"
  )
}
pip_check <- run_command(
  gemelli_python,
  c("-m", "pip", "check"),
  env = gemelli_env_vars
)
add_check(
  "gemelli-pip-check", "environment",
  paste(pip_check$output, collapse = " "), "No broken requirements found.",
  pip_check$status == 0L && any(grepl("No broken requirements found", pip_check$output, fixed = TRUE))
)
add_check(
  "gemelli-filtered-features", "composition",
  gemelli_summary$filter$features_retained, 11063,
  gemelli_summary$filter$features_retained == 11063L
)
add_check(
  "gemelli-filtered-reads", "composition",
  gemelli_summary$filter$reads_retained, 28148221,
  gemelli_summary$filter$reads_retained == 28148221
)

community <- t(counts)
library_size <- rowSums(community)
rarefaction_depth <- 50000L
set.seed(20260722)
community_rarefied <- vegan::rrarefy(
  community,
  sample = rarefaction_depth
)
rarefaction_audit <- data.frame(
  SampleID = rownames(community_rarefied),
  HabitatClass = metadata[rownames(community_rarefied), "HabitatClass"],
  OriginalReads = as.numeric(library_size[rownames(community_rarefied)]),
  RarefiedReads = rowSums(community_rarefied),
  FeaturesBefore = rowSums(community > 0),
  FeaturesAfter = rowSums(community_rarefied > 0),
  ReadRetention = rowSums(community_rarefied) /
    as.numeric(library_size[rownames(community_rarefied)]),
  stringsAsFactors = FALSE
)
add_check(
  "rarefaction-minimum-depth-gate", "standardization",
  min(library_size), paste0(">=", rarefaction_depth),
  min(library_size) >= rarefaction_depth
)
add_check(
  "rarefaction-sample-count", "standardization",
  nrow(community_rarefied), 26, nrow(community_rarefied) == 26L
)
add_check(
  "rarefaction-equal-depth", "standardization",
  paste(range(rowSums(community_rarefied)), collapse = ","),
  "50000,50000",
  all(rowSums(community_rarefied) == rarefaction_depth)
)
add_check(
  "rarefaction-total-reads", "standardization",
  sum(community_rarefied), 1300000,
  sum(community_rarefied) == 1300000
)
add_check(
  "rarefaction-nonzero-features", "standardization",
  sum(colSums(community_rarefied) > 0), 12443,
  sum(colSums(community_rarefied) > 0) == 12443L
)

bray <- vegan::vegdist(community_rarefied, method = "bray")
jaccard <- vegan::vegdist(
  community_rarefied,
  method = "jaccard",
  binary = TRUE
)
bray_relative <- vegan::vegdist(
  vegan::decostand(community, method = "total"),
  method = "bray"
)
bray_relative_spearman <- stats::cor(
  as.vector(bray),
  as.vector(bray_relative),
  method = "spearman"
)
add_check(
  "bray-relative-sensitivity", "standardization",
  round(bray_relative_spearman, 8), ">=0.99",
  bray_relative_spearman >= 0.99
)

ps_rarefied <- phyloseq::phyloseq(
  phyloseq::otu_table(t(community_rarefied), taxa_are_rows = TRUE),
  phyloseq::phy_tree(tree)
)
ps_rarefied <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps_rarefied) > 0,
  ps_rarefied
)
unweighted_unifrac <- phyloseq::UniFrac(
  ps_rarefied,
  weighted = FALSE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)
weighted_unifrac <- phyloseq::UniFrac(
  ps_rarefied,
  weighted = TRUE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)
add_check(
  "unifrac-pruned-tree-tips", "tree",
  phyloseq::ntaxa(ps_rarefied), 12443,
  phyloseq::ntaxa(ps_rarefied) == 12443L
)

feature_total <- rowSums(counts)
feature_prevalence <- rowMeans(counts > 0)
keep_composition <- feature_total > 10 & feature_prevalence > 0.10
feature_filter_audit <- data.frame(
  FeatureID = feature_ids,
  TotalCount = as.numeric(feature_total),
  Prevalence = as.numeric(feature_prevalence),
  PassTotalCountGT10 = feature_total > 10,
  PassPrevalenceGT10Percent = feature_prevalence > 0.10,
  Retained = keep_composition,
  stringsAsFactors = FALSE
)
add_check(
  "composition-feature-count", "composition",
  sum(keep_composition), 11063, sum(keep_composition) == 11063L
)
add_check(
  "composition-read-count", "composition",
  sum(counts[keep_composition, , drop = FALSE]), 28148221,
  sum(counts[keep_composition, , drop = FALSE]) == 28148221
)

aitchison_distance <- function(pseudocount) {
  log_counts <- log(
    t(counts[keep_composition, , drop = FALSE]) + pseudocount
  )
  clr <- log_counts - rowMeans(log_counts)
  stats::dist(clr, method = "euclidean")
}
aitchison_pc01 <- aitchison_distance(0.1)
aitchison_pc05 <- aitchison_distance(0.5)
aitchison_pc10 <- aitchison_distance(1.0)

robust_frame <- readr::read_tsv(
  file.path(output_dir, "robust-aitchison-distance.tsv"),
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)
robust_ids <- as.character(robust_frame[[1L]])
robust_matrix <- data.matrix(robust_frame[-1L])
rownames(robust_matrix) <- robust_ids
colnames(robust_matrix) <- colnames(robust_frame)[-1L]
add_check(
  "robust-sample-id-set", "composition",
  length(intersect(sample_ids, robust_ids)), length(sample_ids),
  setequal(sample_ids, robust_ids) && setequal(sample_ids, colnames(robust_matrix))
)
robust_aitchison <- stats::as.dist(
  robust_matrix[sample_ids, sample_ids, drop = FALSE]
)

align_distance <- function(distance, ids = sample_ids) {
  matrix <- as.matrix(distance)
  stats::as.dist(matrix[ids, ids, drop = FALSE])
}
distances <- list(
  "Bray-Curtis" = align_distance(bray),
  "Jaccard" = align_distance(jaccard),
  "Unweighted UniFrac" = align_distance(unweighted_unifrac),
  "Weighted UniFrac" = align_distance(weighted_unifrac),
  "Aitchison" = align_distance(aitchison_pc05),
  "Robust Aitchison" = align_distance(robust_aitchison)
)

validate_distance <- function(distance, metric) {
  matrix <- as.matrix(distance)
  checks <- c(
    dimensions = identical(dim(matrix), c(26L, 26L)),
    ids = identical(rownames(matrix), sample_ids) &&
      identical(colnames(matrix), sample_ids),
    finite = all(is.finite(matrix)),
    nonnegative = min(matrix) >= -1e-12,
    symmetric = max(abs(matrix - t(matrix))) <= 1e-10,
    zero_diagonal = max(abs(diag(matrix))) <= 1e-10
  )
  for (check_name in names(checks)) {
    add_check(
      paste0(
        "distance-",
        gsub("[^a-z0-9]+", "-", tolower(metric)),
        "-",
        check_name
      ),
      "distance",
      checks[[check_name]],
      TRUE,
      checks[[check_name]]
    )
  }
  invisible(checks)
}
invisible(mapply(validate_distance, distances, names(distances)))

pair_indices <- which(upper.tri(matrix(0, 26, 26)), arr.ind = TRUE)
pairwise_rows <- list()
distance_summary_rows <- list()
for (metric in names(distances)) {
  matrix <- as.matrix(distances[[metric]])
  values <- matrix[pair_indices]
  sample_a <- rownames(matrix)[pair_indices[, 1L]]
  sample_b <- colnames(matrix)[pair_indices[, 2L]]
  habitat_a <- metadata[sample_a, "HabitatClass"]
  habitat_b <- metadata[sample_b, "HabitatClass"]
  pair_type <- ifelse(habitat_a == habitat_b, "Within habitat", "Between habitats")
  pairwise_rows[[metric]] <- data.frame(
    Metric = metric,
    SampleA = sample_a,
    SampleB = sample_b,
    HabitatA = habitat_a,
    HabitatB = habitat_b,
    PairType = pair_type,
    Distance = values,
    stringsAsFactors = FALSE
  )
  distance_summary_rows[[metric]] <- data.frame(
    Metric = metric,
    SamplePairs = length(values),
    Minimum = min(values),
    Median = stats::median(values),
    Maximum = max(values),
    WithinMedian = stats::median(values[pair_type == "Within habitat"]),
    BetweenMedian = stats::median(values[pair_type == "Between habitats"]),
    BetweenWithinRatio = stats::median(values[pair_type == "Between habitats"]) /
      stats::median(values[pair_type == "Within habitat"]),
    stringsAsFactors = FALSE
  )
}
pairwise_distances <- dplyr::bind_rows(pairwise_rows)
distance_summary <- dplyr::bind_rows(distance_summary_rows)
add_check(
  "distance-metric-count", "distance", length(distances), 6,
  length(distances) == 6L
)
add_check(
  "distance-pair-count", "distance", nrow(pairwise_distances), 1950,
  nrow(pairwise_distances) == 1950L
)

distance_vectors <- do.call(cbind, lapply(distances, as.vector))
distance_cor_matrix <- stats::cor(distance_vectors, method = "spearman")
distance_correlation <- tidyr::expand_grid(
  MetricX = rownames(distance_cor_matrix),
  MetricY = colnames(distance_cor_matrix)
)
distance_correlation$Spearman <- mapply(
  function(x, y) distance_cor_matrix[x, y],
  distance_correlation$MetricX,
  distance_correlation$MetricY
)
jaccard_unweighted_spearman <- unname(
  distance_cor_matrix["Jaccard", "Unweighted UniFrac"]
)
aitchison_robust_spearman <- unname(
  distance_cor_matrix["Aitchison", "Robust Aitchison"]
)
add_check(
  "jaccard-unweighted-spearman", "comparison",
  round(jaccard_unweighted_spearman, 8), ">=0.95",
  jaccard_unweighted_spearman >= 0.95
)
add_check(
  "aitchison-robust-spearman", "comparison",
  round(aitchison_robust_spearman, 8), "0.79 to 0.81",
  aitchison_robust_spearman >= 0.79 && aitchison_robust_spearman <= 0.81
)

pseudocount_distances <- list(
  "0.1" = align_distance(aitchison_pc01),
  "0.5" = align_distance(aitchison_pc05),
  "1.0" = align_distance(aitchison_pc10)
)
aitchison_pseudocount_sensitivity <- data.frame(
  Pseudocount = as.numeric(names(pseudocount_distances)),
  SpearmanToAitchisonPC05 = vapply(
    pseudocount_distances,
    function(x) stats::cor(
      as.vector(x), as.vector(aitchison_pc05), method = "spearman"
    ),
    numeric(1)
  ),
  SpearmanToRobustAitchison = vapply(
    pseudocount_distances,
    function(x) stats::cor(
      as.vector(x), as.vector(robust_aitchison), method = "spearman"
    ),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)
add_check(
  "pseudocount-sensitivity-minimum", "composition",
  round(min(aitchison_pseudocount_sensitivity$SpearmanToAitchisonPC05), 8),
  ">=0.99",
  min(aitchison_pseudocount_sensitivity$SpearmanToAitchisonPC05) >= 0.99
)

pcoa_score_rows <- list()
pcoa_audit_rows <- list()
for (metric in names(distances)) {
  distance <- distances[[metric]]
  raw_pcoa <- suppressWarnings(
    stats::cmdscale(
      distance,
      k = attr(distance, "Size") - 2L,
      eig = TRUE,
      add = FALSE
    )
  )
  eigenvalues <- as.numeric(raw_pcoa$eig)
  negative_count <- sum(eigenvalues < -1e-10)
  negative_sum <- sum(abs(eigenvalues[eigenvalues < 0]))
  positive_sum <- sum(eigenvalues[eigenvalues > 0])
  negative_positive_ratio <- ifelse(
    positive_sum > 0,
    negative_sum / positive_sum,
    NA_real_
  )
  correction <- if (negative_count > 0L) "Lingoes" else "None"
  pcoa_fit <- ape::pcoa(
    distance,
    correction = if (negative_count > 0L) "lingoes" else "none"
  )
  if (negative_count > 0L) {
    coordinates <- pcoa_fit$vectors.cor[, 1:2, drop = FALSE]
    relative_eigen <- pcoa_fit$values$Rel_corr_eig[1:2]
    lingoes_constant <- abs(min(eigenvalues))
  } else {
    coordinates <- pcoa_fit$vectors[, 1:2, drop = FALSE]
    relative_eigen <- pcoa_fit$values$Relative_eig[1:2]
    lingoes_constant <- 0
  }
  coordinates <- coordinates[sample_ids, , drop = FALSE]
  pcoa_score_rows[[metric]] <- data.frame(
    Metric = metric,
    SampleID = rownames(coordinates),
    HabitatClass = metadata[rownames(coordinates), "HabitatClass"],
    Axis1 = coordinates[, 1L],
    Axis2 = coordinates[, 2L],
    Axis1Explained = relative_eigen[[1L]],
    Axis2Explained = relative_eigen[[2L]],
    Correction = correction,
    stringsAsFactors = FALSE
  )
  pcoa_audit_rows[[metric]] <- data.frame(
    Metric = metric,
    PositiveEigenvalues = sum(eigenvalues > 1e-10),
    NegativeEigenvalues = negative_count,
    NearZeroEigenvalues = sum(abs(eigenvalues) <= 1e-10),
    MinimumEigenvalue = min(eigenvalues),
    MaximumEigenvalue = max(eigenvalues),
    NegativePositiveRatio = negative_positive_ratio,
    Correction = correction,
    LingoesConstant = lingoes_constant,
    Axis1Explained = relative_eigen[[1L]],
    Axis2Explained = relative_eigen[[2L]],
    stringsAsFactors = FALSE
  )
}
pcoa_scores <- dplyr::bind_rows(pcoa_score_rows)
pcoa_eigen_audit <- dplyr::bind_rows(pcoa_audit_rows)
weighted_negative_count <- pcoa_eigen_audit$NegativeEigenvalues[
  pcoa_eigen_audit$Metric == "Weighted UniFrac"
]
add_check(
  "pcoa-score-count", "ordination", nrow(pcoa_scores), 156,
  nrow(pcoa_scores) == 156L
)
add_check(
  "pcoa-weighted-negative-eigenvalues", "ordination",
  weighted_negative_count, 3, weighted_negative_count == 3L
)
add_check(
  "pcoa-weighted-lingoes", "ordination",
  pcoa_eigen_audit$Correction[pcoa_eigen_audit$Metric == "Weighted UniFrac"],
  "Lingoes",
  pcoa_eigen_audit$Correction[
    pcoa_eigen_audit$Metric == "Weighted UniFrac"
  ] == "Lingoes"
)
add_check(
  "pcoa-other-no-correction", "ordination",
  sum(
    pcoa_eigen_audit$Metric != "Weighted UniFrac" &
      pcoa_eigen_audit$Correction == "None"
  ),
  5,
  all(
    pcoa_eigen_audit$Correction[
      pcoa_eigen_audit$Metric != "Weighted UniFrac"
    ] == "None"
  )
)

metric_levels <- names(distances)
habitat_levels <- c(
  "Human-associated", "Freshwater", "Marine/estuarine", "Soil", "Mock"
)
habitat_palette <- c(
  "Human-associated" = pal_pub[["blue"]],
  "Freshwater" = pal_pub[["sky"]],
  "Marine/estuarine" = pal_pub[["green"]],
  "Soil" = pal_pub[["orange"]],
  "Mock" = pal_pub[["purple"]]
)
habitat_shapes <- c(
  "Human-associated" = 21,
  "Freshwater" = 22,
  "Marine/estuarine" = 24,
  "Soil" = 23,
  "Mock" = 25
)

choice_grid <- tidyr::expand_grid(
  Metric = metric_levels,
  Signal = c("Abundance", "Occurrence", "Phylogeny", "Log-ratio")
)
choice_grid$Included <- with(
  choice_grid,
  (Metric == "Bray-Curtis" & Signal == "Abundance") |
    (Metric == "Jaccard" & Signal == "Occurrence") |
    (Metric == "Unweighted UniFrac" &
      Signal %in% c("Occurrence", "Phylogeny")) |
    (Metric == "Weighted UniFrac" &
      Signal %in% c("Abundance", "Phylogeny")) |
    (Metric %in% c("Aitchison", "Robust Aitchison") &
      Signal %in% c("Abundance", "Log-ratio"))
)
choice_grid$Label <- ifelse(choice_grid$Included, "YES", "—")
choice_grid$Metric <- factor(choice_grid$Metric, levels = rev(metric_levels))
choice_grid$Signal <- factor(
  choice_grid$Signal,
  levels = c("Abundance", "Occurrence", "Phylogeny", "Log-ratio")
)
choice_plot <- ggplot2::ggplot(
  choice_grid,
  ggplot2::aes(x = Signal, y = Metric, fill = Included)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 1.1) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label),
    colour = "#1A1A1A",
    fontface = "bold",
    size = 3.2
  ) +
  ggplot2::scale_fill_manual(
    values = c(`FALSE` = "#ECEFF3", `TRUE` = pal_pub[["sky"]]),
    guide = "none"
  ) +
  ggplot2::labs(
    title = "Distance metrics encode different biological contrasts",
    subtitle = "Choose the signal before inspecting group separation",
    x = "Information used by the metric",
    y = NULL,
    caption = paste0(
      "Occurrence and UniFrac branches use 50,000 reads per sample.\n",
      "Aitchison branches use filtered unrarefied counts; classic CLR uses a declared pseudocount."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(face = "bold")
  )

correlation_plot_data <- distance_correlation
correlation_plot_data$MetricX <- factor(
  correlation_plot_data$MetricX,
  levels = metric_levels
)
correlation_plot_data$MetricY <- factor(
  correlation_plot_data$MetricY,
  levels = rev(metric_levels)
)
correlation_plot <- ggplot2::ggplot(
  correlation_plot_data,
  ggplot2::aes(x = MetricX, y = MetricY, fill = Spearman)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", Spearman)),
    size = 2.8,
    colour = "#111827"
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c("#F7F7F7", "#B8DCEB", "#0072B2"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  ggplot2::labs(
    title = "Metric choice changes the rank order of sample pairs",
    subtitle = "Spearman correlation across all 325 pairwise distances",
    x = NULL,
    y = NULL,
    fill = "Spearman\nrho",
    caption = "High correlation in this dataset does not make two metric definitions interchangeable."
  ) +
  theme_pub(base_size = 8.5) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      size = 7.5
    ),
    axis.text.y = ggplot2::element_text(size = 7.5),
    legend.position = "right"
  )

pcoa_plot_data <- pcoa_scores
pcoa_plot_data$HabitatClass <- factor(
  pcoa_plot_data$HabitatClass,
  levels = habitat_levels
)
pcoa_labels <- setNames(
  vapply(
    metric_levels,
    function(metric) {
      row <- pcoa_eigen_audit[pcoa_eigen_audit$Metric == metric, ]
      correction_label <- ifelse(
        row$Correction == "Lingoes",
        "\nLingoes corrected",
        ""
      )
      sprintf(
        "%s\nPCo 1 %.1f%% · PCo 2 %.1f%%%s",
        metric,
        100 * row$Axis1Explained,
        100 * row$Axis2Explained,
        correction_label
      )
    },
    character(1)
  ),
  metric_levels
)
pcoa_plot_data$Metric <- factor(pcoa_plot_data$Metric, levels = metric_levels)
pcoa_plot <- ggplot2::ggplot(
  pcoa_plot_data,
  ggplot2::aes(
    x = Axis1,
    y = Axis2,
    fill = HabitatClass,
    shape = HabitatClass
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "#D1D5DB",
    linewidth = 0.3
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    colour = "#D1D5DB",
    linewidth = 0.3
  ) +
  ggplot2::geom_point(
    size = 2.2,
    stroke = 0.45,
    colour = "#1F2937",
    alpha = 0.92
  ) +
  ggplot2::facet_wrap(
    ~Metric,
    scales = "free",
    ncol = 3,
    labeller = ggplot2::as_labeller(pcoa_labels)
  ) +
  ggplot2::scale_fill_manual(values = habitat_palette, drop = FALSE) +
  ggplot2::scale_shape_manual(values = habitat_shapes, drop = FALSE) +
  ggplot2::labs(
    title = "The same samples occupy six different distance geometries",
    subtitle = "GlobalPatterns 16S communities; each point is one sample",
    x = "Principal coordinate 1",
    y = "Principal coordinate 2",
    fill = "Habitat",
    shape = "Habitat",
    caption = "Coordinates visualize distance structure; they are not tests of group separation."
  ) +
  theme_pub(base_size = 8.2) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = "top",
    legend.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(size = 7.1, lineheight = 0.95)
  )

base_pair_type <- pairwise_rows[["Robust Aitchison"]]$PairType
zero_plot_rows <- lapply(
  names(pseudocount_distances),
  function(pc) {
    data.frame(
      Pseudocount = pc,
      ClassicAitchison = as.vector(pseudocount_distances[[pc]]),
      RobustAitchison = as.vector(robust_aitchison),
      PairType = base_pair_type,
      stringsAsFactors = FALSE
    )
  }
)
zero_plot_data <- dplyr::bind_rows(zero_plot_rows)
zero_plot_data$Pseudocount <- factor(
  zero_plot_data$Pseudocount,
  levels = c("0.1", "0.5", "1.0"),
  labels = c("Pseudocount 0.1", "Pseudocount 0.5", "Pseudocount 1.0")
)
zero_annotation <- data.frame(
  Pseudocount = levels(zero_plot_data$Pseudocount),
  Label = sprintf(
    "rho with robust = %.3f",
    aitchison_pseudocount_sensitivity$SpearmanToRobustAitchison
  ),
  stringsAsFactors = FALSE
)
zero_annotation$Pseudocount <- factor(
  zero_annotation$Pseudocount,
  levels = levels(zero_plot_data$Pseudocount)
)
zero_plot <- ggplot2::ggplot(
  zero_plot_data,
  ggplot2::aes(
    x = ClassicAitchison,
    y = RobustAitchison,
    colour = PairType
  )
) +
  ggplot2::geom_point(size = 1.25, alpha = 0.62) +
  ggplot2::geom_smooth(
    ggplot2::aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.55,
    colour = "#374151",
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = zero_annotation,
    ggplot2::aes(x = -Inf, y = Inf, label = Label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.35,
    size = 2.7,
    fontface = "bold",
    colour = "#1F2937"
  ) +
  ggplot2::facet_wrap(~Pseudocount, nrow = 1, scales = "free_x") +
  ggplot2::scale_colour_manual(
    values = c(
      "Within habitat" = pal_pub[["blue"]],
      "Between habitats" = pal_pub[["vermillion"]]
    )
  ) +
  ggplot2::labs(
    title = "Zero handling changes Aitchison geometry",
    subtitle = "Classic pseudocount CLR versus observed-only Robust Aitchison",
    x = "Classic Aitchison distance",
    y = "Robust Aitchison distance",
    colour = "Sample pair",
    caption = paste(
      "Both branches use the same 11,063 filtered features.",
      "Robust Aitchison does not replace observed zeros with the displayed pseudocounts."
    )
  ) +
  theme_pub(base_size = 8.5) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 8)
  )

figure_specification <- data.frame(
  stem = c(
    "21-distance-choice-map",
    "21-distance-correlation",
    "21-pcoa-metric-comparison",
    "21-aitchison-zero-sensitivity"
  ),
  width_mm = c(183, 183, 183, 183),
  height_mm = c(92, 112, 132, 94),
  raster_ppi = rep(600L, 4L),
  stringsAsFactors = FALSE
)
plot_objects <- list(
  choice_plot,
  correlation_plot,
  pcoa_plot,
  zero_plot
)
for (i in seq_len(nrow(figure_specification))) {
  save_pub(
    plot = plot_objects[[i]],
    file_base = file.path(figure_dir, figure_specification$stem[[i]]),
    width = figure_specification$width_mm[[i]],
    height = figure_specification$height_mm[[i]],
    units = "mm",
    dpi = figure_specification$raster_ppi[[i]],
    write_svg = TRUE,
    write_tiff = TRUE,
    base_family = font_family
  )
}

parse_pdf_info <- function(path) {
  result <- run_command("pdfinfo", shQuote(path))
  if (result$status != 0L) stop("pdfinfo failed for ", path, call. = FALSE)
  page_line <- grep("^Pages:", result$output, value = TRUE)
  size_line <- grep("^Page size:", result$output, value = TRUE)
  size_match <- regexec(
    "Page size:[[:space:]]+([0-9.]+)[[:space:]]+x[[:space:]]+([0-9.]+)[[:space:]]+pts",
    size_line
  )
  size_parts <- regmatches(size_line, size_match)[[1L]]
  if (length(page_line) != 1L || length(size_parts) != 3L) {
    stop("Unexpected pdfinfo output for ", path, call. = FALSE)
  }
  data.frame(
    pages = as.integer(sub("^Pages:[[:space:]]+", "", page_line)),
    width_pt = as.numeric(size_parts[[2L]]),
    height_pt = as.numeric(size_parts[[3L]]),
    stringsAsFactors = FALSE
  )
}

parse_svg_info <- function(path) {
  svg_text <- paste(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  root_match <- regexec(
    "<svg[^>]*width='([0-9.]+)pt'[^>]*height='([0-9.]+)pt'",
    svg_text
  )
  root_parts <- regmatches(svg_text, root_match)[[1L]]
  if (length(root_parts) != 3L) {
    stop("Unexpected SVG root dimensions for ", path, call. = FALSE)
  }
  data.frame(
    width_pt = as.numeric(root_parts[[2L]]),
    height_pt = as.numeric(root_parts[[3L]]),
    text_nodes = lengths(
      regmatches(svg_text, gregexpr("<text\\b", svg_text, perl = TRUE))
    ),
    font_family_present = grepl(
      'font-family: "DejaVu Sans"',
      svg_text,
      fixed = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

identify_raster <- function(path) {
  result <- run_command(
    "identify",
    c("-format", shQuote("%m|%w|%h|%x|%y|%U|%C"), shQuote(path))
  )
  if (result$status != 0L || length(result$output) == 0L) {
    stop("identify failed for ", path, call. = FALSE)
  }
  fields <- strsplit(result$output[[1L]], "\\|", fixed = FALSE)[[1L]]
  if (length(fields) != 7L) {
    stop("Unexpected identify output for ", path, call. = FALSE)
  }
  data.frame(
    detected_format = fields[[1L]],
    pixel_width = as.integer(fields[[2L]]),
    pixel_height = as.integer(fields[[3L]]),
    density_x = as.numeric(fields[[4L]]),
    density_y = as.numeric(fields[[5L]]),
    density_unit = fields[[6L]],
    compression = fields[[7L]],
    stringsAsFactors = FALSE
  )
}

file_signature_ok <- function(path, extension) {
  raw <- readBin(path, what = "raw", n = 16L)
  if (extension == "pdf") {
    return(identical(rawToChar(raw[seq_len(4L)]), "%PDF"))
  }
  if (extension == "png") {
    return(identical(
      as.integer(raw[seq_len(8L)]),
      c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)
    ))
  }
  if (extension == "tiff") {
    signature <- as.integer(raw[seq_len(4L)])
    return(
      identical(signature, c(73L, 73L, 42L, 0L)) ||
        identical(signature, c(77L, 77L, 0L, 42L))
    )
  }
  if (extension == "svg") {
    prefix <- rawToChar(raw)
    return(startsWith(prefix, "<?xml") || startsWith(prefix, "<svg"))
  }
  FALSE
}

format_rows <- list()
for (spec_i in seq_len(nrow(figure_specification))) {
  spec <- figure_specification[spec_i, , drop = FALSE]
  for (extension in c("pdf", "svg", "png", "tiff")) {
    path <- file.path(figure_dir, paste0(spec$stem, ".", extension))
    exists <- file.exists(path)
    bytes <- if (exists) file.info(path)$size else 0
    add_check(
      paste0(spec$stem, "-", extension, "-exists"),
      "export", exists, TRUE, exists
    )
    add_check(
      paste0(spec$stem, "-", extension, "-nonempty"),
      "export", bytes, ">1000", bytes > 1000
    )
    signature_ok <- exists && file_signature_ok(path, extension)
    add_check(
      paste0(spec$stem, "-", extension, "-signature"),
      "export", signature_ok, TRUE, signature_ok
    )

    actual_width_mm <- NA_real_
    actual_height_mm <- NA_real_
    pixel_width <- NA_integer_
    pixel_height <- NA_integer_
    effective_ppi_x <- NA_real_
    effective_ppi_y <- NA_real_
    detected_format <- toupper(extension)
    compression <- NA_character_
    pages <- NA_integer_
    text_nodes <- NA_integer_
    font_status <- NA_character_

    if (extension == "pdf") {
      pdf_info <- parse_pdf_info(path)
      actual_width_mm <- pdf_info$width_pt * 25.4 / 72
      actual_height_mm <- pdf_info$height_pt * 25.4 / 72
      pages <- pdf_info$pages
      font_result <- run_command("pdffonts", shQuote(path))
      font_output <- paste(font_result$output, collapse = "\n")
      font_embedded <- font_result$status == 0L &&
        grepl("DejaVuSans", font_output, fixed = TRUE) &&
        grepl("yes[[:space:]]+yes[[:space:]]+yes", font_output)
      font_status <- ifelse(font_embedded, "embedded", "missing")
      add_check(
        paste0(spec$stem, "-pdf-single-page"),
        "vector", pages, 1, pages == 1L
      )
      add_check(
        paste0(spec$stem, "-pdf-font-embedded"),
        "vector", font_status, "embedded", font_embedded
      )
    } else if (extension == "svg") {
      svg_info <- parse_svg_info(path)
      actual_width_mm <- svg_info$width_pt * 25.4 / 72
      actual_height_mm <- svg_info$height_pt * 25.4 / 72
      text_nodes <- svg_info$text_nodes
      font_status <- ifelse(svg_info$font_family_present, "declared", "missing")
      add_check(
        paste0(spec$stem, "-svg-text-nodes"),
        "vector", text_nodes, ">0", text_nodes > 0L
      )
      add_check(
        paste0(spec$stem, "-svg-font-declared"),
        "vector", font_status, "declared", svg_info$font_family_present
      )
    } else {
      raster_info <- identify_raster(path)
      detected_format <- raster_info$detected_format
      pixel_width <- raster_info$pixel_width
      pixel_height <- raster_info$pixel_height
      compression <- raster_info$compression
      actual_width_mm <- spec$width_mm
      actual_height_mm <- spec$height_mm
      effective_ppi_x <- pixel_width / (spec$width_mm / 25.4)
      effective_ppi_y <- pixel_height / (spec$height_mm / 25.4)
      expected_width_px <- floor(spec$width_mm / 25.4 * spec$raster_ppi)
      expected_height_px <- floor(spec$height_mm / 25.4 * spec$raster_ppi)
      add_check(
        paste0(spec$stem, "-", extension, "-pixel-dimensions"),
        "raster", c(pixel_width, pixel_height),
        c(expected_width_px, expected_height_px),
        all(c(pixel_width, pixel_height) == c(expected_width_px, expected_height_px))
      )
      add_check(
        paste0(spec$stem, "-", extension, "-effective-ppi"),
        "raster", round(c(effective_ppi_x, effective_ppi_y), 3),
        paste0(spec$raster_ppi, " +/- 0.5"),
        all(abs(c(effective_ppi_x, effective_ppi_y) - spec$raster_ppi) <= 0.5)
      )
      add_check(
        paste0(spec$stem, "-", extension, "-format"),
        "raster", detected_format, toupper(extension),
        identical(detected_format, toupper(extension))
      )
      if (extension == "tiff") {
        add_check(
          paste0(spec$stem, "-tiff-lzw"),
          "raster", compression, "LZW", identical(compression, "LZW")
        )
      }
    }

    dimension_pass <-
      abs(actual_width_mm - spec$width_mm) <= 0.4 &&
      abs(actual_height_mm - spec$height_mm) <= 0.4
    add_check(
      paste0(spec$stem, "-", extension, "-physical-dimensions"),
      "export", round(c(actual_width_mm, actual_height_mm), 3),
      c(spec$width_mm, spec$height_mm), dimension_pass
    )

    format_rows[[length(format_rows) + 1L]] <- data.frame(
      figure = spec$stem,
      extension = extension,
      path = file.path("figures", paste0(spec$stem, ".", extension)),
      detected_format = detected_format,
      bytes = bytes,
      target_width_mm = spec$width_mm,
      target_height_mm = spec$height_mm,
      actual_width_mm = actual_width_mm,
      actual_height_mm = actual_height_mm,
      pixel_width = pixel_width,
      pixel_height = pixel_height,
      effective_ppi_x = effective_ppi_x,
      effective_ppi_y = effective_ppi_y,
      compression = compression,
      pages = pages,
      text_nodes = text_nodes,
      font_status = font_status,
      dimension_status = ifelse(dimension_pass, "PASS", "FAIL"),
      stringsAsFactors = FALSE
    )
  }
}
format_audit <- dplyr::bind_rows(format_rows)
add_check(
  "primary-figure-count", "export",
  length(unique(format_audit$figure)), 4,
  length(unique(format_audit$figure)) == 4L
)
add_check(
  "primary-format-count", "export", nrow(format_audit), 16,
  nrow(format_audit) == 16L
)
add_check(
  "vector-file-count", "export",
  sum(format_audit$extension %in% c("pdf", "svg")), 8,
  sum(format_audit$extension %in% c("pdf", "svg")) == 8L
)
add_check(
  "raster-file-count", "export",
  sum(format_audit$extension %in% c("png", "tiff")), 8,
  sum(format_audit$extension %in% c("png", "tiff")) == 8L
)
add_check(
  "all-physical-dimensions-pass", "export",
  sum(format_audit$dimension_status == "PASS"), 16,
  all(format_audit$dimension_status == "PASS")
)

validation_checks <- dplyr::bind_rows(check_rows)
checks_total <- nrow(validation_checks)
checks_passed <- sum(validation_checks$status == "PASS")
checks_failed <- checks_total - checks_passed

write_tsv(input_audit, "input-audit.tsv")
write_tsv(package_audit, "package-audit.tsv")
write_tsv(python_package_audit, "python-package-audit.tsv")
write_tsv(command_audit, "command-audit.tsv")
write_tsv(font_audit, "font-audit.tsv")
write_tsv(tree_audit, "tree-audit.tsv")
write_tsv(rarefaction_audit, "rarefaction-audit.tsv")
write_tsv(feature_filter_audit, "feature-filter-audit.tsv")
write_tsv(distance_summary, "distance-summary.tsv")
write_tsv(distance_correlation, "distance-correlation.tsv")
write_tsv(pcoa_eigen_audit, "pcoa-eigen-audit.tsv")
write_tsv(pcoa_scores, "pcoa-scores.tsv")
write_tsv(pairwise_distances, "pairwise-distances.tsv")
write_tsv(
  aitchison_pseudocount_sensitivity,
  "aitchison-pseudocount-sensitivity.tsv"
)
write_tsv(gemelli_determinism_audit, "gemelli-determinism-audit.tsv")
write_tsv(figure_specification, "figure-specification.tsv")
write_tsv(format_audit, "format-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

summary_payload <- list(
  article = 21L,
  title = "Beta distance metric, tree and zero-handling audit",
  dataset = list(
    source = "phyloseq 1.48.0 GlobalPatterns",
    primary_study_doi = source_summary$primary_study$doi,
    source_features = source_summary$dimensions$source_features,
    zero_sum_features_removed = source_summary$dimensions$zero_sum_features_removed,
    features = nrow(counts),
    samples = ncol(counts),
    reads = sum(counts),
    min_depth = min(library_size),
    median_depth = unname(stats::median(library_size)),
    max_depth = max(library_size),
    zero_fraction = mean(counts == 0),
    habitat_classes = length(unique(metadata$HabitatClass))
  ),
  tree = list(
    rooted = ape::is.rooted(tree),
    tips = length(tree$tip.label),
    internal_nodes = tree$Nnode,
    branch_lengths_present = !is.null(tree$edge.length),
    negative_branch_lengths = sum(tree$edge.length < 0),
    tip_set_matches_counts = tree_tip_match
  ),
  standardization = list(
    seed = 20260722L,
    rarefaction_depth = rarefaction_depth,
    samples_retained = nrow(community_rarefied),
    rarefied_reads = sum(community_rarefied),
    nonzero_features_after_rarefaction = sum(colSums(community_rarefied) > 0),
    bray_relative_abundance_spearman = bray_relative_spearman
  ),
  composition = list(
    min_feature_count_exclusive = 10L,
    min_feature_prevalence_percent_exclusive = 10,
    features_retained = sum(keep_composition),
    reads_retained = sum(counts[keep_composition, , drop = FALSE]),
    classic_aitchison_default_pseudocount = 0.5,
    pseudocounts_audited = c(0.1, 0.5, 1.0)
  ),
  analysis = list(
    distance_metrics = length(distances),
    pairwise_distances_per_metric = choose(ncol(counts), 2L),
    pairwise_distance_rows = nrow(pairwise_distances)
  ),
  comparisons = list(
    jaccard_unweighted_spearman = jaccard_unweighted_spearman,
    aitchison_robust_spearman = aitchison_robust_spearman,
    pseudocount_minimum_spearman_to_default = min(
      aitchison_pseudocount_sensitivity$SpearmanToAitchisonPC05
    )
  ),
  robust_aitchison = list(
    method_family = "DEICODE rclr/RPCA",
    maintained_implementation = "gemelli 0.0.13",
    python = gemelli_summary$versions$python,
    features_retained = gemelli_summary$filter$features_retained,
    samples_retained = gemelli_summary$filter$samples_retained,
    reads_retained = gemelli_summary$filter$reads_retained,
    n_components = gemelli_summary$model$n_components,
    max_iterations = gemelli_summary$model$max_iterations,
    seed = gemelli_summary$seed,
    arpack_seed = gemelli_summary$model$arpack_seed,
    canonical_round_digits = gemelli_summary$model$canonical_round_digits,
    deterministic_outputs_matched = sum(
      gemelli_determinism_audit$status == "PASS"
    )
  ),
  pcoa = list(
    weighted_unifrac_negative_eigenvalues = weighted_negative_count,
    weighted_unifrac_negative_positive_ratio = pcoa_eigen_audit$NegativePositiveRatio[
      pcoa_eigen_audit$Metric == "Weighted UniFrac"
    ],
    weighted_unifrac_correction = pcoa_eigen_audit$Correction[
      pcoa_eigen_audit$Metric == "Weighted UniFrac"
    ],
    metrics_without_correction = sum(pcoa_eigen_audit$Correction == "None")
  ),
  graphics = list(
    font_family = font_family,
    primary_figures = length(unique(format_audit$figure)),
    vector_files = sum(format_audit$extension %in% c("pdf", "svg")),
    raster_files = sum(format_audit$extension %in% c("png", "tiff")),
    format_files = nrow(format_audit),
    raster_ppi = 600L,
    svg_files_with_text = sum(
      format_audit$extension == "svg" & format_audit$text_nodes > 0
    ),
    pdf_files_with_embedded_font = sum(
      format_audit$extension == "pdf" &
        format_audit$font_status == "embedded"
    ),
    lzw_tiff_files = sum(
      format_audit$extension == "tiff" &
        format_audit$compression == "LZW"
    ),
    dimension_checks_passed = sum(format_audit$dimension_status == "PASS")
  ),
  r_versions = as.list(observed_versions),
  python_versions = as.list(observed_python_versions),
  checks_total = checks_total,
  checks_passed = checks_passed,
  checks_failed = checks_failed
)
jsonlite::write_json(
  summary_payload,
  file.path(output_dir, "beta-distances-summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 12
)

session_lines <- sanitize_text(capture.output(utils::sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"))
log_lines <- c(
  "Article 21 beta-distances validation",
  paste0("R version: ", r_version),
  paste0(
    "Dataset: ", nrow(counts), " features x ", ncol(counts),
    " samples; reads=", sum(counts)
  ),
  paste0(
    "Tree: rooted=", ape::is.rooted(tree),
    "; tips=", length(tree$tip.label),
    "; exact tip match=", tree_tip_match
  ),
  paste0(
    "Rarefaction: seed=20260722; depth=", rarefaction_depth,
    "; reads=", sum(community_rarefied),
    "; nonzero features=", sum(colSums(community_rarefied) > 0)
  ),
  paste0(
    "Composition filter: features=", sum(keep_composition),
    "; reads=", sum(counts[keep_composition, , drop = FALSE])
  ),
  paste0(
    "Jaccard / unweighted UniFrac Spearman: ",
    round(jaccard_unweighted_spearman, 8)
  ),
  paste0(
    "Aitchison / Robust Aitchison Spearman: ",
    round(aitchison_robust_spearman, 8)
  ),
  paste0(
    "Weighted UniFrac PCoA: negative eigenvalues=", weighted_negative_count,
    "; correction=Lingoes"
  ),
  paste0(
    "gemelli primary: ",
    paste(gemelli_primary_run$output, collapse = " | ")
  ),
  paste0(
    "gemelli replay: ",
    paste(gemelli_replay_run$output, collapse = " | ")
  ),
  paste0(
    "gemelli deterministic outputs: ",
    sum(gemelli_determinism_audit$status == "PASS"),
    "/", nrow(gemelli_determinism_audit)
  ),
  paste0("Primary exports: ", nrow(format_audit), " files"),
  paste0("Checks passed: ", checks_passed, "/", checks_total),
  "",
  paste(
    validation_checks$check_id,
    validation_checks$status,
    validation_checks$observed,
    sep = "\t"
  )
)
writeLines(sanitize_text(log_lines), file.path(output_dir, "validation.log"))
unlink(replay_dir, recursive = TRUE, force = TRUE)

if (checks_failed > 0L) {
  failed <- validation_checks$check_id[validation_checks$status == "FAIL"]
  stop(
    "Article 21 validation failed: ",
    paste(failed, collapse = ", "),
    call. = FALSE
  )
}

message(
  "Article 21 validation passed: ",
  checks_passed,
  "/",
  checks_total,
  " checks."
)
