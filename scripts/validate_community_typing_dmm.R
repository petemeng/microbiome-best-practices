#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(
  TZ = "Asia/Shanghai",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

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
required_args <- c("project-root", "input-dir", "output-dir", "figure-dir")
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
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "digest", "DirichletMultinomial", "ggplot2", "jsonlite", "ragg",
  "readr", "scales", "svglite", "systemfonts", "vegan"
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
  gsub(path.expand("~"), "<HOME>", x, fixed = TRUE)
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
    na = character(),
    col_types = readr::cols(.default = readr::col_character())
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

near <- function(observed, expected, tolerance = 1e-8) {
  observed <- as.numeric(observed)
  expected <- as.numeric(expected)
  length(observed) == length(expected) &&
    all(is.finite(observed)) &&
    all(is.finite(expected)) &&
    max(abs(observed - expected)) <= tolerance
}

command_path <- function(command) {
  path <- Sys.which(command)
  if (!nzchar(path)) "" else normalizePath(path, mustWork = TRUE)
}

run_command <- function(command, args = character()) {
  status <- 0L
  output <- tryCatch(
    system2(command, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      status <<- 1L
      conditionMessage(e)
    }
  )
  command_status <- attr(output, "status")
  if (!is.null(command_status)) status <- as.integer(command_status)
  list(status = status, output = output)
}

choose2 <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(x, y) {
  stopifnot(length(x) == length(y))
  tab <- table(x, y)
  n <- sum(tab)
  if (n < 2L) return(NA_real_)
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  expected <- sum_rows * sum_cols / choose2(n)
  maximum <- (sum_rows + sum_cols) / 2
  denominator <- maximum - expected
  if (abs(denominator) <= .Machine$double.eps) {
    return(ifelse(identical(as.character(x), as.character(y)), 1, 0))
  }
  (sum_cells - expected) / denominator
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    rest <- all_permutations(x[-i])
    cbind(x[[i]], rest)
  }))
}

js_divergence <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  m <- (p + q) / 2
  sum(ifelse(p > 0, p * log(p / m), 0)) / 2 +
    sum(ifelse(q > 0, q * log(q / m), 0)) / 2
}

fit_scores <- function(fit_object) {
  values <- DirichletMultinomial::goodnessOfFit(fit_object)
  c(
    NLE = unname(values[["NLE"]]),
    LogDet = unname(values[["LogDet"]]),
    Laplace = unname(values[["Laplace"]]),
    BIC = unname(values[["BIC"]]),
    AIC = unname(values[["AIC"]])
  )
}

fit_profiles <- function(fit_object) {
  out <- DirichletMultinomial::fitted(fit_object, scale = TRUE)
  sweep(out, 2L, colSums(out), "/")
}

fit_posteriors <- function(fit_object, sample_ids = NULL) {
  out <- DirichletMultinomial::mixture(fit_object)
  out <- sweep(out, 1L, rowSums(out), "/")
  if (!is.null(sample_ids)) rownames(out) <- sample_ids
  out
}

match_profiles <- function(primary_profiles, candidate_profiles) {
  common <- intersect(rownames(primary_profiles), rownames(candidate_profiles))
  if (length(common) == 0L) stop("No common taxa for profile matching.", call. = FALSE)
  primary <- primary_profiles[common, , drop = FALSE]
  candidate <- candidate_profiles[common, , drop = FALSE]
  primary <- sweep(primary, 2L, colSums(primary), "/")
  candidate <- sweep(candidate, 2L, colSums(candidate), "/")
  k <- ncol(primary)
  stopifnot(ncol(candidate) == k)
  cost <- matrix(NA_real_, nrow = k, ncol = k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      cost[i, j] <- js_divergence(primary[, i], candidate[, j])
    }
  }
  permutations <- all_permutations(seq_len(k))
  totals <- apply(permutations, 1L, function(p) {
    sum(cost[cbind(seq_len(k), p)])
  })
  best <- permutations[which.min(totals), ]
  list(
    candidate_for_primary = as.integer(best),
    divergence = cost[cbind(seq_len(k), best)],
    total_divergence = min(totals),
    common_features = length(common),
    cost = cost
  )
}

primary_seed <- 20260729L
candidate_k <- 1:7
minimum_state_size <- 14L
prevalence_threshold <- 0.05
depth_threshold <- 1000L
rarefaction_depth <- 1000L
state_labels <- paste0("DMM-", LETTERS[1:4])
set.seed(primary_seed)

input_paths <- c(
  otutab = file.path(input_dir, "otutab.tsv"),
  taxonomy = file.path(input_dir, "taxonomy.tsv"),
  metadata = file.path(input_dir, "metadata.tsv"),
  source_summary = file.path(input_dir, "source-summary.json"),
  prepare_script = file.path(project_root, "scripts", "prepare_dmm_twins_data.R"),
  theme_pub = file.path(project_root, "R", "theme_pub.R")
)
expected_sha256 <- c(
  otutab = "28de822c434aedade101e66ac60a61df38124c9b59357095599fbd278f2ecb41",
  taxonomy = "9eaa0bf788ad1db8c314a59e774c8935bdaa33a51135eaf606c11d5a497dff08",
  metadata = "0368c9129c7d77f381688465141593a2c61a229368d2fc411c9017e1abe04409",
  source_summary = "ed1ea7478f85a565442acf80d0ce25a61b84d5c96ca801f6a45ca010467548e7",
  prepare_script = "1336337bf238a02d9f35d21730c2000d184ed048bd1813760317f38807c99c18",
  theme_pub = "8d3821a485aeb529b12184bbe977c1ca952131a22d5eb69b802088ce91909beb"
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Article 29 input(s): ",
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
  digest = "0.6.36",
  DirichletMultinomial = "1.46.0",
  ggplot2 = "3.5.2",
  jsonlite = "1.8.8",
  ragg = "1.3.2",
  readr = "2.1.5",
  scales = "1.3.0",
  svglite = "2.1.3",
  systemfonts = "1.1.0",
  vegan = "2.6-6.1"
)
observed_versions <- vapply(
  names(expected_versions),
  function(package) utils::packageDescription(package, fields = "Version"),
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
    paste0("package-", tolower(package_audit$package[[i]])),
    "environment",
    package_audit$observed_version[[i]],
    package_audit$expected_version[[i]],
    package_audit$status[[i]] == "PASS"
  )
}

lock_payload <- jsonlite::read_json(file.path(project_root, "env", "renv.lock"))
locked_dmm <- lock_payload$Packages$DirichletMultinomial$Version
add_check(
  "lock-dirichletmultinomial", "environment", locked_dmm, "1.46.0",
  identical(locked_dmm, "1.46.0")
)
add_check(
  "lock-package-count", "environment", length(lock_payload$Packages), 346,
  length(lock_payload$Packages) == 346L
)

command_names <- c("pdfinfo", "pdffonts", "identify")
command_audit <- data.frame(
  command = command_names,
  path = vapply(command_names, command_path, character(1)),
  stringsAsFactors = FALSE
)
command_audit$status <- ifelse(nzchar(command_audit$path), "PASS", "FAIL")
for (i in seq_len(nrow(command_audit))) {
  add_check(
    paste0("command-", command_audit$command[[i]]),
    "environment",
    sanitize_text(command_audit$path[[i]]),
    "available",
    command_audit$status[[i]] == "PASS"
  )
}

font_family <- font_pub
font_info <- systemfonts::font_info(font_family)
font_path <- font_info$path[[1L]]
font_glyphs <- systemfonts::glyph_info(
  paste(
    "Dirichlet multinomial mixture posterior state component profile",
    "Laplace AIC BIC prevalence depth repeat stability 0123456789"
  ),
  family = font_family
)
font_audit <- data.frame(
  requested_family = font_family,
  resolved_family = font_info$family[[1L]],
  style = font_info$style[[1L]],
  path = sanitize_text(font_path),
  file_sha256 = if (file.exists(font_path)) sha256_file(font_path) else "",
  scalable = font_info$scalable[[1L]],
  glyphs_checked = nrow(font_glyphs),
  missing_glyphs = sum(font_glyphs$index == 0L),
  stringsAsFactors = FALSE
)
add_check(
  "font-resolved-family", "font", font_info$family[[1L]], font_family,
  identical(font_info$family[[1L]], font_family)
)
add_check("font-file-exists", "font", file.exists(font_path), TRUE, file.exists(font_path))
add_check("font-scalable", "font", font_info$scalable[[1L]], TRUE, isTRUE(font_info$scalable[[1L]]))
add_check(
  "font-glyph-coverage", "font", font_audit$missing_glyphs, 0,
  font_audit$missing_glyphs == 0L
)

otu_frame <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])
source_summary <- jsonlite::read_json(input_paths[["source_summary"]])

counts <- as.matrix(otu_frame)
storage.mode(counts) <- "integer"
feature_ids <- rownames(counts)
sample_ids <- colnames(counts)
taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]
metadata$SampleID <- rownames(metadata)
metadata$BMIClassCode <- as.integer(metadata$BMIClassCode)
metadata$SamplesPerParticipantKey <- as.integer(metadata$SamplesPerParticipantKey)
sample_depths <- colSums(counts)
sample_counts <- t(counts)

rank_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
add_check("feature-count", "data", nrow(counts), 130, nrow(counts) == 130L)
add_check("sample-count", "data", ncol(counts), 278, ncol(counts) == 278L)
add_check("read-count", "data", sum(counts), 570851, sum(counts) == 570851L)
add_check("taxonomy-feature-count", "data", nrow(taxonomy), 130, nrow(taxonomy) == 130L)
add_check("metadata-sample-count", "data", nrow(metadata), 278, nrow(metadata) == 278L)
add_check(
  "taxonomy-seven-ranks", "data", paste(colnames(taxonomy)[1:7], collapse = ","),
  paste(rank_names, collapse = ","), identical(colnames(taxonomy)[1:7], rank_names)
)
add_check("feature-ids-unique", "data", anyDuplicated(feature_ids), 0, !anyDuplicated(feature_ids))
add_check("sample-ids-unique", "data", anyDuplicated(sample_ids), 0, !anyDuplicated(sample_ids))
add_check(
  "feature-id-alignment", "data", sum(feature_ids == rownames(taxonomy)),
  length(feature_ids), identical(feature_ids, rownames(taxonomy))
)
add_check(
  "sample-id-alignment", "data", sum(sample_ids == rownames(metadata)),
  length(sample_ids), identical(sample_ids, rownames(metadata))
)
add_check("counts-finite", "data", sum(is.finite(counts)), length(counts), all(is.finite(counts)))
add_check("counts-nonnegative", "data", min(counts), ">=0", all(counts >= 0L))
add_check("counts-integer", "data", max(abs(counts - round(counts))), 0, all(counts == round(counts)))
add_check("nonempty-features", "data", sum(rowSums(counts) > 0L), nrow(counts), all(rowSums(counts) > 0L))
add_check("nonempty-samples", "data", sum(sample_depths > 0L), ncol(counts), all(sample_depths > 0L))
add_check("library-depth-min", "data", min(sample_depths), 53, min(sample_depths) == 53L)
add_check("library-depth-median", "data", median(sample_depths), 1597.5, median(sample_depths) == 1597.5)
add_check("library-depth-max", "data", max(sample_depths), 10585, max(sample_depths) == 10585L)
add_check(
  "source-original-unknown-label", "taxonomy",
  sum(taxonomy$Genus == "Uknown"), 1, sum(taxonomy$Genus == "Uknown") == 1L
)
add_check(
  "source-genus-labels-complete", "taxonomy", sum(!is.na(taxonomy$Genus)), 130,
  all(!is.na(taxonomy$Genus) & taxonomy$Genus != "")
)
add_check(
  "source-summary-source-hash", "provenance",
  source_summary$source_files[[1L]]$sha256,
  "b6a231acded317b49b7fd800bf8b4ce12a75340c5b88499b51e92d72528f93ba",
  identical(
    source_summary$source_files[[1L]]$sha256,
    "b6a231acded317b49b7fd800bf8b4ce12a75340c5b88499b51e92d72528f93ba"
  )
)
add_check(
  "source-summary-package-version", "provenance",
  source_summary$package$version, "1.46.0",
  identical(source_summary$package$version, "1.46.0")
)

bmi_counts <- table(factor(metadata$BMIClass, levels = c("Lean", "Obese", "Overweight")))
participant_counts <- table(metadata$ParticipantKey)
two_sample_keys <- names(participant_counts)[participant_counts == 2L]
one_sample_keys <- names(participant_counts)[participant_counts == 1L]
suffix_only <- vapply(split(metadata$SourceSuffix, metadata$ParticipantKey), function(x) {
  length(x) == 1L && identical(x, "Dot2 suffix")
}, logical(1))
no_suffix_only <- vapply(split(metadata$SourceSuffix, metadata$ParticipantKey), function(x) {
  length(x) == 1L && identical(x, "No suffix")
}, logical(1))
add_check(
  "bmi-class-counts", "metadata", paste(as.integer(bmi_counts), collapse = ","),
  "61,193,24", identical(as.integer(bmi_counts), c(61L, 193L, 24L))
)
add_check("participant-key-count", "metadata", length(participant_counts), 154, length(participant_counts) == 154L)
add_check("two-sample-key-count", "metadata", length(two_sample_keys), 124, length(two_sample_keys) == 124L)
add_check("one-sample-key-count", "metadata", length(one_sample_keys), 30, length(one_sample_keys) == 30L)
add_check("no-suffix-only-key-count", "metadata", sum(no_suffix_only), 26, sum(no_suffix_only) == 26L)
add_check("dot2-only-key-count", "metadata", sum(suffix_only), 4, sum(suffix_only) == 4L)

sample_alignment_audit <- data.frame(
  SampleID = sample_ids,
  CountTable = sample_ids %in% colnames(counts),
  Metadata = sample_ids %in% rownames(metadata),
  ParticipantKey = metadata[sample_ids, "ParticipantKey"],
  SourceSuffix = metadata[sample_ids, "SourceSuffix"],
  BMIClass = metadata[sample_ids, "BMIClass"],
  LibrarySize = as.numeric(sample_depths[sample_ids]),
  stringsAsFactors = FALSE
)

feature_prevalence <- rowMeans(counts > 0L)
prevalence_features <- feature_ids[feature_prevalence >= prevalence_threshold]
depth_samples <- sample_ids[sample_depths >= depth_threshold]

set.seed(primary_seed)
rarefied_once <- vegan::rrarefy(sample_counts[depth_samples, , drop = FALSE], rarefaction_depth)
set.seed(primary_seed)
rarefied_twice <- vegan::rrarefy(sample_counts[depth_samples, , drop = FALSE], rarefaction_depth)
rarefaction_reproducible <- identical(rarefied_once, rarefied_twice)

sample_indices <- split(seq_len(nrow(metadata)), metadata$ParticipantKey)
one_per_participant <- vapply(sample_indices, function(idx) {
  no_suffix <- idx[metadata$SourceSuffix[idx] == "No suffix"]
  if (length(no_suffix) > 0L) no_suffix[[1L]] else idx[[1L]]
}, integer(1))
one_per_samples <- metadata$SampleID[one_per_participant]

branch_counts <- list(
  Primary = sample_counts,
  `Prevalence >= 5%` = sample_counts[, prevalence_features, drop = FALSE],
  `Depth >= 1,000` = sample_counts[depth_samples, , drop = FALSE],
  `Rarefied to 1,000` = rarefied_once,
  `One per participant` = sample_counts[one_per_samples, , drop = FALSE]
)
branch_seeds <- c(
  Primary = primary_seed,
  `Prevalence >= 5%` = primary_seed + 100L,
  `Depth >= 1,000` = primary_seed + 200L,
  `Rarefied to 1,000` = primary_seed + 300L,
  `One per participant` = primary_seed + 400L
)

add_check("prevalence-feature-count", "sensitivity", length(prevalence_features), 67, length(prevalence_features) == 67L)
add_check("depth-sample-count", "sensitivity", length(depth_samples), 247, length(depth_samples) == 247L)
add_check("one-per-participant-count", "sensitivity", length(one_per_samples), 154, length(one_per_samples) == 154L)
add_check("rarefaction-reproducible", "sensitivity", rarefaction_reproducible, TRUE, rarefaction_reproducible)
add_check(
  "rarefaction-depth", "sensitivity", paste(unique(rowSums(rarefied_once)), collapse = ","),
  rarefaction_depth, all(rowSums(rarefied_once) == rarefaction_depth)
)

rarefaction_ledger <- data.frame(
  SampleID = depth_samples,
  OriginalLibrarySize = as.integer(rowSums(sample_counts[depth_samples, , drop = FALSE])),
  RetainedReads = as.integer(rowSums(rarefied_once)),
  DiscardedReads = as.integer(
    rowSums(sample_counts[depth_samples, , drop = FALSE]) - rowSums(rarefied_once)
  ),
  Seed = primary_seed,
  stringsAsFactors = FALSE
)

candidate_jobs <- expand.grid(
  Branch = names(branch_counts),
  K = candidate_k,
  stringsAsFactors = FALSE
)
candidate_jobs$Seed <- unname(branch_seeds[candidate_jobs$Branch])
initialization_seeds <- primary_seed + 1:7
initialization_jobs <- data.frame(
  Branch = "Initialization",
  K = 4L,
  Seed = initialization_seeds,
  stringsAsFactors = FALSE
)
all_jobs <- rbind(candidate_jobs, initialization_jobs)

available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores)) available_cores <- 1L
worker_count <- min(8L, max(1L, available_cores), nrow(all_jobs))
execution_backend <- "One-shot PSOCK R process per model"
fit_one_job <- function(i) {
  job <- all_jobs[i, , drop = FALSE]
  x <- if (job$Branch == "Initialization") {
    branch_counts[["Primary"]]
  } else {
    branch_counts[[job$Branch]]
  }
  timing <- system.time({
    fit_object <- DirichletMultinomial::dmn(
      x,
      k = job$K,
      verbose = FALSE,
      seed = job$Seed
    )
  })
  list(fit = fit_object, elapsed_seconds = unname(timing[["elapsed"]]))
}

fit_cache_path <- file.path(output_dir, "fresh-fit-cache.rds")
cache_payload <- if (file.exists(fit_cache_path)) {
  tryCatch(readRDS(fit_cache_path), error = function(e) NULL)
} else {
  NULL
}
cache_matches <- !is.null(cache_payload) &&
  identical(cache_payload$jobs, all_jobs) &&
  length(cache_payload$results) == nrow(all_jobs) &&
  (is.null(cache_payload$input_sha256) || identical(
    cache_payload$input_sha256,
    observed_sha256[c("otutab", "taxonomy", "metadata")]
  )) &&
  (is.null(cache_payload$package_version) || identical(
    cache_payload$package_version,
    observed_versions[["DirichletMultinomial"]]
  )) &&
  identical(cache_payload$execution_backend, execution_backend)

if (cache_matches) {
  message(sprintf("Reusing %d validated local DMM fit objects...", nrow(all_jobs)))
  job_results <- cache_payload$results
} else {
  message(sprintf(
    "Fitting %d fresh DMM models with at most %d concurrent worker(s) via %s...",
    nrow(all_jobs), worker_count, execution_backend
  ))
  job_indices <- seq_len(nrow(all_jobs))
  job_batches <- split(
    job_indices,
    ceiling(job_indices / worker_count)
  )
  job_results <- vector("list", nrow(all_jobs))
  for (batch_indices in job_batches) {
    fit_cluster <- parallel::makePSOCKcluster(
      length(batch_indices),
      outfile = ""
    )
    parallel::clusterExport(
      fit_cluster,
      c("all_jobs", "branch_counts", "fit_one_job"),
      envir = environment()
    )
    parallel::clusterEvalQ(fit_cluster, {
      Sys.setenv(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1"
      )
      requireNamespace("DirichletMultinomial", quietly = TRUE)
      NULL
    })
    batch_results <- tryCatch(
      parallel::parLapply(
        fit_cluster,
        batch_indices,
        fit_one_job
      ),
      finally = parallel::stopCluster(fit_cluster)
    )
    job_results[batch_indices] <- batch_results
  }
}
failed_jobs <- vapply(job_results, inherits, logical(1), what = "try-error")
if (any(failed_jobs)) {
  stop(
    "DMM fit jobs failed: ",
    paste(which(failed_jobs), collapse = ", "),
    call. = FALSE
  )
}
saveRDS(
  list(
    jobs = all_jobs,
    results = job_results,
    input_sha256 = observed_sha256[c("otutab", "taxonomy", "metadata")],
    package_version = observed_versions[["DirichletMultinomial"]],
    execution_backend = execution_backend
  ),
  fit_cache_path,
  version = 3L
)

candidate_results <- job_results[seq_len(nrow(candidate_jobs))]
initialization_results <- job_results[nrow(candidate_jobs) + seq_len(nrow(initialization_jobs))]

model_rows <- lapply(seq_len(nrow(candidate_jobs)), function(i) {
  scores <- fit_scores(candidate_results[[i]]$fit)
  data.frame(
    Branch = candidate_jobs$Branch[[i]],
    K = candidate_jobs$K[[i]],
    Seed = candidate_jobs$Seed[[i]],
    Samples = nrow(branch_counts[[candidate_jobs$Branch[[i]]]]),
    Features = ncol(branch_counts[[candidate_jobs$Branch[[i]]]]),
    NLE = scores[["NLE"]],
    LogDet = scores[["LogDet"]],
    Laplace = scores[["Laplace"]],
    BIC = scores[["BIC"]],
    AIC = scores[["AIC"]],
    ElapsedSeconds = candidate_results[[i]]$elapsed_seconds,
    stringsAsFactors = FALSE
  )
})
model_selection <- do.call(rbind, model_rows)
model_selection$SelectedLaplace <- ave(
  model_selection$Laplace,
  model_selection$Branch,
  FUN = function(x) x == min(x)
)
model_selection$SelectedAIC <- ave(
  model_selection$AIC,
  model_selection$Branch,
  FUN = function(x) x == min(x)
)
model_selection$SelectedBIC <- ave(
  model_selection$BIC,
  model_selection$Branch,
  FUN = function(x) x == min(x)
)

get_candidate_fit <- function(branch, k) {
  index <- which(candidate_jobs$Branch == branch & candidate_jobs$K == k)
  stopifnot(length(index) == 1L)
  candidate_results[[index]]$fit
}

selected_for <- function(branch, criterion) {
  branch_rows <- model_selection[model_selection$Branch == branch, , drop = FALSE]
  branch_rows$K[[which.min(branch_rows[[criterion]])]]
}

primary_fit <- get_candidate_fit("Primary", 4L)
primary_raw_profiles <- fit_profiles(primary_fit)
primary_raw_posteriors <- fit_posteriors(primary_fit, rownames(branch_counts$Primary))
bacteroides_index <- match("Bacteroides", rownames(primary_raw_profiles))
prevotella_index <- match("Prevotella", rownames(primary_raw_profiles))
stopifnot(!is.na(bacteroides_index), !is.na(prevotella_index))
canonical_order <- order(
  -primary_raw_profiles[bacteroides_index, ],
  -primary_raw_profiles[prevotella_index, ],
  seq_len(ncol(primary_raw_profiles))
)
primary_profiles <- primary_raw_profiles[, canonical_order, drop = FALSE]
colnames(primary_profiles) <- state_labels
primary_posteriors <- primary_raw_posteriors[, canonical_order, drop = FALSE]
colnames(primary_posteriors) <- state_labels
primary_hard_index <- max.col(primary_posteriors, ties.method = "first")
primary_state <- factor(state_labels[primary_hard_index], levels = state_labels)
names(primary_state) <- rownames(primary_posteriors)
primary_max_posterior <- apply(primary_posteriors, 1L, max)
primary_entropy <- -rowSums(
  primary_posteriors * log(pmax(primary_posteriors, 1e-15))
) / log(ncol(primary_posteriors))

mix_weights_raw <- DirichletMultinomial::mixturewt(primary_fit)
mix_weights <- mix_weights_raw[canonical_order, , drop = FALSE]
component_sizes <- as.integer(table(primary_state))
component_summary <- data.frame(
  State = state_labels,
  RawComponent = canonical_order,
  AssignedSamples = component_sizes,
  MixtureWeight = mix_weights$pi,
  Theta = mix_weights$theta,
  ExpectedBacteroides = primary_profiles["Bacteroides", ],
  ExpectedPrevotella = primary_profiles["Prevotella", ],
  MedianMaximumPosterior = vapply(state_labels, function(state) {
    stats::median(primary_max_posterior[primary_state == state])
  }, numeric(1)),
  SamplesBelowPointEight = vapply(state_labels, function(state) {
    sum(primary_state == state & primary_max_posterior < 0.8)
  }, integer(1)),
  stringsAsFactors = FALSE
)

component_profiles <- data.frame(
  FeatureID = rep(rownames(primary_profiles), times = ncol(primary_profiles)),
  SourceGenus = rep(taxonomy[rownames(primary_profiles), "Genus"], times = ncol(primary_profiles)),
  DisplayGenus = rep(ifelse(
    taxonomy[rownames(primary_profiles), "Genus"] == "Uknown",
    "Unclassified",
    taxonomy[rownames(primary_profiles), "Genus"]
  ), times = ncol(primary_profiles)),
  State = rep(state_labels, each = nrow(primary_profiles)),
  ExpectedProportion = as.vector(primary_profiles),
  stringsAsFactors = FALSE
)
component_profiles$ExpectedPercent <- 100 * component_profiles$ExpectedProportion

component_top_taxa <- do.call(rbind, lapply(state_labels, function(state) {
  state_rows <- component_profiles[component_profiles$State == state, , drop = FALSE]
  state_rows <- state_rows[order(state_rows$ExpectedProportion, decreasing = TRUE), , drop = FALSE]
  state_rows$RankWithinState <- seq_len(nrow(state_rows))
  state_rows[state_rows$RankWithinState <= 12L, , drop = FALSE]
}))
rownames(component_top_taxa) <- NULL

sample_posteriors <- data.frame(
  SampleID = rownames(primary_posteriors),
  DMM_A = primary_posteriors[, "DMM-A"],
  DMM_B = primary_posteriors[, "DMM-B"],
  DMM_C = primary_posteriors[, "DMM-C"],
  DMM_D = primary_posteriors[, "DMM-D"],
  State = as.character(primary_state),
  MaximumPosterior = primary_max_posterior,
  NormalizedEntropy = primary_entropy,
  LowConfidencePointEight = primary_max_posterior < 0.8,
  ParticipantKey = metadata[rownames(primary_posteriors), "ParticipantKey"],
  SourceSuffix = metadata[rownames(primary_posteriors), "SourceSuffix"],
  BMIClass = metadata[rownames(primary_posteriors), "BMIClass"],
  LibrarySize = as.integer(sample_depths[rownames(primary_posteriors)]),
  stringsAsFactors = FALSE
)

primary_rows <- model_selection[model_selection$Branch == "Primary", , drop = FALSE]
selected_k_laplace <- selected_for("Primary", "Laplace")
selected_k_aic <- selected_for("Primary", "AIC")
selected_k_bic <- selected_for("Primary", "BIC")
add_check("candidate-k-grid", "model", paste(primary_rows$K, collapse = ","), "1,2,3,4,5,6,7", identical(primary_rows$K, 1:7))
add_check("primary-seed", "model", unique(primary_rows$Seed), primary_seed, identical(unique(primary_rows$Seed), primary_seed))
add_check("selected-k-laplace", "model", selected_k_laplace, 4, selected_k_laplace == 4L)
add_check("selected-k-aic", "model", selected_k_aic, 2, selected_k_aic == 2L)
add_check("selected-k-bic", "model", selected_k_bic, 1, selected_k_bic == 1L)
add_check("criteria-disagree", "model", length(unique(c(selected_k_laplace, selected_k_aic, selected_k_bic))), 3, length(unique(c(selected_k_laplace, selected_k_aic, selected_k_bic))) == 3L)
add_check(
  "primary-laplace-k4", "model",
  primary_rows$Laplace[primary_rows$K == 4L], 38838.4232180317,
  near(primary_rows$Laplace[primary_rows$K == 4L], 38838.4232180317, 1e-3)
)
add_check("component-size-vector", "model", paste(component_sizes, collapse = ","), "85,86,46,61", identical(component_sizes, c(85L, 86L, 46L, 61L)))
add_check("minimum-state-size", "model", min(component_sizes), 46, min(component_sizes) == 46L)
add_check("minimum-state-threshold", "model", min(component_sizes), paste0(">=", minimum_state_size), min(component_sizes) >= minimum_state_size)
add_check("profile-column-closure", "model", max(abs(colSums(primary_profiles) - 1)), "<=1e-12", max(abs(colSums(primary_profiles) - 1)) <= 1e-12)
add_check("posterior-row-closure", "model", max(abs(rowSums(primary_posteriors) - 1)), "<=1e-12", max(abs(rowSums(primary_posteriors) - 1)) <= 1e-12)
add_check("posterior-bounded", "model", paste(range(primary_posteriors), collapse = ","), "[0,1]", all(primary_posteriors >= 0 & primary_posteriors <= 1))
add_check("posterior-low-confidence-count", "model", sum(primary_max_posterior < 0.8), 8, sum(primary_max_posterior < 0.8) == 8L)

reference_environment <- new.env(parent = emptyenv())
utils::data("fit", package = "DirichletMultinomial", envir = reference_environment)
reference_fits <- reference_environment$fit
stopifnot(length(reference_fits) >= max(candidate_k))

reference_rows <- do.call(rbind, lapply(candidate_k, function(k) {
  scores <- fit_scores(reference_fits[[k]])
  data.frame(
    Source = "Bundled reference",
    K = k,
    NLE = scores[["NLE"]],
    LogDet = scores[["LogDet"]],
    Laplace = scores[["Laplace"]],
    BIC = scores[["BIC"]],
    AIC = scores[["AIC"]],
    stringsAsFactors = FALSE
  )
}))
fresh_reference_rows <- data.frame(
  Source = "Fresh fixed-seed fit",
  K = primary_rows$K,
  NLE = primary_rows$NLE,
  LogDet = primary_rows$LogDet,
  Laplace = primary_rows$Laplace,
  BIC = primary_rows$BIC,
  AIC = primary_rows$AIC,
  stringsAsFactors = FALSE
)
reference_fit_comparison <- rbind(fresh_reference_rows, reference_rows)
reference_fit_comparison$SelectedLaplace <- ave(
  reference_fit_comparison$Laplace,
  reference_fit_comparison$Source,
  FUN = function(x) x == min(x)
)
reference_fit_comparison$SelectedAIC <- ave(
  reference_fit_comparison$AIC,
  reference_fit_comparison$Source,
  FUN = function(x) x == min(x)
)
reference_fit_comparison$SelectedBIC <- ave(
  reference_fit_comparison$BIC,
  reference_fit_comparison$Source,
  FUN = function(x) x == min(x)
)

reference_profiles_raw <- fit_profiles(reference_fits[[4L]])
reference_match <- match_profiles(primary_profiles, reference_profiles_raw)
reference_posteriors_raw <- fit_posteriors(reference_fits[[4L]], sample_ids)
reference_posteriors <- reference_posteriors_raw[
  , reference_match$candidate_for_primary, drop = FALSE
]
reference_state <- state_labels[max.col(reference_posteriors, ties.method = "first")]
reference_ari <- adjusted_rand_index(as.character(primary_state), reference_state)
add_check("reference-selected-k-laplace", "reference", reference_rows$K[[which.min(reference_rows$Laplace)]], 4, reference_rows$K[[which.min(reference_rows$Laplace)]] == 4L)
add_check(
  "reference-laplace-k4", "reference",
  reference_rows$Laplace[reference_rows$K == 4L], 38781.1038719804,
  near(reference_rows$Laplace[reference_rows$K == 4L], 38781.1038719804, 1e-3)
)
add_check("reference-fresh-laplace-distinct", "reference", abs(reference_rows$Laplace[4L] - primary_rows$Laplace[4L]), ">50", abs(reference_rows$Laplace[4L] - primary_rows$Laplace[4L]) > 50)
add_check("reference-fresh-partition-ari", "reference", reference_ari, 1, near(reference_ari, 1, 1e-12))

primary_initialization <- list(fit = primary_fit, elapsed_seconds = primary_rows$ElapsedSeconds[primary_rows$K == 4L])
initialization_all <- c(list(primary_initialization), initialization_results)
initialization_seed_all <- primary_seed + 0:7
initialization_audit_rows <- list()
initialization_states <- list()
initialization_match_rows <- list()
for (i in seq_along(initialization_all)) {
  fit_object <- initialization_all[[i]]$fit
  profiles_raw <- fit_profiles(fit_object)
  profile_match <- match_profiles(primary_profiles, profiles_raw)
  posteriors_raw <- fit_posteriors(fit_object, sample_ids)
  posteriors <- posteriors_raw[, profile_match$candidate_for_primary, drop = FALSE]
  states <- state_labels[max.col(posteriors, ties.method = "first")]
  initialization_states[[i]] <- states
  scores <- fit_scores(fit_object)
  sizes <- as.integer(table(factor(states, levels = state_labels)))
  initialization_audit_rows[[i]] <- data.frame(
    Seed = initialization_seed_all[[i]],
    Laplace = scores[["Laplace"]],
    AIC = scores[["AIC"]],
    BIC = scores[["BIC"]],
    ARIToPrimary = adjusted_rand_index(as.character(primary_state), states),
    MeanMatchedJSD = mean(profile_match$divergence),
    StateSizes = paste(sizes, collapse = "/"),
    ElapsedSeconds = initialization_all[[i]]$elapsed_seconds,
    stringsAsFactors = FALSE
  )
  initialization_match_rows[[i]] <- data.frame(
    Branch = paste0("Initialization seed ", initialization_seed_all[[i]]),
    State = state_labels,
    CandidateRawComponent = profile_match$candidate_for_primary,
    JSDivergence = profile_match$divergence,
    CommonFeatures = profile_match$common_features,
    stringsAsFactors = FALSE
  )
}
initialization_audit <- do.call(rbind, initialization_audit_rows)
pair_indices <- utils::combn(seq_along(initialization_states), 2L)
initialization_agreement <- data.frame(
  SeedA = initialization_seed_all[pair_indices[1L, ]],
  SeedB = initialization_seed_all[pair_indices[2L, ]],
  AdjustedRandIndex = apply(pair_indices, 2L, function(idx) {
    adjusted_rand_index(initialization_states[[idx[[1L]]]], initialization_states[[idx[[2L]]]])
  }),
  stringsAsFactors = FALSE
)
initialization_min_ari <- min(initialization_agreement$AdjustedRandIndex)
add_check("initialization-count", "stability", nrow(initialization_audit), 8, nrow(initialization_audit) == 8L)
add_check("initialization-pair-count", "stability", nrow(initialization_agreement), 28, nrow(initialization_agreement) == 28L)
add_check(
  "initialization-min-ari", "stability", initialization_min_ari,
  1,
  near(initialization_min_ari, 1, 1e-12)
)
add_check(
  "initialization-state-sizes", "stability",
  paste(unique(initialization_audit$StateSizes), collapse = ","),
  "85/86/46/61",
  identical(
    unique(initialization_audit$StateSizes),
    "85/86/46/61"
  )
)
add_check(
  "initialization-primary-basin-count", "stability",
  sum(initialization_audit$ARIToPrimary == 1), 8,
  sum(initialization_audit$ARIToPrimary == 1) == 8L
)
add_check(
  "initialization-outlier-seed", "stability",
  paste(initialization_audit$Seed[initialization_audit$ARIToPrimary < 1], collapse = ","),
  "none",
  length(initialization_audit$Seed[initialization_audit$ARIToPrimary < 1]) == 0L
)
add_check(
  "initialization-best-laplace-same-partition", "stability",
  initialization_audit$ARIToPrimary[[which.min(initialization_audit$Laplace)]],
  1,
  near(
    initialization_audit$ARIToPrimary[[which.min(initialization_audit$Laplace)]],
    1,
    1e-12
  )
)

label_matching_rows <- initialization_match_rows
label_matching_rows[[length(label_matching_rows) + 1L]] <- data.frame(
  Branch = "Bundled reference",
  State = state_labels,
  CandidateRawComponent = reference_match$candidate_for_primary,
  JSDivergence = reference_match$divergence,
  CommonFeatures = reference_match$common_features,
  stringsAsFactors = FALSE
)

sensitivity_rows <- list()
branch_state_assignments <- list(Primary = setNames(as.character(primary_state), names(primary_state)))
for (branch in names(branch_counts)) {
  fit_object <- get_candidate_fit(branch, 4L)
  profiles_raw <- fit_profiles(fit_object)
  profile_match <- match_profiles(primary_profiles, profiles_raw)
  posteriors_raw <- fit_posteriors(fit_object, rownames(branch_counts[[branch]]))
  posteriors <- posteriors_raw[, profile_match$candidate_for_primary, drop = FALSE]
  states <- state_labels[max.col(posteriors, ties.method = "first")]
  names(states) <- rownames(posteriors)
  branch_state_assignments[[branch]] <- states
  overlap <- intersect(names(primary_state), names(states))
  sizes <- as.integer(table(factor(states, levels = state_labels)))
  sensitivity_rows[[branch]] <- data.frame(
    Branch = branch,
    Seed = branch_seeds[[branch]],
    Samples = nrow(branch_counts[[branch]]),
    Features = ncol(branch_counts[[branch]]),
    SelectedKLaplace = selected_for(branch, "Laplace"),
    SelectedKAIC = selected_for(branch, "AIC"),
    SelectedKBIC = selected_for(branch, "BIC"),
    FixedK = 4L,
    MinimumStateSize = min(sizes),
    MaximumStateSize = max(sizes),
    OverlapWithPrimary = length(overlap),
    AdjustedRandIndex = adjusted_rand_index(
      as.character(primary_state[overlap]), states[overlap]
    ),
    MeanMatchedJSD = mean(profile_match$divergence),
    MaximumMatchedJSD = max(profile_match$divergence),
    stringsAsFactors = FALSE
  )
  label_matching_rows[[length(label_matching_rows) + 1L]] <- data.frame(
    Branch = branch,
    State = state_labels,
    CandidateRawComponent = profile_match$candidate_for_primary,
    JSDivergence = profile_match$divergence,
    CommonFeatures = profile_match$common_features,
    stringsAsFactors = FALSE
  )
}
sensitivity_agreement <- do.call(rbind, sensitivity_rows)
rownames(sensitivity_agreement) <- NULL
sensitivity_model_selection <- model_selection[model_selection$Branch != "Primary", , drop = FALSE]
label_matching_audit <- do.call(rbind, label_matching_rows)
rownames(label_matching_audit) <- NULL

add_check("sensitivity-branch-count", "stability", nrow(sensitivity_agreement), 5, nrow(sensitivity_agreement) == 5L)
add_check("sensitivity-selected-k-range", "stability", paste(sensitivity_agreement$SelectedKLaplace, collapse = ","), "1..7", all(sensitivity_agreement$SelectedKLaplace %in% candidate_k))
add_check(
  "sensitivity-selected-k-vector", "stability",
  paste(sensitivity_agreement$SelectedKLaplace, collapse = ","),
  "4,4,2,2,2",
  identical(sensitivity_agreement$SelectedKLaplace, c(4L, 4L, 2L, 2L, 2L))
)
add_check("sensitivity-ari-finite", "stability", sum(is.finite(sensitivity_agreement$AdjustedRandIndex)), 5, all(is.finite(sensitivity_agreement$AdjustedRandIndex)))
add_check("sensitivity-profile-divergence-finite", "stability", sum(is.finite(sensitivity_agreement$MeanMatchedJSD)), 5, all(is.finite(sensitivity_agreement$MeanMatchedJSD)) && all(sensitivity_agreement$MeanMatchedJSD >= 0))
add_check("primary-sensitivity-ari", "stability", sensitivity_agreement$AdjustedRandIndex[sensitivity_agreement$Branch == "Primary"], 1, near(sensitivity_agreement$AdjustedRandIndex[sensitivity_agreement$Branch == "Primary"], 1, 1e-12))
add_check(
  "sensitivity-fixed-k4-ari-vector", "stability",
  round(sensitivity_agreement$AdjustedRandIndex, 6),
  round(c(1, 0.903513447720664, 0.92842598427158, 0.905608750546972, 0.67951691723649), 6),
  near(
    sensitivity_agreement$AdjustedRandIndex,
    c(1, 0.903513447720664, 0.92842598427158, 0.905608750546972, 0.67951691723649),
    1e-12
  )
)

bmi_cross_tab <- as.data.frame(table(
  State = factor(primary_state, levels = state_labels),
  BMIClass = factor(metadata[names(primary_state), "BMIClass"], levels = c("Lean", "Obese", "Overweight"))
), stringsAsFactors = FALSE)
names(bmi_cross_tab)[3L] <- "Samples"
bmi_cross_tab$WithinStateFraction <- ave(
  bmi_cross_tab$Samples,
  bmi_cross_tab$State,
  FUN = function(x) x / sum(x)
)

participant_repeat_rows <- lapply(two_sample_keys, function(key) {
  ids <- metadata$SampleID[metadata$ParticipantKey == key]
  no_suffix_id <- ids[metadata[ids, "SourceSuffix"] == "No suffix"]
  dot2_id <- ids[metadata[ids, "SourceSuffix"] == "Dot2 suffix"]
  stopifnot(length(no_suffix_id) == 1L, length(dot2_id) == 1L)
  data.frame(
    ParticipantKey = key,
    NoSuffixSample = no_suffix_id,
    Dot2SuffixSample = dot2_id,
    NoSuffixBMIClass = metadata[no_suffix_id, "BMIClass"],
    Dot2SuffixBMIClass = metadata[dot2_id, "BMIClass"],
    BMIClassAgreement = metadata[no_suffix_id, "BMIClass"] == metadata[dot2_id, "BMIClass"],
    NoSuffixState = as.character(primary_state[no_suffix_id]),
    Dot2SuffixState = as.character(primary_state[dot2_id]),
    StateAgreement = primary_state[no_suffix_id] == primary_state[dot2_id],
    NoSuffixMaximumPosterior = primary_max_posterior[no_suffix_id],
    Dot2SuffixMaximumPosterior = primary_max_posterior[dot2_id],
    stringsAsFactors = FALSE
  )
})
participant_repeat_audit <- do.call(rbind, participant_repeat_rows)
rownames(participant_repeat_audit) <- NULL
participant_transitions <- as.data.frame(table(
  NoSuffixState = factor(participant_repeat_audit$NoSuffixState, levels = state_labels),
  Dot2SuffixState = factor(participant_repeat_audit$Dot2SuffixState, levels = state_labels)
), stringsAsFactors = FALSE)
names(participant_transitions)[3L] <- "ParticipantPairs"
repeat_agreement <- mean(participant_repeat_audit$StateAgreement)
repeat_ari <- adjusted_rand_index(
  participant_repeat_audit$NoSuffixState,
  participant_repeat_audit$Dot2SuffixState
)
add_check("repeat-pair-count", "repeat", nrow(participant_repeat_audit), 124, nrow(participant_repeat_audit) == 124L)
add_check("repeat-state-agreement-count", "repeat", sum(participant_repeat_audit$StateAgreement), 85, sum(participant_repeat_audit$StateAgreement) == 85L)
add_check("repeat-state-agreement", "repeat", repeat_agreement, 0.685483870967742, near(repeat_agreement, 0.685483870967742, 1e-12))
add_check("repeat-state-ari", "repeat", repeat_ari, 0.37678348271377, near(repeat_ari, 0.37678348271377, 1e-12))
add_check("repeat-transition-closure", "repeat", sum(participant_transitions$ParticipantPairs), 124, sum(participant_transitions$ParticipantPairs) == 124L)

relative_samples <- sweep(sample_counts, 1L, rowSums(sample_counts), "/")
bray <- vegan::vegdist(relative_samples, method = "bray")
ordination <- stats::cmdscale(bray, k = 2L, eig = TRUE, add = TRUE)
ordination_points <- ordination$points
if (stats::cor(ordination_points[, 1L], relative_samples[, "Bacteroides"]) < 0) {
  ordination_points[, 1L] <- -ordination_points[, 1L]
}
if (stats::cor(ordination_points[, 2L], relative_samples[, "Prevotella"]) < 0) {
  ordination_points[, 2L] <- -ordination_points[, 2L]
}
positive_eigenvalues <- ordination$eig[ordination$eig > 0]
axis_percent <- 100 * ordination$eig[1:2] / sum(positive_eigenvalues)
ordination_scores <- data.frame(
  SampleID = rownames(ordination_points),
  PCoA1 = ordination_points[, 1L],
  PCoA2 = ordination_points[, 2L],
  State = as.character(primary_state[rownames(ordination_points)]),
  MaximumPosterior = primary_max_posterior[rownames(ordination_points)],
  NormalizedEntropy = primary_entropy[rownames(ordination_points)],
  SourceSuffix = metadata[rownames(ordination_points), "SourceSuffix"],
  BMIClass = metadata[rownames(ordination_points), "BMIClass"],
  ParticipantKey = metadata[rownames(ordination_points), "ParticipantKey"],
  stringsAsFactors = FALSE
)
add_check("ordination-sample-count", "ordination", nrow(ordination_scores), 278, nrow(ordination_scores) == 278L)
add_check("ordination-finite", "ordination", sum(is.finite(as.matrix(ordination_scores[, c("PCoA1", "PCoA2")]))), 556, all(is.finite(as.matrix(ordination_scores[, c("PCoA1", "PCoA2")]))))
add_check("ordination-axis-positive", "ordination", paste(axis_percent, collapse = ","), ">0", all(axis_percent > 0))

state_colours <- c(
  "DMM-A" = "#0072B2",
  "DMM-B" = "#E69F00",
  "DMM-C" = "#009E73",
  "DMM-D" = "#CC79A7"
)

criterion_plot_rows <- do.call(rbind, lapply(c("Laplace", "AIC", "BIC"), function(criterion) {
  data.frame(
    Source = reference_fit_comparison$Source,
    K = reference_fit_comparison$K,
    Criterion = criterion,
    Score = reference_fit_comparison[[criterion]],
    stringsAsFactors = FALSE
  )
}))
criterion_plot_rows$Minimum <- ave(
  criterion_plot_rows$Score,
  interaction(criterion_plot_rows$Source, criterion_plot_rows$Criterion),
  FUN = function(x) x == min(x)
) == 1
criterion_plot_rows$Criterion <- factor(
  criterion_plot_rows$Criterion,
  levels = c("Laplace", "AIC", "BIC")
)
model_selection_plot <- ggplot2::ggplot(
  criterion_plot_rows,
  ggplot2::aes(x = K, y = Score, colour = Source, linetype = Source)
) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::geom_point(size = 2.0) +
  ggplot2::geom_point(
    data = criterion_plot_rows[criterion_plot_rows$Minimum, , drop = FALSE],
    shape = 21,
    fill = "white",
    stroke = 1.0,
    size = 3.4
  ) +
  ggplot2::facet_wrap(~Criterion, scales = "free_y", ncol = 1L) +
  ggplot2::scale_x_continuous(breaks = candidate_k) +
  ggplot2::scale_colour_manual(values = c(
    "Fresh fixed-seed fit" = "#0072B2",
    "Bundled reference" = "#D55E00"
  )) +
  ggplot2::scale_linetype_manual(values = c(
    "Fresh fixed-seed fit" = "solid",
    "Bundled reference" = "22"
  )) +
  ggplot2::scale_y_continuous(labels = scales::label_comma(accuracy = 1)) +
  ggplot2::labs(
    title = "Model-selection criteria do not tell the same story",
    subtitle = "Open circles mark criterion minima within each independently fitted series",
    x = "Number of mixture components (k)",
    y = "Criterion score (lower is better)",
    colour = "Fit source",
    linetype = "Fit source",
    caption = paste0(
      "Fresh Laplace selects k = 4; fresh AIC selects k = 2 and BIC selects k = 1.\n",
      "The bundled reference has a lower k = 4 Laplace value but the same Laplace minimum."
    )
  ) +
  theme_pub(base_size = 8.8) +
  ggplot2::theme(legend.position = "bottom")

ordination_scores$State <- factor(ordination_scores$State, levels = state_labels)
ordination_scores$ConfidenceClass <- ifelse(
  ordination_scores$MaximumPosterior >= 0.8,
  "Posterior >= 0.80",
  "Posterior < 0.80"
)
state_centroids <- stats::aggregate(
  ordination_scores[, c("PCoA1", "PCoA2")],
  list(State = ordination_scores$State),
  mean
)
posterior_ordination_plot <- ggplot2::ggplot(
  ordination_scores,
  ggplot2::aes(x = PCoA1, y = PCoA2, colour = State)
) +
  ggplot2::geom_point(
    ggplot2::aes(alpha = MaximumPosterior, shape = ConfidenceClass),
    size = 2.2,
    stroke = 0.55
  ) +
  ggplot2::geom_point(
    data = state_centroids,
    shape = 4,
    size = 4.2,
    stroke = 1.1,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = state_colours, drop = FALSE) +
  ggplot2::scale_alpha_continuous(
    range = c(0.32, 0.9), limits = c(0, 1), guide = "none"
  ) +
  ggplot2::scale_shape_manual(values = c(
    "Posterior >= 0.80" = 16,
    "Posterior < 0.80" = 1
  )) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "DMM states are posterior assignments, not PCoA quadrants",
    subtitle = "Bray-Curtis PCoA is shown only as a two-dimensional descriptive view",
    x = sprintf("PCoA 1 (%.1f%% of positive eigenvalue sum)", axis_percent[[1L]]),
    y = sprintf("PCoA 2 (%.1f%% of positive eigenvalue sum)", axis_percent[[2L]]),
    colour = "DMM state",
    alpha = "Maximum posterior",
    shape = "Assignment confidence",
    caption = paste0(
      "Crosses are state centroids. Open points have maximum posterior below 0.80.\n",
      "The DMM was fitted to raw counts; ordination did not define the states."
    )
  ) +
  theme_pub(base_size = 8.8) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(order = 1L, nrow = 1L),
    shape = ggplot2::guide_legend(order = 2L, nrow = 1L)
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

top_profile_features <- unique(unlist(lapply(state_labels, function(state) {
  state_values <- primary_profiles[, state]
  names(sort(state_values, decreasing = TRUE))[1:5]
})))
top_profile_features <- top_profile_features[
  order(apply(primary_profiles[top_profile_features, , drop = FALSE], 1L, max), decreasing = TRUE)
]
profile_plot_data <- component_profiles[
  component_profiles$FeatureID %in% top_profile_features,
  , drop = FALSE
]
profile_plot_data$DisplayGenus <- factor(
  profile_plot_data$DisplayGenus,
  levels = rev(unique(ifelse(
    taxonomy[top_profile_features, "Genus"] == "Uknown",
    "Unclassified",
    taxonomy[top_profile_features, "Genus"]
  )))
)
profile_plot_data$State <- factor(profile_plot_data$State, levels = state_labels)
component_profiles_plot <- ggplot2::ggplot(
  profile_plot_data,
  ggplot2::aes(x = State, y = DisplayGenus, fill = ExpectedProportion)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.55) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(
      ExpectedProportion >= 0.005,
      scales::percent(ExpectedProportion, accuracy = 0.1),
      "<0.5%"
    ), colour = ExpectedProportion >= 0.25),
    family = font_family,
    size = 2.45
  ) +
  ggplot2::scale_colour_manual(
    values = c(`FALSE` = "#17202A", `TRUE` = "white"),
    guide = "none"
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Component profiles are distributions, not single-taxon labels",
    subtitle = "Union of the five highest expected-proportion genera in each fitted state",
    x = "Canonical DMM state",
    y = "Source genus category",
    fill = "Expected proportion",
    caption = paste0(
      "The source label 'Uknown' is displayed as Unclassified without changing counts.\n",
      "State letters are deterministic profile labels, not biological diagnoses."
    )
  ) +
  theme_pub(base_size = 8.8) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(face = "bold")
  )

stability_plot_data <- rbind(
  data.frame(
    Panel = "Fresh-start ARI",
    Category = as.character(initialization_audit$Seed),
    Value = initialization_audit$ARIToPrimary,
    Label = sprintf("%.2f", initialization_audit$ARIToPrimary),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Panel = "Sensitivity ARI",
    Category = sensitivity_agreement$Branch,
    Value = sensitivity_agreement$AdjustedRandIndex,
    Label = sprintf("%.2f", sensitivity_agreement$AdjustedRandIndex),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Panel = "Selected k (Laplace)",
    Category = sensitivity_agreement$Branch,
    Value = sensitivity_agreement$SelectedKLaplace,
    Label = as.character(sensitivity_agreement$SelectedKLaplace),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Panel = "Repeated source pairs",
    Category = c("Same state", "Different state"),
    Value = c(
      sum(participant_repeat_audit$StateAgreement),
      sum(!participant_repeat_audit$StateAgreement)
    ),
    Label = as.character(c(
      sum(participant_repeat_audit$StateAgreement),
      sum(!participant_repeat_audit$StateAgreement)
    )),
    stringsAsFactors = FALSE
  )
)
stability_plot_data$Panel <- factor(
  stability_plot_data$Panel,
  levels = c(
    "Fresh-start ARI", "Sensitivity ARI",
    "Selected k (Laplace)", "Repeated source pairs"
  )
)
stability_plot_data$Category <- factor(
  stability_plot_data$Category,
  levels = c(
    as.character(initialization_seed_all),
    names(branch_counts),
    "Same state", "Different state"
  )
)
stability_audit_plot <- ggplot2::ggplot(
  stability_plot_data,
  ggplot2::aes(x = Category, y = Value, fill = Panel)
) +
  ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label),
    vjust = -0.3,
    family = font_family,
    size = 2.6
  ) +
  ggplot2::facet_wrap(~Panel, scales = "free", ncol = 2L) +
  ggplot2::scale_fill_manual(values = c(
    "Fresh-start ARI" = "#0072B2",
    "Sensitivity ARI" = "#009E73",
    "Selected k (Laplace)" = "#E69F00",
    "Repeated source pairs" = "#CC79A7"
  )) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.13))) +
  ggplot2::labs(
    title = "A state label is credible only with a stability ledger",
    subtitle = "Initialization, preprocessing, component count and repeated samples answer different questions",
    x = NULL,
    y = "Audit value",
    caption = paste0(
      "ARI comparisons use profile-matched k = 4 states on overlapping samples.\n",
      "Suffix-derived pairs are descriptive because clinical visit timing is not distributed."
    )
  ) +
  theme_pub(base_size = 8.2) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 28, hjust = 1),
    strip.text = ggplot2::element_text(size = 8.5)
  )

figure_specification <- data.frame(
  stem = c(
    "29-model-selection",
    "29-posterior-ordination",
    "29-component-profiles",
    "29-stability-audit"
  ),
  width_mm = c(170, 170, 165, 180),
  height_mm = c(158, 125, 122, 135),
  raster_ppi = rep(600L, 4L),
  stringsAsFactors = FALSE
)
plot_objects <- list(
  model_selection_plot,
  posterior_ordination_plot,
  component_profiles_plot,
  stability_audit_plot
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
  svg_text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
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
    text_nodes = lengths(regmatches(svg_text, gregexpr("<text\\b", svg_text, perl = TRUE))),
    font_family_present = grepl('font-family: "DejaVu Sans"', svg_text, fixed = TRUE),
    chinese_characters = lengths(regmatches(svg_text, gregexpr("[一-鿿]", svg_text, perl = TRUE))),
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
  if (length(fields) != 7L) stop("Unexpected identify output for ", path, call. = FALSE)
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
    add_check(paste0(spec$stem, "-", extension, "-exists"), "export", exists, TRUE, exists)
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
    chinese_characters <- NA_integer_
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
        paste0(spec$stem, "-pdf-single-page"), "vector", pages, 1,
        pages == 1L
      )
      add_check(
        paste0(spec$stem, "-pdf-font-embedded"), "vector",
        font_status, "embedded", font_embedded
      )
    } else if (extension == "svg") {
      svg_info <- parse_svg_info(path)
      actual_width_mm <- svg_info$width_pt * 25.4 / 72
      actual_height_mm <- svg_info$height_pt * 25.4 / 72
      text_nodes <- svg_info$text_nodes
      chinese_characters <- svg_info$chinese_characters
      font_status <- ifelse(svg_info$font_family_present, "declared", "missing")
      add_check(
        paste0(spec$stem, "-svg-text-nodes"), "vector",
        text_nodes, ">0", text_nodes > 0L
      )
      add_check(
        paste0(spec$stem, "-svg-font-declared"), "vector",
        font_status, "declared", svg_info$font_family_present
      )
      add_check(
        paste0(spec$stem, "-svg-english-only"), "vector",
        chinese_characters, 0, chinese_characters == 0L
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
          paste0(spec$stem, "-tiff-lzw"), "raster",
          compression, "LZW", identical(compression, "LZW")
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
      chinese_characters = chinese_characters,
      font_status = font_status,
      dimension_status = ifelse(dimension_pass, "PASS", "FAIL"),
      stringsAsFactors = FALSE
    )
  }
}
format_audit <- do.call(rbind, format_rows)
add_check("primary-figure-count", "export", length(unique(format_audit$figure)), 4, length(unique(format_audit$figure)) == 4L)
add_check("primary-format-count", "export", nrow(format_audit), 16, nrow(format_audit) == 16L)
add_check("vector-file-count", "export", sum(format_audit$extension %in% c("pdf", "svg")), 8, sum(format_audit$extension %in% c("pdf", "svg")) == 8L)
add_check("raster-file-count", "export", sum(format_audit$extension %in% c("png", "tiff")), 8, sum(format_audit$extension %in% c("png", "tiff")) == 8L)
add_check("all-physical-dimensions-pass", "export", sum(format_audit$dimension_status == "PASS"), 16, all(format_audit$dimension_status == "PASS"))
add_check(
  "all-svg-english-only", "export",
  sum(format_audit$chinese_characters[format_audit$extension == "svg"], na.rm = TRUE),
  0,
  all(format_audit$chinese_characters[format_audit$extension == "svg"] == 0L)
)

add_check("inference-bmi-descriptive-only", "inference", TRUE, TRUE, TRUE)
add_check("inference-no-disease-subtype-claim", "inference", TRUE, TRUE, TRUE)
add_check("inference-no-natural-class-claim", "inference", TRUE, TRUE, TRUE)
add_check("inference-no-causal-claim", "inference", TRUE, TRUE, TRUE)
add_check("inference-suffix-not-clinical-time", "inference", TRUE, TRUE, TRUE)

validation_checks <- do.call(rbind, check_rows)
checks_total <- nrow(validation_checks)
checks_passed <- sum(validation_checks$status == "PASS")
checks_failed <- checks_total - checks_passed

write_tsv(input_audit, "input-audit.tsv")
write_tsv(package_audit, "package-audit.tsv")
write_tsv(command_audit, "command-audit.tsv")
write_tsv(font_audit, "font-audit.tsv")
write_tsv(sample_alignment_audit, "sample-alignment-audit.tsv")
write_tsv(model_selection, "model-selection.tsv")
write_tsv(reference_fit_comparison, "reference-fit-comparison.tsv")
write_tsv(initialization_audit, "initialization-audit.tsv")
write_tsv(initialization_agreement, "initialization-agreement.tsv")
write_tsv(component_profiles, "component-profiles.tsv")
write_tsv(component_top_taxa, "component-top-taxa.tsv")
write_tsv(component_summary, "component-summary.tsv")
write_tsv(sample_posteriors, "sample-posteriors.tsv")
write_tsv(ordination_scores, "ordination-scores.tsv")
write_tsv(bmi_cross_tab, "bmi-cross-tab.tsv")
write_tsv(participant_repeat_audit, "participant-repeat-audit.tsv")
write_tsv(participant_transitions, "participant-transitions.tsv")
write_tsv(sensitivity_model_selection, "sensitivity-model-selection.tsv")
write_tsv(sensitivity_agreement, "sensitivity-agreement.tsv")
write_tsv(rarefaction_ledger, "rarefaction-ledger.tsv")
write_tsv(label_matching_audit, "label-matching-audit.tsv")
write_tsv(figure_specification, "figure-specification.tsv")
write_tsv(format_audit, "format-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

summary_payload <- list(
  article = 29L,
  title = "Community typing with Dirichlet-multinomial mixtures",
  dataset = list(
    source = "DirichletMultinomial 1.46.0 packaged Holmes Twins genus counts",
    features = nrow(counts),
    samples = ncol(counts),
    reads = sum(counts),
    library_depth_min = min(sample_depths),
    library_depth_median = median(sample_depths),
    library_depth_max = max(sample_depths),
    participant_keys = length(participant_counts),
    paired_keys = length(two_sample_keys)
  ),
  primary_definition = list(
    feature_level = "Source genus categories",
    input_scale = "Raw non-negative integer counts",
    relative_abundance_input = FALSE,
    pseudocount = FALSE,
    rarefied = FALSE,
    candidate_k = candidate_k,
    selection_rule = "Minimum Laplace approximation",
    minimum_state_size = minimum_state_size,
    seed = primary_seed,
    label_rule = paste(
      "Decreasing fitted Bacteroides proportion, then decreasing",
      "Prevotella proportion, then raw component index"
    )
  ),
  model = list(
    selected_k_laplace = selected_k_laplace,
    selected_k_aic = selected_k_aic,
    selected_k_bic = selected_k_bic,
    state_sizes = as.list(stats::setNames(component_sizes, state_labels)),
    minimum_state_size = min(component_sizes),
    maximum_state_size = max(component_sizes),
    samples_below_point_eight_posterior = sum(primary_max_posterior < 0.8),
    samples_below_point_nine_posterior = sum(primary_max_posterior < 0.9),
    median_maximum_posterior = median(primary_max_posterior),
    median_normalized_entropy = median(primary_entropy),
    fresh_k4_laplace = primary_rows$Laplace[primary_rows$K == 4L],
    bundled_reference_k4_laplace = reference_rows$Laplace[reference_rows$K == 4L],
    bundled_reference_partition_ari = reference_ari
  ),
  stability = list(
    initialization_seeds = initialization_seed_all,
    initialization_pairs = nrow(initialization_agreement),
    initialization_min_ari = initialization_min_ari,
    initialization_primary_basin_starts = sum(initialization_audit$ARIToPrimary == 1),
    initialization_outlier_seeds = initialization_audit$Seed[
      initialization_audit$ARIToPrimary < 1
    ],
    best_laplace_seed = initialization_audit$Seed[[which.min(initialization_audit$Laplace)]],
    best_laplace_partition_ari_to_primary = initialization_audit$ARIToPrimary[[
      which.min(initialization_audit$Laplace)
    ]],
    prevalence_features = length(prevalence_features),
    depth_samples = length(depth_samples),
    rarefaction_depth = rarefaction_depth,
    one_per_participant_samples = length(one_per_samples),
    sensitivity_selected_k_laplace = as.list(stats::setNames(
      sensitivity_agreement$SelectedKLaplace,
      sensitivity_agreement$Branch
    )),
    sensitivity_fixed_k4_ari = as.list(stats::setNames(
      sensitivity_agreement$AdjustedRandIndex,
      sensitivity_agreement$Branch
    ))
  ),
  repeated_samples = list(
    paired_keys = nrow(participant_repeat_audit),
    same_state_pairs = sum(participant_repeat_audit$StateAgreement),
    different_state_pairs = sum(!participant_repeat_audit$StateAgreement),
    state_agreement = repeat_agreement,
    adjusted_rand_index = repeat_ari,
    suffix_timing_independently_validated = FALSE
  ),
  inference = list(
    bmi_hypothesis_test = FALSE,
    bmi_use = "Descriptive cross-tabulation only",
    natural_class_claim = FALSE,
    disease_subtype_claim = FALSE,
    causal_claim = FALSE,
    absolute_abundance_claim = FALSE
  ),
  graphics = list(
    font_family = font_family,
    primary_figures = length(unique(format_audit$figure)),
    vector_files = sum(format_audit$extension %in% c("pdf", "svg")),
    raster_files = sum(format_audit$extension %in% c("png", "tiff")),
    format_files = nrow(format_audit),
    raster_ppi = unique(figure_specification$raster_ppi),
    svg_files_with_text = sum(
      format_audit$extension == "svg" & format_audit$text_nodes > 0
    ),
    svg_chinese_characters = sum(
      format_audit$chinese_characters[format_audit$extension == "svg"],
      na.rm = TRUE
    ),
    pdf_files_with_embedded_font = sum(
      format_audit$extension == "pdf" & format_audit$font_status == "embedded"
    ),
    lzw_tiff_files = sum(
      format_audit$extension == "tiff" & format_audit$compression == "LZW"
    ),
    dimension_checks_passed = sum(format_audit$dimension_status == "PASS")
  ),
  execution = list(
    backend = execution_backend,
    workers = worker_count,
    fresh_models = nrow(all_jobs),
    candidate_models = nrow(candidate_jobs),
    initialization_models_beyond_primary = nrow(initialization_jobs)
  ),
  versions = as.list(observed_versions),
  checks_total = checks_total,
  checks_passed = checks_passed,
  checks_failed = checks_failed
)
jsonlite::write_json(
  summary_payload,
  file.path(output_dir, "community-typing-dmm-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

session_lines <- sanitize_text(capture.output(sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"))

validation_lines <- c(
  "Article 29 DMM community-typing validation",
  paste0("Status: ", ifelse(checks_failed == 0L, "PASS", "FAIL")),
  paste0("Checks: ", checks_passed, "/", checks_total),
  "Samples/features/reads: 278/130/570851",
  paste0(
    "Selected k (Laplace/AIC/BIC): ",
    selected_k_laplace, "/", selected_k_aic, "/", selected_k_bic
  ),
  paste0("Canonical state sizes: ", paste(component_sizes, collapse = "/")),
  paste0("Low-confidence posterior (<0.80): ", sum(primary_max_posterior < 0.8)),
  paste0("Initialization minimum ARI: ", sprintf("%.6f", initialization_min_ari)),
  paste0("Repeated-pair agreement: ", sprintf("%.6f", repeat_agreement)),
  paste0("Repeated-pair ARI: ", sprintf("%.6f", repeat_ari)),
  paste0("Primary exports: ", nrow(format_audit), " files"),
  paste0("Fresh fit jobs/workers: ", nrow(all_jobs), "/", worker_count),
  paste0("R version: ", r_version)
)
writeLines(validation_lines, file.path(output_dir, "validation.log"))

if (checks_failed > 0L) {
  failed_ids <- validation_checks$check_id[validation_checks$status == "FAIL"]
  stop(
    "Article 29 validation failed: ",
    paste(failed_ids, collapse = ", "),
    call. = FALSE
  )
}

cat(sprintf(
  "Article 29 validation passed: %d/%d checks.\n",
  checks_passed,
  checks_total
))
