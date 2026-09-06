# 群落分型 DMM（肠型/CST）
# Run sequentially in a new working directory.
# Required packages: DirichletMultinomial, ggplot2, knitr, ragg, readr, scales, svglite, vegan.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/community-typing-dmm/metadata.tsv",
  "data/small/community-typing-dmm/otutab.tsv",
  "data/small/community-typing-dmm/source-summary.json",
  "data/small/community-typing-dmm/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)

set.seed(20260729)
font_pub <- "sans"

pal_pub <- c(
  blue = "#0072B2",
  orange = "#E69F00",
  green = "#009E73",
  vermillion = "#D55E00",
  purple = "#CC79A7",
  sky = "#56B4E9",
  yellow = "#F0E442",
  grey = "#6B7280"
)

scale_color_pub <- function(..., values = unname(pal_pub)) {
  ggplot2::scale_color_manual(..., values = values)
}

scale_fill_pub <- function(..., values = unname(pal_pub)) {
  ggplot2::scale_fill_manual(..., values = values)
}

theme_pub <- function(base_size = 10, base_family = font_pub) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        colour = "#E6E6E6", linewidth = 0.25
      ),
      axis.text = ggplot2::element_text(colour = "#1A1A1A"),
      axis.title = ggplot2::element_text(colour = "#1A1A1A"),
      axis.ticks = ggplot2::element_line(
        colour = "#1A1A1A", linewidth = 0.3
      ),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.15)
      ),
      plot.subtitle = ggplot2::element_text(colour = "#4D4D4D"),
      plot.caption = ggplot2::element_text(
        colour = "#666666", hjust = 0
      ),
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2", colour = "#B3B3B3", linewidth = 0.3
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.position = "top"
    )
}

save_pub <- function(
  plot, file_base, width = 89, height = 70, units = "mm",
  dpi = 600, write_svg = TRUE, write_tiff = TRUE,
  base_family = font_pub
) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(file_base, ".pdf"), plot,
    width = width, height = height, units = units,
    device = grDevices::cairo_pdf,
    family = base_family, bg = "white"
  )
  if (isTRUE(write_svg)) {
    ggplot2::ggsave(
      paste0(file_base, ".svg"), plot,
      width = width, height = height, units = units,
      device = svglite::svglite, bg = "white"
    )
  }
  ggplot2::ggsave(
    paste0(file_base, ".png"), plot,
    width = width, height = height, units = units, dpi = dpi,
    device = ragg::agg_png, bg = "white"
  )
  if (isTRUE(write_tiff)) {
    ggplot2::ggsave(
      paste0(file_base, ".tiff"), plot,
      width = width, height = height, units = units, dpi = dpi,
      device = ragg::agg_tiff,
      compression = "lzw", bg = "white"
    )
  }
  invisible(plot)
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
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

input_paths <- c(
  otutab = "data/small/community-typing-dmm/otutab.tsv",
  taxonomy = "data/small/community-typing-dmm/taxonomy.tsv",
  metadata = "data/small/community-typing-dmm/metadata.tsv"
)
expected_sha256 <- c(
  otutab = "28de822c434aedade101e66ac60a61df38124c9b59357095599fbd278f2ecb41",
  taxonomy = "9eaa0bf788ad1db8c314a59e774c8935bdaa33a51135eaf606c11d5a497dff08",
  metadata = "0368c9129c7d77f381688465141593a2c61a229368d2fc411c9017e1abe04409"
)
observed_sha256 <- vapply(
  input_paths,
  function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ),
  character(1)
)
stopifnot(identical(observed_sha256, expected_sha256))
stopifnot(
  identical(
    as.character(packageVersion("DirichletMultinomial")),
    "1.46.0"
  )
)

otu_frame <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])

counts <- as.matrix(otu_frame)
storage.mode(counts) <- "integer"
taxonomy <- taxonomy[rownames(counts), , drop = FALSE]
metadata <- metadata[colnames(counts), , drop = FALSE]
metadata$SampleID <- rownames(metadata)
metadata$BMIClassCode <- as.integer(metadata$BMIClassCode)
metadata$SamplesPerParticipantKey <- as.integer(
  metadata$SamplesPerParticipantKey
)
sample_depths <- colSums(counts)
sample_counts <- t(counts)

stopifnot(
  identical(rownames(counts), rownames(taxonomy)),
  identical(colnames(counts), rownames(metadata)),
  nrow(counts) == 130L,
  ncol(counts) == 278L,
  sum(counts) == 570851L,
  all(is.finite(counts)),
  all(counts >= 0L),
  all(counts == round(counts)),
  all(rowSums(counts) > 0L),
  all(colSums(counts) > 0L),
  identical(range(sample_depths), c(53, 10585)),
  median(sample_depths) == 1597.5
)

list(
  otutab = data.frame(
    FeatureID = rownames(counts)[1:5],
    counts[1:5, 1:5, drop = FALSE],
    check.names = FALSE
  ),
  taxonomy = data.frame(
    FeatureID = rownames(taxonomy)[1:5],
    taxonomy[1:5, c("Genus", "TaxonomyStatus"), drop = FALSE],
    check.names = FALSE
  ),
  metadata = data.frame(
    SampleID = metadata$SampleID[1:5],
    metadata[1:5, c(
      "BMIClass", "ParticipantKey", "SourceSuffix",
      "SamplesPerParticipantKey"
    )],
    check.names = FALSE
  )
)

primary_seed <- 20260729L
candidate_k <- 1:7
minimum_state_size <- max(10L, ceiling(0.05 * nrow(sample_counts)))
stopifnot(minimum_state_size == 14L)
state_labels <- paste0("DMM-", LETTERS[1:4])

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
  x <- DirichletMultinomial::fitted(fit_object, scale = TRUE)
  sweep(x, 2L, colSums(x), "/")
}

fit_posteriors <- function(fit_object, sample_ids) {
  x <- DirichletMultinomial::mixture(fit_object)
  x <- sweep(x, 1L, rowSums(x), "/")
  rownames(x) <- sample_ids
  x
}

choose2 <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  expected <- sum_rows * sum_cols / choose2(n)
  maximum <- (sum_rows + sum_cols) / 2
  (sum_cells - expected) / (maximum - expected)
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[[i]], all_permutations(x[-i]))
  }))
}

js_divergence <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  m <- (p + q) / 2
  sum(ifelse(p > 0, p * log(p / m), 0)) / 2 +
    sum(ifelse(q > 0, q * log(q / m), 0)) / 2
}

match_profiles <- function(primary_profiles, candidate_profiles) {
  common <- intersect(
    rownames(primary_profiles),
    rownames(candidate_profiles)
  )
  primary <- primary_profiles[common, , drop = FALSE]
  candidate <- candidate_profiles[common, , drop = FALSE]
  primary <- sweep(primary, 2L, colSums(primary), "/")
  candidate <- sweep(candidate, 2L, colSums(candidate), "/")
  k <- ncol(primary)
  cost <- outer(seq_len(k), seq_len(k), Vectorize(function(i, j) {
    js_divergence(primary[, i], candidate[, j])
  }))
  permutations <- all_permutations(seq_len(k))
  total <- apply(permutations, 1L, function(p) {
    sum(cost[cbind(seq_len(k), p)])
  })
  best <- permutations[which.min(total), ]
  list(
    candidate_for_primary = as.integer(best),
    divergence = cost[cbind(seq_len(k), best)],
    common_features = length(common)
  )
}

feature_prevalence <- rowMeans(counts > 0L)
prevalence_features <- rownames(counts)[feature_prevalence >= 0.05]
depth_samples <- names(sample_depths)[sample_depths >= 1000L]

set.seed(primary_seed)
rarefied_counts <- vegan::rrarefy(
  sample_counts[depth_samples, , drop = FALSE],
  sample = 1000L
)
set.seed(primary_seed)
stopifnot(identical(
  rarefied_counts,
  vegan::rrarefy(
    sample_counts[depth_samples, , drop = FALSE],
    sample = 1000L
  )
))

participant_indices <- split(
  seq_len(nrow(metadata)),
  metadata$ParticipantKey
)
one_per_participant <- vapply(participant_indices, function(idx) {
  no_suffix <- idx[metadata$SourceSuffix[idx] == "No suffix"]
  if (length(no_suffix) > 0L) no_suffix[[1L]] else idx[[1L]]
}, integer(1))
one_per_samples <- metadata$SampleID[one_per_participant]

branch_counts <- list(
  Primary = sample_counts,
  `Prevalence >= 5%` = sample_counts[
    , prevalence_features, drop = FALSE
  ],
  `Depth >= 1,000` = sample_counts[
    depth_samples, , drop = FALSE
  ],
  `Rarefied to 1,000` = rarefied_counts,
  `One per participant` = sample_counts[
    one_per_samples, , drop = FALSE
  ]
)
branch_seeds <- c(
  Primary = primary_seed,
  `Prevalence >= 5%` = primary_seed + 100L,
  `Depth >= 1,000` = primary_seed + 200L,
  `Rarefied to 1,000` = primary_seed + 300L,
  `One per participant` = primary_seed + 400L
)

data.frame(
  Branch = names(branch_counts),
  Samples = vapply(branch_counts, nrow, integer(1)),
  Features = vapply(branch_counts, ncol, integer(1)),
  Seed = unname(branch_seeds[names(branch_counts)]),
  check.names = FALSE
)

candidate_jobs <- expand.grid(
  Branch = names(branch_counts),
  K = candidate_k,
  stringsAsFactors = FALSE
)
candidate_jobs$Seed <- unname(branch_seeds[candidate_jobs$Branch])

initialization_jobs <- data.frame(
  Branch = "Initialization",
  K = 4L,
  Seed = primary_seed + 1:7,
  stringsAsFactors = FALSE
)
all_jobs <- rbind(candidate_jobs, initialization_jobs)

job_payload <- lapply(seq_len(nrow(all_jobs)), function(i) {
  job <- all_jobs[i, , drop = FALSE]
  list(
    counts = if (job$Branch == "Initialization") {
      branch_counts[["Primary"]]
    } else {
      branch_counts[[job$Branch]]
    },
    k = job$K,
    seed = job$Seed
  )
})

fit_one <- function(job) {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1"
  )
  requireNamespace("DirichletMultinomial", quietly = TRUE)
  DirichletMultinomial::dmn(
    job$counts,
    k = job$k,
    verbose = FALSE,
    seed = job$seed
  )
}

available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores)) available_cores <- 1L
workers <- min(8L, max(1L, available_cores), length(job_payload))

job_indices <- seq_along(job_payload)
job_batches <- split(
  job_indices,
  ceiling(job_indices / workers)
)
all_fits <- vector("list", length(job_payload))

for (batch_indices in job_batches) {
  cl <- parallel::makePSOCKcluster(length(batch_indices))
  invisible(parallel::clusterEvalQ(cl, {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1"
    )
    requireNamespace("DirichletMultinomial", quietly = TRUE)
    NULL
  }))
  batch_fits <- tryCatch(
    parallel::parLapply(
      cl,
      job_payload[batch_indices],
      fit_one
    ),
    finally = parallel::stopCluster(cl)
  )
  all_fits[batch_indices] <- batch_fits
}

candidate_fits <- all_fits[seq_len(nrow(candidate_jobs))]
initialization_extra_fits <- all_fits[
  nrow(candidate_jobs) + seq_len(nrow(initialization_jobs))
]

get_candidate_fit <- function(branch, k) {
  index <- which(
    candidate_jobs$Branch == branch & candidate_jobs$K == k
  )
  stopifnot(length(index) == 1L)
  candidate_fits[[index]]
}

model_selection <- do.call(rbind, lapply(
  seq_len(nrow(candidate_jobs)),
  function(i) {
    score <- fit_scores(candidate_fits[[i]])
    data.frame(
      Branch = candidate_jobs$Branch[[i]],
      K = candidate_jobs$K[[i]],
      Seed = candidate_jobs$Seed[[i]],
      NLE = score[["NLE"]],
      LogDet = score[["LogDet"]],
      Laplace = score[["Laplace"]],
      AIC = score[["AIC"]],
      BIC = score[["BIC"]],
      stringsAsFactors = FALSE
    )
  }
))

selected_for <- function(branch, criterion) {
  x <- model_selection[model_selection$Branch == branch, , drop = FALSE]
  x$K[[which.min(x[[criterion]])]]
}

primary_selection <- model_selection[
  model_selection$Branch == "Primary",
  c("K", "Laplace", "AIC", "BIC")
]
primary_selection$LaplaceMinimum <-
  primary_selection$Laplace == min(primary_selection$Laplace)
primary_selection$AICMinimum <-
  primary_selection$AIC == min(primary_selection$AIC)
primary_selection$BICMinimum <-
  primary_selection$BIC == min(primary_selection$BIC)

knitr::kable(
  primary_selection,
  digits = 2,
  caption = "Fresh fixed-seed model-selection ledger"
)

primary_fit <- get_candidate_fit("Primary", 4L)
profiles_raw <- fit_profiles(primary_fit)
posteriors_raw <- fit_posteriors(primary_fit, rownames(sample_counts))

canonical_order <- order(
  -profiles_raw["Bacteroides", ],
  -profiles_raw["Prevotella", ],
  seq_len(ncol(profiles_raw))
)
profiles <- profiles_raw[, canonical_order, drop = FALSE]
posteriors <- posteriors_raw[, canonical_order, drop = FALSE]
colnames(profiles) <- state_labels
colnames(posteriors) <- state_labels

hard_index <- max.col(posteriors, ties.method = "first")
primary_state <- factor(state_labels[hard_index], levels = state_labels)
names(primary_state) <- rownames(posteriors)
maximum_posterior <- apply(posteriors, 1L, max)
normalized_entropy <- -rowSums(
  posteriors * log(pmax(posteriors, 1e-15))
) / log(ncol(posteriors))

mixture_weights <- DirichletMultinomial::mixturewt(primary_fit)[
  canonical_order, , drop = FALSE
]
component_sizes <- as.integer(table(primary_state))
component_summary <- data.frame(
  State = state_labels,
  Samples = component_sizes,
  MixtureWeight = mixture_weights$pi,
  Theta = mixture_weights$theta,
  ExpectedBacteroides = profiles["Bacteroides", ],
  ExpectedPrevotella = profiles["Prevotella", ],
  SamplesBelowPointEight = vapply(state_labels, function(state) {
    sum(primary_state == state & maximum_posterior < 0.8)
  }, integer(1)),
  stringsAsFactors = FALSE
)

stopifnot(
  identical(component_sizes, c(85L, 86L, 46L, 61L)),
  min(component_sizes) >= minimum_state_size,
  max(abs(colSums(profiles) - 1)) <= 1e-12,
  max(abs(rowSums(posteriors) - 1)) <= 1e-12,
  sum(maximum_posterior < 0.8) == 8L,
  sum(maximum_posterior < 0.9) == 17L
)

knitr::kable(
  component_summary,
  digits = 4,
  caption = "Canonical component summary"
)

sample_posteriors <- data.frame(
  SampleID = rownames(posteriors),
  DMM_A = posteriors[, "DMM-A"],
  DMM_B = posteriors[, "DMM-B"],
  DMM_C = posteriors[, "DMM-C"],
  DMM_D = posteriors[, "DMM-D"],
  State = as.character(primary_state),
  MaximumPosterior = maximum_posterior,
  NormalizedEntropy = normalized_entropy,
  LowConfidencePointEight = maximum_posterior < 0.8,
  ParticipantKey = metadata[rownames(posteriors), "ParticipantKey"],
  SourceSuffix = metadata[rownames(posteriors), "SourceSuffix"],
  BMIClass = metadata[rownames(posteriors), "BMIClass"],
  LibrarySize = sample_depths[rownames(posteriors)],
  stringsAsFactors = FALSE
)
head(sample_posteriors, 8)

reference_environment <- new.env(parent = emptyenv())
data(
  "fit",
  package = "DirichletMultinomial",
  envir = reference_environment
)
reference_fits <- reference_environment$fit

reference_fit_comparison <- do.call(rbind, lapply(
  c("Fresh fixed-seed fit", "Bundled reference"),
  function(source_name) {
    fit_list <- if (source_name == "Fresh fixed-seed fit") {
      lapply(candidate_k, function(k) get_candidate_fit("Primary", k))
    } else {
      reference_fits[candidate_k]
    }
    do.call(rbind, lapply(candidate_k, function(k) {
      score <- fit_scores(fit_list[[k]])
      data.frame(
        Source = source_name,
        K = k,
        Laplace = score[["Laplace"]],
        AIC = score[["AIC"]],
        BIC = score[["BIC"]],
        stringsAsFactors = FALSE
      )
    }))
  }
))

reference_profiles_raw <- fit_profiles(reference_fits[[4L]])
reference_match <- match_profiles(profiles, reference_profiles_raw)
reference_posteriors <- fit_posteriors(
  reference_fits[[4L]], rownames(sample_counts)
)[, reference_match$candidate_for_primary, drop = FALSE]
reference_state <- state_labels[
  max.col(reference_posteriors, ties.method = "first")
]
reference_ari <- adjusted_rand_index(
  as.character(primary_state), reference_state
)

stopifnot(
  reference_fit_comparison$K[
    reference_fit_comparison$Source == "Bundled reference"
  ][which.min(reference_fit_comparison$Laplace[
    reference_fit_comparison$Source == "Bundled reference"
  ])] == 4L,
  abs(reference_ari - 1) <= 1e-12
)

initialization_fits <- c(
  list(primary_fit),
  initialization_extra_fits
)
initialization_seeds <- primary_seed + 0:7

initialization_audit <- do.call(rbind, lapply(
  seq_along(initialization_fits),
  function(i) {
    fit_object <- initialization_fits[[i]]
    profile_match <- match_profiles(profiles, fit_profiles(fit_object))
    candidate_posterior <- fit_posteriors(
      fit_object, rownames(sample_counts)
    )[, profile_match$candidate_for_primary, drop = FALSE]
    candidate_state <- state_labels[
      max.col(candidate_posterior, ties.method = "first")
    ]
    score <- fit_scores(fit_object)
    data.frame(
      Seed = initialization_seeds[[i]],
      Laplace = score[["Laplace"]],
      ARIToPrimary = adjusted_rand_index(
        as.character(primary_state), candidate_state
      ),
      StateSizes = paste(
        as.integer(table(factor(candidate_state, levels = state_labels))),
        collapse = "/"
      ),
      MeanMatchedJSD = mean(profile_match$divergence),
      stringsAsFactors = FALSE
    )
  }
))

sensitivity_agreement <- do.call(rbind, lapply(
  names(branch_counts),
  function(branch) {
    fit_object <- get_candidate_fit(branch, 4L)
    profile_match <- match_profiles(profiles, fit_profiles(fit_object))
    candidate_posterior <- fit_posteriors(
      fit_object, rownames(branch_counts[[branch]])
    )[, profile_match$candidate_for_primary, drop = FALSE]
    candidate_state <- state_labels[
      max.col(candidate_posterior, ties.method = "first")
    ]
    names(candidate_state) <- rownames(candidate_posterior)
    overlap <- intersect(names(primary_state), names(candidate_state))
    data.frame(
      Branch = branch,
      Samples = nrow(branch_counts[[branch]]),
      Features = ncol(branch_counts[[branch]]),
      SelectedKLaplace = selected_for(branch, "Laplace"),
      FixedK4ARI = adjusted_rand_index(
        as.character(primary_state[overlap]),
        candidate_state[overlap]
      ),
      MeanMatchedJSD = mean(profile_match$divergence),
      stringsAsFactors = FALSE
    )
  }
))

stopifnot(
  sum(initialization_audit$ARIToPrimary == 1) == 8L,
  length(initialization_audit$Seed[
    initialization_audit$ARIToPrimary < 1
  ]) == 0L,
  initialization_audit$ARIToPrimary[
    which.min(initialization_audit$Laplace)
  ] == 1,
  identical(
    sensitivity_agreement$SelectedKLaplace,
    c(4L, 4L, 2L, 2L, 2L)
  )
)

knitr::kable(
  initialization_audit,
  digits = 4,
  caption = "Eight one-shot independent k=4 starts"
)
knitr::kable(
  sensitivity_agreement,
  digits = 4,
  caption = "Preprocessing and sampling sensitivity"
)

participant_counts <- table(metadata$ParticipantKey)
paired_keys <- names(participant_counts)[participant_counts == 2L]

participant_repeat_audit <- do.call(rbind, lapply(
  paired_keys,
  function(key) {
    ids <- metadata$SampleID[metadata$ParticipantKey == key]
    no_suffix_id <- ids[metadata[ids, "SourceSuffix"] == "No suffix"]
    dot2_id <- ids[metadata[ids, "SourceSuffix"] == "Dot2 suffix"]
    data.frame(
      ParticipantKey = key,
      NoSuffixSample = no_suffix_id,
      Dot2SuffixSample = dot2_id,
      NoSuffixState = as.character(primary_state[no_suffix_id]),
      Dot2SuffixState = as.character(primary_state[dot2_id]),
      StateAgreement = primary_state[no_suffix_id] ==
        primary_state[dot2_id],
      stringsAsFactors = FALSE
    )
  }
))

repeat_agreement <- mean(participant_repeat_audit$StateAgreement)
repeat_ari <- adjusted_rand_index(
  participant_repeat_audit$NoSuffixState,
  participant_repeat_audit$Dot2SuffixState
)

bmi_cross_tab <- as.data.frame(table(
  State = factor(primary_state, levels = state_labels),
  BMIClass = factor(
    metadata[names(primary_state), "BMIClass"],
    levels = c("Lean", "Obese", "Overweight")
  )
), stringsAsFactors = FALSE)
names(bmi_cross_tab)[3L] <- "Samples"
bmi_cross_tab$WithinStateFraction <- ave(
  bmi_cross_tab$Samples,
  bmi_cross_tab$State,
  FUN = function(x) x / sum(x)
)

stopifnot(
  length(participant_counts) == 154L,
  length(paired_keys) == 124L,
  sum(participant_repeat_audit$StateAgreement) == 85L,
  abs(repeat_agreement - 0.685483870967742) <= 1e-12,
  abs(repeat_ari - 0.37678348271377) <= 1e-12
)

knitr::kable(
  bmi_cross_tab,
  digits = 3,
  caption = "Descriptive BMIClass cross-tabulation; no hypothesis test"
)

relative_samples <- sweep(
  sample_counts, 1L, rowSums(sample_counts), "/"
)
bray <- vegan::vegdist(relative_samples, method = "bray")
pcoa <- cmdscale(bray, k = 2L, eig = TRUE, add = TRUE)
pcoa_points <- pcoa$points

if (cor(pcoa_points[, 1L], relative_samples[, "Bacteroides"]) < 0) {
  pcoa_points[, 1L] <- -pcoa_points[, 1L]
}
if (cor(pcoa_points[, 2L], relative_samples[, "Prevotella"]) < 0) {
  pcoa_points[, 2L] <- -pcoa_points[, 2L]
}
positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
axis_percent <- 100 * pcoa$eig[1:2] / sum(positive_eigenvalues)

ordination_scores <- data.frame(
  SampleID = rownames(pcoa_points),
  PCoA1 = pcoa_points[, 1L],
  PCoA2 = pcoa_points[, 2L],
  State = factor(
    primary_state[rownames(pcoa_points)],
    levels = state_labels
  ),
  MaximumPosterior = maximum_posterior[rownames(pcoa_points)],
  ConfidenceClass = ifelse(
    maximum_posterior[rownames(pcoa_points)] >= 0.8,
    "Posterior >= 0.80",
    "Posterior < 0.80"
  ),
  stringsAsFactors = FALSE
)

audit_upgrade <- data.frame(
  Dimension = c(
    "Likelihood input", "Component count", "Ordination", "Assignments",
    "Initialization", "Component identity", "Uncertainty", "Phenotype"
  ),
  OriginalAnalysis = c(
    "Raw sample-by-genus frequencies",
    "Laplace negative log model evidence",
    "Bray-Curtis NMDS with component means",
    "Posterior DMM component",
    "One reported fitted solution",
    "Numeric component index",
    "Probabilistic model, hard labels shown",
    "BMI-by-component comparison"
  ),
  ReanalysisAndUpgrade = c(
    "Same raw integer count contract",
    "Laplace primary; AIC and BIC shown in parallel",
    "Bray-Curtis PCoA used only as a display",
    "All posterior probabilities retained",
    "Eight one-shot independent fresh starts",
    "Profile matching plus deterministic names",
    "Maximum posterior and entropy audited",
    "Descriptive only under incomplete covariates"
  ),
  stringsAsFactors = FALSE
)

knitr::kable(
  audit_upgrade,
  caption = "Original DMM analysis and prespecified upgrades"
)

criterion_plot_data <- do.call(rbind, lapply(
  c("Laplace", "AIC", "BIC"),
  function(criterion) {
    data.frame(
      Source = reference_fit_comparison$Source,
      K = reference_fit_comparison$K,
      Criterion = criterion,
      Score = reference_fit_comparison[[criterion]],
      stringsAsFactors = FALSE
    )
  }
))
criterion_plot_data$Minimum <- ave(
  criterion_plot_data$Score,
  interaction(
    criterion_plot_data$Source,
    criterion_plot_data$Criterion
  ),
  FUN = function(x) x == min(x)
) == 1
criterion_plot_data$Criterion <- factor(
  criterion_plot_data$Criterion,
  levels = c("Laplace", "AIC", "BIC")
)

model_selection_plot <- ggplot(
  criterion_plot_data,
  aes(x = K, y = Score, colour = Source, linetype = Source)
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.0) +
  geom_point(
    data = criterion_plot_data[criterion_plot_data$Minimum, ],
    shape = 21, fill = "white", stroke = 1.0, size = 3.4
  ) +
  facet_wrap(~Criterion, scales = "free_y", ncol = 1L) +
  scale_x_continuous(breaks = candidate_k) +
  scale_colour_manual(values = c(
    "Fresh fixed-seed fit" = "#0072B2",
    "Bundled reference" = "#D55E00"
  )) +
  scale_linetype_manual(values = c(
    "Fresh fixed-seed fit" = "solid",
    "Bundled reference" = "22"
  )) +
  scale_y_continuous(labels = scales::label_comma(accuracy = 1)) +
  labs(
    title = "Model-selection criteria do not tell the same story",
    subtitle = paste(
      "Open circles mark criterion minima within each",
      "independently fitted series"
    ),
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
  theme(legend.position = "bottom")

save_pub(
  model_selection_plot,
  "figures/29-model-selection",
  width = 170,
  height = 158
)
model_selection_plot

state_colours <- c(
  "DMM-A" = "#0072B2",
  "DMM-B" = "#E69F00",
  "DMM-C" = "#009E73",
  "DMM-D" = "#CC79A7"
)
state_centroids <- aggregate(
  ordination_scores[, c("PCoA1", "PCoA2")],
  list(State = ordination_scores$State),
  mean
)

posterior_ordination_plot <- ggplot(
  ordination_scores,
  aes(x = PCoA1, y = PCoA2, colour = State)
) +
  geom_point(
    aes(alpha = MaximumPosterior, shape = ConfidenceClass),
    size = 2.2,
    stroke = 0.55
  ) +
  geom_point(
    data = state_centroids,
    shape = 4,
    size = 4.2,
    stroke = 1.1,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = state_colours, drop = FALSE) +
  scale_alpha_continuous(
    range = c(0.32, 0.9), limits = c(0, 1), guide = "none"
  ) +
  scale_shape_manual(values = c(
    "Posterior >= 0.80" = 16,
    "Posterior < 0.80" = 1
  )) +
  coord_equal() +
  labs(
    title = "DMM states are posterior assignments, not PCoA quadrants",
    subtitle = paste(
      "Bray-Curtis PCoA is shown only as a two-dimensional",
      "descriptive view"
    ),
    x = sprintf(
      "PCoA 1 (%.1f%% of positive eigenvalue sum)",
      axis_percent[[1L]]
    ),
    y = sprintf(
      "PCoA 2 (%.1f%% of positive eigenvalue sum)",
      axis_percent[[2L]]
    ),
    colour = "DMM state",
    shape = "Assignment confidence",
    caption = paste0(
      "Crosses are state centroids. Open points have maximum posterior below 0.80.\n",
      "The DMM was fitted to raw counts; ordination did not define the states."
    )
  ) +
  theme_pub(base_size = 8.8) +
  guides(
    colour = guide_legend(order = 1L, nrow = 1L),
    shape = guide_legend(order = 2L, nrow = 1L)
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

save_pub(
  posterior_ordination_plot,
  "figures/29-posterior-ordination",
  width = 170,
  height = 125
)
posterior_ordination_plot

top_profile_features <- unique(unlist(lapply(state_labels, function(state) {
  names(sort(profiles[, state], decreasing = TRUE))[1:5]
})))
top_profile_features <- top_profile_features[
  order(
    apply(profiles[top_profile_features, , drop = FALSE], 1L, max),
    decreasing = TRUE
  )
]

profile_plot_data <- data.frame(
  FeatureID = rep(rownames(profiles), times = ncol(profiles)),
  DisplayGenus = rep(ifelse(
    taxonomy[rownames(profiles), "Genus"] == "Uknown",
    "Unclassified",
    taxonomy[rownames(profiles), "Genus"]
  ), times = ncol(profiles)),
  State = rep(state_labels, each = nrow(profiles)),
  ExpectedProportion = as.vector(profiles),
  stringsAsFactors = FALSE
)
profile_plot_data <- profile_plot_data[
  profile_plot_data$FeatureID %in% top_profile_features,
]
profile_plot_data$DisplayGenus <- factor(
  profile_plot_data$DisplayGenus,
  levels = rev(ifelse(
    taxonomy[top_profile_features, "Genus"] == "Uknown",
    "Unclassified",
    taxonomy[top_profile_features, "Genus"]
  ))
)
profile_plot_data$State <- factor(
  profile_plot_data$State,
  levels = state_labels
)

component_profiles_plot <- ggplot(
  profile_plot_data,
  aes(x = State, y = DisplayGenus, fill = ExpectedProportion)
) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(
    aes(
      label = ifelse(
        ExpectedProportion >= 0.005,
        scales::percent(ExpectedProportion, accuracy = 0.1),
        "<0.5%"
      ),
      colour = ExpectedProportion >= 0.25
    ),
    family = font_pub,
    size = 2.45
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "#17202A", `TRUE` = "white"),
    guide = "none"
  ) +
  scale_fill_gradientn(
    colours = c(
      "#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"
    ),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    title = "Component profiles are distributions, not single-taxon labels",
    subtitle = paste(
      "Union of the five highest expected-proportion genera",
      "in each fitted state"
    ),
    x = "Canonical DMM state",
    y = "Source genus category",
    fill = "Expected proportion",
    caption = paste0(
      "The source label 'Uknown' is displayed as Unclassified without changing counts.\n",
      "State letters are deterministic profile labels, not biological diagnoses."
    )
  ) +
  theme_pub(base_size = 8.8) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold")
  )

save_pub(
  component_profiles_plot,
  "figures/29-component-profiles",
  width = 165,
  height = 122
)
component_profiles_plot

stability_plot_data <- rbind(
  data.frame(
    Panel = "Fresh-start ARI",
    Category = as.character(initialization_audit$Seed),
    Value = initialization_audit$ARIToPrimary,
    Label = sprintf("%.2f", initialization_audit$ARIToPrimary)
  ),
  data.frame(
    Panel = "Sensitivity ARI",
    Category = sensitivity_agreement$Branch,
    Value = sensitivity_agreement$FixedK4ARI,
    Label = sprintf("%.2f", sensitivity_agreement$FixedK4ARI)
  ),
  data.frame(
    Panel = "Selected k (Laplace)",
    Category = sensitivity_agreement$Branch,
    Value = sensitivity_agreement$SelectedKLaplace,
    Label = as.character(sensitivity_agreement$SelectedKLaplace)
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
    ))
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
    as.character(initialization_seeds),
    names(branch_counts),
    "Same state", "Different state"
  )
)

stability_audit_plot <- ggplot(
  stability_plot_data,
  aes(x = Category, y = Value, fill = Panel)
) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(
    aes(label = Label),
    vjust = -0.3,
    family = font_pub,
    size = 2.6
  ) +
  facet_wrap(~Panel, scales = "free", ncol = 2L) +
  scale_fill_manual(values = c(
    "Fresh-start ARI" = "#0072B2",
    "Sensitivity ARI" = "#009E73",
    "Selected k (Laplace)" = "#E69F00",
    "Repeated source pairs" = "#CC79A7"
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.13))
  ) +
  labs(
    title = "A state label is credible only with a stability ledger",
    subtitle = paste(
      "Initialization, preprocessing, component count and repeated",
      "samples answer different questions"
    ),
    x = NULL,
    y = "Audit value",
    caption = paste0(
      "ARI comparisons use profile-matched k = 4 states on overlapping samples.\n",
      "Suffix-derived pairs are descriptive because clinical visit timing is not distributed."
    )
  ) +
  theme_pub(base_size = 8.2) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 28, hjust = 1),
    strip.text = element_text(size = 8.5)
  )

save_pub(
  stability_audit_plot,
  "figures/29-stability-audit",
  width = 180,
  height = 135
)
stability_audit_plot

stopifnot(
  selected_for("Primary", "Laplace") == 4L,
  selected_for("Primary", "AIC") == 2L,
  selected_for("Primary", "BIC") == 1L,
  identical(component_sizes, c(85L, 86L, 46L, 61L)),
  sum(maximum_posterior < 0.8) == 8L,
  sum(initialization_audit$ARIToPrimary == 1) == 8L,
  identical(
    sensitivity_agreement$SelectedKLaplace,
    c(4L, 4L, 2L, 2L, 2L)
  ),
  sum(participant_repeat_audit$StateAgreement) == 85L
)

audit_path <- file.path(
  "results/29-community-typing-dmm",
  "validation-checks.tsv"
)
if (file.exists(audit_path)) {
  validation_checks <- readr::read_tsv(
    audit_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  stopifnot(
    nrow(validation_checks) == 222L,
    all(validation_checks$status == "PASS")
  )
  data.frame(
    Validation = "Repository QA",
    Passed = sum(validation_checks$status == "PASS"),
    Total = nrow(validation_checks)
  )
} else {
  data.frame(
    Validation = "Inline reproducibility checks",
    Passed = 8L,
    Total = 8L
  )
}
