# 核心菌群（core microbiome）+ 稀有生物圈
# Run sequentially in a new working directory.
# Required packages: ggplot2, ggrepel, ragg, readr, scales, svglite, vegan.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/metadata.tsv",
  "data/small/otutab.tsv",
  "data/small/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)

set.seed(20260728)
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
      plot.caption = ggplot2::element_text(colour = "#666666", hjust = 0),
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
    device = grDevices::cairo_pdf, family = base_family, bg = "white"
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
      device = ragg::agg_tiff, compression = "lzw", bg = "white"
    )
  }
  invisible(plot)
}

clean_taxon <- function(x, rank_name) {
  x <- trimws(as.character(x))
  x <- sub("^[A-Za-z]__", "", x)
  missing <- is.na(x) | x == "" | tolower(x) %in% c(
    "unassigned", "unclassified", "unknown", "na", "nan"
  )
  x[missing] <- paste("Unknown", tolower(rank_name))
  x
}

jaccard <- function(a, b) {
  union_n <- sum(a | b)
  if (union_n == 0L) 1 else sum(a & b) / union_n
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
  otutab = "data/small/otutab.tsv",
  taxonomy = "data/small/taxonomy.tsv",
  metadata = "data/small/metadata.tsv"
)
expected_sha256 <- c(
  otutab = "76fa79c38da889f35978dc86da4641a270961746709ff38049ee5f67e3c6f7a3",
  taxonomy = "725280bb9a0cd9bda7b540022e92af945ceed52527f8b2d220b055b4489f6901",
  metadata = "df24771dccf27607ddf922c6bca2cafa876d946fbe2e09d14b601accce66ba64"
)
observed_sha256 <- vapply(
  input_paths,
  function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ),
  character(1)
)
stopifnot(identical(observed_sha256, expected_sha256))

otutab <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])

counts <- as.matrix(otutab)
storage.mode(counts) <- "integer"
stopifnot(
  identical(rownames(counts), rownames(taxonomy)),
  setequal(colnames(counts), rownames(metadata)),
  all(counts >= 0L),
  all(counts == round(counts)),
  all(rowSums(counts) > 0L),
  all(colSums(counts) > 0L)
)
metadata <- metadata[colnames(counts), , drop = FALSE]

group_codes <- c("IW", "CW", "TW")
group_display <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan wetland"
)
metadata$SampleID <- rownames(metadata)
metadata$GroupDisplay <- unname(group_display[metadata$Group])

sample_depths <- colSums(counts)
relative_feature <- sweep(counts, 2L, sample_depths, "/")
stopifnot(max(abs(colSums(relative_feature) - 1)) <= 1e-12)

data.frame(
  Features = nrow(counts),
  Samples = ncol(counts),
  Reads = sum(counts),
  MinimumDepth = min(sample_depths),
  MaximumDepth = max(sample_depths)
)

otutab[1:5, 1:5, drop = FALSE]
taxonomy[1:6, , drop = FALSE]
metadata[1:6, , drop = FALSE]

primary_detection <- 0.0001
primary_prevalence <- 0.80
primary_rare_threshold <- 0.001
primary_seed <- 20260728L

detected_primary <- relative_feature >= primary_detection
detected_samples <- rowSums(detected_primary)
global_prevalence <- detected_samples / ncol(counts)
global_core <- global_prevalence >= primary_prevalence
feature_equal_mean <- rowMeans(relative_feature)
feature_pooled <- rowSums(counts) / sum(counts)
feature_maximum <- apply(relative_feature, 1L, max)
positive_samples <- rowSums(counts > 0L)

data.frame(
  DetectionThreshold = primary_detection,
  PrevalenceThreshold = primary_prevalence,
  CoreFeatures = sum(global_core),
  CoreAbundanceMass = sum(feature_equal_mean[global_core]),
  CoreBelowPointOnePercentMean = sum(
    global_core & feature_equal_mean < primary_rare_threshold
  )
)

rank_names <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
)
rank_prefix <- c(
  Kingdom = "K", Phylum = "P", Class = "C", Order = "O",
  Family = "F", Genus = "G", Species = "S"
)
taxonomy_clean <- as.data.frame(
  setNames(
    lapply(
      rank_names,
      function(rank_name) clean_taxon(taxonomy[[rank_name]], rank_name)
    ),
    rank_names
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(taxonomy_clean) <- rownames(taxonomy)

known_matrix <- vapply(
  rank_names,
  function(rank_name) {
    taxonomy_clean[[rank_name]] != paste("Unknown", tolower(rank_name))
  },
  logical(nrow(taxonomy_clean))
)
best_index <- apply(known_matrix, 1L, function(x) {
  known <- which(x)
  if (length(known) == 0L) 1L else max(known)
})
best_rank <- rank_names[best_index]
best_taxon <- vapply(
  seq_len(nrow(taxonomy_clean)),
  function(i) taxonomy_clean[[best_rank[[i]]]][[i]],
  character(1)
)
display_label <- paste0(unname(rank_prefix[best_rank]), " · ", best_taxon)
lineage <- apply(taxonomy_clean, 1L, paste, collapse = " > ")

group_core <- list()
group_prevalence <- list()
for (group_code in group_codes) {
  ids <- metadata$SampleID[metadata$Group == group_code]
  group_prevalence[[group_code]] <- rowMeans(
    detected_primary[, ids, drop = FALSE]
  )
  group_core[[group_code]] <-
    group_prevalence[[group_code]] >= primary_prevalence
}

group_core_summary <- data.frame(
  GroupCode = group_codes,
  Group = unname(group_display[group_codes]),
  Samples = vapply(
    group_codes,
    function(group_code) sum(metadata$Group == group_code),
    integer(1)
  ),
  CoreFeatures = vapply(group_core, sum, integer(1)),
  stringsAsFactors = FALSE
)

shared_three <- sum(
  group_core[["IW"]] & group_core[["CW"]] & group_core[["TW"]]
)
group_union <- sum(
  group_core[["IW"]] | group_core[["CW"]] | group_core[["TW"]]
)
data.frame(
  GlobalCore = sum(global_core),
  SharedThreeGroupCore = shared_three,
  GroupCoreUnion = group_union
)
group_core_summary

membership_flags <- data.frame(
  IW = group_core[["IW"]],
  CW = group_core[["CW"]],
  TW = group_core[["TW"]],
  stringsAsFactors = FALSE
)
membership_key <- apply(membership_flags, 1L, function(x) {
  active <- group_codes[as.logical(x[group_codes])]
  if (length(active) == 0L) "None" else paste(active, collapse = " + ")
})
membership_levels <- c(
  "IW", "CW", "TW", "IW + CW", "IW + TW", "CW + TW",
  "IW + CW + TW", "None"
)
membership_counts <- table(factor(membership_key, levels = membership_levels))
core_membership_patterns <- data.frame(
  Pattern = membership_levels,
  IW = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE),
  CW = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE),
  TW = c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
  CoreFeatures = as.integer(membership_counts),
  stringsAsFactors = FALSE
)
core_membership_patterns

detection_thresholds <- c(0, 0.0001, 0.0005, 0.001)
detection_labels <- c("Presence", "0.01%", "0.05%", "0.1%")
prevalence_thresholds <- c(0.50, 0.80, 0.90, 1.00)

threshold_rows <- list()
for (detection_i in seq_along(detection_thresholds)) {
  threshold <- detection_thresholds[[detection_i]]
  detected <- if (threshold == 0) {
    counts > 0L
  } else {
    relative_feature >= threshold
  }
  prevalence <- rowMeans(detected)
  for (prevalence_threshold in prevalence_thresholds) {
    mask <- prevalence >= prevalence_threshold
    threshold_rows[[length(threshold_rows) + 1L]] <- data.frame(
      DetectionThreshold = threshold,
      DetectionLabel = detection_labels[[detection_i]],
      PrevalenceThreshold = prevalence_threshold,
      PrevalenceLabel = scales::percent(prevalence_threshold, accuracy = 1),
      CoreFeatures = sum(mask),
      CoreAbundanceMass = sum(feature_equal_mean[mask]),
      stringsAsFactors = FALSE
    )
  }
}
core_threshold_sensitivity <- do.call(rbind, threshold_rows)
core_threshold_sensitivity

detection_count_audit <- do.call(rbind, lapply(seq_along(sample_depths), function(i) {
  data.frame(
    SampleID = names(sample_depths)[[i]],
    GroupCode = metadata[names(sample_depths)[[i]], "Group"],
    LibrarySize = sample_depths[[i]],
    DetectionThreshold = detection_thresholds,
    DetectionLabel = detection_labels,
    MinimumReads = ifelse(
      detection_thresholds == 0,
      1L,
      ceiling(sample_depths[[i]] * detection_thresholds)
    ),
    stringsAsFactors = FALSE
  )
}))
subset(detection_count_audit, DetectionLabel == "0.01%")[1:8, ]
range(subset(
  detection_count_audit,
  DetectionLabel == "0.01%"
)$MinimumReads)

rarefaction_depth <- 10000L
set.seed(primary_seed)
rarefied_counts <- t(vegan::rrarefy(
  t(counts),
  sample = rarefaction_depth
))
set.seed(primary_seed)
rarefied_repeat <- t(vegan::rrarefy(
  t(counts),
  sample = rarefaction_depth
))
stopifnot(
  identical(rarefied_counts, rarefied_repeat),
  all(colSums(rarefied_counts) == rarefaction_depth)
)

raw_presence_core <- rowMeans(counts > 0L) >= primary_prevalence
rarefied_one_read_prevalence <- rowMeans(rarefied_counts > 0L)
rarefied_two_read_prevalence <- rowMeans(rarefied_counts >= 2L)
rarefied_one_read_core <-
  rarefied_one_read_prevalence >= primary_prevalence
rarefied_two_read_core <-
  rarefied_two_read_prevalence >= primary_prevalence

depth_sensitivity <- data.frame(
  Definition = c(
    "Raw presence", "Raw 0.01%", "Rarefied one read", "Rarefied two reads"
  ),
  CoreFeatures = c(
    sum(raw_presence_core),
    sum(global_core),
    sum(rarefied_one_read_core),
    sum(rarefied_two_read_core)
  ),
  JaccardWithPrimary = c(
    jaccard(raw_presence_core, global_core),
    1,
    jaccard(rarefied_one_read_core, global_core),
    jaccard(rarefied_two_read_core, global_core)
  ),
  stringsAsFactors = FALSE
)
depth_sensitivity

feature_state <- ifelse(
  global_core,
  "Primary core",
  ifelse(
    rowSums(counts) <= 10L & positive_samples <= 2L,
    "Sampling-limited",
    ifelse(
      feature_equal_mean < primary_rare_threshold & feature_maximum >= 0.01,
      "Conditionally rare",
      ifelse(
        feature_maximum < primary_rare_threshold,
        "Persistently rare",
        "Intermediate"
      )
    )
  )
)
state_levels <- c(
  "Primary core", "Conditionally rare", "Persistently rare",
  "Sampling-limited", "Intermediate"
)
feature_state <- factor(feature_state, levels = state_levels)

feature_state_summary <- do.call(rbind, lapply(state_levels, function(state) {
  mask <- feature_state == state
  data.frame(
    FeatureState = state,
    Features = sum(mask),
    FeatureFraction = mean(mask),
    EqualSampleAbundanceMass = sum(feature_equal_mean[mask]),
    stringsAsFactors = FALSE
  )
}))
stopifnot(
  sum(feature_state_summary$Features) == nrow(counts),
  abs(sum(feature_state_summary$EqualSampleAbundanceMass) - 1) <= 1e-12
)
feature_state_summary

rare_thresholds <- c(0.0001, 0.0005, 0.001, 0.005)
rare_threshold_labels <- c("0.01%", "0.05%", "0.1%", "0.5%")
rare_sample_rows <- list()
for (threshold_i in seq_along(rare_thresholds)) {
  threshold <- rare_thresholds[[threshold_i]]
  rare_mass <- colSums(
    relative_feature *
      (relative_feature > 0 & relative_feature < threshold)
  )
  rare_sample_rows[[threshold_i]] <- data.frame(
    SampleID = colnames(counts),
    GroupCode = metadata[colnames(counts), "Group"],
    Group = unname(group_display[metadata[colnames(counts), "Group"]]),
    Type = metadata[colnames(counts), "Type"],
    Saline = metadata[colnames(counts), "Saline"],
    RareThreshold = threshold,
    ThresholdLabel = rare_threshold_labels[[threshold_i]],
    RareMass = rare_mass,
    stringsAsFactors = FALSE
  )
}
rare_mass_all_thresholds <- do.call(rbind, rare_sample_rows)
rare_mass_by_sample <- subset(
  rare_mass_all_thresholds,
  RareThreshold == primary_rare_threshold
)

rare_threshold_sensitivity <- do.call(
  rbind,
  lapply(seq_along(rare_thresholds), function(i) {
    x <- subset(
      rare_mass_all_thresholds,
      RareThreshold == rare_thresholds[[i]]
    )
    rows <- list(data.frame(
      GroupCode = "ALL",
      Group = "All samples",
      RareThreshold = rare_thresholds[[i]],
      ThresholdLabel = rare_threshold_labels[[i]],
      Samples = nrow(x),
      MeanRareMass = mean(x$RareMass),
      MinimumRareMass = min(x$RareMass),
      MaximumRareMass = max(x$RareMass),
      stringsAsFactors = FALSE
    ))
    for (group_code in group_codes) {
      values <- x$RareMass[x$GroupCode == group_code]
      rows[[length(rows) + 1L]] <- data.frame(
        GroupCode = group_code,
        Group = unname(group_display[group_code]),
        RareThreshold = rare_thresholds[[i]],
        ThresholdLabel = rare_threshold_labels[[i]],
        Samples = length(values),
        MeanRareMass = mean(values),
        MinimumRareMass = min(values),
        MaximumRareMass = max(values),
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  })
)
subset(rare_threshold_sensitivity, GroupCode == "ALL")

result_dir <- "results/28-core-rare-biosphere"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

feature_occupancy_abundance <- data.frame(
  FeatureID = rownames(counts),
  FeatureState = as.character(feature_state),
  TotalReads = rowSums(counts),
  PositiveSamples = positive_samples,
  RawPresence = positive_samples / ncol(counts),
  DetectedSamples = detected_samples,
  Prevalence = global_prevalence,
  EqualSampleMean = feature_equal_mean,
  PooledReadFraction = feature_pooled,
  MaximumSampleAbundance = feature_maximum,
  PrimaryCore = global_core,
  IWPrevalence = group_prevalence[["IW"]],
  CWPrevalence = group_prevalence[["CW"]],
  TWPrevalence = group_prevalence[["TW"]],
  BestRank = best_rank,
  BestTaxon = best_taxon,
  DisplayLabel = display_label,
  Lineage = lineage,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

core_feature_taxonomy <- data.frame(
  FeatureID = rownames(counts)[global_core],
  taxonomy_clean[global_core, , drop = FALSE],
  DetectedSamples = detected_samples[global_core],
  Prevalence = global_prevalence[global_core],
  EqualSampleMean = feature_equal_mean[global_core],
  MaximumSampleAbundance = feature_maximum[global_core],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
core_feature_taxonomy <- core_feature_taxonomy[
  order(core_feature_taxonomy$EqualSampleMean, decreasing = TRUE),
  , drop = FALSE
]

rarefaction_ledger <- data.frame(
  SampleID = colnames(counts),
  GroupCode = metadata[colnames(counts), "Group"],
  OriginalLibrarySize = sample_depths,
  RarefactionDepth = rarefaction_depth,
  RarefiedLibrarySize = colSums(rarefied_counts),
  ReadsDiscarded = sample_depths - colSums(rarefied_counts),
  Seed = primary_seed,
  stringsAsFactors = FALSE
)

readr::write_tsv(
  feature_occupancy_abundance,
  file.path(result_dir, "article-feature-occupancy-abundance.tsv")
)
readr::write_tsv(
  core_feature_taxonomy,
  file.path(result_dir, "article-core-feature-taxonomy.tsv")
)
readr::write_tsv(
  group_core_summary,
  file.path(result_dir, "article-group-core-summary.tsv")
)
readr::write_tsv(
  core_membership_patterns,
  file.path(result_dir, "article-core-membership-patterns.tsv")
)
readr::write_tsv(
  core_threshold_sensitivity,
  file.path(result_dir, "article-core-threshold-sensitivity.tsv")
)
readr::write_tsv(
  detection_count_audit,
  file.path(result_dir, "article-detection-count-audit.tsv")
)
readr::write_tsv(
  rarefaction_ledger,
  file.path(result_dir, "article-rarefaction-ledger.tsv")
)
readr::write_tsv(
  depth_sensitivity,
  file.path(result_dir, "article-depth-sensitivity.tsv")
)
readr::write_tsv(
  feature_state_summary,
  file.path(result_dir, "article-feature-state-summary.tsv")
)
readr::write_tsv(
  rare_mass_by_sample,
  file.path(result_dir, "article-rare-mass-by-sample.tsv")
)
readr::write_tsv(
  rare_threshold_sensitivity,
  file.path(result_dir, "article-rare-threshold-sensitivity.tsv")
)

state_colours <- c(
  `Primary core` = "#0072B2",
  `Conditionally rare` = "#D55E00",
  `Persistently rare` = "#56B4E9",
  `Sampling-limited` = "#9CA3AF",
  Intermediate = "#CC79A7"
)
occupancy_plot_data <- feature_occupancy_abundance
occupancy_plot_data$FeatureState <- factor(
  occupancy_plot_data$FeatureState,
  levels = state_levels
)
core_label_ids <- head(
  rownames(counts)[global_core][order(
    feature_equal_mean[global_core], decreasing = TRUE
  )],
  5L
)
conditional_mask <- feature_state == "Conditionally rare"
conditional_label_ids <- head(
  rownames(counts)[conditional_mask][order(
    feature_maximum[conditional_mask], decreasing = TRUE
  )],
  4L
)
occupancy_labels <- subset(
  occupancy_plot_data,
  FeatureID %in% c(core_label_ids, conditional_label_ids)
)
occupancy_labels$PointLabel <- paste0(
  occupancy_labels$DisplayLabel,
  " (",
  occupancy_labels$FeatureID,
  ")"
)

occupancy_plot <- ggplot2::ggplot(
  occupancy_plot_data,
  ggplot2::aes(
    x = Prevalence,
    y = EqualSampleMean,
    colour = FeatureState
  )
) +
  ggplot2::geom_hline(
    yintercept = primary_rare_threshold,
    colour = "#6B7280",
    linewidth = 0.45,
    linetype = "22"
  ) +
  ggplot2::geom_vline(
    xintercept = primary_prevalence,
    colour = "#1F2937",
    linewidth = 0.55,
    linetype = "22"
  ) +
  ggplot2::geom_point(size = 1.15, alpha = 0.58) +
  ggrepel::geom_text_repel(
    data = occupancy_labels,
    ggplot2::aes(label = PointLabel),
    seed = primary_seed,
    family = font_pub,
    size = 2.25,
    colour = "#17202A",
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = state_colours, drop = FALSE) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1.01),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.005, 0.01))
  ) +
  ggplot2::scale_y_log10(
    labels = scales::label_percent(accuracy = 0.001),
    expand = ggplot2::expansion(mult = c(0.04, 0.16))
  ) +
  ggplot2::annotation_logticks(sides = "l", linewidth = 0.25) +
  ggplot2::labs(
    title = "Core membership depends on an explicit abundance–occupancy rule",
    subtitle = "Primary core: at least 0.01% abundance in at least 80% of 90 samples",
    x = "Prevalence at the 0.01% detection threshold",
    y = "Equal-sample mean relative abundance",
    colour = "Feature state",
    caption = paste(
      "Each point is one feature/OTU; all-read denominators and equal sample weights are retained.",
      "The horizontal guide marks 0.1% mean abundance, not a core boundary."
    )
  ) +
  theme_pub(base_size = 8.5) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

save_pub(
  occupancy_plot,
  "figures/28-occupancy-abundance",
  width = 183,
  height = 132
)
occupancy_plot

membership_plot_data <- subset(
  core_membership_patterns,
  Pattern != "None"
)
membership_plot_data <- membership_plot_data[
  order(membership_plot_data$CoreFeatures, decreasing = TRUE),
  , drop = FALSE
]
membership_plot_data$X <- seq_len(nrow(membership_plot_data))
membership_plot_data$PatternLabel <- gsub(
  " \\+ ",
  " · ",
  membership_plot_data$Pattern
)
membership_y <- c(IW = -18, CW = -39, TW = -60)
membership_grid <- do.call(rbind, lapply(
  seq_len(nrow(membership_plot_data)),
  function(i) {
    data.frame(
      X = membership_plot_data$X[[i]],
      GroupCode = group_codes,
      Y = unname(membership_y[group_codes]),
      Active = as.logical(membership_plot_data[i, group_codes]),
      stringsAsFactors = FALSE
    )
  }
))
membership_segments <- do.call(rbind, lapply(
  seq_len(nrow(membership_plot_data)),
  function(i) {
    active_y <- membership_y[group_codes[
      as.logical(membership_plot_data[i, group_codes])
    ]]
    data.frame(
      X = membership_plot_data$X[[i]],
      YMin = min(active_y),
      YMax = max(active_y),
      stringsAsFactors = FALSE
    )
  }
))

membership_plot <- ggplot2::ggplot(
  membership_plot_data,
  ggplot2::aes(x = X)
) +
  ggplot2::geom_col(
    ggplot2::aes(y = CoreFeatures),
    width = 0.68,
    fill = "#0072B2"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = CoreFeatures, label = CoreFeatures),
    vjust = -0.35,
    family = font_pub,
    size = 2.8
  ) +
  ggplot2::geom_segment(
    data = membership_segments,
    ggplot2::aes(x = X, xend = X, y = YMin, yend = YMax),
    inherit.aes = FALSE,
    linewidth = 0.75,
    colour = "#374151"
  ) +
  ggplot2::geom_point(
    data = membership_grid,
    ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    size = 3.1,
    shape = 21,
    stroke = 0.45,
    colour = "#B8C0CC",
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = subset(membership_grid, Active),
    ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    size = 3.1,
    colour = "#17202A"
  ) +
  ggplot2::annotate(
    "text",
    x = 0.38,
    y = unname(membership_y),
    label = group_codes,
    hjust = 0,
    family = font_pub,
    size = 2.7,
    fontface = "bold"
  ) +
  ggplot2::scale_x_continuous(
    breaks = membership_plot_data$X,
    labels = membership_plot_data$PatternLabel,
    limits = c(0.3, nrow(membership_plot_data) + 0.55),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(0, 50, 100, 150),
    labels = c("0", "50", "100", "150"),
    limits = c(-72, 190),
    expand = ggplot2::expansion(mult = c(0, 0.015))
  ) +
  ggplot2::labs(
    title = "Group cores overlap, but they are not the global core",
    subtitle = "The same 0.01% detection and 80% prevalence rules are applied within each 30-sample group",
    x = "Group-core membership pattern",
    y = "Features",
    caption = paste(
      "The seven non-empty intersections contain 425 features; 59 are core in all three groups.",
      "IW, inland; CW, coastal; TW, Tibetan wetland."
    )
  ) +
  theme_pub(base_size = 8.5) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 24, hjust = 1),
    axis.ticks.y = ggplot2::element_blank()
  )

save_pub(
  membership_plot,
  "figures/28-group-core-membership",
  width = 165,
  height = 118
)
membership_plot

fixed_depth_rows <- list()
for (prevalence_threshold in prevalence_thresholds) {
  fixed_depth_rows[[length(fixed_depth_rows) + 1L]] <- data.frame(
    DetectionLabel = c("One read", "Two reads"),
    PrevalenceThreshold = prevalence_threshold,
    PrevalenceLabel = scales::percent(prevalence_threshold, accuracy = 1),
    CoreFeatures = c(
      sum(rarefied_one_read_prevalence >= prevalence_threshold),
      sum(rarefied_two_read_prevalence >= prevalence_threshold)
    ),
    stringsAsFactors = FALSE
  )
}
fixed_depth_sensitivity <- do.call(rbind, fixed_depth_rows)
threshold_plot_data <- rbind(
  data.frame(
    Panel = "Raw libraries",
    DetectionLabel = core_threshold_sensitivity$DetectionLabel,
    PrevalenceLabel = core_threshold_sensitivity$PrevalenceLabel,
    CoreFeatures = core_threshold_sensitivity$CoreFeatures,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Panel = "10,000-read branch",
    DetectionLabel = fixed_depth_sensitivity$DetectionLabel,
    PrevalenceLabel = fixed_depth_sensitivity$PrevalenceLabel,
    CoreFeatures = fixed_depth_sensitivity$CoreFeatures,
    stringsAsFactors = FALSE
  )
)
threshold_plot_data$DetectionLabel <- factor(
  threshold_plot_data$DetectionLabel,
  levels = c(detection_labels, "One read", "Two reads")
)
threshold_plot_data$PrevalenceLabel <- factor(
  threshold_plot_data$PrevalenceLabel,
  levels = scales::percent(prevalence_thresholds, accuracy = 1)
)
threshold_plot_data$Panel <- factor(
  threshold_plot_data$Panel,
  levels = c("Raw libraries", "10,000-read branch")
)
primary_tiles <- subset(
  threshold_plot_data,
  (Panel == "Raw libraries" &
    DetectionLabel == "0.01%" & PrevalenceLabel == "80%") |
    (Panel == "10,000-read branch" &
      DetectionLabel == "One read" & PrevalenceLabel == "80%")
)
primary_rarefied_jaccard <- jaccard(
  global_core,
  rarefied_one_read_core
)

threshold_plot <- ggplot2::ggplot(
  threshold_plot_data,
  ggplot2::aes(
    x = DetectionLabel,
    y = PrevalenceLabel,
    fill = CoreFeatures
  )
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
  ggplot2::geom_tile(
    data = primary_tiles,
    fill = NA,
    colour = "#D55E00",
    linewidth = 1.15
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(CoreFeatures)),
    family = font_pub,
    size = 3
  ) +
  ggplot2::facet_grid(
    . ~ Panel,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c(
      "#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"
    ),
    trans = "sqrt",
    labels = scales::label_comma()
  ) +
  ggplot2::labs(
    title = "Core size changes with detection, prevalence and sequencing depth",
    subtitle = "Raw thresholds use all reads; the fixed-depth branch uses seed 20260728 and is a sensitivity analysis",
    x = "Detection rule",
    y = "Minimum prevalence",
    fill = "Core features",
    caption = paste0(
      "Orange outlines compare the primary raw core (97) with the rarefied one-read core (167) at 80%; ",
      "Jaccard = ", sprintf("%.3f", primary_rarefied_jaccard), "."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 24, hjust = 1),
    strip.text = ggplot2::element_text(size = 9)
  )

save_pub(
  threshold_plot,
  "figures/28-threshold-depth-sensitivity",
  width = 170,
  height = 112
)
threshold_plot

group_colours <- c(
  `Inland wetland` = "#E69F00",
  `Coastal wetland` = "#0072B2",
  `Tibetan wetland` = "#009E73"
)
rare_plot_data <- rare_mass_all_thresholds
rare_plot_data$ThresholdLabel <- factor(
  rare_plot_data$ThresholdLabel,
  levels = rare_threshold_labels
)
rare_plot_data$Group <- factor(
  rare_plot_data$Group,
  levels = unname(group_display[group_codes])
)
rare_group_means <- subset(
  rare_threshold_sensitivity,
  GroupCode %in% group_codes
)
rare_group_means$ThresholdLabel <- factor(
  rare_group_means$ThresholdLabel,
  levels = rare_threshold_labels
)
rare_group_means$Group <- factor(
  rare_group_means$Group,
  levels = unname(group_display[group_codes])
)

rare_plot <- ggplot2::ggplot(
  rare_plot_data,
  ggplot2::aes(
    x = ThresholdLabel,
    y = RareMass,
    fill = Group,
    colour = Group
  )
) +
  ggplot2::annotate(
    "rect",
    xmin = 2.5,
    xmax = 3.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "#F0E442",
    alpha = 0.10
  ) +
  ggplot2::geom_boxplot(
    width = 0.68,
    outlier.shape = NA,
    alpha = 0.35,
    position = ggplot2::position_dodge(width = 0.78),
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    size = 0.85,
    alpha = 0.30,
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.11,
      dodge.width = 0.78,
      seed = primary_seed
    )
  ) +
  ggplot2::geom_line(
    data = rare_group_means,
    ggplot2::aes(y = MeanRareMass, group = Group),
    position = ggplot2::position_dodge(width = 0.78),
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    data = rare_group_means,
    ggplot2::aes(y = MeanRareMass),
    position = ggplot2::position_dodge(width = 0.78),
    size = 2.2,
    shape = 21,
    stroke = 0.7
  ) +
  ggplot2::scale_fill_manual(values = group_colours, drop = FALSE) +
  ggplot2::scale_colour_manual(values = group_colours, drop = FALSE) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 0.8, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 0.93)) +
  ggplot2::labs(
    title = "The measured rare-biosphere mass is threshold dependent",
    subtitle = "Local rare mass sums features below the selected abundance threshold within each sample",
    x = "Local rare-abundance threshold",
    y = "Relative-abundance mass below threshold",
    fill = "Wetland group",
    colour = "Wetland group",
    caption = paste(
      "The highlighted 0.1% rule gives a 44.2% equal-sample mean (range 23.3–59.4%).",
      "Points are samples; lines connect descriptive group means without hypothesis tests."
    )
  ) +
  theme_pub(base_size = 8.8) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.major.x = ggplot2::element_blank()
  )

save_pub(
  rare_plot,
  "figures/28-rare-biosphere-mass",
  width = 170,
  height = 116
)
rare_plot

expected_threshold_counts <- c(
  1499L, 334L, 166L, 27L,
  404L, 97L, 53L, 10L,
  48L, 9L, 5L, 0L,
  16L, 0L, 0L, 0L
)
stopifnot(
  sum(global_core) == 97L,
  abs(sum(feature_equal_mean[global_core]) - 0.170361211031733) <= 1e-12,
  identical(unname(group_core_summary$CoreFeatures), c(166L, 160L, 278L)),
  shared_three == 59L,
  group_union == 425L,
  identical(core_threshold_sensitivity$CoreFeatures, expected_threshold_counts),
  identical(depth_sensitivity$CoreFeatures, c(334L, 97L, 167L, 52L)),
  abs(primary_rarefied_jaccard - 0.580838323353293) <= 1e-12,
  identical(feature_state_summary$Features, c(97L, 200L, 8407L, 1155L, 3769L)),
  abs(mean(rare_mass_by_sample$RareMass) - 0.441501476980241) <= 1e-12
)

audit_path <- file.path(result_dir, "validation-checks.tsv")
if (file.exists(audit_path)) {
  validation_checks <- readr::read_tsv(
    audit_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  stopifnot(
    nrow(validation_checks) == 320L,
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
    Passed = 10L,
    Total = 10L
  )
}
