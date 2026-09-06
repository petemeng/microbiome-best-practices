# 方差解构：envfit/Mantel/VPA 变异分解
# Run sequentially in a new working directory.
# Required packages: ggplot2, jsonlite, permute, ragg, scales, svglite, vegan.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/environment.tsv",
  "data/small/metadata.tsv",
  "data/small/otutab.tsv",
  "data/small/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)
library(vegan)

set.seed(20260725)
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

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    na = character(),
    col_types = readr::cols(.default = readr::col_character())
  )
  ids <- as.character(x[[1]])
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1], check.names = FALSE)
  rownames(out) <- ids
  out
}

input_paths <- c(
  otutab = "data/small/otutab.tsv",
  taxonomy = "data/small/taxonomy.tsv",
  metadata = "data/small/metadata.tsv",
  environment = "data/small/environment.tsv"
)
expected_sha256 <- c(
  otutab = "76fa79c38da889f35978dc86da4641a270961746709ff38049ee5f67e3c6f7a3",
  taxonomy = "725280bb9a0cd9bda7b540022e92af945ceed52527f8b2d220b055b4489f6901",
  metadata = "df24771dccf27607ddf922c6bca2cafa876d946fbe2e09d14b601accce66ba64",
  environment = "45814f6724e8abd901b8063ba227c42716777cd0749ac130e9cae99f0477f7a8"
)
observed_sha256 <- vapply(
  input_paths,
  digest::digest,
  character(1),
  algo = "sha256",
  serialize = FALSE,
  file = TRUE
)
stopifnot(identical(unname(observed_sha256), unname(expected_sha256)))

otu_frame <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])
environment_frame <- read_keyed_tsv(input_paths[["environment"]])

counts_character <- as.matrix(otu_frame)
counts <- matrix(
  as.numeric(counts_character),
  nrow = nrow(counts_character),
  ncol = ncol(counts_character),
  dimnames = dimnames(counts_character)
)
environment_character <- as.matrix(environment_frame)
environment_all <- matrix(
  as.numeric(environment_character),
  nrow = nrow(environment_character),
  ncol = ncol(environment_character),
  dimnames = dimnames(environment_character)
)

sample_ids <- colnames(counts)
feature_ids <- rownames(counts)
environment_variables <- c(
  "Latitude", "Longitude", "Altitude", "Temperature", "Precipitation",
  "TOC", "NH4", "NO3", "pH", "Conductivity", "TN"
)
stopifnot(
  setequal(feature_ids, rownames(taxonomy)),
  setequal(sample_ids, rownames(metadata)),
  all(sample_ids %in% rownames(environment_all)),
  all(environment_variables %in% colnames(environment_all)),
  all(is.finite(counts)),
  all(counts >= 0),
  all(abs(counts - round(counts)) < .Machine$double.eps^0.5)
)

taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]
environment <- as.data.frame(
  environment_all[sample_ids, environment_variables, drop = FALSE],
  check.names = FALSE
)

group_labels <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan Plateau"
)
metadata$GroupCode <- metadata$Group
metadata$Group <- factor(
  unname(group_labels[metadata$GroupCode]),
  levels = unname(group_labels)
)
metadata$Saline <- factor(metadata$Saline)
metadata$Type <- factor(metadata$Type)

dim(counts)
sum(counts)
head(counts[, 1:5])
head(taxonomy)
head(metadata)
head(environment)

design_cross <- with(metadata, table(Group, Saline))
design_cross

group_determines_saline <- all(rowSums(design_cross > 0) == 1L)
repeat_columns <- intersect(
  c("SubjectID", "PairID", "Block", "SiteID", "Time", "Batch"),
  colnames(metadata)
)
stopifnot(group_determines_saline, length(repeat_columns) == 0L)

correlation_matrix <- stats::cor(environment, method = "pearson")
correlation_index <- expand.grid(
  Variable1 = rownames(correlation_matrix),
  Variable2 = colnames(correlation_matrix),
  stringsAsFactors = FALSE
)
correlation_index$PearsonR <- as.vector(correlation_matrix)
correlation_index$AbsR <- abs(correlation_index$PearsonR)
variable_order <- setNames(seq_along(environment_variables), environment_variables)
correlation_index$UniquePair <-
  variable_order[correlation_index$Variable1] <
  variable_order[correlation_index$Variable2]
high_correlation_pairs <- subset(
  correlation_index,
  UniquePair & AbsR >= 0.80,
  select = c(Variable1, Variable2, PearsonR, AbsR)
)
high_correlation_pairs

primary_seed <- 20260725L
permutation_count <- 999L
min_prevalence <- 0.05
minimum_samples <- ceiling(min_prevalence * ncol(counts))

primary_keep <- rowSums(counts > 0) >= minimum_samples
primary_counts <- counts[primary_keep, , drop = FALSE]
community_counts <- t(primary_counts)
community_relative <- sweep(
  community_counts,
  1,
  rowSums(community_counts),
  "/"
)
community_hellinger <- vegan::decostand(
  community_counts,
  method = "hellinger"
)
bray_distance <- vegan::vegdist(community_relative, method = "bray")

permutation_control <- permute::how(
  within = permute::Within(type = "free"),
  nperm = permutation_count
)
set.seed(primary_seed)
permutation_matrix <- permute::shuffleSet(
  n = length(sample_ids),
  nset = permutation_count,
  control = permutation_control
)

stopifnot(
  sum(primary_keep) == 11113L,
  minimum_samples == 5L,
  max(abs(rowSums(community_relative) - 1)) < 1e-12,
  all(is.finite(community_hellinger)),
  nrow(unique(as.data.frame(permutation_matrix))) == 999L
)
data.frame(
  Samples = nrow(community_counts),
  FeaturesBefore = nrow(counts),
  FeaturesAfter = nrow(primary_counts),
  MinimumSamples = minimum_samples,
  Permutations = nrow(permutation_matrix),
  MinimumP = 1 / (permutation_count + 1)
)

pcoa <- vegan::wcmdscale(
  bray_distance,
  k = 2,
  eig = TRUE,
  add = "lingoes",
  x.ret = TRUE
)
pcoa_scores <- data.frame(
  SampleID = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  Group = metadata[rownames(pcoa$points), "Group"],
  stringsAsFactors = FALSE
)
positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
pcoa_axis_percent <- 100 * pcoa$eig[1:2] / sum(positive_eigenvalues)

environment_z <- as.data.frame(scale(environment))
set.seed(primary_seed)
envfit_result <- vegan::envfit(
  pcoa$points,
  environment_z,
  permutations = permutation_matrix,
  w = rep(1, nrow(pcoa$points)),
  na.rm = FALSE
)
envfit_scaled_scores <- vegan::scores(envfit_result, display = "vectors")
envfit_results <- data.frame(
  Variable = rownames(envfit_result$vectors$arrows),
  Axis1Direction = envfit_result$vectors$arrows[, 1],
  Axis2Direction = envfit_result$vectors$arrows[, 2],
  PlotAxis1 = envfit_scaled_scores[, 1],
  PlotAxis2 = envfit_scaled_scores[, 2],
  TwoAxisR2 = as.numeric(envfit_result$vectors$r),
  PValue = as.numeric(envfit_result$vectors$pvals),
  stringsAsFactors = FALSE
)
envfit_results$PAdjustedBH <- stats::p.adjust(
  envfit_results$PValue,
  method = "BH"
)
envfit_results$RejectBH05 <- envfit_results$PAdjustedBH < 0.05
envfit_results <- envfit_results[
  order(envfit_results$PAdjustedBH, -envfit_results$TwoAxisR2),
  ,
  drop = FALSE
]
rownames(envfit_results) <- NULL

envfit_results[, c(
  "Variable", "TwoAxisR2", "PValue", "PAdjustedBH", "RejectBH05"
)]

edaphic_names <- c("pH", "TOC", "NH4", "NO3", "Conductivity")
climate_names <- c("Temperature", "Precipitation")
edaphic_z <- as.data.frame(scale(environment[, edaphic_names, drop = FALSE]))
climate_z <- as.data.frame(scale(environment[, climate_names, drop = FALSE]))

haversine_dist <- function(latitude, longitude, ids = names(latitude)) {
  rad <- pi / 180
  lat <- latitude * rad
  lon <- longitude * rad
  n <- length(lat)
  out <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) {
    j <- (i + 1L):n
    dlat <- lat[j] - lat[i]
    dlon <- lon[j] - lon[i]
    a <- sin(dlat / 2)^2 +
      cos(lat[i]) * cos(lat[j]) * sin(dlon / 2)^2
    out[i, j] <- 6371.0088 * 2 * atan2(
      sqrt(a),
      sqrt(pmax(0, 1 - a))
    )
    out[j, i] <- out[i, j]
  }
  dimnames(out) <- list(ids, ids)
  stats::as.dist(out)
}

environment_distances <- list(
  Edaphic = stats::dist(edaphic_z, method = "euclidean"),
  Climate = stats::dist(climate_z, method = "euclidean"),
  Geographic = haversine_dist(
    environment$Latitude,
    environment$Longitude,
    sample_ids
  )
)

mantel_results <- do.call(
  rbind,
  lapply(seq_along(environment_distances), function(i) {
    set.seed(primary_seed)
    fit <- vegan::mantel(
      bray_distance,
      environment_distances[[i]],
      method = "spearman",
      permutations = permutation_matrix,
      na.rm = FALSE
    )
    data.frame(
      DistanceBlock = names(environment_distances)[i],
      Statistic = as.numeric(fit$statistic),
      PValue = as.numeric(fit$signif),
      Permutations = permutation_count,
      stringsAsFactors = FALSE
    )
  })
)
mantel_results$PAdjustedHolm <- stats::p.adjust(
  mantel_results$PValue,
  method = "holm"
)
mantel_results$RejectHolm05 <- mantel_results$PAdjustedHolm < 0.05
mantel_results

set.seed(primary_seed)
partial_mantel <- vegan::mantel.partial(
  bray_distance,
  environment_distances$Edaphic,
  environment_distances$Geographic,
  method = "pearson",
  permutations = permutation_matrix,
  na.rm = FALSE
)
partial_mantel_sensitivity <- data.frame(
  Analysis = "Edaphic distance | geographic distance",
  Method = "Pearson partial Mantel",
  Statistic = as.numeric(partial_mantel$statistic),
  PValue = as.numeric(partial_mantel$signif),
  Permutations = permutation_count,
  PrimaryInference = FALSE,
  stringsAsFactors = FALSE
)
partial_mantel_sensitivity

group_df <- data.frame(Group = metadata$Group, row.names = sample_ids)
variation_partition_object <- vegan::varpart(
  community_hellinger,
  edaphic_z,
  group_df
)
individual_fractions <- variation_partition_object$part$indfract
variation_partition <- data.frame(
  FractionCode = c("[a]", "[b]", "[c]", "[d]"),
  Fraction = c(
    "Pure edaphic | wetland group",
    "Pure wetland group | edaphic",
    "Shared edaphic + wetland group",
    "Unexplained"
  ),
  AdjustedR2 = as.numeric(individual_fractions$Adj.R.squared),
  Testable = as.logical(individual_fractions$Testable),
  stringsAsFactors = FALSE
)

model_data <- cbind(edaphic_z, Group = metadata$Group, Saline = metadata$Saline)
combined_rda <- vegan::rda(
  community_hellinger ~ pH + TOC + NH4 + NO3 + Conductivity + Group,
  data = model_data
)
pure_edaphic_rda <- vegan::rda(
  community_hellinger ~ pH + TOC + NH4 + NO3 + Conductivity +
    Condition(Group),
  data = model_data
)
pure_group_rda <- vegan::rda(
  community_hellinger ~ Group +
    Condition(pH + TOC + NH4 + NO3 + Conductivity),
  data = model_data
)

combined_r2 <- vegan::RsquareAdj(combined_rda)
set.seed(primary_seed)
pure_edaphic_test <- anova(
  pure_edaphic_rda,
  permutations = permutation_matrix
)
set.seed(primary_seed)
pure_group_test <- anova(
  pure_group_rda,
  permutations = permutation_matrix
)
variation_fraction_tests <- data.frame(
  Fraction = variation_partition$Fraction[1:2],
  AdjustedR2 = c(
    vegan::RsquareAdj(pure_edaphic_rda)$adj.r.squared,
    vegan::RsquareAdj(pure_group_rda)$adj.r.squared
  ),
  FValue = c(pure_edaphic_test$F[1], pure_group_test$F[1]),
  PValue = c(
    pure_edaphic_test$`Pr(>F)`[1],
    pure_group_test$`Pr(>F)`[1]
  ),
  stringsAsFactors = FALSE
)
variation_fraction_tests$PAdjustedHolm <- stats::p.adjust(
  variation_fraction_tests$PValue,
  method = "holm"
)
variation_fraction_tests$RejectHolm05 <-
  variation_fraction_tests$PAdjustedHolm < 0.05

variation_partition
data.frame(
  CombinedRawR2 = combined_r2$r.squared,
  CombinedAdjustedR2 = combined_r2$adj.r.squared
)
variation_fraction_tests

combined_vif <- vegan::vif.cca(combined_rda)
vif_audit <- data.frame(
  Term = names(combined_vif),
  VIF = as.numeric(combined_vif),
  Status = ifelse(combined_vif < 5, "PASS", "REVIEW"),
  stringsAsFactors = FALSE
)

edaphic_unconditioned_rda <- vegan::rda(
  community_hellinger ~ pH + TOC + NH4 + NO3 + Conductivity,
  data = model_data
)
edaphic_condition_saline_rda <- vegan::rda(
  community_hellinger ~ pH + TOC + NH4 + NO3 + Conductivity +
    Condition(Saline),
  data = model_data
)
conditioning_sensitivity <- data.frame(
  Conditioning = c("None", "Binary salinity", "Three-level wetland group"),
  AdjustedR2 = c(
    vegan::RsquareAdj(edaphic_unconditioned_rda)$adj.r.squared,
    vegan::RsquareAdj(edaphic_condition_saline_rda)$adj.r.squared,
    vegan::RsquareAdj(pure_edaphic_rda)$adj.r.squared
  ),
  Primary = c(FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

prevalence_thresholds <- c(0, 0.05, 0.10, 0.20)
prevalence_sensitivity <- do.call(
  rbind,
  lapply(prevalence_thresholds, function(threshold) {
    keep <- if (threshold == 0) {
      rep(TRUE, nrow(counts))
    } else {
      rowSums(counts > 0) >= ceiling(threshold * ncol(counts))
    }
    branch_counts <- t(counts[keep, , drop = FALSE])
    branch_relative <- sweep(
      branch_counts,
      1,
      rowSums(branch_counts),
      "/"
    )
    branch_hellinger <- vegan::decostand(
      branch_counts,
      method = "hellinger"
    )
    branch_vpa <- vegan::varpart(branch_hellinger, edaphic_z, group_df)
    fractions <- as.numeric(branch_vpa$part$indfract$Adj.R.squared)
    set.seed(primary_seed)
    branch_mantel <- vegan::mantel(
      vegan::vegdist(branch_relative, method = "bray"),
      environment_distances$Edaphic,
      method = "spearman",
      permutations = permutation_matrix
    )
    data.frame(
      PrevalenceThreshold = threshold,
      MinimumSamples = ifelse(
        threshold == 0,
        1L,
        ceiling(threshold * ncol(counts))
      ),
      Features = sum(keep),
      PureEdaphicAdjustedR2 = fractions[1],
      PureGroupAdjustedR2 = fractions[2],
      SharedAdjustedR2 = fractions[3],
      UnexplainedAdjustedR2 = fractions[4],
      EdaphicMantelRho = as.numeric(branch_mantel$statistic),
      EdaphicMantelP = as.numeric(branch_mantel$signif),
      Primary = threshold == min_prevalence,
      stringsAsFactors = FALSE
    )
  })
)

vif_audit
conditioning_sensitivity
prevalence_sensitivity

correlation_index$Diagonal <-
  correlation_index$Variable1 == correlation_index$Variable2
correlation_index$HighCollinearity <-
  correlation_index$UniquePair & correlation_index$AbsR >= 0.80
correlation_plot_data <- correlation_index[
  variable_order[correlation_index$Variable1] <=
    variable_order[correlation_index$Variable2],
  ,
  drop = FALSE
]
correlation_plot_data$Variable1 <- factor(
  correlation_plot_data$Variable1,
  levels = rev(environment_variables)
)
correlation_plot_data$Variable2 <- factor(
  correlation_plot_data$Variable2,
  levels = environment_variables
)
correlation_plot_data$Label <- ifelse(
  correlation_plot_data$Diagonal,
  "1.00",
  ifelse(
    correlation_plot_data$AbsR >= 0.50,
    sprintf("%.2f", correlation_plot_data$PearsonR),
    ""
  )
)

correlation_plot <- ggplot2::ggplot(
  correlation_plot_data,
  ggplot2::aes(Variable2, Variable1, fill = PearsonR)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::geom_tile(
    data = subset(correlation_plot_data, HighCollinearity),
    fill = NA,
    colour = "#111111",
    linewidth = 0.85
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label),
    size = 2.5,
    colour = "#111111",
    family = font_pub
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Environmental covariates are not independent",
    subtitle = "Black outlines mark |r| >= 0.80; labels are shown for |r| >= 0.50",
    x = NULL,
    y = NULL,
    caption = paste(
      "All 11 variables enter the descriptive envfit family.",
      "TN is excluded from the primary edaphic set because it is nearly redundant with TOC.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 8.4) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 42, hjust = 1, size = 7.4),
    axis.text.y = ggplot2::element_text(size = 7.4),
    legend.position = "right"
  )

save_pub(
  correlation_plot,
  "figures/25-environment-correlation",
  width = 183,
  height = 126
)
correlation_plot

group_palette <- c(
  "Inland wetland" = "#0072B2",
  "Coastal wetland" = "#D55E00",
  "Tibetan Plateau" = "#009E73"
)
group_shapes <- c(
  "Inland wetland" = 16,
  "Coastal wetland" = 17,
  "Tibetan Plateau" = 15
)
envfit_plot_arrows <- subset(envfit_results, RejectBH05)
envfit_plot_arrows <- head(
  envfit_plot_arrows[order(-envfit_plot_arrows$TwoAxisR2), ],
  7
)
arrow_multiplier <- min(
  diff(range(pcoa_scores$PCoA1)) * 0.40 /
    max(abs(envfit_plot_arrows$PlotAxis1)),
  diff(range(pcoa_scores$PCoA2)) * 0.40 /
    max(abs(envfit_plot_arrows$PlotAxis2))
)
envfit_plot_arrows$ArrowX <-
  envfit_plot_arrows$PlotAxis1 * arrow_multiplier
envfit_plot_arrows$ArrowY <-
  envfit_plot_arrows$PlotAxis2 * arrow_multiplier

envfit_label_offsets <- data.frame(
  Variable = c(
    "Temperature", "Latitude", "Altitude", "Longitude",
    "Conductivity", "Precipitation", "pH"
  ),
  OffsetX = c(0.025, -0.015, -0.025, 0.005, -0.060, -0.035, 0.020),
  OffsetY = c(0.030, -0.045, 0.035, -0.030, 0.035, -0.055, 0.040),
  Hjust = c(0, 1, 0.5, 0.5, 1, 0.5, 0),
  stringsAsFactors = FALSE
)
offset_index <- match(envfit_plot_arrows$Variable, envfit_label_offsets$Variable)
envfit_plot_arrows$LabelX <-
  envfit_plot_arrows$ArrowX + envfit_label_offsets$OffsetX[offset_index]
envfit_plot_arrows$LabelY <-
  envfit_plot_arrows$ArrowY + envfit_label_offsets$OffsetY[offset_index]
envfit_plot_arrows$LabelHjust <- envfit_label_offsets$Hjust[offset_index]

envfit_plot <- ggplot2::ggplot(
  pcoa_scores,
  ggplot2::aes(PCoA1, PCoA2, colour = Group, shape = Group)
) +
  ggplot2::stat_ellipse(
    type = "t",
    level = 0.95,
    linewidth = 0.55,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(size = 2.4, alpha = 0.83) +
  ggplot2::geom_segment(
    data = envfit_plot_arrows,
    ggplot2::aes(x = 0, y = 0, xend = ArrowX, yend = ArrowY),
    inherit.aes = FALSE,
    colour = "#222222",
    linewidth = 0.55,
    arrow = grid::arrow(length = grid::unit(2.1, "mm"), type = "closed")
  ) +
  ggplot2::geom_segment(
    data = envfit_plot_arrows,
    ggplot2::aes(x = ArrowX, y = ArrowY, xend = LabelX, yend = LabelY),
    inherit.aes = FALSE,
    colour = "#777777",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    data = envfit_plot_arrows,
    ggplot2::aes(
      x = LabelX,
      y = LabelY,
      label = Variable,
      hjust = LabelHjust
    ),
    inherit.aes = FALSE,
    colour = "#111111",
    family = font_pub,
    fontface = "bold",
    size = 2.8
  ) +
  ggplot2::scale_colour_manual(values = group_palette, name = "Wetland group") +
  ggplot2::scale_shape_manual(values = group_shapes, name = "Wetland group") +
  ggplot2::coord_equal(clip = "off") +
  ggplot2::labs(
    title = "envfit overlays external gradients on an unconstrained PCoA",
    subtitle = paste(
      "Top 7 of",
      sum(envfit_results$RejectBH05),
      "BH-significant variables; all 11 tests remain in the audit table"
    ),
    x = sprintf("PCoA1 (%.1f%% positive inertia)", pcoa_axis_percent[1]),
    y = sprintf("PCoA2 (%.1f%% positive inertia)", pcoa_axis_percent[2]),
    caption = paste(
      "Arrow direction and length describe association with this two-axis display.",
      "They are not causal effects or full-community explained fractions.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    legend.position = "top",
    legend.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )

save_pub(
  envfit_plot,
  "figures/25-envfit-pcoa",
  width = 183,
  height = 122
)
envfit_plot

format_p <- function(p) {
  ifelse(p <= 0.001, "<=0.001", sprintf("%.3f", p))
}
mantel_plot_data <- data.frame(
  Label = c(
    "Geographic distance",
    "Climate distance",
    "Edaphic distance",
    "Edaphic | geographic"
  ),
  Statistic = c(
    mantel_results$Statistic[mantel_results$DistanceBlock == "Geographic"],
    mantel_results$Statistic[mantel_results$DistanceBlock == "Climate"],
    mantel_results$Statistic[mantel_results$DistanceBlock == "Edaphic"],
    partial_mantel_sensitivity$Statistic
  ),
  PValue = c(
    mantel_results$PAdjustedHolm[mantel_results$DistanceBlock == "Geographic"],
    mantel_results$PAdjustedHolm[mantel_results$DistanceBlock == "Climate"],
    mantel_results$PAdjustedHolm[mantel_results$DistanceBlock == "Edaphic"],
    partial_mantel_sensitivity$PValue
  ),
  Family = c(rep("Primary Mantel family", 3), "Partial Mantel sensitivity"),
  stringsAsFactors = FALSE
)
mantel_plot_data$Label <- factor(
  mantel_plot_data$Label,
  levels = rev(mantel_plot_data$Label)
)
mantel_plot_data$Annotation <- sprintf(
  "r = %.3f; P %s",
  mantel_plot_data$Statistic,
  format_p(mantel_plot_data$PValue)
)

mantel_plot <- ggplot2::ggplot(
  mantel_plot_data,
  ggplot2::aes(Statistic, Label, colour = Family, shape = Family)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = Statistic, yend = Label),
    linewidth = 0.8,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(size = 3.2) +
  ggplot2::geom_text(
    ggplot2::aes(label = Annotation),
    hjust = -0.08,
    colour = "#111111",
    family = font_pub,
    size = 3.1,
    show.legend = FALSE
  ) +
  ggplot2::geom_vline(xintercept = 0, colour = "#555555", linewidth = 0.4) +
  ggplot2::scale_colour_manual(
    values = c(
      "Primary Mantel family" = "#0072B2",
      "Partial Mantel sensitivity" = "#D55E00"
    ),
    name = NULL
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Primary Mantel family" = 16,
      "Partial Mantel sensitivity" = 17
    ),
    name = NULL
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 0.72), clip = "off") +
  ggplot2::labs(
    title = "Mantel statistics compare two distance structures",
    subtitle = "Primary Spearman tests use Holm correction across three prespecified blocks",
    x = "Matrix correlation",
    y = NULL,
    caption = paste(
      "The partial Mantel row uses Pearson correlation and is a sensitivity analysis only.",
      "Distance correlation does not identify an independent environmental cause.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "top",
    plot.margin = ggplot2::margin(5.5, 38, 5.5, 5.5)
  )

save_pub(
  mantel_plot,
  "figures/25-mantel-distance-blocks",
  width = 183,
  height = 96
)
mantel_plot

variation_plot_data <- variation_partition
variation_plot_data$Fraction <- factor(
  variation_plot_data$Fraction,
  levels = rev(variation_plot_data$Fraction)
)
variation_palette <- c(
  "Pure edaphic | wetland group" = "#0072B2",
  "Pure wetland group | edaphic" = "#D55E00",
  "Shared edaphic + wetland group" = "#CC79A7",
  "Unexplained" = "#9CA3AF"
)

variation_plot <- ggplot2::ggplot(
  variation_plot_data,
  ggplot2::aes(AdjustedR2, Fraction, fill = Fraction)
) +
  ggplot2::geom_col(width = 0.64, colour = "#333333", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.4) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = scales::percent(AdjustedR2, accuracy = 0.1),
      hjust = ifelse(AdjustedR2 >= 0, -0.08, 1.08)
    ),
    family = font_pub,
    size = 3.2
  ) +
  ggplot2::scale_fill_manual(values = variation_palette, guide = "none") +
  ggplot2::scale_x_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.02, 0.13))
  ) +
  ggplot2::labs(
    title = "Adjusted variation partition separates unique and shared fractions",
    subtitle = sprintf(
      "Combined adjusted R2 = %.3f; both pure fractions Holm P <= 0.002",
      combined_r2$adj.r.squared
    ),
    x = "Adjusted fraction of Hellinger community variation",
    y = NULL,
    caption = paste(
      "Shared variation is not directly testable or attributable to one predictor set.",
      "Wetland group is a broad context variable and is completely confounded with salinity here.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5.5, 18, 5.5, 5.5)
  )

save_pub(
  variation_plot,
  "figures/25-variation-partition",
  width = 183,
  height = 100
)
variation_plot

audit_path <-
  "results/25-environment-variance/environment-variance-summary.json"

if (file.exists(audit_path)) {
  audit_summary <- jsonlite::read_json(audit_path)
  stopifnot(
    audit_summary$checks_total == 221L,
    audit_summary$checks_passed == 221L,
    audit_summary$checks_failed == 0L,
    audit_summary$dataset$samples == 90L,
    audit_summary$standardization$retained_features == 11113L,
    audit_summary$graphics$primary_figures == 4L,
    audit_summary$graphics$format_files == 16L
  )
  audit_summary[c(
    "dataset", "standardization", "permutations", "mantel",
    "variation_partition", "graphics", "checks_total", "checks_failed"
  )]
} else {
  data.frame(
    Audit = "Optional repository-level verification",
    Status = "NOT RUN",
    NextStep = "Run the Rscript command shown immediately above"
  )
}
