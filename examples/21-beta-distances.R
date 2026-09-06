# Beta 距离度量：Bray/Jaccard/UniFrac/Aitchison/Robust Aitchison（DEICODE）
# Run sequentially in a new working directory.
# Required packages: ape, digest, dplyr, ggplot2, jsonlite, phyloseq, ragg, readr, svglite, tidyr, vegan.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "env/gemelli.yml",
  "scripts/run_gemelli_rpca.py",
  "data/small/beta-distances/metadata.tsv",
  "data/small/beta-distances/otutab.tsv",
  "data/small/beta-distances/rooted-tree.nwk.gz",
  "data/small/beta-distances/source-summary.json",
  "data/small/beta-distances/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ape)
library(digest)
library(dplyr)
library(ggplot2)
library(jsonlite)
library(phyloseq)
library(readr)
library(tidyr)
library(vegan)

set.seed(20260722)

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
  plot,
  file_base,
  width = 89,
  height = 70,
  units = "mm",
  dpi = 600,
  write_svg = TRUE,
  write_tiff = TRUE,
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
    na = character()
  )
  ids <- as.character(x[[1L]])
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

input_paths <- c(
  otutab = "data/small/beta-distances/otutab.tsv",
  taxonomy = "data/small/beta-distances/taxonomy.tsv",
  metadata = "data/small/beta-distances/metadata.tsv",
  tree = "data/small/beta-distances/rooted-tree.nwk.gz"
)
expected_sha256 <- c(
  otutab = "727073a3f26aff6730ad54b55c11e1d52bd5d3d34f99dbf75d125d2e84d9722d",
  taxonomy = "1189e6e2a02474a665283f6b9962ec25716b6c7faf66011a6506e87a69ddf0ff",
  metadata = "c708b3b89edbb26a4aab26b5e575cdefef1eae28b9b7f93cea5f8f2edfc94221",
  tree = "535f512b3018da462a416b8f04ef36369bf1c1e383510494a790508361fc8c15"
)
observed_sha256 <- vapply(
  input_paths,
  function(path) digest::digest(
    path,
    algo = "sha256",
    serialize = FALSE,
    file = TRUE
  ),
  character(1)
)
stopifnot(identical(observed_sha256, expected_sha256))

otutab <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])
tree <- ape::read.tree(gzfile(input_paths[["tree"]]))

otutab <- data.matrix(otutab)
storage.mode(otutab) <- "numeric"
feature_ids <- rownames(otutab)
sample_ids <- colnames(otutab)

stopifnot(
  nrow(otutab) == 18988L,
  ncol(otutab) == 26L,
  sum(otutab) == 28216678,
  setequal(feature_ids, rownames(taxonomy)),
  setequal(sample_ids, rownames(metadata)),
  all(is.finite(otutab)),
  all(otutab >= 0),
  all(abs(otutab - round(otutab)) < .Machine$double.eps^0.5),
  all(rowSums(otutab) > 0),
  all(colSums(otutab) > 0)
)
taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]

dim(otutab)
sum(otutab)
otutab[1:5, 1:5]
taxonomy[1:5, , drop = FALSE]
head(metadata)

tree_audit <- data.frame(
  Metric = c(
    "Rooted",
    "Tips",
    "Unique tips",
    "Internal nodes",
    "Branch lengths present",
    "Finite branch lengths",
    "Negative branch lengths",
    "Tip set matches otutab"
  ),
  Value = c(
    ape::is.rooted(tree),
    length(tree$tip.label),
    length(unique(tree$tip.label)),
    tree$Nnode,
    !is.null(tree$edge.length),
    all(is.finite(tree$edge.length)),
    sum(tree$edge.length < 0),
    setequal(tree$tip.label, feature_ids)
  )
)
stopifnot(
  ape::is.rooted(tree),
  length(tree$tip.label) == 18988L,
  !anyDuplicated(tree$tip.label),
  !is.null(tree$edge.length),
  all(is.finite(tree$edge.length)),
  all(tree$edge.length >= 0),
  setequal(tree$tip.label, feature_ids)
)
tree_audit

community <- t(otutab)
library_size <- rowSums(community)

stopifnot(min(library_size) >= 50000)
set.seed(20260722)
community_rarefied <- vegan::rrarefy(
  community,
  sample = 50000
)
stopifnot(
  nrow(community_rarefied) == 26L,
  all(rowSums(community_rarefied) == 50000L),
  sum(community_rarefied) == 1300000
)

rarefaction_audit <- data.frame(
  SampleID = rownames(community_rarefied),
  HabitatClass = metadata[rownames(community_rarefied), "HabitatClass"],
  OriginalReads = as.numeric(library_size[rownames(community_rarefied)]),
  RarefiedReads = rowSums(community_rarefied),
  FeaturesBefore = rowSums(community > 0),
  FeaturesAfter = rowSums(community_rarefied > 0)
)
head(rarefaction_audit)
c(
  samples = nrow(community_rarefied),
  reads = sum(community_rarefied),
  retained_features = sum(colSums(community_rarefied) > 0)
)

bray <- vegan::vegdist(community_rarefied, method = "bray")
jaccard <- vegan::vegdist(
  community_rarefied,
  method = "jaccard",
  binary = TRUE
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

# Bray sensitivity: rarefied counts versus unrarefied relative abundance
bray_relative <- vegan::vegdist(
  vegan::decostand(community, method = "total"),
  method = "bray"
)
bray_sensitivity_rho <- cor(
  as.vector(bray),
  as.vector(bray_relative),
  method = "spearman"
)
bray_sensitivity_rho

feature_total <- rowSums(otutab)
feature_prevalence <- rowMeans(otutab > 0)
keep_composition <- feature_total > 10 & feature_prevalence > 0.10

feature_filter_audit <- data.frame(
  FeatureID = feature_ids,
  TotalCount = feature_total,
  Prevalence = feature_prevalence,
  Retained = keep_composition
)
stopifnot(
  sum(keep_composition) == 11063L,
  sum(otutab[keep_composition, , drop = FALSE]) == 28148221
)

aitchison_distance <- function(pseudocount) {
  log_counts <- log(
    t(otutab[keep_composition, , drop = FALSE]) + pseudocount
  )
  clr <- log_counts - rowMeans(log_counts)
  dist(clr, method = "euclidean")
}

aitchison_pc01 <- aitchison_distance(0.1)
aitchison_pc05 <- aitchison_distance(0.5)
aitchison_pc10 <- aitchison_distance(1.0)

c(
  retained_features = sum(keep_composition),
  retained_reads = sum(otutab[keep_composition, , drop = FALSE])
)

gemelli_python <- Sys.getenv(
  "GEMELLI_PYTHON",
  "gemelli-env/bin/python"
)
gemelli_runner <- Sys.getenv(
  "GEMELLI_RUNNER",
  "scripts/run_gemelli_rpca.py"
)
stopifnot(file.exists(gemelli_python), file.exists(gemelli_runner))

dir.create("results/21-beta-distances", recursive = TRUE, showWarnings = FALSE)
gemelli_log <- system2(
  gemelli_python,
  args = c(
    shQuote(gemelli_runner),
    "--input", shQuote(input_paths[["otutab"]]),
    "--output-dir", shQuote("results/21-beta-distances"),
    "--n-components", "3",
    "--min-sample-count", "500",
    "--min-feature-count", "10",
    "--min-feature-frequency", "10",
    "--max-iterations", "5"
  ),
  stdout = TRUE,
  stderr = TRUE,
  env = c(
    "PYTHONHASHSEED=20260722",
    "OMP_NUM_THREADS=1",
    "OPENBLAS_NUM_THREADS=1",
    "MKL_NUM_THREADS=1"
  )
)
gemelli_status <- attr(gemelli_log, "status")
if (is.null(gemelli_status)) gemelli_status <- 0L
stopifnot(gemelli_status == 0L)
gemelli_log

robust_frame <- readr::read_tsv(
  "results/21-beta-distances/robust-aitchison-distance.tsv",
  show_col_types = FALSE,
  name_repair = "minimal"
)
robust_ids <- as.character(robust_frame[[1L]])
robust_matrix <- data.matrix(robust_frame[-1L])
rownames(robust_matrix) <- robust_ids
colnames(robust_matrix) <- colnames(robust_frame)[-1L]
stopifnot(
  setequal(sample_ids, rownames(robust_matrix)),
  setequal(sample_ids, colnames(robust_matrix))
)
robust_aitchison <- as.dist(
  robust_matrix[sample_ids, sample_ids, drop = FALSE]
)

align_distance <- function(distance, ids = sample_ids) {
  matrix <- as.matrix(distance)
  as.dist(matrix[ids, ids, drop = FALSE])
}

distances <- list(
  "Bray-Curtis" = align_distance(bray),
  "Jaccard" = align_distance(jaccard),
  "Unweighted UniFrac" = align_distance(unweighted_unifrac),
  "Weighted UniFrac" = align_distance(weighted_unifrac),
  "Aitchison" = align_distance(aitchison_pc05),
  "Robust Aitchison" = align_distance(robust_aitchison)
)

validate_distance <- function(distance) {
  matrix <- as.matrix(distance)
  stopifnot(
    identical(dim(matrix), c(26L, 26L)),
    identical(rownames(matrix), sample_ids),
    identical(colnames(matrix), sample_ids),
    all(is.finite(matrix)),
    min(matrix) >= -1e-12,
    max(abs(matrix - t(matrix))) <= 1e-10,
    max(abs(diag(matrix))) <= 1e-10
  )
  TRUE
}
vapply(distances, validate_distance, logical(1))

pair_indices <- which(upper.tri(matrix(0, 26, 26)), arr.ind = TRUE)
distance_summary <- bind_rows(lapply(names(distances), function(metric) {
  matrix <- as.matrix(distances[[metric]])
  values <- matrix[pair_indices]
  sample_a <- rownames(matrix)[pair_indices[, 1L]]
  sample_b <- colnames(matrix)[pair_indices[, 2L]]
  pair_type <- ifelse(
    metadata[sample_a, "HabitatClass"] == metadata[sample_b, "HabitatClass"],
    "Within habitat",
    "Between habitats"
  )
  data.frame(
    Metric = metric,
    Minimum = min(values),
    Median = median(values),
    Maximum = max(values),
    WithinMedian = median(values[pair_type == "Within habitat"]),
    BetweenMedian = median(values[pair_type == "Between habitats"]),
    BetweenWithinRatio =
      median(values[pair_type == "Between habitats"]) /
      median(values[pair_type == "Within habitat"])
  )
}))
distance_summary

metric_levels <- names(distances)
distance_vectors <- do.call(cbind, lapply(distances, as.vector))
distance_cor_matrix <- cor(distance_vectors, method = "spearman")

distance_correlation <- tidyr::expand_grid(
  MetricX = rownames(distance_cor_matrix),
  MetricY = colnames(distance_cor_matrix)
)
distance_correlation$Spearman <- mapply(
  function(x, y) distance_cor_matrix[x, y],
  distance_correlation$MetricX,
  distance_correlation$MetricY
)

round(distance_cor_matrix, 3)

pseudocount_distances <- list(
  "0.1" = align_distance(aitchison_pc01),
  "0.5" = align_distance(aitchison_pc05),
  "1.0" = align_distance(aitchison_pc10)
)
aitchison_pseudocount_sensitivity <- data.frame(
  Pseudocount = as.numeric(names(pseudocount_distances)),
  SpearmanToAitchisonPC05 = vapply(
    pseudocount_distances,
    function(x) cor(
      as.vector(x),
      as.vector(aitchison_pc05),
      method = "spearman"
    ),
    numeric(1)
  ),
  SpearmanToRobustAitchison = vapply(
    pseudocount_distances,
    function(x) cor(
      as.vector(x),
      as.vector(robust_aitchison),
      method = "spearman"
    ),
    numeric(1)
  )
)
aitchison_pseudocount_sensitivity

pcoa_score_rows <- list()
pcoa_audit_rows <- list()

for (metric in names(distances)) {
  distance <- distances[[metric]]
  raw_pcoa <- suppressWarnings(
    cmdscale(
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
  correction <- if (negative_count > 0L) "Lingoes" else "None"

  fit <- ape::pcoa(
    distance,
    correction = if (negative_count > 0L) "lingoes" else "none"
  )
  if (negative_count > 0L) {
    coordinates <- fit$vectors.cor[, 1:2, drop = FALSE]
    relative_eigen <- fit$values$Rel_corr_eig[1:2]
  } else {
    coordinates <- fit$vectors[, 1:2, drop = FALSE]
    relative_eigen <- fit$values$Relative_eig[1:2]
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
    Correction = correction
  )
  pcoa_audit_rows[[metric]] <- data.frame(
    Metric = metric,
    PositiveEigenvalues = sum(eigenvalues > 1e-10),
    NegativeEigenvalues = negative_count,
    MinimumEigenvalue = min(eigenvalues),
    NegativePositiveRatio = negative_sum / positive_sum,
    Correction = correction,
    Axis1Explained = relative_eigen[[1L]],
    Axis2Explained = relative_eigen[[2L]]
  )
}

pcoa_scores <- bind_rows(pcoa_score_rows)
pcoa_eigen_audit <- bind_rows(pcoa_audit_rows)
stopifnot(
  nrow(pcoa_scores) == 156L,
  pcoa_eigen_audit$NegativeEigenvalues[
    pcoa_eigen_audit$Metric == "Weighted UniFrac"
  ] == 3L,
  all(
    pcoa_eigen_audit$Correction[
      pcoa_eigen_audit$Metric != "Weighted UniFrac"
    ] == "None"
  )
)
pcoa_eigen_audit

dir.create("results/21-beta-distances", recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(
  distance_summary,
  "results/21-beta-distances/distance-summary.tsv"
)
readr::write_tsv(
  distance_correlation,
  "results/21-beta-distances/distance-correlation.tsv"
)
readr::write_tsv(
  pcoa_eigen_audit,
  "results/21-beta-distances/pcoa-eigen-audit.tsv"
)
readr::write_tsv(
  pcoa_scores,
  "results/21-beta-distances/pcoa-scores.tsv"
)
readr::write_tsv(
  aitchison_pseudocount_sensitivity,
  "results/21-beta-distances/aitchison-pseudocount-sensitivity.tsv"
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

choice_plot <- ggplot(
  choice_grid,
  aes(x = Signal, y = Metric, fill = Included)
) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = Label), fontface = "bold", size = 3.2) +
  scale_fill_manual(
    values = c(`FALSE` = "#ECEFF3", `TRUE` = pal_pub[["sky"]]),
    guide = "none"
  ) +
  labs(
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
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold")
  )

save_pub(
  choice_plot,
  "figures/21-distance-choice-map",
  width = 183,
  height = 92
)
choice_plot

correlation_plot_data <- distance_correlation
correlation_plot_data$MetricX <- factor(
  correlation_plot_data$MetricX,
  levels = metric_levels
)
correlation_plot_data$MetricY <- factor(
  correlation_plot_data$MetricY,
  levels = rev(metric_levels)
)

correlation_plot <- ggplot(
  correlation_plot_data,
  aes(x = MetricX, y = MetricY, fill = Spearman)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Spearman)), size = 2.8) +
  scale_fill_gradientn(
    colours = c("#F7F7F7", "#B8DCEB", "#0072B2"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  labs(
    title = "Metric choice changes the rank order of sample pairs",
    subtitle = "Spearman correlation across all 325 pairwise distances",
    x = NULL,
    y = NULL,
    fill = "Spearman\nrho",
    caption = "High correlation in this dataset does not make two metric definitions interchangeable."
  ) +
  theme_pub(base_size = 8.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "right"
  )

save_pub(
  correlation_plot,
  "figures/21-distance-correlation",
  width = 183,
  height = 112
)
correlation_plot

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

pcoa_labels <- setNames(
  vapply(metric_levels, function(metric) {
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
  }, character(1)),
  metric_levels
)

pcoa_plot_data <- pcoa_scores
pcoa_plot_data$Metric <- factor(pcoa_plot_data$Metric, levels = metric_levels)
pcoa_plot_data$HabitatClass <- factor(
  pcoa_plot_data$HabitatClass,
  levels = habitat_levels
)

pcoa_plot <- ggplot(
  pcoa_plot_data,
  aes(
    x = Axis1,
    y = Axis2,
    fill = HabitatClass,
    shape = HabitatClass
  )
) +
  geom_hline(yintercept = 0, colour = "#D1D5DB", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#D1D5DB", linewidth = 0.3) +
  geom_point(size = 2.2, stroke = 0.45, colour = "#1F2937", alpha = 0.92) +
  facet_wrap(
    ~Metric,
    scales = "free",
    ncol = 3,
    labeller = as_labeller(pcoa_labels)
  ) +
  scale_fill_manual(values = habitat_palette, drop = FALSE) +
  scale_shape_manual(values = habitat_shapes, drop = FALSE) +
  labs(
    title = "The same samples occupy six different distance geometries",
    subtitle = "GlobalPatterns 16S communities; each point is one sample",
    x = "Principal coordinate 1",
    y = "Principal coordinate 2",
    fill = "Habitat",
    shape = "Habitat",
    caption = "Coordinates visualize distance structure; they are not tests of group separation."
  ) +
  theme_pub(base_size = 8.2) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(size = 7.1, lineheight = 0.95)
  )

save_pub(
  pcoa_plot,
  "figures/21-pcoa-metric-comparison",
  width = 183,
  height = 132
)
pcoa_plot

sample_a <- sample_ids[pair_indices[, 1L]]
sample_b <- sample_ids[pair_indices[, 2L]]
pair_type <- ifelse(
  metadata[sample_a, "HabitatClass"] == metadata[sample_b, "HabitatClass"],
  "Within habitat",
  "Between habitats"
)

zero_plot_data <- bind_rows(lapply(names(pseudocount_distances), function(pc) {
  data.frame(
    Pseudocount = pc,
    ClassicAitchison = as.vector(pseudocount_distances[[pc]]),
    RobustAitchison = as.vector(robust_aitchison),
    PairType = pair_type
  )
}))
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
  )
)
zero_annotation$Pseudocount <- factor(
  zero_annotation$Pseudocount,
  levels = levels(zero_plot_data$Pseudocount)
)

zero_plot <- ggplot(
  zero_plot_data,
  aes(
    x = ClassicAitchison,
    y = RobustAitchison,
    colour = PairType
  )
) +
  geom_point(size = 1.25, alpha = 0.62) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.55,
    colour = "#374151",
    show.legend = FALSE
  ) +
  geom_text(
    data = zero_annotation,
    aes(x = -Inf, y = Inf, label = Label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.35,
    size = 2.7,
    fontface = "bold"
  ) +
  facet_wrap(~Pseudocount, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c(
    "Within habitat" = pal_pub[["blue"]],
    "Between habitats" = pal_pub[["vermillion"]]
  )) +
  labs(
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
  theme(legend.position = "top")

save_pub(
  zero_plot,
  "figures/21-aitchison-zero-sensitivity",
  width = 183,
  height = 94
)
zero_plot
