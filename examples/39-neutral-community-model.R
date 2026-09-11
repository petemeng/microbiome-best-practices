# 中性群落模型：Sloan 拟合、迁入率与偏离类群
# Run sequentially in a new working directory.
# Required packages: ggplot2, ggrepel, knitr, patchwork, ragg, scales, svglite, vegan.

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
library(vegan)

set.seed(20260739)
font_pub <- "sans"
pal_pub <- c(
  blue = "#0072B2", orange = "#E69F00", green = "#009E73",
  vermillion = "#D55E00", purple = "#CC79A7", sky = "#56B4E9",
  yellow = "#F0E442", grey = "#6B7280"
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
      axis.ticks = ggplot2::element_line(colour = "#1A1A1A", linewidth = 0.3),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.15)),
      plot.subtitle = ggplot2::element_text(colour = "#4D4D4D"),
      plot.caption = ggplot2::element_text(colour = "#666666", hjust = 0),
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2", colour = "#B3B3B3", linewidth = 0.3
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(), legend.position = "top"
    )
}
save_pub <- function(
  plot, file_base, width = 89, height = 70, units = "mm",
  dpi = 600, write_svg = TRUE, write_tiff = TRUE,
  base_family = font_pub
) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(file_base, ".pdf"), plot, width = width, height = height,
    units = units, device = grDevices::cairo_pdf,
    family = base_family, bg = "white"
  )
  if (write_svg) ggplot2::ggsave(
    paste0(file_base, ".svg"), plot, width = width, height = height,
    units = units, device = svglite::svglite, bg = "white"
  )
  ggplot2::ggsave(
    paste0(file_base, ".png"), plot, width = width, height = height,
    units = units, dpi = dpi, device = ragg::agg_png, bg = "white"
  )
  if (write_tiff) ggplot2::ggsave(
    paste0(file_base, ".tiff"), plot, width = width, height = height,
    units = units, dpi = dpi, device = ragg::agg_tiff,
    compression = "lzw", bg = "white"
  )
  invisible(plot)
}

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path, show_col_types = FALSE, progress = FALSE,
    name_repair = "minimal", na = character()
  )
  ids <- as.character(x[[1L]])
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}
otutab_raw <- read_keyed_tsv("data/small/otutab.tsv")
taxonomy <- read_keyed_tsv("data/small/taxonomy.tsv")
metadata <- read_keyed_tsv("data/small/metadata.tsv")
otutab <- as.matrix(data.frame(
  lapply(otutab_raw, as.integer), check.names = FALSE
))
rownames(otutab) <- rownames(otutab_raw)
metadata$Group <- factor(metadata$Group, levels = c("CW", "IW", "TW"))
otutab <- otutab[, rownames(metadata), drop = FALSE]

clean_taxon <- function(x) {
  x <- sub("^[A-Za-z]__", "", trimws(as.character(x)))
  x[
    is.na(x) | x == "" |
      grepl("unclassified|unknown|unassigned|uncultured|metagenome", x,
            ignore.case = TRUE)
  ] <- NA_character_
  x
}
rank_names <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
)
taxonomy_clean <- as.data.frame(
  lapply(taxonomy[rownames(otutab), rank_names], clean_taxon),
  check.names = FALSE
)
rownames(taxonomy_clean) <- rownames(otutab)
display_taxon <- apply(
  taxonomy_clean[, rev(rank_names), drop = FALSE], 1L,
  function(x) {
    value <- x[!is.na(x)][1L]
    if (length(value)) value else "Unclassified OTU"
  }
)
names(display_taxon) <- rownames(taxonomy_clean)

stopifnot(
  nrow(otutab) == 13628L,
  ncol(otutab) == 90L,
  sum(otutab) == 1619670L,
  identical(colnames(otutab), rownames(metadata)),
  all(colSums(otutab) >= 10000L)
)

head(otutab[, 1:5])
head(taxonomy[, rank_names])
head(metadata)

fit_ncm <- function(counts, detection = NULL) {
  stopifnot(all(counts >= 0), all(colSums(counts) > 0))
  relative <- sweep(counts, 2L, colSums(counts), "/")
  mean_abundance <- rowMeans(relative)
  mean_library <- mean(colSums(counts))
  if (is.null(detection)) detection <- 1 / mean_library
  occupancy <- rowMeans(relative >= detection)
  estimable <- mean_abundance > 0 & mean_abundance < 1

  objective <- function(log10_m) {
    migration <- 10^log10_m
    predicted <- stats::pbeta(
      detection,
      mean_library * migration * mean_abundance[estimable],
      mean_library * migration * (1 - mean_abundance[estimable]),
      lower.tail = FALSE
    )
    sum((occupancy[estimable] - predicted)^2)
  }
  optimum <- stats::optimize(objective, interval = c(-8, 0))
  migration <- 10^optimum$minimum
  predicted <- rep(NA_real_, length(mean_abundance))
  predicted[estimable] <- stats::pbeta(
    detection,
    mean_library * migration * mean_abundance[estimable],
    mean_library * migration * (1 - mean_abundance[estimable]),
    lower.tail = FALSE
  )
  predicted[mean_abundance <= 0] <- 0
  predicted[mean_abundance >= 1] <- 1
  sample_n <- ncol(counts)
  lower <- stats::qbinom(0.025, sample_n, predicted) / sample_n
  upper <- stats::qbinom(0.975, sample_n, predicted) / sample_n
  classification <- ifelse(
    occupancy > upper, "Above neutral",
    ifelse(occupancy < lower, "Below neutral", "Neutral envelope")
  )
  r_squared <- 1 -
    sum((occupancy[estimable] - predicted[estimable])^2) /
    sum((occupancy[estimable] - mean(occupancy[estimable]))^2)

  data.frame(
    FeatureID = rownames(counts),
    MeanAbundance = mean_abundance,
    Occupancy = occupancy,
    Predicted = predicted,
    Lower = lower,
    Upper = upper,
    Classification = classification,
    Migration = migration,
    MeanLibrary = mean_library,
    Detection = detection,
    RSquared = r_squared,
    stringsAsFactors = FALSE
  )
}

pooled_fit <- fit_ncm(otutab)
pooled_fit$DisplayTaxon <- display_taxon[pooled_fit$FeatureID]
pooled_fit$StandardizedDeviation <- with(
  pooled_fit,
  (Occupancy - Predicted) /
    pmax(sqrt(Predicted * (1 - Predicted) / ncol(otutab)), 1 / ncol(otutab))
)
primary_summary <- unique(pooled_fit[
  , c("Migration", "MeanLibrary", "Detection", "RSquared")
])
classification_summary <- as.data.frame(
  table(pooled_fit$Classification), stringsAsFactors = FALSE
)
names(classification_summary) <- c("Classification", "OTUs")

stopifnot(
  nrow(primary_summary) == 1L,
  primary_summary$Migration > 0,
  primary_summary$Migration < 1,
  is.finite(primary_summary$RSquared)
)
primary_summary
classification_summary

bootstrap_migration <- function(counts, replicates, seed) {
  set.seed(seed)
  replicate(replicates, {
    columns <- sample(seq_len(ncol(counts)), replace = TRUE)
    unique(fit_ncm(counts[, columns, drop = FALSE])$Migration)
  })
}
summarise_pool <- function(counts, pool, seed) {
  fit <- fit_ncm(counts)
  migration_boot <- bootstrap_migration(
    counts, replicates = 200L, seed = seed
  )
  data.frame(
    Pool = pool,
    Samples = ncol(counts),
    Migration = unique(fit$Migration),
    Lower = unname(stats::quantile(migration_boot, 0.025)),
    Upper = unname(stats::quantile(migration_boot, 0.975)),
    RSquared = unique(fit$RSquared),
    Detection = unique(fit$Detection),
    stringsAsFactors = FALSE
  )
}

pool_ledger <- rbind(
  summarise_pool(otutab, "All wetlands", 20260739L),
  summarise_pool(
    otutab[, metadata$Group == "CW", drop = FALSE], "CW", 20260740L
  ),
  summarise_pool(
    otutab[, metadata$Group == "IW", drop = FALSE], "IW", 20260741L
  ),
  summarise_pool(
    otutab[, metadata$Group == "TW", drop = FALSE], "TW", 20260742L
  )
)
pool_ledger

fixed_detection_fit <- fit_ncm(otutab, detection = 1e-4)
set.seed(20260739)
rarefied_otutab <- t(vegan::rrarefy(t(otutab), sample = 10000L))
rarefied_fit <- fit_ncm(rarefied_otutab)

classification_sensitivity <- rbind(
  transform(
    as.data.frame(table(pooled_fit$Classification), stringsAsFactors = FALSE),
    Model = "Primary: 1 / mean depth"
  ),
  transform(
    as.data.frame(
      table(fixed_detection_fit$Classification), stringsAsFactors = FALSE
    ),
    Model = "Fixed relative threshold"
  ),
  transform(
    as.data.frame(table(rarefied_fit$Classification), stringsAsFactors = FALSE),
    Model = "Rarefied to 10k reads"
  )
)
names(classification_sensitivity)[1:2] <- c("Classification", "OTUs")
classification_sensitivity$Fraction <- ave(
  classification_sensitivity$OTUs,
  classification_sensitivity$Model,
  FUN = function(x) x / sum(x)
)

classification_agreement <- data.frame(
  Comparison = c(
    "Primary vs fixed threshold", "Primary vs rarefied"
  ),
  Agreement = c(
    mean(pooled_fit$Classification == fixed_detection_fit$Classification),
    mean(pooled_fit$Classification == rarefied_fit$Classification)
  ),
  stringsAsFactors = FALSE
)

detection_grid <- 10^seq(-5, -3, length.out = 17L)
detection_sensitivity <- do.call(rbind, lapply(detection_grid, function(value) {
  fit <- fit_ncm(otutab, detection = value)
  data.frame(
    Detection = value,
    Migration = unique(fit$Migration),
    RSquared = unique(fit$RSquared),
    NeutralFraction = mean(fit$Classification == "Neutral envelope"),
    stringsAsFactors = FALSE
  )
}))

stopifnot(
  all(colSums(rarefied_otutab) == 10000L),
  nrow(classification_sensitivity) == 9L,
  nrow(detection_sensitivity) == 17L,
  all(classification_agreement$Agreement >= 0 &
        classification_agreement$Agreement <= 1)
)
classification_agreement

neutral_model_audit <- data.frame(
  Dimension = c(
    "Regional pool", "Detection", "Sequencing depth", "Model fit",
    "Envelope", "Taxonomic rank", "Group contrast"
  ),
  FragilePractice = c(
    "Pool all habitats without justification",
    "Define presence differently across analyses",
    "Ignore library-size-driven occupancy",
    "Translate R-squared into neutral-process percentage",
    "Call every in-envelope taxon neutral",
    "Compare ASV and genus parameters directly",
    "Treat group-specific m as a randomized effect"
  ),
  UpgradedContract = c(
    "State source pool and refit plausible strata",
    "Prespecify threshold and show sensitivity",
    "Report N and a depth-handling sensitivity",
    "Use R-squared only as abundance-occupancy fit",
    "Use predictive labels as diagnostic evidence",
    "Keep one feature universe per comparison",
    "Report pool, depth and confounding boundaries"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(
  neutral_model_audit,
  caption = "Audit contract for neutral community-model claims"
)

class_colours <- c(
  "Above neutral" = pal_pub[["vermillion"]],
  "Neutral envelope" = pal_pub[["grey"]],
  "Below neutral" = pal_pub[["blue"]]
)
curve_data <- pooled_fit[order(pooled_fit$MeanAbundance), ]
label_data <- pooled_fit[
  order(-abs(pooled_fit$StandardizedDeviation)),
]
label_data <- label_data[!duplicated(label_data$DisplayTaxon), ]
label_data <- head(label_data, 8L)
neutral_curve_plot <- ggplot() +
  geom_ribbon(
    data = curve_data,
    aes(x = MeanAbundance, ymin = Lower, ymax = Upper),
    fill = "#BFC5CC", alpha = 0.45
  ) +
  geom_line(
    data = curve_data,
    aes(x = MeanAbundance, y = Predicted),
    colour = "#202020", linewidth = 0.65
  ) +
  geom_point(
    data = pooled_fit,
    aes(x = MeanAbundance, y = Occupancy, colour = Classification),
    size = 0.8, alpha = 0.42
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(x = MeanAbundance, y = Occupancy, label = DisplayTaxon),
    family = font_pub, size = 2.35, max.overlaps = Inf,
    min.segment.length = 0, box.padding = 0.45, point.padding = 0.15,
    force = 1.5, max.time = 2, show.legend = FALSE
  ) +
  scale_x_log10(labels = scales::label_scientific()) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
  scale_colour_manual(values = class_colours) +
  labs(
    title = "Neutral model reproduces the pooled abundance-occupancy pattern",
    subtitle = sprintf(
      "m = %.3f; R² = %.3f; N = %.0f reads",
      primary_summary$Migration, primary_summary$RSquared,
      primary_summary$MeanLibrary
    ),
    x = "Mean relative abundance in the metacommunity",
    y = "Sample occupancy", colour = "Predictive class",
    caption = "Ribbon: exact-binomial 95% predictive envelope; labels mark the largest standardized deviations."
  ) +
  theme_pub(base_size = 8.4) +
  theme(legend.position = "bottom")
save_pub(
  neutral_curve_plot, "figures/39-neutral-abundance-occupancy",
  width = 180, height = 143
)
neutral_curve_plot

pool_ledger$Pool <- factor(
  pool_ledger$Pool, levels = rev(c("All wetlands", "CW", "IW", "TW"))
)
pool_plot <- ggplot(
  pool_ledger,
  aes(x = Migration, y = Pool, colour = Pool)
) +
  geom_errorbarh(
    aes(xmin = Lower, xmax = Upper), height = 0.15, linewidth = 0.7
  ) +
  geom_point(aes(size = RSquared), shape = 16) +
  scale_x_log10() +
  scale_colour_manual(values = c(
    "All wetlands" = pal_pub[["grey"]], "CW" = pal_pub[["blue"]],
    "IW" = pal_pub[["orange"]], "TW" = pal_pub[["green"]]
  )) +
  scale_size_continuous(range = c(3.2, 6.2)) +
  labs(
    title = "Immigration estimates depend on the defined regional pool",
    subtitle = "Points: fitted m; bars: 95% sample-bootstrap interval (200 replicates)",
    x = "Migration parameter m (log scale)", y = NULL,
    colour = "Regional pool", size = "Model R²",
    caption = "Group-specific fits also change sample number, mean depth, composition and spatial extent."
  ) +
  theme_pub(base_size = 8.6) +
  theme(legend.position = "bottom")
save_pub(
  pool_plot, "figures/39-neutral-pool-migration",
  width = 180, height = 122
)
pool_plot

classification_sensitivity$Model <- factor(
  classification_sensitivity$Model,
  levels = c(
    "Primary: 1 / mean depth", "Fixed relative threshold",
    "Rarefied to 10k reads"
  )
)
classification_sensitivity$Classification <- factor(
  classification_sensitivity$Classification,
  levels = c("Below neutral", "Neutral envelope", "Above neutral")
)
class_sensitivity_plot <- ggplot(
  classification_sensitivity,
  aes(x = Model, y = Fraction, fill = Classification)
) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = class_colours, drop = FALSE) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  labs(
    title = "Predictive labels move when observation rules change",
    subtitle = sprintf(
      "Agreement with primary labels: fixed threshold %.1f%%; rarefied %.1f%%",
      100 * classification_agreement$Agreement[1],
      100 * classification_agreement$Agreement[2]
    ),
    x = NULL, y = "Fraction of OTUs", fill = "Predictive class",
    caption = "Rarefaction is used only as a sensitivity analysis, with seed 20260739."
  ) +
  theme_pub(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 18, hjust = 1),
    legend.position = "bottom"
  )
save_pub(
  class_sensitivity_plot, "figures/39-neutral-classification-sensitivity",
  width = 180, height = 122
)
class_sensitivity_plot

m_panel <- ggplot(
  detection_sensitivity, aes(x = Detection, y = Migration)
) +
  geom_vline(
    xintercept = primary_summary$Detection,
    linetype = 2, colour = pal_pub[["grey"]]
  ) +
  geom_line(colour = pal_pub[["blue"]], linewidth = 0.8) +
  geom_point(colour = pal_pub[["blue"]], size = 1.9) +
  scale_x_log10(labels = scales::label_scientific()) +
  scale_y_log10() +
  labs(
    title = "A. Migration estimate", x = "Detection threshold",
    y = "Migration parameter m"
  ) +
  theme_pub(base_size = 8.0)
r2_panel <- ggplot(
  detection_sensitivity, aes(x = Detection, y = RSquared)
) +
  geom_vline(
    xintercept = primary_summary$Detection,
    linetype = 2, colour = pal_pub[["grey"]]
  ) +
  geom_line(colour = pal_pub[["vermillion"]], linewidth = 0.8) +
  geom_point(colour = pal_pub[["vermillion"]], size = 1.9) +
  scale_x_log10(labels = scales::label_scientific()) +
  labs(
    title = "B. Abundance-occupancy fit", x = "Detection threshold",
    y = "Model R²"
  ) +
  theme_pub(base_size = 8.0)
detection_plot <- patchwork::wrap_plots(m_panel, r2_panel, nrow = 1L) +
  patchwork::plot_annotation(
    title = "Observation limits are part of the neutral-model estimand",
    caption = "Dashed lines mark the primary threshold (1 / mean library size)."
  )
save_pub(
  detection_plot, "figures/39-neutral-detection-sensitivity",
  width = 180, height = 122
)
detection_plot

result_dir <- "results/39-neutral-community-model"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  pooled_fit,
  file.path(result_dir, "otu-neutral-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  pool_ledger,
  file.path(result_dir, "pool-migration-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  classification_sensitivity,
  file.path(result_dir, "classification-sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  detection_sensitivity,
  file.path(result_dir, "detection-sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

expected_bases <- c(
  "figures/39-neutral-abundance-occupancy",
  "figures/39-neutral-pool-migration",
  "figures/39-neutral-classification-sensitivity",
  "figures/39-neutral-detection-sensitivity"
)
expected_files <- as.vector(outer(
  expected_bases, c(".pdf", ".svg", ".png", ".tiff"), paste0
))
stopifnot(
  all(file.exists(expected_files)),
  nrow(pooled_fit) == 13628L,
  nrow(pool_ledger) == 4L,
  all(pool_ledger$Migration > 0 & pool_ledger$Migration < 1),
  unique(pooled_fit$RSquared) > 0.50,
  all(colSums(rarefied_otutab) == 10000L),
  setequal(
    unique(pooled_fit$Classification),
    c("Above neutral", "Below neutral", "Neutral envelope")
  )
)
