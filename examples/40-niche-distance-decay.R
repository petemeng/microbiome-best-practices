# 生态位宽度、物种丰度分布与距离衰减
# Run sequentially in a new working directory.
# Required packages: ggplot2, ggrepel, knitr, patchwork, permute, ragg, scales, svglite, vegan.

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

set.seed(20260740)
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
environment <- read_keyed_tsv("data/small/environment.tsv")
otutab <- as.matrix(data.frame(
  lapply(otutab_raw, as.integer), check.names = FALSE
))
rownames(otutab) <- rownames(otutab_raw)
metadata$Group <- factor(metadata$Group, levels = c("CW", "IW", "TW"))
otutab <- otutab[, rownames(metadata), drop = FALSE]
environment <- environment[rownames(metadata), , drop = FALSE]

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
  identical(rownames(environment), rownames(metadata)),
  all(stats::complete.cases(environment[, c("Latitude", "Longitude")]))
)

head(otutab[, 1:5])
head(taxonomy[, rank_names])
head(metadata)
head(environment[, c("Latitude", "Longitude")])

relative <- sweep(otutab, 2L, colSums(otutab), "/")
prevalence <- rowMeans(otutab > 0)
total_reads <- rowSums(otutab)
breadth_keep <- prevalence >= 0.05 & total_reads >= 100L
breadth_relative <- relative[breadth_keep, , drop = FALSE]

sample_share <- sweep(
  breadth_relative, 1L, rowSums(breadth_relative), "/"
)
levins_b <- 1 / rowSums(sample_share^2)
standardized_b <- (levins_b - 1) / (ncol(otutab) - 1)
breadth_quantiles <- stats::quantile(
  standardized_b, c(0.10, 0.90), names = FALSE
)
breadth_ledger <- data.frame(
  FeatureID = names(standardized_b),
  DisplayTaxon = display_taxon[names(standardized_b)],
  Prevalence = prevalence[names(standardized_b)],
  TotalReads = total_reads[names(standardized_b)],
  LevinsB = unname(levins_b),
  StandardizedBreadth = unname(standardized_b),
  BreadthClass = ifelse(
    standardized_b <= breadth_quantiles[[1L]], "Narrow decile",
    ifelse(
      standardized_b >= breadth_quantiles[[2L]],
      "Broad decile", "Intermediate"
    )
  ),
  stringsAsFactors = FALSE
)
breadth_occupancy_rho <- stats::cor(
  breadth_ledger$StandardizedBreadth,
  breadth_ledger$Prevalence,
  method = "spearman"
)

stopifnot(
  nrow(breadth_ledger) == 3208L,
  all(breadth_ledger$StandardizedBreadth >= 0 &
        breadth_ledger$StandardizedBreadth <= 1)
)
summary(breadth_ledger$StandardizedBreadth)

rad_models <- c("Null", "Preemption", "Lognormal", "Zipf", "Mandelbrot")
rad_aic <- matrix(
  NA_real_, nrow = ncol(otutab), ncol = length(rad_models),
  dimnames = list(colnames(otutab), rad_models)
)
rad_warning_ledger <- data.frame(
  SampleID = character(), Message = character(), stringsAsFactors = FALSE
)

for (sample_id in colnames(otutab)) {
  abundance <- otutab[, sample_id]
  abundance <- abundance[abundance > 0]
  fit <- withCallingHandlers(
    vegan::radfit(abundance),
    warning = function(warning_condition) {
      warning_message <- conditionMessage(warning_condition)
      warning_message <- if (
        grepl("did not converge|算法没有聚合", warning_message)
      ) {
        "GLM did not converge"
      } else {
        warning_message
      }
      rad_warning_ledger <<- rbind(
        rad_warning_ledger,
        data.frame(
          SampleID = sample_id,
          Message = warning_message,
          stringsAsFactors = FALSE
        )
      )
      invokeRestart("muffleWarning")
    }
  )
  rad_aic[sample_id, ] <- vapply(
    fit$models, stats::AIC, numeric(1)
  )[rad_models]
}
rad_winner <- rad_models[max.col(-rad_aic, ties.method = "first")]
rad_delta_aic <- sweep(rad_aic, 1L, apply(rad_aic, 1L, min), "-")
rad_aic_ledger <- data.frame(
  SampleID = rep(rownames(rad_aic), times = ncol(rad_aic)),
  Group = rep(as.character(metadata[rownames(rad_aic), "Group"]),
              times = ncol(rad_aic)),
  Model = rep(colnames(rad_aic), each = nrow(rad_aic)),
  AIC = as.vector(rad_aic),
  DeltaAIC = as.vector(rad_delta_aic),
  stringsAsFactors = FALSE
)
rad_model_summary <- data.frame(
  Model = rad_models,
  WinningSamples = as.integer(table(factor(rad_winner, levels = rad_models))),
  MedianDeltaAIC = apply(rad_delta_aic, 2L, stats::median),
  stringsAsFactors = FALSE
)

sample_richness <- colSums(otutab > 0)
rad_warning_samples <- unique(rad_warning_ledger$SampleID)
representative_candidates <- setdiff(
  names(sample_richness), rad_warning_samples
)
representative_sample <- representative_candidates[
  which.min(abs(
    sample_richness[representative_candidates] -
      stats::median(sample_richness[representative_candidates])
  ))
]
representative_abundance <- otutab[, representative_sample]
representative_abundance <- representative_abundance[
  representative_abundance > 0
]
representative_fit <- suppressWarnings(
  vegan::radfit(representative_abundance)
)
representative_rad <- do.call(rbind, lapply(rad_models, function(model) {
  data.frame(
    Rank = seq_along(representative_fit$models[[model]]$fitted.values),
    Abundance = as.numeric(
      representative_fit$models[[model]]$fitted.values
    ),
    Model = model,
    stringsAsFactors = FALSE
  )
}))
representative_observed <- data.frame(
  Rank = seq_along(representative_fit$y),
  Abundance = as.numeric(representative_fit$y)
)

stopifnot(
  nrow(rad_aic_ledger) == 450L,
  sum(rad_model_summary$WinningSamples) == 90L,
  all(is.finite(rad_aic_ledger$AIC)),
  representative_sample %in% rownames(metadata)
)
rad_model_summary
rad_warning_ledger

haversine_km <- function(lat1, lon1, lat2, lon2) {
  earth_radius <- 6371.0088
  phi1 <- lat1 * pi / 180
  phi2 <- lat2 * pi / 180
  delta_phi <- (lat2 - lat1) * pi / 180
  delta_lambda <- (lon2 - lon1) * pi / 180
  a <- sin(delta_phi / 2)^2 +
    cos(phi1) * cos(phi2) * sin(delta_lambda / 2)^2
  2 * earth_radius * asin(pmin(1, sqrt(a)))
}
sample_ids <- rownames(metadata)
sample_index <- seq_along(sample_ids)
geographic_matrix <- outer(
  sample_index, sample_index,
  Vectorize(function(i, j) {
    haversine_km(
      environment$Latitude[[i]], environment$Longitude[[i]],
      environment$Latitude[[j]], environment$Longitude[[j]]
    )
  })
)
dimnames(geographic_matrix) <- list(sample_ids, sample_ids)
bray_matrix <- as.matrix(vegan::vegdist(t(relative), method = "bray"))

pair_index <- which(upper.tri(bray_matrix), arr.ind = TRUE)
distance_pairs <- data.frame(
  Sample1 = sample_ids[pair_index[, "row"]],
  Sample2 = sample_ids[pair_index[, "col"]],
  Group1 = as.character(metadata$Group[pair_index[, "row"]]),
  Group2 = as.character(metadata$Group[pair_index[, "col"]]),
  DistanceKm = geographic_matrix[pair_index],
  BrayCurtis = bray_matrix[pair_index],
  stringsAsFactors = FALSE
)
distance_pairs$Similarity <- pmax(1 - distance_pairs$BrayCurtis, 1e-6)
distance_pairs$Scope <- ifelse(
  distance_pairs$Group1 == distance_pairs$Group2,
  distance_pairs$Group1, "Between regions"
)
distance_pairs$LogDistance <- log10(distance_pairs$DistanceKm + 1)
distance_pairs$LogSimilarity <- log10(distance_pairs$Similarity)
within_pairs <- distance_pairs[
  distance_pairs$Scope %in% levels(metadata$Group), , drop = FALSE
]

fit_distance_slope <- function(data, region) {
  fit <- stats::lm(LogSimilarity ~ LogDistance, data = data)
  interval <- stats::confint(fit, "LogDistance", level = 0.95)
  data.frame(
    Region = region,
    Pairs = nrow(data),
    MinDistanceKm = min(data$DistanceKm),
    MaxDistanceKm = max(data$DistanceKm),
    Slope = unname(stats::coef(fit)[["LogDistance"]]),
    Lower = unname(interval[[1L]]),
    Upper = unname(interval[[2L]]),
    RSquared = summary(fit)$r.squared,
    stringsAsFactors = FALSE
  )
}
distance_slope_ledger <- rbind(
  fit_distance_slope(distance_pairs, "Pooled"),
  do.call(rbind, lapply(levels(metadata$Group), function(group_name) {
    fit_distance_slope(
      within_pairs[within_pairs$Scope == group_name, , drop = FALSE],
      group_name
    )
  }))
)

set.seed(20260740)
blocked_permutations <- permute::shuffleSet(
  nrow(metadata), nset = 999L,
  control = permute::how(blocks = metadata$Group)
)
pooled_mantel <- vegan::mantel(
  stats::as.dist(bray_matrix), stats::as.dist(geographic_matrix),
  method = "spearman", permutations = blocked_permutations
)
mantel_ledger <- data.frame(
  Region = "Pooled; permutations within region",
  Samples = nrow(metadata),
  MantelRho = unname(pooled_mantel$statistic),
  PermutationP = pooled_mantel$signif,
  stringsAsFactors = FALSE
)
for (group_index in seq_along(levels(metadata$Group))) {
  group_name <- levels(metadata$Group)[[group_index]]
  rows <- metadata$Group == group_name
  set.seed(20260740L + group_index)
  fit <- vegan::mantel(
    stats::as.dist(bray_matrix[rows, rows]),
    stats::as.dist(geographic_matrix[rows, rows]),
    method = "spearman", permutations = 999L
  )
  mantel_ledger <- rbind(
    mantel_ledger,
    data.frame(
      Region = group_name,
      Samples = sum(rows),
      MantelRho = unname(fit$statistic),
      PermutationP = fit$signif,
      stringsAsFactors = FALSE
    )
  )
}

stopifnot(
  nrow(distance_pairs) == choose(90L, 2L),
  nrow(within_pairs) == 3L * choose(30L, 2L),
  max(distance_pairs$DistanceKm) > 4000,
  all(distance_slope_ledger$Slope < 0),
  all(mantel_ledger$PermutationP >= 1 / 1000 &
        mantel_ledger$PermutationP <= 1)
)
distance_slope_ledger
mantel_ledger

niche_distance_audit <- data.frame(
  Dimension = c(
    "Breadth estimand", "Breadth filter", "RAD replicate", "RAD choice",
    "Geographic distance", "Pair dependence", "Spatial scale", "Process claim"
  ),
  FragilePractice = c(
    "Call sample evenness a physiological niche",
    "Let singletons define specialists",
    "Pool samples into one curve",
    "Select a model without diagnostics",
    "Use Euclidean degrees as kilometres",
    "Run ordinary p-values on all pairs",
    "Mix within- and between-region pairs",
    "Equate decay with dispersal limitation"
  ),
  UpgradedContract = c(
    "Use sample-use breadth and state its limit",
    "Prespecify prevalence and abundance gates",
    "Fit each sample and retain winner variation",
    "Report AIC ledger and convergence warnings",
    "Use Haversine great-circle distance",
    "Use slopes descriptively and permute labels",
    "Lead with stratified curves and pooled audit",
    "Require environment and process-aware evidence"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(
  niche_distance_audit,
  caption = "Audit contract for macroecological pattern claims"
)

breadth_colours <- c(
  "Narrow decile" = pal_pub[["blue"]],
  "Intermediate" = pal_pub[["grey"]],
  "Broad decile" = pal_pub[["vermillion"]]
)
narrow_label_pool <- head(
  breadth_ledger[order(breadth_ledger$StandardizedBreadth), ], 40L
)
narrow_label_pool <- narrow_label_pool[
  !duplicated(narrow_label_pool$DisplayTaxon),
]
narrow_labels <- head(
  narrow_label_pool[order(nchar(narrow_label_pool$DisplayTaxon)), ], 3L
)
broad_labels <- head(
  breadth_ledger[order(-breadth_ledger$StandardizedBreadth), ], 6L
)
label_breadth <- rbind(narrow_labels, broad_labels)
breadth_plot <- ggplot(
  breadth_ledger,
  aes(x = Prevalence, y = StandardizedBreadth, colour = BreadthClass)
) +
  geom_point(size = 1.1, alpha = 0.48) +
  ggrepel::geom_text_repel(
    data = label_breadth,
    aes(label = DisplayTaxon),
    family = font_pub, size = 2.35, max.overlaps = Inf,
    min.segment.length = 0, show.legend = FALSE
  ) +
  scale_colour_manual(values = breadth_colours) +
  scale_x_continuous(
    labels = scales::label_percent(), limits = c(0.04, 1.03),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, max(breadth_ledger$StandardizedBreadth) * 1.08),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Sample-use breadth is constrained by occupancy",
    subtitle = sprintf(
      "%d OTUs; Spearman rho = %.3f",
      nrow(breadth_ledger), breadth_occupancy_rho
    ),
    x = "Sample occupancy", y = "Standardized Levins breadth",
    colour = "Breadth class",
    caption = "Deciles are descriptive within this filtered feature universe; breadth is not a physiological tolerance measurement."
  ) +
  theme_pub(base_size = 8.4) +
  theme(legend.position = "bottom")
save_pub(
  breadth_plot, "figures/40-niche-breadth", width = 180, height = 140
)
breadth_plot

rad_colours <- c(
  "Null" = pal_pub[["grey"]], "Preemption" = pal_pub[["orange"]],
  "Lognormal" = pal_pub[["green"]], "Zipf" = pal_pub[["blue"]],
  "Mandelbrot" = pal_pub[["vermillion"]]
)
rad_curve_panel <- ggplot() +
  geom_point(
    data = representative_observed,
    aes(x = Rank, y = Abundance), size = 0.75, alpha = 0.45
  ) +
  geom_line(
    data = representative_rad,
    aes(x = Rank, y = Abundance, colour = Model), linewidth = 0.65
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = rad_colours) +
  labs(
    title = "A. Representative rank-abundance fit",
    subtitle = sprintf(
      "%s (%s); richness = %d",
      representative_sample, metadata[representative_sample, "Group"],
      sample_richness[[representative_sample]]
    ),
    x = "Abundance rank", y = "Read count", colour = "RAD model"
  ) +
  theme_pub(base_size = 7.8) +
  theme(legend.position = "bottom")
rad_winner_panel <- ggplot(
  rad_model_summary,
  aes(x = reorder(Model, WinningSamples), y = WinningSamples, fill = Model)
) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.3) +
  geom_text(
    aes(label = WinningSamples), hjust = -0.15,
    family = font_pub, size = 2.7
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = rad_colours) +
  scale_y_continuous(
    limits = c(0, max(rad_model_summary$WinningSamples) * 1.14),
    expand = c(0, 0)
  ) +
  labs(
    title = "B. Best AIC model across samples",
    subtitle = sprintf("%d captured fitting warnings", nrow(rad_warning_ledger)),
    x = NULL, y = "Samples"
  ) +
  theme_pub(base_size = 7.8) +
  theme(legend.position = "none")
rad_plot <- patchwork::wrap_plots(
  rad_curve_panel, rad_winner_panel, nrow = 1L, widths = c(1.25, 0.75)
) +
  patchwork::plot_annotation(
    title = "RAD models summarize dominance, not unique mechanisms",
    caption = "Model winners are relative to the five-model candidate set; each biological sample was fitted independently."
  )
save_pub(
  rad_plot, "figures/40-rank-abundance-models",
  width = 180, height = 132
)
rad_plot

region_colours <- c(
  "CW" = pal_pub[["blue"]],
  "IW" = pal_pub[["orange"]],
  "TW" = pal_pub[["green"]]
)
within_pairs$Scope <- factor(
  within_pairs$Scope, levels = c("CW", "IW", "TW")
)
distance_decay_plot <- ggplot(
  within_pairs,
  aes(x = LogDistance, y = LogSimilarity, colour = Scope)
) +
  geom_point(size = 0.9, alpha = 0.32) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  facet_wrap(~Scope, scales = "free_x") +
  scale_colour_manual(values = region_colours) +
  labs(
    title = "Distance-decay differs across wetland regions",
    subtitle = "Bray-Curtis similarity regressed on Haversine distance within each region",
    x = "log10(Geographic distance [km] + 1)",
    y = "log10(Bray-Curtis similarity)", colour = "Region",
    caption = "Pairwise points are not independent; slopes are effect-size descriptions and permutation evidence is reported separately."
  ) +
  theme_pub(base_size = 8.2) +
  theme(legend.position = "none")
save_pub(
  distance_decay_plot, "figures/40-distance-decay",
  width = 180, height = 137
)
distance_decay_plot

plot_slopes <- distance_slope_ledger[
  distance_slope_ledger$Region != "Pooled",
]
slope_panel <- ggplot(
  plot_slopes,
  aes(x = Slope, y = Region, colour = Region)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#555555") +
  geom_errorbarh(
    aes(xmin = Lower, xmax = Upper), height = 0.15, linewidth = 0.75
  ) +
  geom_point(size = 3) +
  scale_colour_manual(values = region_colours) +
  labs(
    title = "A. Within-region decay slope",
    subtitle = "Estimate and ordinary-regression 95% CI",
    x = "Log-log slope", y = NULL
  ) +
  theme_pub(base_size = 8.0) +
  theme(legend.position = "none")
plot_mantel <- mantel_ledger[mantel_ledger$Region %in% c("CW", "IW", "TW"), ]
mantel_panel <- ggplot(
  plot_mantel,
  aes(x = MantelRho, y = Region, colour = Region)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#555555") +
  geom_point(size = 3) +
  geom_text(
    aes(label = sprintf("p = %.3f", PermutationP)),
    nudge_x = 0.055, hjust = 0, family = font_pub, size = 2.5
  ) +
  scale_colour_manual(values = region_colours) +
  scale_x_continuous(
    limits = c(min(0, min(plot_mantel$MantelRho) - 0.05),
               max(plot_mantel$MantelRho) + 0.18)
  ) +
  labs(
    title = "B. Mantel rank association",
    subtitle = "999 sample-label permutations per region",
    x = "Spearman Mantel rho", y = NULL
  ) +
  theme_pub(base_size = 8.0) +
  theme(legend.position = "none")
distance_audit_plot <- patchwork::wrap_plots(
  slope_panel, mantel_panel, nrow = 1L
) +
  patchwork::plot_annotation(
    title = "Effect size and permutation evidence answer different questions",
    caption = "Slope CIs ignore pair dependence; Mantel does not separate environment from dispersal.",
    theme = theme(
      plot.caption = element_text(
        family = font_pub, colour = "#666666", hjust = 0,
        margin = margin(t = 5, r = 4, b = 2, l = 8)
      )
    )
  )
save_pub(
  distance_audit_plot, "figures/40-distance-decay-audit",
  width = 180, height = 125
)
distance_audit_plot

result_dir <- "results/40-niche-distance-decay"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  breadth_ledger,
  file.path(result_dir, "niche-breadth-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  rad_aic_ledger,
  file.path(result_dir, "rad-model-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  rad_warning_ledger,
  file.path(result_dir, "rad-warning-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  distance_slope_ledger,
  file.path(result_dir, "distance-decay-slopes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  mantel_ledger,
  file.path(result_dir, "distance-mantel-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

expected_bases <- c(
  "figures/40-niche-breadth",
  "figures/40-rank-abundance-models",
  "figures/40-distance-decay",
  "figures/40-distance-decay-audit"
)
expected_files <- as.vector(outer(
  expected_bases, c(".pdf", ".svg", ".png", ".tiff"), paste0
))
stopifnot(
  all(file.exists(expected_files)),
  nrow(breadth_ledger) == 3208L,
  nrow(rad_aic_ledger) == 450L,
  sum(rad_model_summary$WinningSamples) == 90L,
  nrow(distance_pairs) == 4005L,
  nrow(within_pairs) == 1305L,
  nrow(mantel_ledger) == 4L,
  all(plot_slopes$Slope < 0)
)
