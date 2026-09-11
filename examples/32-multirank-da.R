# 多分类级差异（差异科/属/种，collapse 逐级）
# Run sequentially in a new working directory.
# Required packages: ALDEx2, ggplot2, ggrepel, knitr, patchwork, ragg, scales, svglite.

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

set.seed(20260732)
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
      plot.title = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.15)
      ),
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
metadata_all <- read_keyed_tsv("data/small/metadata.tsv")
otutab_all <- as.matrix(data.frame(
  lapply(otutab_raw, as.integer), check.names = FALSE
))
rownames(otutab_all) <- rownames(otutab_raw)
sample_ids <- rownames(metadata_all)[metadata_all$Group %in% c("IW", "CW")]
metadata <- metadata_all[sample_ids, , drop = FALSE]
metadata$Group <- relevel(factor(metadata$Group), ref = "IW")
otutab <- otutab_all[, sample_ids, drop = FALSE]

stopifnot(
  nrow(otutab_all) == 13628L,
  ncol(otutab_all) == 90L,
  sum(otutab_all) == 1619670L,
  ncol(otutab) == 60L,
  all(otutab >= 0L),
  all(otutab == round(otutab))
)

head(otutab[, 1:5])
head(taxonomy)
head(metadata)

rank_names <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
)
clean_taxon <- function(x) {
  x <- sub("^[A-Za-z]__", "", trimws(as.character(x)))
  x[
    is.na(x) | x == "" |
      grepl(
        "unclassified|unknown|unassigned|uncultured|metagenome",
        x, ignore.case = TRUE
      )
  ] <- NA_character_
  x
}
taxonomy_clean <- as.data.frame(
  lapply(taxonomy[rownames(otutab_all), rank_names], clean_taxon),
  check.names = FALSE
)

aggregate_rank <- function(rank, sample_ids = colnames(otutab)) {
  rank_index <- match(rank, rank_names)
  stopifnot(!is.na(rank_index))
  use_ranks <- rank_names[seq_len(rank_index)]
  known <- !is.na(taxonomy_clean[[rank]])
  lineage <- apply(
    taxonomy_clean[known, use_ranks, drop = FALSE], 1L,
    function(x) paste(ifelse(is.na(x), "?", x), collapse = ";")
  )
  counts <- rowsum(
    otutab_all[known, sample_ids, drop = FALSE], lineage, reorder = FALSE
  )
  lineages <- rownames(counts)
  labels <- vapply(
    strsplit(lineages, ";", fixed = TRUE), tail, character(1), 1L
  )
  feature_ids <- make.unique(paste0(rank, ":", labels))
  rownames(counts) <- feature_ids

  lineage_parts <- strsplit(lineages, ";", fixed = TRUE)
  info <- data.frame(
    FeatureID = feature_ids,
    Rank = rank,
    DisplayTaxon = labels,
    Lineage = lineages,
    ParentTaxon = if (rank_index > 1L) vapply(
      lineage_parts, function(x) x[[rank_index - 1L]], character(1)
    ) else NA_character_,
    Reportable = TRUE,
    stringsAsFactors = FALSE,
    row.names = feature_ids
  )
  residual_id <- paste0("Unclassified_", tolower(rank), "_residual")
  residual_counts <- matrix(
    colSums(otutab_all[!known, sample_ids, drop = FALSE]),
    nrow = 1L,
    dimnames = list(residual_id, sample_ids)
  )
  counts <- rbind(counts, residual_counts)
  residual_id <- tail(rownames(counts), 1L)
  info <- rbind(
    info,
    data.frame(
      FeatureID = residual_id, Rank = rank,
      DisplayTaxon = paste("Unclassified", tolower(rank), "residual"),
      Lineage = paste("Unclassified at", rank, "rank"),
      ParentTaxon = NA_character_, Reportable = FALSE,
      row.names = residual_id
    )
  )
  list(
    counts = counts,
    info = info,
    known_asvs = sum(known),
    known_read_fraction = sum(otutab_all[known, , drop = FALSE]) /
      sum(otutab_all)
  )
}

analyze_rank <- function(rank, seed) {
  aggregated <- aggregate_rank(rank)
  counts_all <- aggregated$counts
  prevalence <- rowMeans(counts_all > 0)
  total_reads <- rowSums(counts_all)
  eligible <- prevalence >= 0.10 & total_reads >= 20
  score <- prevalence * log1p(total_reads)
  selected <- names(sort(score[eligible], decreasing = TRUE))
  selected <- head(selected, 80L)
  counts <- counts_all[selected, , drop = FALSE]
  info <- aggregated$info[selected, , drop = FALSE]
  denominator <- if (rank == "Species") "all" else "iqlr"

  stopifnot(nrow(counts) >= 4L, ncol(counts) == 60L)
  set.seed(seed)
  clr <- ALDEx2::aldex.clr(
    counts,
    as.character(metadata$Group),
    mc.samples = 128,
    denom = denominator,
    verbose = FALSE
  )
  test <- ALDEx2::aldex.ttest(
    clr, paired.test = FALSE, verbose = FALSE
  )
  effect <- ALDEx2::aldex.effect(
    clr, paired.test = FALSE, verbose = FALSE
  )
  results <- data.frame(
    FeatureID = rownames(test),
    Rank = rank,
    DisplayTaxon = info[rownames(test), "DisplayTaxon"],
    ParentTaxon = info[rownames(test), "ParentTaxon"],
    Lineage = info[rownames(test), "Lineage"],
    Reportable = info[rownames(test), "Reportable"],
    Prevalence = prevalence[rownames(test)],
    TotalReads = total_reads[rownames(test)],
    EffectCW = -effect$effect,
    DifferenceCW = -effect$diff.btw,
    PValue = test$we.ep,
    QValue = test$we.eBH,
    stringsAsFactors = FALSE
  )
  results$StatisticallySignificant <- with(
    results,
    Reportable & is.finite(QValue) & QValue < 0.05
  )

  named_tested <- sum(results$Reportable)
  reporting_gate <- aggregated$known_read_fraction >= 0.05 &&
    named_tested >= 20L
  results$PrimaryReportingEligible <- reporting_gate
  results$ReportedSignificant <- results$StatisticallySignificant & reporting_gate
  list(
    results = results,
    summary = data.frame(
      Rank = rank,
      KnownASVs = aggregated$known_asvs,
      KnownReadFraction = aggregated$known_read_fraction,
      NamedTaxaTested = named_tested,
      StatisticalHits = sum(results$StatisticallySignificant),
      ReportedHits = sum(results$ReportedSignificant),
      ReportingGate = ifelse(reporting_gate, "PASS", "FAIL"),
      Denominator = denominator,
      stringsAsFactors = FALSE
    )
  )
}

rank_order <- c("Family", "Genus", "Species")
rank_runs <- Map(
  analyze_rank,
  rank = rank_order,
  seed = 20260732L + seq_along(rank_order) - 1L
)
names(rank_runs) <- rank_order
multirank_results <- do.call(rbind, lapply(rank_runs, `[[`, "results"))
rank_summary <- do.call(rbind, lapply(rank_runs, `[[`, "summary"))
rank_summary$Rank <- factor(rank_summary$Rank, levels = rank_order)

stopifnot(
  nrow(rank_summary) == 3L,
  rank_summary$ReportingGate[rank_summary$Rank == "Family"] == "PASS",
  rank_summary$ReportingGate[rank_summary$Rank == "Genus"] == "PASS",
  rank_summary$ReportingGate[rank_summary$Rank == "Species"] == "FAIL",
  rank_summary$KnownReadFraction[rank_summary$Rank == "Species"] < 0.01,
  all(is.finite(multirank_results$EffectCW)),
  all(multirank_results$QValue >= 0 & multirank_results$QValue <= 1)
)

rank_summary

family_results <- multirank_results[
  multirank_results$Rank == "Family" & multirank_results$Reportable,
]
genus_results <- multirank_results[
  multirank_results$Rank == "Genus" & multirank_results$Reportable,
]

family_effect_by_name <- tapply(
  family_results$EffectCW,
  family_results$DisplayTaxon,
  function(x) x[[which.max(abs(x))]]
)
lineage_audit <- genus_results[
  genus_results$ParentTaxon %in% names(family_effect_by_name),
  c(
    "FeatureID", "DisplayTaxon", "ParentTaxon", "EffectCW",
    "QValue", "ReportedSignificant"
  )
]
names(lineage_audit)[names(lineage_audit) == "EffectCW"] <- "GenusEffectCW"
names(lineage_audit)[names(lineage_audit) == "QValue"] <- "GenusQValue"
lineage_audit$FamilyEffectCW <- family_effect_by_name[lineage_audit$ParentTaxon]
lineage_audit$DirectionAgreement <- sign(lineage_audit$GenusEffectCW) ==
  sign(lineage_audit$FamilyEffectCW)

lineage_audit <- lineage_audit[order(
  lineage_audit$GenusQValue, -abs(lineage_audit$GenusEffectCW)
), ]
head(lineage_audit, 15)

rank_audit <- data.frame(
  Decision = c(
    "Aggregation key", "Unknown reads", "Filtering", "FDR family",
    "Rank eligibility", "Parent-child interpretation"
  ),
  FragilePractice = c(
    "Collapse by the final label only",
    "Drop or redistribute unknown reads",
    "Reuse the genus hit list at every rank",
    "Carry genus q-values to family/species",
    "Assume lower rank is always better",
    "Treat significance as inherited"
  ),
  UpgradedContract = c(
    "Use complete lineage through the target rank",
    "Keep a non-reportable denominator residual",
    "Recompute prevalence and total reads per rank",
    "Refit and adjust within each prespecified rank",
    "Gate on annotation coverage and tested taxa",
    "Audit direction and aggregation; test each node"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(rank_audit, caption = "Multi-rank differential-abundance audit")

rank_plot_data <- rank_summary
rank_plot_data$GateLabel <- paste0(
  rank_plot_data$ReportingGate, "\n",
  rank_plot_data$ReportedHits, " reported hits"
)

coverage_panel <- ggplot(
  rank_plot_data,
  aes(x = Rank, y = 100 * KnownReadFraction, fill = Rank)
) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_hline(yintercept = 5, linetype = "22", colour = "#555555") +
  geom_text(
    aes(label = scales::percent(KnownReadFraction, accuracy = 0.1)),
    vjust = -0.45, family = font_pub, size = 2.8
  ) +
  scale_fill_manual(values = c(
    Family = pal_pub[["blue"]], Genus = pal_pub[["orange"]],
    Species = pal_pub[["grey"]]
  )) +
  scale_y_continuous(
    limits = c(0, max(100 * rank_plot_data$KnownReadFraction) * 1.18),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "A. Annotation coverage",
    x = NULL, y = "Reads with a reportable label (%)",
    caption = "Dashed line: prespecified 5% read-coverage gate."
  ) +
  theme_pub(base_size = 8.2)

testing_panel <- ggplot(
  rank_plot_data,
  aes(x = Rank, y = NamedTaxaTested, fill = ReportingGate)
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = GateLabel), vjust = -0.35,
    family = font_pub, size = 2.7, lineheight = 0.9
  ) +
  scale_fill_manual(values = c(PASS = pal_pub[["green"]], FAIL = pal_pub[["grey"]])) +
  scale_y_continuous(
    limits = c(0, max(rank_plot_data$NamedTaxaTested) * 1.30),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "B. Rank-specific testing and reporting",
    x = NULL, y = "Reportable taxa tested", fill = "Reporting gate",
    caption = "Statistical hits at a failed rank remain exploratory."
  ) +
  theme_pub(base_size = 8.2) +
  theme(legend.position = "bottom")

evidence_cascade_plot <- patchwork::wrap_plots(
  coverage_panel, testing_panel, nrow = 1L
) + patchwork::plot_annotation(
    title = "Taxonomic resolution must pass an evidence gate",
    subtitle = "Same samples and contrast; each rank is independently collapsed, filtered and tested"
  )
save_pub(
  evidence_cascade_plot, "figures/32-rank-evidence-cascade",
  width = 180, height = 112
)
evidence_cascade_plot

effect_candidates <- do.call(rbind, lapply(rank_order, function(rank) {
  z <- multirank_results[
    multirank_results$Rank == rank & multirank_results$Reportable,
  ]
  z <- z[order(z$QValue, -abs(z$EffectCW)), ]
  head(z, 12L)
}))
effect_candidates$DisplayLabel <- paste0(
  effect_candidates$DisplayTaxon, "  [", effect_candidates$Rank, "]"
)
effect_candidates$DisplayLabel <- factor(
  effect_candidates$DisplayLabel,
  levels = rev(unique(effect_candidates$DisplayLabel))
)
effect_candidates$Evidence <- ifelse(
  effect_candidates$ReportedSignificant,
  "Reportable q < 0.05",
  ifelse(
    effect_candidates$StatisticallySignificant,
    "Exploratory q < 0.05",
    "Not significant"
  )
)

multirank_effect_plot <- ggplot(
  effect_candidates,
  aes(x = EffectCW, y = DisplayLabel, colour = Rank, shape = Evidence)
) +
  geom_vline(xintercept = 0, linetype = "22", colour = "#777777") +
  geom_segment(aes(x = 0, xend = EffectCW, yend = DisplayLabel), alpha = 0.35) +
  geom_point(size = 2.5, stroke = 0.75) +
  scale_colour_manual(values = c(
    Family = pal_pub[["blue"]], Genus = pal_pub[["orange"]],
    Species = pal_pub[["grey"]]
  )) +
  scale_shape_manual(values = c(
    "Reportable q < 0.05" = 16,
    "Exploratory q < 0.05" = 1,
    "Not significant" = 4
  )) +
  labs(
    title = "Rank-specific effects cannot outrun annotation evidence",
    subtitle = "ALDEx2 standardized effect; positive values indicate higher abundance in CW",
    x = "ALDEx2 effect (CW versus IW)", y = NULL,
    colour = "Taxonomic rank", shape = "Evidence status",
    caption = "Species points are exploratory because the rank-level annotation gate failed."
  ) +
  theme_pub(base_size = 7.7) +
  theme(legend.position = "bottom", legend.box = "vertical")

save_pub(
  multirank_effect_plot, "figures/32-multirank-effect-map",
  width = 180, height = 152
)
multirank_effect_plot

lineage_plot_data <- head(lineage_audit, 30L)
lineage_plot_data$Evidence <- ifelse(
  lineage_plot_data$ReportedSignificant, "Genus q < 0.05", "Genus not significant"
)

lineage_coherence_plot <- ggplot(
  lineage_plot_data,
  aes(x = FamilyEffectCW, y = GenusEffectCW)
) +
  geom_hline(yintercept = 0, colour = "#888888", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "#888888", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#777777") +
  geom_point(
    aes(colour = DirectionAgreement, shape = Evidence),
    size = 2.3, alpha = 0.85, stroke = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = head(lineage_plot_data, 10L),
    aes(label = DisplayTaxon),
    family = font_pub, size = 2.35, min.segment.length = 0,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    `TRUE` = pal_pub[["green"]], `FALSE` = pal_pub[["vermillion"]]
  ), labels = c(`TRUE` = "Same direction", `FALSE` = "Opposite direction")) +
  scale_shape_manual(values = c("Genus q < 0.05" = 16, "Genus not significant" = 1)) +
  labs(
    title = "Family and genus effects are related, not inherited",
    subtitle = "Thirty genera with the smallest q-values and a tested parent family",
    x = "Parent-family ALDEx2 effect", y = "Genus ALDEx2 effect",
    colour = "Direction audit", shape = "Genus evidence",
    caption = "A family aggregates every included child; one child's significance does not propagate."
  ) +
  theme_pub(base_size = 8.3) +
  theme(legend.position = "bottom", legend.box = "vertical")

save_pub(
  lineage_coherence_plot, "figures/32-family-genus-coherence",
  width = 170, height = 127
)
lineage_coherence_plot

expected_bases <- c(
  "figures/32-rank-evidence-cascade",
  "figures/32-multirank-effect-map",
  "figures/32-family-genus-coherence"
)
expected_files <- as.vector(outer(
  expected_bases, c(".pdf", ".svg", ".png", ".tiff"), paste0
))
stopifnot(
  all(file.exists(expected_files)),
  identical(as.character(rank_summary$Rank), rank_order),
  sum(rank_summary$ReportingGate == "PASS") == 2L,
  sum(multirank_results$ReportedSignificant[
    multirank_results$Rank == "Species"
  ]) == 0L,
  nrow(lineage_audit) >= 20L
)
