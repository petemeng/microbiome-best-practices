# 差异方法横评（ANCOM-BC2/ALDEx2/MaAsLin2/LEfSe/corncob）
# Run sequentially in a new working directory.
# Required packages: ALDEx2, ANCOMBC, Maaslin2, S4Vectors, SummarizedExperiment, corncob, ggplot2, knitr, lefser, phyloseq, ragg, svglite.

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

set.seed(20260731)
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
  lapply(taxonomy[rownames(otutab), rank_names], clean_taxon),
  check.names = FALSE
)
known_genus <- !is.na(taxonomy_clean$Genus)
lineage_key <- apply(
  taxonomy_clean[known_genus, rank_names[1:6], drop = FALSE], 1L,
  function(x) paste(ifelse(is.na(x), "?", x), collapse = ";")
)
genus_counts <- rowsum(
  otutab[known_genus, , drop = FALSE], lineage_key, reorder = FALSE
)
lineages <- rownames(genus_counts)
labels <- vapply(strsplit(lineages, ";", fixed = TRUE), tail, character(1), 1L)
feature_ids <- make.unique(paste0("Genus:", labels))
rownames(genus_counts) <- feature_ids
feature_info <- data.frame(
  FeatureID = feature_ids, DisplayTaxon = labels, Lineage = lineages,
  Reportable = TRUE, stringsAsFactors = FALSE, row.names = feature_ids
)
genus_counts <- rbind(
  genus_counts,
  Unclassified_residual = colSums(otutab[!known_genus, , drop = FALSE])
)
feature_info <- rbind(
  feature_info,
  data.frame(
    FeatureID = "Unclassified_residual",
    DisplayTaxon = "Unclassified residual",
    Lineage = "Unclassified at genus rank", Reportable = FALSE,
    row.names = "Unclassified_residual"
  )
)

prevalence <- rowMeans(genus_counts > 0)
total_reads <- rowSums(genus_counts)
eligible <- prevalence >= 0.10 & total_reads >= 20
selection_score <- prevalence * log1p(total_reads)
selected_ids <- names(sort(selection_score[eligible], decreasing = TRUE))[1:80]
counts <- genus_counts[selected_ids, , drop = FALSE]
feature_info <- feature_info[selected_ids, , drop = FALSE]

stopifnot(
  ncol(counts) == 60L,
  nrow(counts) == 80L,
  sum(metadata$Group == "IW") == 30L,
  sum(metadata$Group == "CW") == 30L,
  "Unclassified_residual" %in% rownames(counts),
  sum(feature_info$Reportable) == 79L
)

data.frame(
  Samples = ncol(counts), TestedRows = nrow(counts),
  ReportableGenera = sum(feature_info$Reportable),
  MinPrevalence = min(prevalence[selected_ids]),
  ReadFractionRetained = sum(counts) / sum(otutab)
)

taxonomy_model <- matrix(
  feature_info$DisplayTaxon,
  ncol = 1L,
  dimnames = list(rownames(feature_info), "Genus")
)
ps <- phyloseq::phyloseq(
  phyloseq::otu_table(counts, taxa_are_rows = TRUE),
  phyloseq::tax_table(taxonomy_model),
  phyloseq::sample_data(metadata)
)

descriptive_direction <- rowMeans(
  log2(counts[, metadata$Group == "CW", drop = FALSE] + 0.5)
) - rowMeans(
  log2(counts[, metadata$Group == "IW", drop = FALSE] + 0.5)
)

set.seed(20260731)
ancom_warnings <- character()
ancom_fit <- withCallingHandlers(
  ANCOMBC::ancombc2(
    data = ps,
    fix_formula = "Group",
    p_adj_method = "BH",
    pseudo_sens = FALSE,
    prv_cut = 0,
    lib_cut = 0,
    group = "Group",
    struc_zero = TRUE,
    neg_lb = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    n_cl = 1,
    verbose = FALSE
  ),
  warning = function(w) {
    ancom_warnings <<- c(ancom_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
ancom_res <- ancom_fit$res
ancom_results <- data.frame(
  FeatureID = ancom_res$taxon,
  Method = "ANCOM-BC2",
  EffectCW = ancom_res$lfc_GroupCW,
  StandardError = ancom_res$se_GroupCW,
  PValue = ancom_res$p_GroupCW,
  QValue = ancom_res$q_GroupCW,
  stringsAsFactors = FALSE
)
ancom_results$Significant <- with(
  ancom_results,
  is.finite(EffectCW) & is.finite(QValue) & QValue < 0.05
)
ancom_diagnostics <- data.frame(
  Rows = nrow(ancom_results),
  FiniteGroupEffects = sum(
    is.finite(ancom_results$EffectCW) & is.finite(ancom_results$StandardError)
  ),
  CapturedWarnings = length(ancom_warnings),
  stringsAsFactors = FALSE
)
ancom_diagnostics

set.seed(20260731)
aldex_clr <- ALDEx2::aldex.clr(
  counts,
  as.character(metadata$Group),
  mc.samples = 128,
  denom = "iqlr",
  verbose = FALSE
)
aldex_test <- ALDEx2::aldex.ttest(
  aldex_clr, paired.test = FALSE, verbose = FALSE
)
aldex_effect <- ALDEx2::aldex.effect(
  aldex_clr, paired.test = FALSE, verbose = FALSE
)
aldex_results <- data.frame(
  FeatureID = rownames(aldex_test),
  Method = "ALDEx2",
  EffectCW = -aldex_effect$effect,
  StandardError = NA_real_,
  PValue = aldex_test$we.ep,
  QValue = aldex_test$we.eBH,
  stringsAsFactors = FALSE
)
aldex_results$Significant <- with(
  aldex_results,
  is.finite(EffectCW) & is.finite(QValue) & QValue < 0.05
)

maaslin_output <- file.path(
  tempdir(), paste0("article31-maaslin2-", Sys.getpid())
)
dir.create(maaslin_output, recursive = TRUE, showWarnings = FALSE)
maaslin_console <- utils::capture.output({
  maaslin_fit <- Maaslin2::Maaslin2(
    input_data = as.data.frame(t(counts), check.names = FALSE),
    input_metadata = data.frame(
      Group = metadata$Group, row.names = rownames(metadata)
    ),
    output = maaslin_output,
    min_abundance = 0,
    min_prevalence = 0,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    fixed_effects = "Group",
    correction = "BH",
    standardize = FALSE,
    cores = 1,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    save_models = FALSE,
    max_significance = 1,
    reference = "Group,IW"
  )
}, type = "output")
maaslin_table <- maaslin_fit$results
maaslin_table <- maaslin_table[
  maaslin_table$metadata == "Group" & maaslin_table$value == "CW",
]
maaslin_feature_map <- setNames(
  rownames(counts),
  make.names(rownames(counts), unique = TRUE)
)
maaslin_results <- data.frame(
  FeatureID = unname(maaslin_feature_map[maaslin_table$feature]),
  Method = "MaAsLin2",
  EffectCW = maaslin_table$coef,
  StandardError = maaslin_table$stderr,
  PValue = maaslin_table$pval,
  QValue = maaslin_table$qval,
  stringsAsFactors = FALSE
)
maaslin_results$Significant <- with(
  maaslin_results,
  is.finite(EffectCW) & is.finite(QValue) & QValue < 0.05
)

relative_per_million <- sweep(counts, 2L, colSums(counts), "/") * 1e6
lefse_input <- SummarizedExperiment::SummarizedExperiment(
  assays = list(relative_abundance = relative_per_million),
  colData = S4Vectors::DataFrame(
    GROUP = metadata$Group, row.names = rownames(metadata)
  )
)
set.seed(20260731)
lefse_fit <- lefser::lefser(
  lefse_input,
  kruskal.threshold = 1,
  wilcox.threshold = 1,
  lda.threshold = 0,
  groupCol = "GROUP",
  checkAbundances = TRUE
)

lefse_p <- vapply(rownames(counts), function(feature_id) {
  stats::wilcox.test(
    relative_per_million[feature_id, metadata$Group == "CW"],
    relative_per_million[feature_id, metadata$Group == "IW"],
    exact = FALSE
  )$p.value
}, numeric(1))
lefse_q <- stats::p.adjust(lefse_p, method = "BH")
lefse_scores <- setNames(lefse_fit$scores, lefse_fit$Names)
lefse_results <- data.frame(
  FeatureID = rownames(counts),
  Method = "LEfSe + BH",
  EffectCW = sign(descriptive_direction) * abs(lefse_scores[rownames(counts)]),
  StandardError = NA_real_,
  PValue = unname(lefse_p[rownames(counts)]),
  QValue = unname(lefse_q[rownames(counts)]),
  stringsAsFactors = FALSE
)
lefse_results$Significant <- with(
  lefse_results,
  is.finite(EffectCW) & is.finite(QValue) &
    QValue < 0.05 & abs(EffectCW) >= 2
)

set.seed(20260731)
corncob_fit <- corncob::differentialTest(
  formula = ~ Group,
  phi.formula = ~ Group,
  formula_null = ~ 1,
  phi.formula_null = ~ Group,
  data = ps,
  test = "Wald",
  boot = FALSE,
  fdr_cutoff = 1,
  full_output = TRUE,
  verbose = FALSE
)
corncob_rows <- lapply(seq_along(corncob_fit$p), function(i) {
  feature_id <- names(corncob_fit$p)[[i]]
  model_summary <- corncob_fit$all_models[[i]]
  coefficient_table <- if (is.null(model_summary)) NULL else coef(model_summary)
  coefficient_name <- "mu.GroupCW"
  has_coefficient <- !is.null(coefficient_table) &&
    coefficient_name %in% rownames(coefficient_table)
  data.frame(
    FeatureID = feature_id,
    Method = "corncob",
    EffectCW = if (has_coefficient) coefficient_table[coefficient_name, "Estimate"] else NA_real_,
    StandardError = if (has_coefficient) coefficient_table[coefficient_name, "Std. Error"] else NA_real_,
    PValue = unname(corncob_fit$p[[i]]),
    QValue = unname(corncob_fit$p_fdr[[i]]),
    stringsAsFactors = FALSE
  )
})
corncob_results <- do.call(rbind, corncob_rows)
corncob_results$Significant <- with(
  corncob_results,
  is.finite(EffectCW) & is.finite(QValue) & QValue < 0.05
)

method_order <- c("ANCOM-BC2", "ALDEx2", "MaAsLin2", "LEfSe + BH", "corncob")
method_results_all <- rbind(
  ancom_results, aldex_results, maaslin_results, lefse_results, corncob_results
)
method_results_all$DisplayTaxon <- feature_info[
  method_results_all$FeatureID, "DisplayTaxon"
]
method_results_all$Reportable <- feature_info[
  method_results_all$FeatureID, "Reportable"
]
method_results_all$Method <- factor(method_results_all$Method, levels = method_order)
method_results <- method_results_all[
  !is.na(method_results_all$Reportable) & method_results_all$Reportable,
]

direction_check <- do.call(rbind, lapply(method_order, function(method_name) {
  z <- method_results[method_results$Method == method_name, ]
  common <- intersect(z$FeatureID, names(descriptive_direction))
  data.frame(
    Method = method_name,
    FiniteRows = sum(is.finite(z$EffectCW)),
    SpearmanWithDescriptiveDirection = suppressWarnings(stats::cor(
      z$EffectCW[match(common, z$FeatureID)],
      descriptive_direction[common],
      method = "spearman", use = "complete.obs"
    )),
    stringsAsFactors = FALSE
  )
}))

hit_sets <- setNames(lapply(method_order, function(method_name) {
  z <- method_results[
    method_results$Method == method_name & method_results$Significant,
  ]
  unique(z$FeatureID)
}), method_order)

hit_summary <- data.frame(
  Method = factor(method_order, levels = method_order),
  Tested = vapply(method_order, function(x) {
    sum(method_results$Method == x & is.finite(method_results$QValue))
  }, integer(1)),
  Significant = lengths(hit_sets),
  stringsAsFactors = FALSE
)

jaccard_matrix <- outer(method_order, method_order, Vectorize(function(a, b) {
  union_set <- union(hit_sets[[a]], hit_sets[[b]])
  if (length(union_set) == 0L) 1 else
    length(intersect(hit_sets[[a]], hit_sets[[b]])) / length(union_set)
}))
dimnames(jaccard_matrix) <- list(method_order, method_order)

feature_support <- aggregate(
  Significant ~ FeatureID + DisplayTaxon,
  data = method_results,
  FUN = sum
)
names(feature_support)[names(feature_support) == "Significant"] <- "MethodSupport"
feature_support <- feature_support[order(
  -feature_support$MethodSupport, feature_support$DisplayTaxon
), ]

stopifnot(
  nrow(method_results_all) == 400L,
  all(direction_check$FiniteRows >= 75L),
  all(direction_check$SpearmanWithDescriptiveDirection > 0.5),
  all(diag(jaccard_matrix) == 1),
  all(jaccard_matrix >= 0 & jaccard_matrix <= 1),
  ancom_diagnostics$FiniteGroupEffects >= 76L
)

hit_summary
direction_check
head(feature_support, 15)

method_audit <- data.frame(
  Method = method_order,
  PrimaryInput = c(
    "Counts", "Counts", "Counts; internal TSS + LOG",
    "Relative abundance", "Taxon count / library total"
  ),
  ReportedEffect = c(
    "Bias-corrected log abundance coefficient",
    "Standardized CLR effect",
    "Linear-model coefficient after transform",
    "Signed LDA score",
    "Beta-binomial mean logit coefficient"
  ),
  UpgradeApplied = c(
    "Finite coefficient and warning audit",
    "128 Monte Carlo samples; IQLR reference",
    "Reference level, normalization and transform fixed",
    "All-feature BH plus |LDA| >= 2",
    "Mean and dispersion formulas reported separately"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(method_audit, caption = "Method-specific estimands and audit upgrades")

hit_count_plot <- ggplot(hit_summary, aes(x = Method, y = Significant, fill = Method)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(
    aes(label = paste0(Significant, " / ", Tested)),
    vjust = -0.45, family = font_pub, size = 3
  ) +
  scale_fill_manual(values = c(
    "ANCOM-BC2" = pal_pub[["blue"]],
    "ALDEx2" = pal_pub[["orange"]],
    "MaAsLin2" = pal_pub[["green"]],
    "LEfSe + BH" = pal_pub[["purple"]],
    "corncob" = pal_pub[["vermillion"]]
  )) +
  scale_y_continuous(
    limits = c(0, max(hit_summary$Significant) * 1.18 + 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Differential-abundance methods return different hit counts",
    subtitle = "Same 60 samples, 79 reportable genera, CW versus IW, BH q < 0.05",
    x = NULL, y = "Significant reportable genera",
    caption = "Labels show significant / finite-tested genera; LEfSe additionally requires |LDA| >= 2."
  ) +
  theme_pub(base_size = 8.8) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

save_pub(hit_count_plot, "figures/31-da-hit-counts", width = 165, height = 107)
hit_count_plot

jaccard_plot_data <- as.data.frame(as.table(jaccard_matrix))
names(jaccard_plot_data) <- c("Method1", "Method2", "Jaccard")
jaccard_plot <- ggplot(
  jaccard_plot_data,
  aes(x = Method1, y = Method2, fill = Jaccard)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Jaccard)), family = font_pub, size = 3) +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#9ECAE1", "#3182BD", "#08306B"),
    limits = c(0, 1)
  ) +
  coord_equal() +
  labs(
    title = "Hit-set agreement is method dependent",
    subtitle = "Jaccard = intersection / union",
    x = NULL, y = NULL, fill = "Jaccard",
    caption = "A low value is a disagreement signal to audit, not permission to take the union."
  ) +
  theme_pub(base_size = 8.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )

save_pub(jaccard_plot, "figures/31-da-jaccard", width = 152, height = 132)
jaccard_plot

top_support <- head(feature_support$FeatureID, 20L)
evidence_plot_data <- method_results[
  method_results$FeatureID %in% top_support,
]
evidence_plot_data$DisplayTaxon <- factor(
  evidence_plot_data$DisplayTaxon,
  levels = rev(feature_support$DisplayTaxon[
    match(top_support, feature_support$FeatureID)
  ])
)
evidence_plot_data$Direction <- ifelse(
  evidence_plot_data$EffectCW >= 0, "Higher in CW", "Higher in IW"
)
evidence_plot_data$Evidence <- ifelse(
  evidence_plot_data$Significant, "q < 0.05", "Not significant"
)

evidence_map_plot <- ggplot(
  evidence_plot_data,
  aes(x = Method, y = DisplayTaxon)
) +
  geom_point(
    aes(colour = Direction, shape = Evidence, size = -log10(pmax(QValue, 1e-12))),
    alpha = 0.86, stroke = 0.7
  ) +
  scale_colour_manual(values = c(
    "Higher in CW" = pal_pub[["vermillion"]],
    "Higher in IW" = pal_pub[["blue"]]
  )) +
  scale_shape_manual(values = c("q < 0.05" = 16, "Not significant" = 1)) +
  scale_size_continuous(range = c(1.3, 5.0), name = expression(-log[10](q))) +
  labs(
    title = "Consensus requires taxon-level direction and evidence",
    subtitle = "Twenty genera with the largest number of supporting methods",
    x = NULL, y = "Genus", colour = "Direction", shape = "Evidence",
    caption = "Effect magnitudes are not compared because the five methods use different scales."
  ) +
  theme_pub(base_size = 8.0) +
  theme(
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.25),
    axis.text.x = element_text(angle = 28, hjust = 1),
    legend.position = "bottom", legend.box = "vertical"
  )

save_pub(evidence_map_plot, "figures/31-da-evidence-map", width = 180, height = 158)
evidence_map_plot

expected_bases <- c(
  "figures/31-da-hit-counts", "figures/31-da-jaccard",
  "figures/31-da-evidence-map"
)
expected_files <- as.vector(outer(
  expected_bases, c(".pdf", ".svg", ".png", ".tiff"), paste0
))
stopifnot(
  all(file.exists(expected_files)),
  nrow(hit_summary) == 5L,
  nrow(direction_check) == 5L,
  max(feature_support$MethodSupport) >= 1L
)
