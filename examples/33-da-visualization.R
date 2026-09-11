# 差异可视化：火山/LEfSe cladogram/Manhattan/森林图
# Run sequentially in a new working directory.
# Required packages: ANCOMBC, S4Vectors, SummarizedExperiment, ggplot2, ggraph, ggrepel, igraph, knitr, lefser, phyloseq, ragg, readr, svglite.

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

set.seed(20260733)
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
lineage_matrix <- do.call(rbind, strsplit(lineages, ";", fixed = TRUE))
colnames(lineage_matrix) <- rank_names[1:6]
labels <- lineage_matrix[, "Genus"]
feature_ids <- make.unique(paste0("Genus:", labels))
rownames(genus_counts) <- feature_ids
feature_info <- data.frame(
  FeatureID = feature_ids,
  as.data.frame(lineage_matrix, stringsAsFactors = FALSE),
  DisplayTaxon = labels,
  Lineage = lineages,
  Reportable = TRUE,
  stringsAsFactors = FALSE,
  row.names = feature_ids,
  check.names = FALSE
)
genus_counts <- rbind(
  genus_counts,
  Unclassified_residual = colSums(otutab[!known_genus, , drop = FALSE])
)
feature_info <- rbind(
  feature_info,
  data.frame(
    FeatureID = "Unclassified_residual",
    Kingdom = NA, Phylum = NA, Class = NA, Order = NA, Family = NA, Genus = NA,
    DisplayTaxon = "Unclassified residual",
    Lineage = "Unclassified at genus rank", Reportable = FALSE,
    row.names = "Unclassified_residual", check.names = FALSE
  )
)

prevalence <- rowMeans(genus_counts > 0)
total_reads <- rowSums(genus_counts)
eligible <- prevalence >= 0.10 & total_reads >= 20
score <- prevalence * log1p(total_reads)
selected_ids <- names(sort(score[eligible], decreasing = TRUE))[1:80]
counts <- genus_counts[selected_ids, , drop = FALSE]
feature_info <- feature_info[selected_ids, , drop = FALSE]

stopifnot(
  nrow(counts) == 80L, ncol(counts) == 60L,
  sum(feature_info$Reportable) == 79L,
  "Unclassified_residual" %in% rownames(counts)
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

set.seed(20260733)
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
master_results <- data.frame(
  FeatureID = ancom_res$taxon,
  EffectCW = ancom_res$lfc_GroupCW,
  StandardError = ancom_res$se_GroupCW,
  PValue = ancom_res$p_GroupCW,
  QValue = ancom_res$q_GroupCW,
  stringsAsFactors = FALSE
)
master_results <- cbind(
  master_results,
  feature_info[master_results$FeatureID, c(
    "Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
    "DisplayTaxon", "Lineage", "Reportable"
  )]
)
master_results$CI95Lower <- master_results$EffectCW -
  1.96 * master_results$StandardError
master_results$CI95Upper <- master_results$EffectCW +
  1.96 * master_results$StandardError
master_results$StatisticallySignificant <- with(
  master_results,
  Reportable & is.finite(QValue) & QValue < 0.05
)
practical_threshold <- log(2)
master_results$Highlighted <- with(
  master_results,
  StatisticallySignificant & abs(EffectCW) >= practical_threshold
)
master_results$Direction <- ifelse(
  master_results$EffectCW >= 0, "Higher in CW", "Higher in IW"
)

stopifnot(
  nrow(master_results) == 80L,
  sum(master_results$Reportable) == 79L,
  sum(is.finite(master_results$EffectCW) & is.finite(master_results$StandardError)) >= 76L,
  all(master_results$CI95Lower <= master_results$CI95Upper, na.rm = TRUE)
)

head(master_results[order(master_results$QValue), c(
  "FeatureID", "DisplayTaxon", "Phylum", "EffectCW",
  "CI95Lower", "CI95Upper", "QValue", "Highlighted"
)], 12)

relative_per_million <- sweep(counts, 2L, colSums(counts), "/") * 1e6
lefse_input <- SummarizedExperiment::SummarizedExperiment(
  assays = list(relative_abundance = relative_per_million),
  colData = S4Vectors::DataFrame(
    GROUP = metadata$Group, row.names = rownames(metadata)
  )
)
set.seed(20260733)
lefse_fit <- lefser::lefser(
  lefse_input,
  kruskal.threshold = 1,
  wilcox.threshold = 1,
  lda.threshold = 0,
  groupCol = "GROUP",
  checkAbundances = TRUE
)
lefse_scores <- setNames(lefse_fit$scores, lefse_fit$Names)
lefse_p <- vapply(rownames(counts), function(feature_id) {
  stats::wilcox.test(
    relative_per_million[feature_id, metadata$Group == "CW"],
    relative_per_million[feature_id, metadata$Group == "IW"],
    exact = FALSE
  )$p.value
}, numeric(1))
lefse_q <- stats::p.adjust(lefse_p, method = "BH")
descriptive_direction <- rowMeans(
  log2(counts[, metadata$Group == "CW", drop = FALSE] + 0.5)
) - rowMeans(
  log2(counts[, metadata$Group == "IW", drop = FALSE] + 0.5)
)
lefse_results <- data.frame(
  FeatureID = rownames(counts),
  LDAScore = sign(descriptive_direction) * abs(lefse_scores[rownames(counts)]),
  QValue = unname(lefse_q[rownames(counts)]),
  Reportable = feature_info[rownames(counts), "Reportable"],
  stringsAsFactors = FALSE
)
lefse_results$CladogramEligible <- with(
  lefse_results,
  Reportable & is.finite(LDAScore) & QValue < 0.05 & abs(LDAScore) >= 2
)

leaf_candidates <- lefse_results[lefse_results$CladogramEligible, ]
leaf_candidates <- leaf_candidates[order(-abs(leaf_candidates$LDAScore)), ]
leaf_candidates <- head(leaf_candidates, 24L)
leaf_info <- feature_info[leaf_candidates$FeatureID, , drop = FALSE]
leaf_direction <- setNames(
  ifelse(leaf_candidates$LDAScore >= 0, "Higher in CW", "Higher in IW"),
  leaf_candidates$FeatureID
)

edge_rows <- list()
node_rows <- list(data.frame(
  name = "Root", Label = "Bacteria", Rank = "Root",
  Evidence = "Context", stringsAsFactors = FALSE
))
for (feature_id in rownames(leaf_info)) {
  parts <- unlist(leaf_info[feature_id, rank_names[1:6]], use.names = FALSE)
  parent <- "Root"
  for (j in seq_along(parts)) {
    node_id <- paste(parts[seq_len(j)], collapse = ";")
    edge_rows[[length(edge_rows) + 1L]] <- data.frame(
      from = parent, to = node_id, stringsAsFactors = FALSE
    )
    node_rows[[length(node_rows) + 1L]] <- data.frame(
      name = node_id,
      Label = parts[[j]],
      Rank = rank_names[[j]],
      Evidence = if (j == 6L) leaf_direction[[feature_id]] else "Context",
      stringsAsFactors = FALSE
    )
    parent <- node_id
  }
}
cladogram_edges <- unique(do.call(rbind, edge_rows))
cladogram_nodes <- unique(do.call(rbind, node_rows))
cladogram_graph <- igraph::graph_from_data_frame(
  cladogram_edges, directed = TRUE, vertices = cladogram_nodes
)

stopifnot(
  igraph::is_dag(cladogram_graph),
  nrow(leaf_candidates) >= 8L,
  all(leaf_candidates$FeatureID %in% rownames(feature_info))
)

visual_audit <- data.frame(
  View = c("Volcano", "Cladogram", "Manhattan", "Forest"),
  PrimaryQuestion = c(
    "Where are direction and q-value extremes?",
    "Which taxonomic lineages contain supported leaves?",
    "Are signals concentrated in taxonomic blocks?",
    "How large and precise are selected effects?"
  ),
  RequiredFields = c(
    "effect, q, label policy",
    "full lineage, LDA, all-feature BH",
    "phylum, stable order, q",
    "effect, SE or CI, contrast"
  ),
  ProhibitedInference = c(
    "Distance from origin proves biology",
    "Internal nodes inherit significance",
    "Adjacent taxa are genetically linked",
    "CI proves absolute abundance or causality"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(visual_audit, caption = "One master table, four visual tasks")

volcano_data <- master_results[master_results$Reportable, ]
volcano_data$PlotQ <- pmax(volcano_data$QValue, 1e-12)
volcano_data$Status <- ifelse(
  volcano_data$Highlighted, volcano_data$Direction, "Not highlighted"
)
label_data <- head(
  volcano_data[order(volcano_data$QValue, -abs(volcano_data$EffectCW)), ],
  12L
)

volcano_plot <- ggplot(
  volcano_data,
  aes(x = EffectCW, y = -log10(PlotQ), colour = Status)
) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", colour = "#666666") +
  geom_vline(
    xintercept = c(-practical_threshold, practical_threshold),
    linetype = "22", colour = "#999999"
  ) +
  geom_point(size = 2.0, alpha = 0.82) +
  ggrepel::geom_text_repel(
    data = label_data, aes(label = DisplayTaxon),
    family = font_pub, size = 2.5, min.segment.length = 0,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Higher in CW" = pal_pub[["vermillion"]],
    "Higher in IW" = pal_pub[["blue"]],
    "Not highlighted" = "#B8B8B8"
  )) +
  labs(
    title = "Volcano plot: evidence and practical effect are separate",
    subtitle = "ANCOM-BC2, CW versus IW; labels are the 12 smallest-q genera",
    x = "Bias-corrected log fold change (positive = higher in CW)",
    y = expression(-log[10](BH~q)), colour = "Display status",
    caption = "Horizontal line: q = 0.05; vertical lines: |log fold change| = log(2)."
  ) +
  theme_pub(base_size = 8.5) +
  theme(legend.position = "bottom")

save_pub(volcano_plot, "figures/33-da-volcano", width = 170, height = 127)
volcano_plot

cladogram_plot <- ggraph::ggraph(
  cladogram_graph, layout = "dendrogram", circular = TRUE
) +
  ggraph::geom_edge_diagonal(colour = "#C7C7C7", alpha = 0.72, linewidth = 0.45) +
  ggraph::geom_node_point(
    aes(colour = Evidence, size = Rank == "Genus"), alpha = 0.92
  ) +
  ggraph::geom_node_text(
    aes(label = ifelse(Rank == "Genus", Label, "")),
    repel = TRUE, family = font_pub, size = 2.35,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Higher in CW" = pal_pub[["vermillion"]],
    "Higher in IW" = pal_pub[["blue"]],
    "Context" = "#A8A8A8"
  )) +
  scale_size_manual(values = c(`FALSE` = 1.4, `TRUE` = 3.1), guide = "none") +
  labs(
    title = "Audited LEfSe cladogram",
    subtitle = "Leaves require all-feature BH q < 0.05 and |LDA| >= 2",
    colour = "Genus direction",
    caption = paste0(
      "Internal nodes provide taxonomy context and are not called significant. ",
      "This is a taxonomy dendrogram, not a branch-length phylogeny."
    )
  ) +
  theme_void(base_family = font_pub) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(colour = "#4D4D4D", size = 9),
    plot.caption = element_text(colour = "#666666", hjust = 0, size = 7),
    legend.position = "bottom"
  )

save_pub(cladogram_plot, "figures/33-da-cladogram", width = 180, height = 180)
cladogram_plot

manhattan_data <- master_results[master_results$Reportable, ]
manhattan_data$PhylumDisplay <- ifelse(
  is.na(manhattan_data$Phylum) | manhattan_data$Phylum == "?",
  "Unresolved phylum", manhattan_data$Phylum
)
phylum_order <- names(sort(tapply(
  manhattan_data$QValue, manhattan_data$PhylumDisplay, min, na.rm = TRUE
)))
manhattan_data <- manhattan_data[order(
  match(manhattan_data$PhylumDisplay, phylum_order),
  manhattan_data$QValue, manhattan_data$DisplayTaxon
), ]
manhattan_data$Index <- seq_len(nrow(manhattan_data))
manhattan_data$PlotQ <- pmax(manhattan_data$QValue, 1e-12)
phylum_centres <- aggregate(Index ~ PhylumDisplay, manhattan_data, mean)
phylum_centres <- phylum_centres[match(phylum_order, phylum_centres$PhylumDisplay), ]
phylum_palette <- setNames(
  rep(c(
    pal_pub[["blue"]], pal_pub[["orange"]], pal_pub[["green"]],
    pal_pub[["purple"]], pal_pub[["sky"]], pal_pub[["vermillion"]],
    "#7A7A7A", "#B28E00"
  ), length.out = length(phylum_order)),
  phylum_order
)

manhattan_plot <- ggplot(
  manhattan_data,
  aes(x = Index, y = -log10(PlotQ), colour = PhylumDisplay)
) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", colour = "#666666") +
  geom_point(aes(shape = Highlighted), size = 2.0, alpha = 0.86, stroke = 0.65) +
  scale_colour_manual(values = phylum_palette) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  scale_x_continuous(
    breaks = phylum_centres$Index,
    labels = phylum_centres$PhylumDisplay,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Taxonomic Manhattan plot retains the full tested denominator",
    subtitle = "Open points are not highlighted; filled points pass q and practical-effect thresholds",
    x = "Phylum block", y = expression(-log[10](BH~q)),
    colour = "Phylum",
    caption = "Horizontal order is taxonomic bookkeeping, not a genomic coordinate or phylogenetic distance."
  ) +
  theme_pub(base_size = 8.0) +
  theme(
    axis.text.x = element_text(angle = 32, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

save_pub(manhattan_plot, "figures/33-da-manhattan", width = 180, height = 117)
manhattan_plot

forest_data <- master_results[
  master_results$Reportable & master_results$StatisticallySignificant &
    is.finite(master_results$CI95Lower) & is.finite(master_results$CI95Upper),
]
forest_data <- head(
  forest_data[order(forest_data$QValue, -abs(forest_data$EffectCW)), ],
  16L
)
forest_data$DisplayTaxon <- factor(
  forest_data$DisplayTaxon,
  levels = rev(forest_data$DisplayTaxon)
)

forest_plot <- ggplot(
  forest_data,
  aes(x = EffectCW, y = DisplayTaxon, colour = Direction)
) +
  geom_vline(xintercept = 0, linetype = "22", colour = "#666666") +
  geom_errorbarh(
    aes(xmin = CI95Lower, xmax = CI95Upper),
    height = 0.22, linewidth = 0.65
  ) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(
    "Higher in CW" = pal_pub[["vermillion"]],
    "Higher in IW" = pal_pub[["blue"]]
  )) +
  labs(
    title = "Forest plot makes effect precision visible",
    subtitle = "Sixteen smallest-q reportable genera from the ANCOM-BC2 master table",
    x = "Bias-corrected log fold change with Wald 95% CI",
    y = "Genus", colour = "Direction",
    caption = "Positive estimates indicate higher abundance in CW relative to IW; intervals are model based."
  ) +
  theme_pub(base_size = 8.2) +
  theme(legend.position = "bottom")

save_pub(forest_plot, "figures/33-da-forest", width = 170, height = 147)
forest_plot

expected_bases <- c(
  "figures/33-da-volcano", "figures/33-da-cladogram",
  "figures/33-da-manhattan", "figures/33-da-forest"
)
expected_files <- as.vector(outer(
  expected_bases, c(".pdf", ".svg", ".png", ".tiff"), paste0
))
stopifnot(
  all(file.exists(expected_files)),
  nrow(volcano_data) == 79L,
  nrow(manhattan_data) == 79L,
  nrow(forest_data) == 16L,
  nrow(leaf_candidates) >= 8L,
  all(master_results$FeatureID %in% rownames(feature_info))
)
