# 共现网络基础：SparCC/SPIEC-EASI/ggClusterNet
# Run sequentially in a new working directory.
# Required packages: SpiecEasi, ggplot2, ggrepel, igraph, knitr, patchwork, ragg, svglite.

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
library(SpiecEasi)

set.seed(20260735)
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
metadata <- read_keyed_tsv("data/small/metadata.tsv")
otutab <- as.matrix(data.frame(
  lapply(otutab_raw, as.integer), check.names = FALSE
))
rownames(otutab) <- rownames(otutab_raw)
metadata$Group <- factor(metadata$Group, levels = c("CW", "IW", "TW"))

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
feature_info_all <- data.frame(
  FeatureID = feature_ids,
  as.data.frame(lineage_matrix, stringsAsFactors = FALSE),
  DisplayTaxon = labels,
  Lineage = lineages,
  stringsAsFactors = FALSE,
  row.names = feature_ids,
  check.names = FALSE
)

prevalence <- rowMeans(genus_counts > 0)
total_reads <- rowSums(genus_counts)
eligible <- prevalence >= 0.20 & total_reads >= 100
selection_score <- prevalence * log1p(total_reads)
selected_ids <- names(sort(selection_score[eligible], decreasing = TRUE))[1:40]
counts <- genus_counts[selected_ids, , drop = FALSE]
feature_info <- feature_info_all[selected_ids, , drop = FALSE]
sample_by_taxon <- t(counts)

stopifnot(
  nrow(otutab) == 13628L,
  ncol(otutab) == 90L,
  sum(otutab) == 1619670L,
  nrow(sample_by_taxon) == 90L,
  ncol(sample_by_taxon) == 40L,
  all(rowMeans(counts > 0) >= 0.20),
  all(rowSums(counts) >= 100L),
  !any(grepl("Unclassified", rownames(counts), ignore.case = TRUE))
)

data.frame(
  Samples = nrow(sample_by_taxon),
  Genera = ncol(sample_by_taxon),
  CW = sum(metadata$Group == "CW"),
  IW = sum(metadata$Group == "IW"),
  TW = sum(metadata$Group == "TW"),
  RetainedNamedReadFraction = sum(counts) / sum(otutab)
)

sparcc_replicates <- 200L
set.seed(20260735)
sparcc_boot <- SpiecEasi::sparccboot(
  sample_by_taxon,
  R = sparcc_replicates,
  ncpus = 1
)
sparcc_inference <- SpiecEasi::pval.sparccboot(
  sparcc_boot, sided = "both"
)

# `pval.sparccboot()` returns the reference coefficients in `cors`.
# Reusing them keeps the effect estimates and empirical p-values tied to the
# same fitted object; a second `sparcc()` call can follow a different internal
# initialization path and should not be mixed with this bootstrap result.
sparcc_correlation <- matrix(
  0,
  nrow = length(selected_ids),
  ncol = length(selected_ids),
  dimnames = list(selected_ids, selected_ids)
)
sparcc_correlation[upper.tri(sparcc_correlation)] <- sparcc_inference$cors
sparcc_correlation <- sparcc_correlation + t(sparcc_correlation)
diag(sparcc_correlation) <- 1
upper_index <- which(upper.tri(sparcc_correlation), arr.ind = TRUE)

# Add-one correction prevents zero empirical p-values.
sparcc_unstable <- !is.finite(sparcc_inference$pvals)
sparcc_raw_p <- sparcc_inference$pvals
sparcc_raw_p[sparcc_unstable] <- 1
sparcc_p <- (
  sparcc_raw_p * sparcc_replicates + 1
) / (sparcc_replicates + 1)
sparcc_q <- stats::p.adjust(sparcc_p, method = "BH")
sparcc_all_pairs <- data.frame(
  From = rownames(sparcc_correlation)[upper_index[, "row"]],
  To = colnames(sparcc_correlation)[upper_index[, "col"]],
  Association = sparcc_inference$cors,
  EmpiricalP = sparcc_p,
  QValue = sparcc_q,
  BootstrapStable = !sparcc_unstable,
  stringsAsFactors = FALSE
)
sparcc_all_pairs$Selected <- with(
  sparcc_all_pairs,
  BootstrapStable & abs(Association) >= 0.30 & QValue < 0.05
)
sparcc_edges <- sparcc_all_pairs[sparcc_all_pairs$Selected, ]
sparcc_edges$Method <- "SparCC"

stopifnot(
  nrow(sparcc_all_pairs) == choose(40L, 2L),
  all(is.finite(sparcc_all_pairs$Association)),
  all(sparcc_all_pairs$EmpiricalP >= 1 / 201),
  all(sparcc_all_pairs$QValue >= 0 & sparcc_all_pairs$QValue <= 1),
  nrow(sparcc_edges) >= 10L
)

set.seed(20260735)
spiec_fit <- SpiecEasi::spiec.easi(
  sample_by_taxon,
  method = "mb",
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  pulsar.params = list(
    rep.num = 50,
    seed = 20260735L
  ),
  verbose = FALSE
)
spiec_adjacency <- as.matrix(SpiecEasi::getRefit(spiec_fit))
spiec_beta <- as.matrix(SpiecEasi::symBeta(
  SpiecEasi::getOptBeta(spiec_fit), mode = "maxabs"
))
spiec_edge_stability <- as.matrix(SpiecEasi::getOptMerge(spiec_fit))
dimnames(spiec_adjacency) <- list(selected_ids, selected_ids)
dimnames(spiec_beta) <- list(selected_ids, selected_ids)
dimnames(spiec_edge_stability) <- list(selected_ids, selected_ids)

spiec_index <- which(
  upper.tri(spiec_adjacency) & spiec_adjacency != 0,
  arr.ind = TRUE
)
spiec_edges <- data.frame(
  From = rownames(spiec_adjacency)[spiec_index[, "row"]],
  To = colnames(spiec_adjacency)[spiec_index[, "col"]],
  Association = spiec_beta[spiec_index],
  EdgeStability = spiec_edge_stability[spiec_index],
  Method = "SPIEC-EASI",
  stringsAsFactors = FALSE
)
spiec_instability <- as.numeric(SpiecEasi::getStability(spiec_fit))
spiec_lambda <- as.numeric(SpiecEasi::getOptLambda(spiec_fit))

stopifnot(
  nrow(spiec_adjacency) == 40L,
  nrow(spiec_edges) >= 5L,
  spiec_instability <= 0.05,
  all(spiec_edges$EdgeStability >= 0 & spiec_edges$EdgeStability <= 1),
  all(is.finite(spiec_edges$Association))
)

edge_key <- function(from, to) {
  paste(pmin(from, to), pmax(from, to), sep = " || ")
}
sparcc_edges$EdgeKey <- edge_key(sparcc_edges$From, sparcc_edges$To)
spiec_edges$EdgeKey <- edge_key(spiec_edges$From, spiec_edges$To)
all_edge_keys <- union(sparcc_edges$EdgeKey, spiec_edges$EdgeKey)

edge_overlap <- data.frame(
  EdgeKey = all_edge_keys,
  InSparCC = all_edge_keys %in% sparcc_edges$EdgeKey,
  InSPIECEASI = all_edge_keys %in% spiec_edges$EdgeKey,
  stringsAsFactors = FALSE
)
edge_overlap$Category <- with(
  edge_overlap,
  ifelse(
    InSparCC & InSPIECEASI, "Both methods",
    ifelse(InSparCC, "SparCC only", "SPIEC-EASI only")
  )
)
network_jaccard <- sum(edge_overlap$InSparCC & edge_overlap$InSPIECEASI) /
  nrow(edge_overlap)

common_keys <- intersect(sparcc_edges$EdgeKey, spiec_edges$EdgeKey)
common_edges <- data.frame(
  EdgeKey = common_keys,
  SparCC = sparcc_edges$Association[match(common_keys, sparcc_edges$EdgeKey)],
  SPIECEASI = spiec_edges$Association[match(common_keys, spiec_edges$EdgeKey)],
  stringsAsFactors = FALSE
)
common_edges$SignAgreement <- sign(common_edges$SparCC) ==
  sign(common_edges$SPIECEASI)

network_method_summary <- data.frame(
  Metric = c(
    "SparCC selected edges", "SPIEC-EASI selected edges",
    "Shared edges", "Edge-set Jaccard", "Shared-edge sign agreement",
    "SPIEC-EASI StARS instability", "SPIEC-EASI selected lambda"
  ),
  Value = c(
    nrow(sparcc_edges), nrow(spiec_edges), length(common_keys),
    network_jaccard,
    if (nrow(common_edges)) mean(common_edges$SignAgreement) else NA_real_,
    spiec_instability, spiec_lambda
  ),
  stringsAsFactors = FALSE
)
network_method_summary

log_counts <- log2(sample_by_taxon + 0.5)
clr_values <- log_counts - rowMeans(log_counts)
group_centered_clr <- clr_values
for (group_name in levels(metadata$Group)) {
  rows <- metadata$Group == group_name
  group_centered_clr[rows, ] <- sweep(
    clr_values[rows, , drop = FALSE],
    2L,
    colMeans(clr_values[rows, , drop = FALSE]),
    "-"
  )
}
group_centered_correlation <- stats::cor(
  group_centered_clr, method = "spearman"
)

sparcc_edges$GroupCenteredCLR <- mapply(
  function(from, to) group_centered_correlation[from, to],
  sparcc_edges$From, sparcc_edges$To
)
sparcc_edges$GroupRobust <- with(
  sparcc_edges,
  sign(Association) == sign(GroupCenteredCLR) & abs(GroupCenteredCLR) >= 0.20
)
group_audit_summary <- data.frame(
  SparCCEdges = nrow(sparcc_edges),
  SameDirectionAfterGroupCentering = sum(
    sign(sparcc_edges$Association) == sign(sparcc_edges$GroupCenteredCLR)
  ),
  GroupRobustEdges = sum(sparcc_edges$GroupRobust),
  GroupRobustFraction = mean(sparcc_edges$GroupRobust),
  stringsAsFactors = FALSE
)
group_audit_summary

network_graph <- igraph::graph_from_data_frame(
  spiec_edges[, c("From", "To")],
  directed = FALSE,
  vertices = data.frame(
    name = selected_ids,
    DisplayTaxon = feature_info[selected_ids, "DisplayTaxon"],
    Phylum = feature_info[selected_ids, "Phylum"],
    stringsAsFactors = FALSE
  )
)
set.seed(20260735)
network_layout <- igraph::layout_with_fr(
  network_graph, niter = 1000, grid = "nogrid"
)
node_plot_data <- data.frame(
  FeatureID = igraph::V(network_graph)$name,
  X = network_layout[, 1L],
  Y = network_layout[, 2L],
  DisplayTaxon = igraph::vertex_attr(network_graph, "DisplayTaxon"),
  Phylum = igraph::vertex_attr(network_graph, "Phylum"),
  Degree = igraph::degree(network_graph),
  stringsAsFactors = FALSE
)
top_phyla <- names(head(sort(table(node_plot_data$Phylum), decreasing = TRUE), 5L))
node_plot_data$PhylumDisplay <- ifelse(
  node_plot_data$Phylum %in% top_phyla,
  node_plot_data$Phylum, "Other phyla"
)
label_nodes <- head(
  node_plot_data[order(-node_plot_data$Degree, node_plot_data$DisplayTaxon), ],
  10L
)

edge_plot_data <- spiec_edges
edge_plot_data$X <- node_plot_data$X[match(edge_plot_data$From, node_plot_data$FeatureID)]
edge_plot_data$Y <- node_plot_data$Y[match(edge_plot_data$From, node_plot_data$FeatureID)]
edge_plot_data$XEnd <- node_plot_data$X[match(edge_plot_data$To, node_plot_data$FeatureID)]
edge_plot_data$YEnd <- node_plot_data$Y[match(edge_plot_data$To, node_plot_data$FeatureID)]
edge_plot_data$Direction <- ifelse(
  edge_plot_data$Association >= 0, "Positive", "Negative"
)

stopifnot(
  all(is.finite(node_plot_data$X)),
  all(is.finite(edge_plot_data$X)),
  sum(node_plot_data$Degree) == 2L * nrow(spiec_edges)
)

network_audit <- data.frame(
  Dimension = c(
    "Composition", "Feature filter", "Uncertainty", "Covariates",
    "Method", "Topology", "Biological claim"
  ),
  FragilePractice = c(
    "Pearson on percentages",
    "Choose nodes after seeing edges",
    "Threshold r without resampling",
    "Pool groups without audit",
    "Treat every method as interchangeable",
    "Call high degree a keystone",
    "Call an edge an interaction"
  ),
  UpgradedContract = c(
    "Use composition-aware inference",
    "Group-blind prevalence/abundance contract",
    "Permutation/bootstrapping and StARS ledger",
    "Group-centered sensitivity or covariate model",
    "Compare marginal and conditional estimands",
    "Report degree with selection uncertainty",
    "Validate with perturbation, culture or mechanism"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(network_audit, caption = "Microbial association-network audit")

taxon_order <- hclust(as.dist(1 - abs(sparcc_correlation)))$order
ordered_ids <- selected_ids[taxon_order]
matrix_plot_data <- as.data.frame(as.table(sparcc_correlation))
names(matrix_plot_data) <- c("From", "To", "Correlation")
matrix_plot_data$From <- factor(matrix_plot_data$From, levels = rev(ordered_ids))
matrix_plot_data$To <- factor(matrix_plot_data$To, levels = ordered_ids)
matrix_plot_data$SPIECEdge <- edge_key(
  as.character(matrix_plot_data$From), as.character(matrix_plot_data$To)
) %in% spiec_edges$EdgeKey &
  as.character(matrix_plot_data$From) != as.character(matrix_plot_data$To)

display_labels <- setNames(feature_info$DisplayTaxon, rownames(feature_info))
network_matrix_plot <- ggplot(
  matrix_plot_data,
  aes(x = To, y = From, fill = Correlation)
) +
  geom_tile() +
  geom_point(
    data = matrix_plot_data[matrix_plot_data$SPIECEdge, ],
    shape = 21, fill = NA, colour = "#111111", stroke = 0.45, size = 1.25,
    inherit.aes = FALSE,
    aes(x = To, y = From)
  ) +
  scale_fill_gradient2(
    low = pal_pub[["blue"]], mid = "white", high = pal_pub[["vermillion"]],
    midpoint = 0, limits = c(-1, 1)
  ) +
  scale_x_discrete(labels = display_labels) +
  scale_y_discrete(labels = display_labels) +
  coord_equal() +
  labs(
    title = "Marginal and conditional networks are different views",
    subtitle = "SparCC correlation fill; black rings mark SPIEC-EASI selected edges",
    x = NULL, y = NULL, fill = "SparCC r",
    caption = "Rows and columns were clustered by 1 - |SparCC r|; ordering is descriptive."
  ) +
  theme_pub(base_size = 6.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, size = 5.7),
    axis.text.y = element_text(size = 5.7),
    legend.position = "bottom"
  )

save_pub(network_matrix_plot, "figures/35-network-matrix", width = 180, height = 165)
network_matrix_plot

phylum_colours <- setNames(
  c(
    pal_pub[["blue"]], pal_pub[["orange"]], pal_pub[["green"]],
    pal_pub[["purple"]], pal_pub[["sky"]], "#8A8A8A"
  )[seq_along(unique(node_plot_data$PhylumDisplay))],
  unique(node_plot_data$PhylumDisplay)
)

spiec_network_plot <- ggplot() +
  geom_segment(
    data = edge_plot_data,
    aes(
      x = X, y = Y, xend = XEnd, yend = YEnd,
      colour = Direction,
      linewidth = abs(Association),
      alpha = EdgeStability
    ),
    lineend = "round"
  ) +
  geom_point(
    data = node_plot_data,
    aes(x = X, y = Y, fill = PhylumDisplay, size = Degree),
    shape = 21, colour = "white", stroke = 0.45
  ) +
  ggrepel::geom_text_repel(
    data = label_nodes,
    aes(x = X, y = Y, label = DisplayTaxon),
    family = font_pub, size = 2.45, min.segment.length = 0,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    Positive = pal_pub[["vermillion"]], Negative = pal_pub[["blue"]]
  )) +
  scale_fill_manual(values = phylum_colours) +
  scale_linewidth_continuous(range = c(0.35, 1.55), guide = "none") +
  scale_alpha_continuous(range = c(0.35, 0.92), limits = c(0, 1), guide = "none") +
  scale_size_continuous(range = c(2.2, 7.2)) +
  coord_equal() +
  labs(
    title = "SPIEC-EASI sparse conditional-association network",
    subtitle = sprintf(
      "MB neighborhood selection; %d edges; StARS instability %.3f",
      nrow(spiec_edges), spiec_instability
    ),
    x = NULL, y = NULL, colour = "Association sign",
    fill = "Phylum", size = "Degree",
    caption = paste0(
      "Edge width is |symmetrized MB coefficient| and alpha is selection stability. ",
      "Labels mark the ten highest-degree nodes."
    )
  ) +
  theme_void(base_family = font_pub) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(colour = "#4D4D4D", size = 9),
    plot.caption = element_text(colour = "#666666", hjust = 0, size = 7),
    legend.position = "bottom", legend.box = "vertical"
  )

save_pub(spiec_network_plot, "figures/35-spiec-network", width = 180, height = 147)
spiec_network_plot

overlap_counts <- as.data.frame(table(edge_overlap$Category), stringsAsFactors = FALSE)
names(overlap_counts) <- c("Category", "Edges")
overlap_counts$Category <- factor(
  overlap_counts$Category,
  levels = c("SparCC only", "Both methods", "SPIEC-EASI only")
)

overlap_panel <- ggplot(
  overlap_counts,
  aes(x = Category, y = Edges, fill = Category)
) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = Edges), vjust = -0.35, family = font_pub, size = 3) +
  scale_fill_manual(values = c(
    "SparCC only" = pal_pub[["orange"]],
    "Both methods" = pal_pub[["purple"]],
    "SPIEC-EASI only" = pal_pub[["green"]]
  )) +
  scale_y_continuous(
    limits = c(0, max(overlap_counts$Edges) * 1.20 + 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "A. Method-specific edge sets",
    subtitle = sprintf("Edge-set Jaccard = %.2f", network_jaccard),
    x = NULL, y = "Selected edges"
  ) +
  theme_pub(base_size = 8.1) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))

group_counts <- data.frame(
  Category = c("Group-robust", "Group-sensitive"),
  Edges = c(
    sum(sparcc_edges$GroupRobust),
    sum(!sparcc_edges$GroupRobust)
  )
)
group_panel <- ggplot(
  group_counts,
  aes(x = Category, y = Edges, fill = Category)
) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = Edges), vjust = -0.35, family = font_pub, size = 3) +
  scale_fill_manual(values = c(
    "Group-robust" = pal_pub[["green"]],
    "Group-sensitive" = pal_pub[["grey"]]
  )) +
  scale_y_continuous(
    limits = c(0, max(group_counts$Edges) * 1.20 + 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "B. SparCC group-centering audit",
    subtitle = "Same sign and |group-centered CLR Spearman| >= 0.20",
    x = NULL, y = "SparCC edges"
  ) +
  theme_pub(base_size = 8.1) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))

edge_audit_plot <- patchwork::wrap_plots(overlap_panel, group_panel, nrow = 1L) +
  patchwork::plot_annotation(
    title = "A network figure needs an edge audit ledger",
    subtitle = "Method overlap and group sensitivity answer different questions",
    caption = "Neither panel validates ecological interaction or causality."
  )
save_pub(edge_audit_plot, "figures/35-network-edge-audit", width = 180, height = 115)
edge_audit_plot

expected_bases <- c(
  "figures/35-network-matrix", "figures/35-spiec-network",
  "figures/35-network-edge-audit"
)
expected_files <- as.vector(outer(
  expected_bases, c(".pdf", ".svg", ".png", ".tiff"), paste0
))
stopifnot(
  all(file.exists(expected_files)),
  sparcc_replicates == 200L,
  nrow(sparcc_all_pairs) == 780L,
  nrow(spiec_edges) == sum(spiec_adjacency[upper.tri(spiec_adjacency)] != 0),
  spiec_instability <= 0.05,
  group_audit_summary$GroupRobustEdges <= nrow(sparcc_edges),
  sum(node_plot_data$Degree) == 2L * nrow(spiec_edges)
)
