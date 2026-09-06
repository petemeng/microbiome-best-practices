# 网络进阶：鲁棒性、拓扑角色与组间比较
# Run sequentially in a new working directory.
# Required packages: SpiecEasi, ggplot2, ggrepel, igraph, knitr, patchwork, ragg, scales, svglite.

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

set.seed(20260736)
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
lineage_matrix <- do.call(
  rbind, strsplit(rownames(genus_counts), ";", fixed = TRUE)
)
colnames(lineage_matrix) <- rank_names[1:6]
display_taxon <- lineage_matrix[, "Genus"]
feature_ids <- make.unique(paste0("Genus:", display_taxon))
rownames(genus_counts) <- feature_ids
feature_info_all <- data.frame(
  FeatureID = feature_ids,
  as.data.frame(lineage_matrix, stringsAsFactors = FALSE),
  DisplayTaxon = display_taxon,
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
  identical(rownames(sample_by_taxon), rownames(metadata)),
  ncol(sample_by_taxon) == 40L,
  all(rowMeans(counts > 0) >= 0.20),
  all(rowSums(counts) >= 100L)
)

data.frame(
  Samples = nrow(sample_by_taxon),
  Nodes = ncol(sample_by_taxon),
  CW = sum(metadata$Group == "CW"),
  IW = sum(metadata$Group == "IW"),
  TW = sum(metadata$Group == "TW"),
  RetainedReadFraction = sum(counts) / sum(otutab)
)

set.seed(20260736)
sparcc_cor <- SpiecEasi::sparcc(sample_by_taxon)$Cor
dimnames(sparcc_cor) <- list(selected_ids, selected_ids)
pair_index <- which(upper.tri(sparcc_cor), arr.ind = TRUE)
sparcc_pairs <- data.frame(
  From = rownames(sparcc_cor)[pair_index[, "row"]],
  To = colnames(sparcc_cor)[pair_index[, "col"]],
  Association = sparcc_cor[pair_index],
  stringsAsFactors = FALSE
)
sparcc_pairs$Selected <- abs(sparcc_pairs$Association) >= 0.30

# Non-parametric sample bootstrap audits threshold and sign stability.
sparcc_bootstrap_replicates <- 100L
set.seed(20260736)
bootstrap_indices <- replicate(
  sparcc_bootstrap_replicates,
  sample(seq_len(nrow(sample_by_taxon)), replace = TRUE),
  simplify = FALSE
)
bootstrap_association <- vapply(
  bootstrap_indices,
  function(index) {
    fit <- SpiecEasi::sparcc(sample_by_taxon[index, , drop = FALSE])$Cor
    fit[upper.tri(fit)]
  },
  numeric(nrow(sparcc_pairs))
)
sparcc_pairs$SelectionFrequency <- rowMeans(abs(bootstrap_association) >= 0.30)
sparcc_pairs$SignConsistency <- rowMeans(
  sign(bootstrap_association) == sign(sparcc_pairs$Association)
)
sparcc_edges <- sparcc_pairs[sparcc_pairs$Selected, , drop = FALSE]

set.seed(20260736)
spiec_fit <- SpiecEasi::spiec.easi(
  sample_by_taxon,
  method = "mb",
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  pulsar.params = list(rep.num = 50, seed = 20260736L),
  verbose = FALSE
)
spiec_adjacency <- as.matrix(SpiecEasi::getRefit(spiec_fit))
spiec_beta <- as.matrix(SpiecEasi::symBeta(
  SpiecEasi::getOptBeta(spiec_fit), mode = "maxabs"
))
dimnames(spiec_adjacency) <- dimnames(spiec_beta) <- list(
  selected_ids, selected_ids
)
spiec_index <- which(
  upper.tri(spiec_adjacency) & spiec_adjacency != 0,
  arr.ind = TRUE
)
spiec_edges <- data.frame(
  From = rownames(spiec_adjacency)[spiec_index[, "row"]],
  To = colnames(spiec_adjacency)[spiec_index[, "col"]],
  Association = spiec_beta[spiec_index],
  stringsAsFactors = FALSE
)
spiec_instability <- as.numeric(SpiecEasi::getStability(spiec_fit))

make_graph <- function(edges) {
  graph <- igraph::graph_from_data_frame(
    edges[, c("From", "To", "Association")],
    directed = FALSE,
    vertices = data.frame(
      name = selected_ids,
      DisplayTaxon = feature_info[selected_ids, "DisplayTaxon"],
      stringsAsFactors = FALSE
    )
  )
  graph <- igraph::set_edge_attr(
    graph, "weight_abs", value = abs(igraph::E(graph)$Association)
  )
  graph
}
sparcc_graph <- make_graph(sparcc_edges)
spiec_graph <- make_graph(spiec_edges)

stopifnot(
  nrow(sparcc_pairs) == choose(40L, 2L),
  nrow(sparcc_edges) >= 20L,
  nrow(spiec_edges) >= 5L,
  spiec_instability <= 0.05,
  all(sparcc_pairs$SelectionFrequency >= 0 & sparcc_pairs$SelectionFrequency <= 1),
  all(sparcc_pairs$SignConsistency >= 0 & sparcc_pairs$SignConsistency <= 1)
)

calculate_roles <- function(graph, method) {
  set.seed(20260736)
  modules <- igraph::cluster_louvain(
    graph, weights = igraph::E(graph)$weight_abs
  )
  membership <- igraph::membership(modules)
  weight_matrix <- as.matrix(igraph::as_adjacency_matrix(
    graph, attr = "weight_abs", sparse = FALSE
  ))
  zi <- pi <- numeric(igraph::vcount(graph))
  names(zi) <- names(pi) <- igraph::V(graph)$name
  module_ids <- sort(unique(membership))

  for (node in igraph::V(graph)$name) {
    own_module <- membership[[node]]
    own_nodes <- names(membership)[membership == own_module]
    within_degree <- rowSums(
      weight_matrix[own_nodes, own_nodes, drop = FALSE]
    )
    zi[[node]] <- if (length(own_nodes) > 1L && stats::sd(within_degree) > 0) {
      (sum(weight_matrix[node, own_nodes]) - mean(within_degree)) /
        stats::sd(within_degree)
    } else {
      0
    }
    total_degree <- sum(weight_matrix[node, ])
    by_module <- vapply(
      module_ids,
      function(module_id) sum(
        weight_matrix[node, names(membership)[membership == module_id]]
      ),
      numeric(1)
    )
    pi[[node]] <- if (total_degree > 0) {
      1 - sum((by_module / total_degree)^2)
    } else {
      0
    }
  }

  out <- data.frame(
    FeatureID = names(zi),
    DisplayTaxon = igraph::vertex_attr(
      graph, "DisplayTaxon", index = names(zi)
    ),
    Method = method,
    Module = unname(membership[names(zi)]),
    Degree = unname(igraph::degree(graph)[names(zi)]),
    Strength = unname(igraph::strength(
      graph, vids = names(zi), weights = igraph::E(graph)$weight_abs
    )),
    Zi = unname(zi),
    Pi = unname(pi),
    stringsAsFactors = FALSE
  )
  out$Role <- with(
    out,
    ifelse(
      Zi > 2.5 & Pi > 0.62, "Network hub",
      ifelse(
        Zi > 2.5, "Module hub",
        ifelse(Pi > 0.62, "Connector", "Peripheral")
      )
    )
  )
  attr(out, "modularity") <- igraph::modularity(modules)
  out
}

sparcc_roles <- calculate_roles(sparcc_graph, "SparCC")
spiec_roles <- calculate_roles(spiec_graph, "SPIEC-EASI")
all_roles <- rbind(sparcc_roles, spiec_roles)

incident_stability <- vapply(
  sparcc_roles$FeatureID,
  function(node) {
    values <- sparcc_edges$SelectionFrequency[
      sparcc_edges$From == node | sparcc_edges$To == node
    ]
    if (length(values)) stats::median(values) else 0
  },
  numeric(1)
)
sparcc_roles$MedianIncidentStability <- incident_stability
sparcc_roles$TopologicalCandidate <- with(
  sparcc_roles,
  Role != "Peripheral" & MedianIncidentStability >= 0.70
)

role_summary <- as.data.frame(
  table(all_roles$Method, all_roles$Role), stringsAsFactors = FALSE
)
names(role_summary) <- c("Method", "Role", "Nodes")
role_summary <- role_summary[role_summary$Nodes > 0, ]
role_summary

giant_component <- function(graph) {
  component <- igraph::components(graph)
  igraph::induced_subgraph(
    graph,
    igraph::V(graph)[component$membership == which.max(component$csize)]
  )
}
largest_component_fraction <- function(graph, initial_n) {
  if (igraph::vcount(graph) == 0L) return(0)
  max(igraph::components(graph)$csize) / initial_n
}
adaptive_attack <- function(graph, strategy) {
  graph <- giant_component(graph)
  initial_n <- igraph::vcount(graph)
  result <- data.frame(
    RemovedFraction = seq(0, 1, length.out = initial_n + 1L),
    LCCFraction = numeric(initial_n + 1L),
    Strategy = strategy,
    stringsAsFactors = FALSE
  )
  result$LCCFraction[[1L]] <- 1
  for (step in seq_len(initial_n)) {
    score <- if (strategy == "Degree attack") {
      igraph::degree(graph)
    } else {
      igraph::betweenness(
        graph, directed = FALSE, normalized = FALSE, weights = NA
      )
    }
    victim <- names(sort(score, decreasing = TRUE))[1L]
    graph <- igraph::delete_vertices(graph, victim)
    result$LCCFraction[[step + 1L]] <- largest_component_fraction(
      graph, initial_n
    )
  }
  result
}
random_attacks <- function(graph, replicates = 500L) {
  graph <- giant_component(graph)
  initial_n <- igraph::vcount(graph)
  set.seed(20260736)
  curves <- vapply(
    seq_len(replicates),
    function(replicate_id) {
      order <- sample(igraph::V(graph)$name)
      working <- graph
      curve <- numeric(initial_n + 1L)
      curve[[1L]] <- 1
      for (step in seq_along(order)) {
        working <- igraph::delete_vertices(working, order[[step]])
        curve[[step + 1L]] <- largest_component_fraction(
          working, initial_n
        )
      }
      curve
    },
    numeric(initial_n + 1L)
  )
  data.frame(
    RemovedFraction = seq(0, 1, length.out = initial_n + 1L),
    Median = apply(curves, 1L, stats::median),
    Lower = apply(curves, 1L, stats::quantile, probs = 0.025),
    Upper = apply(curves, 1L, stats::quantile, probs = 0.975),
    stringsAsFactors = FALSE
  )
}
curve_auc <- function(y) mean((head(y, -1L) + tail(y, -1L)) / 2)

sparcc_degree_attack <- adaptive_attack(sparcc_graph, "Degree attack")
sparcc_betweenness_attack <- adaptive_attack(
  sparcc_graph, "Betweenness attack"
)
sparcc_random_attack <- random_attacks(sparcc_graph, replicates = 500L)

network_attack_summary <- function(graph, method) {
  degree_curve <- adaptive_attack(graph, "Degree attack")
  betweenness_curve <- adaptive_attack(graph, "Betweenness attack")
  random_curve <- random_attacks(graph, replicates = 500L)
  data.frame(
    Method = method,
    InitialGiantNodes = igraph::vcount(giant_component(graph)),
    RandomAUC = curve_auc(random_curve$Median),
    DegreeAttackAUC = curve_auc(degree_curve$LCCFraction),
    BetweennessAttackAUC = curve_auc(betweenness_curve$LCCFraction),
    stringsAsFactors = FALSE
  )
}
attack_summary <- rbind(
  network_attack_summary(sparcc_graph, "SparCC"),
  network_attack_summary(spiec_graph, "SPIEC-EASI")
)
attack_summary

clr <- log2(sample_by_taxon + 0.5)
clr <- clr - rowMeans(clr)
edge_keys <- apply(
  utils::combn(colnames(clr), 2L), 2L,
  function(x) paste(x, collapse = " || ")
)
group_edge_count <- nrow(spiec_edges)

select_group_edges <- function(rows, edge_count = group_edge_count) {
  correlation <- stats::cor(
    clr[rows, , drop = FALSE], method = "spearman"
  )
  values <- correlation[upper.tri(correlation)]
  edge_keys[order(abs(values), decreasing = TRUE)[seq_len(edge_count)]]
}
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

group_names <- levels(metadata$Group)
group_edge_sets <- setNames(
  lapply(group_names, function(group_name) {
    select_group_edges(metadata$Group == group_name)
  }),
  group_names
)
group_pairs <- utils::combn(group_names, 2L)
group_jaccard <- data.frame(
  Comparison = apply(group_pairs, 2L, paste, collapse = " vs "),
  Jaccard = apply(group_pairs, 2L, function(x) {
    jaccard(group_edge_sets[[x[[1L]]]], group_edge_sets[[x[[2L]]]])
  }),
  stringsAsFactors = FALSE
)
observed_mean_jaccard <- mean(group_jaccard$Jaccard)

set.seed(20260736)
permuted_mean_jaccard <- replicate(499L, {
  permuted_group <- sample(metadata$Group)
  permuted_sets <- lapply(group_names, function(group_name) {
    select_group_edges(permuted_group == group_name)
  })
  mean(apply(utils::combn(seq_along(group_names), 2L), 2L, function(x) {
    jaccard(permuted_sets[[x[[1L]]]], permuted_sets[[x[[2L]]]])
  }))
})
rewiring_p <- (
  1 + sum(permuted_mean_jaccard <= observed_mean_jaccard)
) / 500

set.seed(20260736)
group_bootstrap_frequency <- lapply(group_names, function(group_name) {
  rows <- which(metadata$Group == group_name)
  selected <- replicate(
    500L,
    select_group_edges(sample(rows, replace = TRUE)),
    simplify = FALSE
  )
  frequency <- table(unlist(selected, use.names = FALSE)) / 500
  data.frame(
    Group = group_name,
    EdgeKey = names(frequency),
    SelectionFrequency = as.numeric(frequency),
    stringsAsFactors = FALSE
  )
})
group_bootstrap_frequency <- do.call(rbind, group_bootstrap_frequency)
stable_group_edges <- aggregate(
  SelectionFrequency ~ Group,
  data = group_bootstrap_frequency,
  FUN = function(x) sum(x >= 0.70)
)
names(stable_group_edges)[[2L]] <- "StableEdges"

stopifnot(
  all(vapply(group_edge_sets, length, integer(1)) == group_edge_count),
  rewiring_p >= 1 / 500 && rewiring_p <= 1,
  all(stable_group_edges$StableEdges <= group_edge_count)
)

data.frame(
  MeanPairwiseJaccard = observed_mean_jaccard,
  PermutationP = rewiring_p,
  NullMedian = stats::median(permuted_mean_jaccard),
  NullLower = stats::quantile(permuted_mean_jaccard, 0.025),
  NullUpper = stats::quantile(permuted_mean_jaccard, 0.975)
)

network_robustness_audit <- data.frame(
  Dimension = c(
    "Node universe", "Edge density", "Edge uncertainty", "Role threshold",
    "Attack curve", "Group comparison", "Keystone claim"
  ),
  FragilePractice = c(
    "Select taxa separately by group",
    "Compare raw topology at different densities",
    "Treat one thresholded edge set as fixed",
    "Use Zi/Pi labels as biological truth",
    "Equate graph fragmentation with resilience",
    "Read rewiring as treatment effect",
    "Call the highest-degree node a keystone"
  ),
  UpgradedContract = c(
    "One group-blind 40-node set",
    "Same samples per group and same edge count",
    "Sample bootstrap selection frequency",
    "Report method-dependent role ledger",
    "State graph-only estimand and attack rule",
    "Permutation plus design/confounding boundary",
    "Candidate label plus perturbation requirement"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(
  network_robustness_audit,
  caption = "Audit contract for advanced microbial-network claims"
)

role_colours <- c(
  "Peripheral" = pal_pub[["grey"]],
  "Connector" = pal_pub[["orange"]],
  "Module hub" = pal_pub[["blue"]],
  "Network hub" = pal_pub[["vermillion"]]
)
label_roles <- sparcc_roles[
  sparcc_roles$Role != "Peripheral" |
    rank(-sparcc_roles$Zi, ties.method = "first") <= 4L,
]
role_plot <- ggplot(
  sparcc_roles,
  aes(x = Pi, y = Zi, fill = Role, size = Degree)
) +
  geom_vline(xintercept = 0.62, linetype = 2, colour = "#555555") +
  geom_hline(yintercept = 2.5, linetype = 2, colour = "#555555") +
  geom_point(shape = 21, colour = "white", stroke = 0.45, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = label_roles,
    aes(label = DisplayTaxon),
    family = font_pub, size = 2.6, min.segment.length = 0,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_fill_manual(values = role_colours, drop = FALSE) +
  scale_size_continuous(range = c(2.4, 7.5)) +
  coord_cartesian(xlim = c(0, max(0.72, max(sparcc_roles$Pi) * 1.08))) +
  labs(
    title = "Zi-Pi roles describe one inferred topology",
    subtitle = sprintf(
      "SparCC |r| >= 0.30; %d nodes, %d edges, %d modules",
      igraph::vcount(sparcc_graph), igraph::ecount(sparcc_graph),
      length(unique(sparcc_roles$Module))
    ),
    x = "Participation coefficient (Pi)",
    y = "Within-module degree z-score (Zi)",
    fill = "Topological role", size = "Degree",
    caption = "Dashed lines are the classic Pi = 0.62 and Zi = 2.5 role thresholds; roles are not biological validation."
  ) +
  theme_pub(base_size = 8.5) +
  theme(legend.position = "bottom")
save_pub(role_plot, "figures/36-role-cartography", width = 180, height = 145)
role_plot

attack_lines <- rbind(
  transform(
    sparcc_degree_attack[, c("RemovedFraction", "LCCFraction")],
    Strategy = "Degree attack"
  ),
  transform(
    sparcc_betweenness_attack[, c("RemovedFraction", "LCCFraction")],
    Strategy = "Betweenness attack"
  )
)
attack_colours <- c(
  "Random failure" = pal_pub[["grey"]],
  "Degree attack" = pal_pub[["vermillion"]],
  "Betweenness attack" = pal_pub[["blue"]]
)
attack_plot <- ggplot() +
  geom_ribbon(
    data = sparcc_random_attack,
    aes(x = RemovedFraction, ymin = Lower, ymax = Upper),
    fill = pal_pub[["grey"]], alpha = 0.20
  ) +
  geom_line(
    data = sparcc_random_attack,
    aes(x = RemovedFraction, y = Median, colour = "Random failure"),
    linewidth = 0.9
  ) +
  geom_line(
    data = attack_lines,
    aes(x = RemovedFraction, y = LCCFraction, colour = Strategy),
    linewidth = 0.95
  ) +
  scale_colour_manual(values = attack_colours) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  labs(
    title = "Targeted node removal fragments the inferred graph faster",
    subtitle = sprintf(
      "AUC: random %.3f; degree %.3f; betweenness %.3f",
      attack_summary$RandomAUC[attack_summary$Method == "SparCC"],
      attack_summary$DegreeAttackAUC[attack_summary$Method == "SparCC"],
      attack_summary$BetweennessAttackAUC[attack_summary$Method == "SparCC"]
    ),
    x = "Nodes removed", y = "Largest connected component / initial nodes",
    colour = "Removal rule",
    caption = "Ribbon: 95% interval across 500 random removal orders. Targeted scores are recalculated after every removal."
  ) +
  theme_pub(base_size = 8.8) +
  theme(legend.position = "bottom")
save_pub(
  attack_plot, "figures/36-attack-robustness", width = 180, height = 122
)
attack_plot

group_jaccard$Comparison <- factor(
  group_jaccard$Comparison, levels = group_jaccard$Comparison
)
jaccard_panel <- ggplot(
  group_jaccard,
  aes(x = Comparison, y = Jaccard, fill = Comparison)
) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.3f", Jaccard)),
    vjust = -0.35, family = font_pub, size = 3
  ) +
  scale_fill_manual(values = unname(pal_pub[c("blue", "orange", "green")])) +
  scale_y_continuous(limits = c(0, max(group_jaccard$Jaccard) * 1.28)) +
  labs(
    title = "A. Density-matched edge overlap",
    subtitle = sprintf("%d nodes and %d edges in every group", 40L, group_edge_count),
    x = NULL, y = "Edge-set Jaccard"
  ) +
  theme_pub(base_size = 8.1) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))

permutation_data <- data.frame(MeanJaccard = permuted_mean_jaccard)
null_panel <- ggplot(permutation_data, aes(x = MeanJaccard)) +
  geom_histogram(
    bins = 28, fill = pal_pub[["sky"]], colour = "white", linewidth = 0.25
  ) +
  geom_vline(
    xintercept = observed_mean_jaccard,
    colour = pal_pub[["vermillion"]], linewidth = 1.05
  ) +
  annotate(
    "text", x = observed_mean_jaccard, y = Inf,
    label = sprintf("Observed = %.3f\nP = %.3f", observed_mean_jaccard, rewiring_p),
    hjust = -0.08, vjust = 1.25, family = font_pub, size = 2.8
  ) +
  labs(
    title = "B. Group-label permutation null",
    subtitle = "Lower overlap means stronger apparent rewiring",
    x = "Mean pairwise edge Jaccard", y = "Permutations"
  ) +
  theme_pub(base_size = 8.1)

group_rewiring_plot <- patchwork::wrap_plots(
  jaccard_panel, null_panel, nrow = 1L, widths = c(0.9, 1.1)
) +
  patchwork::plot_annotation(
    title = "Group-associated rewiring survives a matched comparison",
    caption = "The observational groups remain confounded with geography, salinity and other environmental variables."
  )
save_pub(
  group_rewiring_plot, "figures/36-group-rewiring",
  width = 180, height = 122
)
group_rewiring_plot

topology_metrics <- data.frame(
  Method = rep(c("SparCC", "SPIEC-EASI"), each = 4L),
  Metric = rep(
    c("Edges", "Giant component", "Modules", "Non-peripheral roles"),
    times = 2L
  ),
  Value = c(
    igraph::ecount(sparcc_graph),
    max(igraph::components(sparcc_graph)$csize),
    length(unique(sparcc_roles$Module)),
    sum(sparcc_roles$Role != "Peripheral"),
    igraph::ecount(spiec_graph),
    max(igraph::components(spiec_graph)$csize),
    length(unique(spiec_roles$Module)),
    sum(spiec_roles$Role != "Peripheral")
  ),
  stringsAsFactors = FALSE
)
topology_metrics$Metric <- factor(
  topology_metrics$Metric,
  levels = c("Edges", "Giant component", "Modules", "Non-peripheral roles")
)
topology_panel <- ggplot(
  topology_metrics,
  aes(x = Method, y = Value, fill = Method)
) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(
    aes(label = Value), vjust = -0.35, family = font_pub, size = 2.8
  ) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1L) +
  scale_fill_manual(values = c(
    "SparCC" = pal_pub[["orange"]],
    "SPIEC-EASI" = pal_pub[["purple"]]
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Topology is an inference-method-dependent result",
    subtitle = sprintf(
      "Same 90 samples and 40 nodes; SPIEC-EASI StARS instability = %.3f",
      spiec_instability
    ),
    x = NULL, y = NULL,
    caption = "Marginal and sparse conditional networks have different estimands; neither is a validated interaction map."
  ) +
  theme_pub(base_size = 7.8) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.spacing.x = grid::unit(4, "mm")
  )
save_pub(
  topology_panel, "figures/36-topology-sensitivity",
  width = 180, height = 132
)
topology_panel

result_dir <- "results/36-network-robustness"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  sparcc_roles,
  file.path(result_dir, "sparcc-node-roles.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  attack_summary,
  file.path(result_dir, "attack-summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  group_jaccard,
  file.path(result_dir, "group-edge-jaccard.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  group_bootstrap_frequency,
  file.path(result_dir, "group-edge-bootstrap.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

expected_bases <- c(
  "figures/36-role-cartography",
  "figures/36-attack-robustness",
  "figures/36-group-rewiring",
  "figures/36-topology-sensitivity"
)
expected_files <- as.vector(outer(
  expected_bases, c(".pdf", ".svg", ".png", ".tiff"), paste0
))
stopifnot(
  all(file.exists(expected_files)),
  sparcc_bootstrap_replicates == 100L,
  nrow(sparcc_pairs) == 780L,
  all(vapply(group_edge_sets, length, integer(1)) == group_edge_count),
  spiec_instability <= 0.05,
  rewiring_p >= 1 / 500,
  all(is.finite(attack_summary$RandomAUC))
)
