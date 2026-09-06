# βNTI + RCbray：确定性与随机过程的条件性证据
# Run sequentially in a new working directory.
# Required packages: ape, ggplot2, iCAMP, knitr, patchwork, ragg, scales, svglite.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/environment.tsv",
  "data/small/metadata.tsv",
  "data/small/otutab.tsv",
  "data/small/rooted-tree.nwk.gz",
  "data/small/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)
library(ape)
library(iCAMP)

set.seed(20260738)
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
full_tree <- ape::read.tree(gzfile("data/small/rooted-tree.nwk.gz"))

tree_tip_contract <- data.frame(
  TreeTips = ape::Ntip(full_tree),
  CountFeatures = nrow(otutab),
  CountFeaturesMissingFromTree = sum(!rownames(otutab) %in% full_tree$tip.label),
  TreeTipsOutsideCountTable = sum(!full_tree$tip.label %in% rownames(otutab)),
  Rooted = ape::is.rooted(full_tree),
  NonPositiveBranches = sum(
    !is.finite(full_tree$edge.length) | full_tree$edge.length <= 0
  )
)
stopifnot(
  nrow(otutab) == 13628L,
  ncol(otutab) == 90L,
  sum(otutab) == 1619670L,
  nrow(taxonomy) == nrow(otutab),
  tree_tip_contract$TreeTips == 14096L,
  tree_tip_contract$CountFeaturesMissingFromTree == 0L,
  tree_tip_contract$TreeTipsOutsideCountTable == 468L,
  tree_tip_contract$Rooted
)
tree_tip_contract

primary_keep <- rowMeans(otutab > 0) >= 0.10 & rowSums(otutab) >= 500
primary_counts <- otutab[primary_keep, , drop = FALSE]
community <- t(primary_counts)
pruned_tree <- ape::keep.tip(full_tree, colnames(community))
phylogenetic_distance <- ape::cophenetic.phylo(pruned_tree)
phylogenetic_distance <- phylogenetic_distance[
  colnames(community), colnames(community), drop = FALSE
]
sample_retention <- data.frame(
  SampleID = colnames(otutab),
  Group = metadata[colnames(otutab), "Group"],
  RetainedFraction = colSums(primary_counts) / colSums(otutab),
  stringsAsFactors = FALSE
)

stopifnot(
  ncol(community) == 585L,
  ape::Ntip(pruned_tree) == ncol(community),
  setequal(pruned_tree$tip.label, colnames(community)),
  identical(rownames(phylogenetic_distance), colnames(community)),
  all(is.finite(phylogenetic_distance)),
  all(phylogenetic_distance >= 0)
)

data.frame(
  RetainedOTUs = ncol(community),
  RetainedReadFraction = sum(primary_counts) / sum(otutab),
  MinimumSampleRetention = min(sample_retention$RetainedFraction),
  MaximumSampleRetention = max(sample_retention$RetainedFraction),
  PrunedTreeTips = ape::Ntip(pruned_tree)
)

relative_abundance <- sweep(otutab, 2L, colSums(otutab), "/")
retained_relative <- relative_abundance[primary_keep, , drop = FALSE]
taxon_site_weights <- sweep(
  retained_relative, 1L, rowSums(retained_relative), "/"
)
distance_index <- which(upper.tri(phylogenetic_distance), arr.ind = TRUE)
pair_distance <- phylogenetic_distance[distance_index]
short_range <- pair_distance <= stats::quantile(pair_distance, 0.10)
short_index <- distance_index[short_range, , drop = FALSE]
short_distance <- pair_distance[short_range]
signal_traits <- c("pH", "Conductivity", "TOC")

set.seed(20260738)
signal_results <- lapply(seq_along(signal_traits), function(trait_index) {
  trait <- signal_traits[[trait_index]]
  niche_optimum <- as.vector(
    taxon_site_weights %*% as.numeric(environment[[trait]])
  )
  niche_difference <- abs(
    niche_optimum[short_index[, "row"]] -
      niche_optimum[short_index[, "col"]]
  )
  observed <- stats::cor(
    short_distance, niche_difference, method = "spearman"
  )
  set.seed(20260738L + trait_index)
  null_rho <- replicate(499L, {
    permuted <- sample(niche_optimum)
    stats::cor(
      short_distance,
      abs(
        permuted[short_index[, "row"]] -
          permuted[short_index[, "col"]]
      ),
      method = "spearman"
    )
  })
  data.frame(
    Trait = trait,
    ObservedRho = observed,
    NullLower = stats::quantile(null_rho, 0.025),
    NullMedian = stats::median(null_rho),
    NullUpper = stats::quantile(null_rho, 0.975),
    PValue = (1 + sum(null_rho >= observed)) / 500,
    ShortRangePairs = length(short_distance),
    stringsAsFactors = FALSE
  )
})
signal_results <- do.call(rbind, signal_results)
signal_results$QValue <- stats::p.adjust(signal_results$PValue, method = "BH")
signal_results$Supported <- signal_results$ObservedRho > 0 &
  signal_results$QValue < 0.05
phylogenetic_signal_gate <- if (
  sum(signal_results$Supported) >= 2L
) "CONDITIONAL PASS" else "FAIL"
signal_results

seed_icamp_function <- function(fun, target, seed) {
  function_text <- paste(deparse(fun), collapse = "\n")
  replacement <- paste0(
    "parallel::clusterSetRNGStream(c1, iseed = ",
    as.integer(seed), "L)\n    ", target
  )
  stopifnot(grepl(target, function_text, fixed = TRUE))
  function_text <- sub(
    target, replacement, function_text, fixed = TRUE
  )
  seeded_function <- eval(
    parse(text = function_text), envir = environment(fun)
  )
  environment(seeded_function) <- environment(fun)
  seeded_function
}
seeded_bnti <- seed_icamp_function(
  iCAMP::bNTI.cm,
  "bMNTD.rand <- parallel::parLapply",
  20260738L
)
seeded_rcbray <- seed_icamp_function(
  iCAMP::RC.cm,
  "BC.rd <- parallel::parLapply",
  20261738L
)
n_workers <- max(
  1L,
  min(8L, parallel::detectCores(logical = FALSE))
)

# Small double-run probe: same seed and worker count must be byte-identical.
probe_community <- community[, head(colnames(community), 40L), drop = FALSE]
probe_distance <- phylogenetic_distance[
  colnames(probe_community), colnames(probe_community), drop = FALSE
]
probe_one <- seeded_bnti(
  comm = probe_community,
  dis = probe_distance,
  nworker = n_workers,
  weighted = TRUE,
  rand = 9L,
  output.bMNTD = FALSE,
  sig.index = "SES"
)$index
probe_two <- seeded_bnti(
  comm = probe_community,
  dis = probe_distance,
  nworker = n_workers,
  weighted = TRUE,
  rand = 9L,
  output.bMNTD = FALSE,
  sig.index = "SES"
)$index
stopifnot(identical(probe_one, probe_two))

primary_randomizations <- 999L
bnti_fit <- seeded_bnti(
  comm = community,
  dis = phylogenetic_distance,
  nworker = n_workers,
  weighted = TRUE,
  rand = primary_randomizations,
  output.bMNTD = TRUE,
  sig.index = "SES"
)
rcbray_fit <- seeded_rcbray(
  comm = community,
  nworker = n_workers,
  weighted = TRUE,
  rand = primary_randomizations,
  sig.index = "RC",
  output.bray = TRUE,
  silent = TRUE
)
bnti_matrix <- bnti_fit$index
rcbray_matrix <- rcbray_fit$index

pair_index <- which(lower.tri(bnti_matrix), arr.ind = TRUE)
assembly_pairs <- data.frame(
  Sample1 = rownames(bnti_matrix)[pair_index[, "row"]],
  Sample2 = colnames(bnti_matrix)[pair_index[, "col"]],
  BetaNTI = bnti_matrix[pair_index],
  RCbray = rcbray_matrix[pair_index],
  BrayCurtis = rcbray_fit$BC.obs[pair_index],
  stringsAsFactors = FALSE
)
assembly_pairs$Group1 <- metadata[assembly_pairs$Sample1, "Group"]
assembly_pairs$Group2 <- metadata[assembly_pairs$Sample2, "Group"]
assembly_pairs$PairScope <- ifelse(
  assembly_pairs$Group1 == assembly_pairs$Group2,
  paste("Within", assembly_pairs$Group1),
  "Between groups"
)
classify_process <- function(beta_nti, rcbray) {
  ifelse(
    beta_nti > 1.96, "Heterogeneous selection",
    ifelse(
      beta_nti < -1.96, "Homogeneous selection",
      ifelse(
        rcbray > 0.95, "Dispersal limitation",
        ifelse(
          rcbray < -0.95, "Homogenizing dispersal", "Undominated"
        )
      )
    )
  )
}
assembly_pairs$Process <- classify_process(
  assembly_pairs$BetaNTI, assembly_pairs$RCbray
)
process_levels <- c(
  "Heterogeneous selection", "Homogeneous selection",
  "Dispersal limitation", "Homogenizing dispersal", "Undominated"
)
assembly_pairs$Process <- factor(
  assembly_pairs$Process, levels = process_levels
)
assembly_pairs$PairScope <- factor(
  assembly_pairs$PairScope,
  levels = c("Within CW", "Within IW", "Within TW", "Between groups")
)
process_fraction <- as.data.frame(
  prop.table(table(assembly_pairs$PairScope, assembly_pairs$Process), 1L),
  stringsAsFactors = FALSE
)
names(process_fraction) <- c("PairScope", "Process", "Fraction")

stopifnot(
  nrow(assembly_pairs) == choose(90L, 2L),
  all(is.finite(assembly_pairs$BetaNTI)),
  all(assembly_pairs$RCbray >= -1 & assembly_pairs$RCbray <= 1),
  abs(sum(prop.table(table(assembly_pairs$Process))) - 1) < 1e-12
)
as.data.frame(prop.table(table(assembly_pairs$Process)))

strict_keep <- rowMeans(otutab > 0) >= 0.10 & rowSums(otutab) >= 1000
strict_community <- t(otutab[strict_keep, , drop = FALSE])
strict_tree <- ape::keep.tip(full_tree, colnames(strict_community))
strict_distance <- ape::cophenetic.phylo(strict_tree)
strict_distance <- strict_distance[
  colnames(strict_community), colnames(strict_community), drop = FALSE
]
seeded_bnti_strict <- seed_icamp_function(
  iCAMP::bNTI.cm,
  "bMNTD.rand <- parallel::parLapply",
  202607381L
)
seeded_rcbray_strict <- seed_icamp_function(
  iCAMP::RC.cm,
  "BC.rd <- parallel::parLapply",
  202617381L
)
strict_randomizations <- 199L
strict_bnti <- seeded_bnti_strict(
  comm = strict_community,
  dis = strict_distance,
  nworker = n_workers,
  weighted = TRUE,
  rand = strict_randomizations,
  output.bMNTD = FALSE,
  sig.index = "SES"
)$index
strict_rcbray <- seeded_rcbray_strict(
  comm = strict_community,
  nworker = n_workers,
  weighted = TRUE,
  rand = strict_randomizations,
  sig.index = "RC",
  silent = TRUE
)$index
strict_process <- classify_process(
  strict_bnti[pair_index], strict_rcbray[pair_index]
)
filter_process_agreement <- mean(
  strict_process == as.character(assembly_pairs$Process)
)
filter_sensitivity <- rbind(
  data.frame(
    Filter = "Total reads >= 500",
    OTUs = sum(primary_keep),
    ReadFraction = sum(otutab[primary_keep, ]) / sum(otutab),
    Randomizations = primary_randomizations,
    Process = levels(assembly_pairs$Process),
    Fraction = as.numeric(prop.table(table(assembly_pairs$Process))),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Filter = "Total reads >= 1000",
    OTUs = sum(strict_keep),
    ReadFraction = sum(otutab[strict_keep, ]) / sum(otutab),
    Randomizations = strict_randomizations,
    Process = process_levels,
    Fraction = as.numeric(prop.table(table(factor(
      strict_process, levels = process_levels
    )))),
    stringsAsFactors = FALSE
  )
)

data.frame(
  PrimaryOTUs = sum(primary_keep),
  StrictOTUs = sum(strict_keep),
  PrimaryReadFraction = sum(otutab[primary_keep, ]) / sum(otutab),
  StrictReadFraction = sum(otutab[strict_keep, ]) / sum(otutab),
  PairwiseProcessAgreement = filter_process_agreement
)

bnti_audit <- data.frame(
  Dimension = c(
    "Tree", "Phylogenetic signal", "Regional pool", "Feature filter",
    "Randomization", "Sequential gate", "Pair dependence", "Claim"
  ),
  FragilePractice = c(
    "Taxonomy tree or unmatched tips",
    "Assume close relatives share niches",
    "Mix unrelated habitats without definition",
    "Hide retained read fraction",
    "Parallel workers without seeded streams",
    "Use RCbray even when |bNTI| > 1.96",
    "Treat all sample pairs as independent",
    "Call fractions measured mechanisms"
  ),
  UpgradedContract = c(
    "Rooted branch-length tree plus exact tip audit",
    "Short-range niche-signal gate",
    "Declare metacommunity and pair scope",
    "Primary plus strict filter ledger",
    "999 nulls with deterministic cluster RNG",
    "Selection gate before RCbray",
    "Descriptive fractions or sample-level resampling",
    "Conditional process evidence"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(bnti_audit, caption = "βNTI/RCbray interpretation audit")

tree_counts <- data.frame(
  Stage = factor(
    c("Source tree", "Count table", "Primary universe"),
    levels = c("Source tree", "Count table", "Primary universe")
  ),
  Tips = c(ape::Ntip(full_tree), nrow(otutab), sum(primary_keep))
)
tree_panel <- ggplot(tree_counts, aes(x = Stage, y = Tips, fill = Stage)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(
    aes(label = scales::comma(Tips)),
    vjust = -0.35, family = font_pub, size = 2.8
  ) +
  scale_fill_manual(values = unname(pal_pub[c("grey", "blue", "orange")])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "A. Tree tips are explicitly reconciled",
    x = NULL, y = "Tips / features"
  ) +
  theme_pub(base_size = 8.0) +
  theme(axis.text.x = element_text(angle = 24, hjust = 1))

retention_panel <- ggplot(
  sample_retention,
  aes(x = Group, y = RetainedFraction, fill = Group)
) +
  geom_violin(trim = FALSE, alpha = 0.40, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, colour = "#333333") +
  geom_jitter(
    width = 0.10, shape = 21, size = 1.5,
    colour = "white", stroke = 0.25, alpha = 0.75
  ) +
  scale_fill_manual(values = c(
    CW = pal_pub[["blue"]], IW = pal_pub[["orange"]], TW = pal_pub[["green"]]
  )) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    title = "B. Dominant-OTU read retention",
    subtitle = sprintf("Pooled coverage = %.1f%%", 100 * sum(primary_counts) / sum(otutab)),
    x = "Wetland group", y = "Retained sample reads"
  ) +
  theme_pub(base_size = 8.0) +
  theme(legend.position = "none")

tree_filter_plot <- patchwork::wrap_plots(
  tree_panel, retention_panel, nrow = 1L
) +
  patchwork::plot_annotation(
    title = "Phylogenetic inference begins with a feature-universe audit",
    caption = "The primary result describes 585 dominant OTUs, not the full rare biosphere."
  )
save_pub(tree_filter_plot, "figures/38-tree-filter", width = 180, height = 122)
tree_filter_plot

process_colours <- c(
  "Heterogeneous selection" = pal_pub[["vermillion"]],
  "Homogeneous selection" = pal_pub[["blue"]],
  "Dispersal limitation" = pal_pub[["orange"]],
  "Homogenizing dispersal" = pal_pub[["green"]],
  "Undominated" = pal_pub[["grey"]]
)
bnti_rc_plot <- ggplot(
  assembly_pairs,
  aes(x = BetaNTI, y = RCbray, colour = Process)
) +
  geom_hline(yintercept = c(-0.95, 0.95), linetype = 2, colour = "#777777") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 2, colour = "#777777") +
  geom_point(size = 1.25, alpha = 0.42) +
  scale_colour_manual(values = process_colours, drop = FALSE) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Selection is evaluated before RCbray",
    subtitle = sprintf(
      "%s; %d pairwise comparisons; %d null randomizations",
      phylogenetic_signal_gate, nrow(assembly_pairs), primary_randomizations
    ),
    x = expression(beta * "NTI"), y = expression("RC"[Bray]),
    colour = "Process",
    caption = "RCbray labels apply only to pairs with |βNTI| <= 1.96; points are not independent observations."
  ) +
  theme_pub(base_size = 8.5) +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"))
save_pub(bnti_rc_plot, "figures/38-bnti-rcbray", width = 180, height = 148)
bnti_rc_plot

process_fraction$Process <- factor(
  process_fraction$Process, levels = rev(process_levels)
)
process_fraction_plot <- ggplot(
  process_fraction,
  aes(x = PairScope, y = Fraction, fill = Process)
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = process_colours, drop = FALSE) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  labs(
    title = "Dispersal limitation dominates this filtered pairwise ledger",
    subtitle = "Fractions summarize dependent sample pairs and are not tested with an ordinary contingency table",
    x = NULL, y = "Pairwise process fraction", fill = "Process",
    caption = "Between-group and within-group results share the same dominant-OTU regional pool."
  ) +
  theme_pub(base_size = 8.3) +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1),
    legend.position = "bottom", legend.box = "vertical"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"))
save_pub(
  process_fraction_plot, "figures/38-process-fractions",
  width = 180, height = 132
)
process_fraction_plot

signal_results$Trait <- factor(
  signal_results$Trait, levels = rev(signal_traits)
)
signal_panel <- ggplot(
  signal_results,
  aes(x = ObservedRho, y = Trait, colour = Supported)
) +
  geom_errorbarh(
    aes(xmin = NullLower, xmax = NullUpper),
    height = 0.18, colour = pal_pub[["grey"]]
  ) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_point(size = 3) +
  geom_text(
    aes(
      label = sprintf("q=%.3f", QValue),
      hjust = ifelse(ObservedRho > 0.12, 1.1, -0.1)
    ),
    nudge_y = 0.23, family = font_pub, size = 2.5,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "FALSE" = pal_pub[["grey"]], "TRUE" = pal_pub[["vermillion"]]
  )) +
  labs(
    title = "A. Short-range niche signal gate",
    subtitle = "Observed rho with 95% label-permutation interval",
    x = "Phylogenetic distance vs niche difference rho", y = NULL,
    colour = "BH supported"
  ) +
  theme_pub(base_size = 8.0) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(-0.06, 0.185), clip = "off")

filter_sensitivity$Filter <- factor(
  filter_sensitivity$Filter,
  levels = c("Total reads >= 500", "Total reads >= 1000")
)
filter_panel <- ggplot(
  filter_sensitivity,
  aes(x = Filter, y = Fraction, fill = Process)
) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = process_colours, drop = FALSE) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  labs(
    title = "B. Feature-filter sensitivity",
    subtitle = sprintf("Pairwise label agreement = %.1f%%", 100 * filter_process_agreement),
    x = NULL, y = "Overall process fraction", fill = "Process"
  ) +
  theme_pub(base_size = 8.0) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "none"
  )

signal_filter_plot <- patchwork::wrap_plots(
  signal_panel, filter_panel, nrow = 1L, widths = c(1.05, 0.95)
) +
  patchwork::plot_annotation(
    title = "Assembly-process labels need assumption and preprocessing gates",
    caption = "The primary analysis covers 47.4% of reads; rare-biosphere processes are outside its estimand."
  )
save_pub(
  signal_filter_plot, "figures/38-phylogenetic-signal",
  width = 180, height = 127
)
signal_filter_plot

result_dir <- "results/38-community-assembly-bnti"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  tree_tip_contract,
  file.path(result_dir, "tree-tip-contract.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  signal_results,
  file.path(result_dir, "phylogenetic-signal.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  assembly_pairs,
  file.path(result_dir, "pairwise-process-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  filter_sensitivity,
  file.path(result_dir, "filter-sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

expected_bases <- c(
  "figures/38-tree-filter",
  "figures/38-bnti-rcbray",
  "figures/38-process-fractions",
  "figures/38-phylogenetic-signal"
)
expected_files <- as.vector(outer(
  expected_bases, c(".pdf", ".svg", ".png", ".tiff"), paste0
))
stopifnot(
  all(file.exists(expected_files)),
  as.character(utils::packageVersion("iCAMP")) == "1.5.12",
  primary_randomizations == 999L,
  strict_randomizations == 199L,
  nrow(assembly_pairs) == 4005L,
  sum(primary_keep) == 585L,
  sum(strict_keep) == 220L,
  identical(probe_one, probe_two),
  phylogenetic_signal_gate == "CONDITIONAL PASS"
)
