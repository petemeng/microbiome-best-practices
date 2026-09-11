# CoDA 成分数据原理：为什么不能比原始丰度
# Run sequentially in a new working directory.
# Required packages: ggplot2, ggrepel, knitr, ragg, scales, svglite.

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

set.seed(20260730)
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

stopifnot(
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata)),
  nrow(otutab) == 13628L,
  ncol(otutab) == 90L,
  sum(otutab) == 1619670L,
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
  lapply(taxonomy[rownames(otutab), rank_names], clean_taxon),
  check.names = FALSE
)

known_genus <- !is.na(taxonomy_clean$Genus)
lineage_key <- apply(
  taxonomy_clean[known_genus, rank_names[1:6], drop = FALSE],
  1L,
  function(x) paste(ifelse(is.na(x), "?", x), collapse = ";")
)
genus_counts <- rowsum(
  otutab[known_genus, , drop = FALSE], lineage_key, reorder = FALSE
)
genus_lineages <- rownames(genus_counts)
genus_labels <- vapply(
  strsplit(genus_lineages, ";", fixed = TRUE), tail, character(1), 1L
)
genus_ids <- make.unique(paste0("Genus:", genus_labels))
rownames(genus_counts) <- genus_ids

genus_info <- data.frame(
  FeatureID = genus_ids,
  DisplayTaxon = genus_labels,
  Lineage = genus_lineages,
  Reportable = TRUE,
  stringsAsFactors = FALSE,
  row.names = genus_ids
)
genus_counts <- rbind(
  genus_counts,
  Unclassified_residual = colSums(otutab[!known_genus, , drop = FALSE])
)
genus_info <- rbind(
  genus_info,
  data.frame(
    FeatureID = "Unclassified_residual",
    DisplayTaxon = "Unclassified residual",
    Lineage = "Unclassified at genus rank",
    Reportable = FALSE,
    row.names = "Unclassified_residual"
  )
)

data.frame(
  Samples = ncol(genus_counts),
  NamedGenera = sum(genus_info$Reportable),
  NamedReadFraction = sum(genus_counts[genus_info$Reportable, ]) / sum(otutab),
  LibraryMin = min(colSums(otutab)),
  LibraryMedian = median(colSums(otutab)),
  LibraryMax = max(colSums(otutab))
)

named_ids <- rownames(genus_info)[genus_info$Reportable]
target_id <- named_ids[which.max(rowSums(genus_counts[named_ids, , drop = FALSE]))]
target_label <- genus_info[target_id, "DisplayTaxon"]

# 90 个真实样本的平均 count composition；只把目标属乘以 4。
absolute_baseline <- rowMeans(genus_counts)
absolute_after <- absolute_baseline
absolute_after[target_id] <- 4 * absolute_after[target_id]

relative_baseline <- absolute_baseline / sum(absolute_baseline)
relative_after <- absolute_after / sum(absolute_after)
closure_audit <- data.frame(
  FeatureID = names(absolute_baseline),
  DisplayTaxon = genus_info[names(absolute_baseline), "DisplayTaxon"],
  AbsoluteFoldChange = absolute_after / absolute_baseline,
  RelativeFoldChange = relative_after / relative_baseline,
  Status = ifelse(names(absolute_baseline) == target_id, "Changed taxon", "Unchanged taxon"),
  stringsAsFactors = FALSE
)
closure_audit <- closure_audit[is.finite(closure_audit$RelativeFoldChange), ]

unchanged_relative_fold <- unique(round(
  closure_audit$RelativeFoldChange[closure_audit$Status == "Unchanged taxon"],
  12
))
stopifnot(
  length(unchanged_relative_fold) == 1L,
  unchanged_relative_fold < 1,
  all(closure_audit$AbsoluteFoldChange[
    closure_audit$Status == "Unchanged taxon"
  ] == 1)
)

closure_audit[order(closure_audit$RelativeFoldChange), ][1:8, ]

library_size <- colSums(otutab)
relative_all <- sweep(genus_counts, 2L, library_size, "/")
pseudocount <- 0.5
log_counts <- log2(genus_counts + pseudocount)
clr_all <- sweep(log_counts, 2L, colMeans(log_counts), "-")

measurement_views <- rbind(
  data.frame(
    SampleID = colnames(genus_counts), Group = metadata$Group,
    Scale = "Raw reads", Value = as.numeric(genus_counts[target_id, ])
  ),
  data.frame(
    SampleID = colnames(genus_counts), Group = metadata$Group,
    Scale = "Relative abundance (%)",
    Value = 100 * as.numeric(relative_all[target_id, ])
  ),
  data.frame(
    SampleID = colnames(genus_counts), Group = metadata$Group,
    Scale = "CLR abundance", Value = as.numeric(clr_all[target_id, ])
  )
)
measurement_views$Scale <- factor(
  measurement_views$Scale,
  levels = c("Raw reads", "Relative abundance (%)", "CLR abundance")
)

# Scale invariance: doubling every count in a sample leaves its composition unchanged.
sample_one <- genus_counts[, 1L]
stopifnot(all.equal(
  sample_one / sum(sample_one),
  (2 * sample_one) / sum(2 * sample_one),
  tolerance = 1e-14
))

reference_candidates <- setdiff(
  named_ids[order(rowSums(genus_counts[named_ids, , drop = FALSE]), decreasing = TRUE)],
  target_id
)[1:3]

reference_effects <- do.call(rbind, lapply(reference_candidates, function(ref_id) {
  alr <- log2(
    (genus_counts[target_id, ] + pseudocount) /
      (genus_counts[ref_id, ] + pseudocount)
  )
  do.call(rbind, lapply(c("CW", "TW"), function(group_name) {
    data.frame(
      Contrast = paste0(group_name, " - IW"),
      Reference = genus_info[ref_id, "DisplayTaxon"],
      MedianLog2RatioDifference = median(alr[metadata$Group == group_name]) -
        median(alr[metadata$Group == "IW"]),
      stringsAsFactors = FALSE
    )
  }))
}))
reference_effects

coda_audit <- data.frame(
  InputOrClaim = c(
    "Raw reads", "Relative abundance", "CLR/ALR", "Zeros",
    "Taxonomy collapse", "Absolute claim"
  ),
  NaiveInterpretation = c(
    "Read count equals microbial load",
    "Percentage change equals absolute change",
    "Transformation removes all bias",
    "Add one without reporting it",
    "Merge every unknown label",
    "Infer cells per gram from 16S alone"
  ),
  UpgradedContract = c(
    "Model library size and method-specific sampling",
    "Report closure and denominator dependence",
    "State the reference frame and feature set",
    "Prespecify and sensitivity-test replacement",
    "Aggregate by full lineage; residual is not a taxon",
    "Collect QMP, qPCR/ddPCR, flow cytometry or spike-in data"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(coda_audit, caption = "Compositional audit ledger")

label_ids <- unique(c(
  target_id,
  head(
    names(sort(relative_baseline[named_ids], decreasing = TRUE)),
    8L
  )
))
closure_plot_data <- data.frame(
  FeatureID = names(relative_baseline),
  Before = as.numeric(relative_baseline),
  After = as.numeric(relative_after),
  Status = ifelse(names(relative_baseline) == target_id, "Changed taxon", "Unchanged taxon"),
  Label = ifelse(
    names(relative_baseline) %in% label_ids,
    genus_info[names(relative_baseline), "DisplayTaxon"],
    NA_character_
  )
)

closure_plot <- ggplot(
  closure_plot_data,
  aes(x = Before, y = After, colour = Status)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#777777") +
  geom_point(alpha = 0.72, size = 1.6) +
  ggrepel::geom_text_repel(
    aes(label = Label), family = font_pub, size = 2.5,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_x_log10(labels = scales::label_percent(accuracy = 0.01)) +
  scale_y_log10(labels = scales::label_percent(accuracy = 0.01)) +
  scale_colour_manual(values = c(
    "Changed taxon" = pal_pub[["vermillion"]],
    "Unchanged taxon" = pal_pub[["blue"]]
  )) +
  coord_equal() +
  labs(
    title = "Closure creates apparent changes in unchanged taxa",
    subtitle = paste0(
      target_label, " absolute abundance was multiplied by four; all other taxa were fixed"
    ),
    x = "Relative abundance before perturbation",
    y = "Relative abundance after perturbation",
    colour = "Counterfactual status",
    caption = "A point below the dashed identity line lost relative share, not absolute abundance."
  ) +
  theme_pub(base_size = 8.8) +
  theme(legend.position = "bottom")

save_pub(closure_plot, "figures/30-closure-artifact", width = 170, height = 122)
closure_plot

three_scales_plot <- ggplot(
  measurement_views,
  aes(x = Group, y = Value, fill = Group)
) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.68) +
  geom_point(
    position = position_jitter(width = 0.09, height = 0, seed = 20260730),
    shape = 21, colour = "white", stroke = 0.25, size = 1.45, alpha = 0.75
  ) +
  facet_wrap(~Scale, scales = "free_y", nrow = 1L) +
  scale_fill_manual(values = c(
    CW = pal_pub[["blue"]], IW = pal_pub[["orange"]], TW = pal_pub[["green"]]
  )) +
  labs(
    title = paste0(target_label, " under three measurement scales"),
    subtitle = "Raw reads retain library size; relative and CLR values use different references",
    x = "Wetland group", y = "Observed value", fill = "Group",
    caption = "CLR used a 0.5-count replacement across the complete genus table including the residual."
  ) +
  theme_pub(base_size = 8.5) +
  theme(legend.position = "none")

save_pub(three_scales_plot, "figures/30-measurement-scales", width = 180, height = 116)
three_scales_plot

reference_plot <- ggplot(
  reference_effects,
  aes(x = MedianLog2RatioDifference, y = Reference, colour = Contrast)
) +
  geom_vline(xintercept = 0, linetype = "22", colour = "#777777") +
  geom_segment(
    aes(x = 0, xend = MedianLog2RatioDifference, yend = Reference),
    linewidth = 0.65, alpha = 0.55
  ) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = c(
    "CW - IW" = pal_pub[["blue"]], "TW - IW" = pal_pub[["vermillion"]]
  )) +
  labs(
    title = "The reference frame defines the reported effect",
    subtitle = paste0("Target: ", target_label, "; positive values favour the first group"),
    x = "Difference in median log2(target / reference)",
    y = "Reference genus", colour = "Contrast",
    caption = "References were the three most abundant named genera excluding the target."
  ) +
  theme_pub(base_size = 8.8) +
  theme(legend.position = "bottom")

save_pub(reference_plot, "figures/30-reference-frame", width = 170, height = 107)
reference_plot

stopifnot(
  target_id %in% rownames(genus_counts),
  length(reference_candidates) == 3L,
  nrow(reference_effects) == 6L,
  all(is.finite(reference_effects$MedianLog2RatioDifference)),
  all(file.exists(paste0(
    rep(c(
      "figures/30-closure-artifact", "figures/30-measurement-scales",
      "figures/30-reference-frame"
    ), each = 4L),
    rep(c(".pdf", ".svg", ".png", ".tiff"), times = 3L)
  )))
)
