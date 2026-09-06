# FAPROTAX、Tax4Fun2 与 BugBase：环境功能类群
# Run sequentially in a new working directory.
# Required packages: ggplot2, ggrepel, microeco, ragg, readr, scales, svglite.

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
library(microeco)

set.seed(20260742)
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
      panel.grid.major = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.25),
      axis.text = ggplot2::element_text(colour = "#1A1A1A"),
      axis.title = ggplot2::element_text(colour = "#1A1A1A"),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.15)),
      plot.subtitle = ggplot2::element_text(colour = "#4D4D4D"),
      plot.caption = ggplot2::element_text(colour = "#666666", hjust = 0),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = "#B3B3B3"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(), legend.position = "top"
    )
}
save_pub <- function(plot, file_base, width = 89, height = 70, units = "mm", dpi = 600) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
                  units = units, device = grDevices::cairo_pdf, family = font_pub, bg = "white")
  ggplot2::ggsave(paste0(file_base, ".svg"), plot, width = width, height = height,
                  units = units, device = svglite::svglite, bg = "white")
  ggplot2::ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
                  units = units, dpi = dpi, device = ragg::agg_png, bg = "white")
  ggplot2::ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
                  units = units, dpi = dpi, device = ragg::agg_tiff,
                  compression = "lzw", bg = "white")
  invisible(plot)
}

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE,
                       name_repair = "minimal", na = character())
  ids <- as.character(x[[1L]])
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}
otutab_raw <- read_keyed_tsv("data/small/otutab.tsv")
taxonomy <- read_keyed_tsv("data/small/taxonomy.tsv")
metadata <- read_keyed_tsv("data/small/metadata.tsv")
otutab <- as.matrix(data.frame(lapply(otutab_raw, as.integer), check.names = FALSE))
rownames(otutab) <- rownames(otutab_raw)
metadata$Group <- factor(metadata$Group, levels = c("CW", "IW", "TW"))
otutab <- otutab[, rownames(metadata), drop = FALSE]

stopifnot(
  nrow(otutab) == 13628L, ncol(otutab) == 90L,
  sum(otutab) == 1619670L,
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata))
)

dataset <- microeco::microtable$new(
  otu_table = as.data.frame(otutab),
  sample_table = metadata,
  tax_table = taxonomy,
  auto_tidy = TRUE
)
head(otutab[, 1:5])
head(taxonomy)
head(metadata)

faprotax <- microeco::trans_func$new(dataset)
faprotax$cal_func(prok_database = "FAPROTAX")
function_map <- as.matrix(faprotax$res_func)
function_map <- function_map[rownames(otutab), , drop = FALSE]
mapped_taxon <- rowSums(function_map) > 0
function_counts <- t(function_map) %*% otutab
coverage <- colSums(otutab[mapped_taxon, , drop = FALSE]) / colSums(otutab)

coverage_ledger <- data.frame(
  SampleID = names(coverage), Group = metadata[names(coverage), "Group"],
  AnnotationCoverage = unname(coverage), stringsAsFactors = FALSE
)
mapping_summary <- data.frame(
  Features = nrow(function_map), Functions = ncol(function_map),
  MappedFeatures = sum(mapped_taxon), MappedReads = sum(otutab[mapped_taxon, ]),
  TotalReads = sum(otutab), OverallReadCoverage = sum(otutab[mapped_taxon, ]) / sum(otutab)
)

stopifnot(
  all(function_map %in% c(0, 1)), ncol(function_map) == 93L,
  sum(mapped_taxon) == 3992L,
  sum(otutab[mapped_taxon, ]) == 562755L,
  all(coverage >= 0 & coverage <= 1)
)
mapping_summary

assignment_share <- sweep(function_counts, 2L, colSums(function_counts), "/")
clr_counts <- log(function_counts + 0.5)
clr_counts <- sweep(clr_counts, 2L, colMeans(clr_counts), "-")

test_one_function <- function(feature) {
  values <- clr_counts[feature, ]
  test <- stats::kruskal.test(values ~ metadata$Group)
  epsilon2 <- max(0, (unname(test$statistic) - nlevels(metadata$Group) + 1) /
                       (length(values) - nlevels(metadata$Group)))
  medians <- tapply(assignment_share[feature, ], metadata$Group, median)
  data.frame(
    Function = feature, PValue = test$p.value, EpsilonSquared = epsilon2,
    CWMedian = medians[["CW"]], IWMedian = medians[["IW"]],
    TWMedian = medians[["TW"]], stringsAsFactors = FALSE
  )
}
function_tests <- do.call(rbind, lapply(rownames(function_counts), test_one_function))
function_tests$AdjustedP <- stats::p.adjust(function_tests$PValue, method = "BH")
function_tests <- function_tests[order(-function_tests$EpsilonSquared, function_tests$AdjustedP), ]
function_tests$FunctionLabel <- gsub("_", " ", function_tests$Function, fixed = TRUE)

stopifnot(nrow(function_tests) == 93L, all(function_tests$AdjustedP >= function_tests$PValue))
head(function_tests, 12)

taxonomy_genus <- taxonomy
taxonomy_genus[, setdiff(colnames(taxonomy_genus), "Genus")] <- ""
dataset_genus <- microeco::microtable$new(
  otu_table = as.data.frame(otutab), sample_table = metadata,
  tax_table = taxonomy_genus, auto_tidy = TRUE
)
faprotax_genus <- microeco::trans_func$new(dataset_genus)
faprotax_genus$for_what <- "prok"
faprotax_genus$cal_func(prok_database = "FAPROTAX")
map_genus <- as.matrix(faprotax_genus$res_func)[rownames(otutab), , drop = FALSE]
counts_genus <- t(map_genus) %*% otutab
share_genus <- sweep(counts_genus, 2L, pmax(colSums(counts_genus), 1), "/")

function_sensitivity <- data.frame(
  Function = rownames(function_counts),
  FullLineageMean = rowMeans(assignment_share),
  GenusOnlyMean = rowMeans(share_genus),
  SampleSpearman = vapply(rownames(function_counts), function(feature) {
    if (stats::sd(assignment_share[feature, ]) == 0 || stats::sd(share_genus[feature, ]) == 0) return(NA_real_)
    stats::cor(assignment_share[feature, ], share_genus[feature, ], method = "spearman")
  }, numeric(1)),
  stringsAsFactors = FALSE
)
function_sensitivity$MeanAbsoluteDifference <- abs(
  function_sensitivity$FullLineageMean - function_sensitivity$GenusOnlyMean
)
function_sensitivity$FunctionLabel <- gsub("_", " ", function_sensitivity$Function, fixed = TRUE)

resolution_coverage <- data.frame(
  Mapping = c("Full lineage", "Genus only"),
  ReadCoverage = c(sum(otutab[mapped_taxon, ]) / sum(otutab),
                   sum(otutab[rowSums(map_genus) > 0, ]) / sum(otutab))
)
resolution_coverage

p_coverage <- ggplot(coverage_ledger, aes(Group, AnnotationCoverage, colour = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "#4D4D4D") +
  geom_jitter(width = 0.12, height = 0, alpha = 0.65, size = 1.5) +
  scale_color_pub(values = c(CW = pal_pub[["blue"]], IW = pal_pub[["orange"]], TW = pal_pub[["green"]])) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(title = "FAPROTAX annotation coverage", x = NULL,
       y = "Annotated read fraction",
       caption = "Coverage precedes multi-function duplication.") +
  theme_pub() + theme(legend.position = "none")
save_pub(p_coverage, "figures/42-function-coverage", 110, 82)
p_coverage

top_functions <- names(sort(rowMeans(assignment_share), decreasing = TRUE))[1:10]
composition_long <- do.call(rbind, lapply(top_functions, function(feature) {
  aggregate(assignment_share[feature, ], list(Group = metadata$Group), mean) |>
    transform(Function = gsub("_", " ", feature, fixed = TRUE)) |>
    setNames(c("Group", "AssignmentShare", "Function"))
}))
composition_long$Function <- factor(
  composition_long$Function,
  levels = rev(gsub("_", " ", top_functions, fixed = TRUE))
)
p_composition <- ggplot(composition_long, aes(AssignmentShare, Function, fill = Group)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  scale_fill_pub(values = c(CW = pal_pub[["blue"]], IW = pal_pub[["orange"]], TW = pal_pub[["green"]])) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Most abundant predicted functions", x = "Mean share of function assignments", y = NULL,
       caption = "Categories are non-exclusive assignment shares.") +
  theme_pub()
save_pub(p_composition, "figures/42-functional-composition", 178, 105)
p_composition

heat_functions <- function_tests$Function[1:12]
heat <- do.call(rbind, lapply(heat_functions, function(feature) {
  means <- tapply(clr_counts[feature, ], metadata$Group, mean)
  z <- as.numeric(scale(means))
  data.frame(Function = gsub("_", " ", feature, fixed = TRUE), Group = names(means), Z = z)
}))
heat$Function <- factor(heat$Function, levels = rev(gsub("_", " ", heat_functions, fixed = TRUE)))
p_heat <- ggplot(heat, aes(Group, Function, fill = Z)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = pal_pub[["blue"]], mid = "white", high = pal_pub[["vermillion"]], midpoint = 0) +
  labs(title = "Predicted-function effect audit", x = NULL, y = NULL, fill = "Row z-score",
       caption = "BH correction covers all 93 functions.") +
  theme_pub()
save_pub(p_heat, "figures/42-functional-heatmap", 135, 105)
p_heat

label_rows <- head(function_sensitivity[order(-function_sensitivity$MeanAbsoluteDifference), ], 8)
p_sensitivity <- ggplot(function_sensitivity, aes(FullLineageMean, GenusOnlyMean)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#6B7280") +
  geom_point(aes(colour = SampleSpearman), size = 2, alpha = 0.8) +
  ggrepel::geom_text_repel(data = label_rows, aes(label = FunctionLabel),
                           size = 2.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_colour_gradient2(low = pal_pub[["vermillion"]], mid = "#F2F2F2",
                         high = pal_pub[["blue"]], midpoint = 0.5, na.value = pal_pub[["grey"]]) +
  scale_x_log10(labels = scales::label_percent(accuracy = 0.01)) +
  scale_y_log10(labels = scales::label_percent(accuracy = 0.01)) +
  labs(title = "Taxonomic-resolution sensitivity", x = "Mean assignment share: full lineage",
       y = "Mean assignment share: genus only", colour = "Sample rho",
       caption = "Large departures from the diagonal identify taxonomy-dependent predictions.") +
  theme_pub()
save_pub(p_sensitivity, "figures/42-functional-sensitivity", 130, 100)
p_sensitivity

dir.create("results/42-functional-guilds", recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(mapping_summary, "results/42-functional-guilds/mapping-summary.tsv")
readr::write_tsv(coverage_ledger, "results/42-functional-guilds/sample-coverage.tsv")
readr::write_tsv(function_tests, "results/42-functional-guilds/function-tests.tsv")
readr::write_tsv(function_sensitivity, "results/42-functional-guilds/taxonomy-resolution-sensitivity.tsv")
readr::write_tsv(resolution_coverage, "results/42-functional-guilds/resolution-coverage.tsv")
