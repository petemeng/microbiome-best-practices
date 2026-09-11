# 多界微生物组：16S＋ITS/18S 联合分析
# Run sequentially in a new working directory.
# Required packages: ggplot2, jsonlite, permute, ragg, readr, scales, svglite, vegan.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "scripts/run_article49_multi_kingdom.R",
  "data/small/multi-kingdom-duran/fungi-otutab.tsv",
  "data/small/multi-kingdom-duran/fungi-taxonomy.tsv",
  "data/small/multi-kingdom-duran/metadata.tsv",
  "data/small/multi-kingdom-duran/oomycete-otutab.tsv",
  "data/small/multi-kingdom-duran/oomycete-taxonomy.tsv",
  "data/small/multi-kingdom-duran/otutab.tsv",
  "data/small/multi-kingdom-duran/source-summary.json",
  "data/small/multi-kingdom-duran/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)

set.seed(20260749)
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
data_dir <- "data/small/multi-kingdom-duran"
metadata <- read_keyed_tsv(file.path(data_dir, "metadata.tsv"))
view_files <- c(
  Bacteria = "otutab.tsv", Fungi = "fungi-otutab.tsv",
  Oomycete = "oomycete-otutab.tsv"
)
taxonomy_files <- c(
  Bacteria = "taxonomy.tsv", Fungi = "fungi-taxonomy.tsv",
  Oomycete = "oomycete-taxonomy.tsv"
)
counts <- lapply(view_files, function(filename) {
  raw <- read_keyed_tsv(file.path(data_dir, filename))
  value <- as.matrix(data.frame(lapply(raw, as.numeric), check.names = FALSE))
  rownames(value) <- rownames(raw)
  value
})
taxonomy <- lapply(taxonomy_files, function(filename) {
  read_keyed_tsv(file.path(data_dir, filename))
})
stopifnot(
  identical(vapply(counts, nrow, integer(1)), c(Bacteria = 772L, Fungi = 1063L, Oomycete = 219L)),
  all(vapply(counts, ncol, integer(1)) == 36L),
  all(vapply(counts, function(x) identical(colnames(x), rownames(metadata)), logical(1))),
  identical(as.integer(table(metadata$Soil, metadata$Compartment)), rep(4L, 9))
)
head(counts$Bacteria[, 1:5])
head(taxonomy$Bacteria)
head(metadata)

run_analysis <- function(command, args) {
  output <- system2(command, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop(paste(tail(output, 20), collapse = "\n"))
  invisible(output)
}
run_analysis(file.path(R.home("bin"), "Rscript"), c("--vanilla", "scripts/run_article49_multi_kingdom.R"))

result_dir <- "results/49-multi-kingdom"
ordination_scores <- readr::read_tsv(
  file.path(result_dir, "ordination-scores.tsv"), show_col_types = FALSE
)
pcoa_variance <- readr::read_tsv(
  file.path(result_dir, "pcoa-variance.tsv"), show_col_types = FALSE
)
procrustes_summary <- readr::read_tsv(
  file.path(result_dir, "procrustes-summary.tsv"), show_col_types = FALSE
)
associations <- readr::read_tsv(
  file.path(result_dir, "cross-kingdom-associations.tsv"), show_col_types = FALSE
)
multi_summary <- jsonlite::read_json(
  file.path(result_dir, "summary.json"), simplifyVector = TRUE
)
stopifnot(
  nrow(ordination_scores) == 108L,
  nrow(procrustes_summary) == 6L,
  nrow(associations) == 1200L,
  multi_summary$raw_cross_kingdom_q_lt_0_05 == 505L,
  multi_summary$adjusted_cross_kingdom_q_lt_0_05 == 2L,
  multi_summary$stable_candidates_q_lt_0_05 == 2L
)
head(associations, 6)

axis_1 <- setNames(
  pcoa_variance$Percent[pcoa_variance$Axis == "Axis1"],
  pcoa_variance$Kingdom[pcoa_variance$Axis == "Axis1"]
)
axis_2 <- setNames(
  pcoa_variance$Percent[pcoa_variance$Axis == "Axis2"],
  pcoa_variance$Kingdom[pcoa_variance$Axis == "Axis2"]
)
facet_labels <- setNames(
  sprintf("%s\nAxis 1 %.1f%% | Axis 2 %.1f%%", names(axis_1), axis_1, axis_2[names(axis_1)]),
  names(axis_1)
)
compartment_colors <- c(
  Soil = unname(pal_pub["grey"]),
  Rhizosphere = unname(pal_pub["orange"]),
  Root = unname(pal_pub["green"])
)
p49_1 <- ggplot(ordination_scores, aes(Axis1, Axis2, colour = Compartment, shape = Soil)) +
  geom_hline(yintercept = 0, colour = "#E5E5E5", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#E5E5E5", linewidth = 0.3) +
  geom_point(size = 2.1, alpha = 0.82) +
  stat_ellipse(aes(group = Compartment), level = 0.80, linewidth = 0.5,
               show.legend = FALSE) +
  facet_wrap(~ Kingdom, scales = "free", labeller = as_labeller(facet_labels)) +
  scale_color_manual(values = compartment_colors) +
  labs(x = "PCoA axis 1", y = "PCoA axis 2", colour = "Compartment", shape = "Soil",
       title = "Each kingdom is ordinated independently") + theme_pub()
save_pub(p49_1, "figures/49-multi-kingdom/49-1-three-kingdom-pcoa", 180, 82)
p49_1

procrustes_summary$Pair <- factor(
  procrustes_summary$Pair,
  levels = rev(c("Bacteria vs Fungi", "Bacteria vs Oomycete", "Fungi vs Oomycete"))
)
p49_2 <- ggplot(procrustes_summary, aes(Correlation, Pair, colour = Model)) +
  geom_segment(aes(x = 0, xend = Correlation, yend = Pair), linewidth = 0.6, alpha = 0.45) +
  geom_point(aes(shape = QValue < 0.05), size = 3) +
  geom_text(aes(label = sprintf("q=%.3f", QValue)), nudge_x = 0.045,
            family = font_pub, size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = c(
    Observed = unname(pal_pub["grey"]),
    `Niche-adjusted` = unname(pal_pub["blue"])
  )) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "Not significant")) +
  scale_x_continuous(limits = c(0, 0.95), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Procrustes correlation", y = NULL, colour = "Model", shape = "Inference",
       title = "Cross-kingdom community concordance") + theme_pub()
save_pub(p49_2, "figures/49-multi-kingdom/49-2-procrustes-audit", 180, 80)
p49_2

associations$Status <- ifelse(
  associations$RawQValue < 0.05 & associations$AdjustedQValue < 0.05, "Both",
  ifelse(associations$RawQValue < 0.05, "Raw only",
         ifelse(associations$AdjustedQValue < 0.05, "Adjusted only", "Neither"))
)
status_colors <- c(
  Both = unname(pal_pub["blue"]),
  `Raw only` = unname(pal_pub["orange"]),
  `Adjusted only` = unname(pal_pub["green"]),
  Neither = "#CDD1D5"
)
p49_3 <- ggplot(associations, aes(RawRho, AdjustedRho, colour = Status)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "#E5E5E5", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#E5E5E5", linewidth = 0.3) +
  geom_point(size = 1.2, alpha = 0.62) +
  scale_color_manual(values = status_colors) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(x = "Raw Spearman rho", y = "Niche-adjusted Spearman rho", colour = "BH status",
       title = "Raw 505 to adjusted 2 FDR hits",
       subtitle = "Effect sign changed for 491 pairs") + theme_pub()
save_pub(p49_3, "figures/49-multi-kingdom/49-3-shared-niche-audit", 89, 82)
p49_3

bootstrapped <- associations[is.finite(associations$SignStability), ]
bootstrapped <- head(bootstrapped, 15)
bootstrapped$Association <- paste0(bootstrapped$Label1, "  |  ", bootstrapped$Label2)
bootstrapped$Association <- factor(
  bootstrapped$Association,
  levels = rev(bootstrapped$Association)
)
bootstrapped$Pass <- bootstrapped$AdjustedQValue < 0.05 & bootstrapped$SignStability >= 0.90
p49_4 <- ggplot(bootstrapped, aes(AdjustedRho, Association, colour = Pass)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = BootstrapLower, xmax = BootstrapUpper), height = 0.15,
                 linewidth = 0.55) +
  geom_point(aes(size = SignStability), alpha = 0.9) +
  scale_color_manual(values = c(
    `TRUE` = unname(pal_pub["blue"]),
    `FALSE` = unname(pal_pub["grey"])
  ),
                     labels = c(`TRUE` = "q < 0.05 and stability >= 90%", `FALSE` = "Exploratory")) +
  scale_size_continuous(range = c(2, 4), labels = scales::percent) +
  labs(x = "Niche-adjusted Spearman rho", y = NULL, colour = "Decision",
       size = "Sign stability", title = "Adjusted association stability") +
  theme_pub(8.5) + theme(legend.position = "right")
save_pub(p49_4, "figures/49-multi-kingdom/49-4-stable-candidates", 180, 112)
p49_4
