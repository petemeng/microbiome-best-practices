# 相关整合：Procrustes/Mantel/HAllA
# Run sequentially in a new working directory.
# Required packages: ggplot2, jsonlite, permute, ragg, readr, svglite, vegan.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "scripts/run_article46_global_concordance.R",
  "scripts/run_article46_halla.py",
  "env/multiomics.yml",
  "data/small/paired-ibd-multiomics/metabolite-annotation.tsv",
  "data/small/paired-ibd-multiomics/metabolites.tsv",
  "data/small/paired-ibd-multiomics/metadata.tsv",
  "data/small/paired-ibd-multiomics/otutab.tsv",
  "data/small/paired-ibd-multiomics/source-summary.json",
  "data/small/paired-ibd-multiomics/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)

set.seed(20260746)
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
data_dir <- "data/small/paired-ibd-multiomics"
otutab_raw <- read_keyed_tsv(file.path(data_dir, "otutab.tsv"))
taxonomy <- read_keyed_tsv(file.path(data_dir, "taxonomy.tsv"))
metadata <- read_keyed_tsv(file.path(data_dir, "metadata.tsv"))
metabolites_raw <- read_keyed_tsv(file.path(data_dir, "metabolites.tsv"))
metabolite_annotation <- read_keyed_tsv(
  file.path(data_dir, "metabolite-annotation.tsv")
)
otutab <- as.matrix(data.frame(lapply(otutab_raw, as.numeric), check.names = FALSE))
metabolites <- as.matrix(data.frame(
  lapply(metabolites_raw, as.numeric), check.names = FALSE
))
rownames(otutab) <- rownames(otutab_raw)
rownames(metabolites) <- rownames(metabolites_raw)
stopifnot(
  nrow(otutab) == 250L, ncol(otutab) == 220L,
  nrow(metabolites) == 277L, ncol(metabolites) == 220L,
  identical(colnames(otutab), rownames(metadata)),
  identical(colnames(metabolites), rownames(metadata)),
  identical(
    as.integer(table(metadata$StudyGroup)[c("CD", "Control", "UC")]),
    c(88L, 56L, 76L)
  )
)
head(otutab[, 1:5])
head(taxonomy)
head(metadata)
head(metabolites[, 1:5])

run_analysis <- function(command, args) {
  output <- system2(command, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop(paste(tail(output, 20), collapse = "\n"))
  invisible(output)
}
run_analysis(file.path(R.home("bin"), "Rscript"), c("--vanilla", "scripts/run_article46_global_concordance.R"))
run_analysis(Sys.getenv("MULTIOMICS_PYTHON", "multiomics-env/bin/python"), c("scripts/run_article46_halla.py"))

result_dir <- "results/46-integration"
global_tests <- readr::read_tsv(
  file.path(result_dir, "global-concordance.tsv"), show_col_types = FALSE
)
procrustes_arrows <- readr::read_tsv(
  file.path(result_dir, "procrustes-arrows.tsv"), show_col_types = FALSE
)
stopifnot(nrow(global_tests) == 4L, all(global_tests$QValue <= 0.001))
global_tests

raw_pairs <- readr::read_tsv(
  file.path(result_dir, "pairwise-raw.tsv"), show_col_types = FALSE
)
adjusted_pairs <- readr::read_tsv(
  file.path(result_dir, "pairwise-group-residual.tsv"), show_col_types = FALSE
)
halla_associations <- readr::read_tsv(
  file.path(result_dir, "halla-v1", "all_associations.txt"),
  show_col_types = FALSE
)
halla_blocks <- readr::read_tsv(
  file.path(result_dir, "halla-v1", "sig_clusters.txt"),
  show_col_types = FALSE
)
halla_summary <- jsonlite::read_json(
  file.path(result_dir, "summary.json"), simplifyVector = TRUE
)
stopifnot(
  nrow(raw_pairs) == 1000L, nrow(adjusted_pairs) == 1000L,
  sum(raw_pairs$QValue < 0.05) == 437L,
  sum(adjusted_pairs$QValue < 0.05) == 350L,
  nrow(halla_blocks) == 70L
)
head(halla_blocks, 5)

group_colors <- c(
  Control = unname(pal_pub["blue"]),
  CD = unname(pal_pub["vermillion"]),
  UC = unname(pal_pub["orange"])
)
p46_1 <- ggplot(procrustes_arrows) +
  geom_segment(
    aes(x = StartAxis1, y = StartAxis2, xend = EndAxis1, yend = EndAxis2,
        colour = StudyGroup),
    linewidth = 0.25, alpha = 0.28
  ) +
  geom_point(aes(StartAxis1, StartAxis2, colour = StudyGroup), size = 0.8, alpha = 0.75) +
  geom_point(aes(EndAxis1, EndAxis2, colour = StudyGroup), shape = 1, size = 0.9) +
  facet_wrap(~ Model) +
  scale_color_manual(values = group_colors) +
  labs(
    x = "Procrustes axis 1", y = "Procrustes axis 2", colour = "Study group",
    title = "Microbiome–metabolome configuration",
    subtitle = "Filled: microbiome; open: rotated metabolome"
  ) +
  coord_equal() + theme_pub()
save_pub(p46_1, "figures/46-integration/46-1-procrustes", 180, 78)
p46_1

global_tests$YPosition <- ifelse(global_tests$Method == "Mantel", 1, 2) +
  ifelse(global_tests$Model == "Group-adjusted", -0.09, 0.09)
p46_2 <- ggplot(global_tests, aes(Statistic, YPosition, colour = Model)) +
  geom_segment(aes(x = 0, xend = Statistic, yend = YPosition), linewidth = 0.6, alpha = 0.5) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.3f", Statistic)), hjust = -0.18,
            size = 2.8, family = font_pub, show.legend = FALSE) +
  scale_color_manual(values = c(
    Observed = unname(pal_pub["grey"]),
    `Group-adjusted` = unname(pal_pub["blue"])
  )) +
  scale_x_continuous(limits = c(0, 0.68), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = c(1, 2), labels = c("Mantel", "Procrustes"),
                     limits = c(0.65, 2.35)) +
  labs(x = "Correlation statistic", y = NULL, colour = "Model",
       title = "Global concordance", subtitle = "999 permutations; all BH q = 0.001") +
  theme_pub()
save_pub(p46_2, "figures/46-integration/46-2-global-tests", 89, 66)
p46_2

pair_audit <- merge(
  raw_pairs[, c("Microbe", "Metabolite", "Rho", "QValue")],
  adjusted_pairs[, c("Microbe", "Metabolite", "Rho", "QValue")],
  by = c("Microbe", "Metabolite"), suffixes = c("Raw", "Adjusted")
)
pair_audit$Status <- ifelse(
  pair_audit$QValueRaw < 0.05 & pair_audit$QValueAdjusted < 0.05, "Both",
  ifelse(pair_audit$QValueRaw < 0.05, "Raw only",
         ifelse(pair_audit$QValueAdjusted < 0.05, "Adjusted only", "Neither"))
)
status_colors <- c(
  Both = unname(pal_pub["blue"]),
  `Raw only` = unname(pal_pub["orange"]),
  `Adjusted only` = unname(pal_pub["green"]),
  Neither = "#C9CDD3"
)
p46_3 <- ggplot(pair_audit, aes(RhoRaw, RhoAdjusted, colour = Status)) +
  geom_abline(slope = 1, intercept = 0, colour = "#777777", linetype = 2, linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
  geom_point(size = 1.2, alpha = 0.65) +
  scale_color_manual(values = status_colors) +
  labs(
    x = "Raw Spearman rho", y = "Group-adjusted Spearman rho", colour = "BH status",
    title = "Group confounding audit",
    subtitle = "437 raw hits; 350 adjusted hits"
  ) +
  coord_equal(xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7)) + theme_pub()
save_pub(p46_3, "figures/46-integration/46-3-confounding-audit", 89, 82)
p46_3

top_halla <- head(adjusted_pairs, 15)
top_halla$Association <- paste0(
  sub(" \\[MB_.*$", "", top_halla$Microbe), "  |  ",
  sub(" \\[MT_.*$", "", top_halla$Metabolite)
)
top_halla$Association <- factor(top_halla$Association, levels = rev(top_halla$Association))
p46_4 <- ggplot(top_halla, aes(Rho, Association, colour = Rho > 0)) +
  geom_segment(aes(x = 0, xend = Rho, yend = Association), linewidth = 0.65) +
  geom_point(aes(size = -log10(QValue)), alpha = 0.9) +
  scale_color_manual(values = c(
    `TRUE` = unname(pal_pub["blue"]),
    `FALSE` = unname(pal_pub["vermillion"])
  ),
                     labels = c(`TRUE` = "Positive", `FALSE` = "Negative")) +
  scale_size_continuous(range = c(2, 5)) +
  labs(x = "Adjusted Spearman rho", y = NULL, colour = "Direction",
       size = expression(-log[10](q)), title = "Strongest adjusted associations") +
  theme_pub(9) + theme(legend.position = "right")
save_pub(p46_4, "figures/46-integration/46-4-halla-associations", 180, 105)
p46_4
