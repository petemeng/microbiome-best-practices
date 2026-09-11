# 微生物来源追踪：SourceTracker / FEAST 与传播证据
# Run sequentially in a new working directory.
# Required packages: FEAST, ggplot2, jsonlite, ragg, readr, scales, svglite.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "scripts/run_article50_feast.R",
  "scripts/run_article50_sourcetracker.py",
  "env/sourcetracker2.yml",
  "data/small/source-tracking-feast/metadata.tsv",
  "data/small/source-tracking-feast/otutab.tsv",
  "data/small/source-tracking-feast/source-summary.json",
  "data/small/source-tracking-feast/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)

set.seed(20260750)
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
data_dir <- "data/small/source-tracking-feast"
otutab_raw <- read_keyed_tsv(file.path(data_dir, "otutab.tsv"))
taxonomy <- read_keyed_tsv(file.path(data_dir, "taxonomy.tsv"))
metadata <- read_keyed_tsv(file.path(data_dir, "metadata.tsv"))
otutab <- as.matrix(data.frame(lapply(otutab_raw, as.integer), check.names = FALSE))
rownames(otutab) <- rownames(otutab_raw)
source_ids <- rownames(metadata)[metadata$SourceSink == "Source"]
sink_ids <- rownames(metadata)[metadata$SourceSink == "Sink"]
stopifnot(
  nrow(otutab) == 1839L, ncol(otutab) == 10L, sum(otutab) == 475154L,
  identical(colnames(otutab), rownames(metadata)),
  length(source_ids) == 9L, length(sink_ids) == 1L,
  identical(
    as.integer(table(metadata[source_ids, "SourceClass"])[c(
      "Adult gut", "Adult skin", "Infant gut", "Soil"
    )]),
    c(3L, 3L, 1L, 2L)
  ),
  all(taxonomy$TaxonomyStatus == "Anonymous taxon in the official FEAST example")
)
head(otutab[, 1:5])
head(taxonomy)
metadata

run_analysis <- function(command, args) {
  output <- system2(command, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop(paste(tail(output, 20), collapse = "\n"))
  invisible(output)
}
run_analysis(file.path(R.home("bin"), "Rscript"), c("--vanilla", "scripts/run_article50_feast.R"))
run_analysis(Sys.getenv("SOURCETRACKER_PYTHON", "sourcetracker-env/bin/python"), c("scripts/run_article50_sourcetracker.py"))

result_dir <- "results/50-source-tracking"
feast_main <- readr::read_tsv(
  file.path(result_dir, "feast-main-class.tsv"), show_col_types = FALSE
)
feast_loo <- readr::read_tsv(
  file.path(result_dir, "feast-loo-predictions.tsv"), show_col_types = FALSE
)
feast_missing <- readr::read_tsv(
  file.path(result_dir, "feast-missing-source.tsv"), show_col_types = FALSE
)
source_tracker_main <- readr::read_tsv(
  file.path(result_dir, "source-tracker-main.tsv"), show_col_types = FALSE
)
source_tracker_loo <- readr::read_tsv(
  file.path(result_dir, "source-tracker-loo-predictions.tsv"), show_col_types = FALSE
)
feast_summary <- jsonlite::read_json(
  file.path(result_dir, "summary-feast.json"), simplifyVector = TRUE
)
source_tracker_summary <- jsonlite::read_json(
  file.path(result_dir, "summary-source-tracker.json"), simplifyVector = TRUE
)
stopifnot(
  abs(sum(feast_main$Contribution) - 1) < 1e-6,
  abs(sum(source_tracker_main$Contribution) - 1) < 1e-6,
  feast_summary$leave_one_out_known_class_accuracy == 0.5,
  source_tracker_summary$leave_one_out_known_class_accuracy == 0.5,
  feast_summary$main_largest_source_class == "Infant gut",
  source_tracker_summary$main_largest_source_class == "Infant gut"
)
rbind(
  data.frame(Method = "FEAST", Unknown = feast_summary$main_unknown,
             LOOAccuracy = feast_summary$leave_one_out_known_class_accuracy),
  data.frame(Method = "SourceTracker2", Unknown = source_tracker_summary$main_unknown,
             LOOAccuracy = source_tracker_summary$leave_one_out_known_class_accuracy)
)

source_colors <- c(
  `Infant gut` = unname(pal_pub["blue"]),
  `Adult gut` = unname(pal_pub["orange"]),
  `Adult skin` = unname(pal_pub["purple"]),
  Soil = unname(pal_pub["green"]),
  Unknown = unname(pal_pub["grey"])
)
feast_main$SourceClass <- factor(
  feast_main$SourceClass,
  levels = names(source_colors)
)
p50_1 <- ggplot(feast_main, aes(SourceClass, Contribution, fill = SourceClass)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = scales::percent(Contribution, accuracy = 0.1)),
            vjust = -0.35, family = font_pub, size = 3) +
  scale_fill_manual(values = source_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.78),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Estimated contribution", fill = "Source class",
       title = "FEAST source attribution",
       subtitle = "Sink: infant gut 1; 5,000-read analysis depth") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                      legend.position = "none")
save_pub(p50_1, "figures/50-source-tracking/50-1-feast-main", 89, 78)
p50_1

method_comparison <- rbind(
  data.frame(Method = "FEAST", SourceClass = as.character(feast_main$SourceClass),
             Contribution = feast_main$Contribution, MonteCarloSD = NA_real_),
  data.frame(Method = "SourceTracker2", SourceClass = source_tracker_main$SourceClass,
             Contribution = source_tracker_main$Contribution,
             MonteCarloSD = source_tracker_main$MonteCarloSD)
)
method_comparison$SourceClass <- factor(method_comparison$SourceClass, levels = names(source_colors))
p50_2 <- ggplot(method_comparison, aes(SourceClass, Contribution, fill = Method)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.66) +
  geom_errorbar(
    aes(ymin = pmax(0, Contribution - MonteCarloSD),
        ymax = pmin(1, Contribution + MonteCarloSD)),
    position = position_dodge(width = 0.74), width = 0.16, na.rm = TRUE
  ) +
  scale_fill_manual(values = c(
    FEAST = unname(pal_pub["blue"]),
    SourceTracker2 = unname(pal_pub["orange"])
  )) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Estimated contribution", fill = "Method",
       title = "Algorithm-dependent allocation",
       subtitle = "SourceTracker2 bars show Monte Carlo SD, not biological uncertainty") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_pub(p50_2, "figures/50-source-tracking/50-2-method-comparison", 180, 82)
p50_2

feast_loo_long <- reshape(
  as.data.frame(feast_loo[, c("SampleID", "TrueClass", "InfantGut", "AdultGut",
                              "AdultSkin", "Soil", "Unknown")]),
  varying = c("InfantGut", "AdultGut", "AdultSkin", "Soil", "Unknown"),
  v.names = "Contribution", timevar = "SourceClass",
  times = c("Infant gut", "Adult gut", "Adult skin", "Soil", "Unknown"),
  direction = "long"
)
feast_loo_long$Method <- "FEAST"
st_loo_long <- reshape(
  as.data.frame(source_tracker_loo[, c("SampleID", "TrueClass", "Infantgut",
                                       "Adultgut", "Adultskin", "Soil", "Unknown")]),
  varying = c("Infantgut", "Adultgut", "Adultskin", "Soil", "Unknown"),
  v.names = "Contribution", timevar = "SourceClass",
  times = c("Infant gut", "Adult gut", "Adult skin", "Soil", "Unknown"),
  direction = "long"
)
st_loo_long$Method <- "SourceTracker2"
loo_long <- rbind(feast_loo_long, st_loo_long)
loo_long$SampleLabel <- paste0(loo_long$SampleID, " (", loo_long$TrueClass, ")")
loo_long$SourceClass <- factor(loo_long$SourceClass, levels = names(source_colors))
p50_3 <- ggplot(loo_long, aes(SourceClass, SampleLabel, fill = Contribution)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  facet_wrap(~ Method, ncol = 1) +
  scale_fill_gradient(low = "#F7FBFF", high = pal_pub["blue"], labels = scales::percent) +
  labs(x = "Candidate source class", y = "Held-out sample (true class)",
       fill = "Contribution", title = "Leave-one-out calibration",
       subtitle = "Known-class top-1 accuracy: 50% for both methods") +
  theme_pub(8.5) + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                         legend.position = "right")
save_pub(p50_3, "figures/50-source-tracking/50-3-loo-calibration", 180, 135)
p50_3

scenario_order <- c("None", "Infant gut", "Adult gut", "Adult skin", "Soil")
feast_missing$RemovedClass <- factor(feast_missing$RemovedClass, levels = scenario_order)
feast_missing$SourceClass <- factor(feast_missing$SourceClass, levels = names(source_colors))
p50_4 <- ggplot(feast_missing, aes(RemovedClass, Contribution, fill = SourceClass)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = source_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0))) +
  labs(x = "Removed source class", y = "Estimated contribution", fill = "Allocated class",
       title = "Missing-source sensitivity",
       subtitle = "Unknown is model-dependent and need not increase when a true source is omitted") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                      legend.position = "right")
save_pub(p50_4, "figures/50-source-tracking/50-4-missing-source", 180, 90)
p50_4
