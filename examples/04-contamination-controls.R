# 污染与对照：blank、kitome、低生物量、decontam
# Run sequentially in a new working directory.
# Required packages: decontam, ggplot2, scales.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/decontam/metadata.tsv",
  "data/small/decontam/otutab.tsv",
  "data/small/decontam/source_summary.json",
  "data/small/decontam/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(decontam)
library(ggplot2)

stopifnot(
  identical(
    as.character(utils::packageVersion("decontam")),
    "1.24.0"
  )
)

set.seed(20260718)

pal_pub <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#F0E442", "#999999"
)

scale_color_pub <- function(...) {
  ggplot2::scale_color_manual(values = pal_pub, ...)
}

scale_fill_pub <- function(...) {
  ggplot2::scale_fill_manual(values = pal_pub, ...)
}

theme_pub <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.7),
        color = NA
      )
    )
}

save_pub <- function(
  plot,
  file_base,
  width = 120,
  height = 95,
  units = "mm",
  dpi = 300
) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(file_base, ".pdf"),
    plot,
    width = width,
    height = height,
    units = units,
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    paste0(file_base, ".png"),
    plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(file_base, ".tiff"),
    plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
  invisible(plot)
}

otutab <- read.delim(
  "data/small/decontam/otutab.tsv",
  row.names = 1,
  check.names = FALSE
)
taxonomy <- read.delim(
  "data/small/decontam/taxonomy.tsv",
  row.names = 1,
  check.names = FALSE
)
metadata <- read.delim(
  "data/small/decontam/metadata.tsv",
  row.names = 1,
  check.names = FALSE
)

otutab <- as.matrix(otutab)
storage.mode(otutab) <- "numeric"
metadata$IsNegative <- as.logical(metadata$IsNegative)

stopifnot(
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata)),
  all(otutab >= 0),
  all(otutab == floor(otutab)),
  all(colSums(otutab) > 0),
  all(c(
    "IsNegative",
    "quant_reading",
    "PlateNumber",
    "Sample_or_Control"
  ) %in% colnames(metadata)),
  all(is.finite(metadata$quant_reading)),
  all(metadata$quant_reading > 0)
)

dim(otutab)
otutab[1:5, 1:5]
taxonomy[1:5, , drop = FALSE]
head(metadata)

control_by_plate <- with(
  metadata,
  table(PlateNumber, IsNegative)
)
control_by_plate
stopifnot(all(control_by_plate[, "TRUE"] == 5))

sample_qc <- data.frame(
  SampleID = rownames(metadata),
  SampleClass = factor(
    ifelse(
      metadata$IsNegative,
      "Negative control",
      "Biological sample"
    ),
    levels = c("Biological sample", "Negative control")
  ),
  Plate = factor(metadata$PlateNumber),
  DNAConcentration = metadata$quant_reading,
  LibrarySize = colSums(otutab),
  stringsAsFactors = FALSE
)

library_summary <- do.call(
  rbind,
  lapply(
    split(sample_qc$LibrarySize, sample_qc$SampleClass),
    function(x) {
      data.frame(
        N = length(x),
        Minimum = min(x),
        Median = stats::median(x),
        Maximum = max(x)
      )
    }
  )
)
library_summary

seqtab <- t(otutab)
is_negative <- metadata[rownames(seqtab), "IsNegative"]
dna_concentration <- metadata[
  rownames(seqtab),
  "quant_reading"
]

stopifnot(
  identical(rownames(seqtab), rownames(metadata)),
  identical(colnames(seqtab), rownames(taxonomy))
)

contam_frequency <- decontam::isContaminant(
  seqtab,
  method = "frequency",
  conc = dna_concentration,
  threshold = 0.10
)

contam_prevalence <- decontam::isContaminant(
  seqtab,
  method = "prevalence",
  neg = is_negative,
  threshold = 0.10
)

contam_combined <- decontam::isContaminant(
  seqtab,
  method = "combined",
  conc = dna_concentration,
  neg = is_negative,
  threshold = 0.10
)

classification_summary <- data.frame(
  Method = c("Frequency", "Prevalence", "Combined"),
  Threshold = 0.10,
  ContaminantASVs = c(
    sum(contam_frequency$contaminant),
    sum(contam_prevalence$contaminant),
    sum(contam_combined$contaminant)
  )
)
classification_summary

method_overlap <- table(
  Frequency = contam_frequency$contaminant,
  Prevalence = contam_prevalence$contaminant
)
method_overlap

contam_batch <- decontam::isContaminant(
  seqtab,
  method = "combined",
  conc = dna_concentration,
  neg = is_negative,
  batch = factor(metadata[rownames(seqtab), "PlateNumber"]),
  batch.combine = "minimum",
  threshold = 0.10
)

batch_comparison <- data.frame(
  GlobalCombined = contam_combined$contaminant,
  BatchCombined = contam_batch$contaminant
)

table(batch_comparison)
data.frame(
  Global = sum(batch_comparison$GlobalCombined),
  BatchAware = sum(batch_comparison$BatchCombined),
  Shared = sum(
    batch_comparison$GlobalCombined &
      batch_comparison$BatchCombined
  )
)

negative_prevalence <- colMeans(
  seqtab[is_negative, , drop = FALSE] > 0
)
biological_prevalence <- colMeans(
  seqtab[!is_negative, , drop = FALSE] > 0
)
total_feature_reads <- colSums(seqtab)

contaminant_audit <- data.frame(
  FeatureID = rownames(contam_combined),
  taxonomy[rownames(contam_combined), , drop = FALSE],
  NegativePrevalence = negative_prevalence,
  BiologicalPrevalence = biological_prevalence,
  TotalReads = total_feature_reads,
  FrequencyP = contam_combined$p.freq,
  PrevalenceP = contam_combined$p.prev,
  CombinedP = contam_combined$p,
  GlobalContaminant = contam_combined$contaminant,
  BatchContaminant = contam_batch$contaminant,
  check.names = FALSE
)

contaminant_ids <- contaminant_audit$FeatureID[
  contaminant_audit$GlobalContaminant
]
library_size <- colSums(otutab)
contaminant_reads <- colSums(
  otutab[contaminant_ids, , drop = FALSE]
)

sample_burden <- data.frame(
  SampleID = colnames(otutab),
  SampleClass = sample_qc$SampleClass,
  Plate = sample_qc$Plate,
  DNAConcentration = sample_qc$DNAConcentration,
  LibrarySize = library_size,
  ContaminantReads = contaminant_reads,
  ContaminantFraction = contaminant_reads / library_size,
  stringsAsFactors = FALSE
)

burden_summary <- do.call(
  rbind,
  lapply(
    split(
      sample_burden$ContaminantFraction,
      sample_burden$SampleClass
    ),
    function(x) {
      data.frame(
        Minimum = min(x),
        Median = stats::median(x),
        Mean = mean(x),
        Maximum = max(x)
      )
    }
  )
)
burden_summary

burden_correlation <- stats::cor.test(
  sample_burden$ContaminantFraction[
    sample_burden$SampleClass == "Biological sample"
  ],
  sample_burden$DNAConcentration[
    sample_burden$SampleClass == "Biological sample"
  ],
  method = "spearman",
  exact = FALSE
)
burden_correlation

biological_ids <- rownames(metadata)[!metadata$IsNegative]
kept_features <- setdiff(rownames(otutab), contaminant_ids)

otutab_filtered <- otutab[
  kept_features,
  biological_ids,
  drop = FALSE
]
taxonomy_filtered <- taxonomy[
  kept_features,
  ,
  drop = FALSE
]
metadata_filtered <- metadata[
  biological_ids,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(otutab_filtered),
    rownames(taxonomy_filtered)
  ),
  identical(
    colnames(otutab_filtered),
    rownames(metadata_filtered)
  ),
  all(colSums(otutab_filtered) > 0)
)

biological_reads_before <- sum(
  otutab[, biological_ids, drop = FALSE]
)
biological_reads_after <- sum(otutab_filtered)
biological_reads_removed <- (
  1 - biological_reads_after / biological_reads_before
)

filter_summary <- data.frame(
  InputFeatures = nrow(otutab),
  RemovedFeatures = length(contaminant_ids),
  OutputFeatures = nrow(otutab_filtered),
  BiologicalSamples = ncol(otutab_filtered),
  BiologicalReadsRemoved = biological_reads_removed
)
filter_summary

write_result_tsv <- function(x, id_name, path) {
  out <- data.frame(
    setNames(list(rownames(x)), id_name),
    x,
    check.names = FALSE
  )
  utils::write.table(
    out,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

result_dir <- "results/04-contamination"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.table(
  contaminant_audit,
  file.path(result_dir, "contaminant-classification.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
utils::write.table(
  sample_burden,
  file.path(result_dir, "sample-contamination-burden.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
write_result_tsv(
  otutab_filtered,
  "FeatureID",
  file.path(result_dir, "otutab.filtered.tsv")
)
write_result_tsv(
  taxonomy_filtered,
  "FeatureID",
  file.path(result_dir, "taxonomy.filtered.tsv")
)
write_result_tsv(
  metadata_filtered,
  "SampleID",
  file.path(result_dir, "metadata.filtered.tsv")
)

sample_colors <- c(
  "Biological sample" = "#0072B2",
  "Negative control" = "#D55E00"
)

p_control_library <- ggplot(
  sample_qc,
  aes(SampleClass, LibrarySize, fill = SampleClass)
) +
  geom_violin(
    width = 0.72,
    trim = FALSE,
    alpha = 0.25,
    color = NA
  ) +
  geom_boxplot(
    width = 0.24,
    outlier.shape = NA,
    alpha = 0.75,
    linewidth = 0.45
  ) +
  geom_point(
    position = position_jitter(
      width = 0.12,
      height = 0,
      seed = 20260718
    ),
    shape = 21,
    size = 1.35,
    stroke = 0.2,
    alpha = 0.48
  ) +
  scale_fill_manual(
    values = sample_colors,
    guide = "none"
  ) +
  scale_y_log10(
    breaks = c(100, 300, 1000, 3000, 10000, 30000),
    labels = scales::label_number(big.mark = ",")
  ) +
  labs(
    title = "Negative controls carry fewer reads but remain informative",
    subtitle = paste0(
      sum(!metadata$IsNegative),
      " biological samples · ",
      sum(metadata$IsNegative),
      " controls · ",
      length(unique(metadata$PlateNumber)),
      " plates"
    ),
    x = NULL,
    y = "Library size (reads; log scale)"
  ) +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey35", size = 9.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

p_control_library
save_pub(
  p_control_library,
  "figures/04-control-library-size",
  width = 162,
  height = 112
)

prevalence_plot_df <- data.frame(
  FeatureID = contaminant_audit$FeatureID,
  NegativePrevalence = contaminant_audit$NegativePrevalence,
  BiologicalPrevalence = contaminant_audit$BiologicalPrevalence,
  TotalReads = contaminant_audit$TotalReads,
  Classification = factor(
    ifelse(
      contaminant_audit$GlobalContaminant,
      "Likely contaminant",
      "Not classified"
    ),
    levels = c("Not classified", "Likely contaminant")
  )
)

classification_colors <- c(
  "Not classified" = "#B8B8B8",
  "Likely contaminant" = "#D55E00"
)

label_features <- intersect(
  c("Seq30", "Seq175", "Seq3"),
  prevalence_plot_df$FeatureID
)

p_contaminant_prevalence <- ggplot(
  prevalence_plot_df,
  aes(
    NegativePrevalence,
    BiologicalPrevalence,
    color = Classification,
    size = TotalReads
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.55
  ) +
  geom_point(alpha = 0.62) +
  geom_text(
    data = subset(
      prevalence_plot_df,
      FeatureID %in% label_features
    ),
    aes(label = FeatureID),
    size = 3.0,
    color = "grey10",
    nudge_y = 0.035,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = classification_colors,
    name = "Combined classifier"
  ) +
  scale_size_continuous(
    trans = "log10",
    range = c(0.8, 4.4),
    breaks = c(10, 1000, 100000, 1000000),
    labels = scales::label_number(big.mark = ","),
    name = "Total reads"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_equal(clip = "off") +
  labs(
    title = "Negative controls reveal a distinct contaminant branch",
    subtitle = paste0(
      sum(contam_combined$contaminant),
      " of ",
      nrow(contam_combined),
      " ASVs classified by combined evidence at threshold = 0.10"
    ),
    x = "Prevalence in negative controls",
    y = "Prevalence in biological samples"
  ) +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey35", size = 9.2),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(8, 10, 8, 8)
  )

p_contaminant_prevalence
save_pub(
  p_contaminant_prevalence,
  "figures/04-contaminant-prevalence",
  width = 178,
  height = 126
)

rho_label <- sprintf(
  "Biological samples: Spearman rho = %.2f",
  unname(burden_correlation$estimate)
)

p_contaminant_burden <- ggplot(
  sample_burden,
  aes(
    DNAConcentration,
    ContaminantFraction,
    color = SampleClass
  )
) +
  geom_point(
    alpha = 0.62,
    size = 1.8
  ) +
  scale_color_manual(
    values = sample_colors,
    name = "Sample class"
  ) +
  scale_x_log10(
    labels = scales::label_number(big.mark = ",")
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(
      base = 10,
      sigma = 0.0001
    ),
    breaks = c(0, 0.0001, 0.001, 0.01, 0.1, 1),
    labels = scales::percent_format(accuracy = 0.01),
    limits = c(0, 1)
  ) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = rho_label,
    hjust = 1.04,
    vjust = 1.25,
    size = 3.1,
    label.size = 0.25,
    fill = scales::alpha("white", 0.85),
    color = "grey15"
  ) +
  labs(
    title = "Contaminant burden rises as DNA concentration falls",
    subtitle = paste0(
      "Burden is the read fraction assigned to ",
      sum(contam_combined$contaminant),
      " globally classified ASVs"
    ),
    x = "DNA concentration (PicoGreen intensity; log scale)",
    y = "Candidate contaminant read fraction"
  ) +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey35", size = 9.2),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

p_contaminant_burden
save_pub(
  p_contaminant_burden,
  "figures/04-contaminant-burden",
  width = 174,
  height = 116
)
