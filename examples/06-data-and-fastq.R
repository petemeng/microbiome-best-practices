# 认识你的数据：测序原理 + FASTQ / 文件格式拆解
# Run sequentially in a new working directory.
# Required packages: dplyr, ggplot2, knitr, scales, tidyr.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/fastq/barcodes.fastq.gz",
  "data/small/fastq/forward.fastq.gz",
  "data/small/fastq/metadata.tsv",
  "data/small/fastq/reverse.fastq.gz",
  "data/small/fastq/source_summary.json"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)
library(dplyr)
library(tidyr)

stopifnot(
  identical(
    as.character(utils::packageVersion("ggplot2")),
    "3.5.2"
  ),
  identical(
    as.character(utils::packageVersion("dplyr")),
    "1.1.4"
  ),
  identical(
    as.character(utils::packageVersion("tidyr")),
    "1.3.1"
  ),
  identical(
    as.character(utils::packageVersion("jsonlite")),
    "1.8.8"
  ),
  identical(
    as.character(utils::packageVersion("digest")),
    "0.6.36"
  )
)

set.seed(20260719)

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
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.3
      ),
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
  dir.create(
    dirname(file_base),
    recursive = TRUE,
    showWarnings = FALSE
  )
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
}

fastq_dir <- "data/small/fastq"

fastq_paths <- c(
  Forward = file.path(fastq_dir, "forward.fastq.gz"),
  Reverse = file.path(fastq_dir, "reverse.fastq.gz"),
  Barcode = file.path(fastq_dir, "barcodes.fastq.gz"),
  Metadata = file.path(fastq_dir, "metadata.tsv")
)

stopifnot(all(file.exists(fastq_paths)))

known_sha256 <- data.frame(
  File = names(fastq_paths),
  ExcerptSHA256 = c(
    "e6fabdbd31db519c719447f1caf2d1cd334c6ab932353e534ea2433ead88d254",
    "d2fa5841d150fead5a73067a267ec531221de11c512a38ec32ada03f1621bf65",
    "9cce1bb9a57975d14b54a5cf07c89b2b7c2095da4a8f0d1bd2720f1349e74057",
    "8d5342f0a80abf4197088bb1ce92f4903c1bd1212f6a9182e0e1c87ebc322bcb"
  ),
  Full1PercentSHA256 = c(
    "392f178b5f15967fafb9332c0493a8b6ecbf4a7857538023813a630c488c2702",
    "5407370313a3974e5b70878153841a61df408e49ced133b3ace44b2db286cfae",
    "c07f517e920d9e3924ce1e87102f057232b57041947523b90d560f68a0ff3c0f",
    "7cff810ad86a621ebc78a16b690255d839ea58e2622ad45a11a838d0f8c5b3bd"
  ),
  stringsAsFactors = FALSE
)

observed_sha256 <- vapply(
  fastq_paths,
  digest::digest,
  FUN.VALUE = character(1),
  algo = "sha256",
  serialize = FALSE,
  file = TRUE
)

known_sha256$ObservedSHA256 <- unname(observed_sha256)
known_sha256$SourceSet <- ifelse(
  known_sha256$ObservedSHA256 == known_sha256$ExcerptSHA256,
  "Bundled 2,000-record excerpt",
  ifelse(
    known_sha256$ObservedSHA256 ==
      known_sha256$Full1PercentSHA256,
    "Official 1% source",
    "Unknown"
  )
)

stopifnot(!any(known_sha256$SourceSet == "Unknown"))

knitr::kable(
  known_sha256[, c("File", "SourceSet", "ObservedSHA256")],
  caption = "Input identity audit"
)

read_fastq <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  lines <- readLines(
    connection,
    warn = FALSE,
    encoding = "bytes"
  )

  if (length(lines) %% 4 != 0) {
    stop(
      basename(path),
      " does not contain a multiple of four lines"
    )
  }

  records <- matrix(
    lines,
    ncol = 4,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c("Header", "Sequence", "Separator", "Quality")
    )
  )
  records <- as.data.frame(
    records,
    stringsAsFactors = FALSE
  )

  if (!all(startsWith(records$Header, "@"))) {
    stop(basename(path), " contains an invalid header")
  }
  if (!all(startsWith(records$Separator, "+"))) {
    stop(basename(path), " contains an invalid separator")
  }
  if (!all(
    nchar(records$Sequence, type = "bytes") ==
      nchar(records$Quality, type = "bytes")
  )) {
    stop(
      basename(path),
      " contains unequal sequence and quality lengths"
    )
  }

  records$ReadID <- sub(
    "^@([^ ]+).*$",
    "\\1",
    records$Header
  )
  records$Mate <- sub(
    "^@[^ ]+ ([12]):.*$",
    "\\1",
    records$Header
  )
  records
}

forward <- read_fastq(fastq_paths[["Forward"]])
reverse <- read_fastq(fastq_paths[["Reverse"]])
barcodes <- read_fastq(fastq_paths[["Barcode"]])

stopifnot(
  nrow(forward) == nrow(reverse),
  nrow(forward) == nrow(barcodes),
  identical(forward$ReadID, reverse$ReadID),
  identical(forward$ReadID, barcodes$ReadID),
  all(forward$Mate == "1"),
  all(reverse$Mate == "2")
)

record_example <- data.frame(
  Line = 1:4,
  Field = c(
    "Header", "Sequence",
    "Separator", "Quality"
  ),
  Value = unname(unlist(
    forward[1, c(
      "Header", "Sequence",
      "Separator", "Quality"
    )]
  )),
  stringsAsFactors = FALSE
)

record_example$Length <- nchar(
  record_example$Value,
  type = "bytes"
)

knitr::kable(
  transform(
    record_example,
    Value = ifelse(
      nchar(Value) > 72,
      paste0(substr(Value, 1, 69), "..."),
      ifelse(
        startsWith(Value, "@"),
        paste0("\\", Value),
        Value
      )
    )
  ),
  caption = "First real FASTQ record"
)

read_q2_metadata <- function(path) {
  lines <- readLines(path, warn = FALSE)
  stopifnot(
    length(lines) >= 3,
    startsWith(lines[2], "#q2:types")
  )
  connection <- textConnection(
    c(lines[1], lines[-c(1, 2)])
  )
  on.exit(close(connection), add = TRUE)
  read.delim(
    connection,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = ""
  )
}

metadata <- read_q2_metadata(
  fastq_paths[["Metadata"]]
)

metadata_summary <- metadata |>
  count(
    Transect = .data[["transect-name"]],
    name = "Samples"
  ) |>
  arrange(Transect)

stopifnot(
  nrow(metadata) == 75,
  n_distinct(metadata[["sample-id"]]) == 75,
  n_distinct(metadata[["barcode-sequence"]]) == 75,
  all(nchar(metadata[["barcode-sequence"]]) == 12)
)

knitr::kable(
  head(metadata, 5),
  caption = "First five Atacama metadata rows"
)

knitr::kable(
  metadata_summary,
  caption = "Samples by transect"
)

decode_phred33 <- function(quality_string) {
  utf8ToInt(quality_string) - 33L
}

quality_matrix <- function(records) {
  lengths <- nchar(
    records$Quality,
    type = "bytes"
  )
  if (length(unique(lengths)) != 1) {
    stop("variable read lengths require padding before matrix conversion")
  }
  do.call(
    rbind,
    lapply(records$Quality, decode_phred33)
  )
}

summarize_cycles <- function(records, read_name) {
  quality <- quality_matrix(records)
  tibble(
    Read = read_name,
    Cycle = seq_len(ncol(quality)),
    Q10 = apply(
      quality,
      2,
      stats::quantile,
      probs = 0.10,
      names = FALSE,
      type = 8
    ),
    MedianQ = apply(
      quality,
      2,
      stats::median
    ),
    Q90 = apply(
      quality,
      2,
      stats::quantile,
      probs = 0.90,
      names = FALSE,
      type = 8
    ),
    Q30Fraction = colMeans(quality >= 30)
  )
}

summarize_file <- function(records, read_name) {
  quality <- quality_matrix(records)
  sequence_length <- nchar(
    records$Sequence,
    type = "bytes"
  )
  tibble(
    Read = read_name,
    Records = nrow(records),
    MinLength = min(sequence_length),
    MedianLength = stats::median(sequence_length),
    MaxLength = max(sequence_length),
    MedianReadMeanPhred = stats::median(
      rowMeans(quality)
    ),
    Q30BaseFraction = mean(quality >= 30),
    NBaseFraction = NA_real_
  )
}

cycle_summary <- bind_rows(
  summarize_cycles(forward, "Forward"),
  summarize_cycles(reverse, "Reverse")
)

file_summary <- bind_rows(
  summarize_file(forward, "Forward"),
  summarize_file(reverse, "Reverse"),
  summarize_file(barcodes, "Barcode")
)

# NBaseFraction 用明确的逐字符计数重算，避免零匹配的特殊值。
count_n <- function(records) {
  n_bases <- sum(
    vapply(
      strsplit(records$Sequence, "", fixed = TRUE),
      function(x) sum(x == "N"),
      integer(1)
    )
  )
  n_bases / sum(
    nchar(records$Sequence, type = "bytes")
  )
}

file_summary$NBaseFraction <- c(
  count_n(forward),
  count_n(reverse),
  count_n(barcodes)
)

knitr::kable(
  transform(
    file_summary,
    MedianReadMeanPhred = round(
      MedianReadMeanPhred,
      2
    ),
    Q30BaseFraction = scales::percent(
      Q30BaseFraction,
      accuracy = 0.1
    ),
    NBaseFraction = scales::percent(
      NBaseFraction,
      accuracy = 0.01
    )
  ),
  caption = "Read-level FASTQ audit"
)

dir.create(
  "results/06-fastq",
  recursive = TRUE,
  showWarnings = FALSE
)

file_audit <- file_summary |>
  left_join(
    known_sha256[, c(
      "File", "ObservedSHA256",
      "SourceSet"
    )],
    by = c("Read" = "File")
  )

utils::write.table(
  file_audit,
  "results/06-fastq/fastq-file-audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  cycle_summary,
  "results/06-fastq/quality-by-cycle.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  record_example,
  "results/06-fastq/fastq-record-example.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  metadata_summary,
  "results/06-fastq/metadata-summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_architecture <- ggplot() +
  annotate(
    "rect",
    xmin = 0, xmax = 12,
    ymin = 0.92, ymax = 1.18,
    fill = pal_pub[5]
  ) +
  annotate(
    "rect",
    xmin = 12, xmax = 88,
    ymin = 0.92, ymax = 1.18,
    fill = "#E6E6E6"
  ) +
  annotate(
    "rect",
    xmin = 88, xmax = 100,
    ymin = 0.92, ymax = 1.18,
    fill = pal_pub[4]
  ) +
  annotate(
    "segment",
    x = 2, xend = 63,
    y = 0.68, yend = 0.68,
    linewidth = 1.2,
    color = pal_pub[1],
    arrow = grid::arrow(
      length = grid::unit(2.5, "mm"),
      type = "closed"
    )
  ) +
  annotate(
    "segment",
    x = 98, xend = 37,
    y = 0.40, yend = 0.40,
    linewidth = 1.2,
    color = pal_pub[2],
    arrow = grid::arrow(
      length = grid::unit(2.5, "mm"),
      type = "closed"
    )
  ) +
  annotate(
    "rect",
    xmin = 37, xmax = 63,
    ymin = 0.28, ymax = 0.80,
    fill = scales::alpha(pal_pub[3], 0.16),
    color = pal_pub[3],
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 6, y = 1.30,
    label = "515F primer",
    size = 3.2
  ) +
  annotate(
    "text",
    x = 50, y = 1.30,
    label = "16S V4 insert",
    size = 3.2
  ) +
  annotate(
    "text",
    x = 94, y = 1.30,
    label = "806R primer",
    size = 3.2
  ) +
  annotate(
    "text",
    x = 18, y = 0.76,
    label = "Read 1",
    color = pal_pub[1],
    fontface = "bold",
    size = 3.2
  ) +
  annotate(
    "text",
    x = 82, y = 0.48,
    label = "Read 2",
    color = pal_pub[2],
    fontface = "bold",
    size = 3.2
  ) +
  annotate(
    "text",
    x = 50, y = 0.18,
    label = "Required overlap",
    color = pal_pub[3],
    size = 3
  ) +
  annotate(
    "rect",
    xmin = 31, xmax = 69,
    ymin = -0.12, ymax = 0.03,
    fill = "#F7F7F7",
    color = "#666666",
    linewidth = 0.4
  ) +
  annotate(
    "text",
    x = 50, y = -0.045,
    label = "Index read assigns the molecule to a sample",
    size = 2.8
  ) +
  coord_cartesian(
    xlim = c(-2, 102),
    ylim = c(-0.18, 1.45),
    clip = "off"
  ) +
  labs(
    title = "Paired-end amplicon architecture",
    subtitle = "Primer removal and overlap are linked decisions"
  ) +
  theme_void(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      color = "#4D4D4D",
      margin = margin(b = 7)
    ),
    plot.margin = margin(8, 10, 6, 10)
  )

save_pub(
  p_architecture,
  "figures/06-read-architecture",
  width = 155,
  height = 88
)

p_architecture

anatomy <- record_example |>
  mutate(
    Display = case_when(
      Field == "Header" ~ substr(Value, 1, 53),
      Field == "Sequence" ~ paste0(
        substr(Value, 1, 48),
        "…  [151 nt]"
      ),
      Field == "Quality" ~ paste0(
        substr(Value, 1, 48),
        "…  [151 chars]"
      ),
      TRUE ~ Value
    ),
    Meaning = c(
      "Read identifier",
      "Nucleotide sequence",
      "Record separator",
      "Phred+33 quality"
    ),
    PlotRow = rev(Line)
  )

p_anatomy <- ggplot(anatomy) +
  geom_tile(
    aes(
      x = 0,
      y = PlotRow,
      fill = Field
    ),
    width = 0.48,
    height = 0.76,
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      x = 0,
      y = PlotRow,
      label = paste0("Line ", Line)
    ),
    color = "white",
    fontface = "bold",
    size = 3.1
  ) +
  geom_text(
    aes(
      x = 0.38,
      y = PlotRow,
      label = Display
    ),
    family = "mono",
    hjust = 0,
    size = 2.75,
    color = "#222222"
  ) +
  geom_text(
    aes(
      x = 4.20,
      y = PlotRow,
      label = Meaning
    ),
    hjust = 0,
    size = 2.9,
    color = "#4D4D4D"
  ) +
  scale_fill_manual(
    values = setNames(
      pal_pub[c(1, 3, 8, 2)],
      c(
        "Header", "Sequence",
        "Separator", "Quality"
      )
    ),
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(-0.28, 5.55),
    ylim = c(0.45, 4.55),
    clip = "off"
  ) +
  labs(
    title = "One FASTQ record has exactly four lines",
    subtitle = "Real Atacama forward read; long sequence and quality strings are abbreviated"
  ) +
  theme_void(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(
      color = "#4D4D4D",
      margin = margin(b = 8)
    ),
    plot.margin = margin(7, 8, 7, 8)
  )

save_pub(
  p_anatomy,
  "figures/06-fastq-anatomy",
  width = 180,
  height = 90
)

p_anatomy

p_quality <- ggplot(
  cycle_summary,
  aes(
    x = Cycle,
    y = MedianQ,
    color = Read,
    fill = Read
  )
) +
  geom_ribbon(
    aes(
      ymin = Q10,
      ymax = Q90
    ),
    alpha = 0.16,
    color = NA
  ) +
  geom_line(linewidth = 0.8) +
  geom_hline(
    yintercept = 30,
    linetype = 2,
    linewidth = 0.45,
    color = "#444444"
  ) +
  annotate(
    "text",
    x = 149,
    y = 30.8,
    label = "Q30",
    hjust = 1,
    vjust = 0,
    size = 3,
    color = "#444444"
  ) +
  scale_color_manual(
    values = c(
      Forward = pal_pub[1],
      Reverse = pal_pub[2]
    )
  ) +
  scale_fill_manual(
    values = c(
      Forward = pal_pub[1],
      Reverse = pal_pub[2]
    )
  ) +
  scale_x_continuous(
    breaks = c(1, 25, 50, 75, 100, 125, 151),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 42),
    breaks = seq(0, 40, 10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Sequencing cycle",
    y = "Phred quality",
    color = NULL,
    fill = NULL,
    title = "Read 2 quality declines near the tail",
    subtitle = "Median and 10th–90th percentiles; n = 2,000 real read pairs"
  ) +
  theme_pub(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(face = "bold")
  )

save_pub(
  p_quality,
  "figures/06-quality-profile",
  width = 135,
  height = 95
)

p_quality

flow <- tibble(
  Step = 1:4,
  Stage = factor(
    c(
      "Raw", "Demultiplexed",
      "Denoised", "Analysis-ready"
    ),
    levels = c(
      "Raw", "Demultiplexed",
      "Denoised", "Analysis-ready"
    )
  ),
  Label = c(
    "Multiplexed\nFASTQ triplet",
    "Per-sample\nR1 + R2",
    "ASV table +\nsequences",
    "otutab +\ntaxonomy +\nmetadata"
  ),
  Detail = c(
    "Forward, reverse, barcode\nshare record order",
    "Manifest maps files\nto sample IDs",
    "Trim, denoise, merge,\nremove chimeras",
    "Stable input contract\nfor R analysis"
  )
)

p_contract <- ggplot(flow) +
  geom_segment(
    data = tibble(
      x = 1:3 + 0.42,
      xend = 2:4 - 0.42,
      y = 1,
      yend = 1
    ),
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    linewidth = 0.7,
    color = "#666666",
    arrow = grid::arrow(
      type = "closed",
      length = grid::unit(2.2, "mm")
    )
  ) +
  geom_rect(
    aes(
      xmin = Step - 0.40,
      xmax = Step + 0.40,
      ymin = 0.82,
      ymax = 1.18,
      fill = Stage
    ),
    color = "white",
    linewidth = 0.7
  ) +
  geom_text(
    aes(
      x = Step,
      y = 1,
      label = Label
    ),
    color = "white",
    fontface = "bold",
    lineheight = 0.95,
    size = 2.8
  ) +
  geom_text(
    aes(
      x = Step,
      y = 0.58,
      label = Detail
    ),
    size = 2.6,
    color = "#4D4D4D",
    lineheight = 0.95
  ) +
  scale_fill_manual(
    values = setNames(
      pal_pub[c(2, 5, 3, 1)],
      levels(flow$Stage)
    ),
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(0.48, 4.52),
    ylim = c(0.34, 1.30),
    clip = "off"
  ) +
  labs(
    title = "From raw reads to the standard triad",
    subtitle = "Metadata anchors sample identity across every transition"
  ) +
  theme_void(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(
      color = "#4D4D4D",
      margin = margin(b = 9)
    ),
    plot.margin = margin(8, 10, 5, 10)
  )

save_pub(
  p_contract,
  "figures/06-file-contract",
  width = 180,
  height = 92
)

p_contract
