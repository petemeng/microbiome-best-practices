# 出版级作图地基：ggplot2 主题、色盲友好配色、字体、PDF/SVG 矢量 + TIFF/PNG 300–600 ppi 导出
# Run sequentially in a new working directory.
# Required packages: colorspace, ggplot2, jsonlite, ragg, readr, scales, svglite, systemfonts.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "scripts/validate_publication_graphics.R",
  "data/small/metadata.tsv",
  "data/small/otutab.tsv",
  "data/small/source_summary.json",
  "data/small/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

set.seed(20260721)

required_packages <- c(
  "colorspace", "digest", "dplyr", "ggplot2", "jsonlite", "ragg",
  "readr", "scales", "svglite", "systemfonts", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
stopifnot(length(missing_packages) == 0L)

font_pub <- "DejaVu Sans"
font_info <- systemfonts::font_info(font_pub)
stopifnot(
  nrow(font_info) >= 1L,
  identical(font_info$family[[1L]], font_pub),
  file.exists(font_info$path[[1L]]),
  isTRUE(font_info$scalable[[1L]])
)

pal_pub <- c(
  blue = "#0072B2",
  orange = "#E69F00",
  green = "#009E73",
  vermillion = "#D55E00",
  purple = "#CC79A7",
  sky = "#56B4E9",
  yellow = "#F0E442",
  grey = "#6B7280"
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
      plot.caption = ggplot2::element_text(
        colour = "#666666", hjust = 0
      ),
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2", colour = "#B3B3B3", linewidth = 0.3
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.position = "top"
    )
}

save_pub <- function(
  plot,
  file_base,
  width = 89,
  height = 70,
  units = "mm",
  dpi = 600,
  write_svg = TRUE,
  write_tiff = TRUE,
  base_family = font_pub
) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    paste0(file_base, ".pdf"), plot,
    width = width, height = height, units = units,
    device = grDevices::cairo_pdf,
    family = base_family,
    bg = "white"
  )

  if (isTRUE(write_svg)) {
    ggplot2::ggsave(
      paste0(file_base, ".svg"), plot,
      width = width, height = height, units = units,
      device = svglite::svglite,
      bg = "white"
    )
  }

  ggplot2::ggsave(
    paste0(file_base, ".png"), plot,
    width = width, height = height, units = units,
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white"
  )

  if (isTRUE(write_tiff)) {
    ggplot2::ggsave(
      paste0(file_base, ".tiff"), plot,
      width = width, height = height, units = units,
      dpi = dpi,
      device = ragg::agg_tiff,
      compression = "lzw",
      bg = "white"
    )
  }

  invisible(plot)
}

input_paths <- c(
  otutab = "data/small/otutab.tsv",
  taxonomy = "data/small/taxonomy.tsv",
  metadata = "data/small/metadata.tsv",
  source_summary = "data/small/source_summary.json"
)
expected_sha256 <- c(
  otutab = "76fa79c38da889f35978dc86da4641a270961746709ff38049ee5f67e3c6f7a3",
  taxonomy = "725280bb9a0cd9bda7b540022e92af945ceed52527f8b2d220b055b4489f6901",
  metadata = "df24771dccf27607ddf922c6bca2cafa876d946fbe2e09d14b601accce66ba64",
  source_summary = "acc15b18d3f85d6d35770d0db7580d91d0a55a838862876536500a8d7c75711b"
)
observed_sha256 <- vapply(
  input_paths,
  digest::digest,
  character(1),
  algo = "sha256",
  serialize = FALSE,
  file = TRUE
)
stopifnot(identical(observed_sha256, expected_sha256))

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    na = character()
  )
  ids <- as.character(x[[1L]])
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

otutab <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])

otu_matrix <- as.matrix(otutab)
storage.mode(otu_matrix) <- "numeric"

stopifnot(
  identical(rownames(otu_matrix), rownames(taxonomy)),
  identical(colnames(otu_matrix), rownames(metadata)),
  all(is.finite(otu_matrix)),
  all(otu_matrix >= 0),
  all(otu_matrix == round(otu_matrix)),
  all(c("Group", "Type", "Saline") %in% colnames(metadata)),
  nrow(otu_matrix) == 13628L,
  ncol(otu_matrix) == 90L,
  sum(otu_matrix) == 1619670
)

dim(otu_matrix)
otu_matrix[1:5, 1:5]
taxonomy[1:5, , drop = FALSE]
head(metadata)

group_labels <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan Plateau"
)
group_palette <- c(
  "Inland wetland" = "#0072B2",
  "Coastal wetland" = "#D55E00",
  "Tibetan Plateau" = "#009E73"
)
group_shapes <- c(
  "Inland wetland" = 21,
  "Coastal wetland" = 22,
  "Tibetan Plateau" = 24
)

sample_plot_data <- data.frame(
  SampleID = colnames(otu_matrix),
  LibrarySize = as.numeric(colSums(otu_matrix)),
  ObservedFeatures = as.numeric(colSums(otu_matrix > 0)),
  Group = metadata[colnames(otu_matrix), "Group"],
  GroupLabel = unname(
    group_labels[metadata[colnames(otu_matrix), "Group"]]
  ),
  Salinity = metadata[colnames(otu_matrix), "Saline"],
  stringsAsFactors = FALSE
)
sample_plot_data$GroupLabel <- factor(
  sample_plot_data$GroupLabel,
  levels = names(group_palette)
)

palette_views <- list(
  Standard = unname(group_palette),
  "Deuteranopia simulation" = colorspace::deutan(
    unname(group_palette), severity = 1
  ),
  "Protanopia simulation" = colorspace::protan(
    unname(group_palette), severity = 1
  )
)

palette_audit <- do.call(
  rbind,
  lapply(names(palette_views), function(view) {
    data.frame(
      View = view,
      Group = names(group_palette),
      Hex = unname(palette_views[[view]]),
      Shape = unname(group_shapes[names(group_palette)]),
      ContrastWhite = as.numeric(
        colorspace::contrast_ratio(palette_views[[view]], "#FFFFFF")
      ),
      ContrastBlack = as.numeric(
        colorspace::contrast_ratio(palette_views[[view]], "#000000")
      )
    )
  })
)

encoding_plot_data <- do.call(
  rbind,
  lapply(names(palette_views), function(view) {
    out <- sample_plot_data
    out$Vision <- view
    view_palette <- stats::setNames(
      palette_views[[view]],
      names(group_palette)
    )
    out$DisplayColor <- unname(
      view_palette[as.character(out$GroupLabel)]
    )
    out
  })
)
encoding_plot_data$Vision <- factor(
  encoding_plot_data$Vision,
  levels = names(palette_views)
)

palette_audit

encoding_plot <- ggplot2::ggplot(
  encoding_plot_data,
  ggplot2::aes(x = ObservedFeatures, y = LibrarySize)
) +
  ggplot2::geom_point(
    ggplot2::aes(fill = DisplayColor, shape = GroupLabel),
    colour = "#1A1A1A",
    stroke = 0.35,
    size = 2.2,
    alpha = 0.85
  ) +
  ggplot2::facet_wrap(~Vision, nrow = 1L) +
  ggplot2::scale_fill_identity(
    name = "Wetland group",
    breaks = unname(group_palette),
    labels = names(group_palette),
    guide = "legend"
  ) +
  ggplot2::scale_shape_manual(
    name = "Wetland group",
    values = group_shapes
  ) +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Colour is not the only group cue",
    subtitle = paste(
      "The same 90 samples retain distinct shapes under",
      "two colour-vision simulations"
    ),
    x = "Observed features",
    y = "Library size (reads)"
  ) +
  theme_pub(base_size = 9.5) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.spacing.x = grid::unit(3, "mm")
  )

save_pub(
  encoding_plot,
  "figures/18-encoding-redundancy-audit",
  width = 183,
  height = 95,
  dpi = 600
)

phylum <- sub(
  "^p__",
  "",
  as.character(taxonomy[rownames(otu_matrix), "Phylum"])
)
phylum[is.na(phylum) | !nzchar(phylum)] <- "Unassigned"

phylum_totals <- sort(
  rowsum(rowSums(otu_matrix), phylum)[, 1L],
  decreasing = TRUE
)
top_phyla <- names(phylum_totals)[seq_len(5L)]
phylum_class <- ifelse(phylum %in% top_phyla, phylum, "Other")
phylum_counts <- rowsum(
  otu_matrix,
  group = phylum_class,
  reorder = FALSE
)
phylum_relative <- sweep(
  phylum_counts,
  2L,
  colSums(phylum_counts),
  "/"
)

phylum_composition <- do.call(
  rbind,
  lapply(names(group_labels), function(group_code) {
    selected <- rownames(metadata)[metadata$Group == group_code]
    data.frame(
      Group = unname(group_labels[[group_code]]),
      Phylum = rownames(phylum_relative),
      MeanRelativeAbundance = rowMeans(
        phylum_relative[, selected, drop = FALSE]
      ),
      SampleN = length(selected),
      stringsAsFactors = FALSE
    )
  })
)

phylum_order <- c(top_phyla, "Other")
phylum_composition$Group <- factor(
  phylum_composition$Group,
  levels = names(group_palette)
)
phylum_composition$Phylum <- factor(
  phylum_composition$Phylum,
  levels = phylum_order
)
stopifnot(
  all(
    abs(
      tapply(
        phylum_composition$MeanRelativeAbundance,
        phylum_composition$Group,
        sum
      ) - 1
    ) < 1e-12
  )
)

phylum_composition

composition_palette <- stats::setNames(
  c(
    pal_pub[["blue"]],
    pal_pub[["orange"]],
    pal_pub[["green"]],
    pal_pub[["vermillion"]],
    pal_pub[["purple"]],
    pal_pub[["grey"]]
  ),
  phylum_order
)

composition_plot <- ggplot2::ggplot(
  phylum_composition,
  ggplot2::aes(
    x = Group,
    y = MeanRelativeAbundance,
    fill = Phylum
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    colour = "white",
    linewidth = 0.25
  ) +
  ggplot2::scale_fill_manual(
    values = composition_palette,
    breaks = phylum_order,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Community composition at final size",
    subtitle = "Mean relative abundance; n = 30 per group",
    x = NULL,
    y = "Mean relative abundance",
    fill = NULL
  ) +
  theme_pub(base_size = 8.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 8),
    legend.key.height = grid::unit(3.3, "mm"),
    panel.grid.major.x = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE)
  )

save_pub(
  composition_plot,
  "figures/18-final-size-community-summary",
  width = 89,
  height = 105,
  dpi = 600
)

decision_boxes <- data.frame(
  xmin = c(1.35, 0.15, 2.55, 0.15, 2.55, 0.65),
  xmax = c(2.65, 1.45, 3.85, 1.45, 3.85, 3.35),
  ymin = c(3.85, 2.25, 2.25, 0.65, 0.65, -0.35),
  ymax = c(4.65, 3.15, 3.15, 1.55, 1.55, 0.25),
  fill = c(
    "#E8F1F8", "#EAF5F1", "#FFF3E6",
    "#DDECF6", "#FCE6DD", "#F2F2F2"
  ),
  label = c(
    "Finished ggplot\nat final dimensions",
    "Lines, points,\nand text",
    "Pixel-based image\nor journal request",
    "PDF / SVG\nvector geometry",
    "TIFF / PNG\nexplicit ppi",
    "Open the exported file; audit width, height, fonts, and labels"
  ),
  stringsAsFactors = FALSE
)
decision_arrows <- data.frame(
  x = c(1.75, 2.25, 0.8, 3.2, 0.8, 3.2),
  xend = c(0.8, 3.2, 0.8, 3.2, 1.55, 2.45),
  y = c(3.85, 3.85, 2.25, 2.25, 0.65, 0.65),
  yend = c(3.15, 3.15, 1.55, 1.55, 0.25, 0.25)
)

decision_plot <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = decision_arrows,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.55,
    colour = "#4D4D4D",
    arrow = grid::arrow(
      length = grid::unit(2.1, "mm"),
      type = "closed"
    )
  ) +
  ggplot2::geom_rect(
    data = decision_boxes,
    ggplot2::aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill
    ),
    colour = "#4D4D4D",
    linewidth = 0.45
  ) +
  ggplot2::geom_text(
    data = decision_boxes,
    ggplot2::aes(
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2,
      label = label
    ),
    family = font_pub,
    colour = "#1A1A1A",
    lineheight = 0.98,
    size = c(3.2, 3, 3, 3.2, 3.2, 2.8),
    fontface = c("bold", "plain", "plain", "bold", "bold", "plain")
  ) +
  ggplot2::scale_fill_identity() +
  ggplot2::coord_cartesian(
    xlim = c(0, 4),
    ylim = c(-0.45, 4.8),
    clip = "off"
  ) +
  ggplot2::labs(
    title = "Choose the file format from the graphic content",
    subtitle = paste(
      "Physical dimensions come first; ppi is a raster-only",
      "sampling contract"
    )
  ) +
  theme_pub(base_size = 10) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "none",
    plot.margin = ggplot2::margin(5, 7, 5, 7)
  )

save_pub(
  decision_plot,
  "figures/18-export-decision-map",
  width = 183,
  height = 90,
  dpi = 600
)

probe_plot <- ggplot2::ggplot(
  sample_plot_data,
  ggplot2::aes(x = ObservedFeatures, y = LibrarySize)
) +
  ggplot2::geom_point(
    ggplot2::aes(fill = GroupLabel, shape = GroupLabel),
    colour = "#1A1A1A",
    stroke = 0.3,
    size = 1.7,
    alpha = 0.85
  ) +
  ggplot2::scale_fill_manual(values = group_palette) +
  ggplot2::scale_shape_manual(values = group_shapes) +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Resolution probe",
    x = "Observed features",
    y = "Library size (reads)"
  ) +
  theme_pub(base_size = 8) +
  ggplot2::theme(legend.position = "none")

probe_ppi <- c(72L, 300L, 600L)
probe_paths <- file.path(
  "results/18-publication-graphics",
  paste0("resolution-probe-", probe_ppi, "ppi.png")
)
dir.create(
  "results/18-publication-graphics",
  recursive = TRUE,
  showWarnings = FALSE
)

for (i in seq_along(probe_ppi)) {
  ggplot2::ggsave(
    probe_paths[[i]],
    plot = probe_plot,
    width = 89,
    height = 70,
    units = "mm",
    dpi = probe_ppi[[i]],
    device = ragg::agg_png,
    bg = "white"
  )
}

resolution_plot_data <- data.frame(
  target_ppi = probe_ppi,
  pixel_width = floor(89 / 25.4 * probe_ppi),
  pixel_height = floor(70 / 25.4 * probe_ppi)
)
resolution_plot_data$megapixels <- with(
  resolution_plot_data,
  pixel_width * pixel_height / 1e6
)
resolution_plot_data$PPI <- factor(
  paste0(resolution_plot_data$target_ppi, " ppi"),
  levels = paste0(probe_ppi, " ppi")
)
resolution_plot_data$PixelLabel <- paste0(
  scales::comma(resolution_plot_data$pixel_width),
  " x ",
  scales::comma(resolution_plot_data$pixel_height),
  " px\n",
  sprintf("%.2f MP", resolution_plot_data$megapixels)
)

resolution_plot <- ggplot2::ggplot(
  resolution_plot_data,
  ggplot2::aes(x = PPI, y = target_ppi, fill = PPI)
) +
  ggplot2::geom_col(
    width = 0.62,
    colour = "#1A1A1A",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = PixelLabel),
    vjust = -0.35,
    family = font_pub,
    size = 3.2,
    lineheight = 0.95
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      pal_pub[["sky"]],
      pal_pub[["blue"]],
      pal_pub[["vermillion"]]
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 700),
    breaks = c(0, 72, 300, 600),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "Physical size stays fixed; pixel budget does not",
    subtitle = "The same plot was exported at 89 x 70 mm",
    x = NULL,
    y = "Pixels per inch (ppi)",
    caption = paste(
      "Increasing metadata cannot restore detail that was absent",
      "from a low-resolution source."
    )
  ) +
  theme_pub(base_size = 10) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major.x = ggplot2::element_blank()
  )

save_pub(
  resolution_plot,
  "figures/18-raster-resolution-audit",
  width = 183,
  height = 92,
  dpi = 600
)

run_analysis <- function(command, args) {
  output <- system2(command, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop(paste(tail(output, 20), collapse = "\n"))
  invisible(output)
}
run_analysis(file.path(R.home("bin"), "Rscript"), c("--vanilla", "scripts/validate_publication_graphics.R", "--project-root", ".", "--input-dir", "data/small", "--output-dir", "results/18-publication-graphics", "--figure-dir", "figures"))

format_audit <- readr::read_tsv(
  "results/18-publication-graphics/format-audit.tsv",
  show_col_types = FALSE
)
validation_summary <- jsonlite::read_json(
  "results/18-publication-graphics/publication-graphics-summary.json",
  simplifyVector = TRUE
)

stopifnot(
  nrow(format_audit) == 16L,
  all(format_audit$dimension_status == "PASS")
)

format_audit[, c(
  "figure",
  "extension",
  "bytes",
  "actual_width_mm",
  "actual_height_mm",
  "pixel_width",
  "pixel_height",
  "effective_ppi_x",
  "compression",
  "font_status",
  "dimension_status"
)]
