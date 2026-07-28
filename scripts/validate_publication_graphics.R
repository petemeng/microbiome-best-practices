#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "Asia/Shanghai")
set.seed(20260721)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("project-root", "input-dir", "output-dir", "figure-dir")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0L) {
  stop(
    "Missing required arguments: ",
    paste(paste0("--", missing_args), collapse = ", "),
    call. = FALSE
  )
}

project_root <- normalizePath(args[["project-root"]], mustWork = TRUE)
input_dir <- normalizePath(args[["input-dir"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
figure_dir <- normalizePath(args[["figure-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "colorspace", "digest", "dplyr", "ggplot2", "jsonlite", "ragg",
  "readr", "scales", "svglite", "systemfonts", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

sanitize_text <- function(x) {
  x <- gsub(project_root, "<PROJECT_ROOT>", x, fixed = TRUE)
  x <- gsub(input_dir, "<INPUT_DIR>", x, fixed = TRUE)
  x <- gsub(output_dir, "<OUTPUT_DIR>", x, fixed = TRUE)
  x <- gsub(figure_dir, "<FIGURE_DIR>", x, fixed = TRUE)
  x <- gsub(path.expand("~"), "<HOME>", x, fixed = TRUE)
  x
}

write_tsv <- function(x, filename) {
  readr::write_tsv(x, file.path(output_dir, filename), na = "")
}

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    na = character()
  )
  ids <- as.character(x[[1L]])
  if (anyDuplicated(ids)) {
    stop("Duplicated identifiers in ", basename(path), call. = FALSE)
  }
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

check_rows <- list()
add_check <- function(id, category, observed, expected, passed) {
  check_rows[[length(check_rows) + 1L]] <<- data.frame(
    check_id = id,
    category = category,
    observed = paste(observed, collapse = ","),
    expected = paste(expected, collapse = ","),
    status = ifelse(isTRUE(passed), "PASS", "FAIL"),
    check.names = FALSE
  )
}

command_path <- function(command) {
  path <- Sys.which(command)
  if (!nzchar(path)) "" else normalizePath(path, mustWork = TRUE)
}

run_command <- function(command, args = character()) {
  status <- 0L
  output <- tryCatch(
    system2(
      command,
      args = args,
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      status <<- 1L
      conditionMessage(e)
    }
  )
  command_status <- attr(output, "status")
  if (!is.null(command_status)) {
    status <- as.integer(command_status)
  }
  list(status = status, output = output)
}

input_paths <- c(
  otutab = file.path(input_dir, "otutab.tsv"),
  taxonomy = file.path(input_dir, "taxonomy.tsv"),
  metadata = file.path(input_dir, "metadata.tsv"),
  source_summary = file.path(input_dir, "source_summary.json")
)
expected_sha256 <- c(
  otutab = "76fa79c38da889f35978dc86da4641a270961746709ff38049ee5f67e3c6f7a3",
  taxonomy = "725280bb9a0cd9bda7b540022e92af945ceed52527f8b2d220b055b4489f6901",
  metadata = "df24771dccf27607ddf922c6bca2cafa876d946fbe2e09d14b601accce66ba64",
  source_summary = "acc15b18d3f85d6d35770d0db7580d91d0a55a838862876536500a8d7c75711b"
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Article 18 input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

observed_sha256 <- vapply(input_paths, sha256_file, character(1))
input_audit <- data.frame(
  asset = names(input_paths),
  file = basename(input_paths),
  bytes = unname(file.info(input_paths)$size),
  observed_sha256 = unname(observed_sha256),
  expected_sha256 = unname(expected_sha256[names(input_paths)]),
  status = ifelse(
    observed_sha256 == expected_sha256[names(input_paths)],
    "PASS",
    "FAIL"
  ),
  check.names = FALSE
)
for (i in seq_len(nrow(input_audit))) {
  add_check(
    paste0("sha256-", input_audit$asset[[i]]),
    "input",
    input_audit$observed_sha256[[i]],
    input_audit$expected_sha256[[i]],
    input_audit$status[[i]] == "PASS"
  )
}

expected_versions <- c(
  colorspace = "2.1.0",
  digest = "0.6.36",
  dplyr = "1.1.4",
  ggplot2 = "3.5.2",
  jsonlite = "1.8.8",
  ragg = "1.3.2",
  readr = "2.1.5",
  scales = "1.3.0",
  svglite = "2.1.3",
  systemfonts = "1.1.0",
  tibble = "3.2.1",
  tidyr = "1.3.1"
)
observed_versions <- vapply(
  names(expected_versions),
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)
package_audit <- data.frame(
  package = names(expected_versions),
  expected_version = unname(expected_versions),
  observed_version = unname(observed_versions),
  status = ifelse(observed_versions == expected_versions, "PASS", "FAIL"),
  check.names = FALSE
)
add_check(
  "r-version",
  "environment",
  paste(R.version$major, R.version$minor, sep = "."),
  "4.4.1",
  identical(paste(R.version$major, R.version$minor, sep = "."), "4.4.1")
)
for (i in seq_len(nrow(package_audit))) {
  add_check(
    paste0("package-", package_audit$package[[i]]),
    "environment",
    package_audit$observed_version[[i]],
    package_audit$expected_version[[i]],
    package_audit$status[[i]] == "PASS"
  )
}

required_commands <- c("pdfinfo", "pdffonts", "identify")
command_audit <- data.frame(
  command = required_commands,
  path = vapply(required_commands, command_path, character(1)),
  stringsAsFactors = FALSE
)
command_audit$status <- ifelse(nzchar(command_audit$path), "PASS", "FAIL")
for (i in seq_len(nrow(command_audit))) {
  add_check(
    paste0("command-", command_audit$command[[i]]),
    "environment",
    command_audit$path[[i]],
    "available",
    command_audit$status[[i]] == "PASS"
  )
}
if (any(command_audit$status != "PASS")) {
  write_tsv(input_audit, "input-audit.tsv")
  write_tsv(package_audit, "package-audit.tsv")
  write_tsv(command_audit, "command-audit.tsv")
  stop("Required graphics-audit command is unavailable.", call. = FALSE)
}

otutab <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])
source_summary <- jsonlite::read_json(
  input_paths[["source_summary"]],
  simplifyVector = TRUE
)

otu_matrix <- as.matrix(otutab)
storage.mode(otu_matrix) <- "numeric"
feature_ids <- rownames(otu_matrix)
sample_ids <- colnames(otu_matrix)

add_check("feature-count", "data", nrow(otu_matrix), 13628, nrow(otu_matrix) == 13628)
add_check("sample-count", "data", ncol(otu_matrix), 90, ncol(otu_matrix) == 90)
add_check("read-count", "data", sum(otu_matrix), 1619670, sum(otu_matrix) == 1619670)
add_check("taxonomy-feature-count", "data", nrow(taxonomy), 13628, nrow(taxonomy) == 13628)
add_check("metadata-sample-count", "data", nrow(metadata), 90, nrow(metadata) == 90)
add_check("feature-ids-unique", "data", anyDuplicated(feature_ids), 0, !anyDuplicated(feature_ids))
add_check("sample-ids-unique", "data", anyDuplicated(sample_ids), 0, !anyDuplicated(sample_ids))
add_check(
  "taxonomy-feature-ids",
  "data",
  length(intersect(feature_ids, rownames(taxonomy))),
  13628,
  setequal(feature_ids, rownames(taxonomy))
)
add_check(
  "metadata-sample-ids",
  "data",
  length(intersect(sample_ids, rownames(metadata))),
  90,
  setequal(sample_ids, rownames(metadata))
)
add_check("counts-finite", "data", all(is.finite(otu_matrix)), TRUE, all(is.finite(otu_matrix)))
add_check("counts-nonnegative", "data", min(otu_matrix), ">=0", min(otu_matrix) >= 0)
add_check(
  "counts-integer",
  "data",
  max(abs(otu_matrix - round(otu_matrix))),
  0,
  all(otu_matrix == round(otu_matrix))
)
add_check(
  "source-summary-shape",
  "data",
  c(source_summary$features, source_summary$samples),
  c(13628, 90),
  identical(
    as.integer(c(source_summary$features, source_summary$samples)),
    c(13628L, 90L)
  )
)
add_check(
  "metadata-columns",
  "data",
  intersect(c("Group", "Type", "Saline"), colnames(metadata)),
  c("Group", "Type", "Saline"),
  all(c("Group", "Type", "Saline") %in% colnames(metadata))
)
group_counts <- table(metadata$Group)
add_check(
  "group-levels",
  "data",
  sort(names(group_counts)),
  c("CW", "IW", "TW"),
  identical(sort(names(group_counts)), c("CW", "IW", "TW"))
)
add_check(
  "group-balance",
  "data",
  as.integer(group_counts[c("CW", "IW", "TW")]),
  c(30, 30, 30),
  identical(as.integer(group_counts[c("CW", "IW", "TW")]), c(30L, 30L, 30L))
)

font_family <- "DejaVu Sans"
font_info <- systemfonts::font_info(font_family)
font_path <- font_info$path[[1L]]
font_glyphs <- systemfonts::glyph_info(
  "Microbiome Wetland Group 0123456789",
  family = font_family
)
font_audit <- data.frame(
  requested_family = font_family,
  resolved_family = font_info$family[[1L]],
  style = font_info$style[[1L]],
  path = sanitize_text(font_path),
  file_sha256 = if (file.exists(font_path)) sha256_file(font_path) else "",
  scalable = font_info$scalable[[1L]],
  glyphs_checked = nrow(font_glyphs),
  missing_glyphs = sum(font_glyphs$index == 0L),
  stringsAsFactors = FALSE
)
add_check(
  "font-resolved-family",
  "font",
  font_info$family[[1L]],
  font_family,
  identical(font_info$family[[1L]], font_family)
)
add_check("font-file-exists", "font", file.exists(font_path), TRUE, file.exists(font_path))
add_check("font-scalable", "font", font_info$scalable[[1L]], TRUE, isTRUE(font_info$scalable[[1L]]))
add_check(
  "font-glyph-coverage",
  "font",
  sum(font_glyphs$index == 0L),
  0,
  all(font_glyphs$index > 0L)
)

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
      view = view,
      group = names(group_palette),
      hex = unname(palette_views[[view]]),
      shape = unname(group_shapes[names(group_palette)]),
      contrast_vs_white = as.numeric(
        colorspace::contrast_ratio(palette_views[[view]], "#FFFFFF")
      ),
      contrast_vs_black = as.numeric(
        colorspace::contrast_ratio(palette_views[[view]], "#000000")
      ),
      stringsAsFactors = FALSE
    )
  })
)
add_check(
  "palette-hex-valid",
  "encoding",
  sum(grepl("^#[0-9A-F]{6}$", palette_audit$hex)),
  nrow(palette_audit),
  all(grepl("^#[0-9A-F]{6}$", palette_audit$hex))
)
add_check(
  "palette-view-count",
  "encoding",
  length(unique(palette_audit$view)),
  3,
  length(unique(palette_audit$view)) == 3L
)
add_check(
  "palette-colors-distinct-per-view",
  "encoding",
  min(vapply(split(palette_audit$hex, palette_audit$view), function(x) length(unique(x)), integer(1))),
  3,
  all(vapply(split(palette_audit$hex, palette_audit$view), function(x) length(unique(x)), integer(1)) == 3L)
)
add_check(
  "group-shapes-distinct",
  "encoding",
  length(unique(group_shapes)),
  3,
  length(unique(group_shapes)) == 3L
)

sample_plot_data <- data.frame(
  SampleID = sample_ids,
  LibrarySize = as.numeric(colSums(otu_matrix)),
  ObservedFeatures = as.numeric(colSums(otu_matrix > 0)),
  Group = metadata[sample_ids, "Group"],
  GroupLabel = unname(group_labels[metadata[sample_ids, "Group"]]),
  Salinity = metadata[sample_ids, "Saline"],
  stringsAsFactors = FALSE
)
sample_plot_data$GroupLabel <- factor(
  sample_plot_data$GroupLabel,
  levels = names(group_palette)
)
add_check(
  "sample-plot-rows",
  "derived-data",
  nrow(sample_plot_data),
  90,
  nrow(sample_plot_data) == 90L
)
add_check(
  "sample-plot-library-range",
  "derived-data",
  range(sample_plot_data$LibrarySize),
  c(10364, 37374),
  identical(as.numeric(range(sample_plot_data$LibrarySize)), c(10364, 37374))
)
add_check(
  "sample-plot-feature-range",
  "derived-data",
  range(sample_plot_data$ObservedFeatures),
  c(1888, 4480),
  identical(as.numeric(range(sample_plot_data$ObservedFeatures)), c(1888, 4480))
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
    out$DisplayColor <- unname(view_palette[as.character(out$GroupLabel)])
    out
  })
)
encoding_plot_data$Vision <- factor(
  encoding_plot_data$Vision,
  levels = names(palette_views)
)

phylum <- sub("^p__", "", as.character(taxonomy[feature_ids, "Phylum"]))
phylum[is.na(phylum) | !nzchar(phylum)] <- "Unassigned"
phylum_totals <- sort(rowsum(rowSums(otu_matrix), phylum)[, 1L], decreasing = TRUE)
top_phyla <- names(phylum_totals)[seq_len(5L)]
phylum_class <- ifelse(phylum %in% top_phyla, phylum, "Other")
phylum_counts <- rowsum(otu_matrix, group = phylum_class, reorder = FALSE)
phylum_relative <- sweep(phylum_counts, 2L, colSums(phylum_counts), "/")
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
composition_sums <- tapply(
  phylum_composition$MeanRelativeAbundance,
  phylum_composition$Group,
  sum
)
add_check(
  "composition-top-category-count",
  "derived-data",
  length(phylum_order),
  6,
  length(phylum_order) == 6L
)
add_check(
  "composition-group-sums",
  "derived-data",
  round(composition_sums, 12),
  1,
  all(abs(composition_sums - 1) < 1e-12)
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

theme_pub <- function(base_size = 10, base_family = font_family) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        colour = "#E6E6E6", linewidth = 0.25
      ),
      axis.text = ggplot2::element_text(colour = "#1A1A1A"),
      axis.title = ggplot2::element_text(colour = "#1A1A1A"),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.15)),
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

save_publication_plot <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(stem, ".pdf"),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = grDevices::cairo_pdf,
    family = font_family,
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".svg"),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = svglite::svglite,
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".png"),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".tiff"),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    device = ragg::agg_tiff,
    compression = "lzw",
    bg = "white"
  )
  invisible(plot)
}

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
    subtitle = "The same 90 samples retain distinct shapes under two colour-vision simulations",
    x = "Observed features",
    y = "Library size (reads)"
  ) +
  theme_pub(base_size = 9.5) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.spacing.x = grid::unit(3, "mm")
  )

composition_palette <- stats::setNames(
  c(
    pal_pub[["blue"]], pal_pub[["orange"]], pal_pub[["green"]],
    pal_pub[["vermillion"]], pal_pub[["purple"]], pal_pub[["grey"]]
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
  ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
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

decision_boxes <- data.frame(
  xmin = c(1.35, 0.15, 2.55, 0.15, 2.55, 0.65),
  xmax = c(2.65, 1.45, 3.85, 1.45, 3.85, 3.35),
  ymin = c(3.85, 2.25, 2.25, 0.65, 0.65, -0.35),
  ymax = c(4.65, 3.15, 3.15, 1.55, 1.55, 0.25),
  fill = c("#E8F1F8", "#EAF5F1", "#FFF3E6", "#DDECF6", "#FCE6DD", "#F2F2F2"),
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
    arrow = grid::arrow(length = grid::unit(2.1, "mm"), type = "closed")
  ) +
  ggplot2::geom_rect(
    data = decision_boxes,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    colour = "#4D4D4D",
    linewidth = 0.45
  ) +
  ggplot2::geom_text(
    data = decision_boxes,
    ggplot2::aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    family = font_family,
    colour = "#1A1A1A",
    lineheight = 0.98,
    size = c(3.2, 3, 3, 3.2, 3.2, 2.8),
    fontface = c("bold", "plain", "plain", "bold", "bold", "plain")
  ) +
  ggplot2::scale_fill_identity() +
  ggplot2::coord_cartesian(xlim = c(0, 4), ylim = c(-0.45, 4.8), clip = "off") +
  ggplot2::labs(
    title = "Choose the file format from the graphic content",
    subtitle = "Physical dimensions come first; ppi is a raster-only sampling contract"
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

figure_specification <- data.frame(
  stem = c(
    "18-encoding-redundancy-audit",
    "18-final-size-community-summary",
    "18-export-decision-map",
    "18-raster-resolution-audit"
  ),
  width_mm = c(183, 89, 183, 183),
  height_mm = c(95, 105, 90, 92),
  raster_ppi = c(600, 600, 600, 600),
  intended_layout = c("two-column", "single-column", "two-column", "two-column"),
  stringsAsFactors = FALSE
)

save_publication_plot(
  encoding_plot,
  file.path(figure_dir, "18-encoding-redundancy-audit"),
  width_mm = 183,
  height_mm = 95,
  dpi = 600
)
save_publication_plot(
  composition_plot,
  file.path(figure_dir, "18-final-size-community-summary"),
  width_mm = 89,
  height_mm = 105,
  dpi = 600
)
save_publication_plot(
  decision_plot,
  file.path(figure_dir, "18-export-decision-map"),
  width_mm = 183,
  height_mm = 90,
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
  output_dir,
  paste0("resolution-probe-", probe_ppi, "ppi.png")
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

identify_raster <- function(path) {
  result <- run_command(
    "identify",
    c(
      "-format",
      shQuote("%m|%w|%h|%x|%y|%U|%C"),
      shQuote(path)
    )
  )
  if (result$status != 0L || length(result$output) == 0L) {
    stop("identify failed for ", path, call. = FALSE)
  }
  fields <- strsplit(result$output[[1L]], "\\|", fixed = FALSE)[[1L]]
  if (length(fields) != 7L) {
    stop("Unexpected identify output for ", path, call. = FALSE)
  }
  data.frame(
    detected_format = fields[[1L]],
    pixel_width = as.integer(fields[[2L]]),
    pixel_height = as.integer(fields[[3L]]),
    density_x = as.numeric(fields[[4L]]),
    density_y = as.numeric(fields[[5L]]),
    density_unit = fields[[6L]],
    compression = fields[[7L]],
    stringsAsFactors = FALSE
  )
}

probe_audit <- do.call(
  rbind,
  lapply(seq_along(probe_paths), function(i) {
    info <- identify_raster(probe_paths[[i]])
    data.frame(
      file = basename(probe_paths[[i]]),
      target_ppi = probe_ppi[[i]],
      width_mm = 89,
      height_mm = 70,
      pixel_width = info$pixel_width,
      pixel_height = info$pixel_height,
      expected_pixel_width = floor(89 / 25.4 * probe_ppi[[i]]),
      expected_pixel_height = floor(70 / 25.4 * probe_ppi[[i]]),
      effective_ppi_x = info$pixel_width / (89 / 25.4),
      effective_ppi_y = info$pixel_height / (70 / 25.4),
      megapixels = info$pixel_width * info$pixel_height / 1e6,
      bytes = file.info(probe_paths[[i]])$size,
      stringsAsFactors = FALSE
    )
  })
)
for (i in seq_len(nrow(probe_audit))) {
  add_check(
    paste0("probe-", probe_audit$target_ppi[[i]], "-pixel-width"),
    "resolution",
    probe_audit$pixel_width[[i]],
    probe_audit$expected_pixel_width[[i]],
    probe_audit$pixel_width[[i]] == probe_audit$expected_pixel_width[[i]]
  )
  add_check(
    paste0("probe-", probe_audit$target_ppi[[i]], "-pixel-height"),
    "resolution",
    probe_audit$pixel_height[[i]],
    probe_audit$expected_pixel_height[[i]],
    probe_audit$pixel_height[[i]] == probe_audit$expected_pixel_height[[i]]
  )
  add_check(
    paste0("probe-", probe_audit$target_ppi[[i]], "-effective-ppi"),
    "resolution",
    round(c(probe_audit$effective_ppi_x[[i]], probe_audit$effective_ppi_y[[i]]), 3),
    paste0(probe_audit$target_ppi[[i]], " +/- 0.5"),
    all(
      abs(
        c(
          probe_audit$effective_ppi_x[[i]],
          probe_audit$effective_ppi_y[[i]]
        ) - probe_audit$target_ppi[[i]]
      ) <= 0.5
    )
  )
}

resolution_plot_data <- probe_audit
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
  ggplot2::geom_col(width = 0.62, colour = "#1A1A1A", linewidth = 0.35) +
  ggplot2::geom_text(
    ggplot2::aes(label = PixelLabel),
    vjust = -0.35,
    family = font_family,
    size = 3.2,
    lineheight = 0.95
  ) +
  ggplot2::scale_fill_manual(
    values = c(pal_pub[["sky"]], pal_pub[["blue"]], pal_pub[["vermillion"]])
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
    caption = "Increasing metadata cannot restore detail that was absent from a low-resolution source."
  ) +
  theme_pub(base_size = 10) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major.x = ggplot2::element_blank()
  )

save_publication_plot(
  resolution_plot,
  file.path(figure_dir, "18-raster-resolution-audit"),
  width_mm = 183,
  height_mm = 92,
  dpi = 600
)

parse_pdf_info <- function(path) {
  result <- run_command("pdfinfo", shQuote(path))
  if (result$status != 0L) {
    stop("pdfinfo failed for ", path, call. = FALSE)
  }
  page_line <- grep("^Pages:", result$output, value = TRUE)
  size_line <- grep("^Page size:", result$output, value = TRUE)
  size_match <- regexec(
    "Page size:[[:space:]]+([0-9.]+)[[:space:]]+x[[:space:]]+([0-9.]+)[[:space:]]+pts",
    size_line
  )
  size_parts <- regmatches(size_line, size_match)[[1L]]
  if (length(page_line) != 1L || length(size_parts) != 3L) {
    stop("Unexpected pdfinfo output for ", path, call. = FALSE)
  }
  data.frame(
    pages = as.integer(sub("^Pages:[[:space:]]+", "", page_line)),
    width_pt = as.numeric(size_parts[[2L]]),
    height_pt = as.numeric(size_parts[[3L]]),
    stringsAsFactors = FALSE
  )
}

parse_svg_info <- function(path) {
  svg_text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  root_match <- regexec(
    "<svg[^>]*width='([0-9.]+)pt'[^>]*height='([0-9.]+)pt'",
    svg_text
  )
  root_parts <- regmatches(svg_text, root_match)[[1L]]
  if (length(root_parts) != 3L) {
    stop("Unexpected SVG root dimensions for ", path, call. = FALSE)
  }
  data.frame(
    width_pt = as.numeric(root_parts[[2L]]),
    height_pt = as.numeric(root_parts[[3L]]),
    text_nodes = lengths(regmatches(svg_text, gregexpr("<text\\b", svg_text, perl = TRUE))),
    font_family_present = grepl('font-family: "DejaVu Sans"', svg_text, fixed = TRUE),
    stringsAsFactors = FALSE
  )
}

file_signature_ok <- function(path, extension) {
  raw <- readBin(path, what = "raw", n = 16L)
  if (extension == "pdf") {
    return(identical(rawToChar(raw[seq_len(4L)]), "%PDF"))
  }
  if (extension == "png") {
    return(identical(as.integer(raw[seq_len(8L)]), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
  }
  if (extension == "tiff") {
    signature <- as.integer(raw[seq_len(4L)])
    return(
      identical(signature, c(73L, 73L, 42L, 0L)) ||
        identical(signature, c(77L, 77L, 0L, 42L))
    )
  }
  if (extension == "svg") {
    prefix <- rawToChar(raw)
    return(startsWith(prefix, "<?xml") || startsWith(prefix, "<svg"))
  }
  FALSE
}

format_rows <- list()
for (spec_i in seq_len(nrow(figure_specification))) {
  spec <- figure_specification[spec_i, , drop = FALSE]
  for (extension in c("pdf", "svg", "png", "tiff")) {
    path <- file.path(figure_dir, paste0(spec$stem, ".", extension))
    exists <- file.exists(path)
    bytes <- if (exists) file.info(path)$size else 0
    add_check(
      paste0(spec$stem, "-", extension, "-exists"),
      "export",
      exists,
      TRUE,
      exists
    )
    add_check(
      paste0(spec$stem, "-", extension, "-nonempty"),
      "export",
      bytes,
      ">1000",
      bytes > 1000
    )
    signature_ok <- exists && file_signature_ok(path, extension)
    add_check(
      paste0(spec$stem, "-", extension, "-signature"),
      "export",
      signature_ok,
      TRUE,
      signature_ok
    )

    actual_width_mm <- NA_real_
    actual_height_mm <- NA_real_
    pixel_width <- NA_integer_
    pixel_height <- NA_integer_
    effective_ppi_x <- NA_real_
    effective_ppi_y <- NA_real_
    detected_format <- toupper(extension)
    compression <- NA_character_
    pages <- NA_integer_
    text_nodes <- NA_integer_
    font_status <- NA_character_

    if (extension == "pdf") {
      pdf_info <- parse_pdf_info(path)
      actual_width_mm <- pdf_info$width_pt * 25.4 / 72
      actual_height_mm <- pdf_info$height_pt * 25.4 / 72
      pages <- pdf_info$pages
      font_result <- run_command("pdffonts", shQuote(path))
      font_output <- paste(font_result$output, collapse = "\n")
      font_embedded <- font_result$status == 0L &&
        grepl("DejaVuSans", font_output, fixed = TRUE) &&
        grepl("yes[[:space:]]+yes[[:space:]]+yes", font_output)
      font_status <- ifelse(font_embedded, "embedded", "missing")
      add_check(
        paste0(spec$stem, "-pdf-single-page"),
        "vector",
        pages,
        1,
        pages == 1L
      )
      add_check(
        paste0(spec$stem, "-pdf-font-embedded"),
        "vector",
        font_status,
        "embedded",
        font_embedded
      )
    } else if (extension == "svg") {
      svg_info <- parse_svg_info(path)
      actual_width_mm <- svg_info$width_pt * 25.4 / 72
      actual_height_mm <- svg_info$height_pt * 25.4 / 72
      text_nodes <- svg_info$text_nodes
      font_status <- ifelse(svg_info$font_family_present, "declared", "missing")
      add_check(
        paste0(spec$stem, "-svg-text-nodes"),
        "vector",
        text_nodes,
        ">0",
        text_nodes > 0L
      )
      add_check(
        paste0(spec$stem, "-svg-font-declared"),
        "vector",
        font_status,
        "declared",
        svg_info$font_family_present
      )
    } else {
      raster_info <- identify_raster(path)
      detected_format <- raster_info$detected_format
      pixel_width <- raster_info$pixel_width
      pixel_height <- raster_info$pixel_height
      compression <- raster_info$compression
      actual_width_mm <- spec$width_mm
      actual_height_mm <- spec$height_mm
      effective_ppi_x <- pixel_width / (spec$width_mm / 25.4)
      effective_ppi_y <- pixel_height / (spec$height_mm / 25.4)
      expected_width_px <- floor(spec$width_mm / 25.4 * spec$raster_ppi)
      expected_height_px <- floor(spec$height_mm / 25.4 * spec$raster_ppi)
      add_check(
        paste0(spec$stem, "-", extension, "-pixel-dimensions"),
        "raster",
        c(pixel_width, pixel_height),
        c(expected_width_px, expected_height_px),
        all(
          c(pixel_width, pixel_height) ==
            c(expected_width_px, expected_height_px)
        )
      )
      add_check(
        paste0(spec$stem, "-", extension, "-effective-ppi"),
        "raster",
        round(c(effective_ppi_x, effective_ppi_y), 3),
        paste0(spec$raster_ppi, " +/- 0.5"),
        all(
          abs(c(effective_ppi_x, effective_ppi_y) - spec$raster_ppi) <= 0.5
        )
      )
      add_check(
        paste0(spec$stem, "-", extension, "-format"),
        "raster",
        detected_format,
        toupper(extension),
        identical(detected_format, toupper(extension))
      )
      if (extension == "tiff") {
        add_check(
          paste0(spec$stem, "-tiff-lzw"),
          "raster",
          compression,
          "LZW",
          identical(compression, "LZW")
        )
      }
    }

    dimension_pass <-
      abs(actual_width_mm - spec$width_mm) <= 0.4 &&
      abs(actual_height_mm - spec$height_mm) <= 0.4
    add_check(
      paste0(spec$stem, "-", extension, "-physical-dimensions"),
      "export",
      round(c(actual_width_mm, actual_height_mm), 3),
      c(spec$width_mm, spec$height_mm),
      dimension_pass
    )

    format_rows[[length(format_rows) + 1L]] <- data.frame(
      figure = spec$stem,
      extension = extension,
      path = file.path("figures", paste0(spec$stem, ".", extension)),
      detected_format = detected_format,
      bytes = bytes,
      target_width_mm = spec$width_mm,
      target_height_mm = spec$height_mm,
      actual_width_mm = actual_width_mm,
      actual_height_mm = actual_height_mm,
      pixel_width = pixel_width,
      pixel_height = pixel_height,
      effective_ppi_x = effective_ppi_x,
      effective_ppi_y = effective_ppi_y,
      compression = compression,
      pages = pages,
      text_nodes = text_nodes,
      font_status = font_status,
      dimension_status = ifelse(dimension_pass, "PASS", "FAIL"),
      stringsAsFactors = FALSE
    )
  }
}
format_audit <- do.call(rbind, format_rows)

add_check(
  "primary-figure-count",
  "export",
  length(unique(format_audit$figure)),
  4,
  length(unique(format_audit$figure)) == 4L
)
add_check(
  "primary-format-count",
  "export",
  nrow(format_audit),
  16,
  nrow(format_audit) == 16L
)
add_check(
  "vector-file-count",
  "export",
  sum(format_audit$extension %in% c("pdf", "svg")),
  8,
  sum(format_audit$extension %in% c("pdf", "svg")) == 8L
)
add_check(
  "raster-file-count",
  "export",
  sum(format_audit$extension %in% c("png", "tiff")),
  8,
  sum(format_audit$extension %in% c("png", "tiff")) == 8L
)
add_check(
  "all-physical-dimensions-pass",
  "export",
  sum(format_audit$dimension_status == "PASS"),
  16,
  all(format_audit$dimension_status == "PASS")
)

validation_checks <- do.call(rbind, check_rows)
checks_total <- nrow(validation_checks)
checks_passed <- sum(validation_checks$status == "PASS")
checks_failed <- checks_total - checks_passed

write_tsv(input_audit, "input-audit.tsv")
write_tsv(package_audit, "package-audit.tsv")
write_tsv(command_audit, "command-audit.tsv")
write_tsv(font_audit, "font-audit.tsv")
write_tsv(palette_audit, "palette-audit.tsv")
write_tsv(sample_plot_data, "sample-plot-data.tsv")
write_tsv(phylum_composition, "phylum-composition.tsv")
write_tsv(figure_specification, "figure-specification.tsv")
write_tsv(probe_audit, "resolution-probe-audit.tsv")
write_tsv(format_audit, "format-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

summary_payload <- list(
  article = 18L,
  title = "Publication graphics foundation",
  dataset = list(
    source = "microeco 2.0.0 wetland 16S example",
    features = nrow(otu_matrix),
    samples = ncol(otu_matrix),
    reads = sum(otu_matrix),
    groups = length(group_counts),
    samples_per_group = as.integer(group_counts[c("CW", "IW", "TW")])
  ),
  graphics = list(
    font_family = font_family,
    palette_colors = length(group_palette),
    colour_vision_views = length(palette_views),
    redundant_shapes = length(unique(group_shapes)),
    primary_figures = length(unique(format_audit$figure)),
    vector_files = sum(format_audit$extension %in% c("pdf", "svg")),
    raster_files = sum(format_audit$extension %in% c("png", "tiff")),
    raster_ppi = 600L,
    single_column_width_mm = 89L,
    two_column_width_mm = 183L,
    svg_files_with_text = sum(
      format_audit$extension == "svg" & format_audit$text_nodes > 0
    ),
    pdf_files_with_embedded_font = sum(
      format_audit$extension == "pdf" & format_audit$font_status == "embedded"
    ),
    lzw_tiff_files = sum(
      format_audit$extension == "tiff" & format_audit$compression == "LZW"
    ),
    dimension_checks_passed = sum(format_audit$dimension_status == "PASS")
  ),
  resolution_probe = list(
    width_mm = 89L,
    height_mm = 70L,
    ppi = as.integer(probe_audit$target_ppi),
    pixel_width = as.integer(probe_audit$pixel_width),
    pixel_height = as.integer(probe_audit$pixel_height)
  ),
  versions = as.list(observed_versions),
  checks_total = checks_total,
  checks_passed = checks_passed,
  checks_failed = checks_failed
)
jsonlite::write_json(
  summary_payload,
  file.path(output_dir, "publication-graphics-summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 10
)

session_lines <- sanitize_text(capture.output(utils::sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"))

log_lines <- c(
  "Article 18 publication-graphics validation",
  paste0("R version: ", paste(R.version$major, R.version$minor, sep = ".")),
  paste0("Dataset: ", nrow(otu_matrix), " features x ", ncol(otu_matrix), " samples"),
  paste0("Reads: ", sum(otu_matrix)),
  paste0("Font: ", font_family, " (", sanitize_text(font_path), ")"),
  paste0("Primary exports: ", nrow(format_audit), " files"),
  paste0("Resolution probes: ", paste(probe_audit$target_ppi, collapse = ", "), " ppi"),
  paste0("Checks passed: ", checks_passed, "/", checks_total),
  "",
  paste(
    validation_checks$check_id,
    validation_checks$status,
    validation_checks$observed,
    sep = "\t"
  )
)
writeLines(sanitize_text(log_lines), file.path(output_dir, "validation.log"))

if (checks_failed > 0L) {
  failed <- validation_checks$check_id[validation_checks$status == "FAIL"]
  stop(
    "Article 18 validation failed: ",
    paste(failed, collapse = ", "),
    call. = FALSE
  )
}

message(
  "Article 18 validation passed: ",
  checks_passed,
  "/",
  checks_total,
  " checks."
)
