# PERMANOVA + betadisper（必配）+ ANOSIM/MRPP
# Run sequentially in a new working directory.
# Required packages: digest, ggplot2, jsonlite, permute, ragg, readr, scales, svglite, vegan.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/permanova-dispersion/metadata.tsv",
  "data/small/permanova-dispersion/otutab.tsv",
  "data/small/permanova-dispersion/source-summary.json",
  "data/small/permanova-dispersion/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(digest)
library(ggplot2)
library(jsonlite)
library(permute)
library(readr)
library(scales)
library(vegan)

set.seed(20260723)

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
    paste0(file_base, ".pdf"), plot,
    width = width, height = height, units = units,
    device = grDevices::cairo_pdf, family = base_family, bg = "white"
  )
  if (isTRUE(write_svg)) {
    ggplot2::ggsave(
      paste0(file_base, ".svg"), plot,
      width = width, height = height, units = units,
      device = svglite::svglite, bg = "white"
    )
  }
  ggplot2::ggsave(
    paste0(file_base, ".png"), plot,
    width = width, height = height, units = units, dpi = dpi,
    device = ragg::agg_png, bg = "white"
  )
  if (isTRUE(write_tiff)) {
    ggplot2::ggsave(
      paste0(file_base, ".tiff"), plot,
      width = width, height = height, units = units, dpi = dpi,
      device = ragg::agg_tiff, compression = "lzw", bg = "white"
    )
  }
  invisible(plot)
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
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
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

input_paths <- c(
  otutab = "data/small/permanova-dispersion/otutab.tsv",
  taxonomy = "data/small/permanova-dispersion/taxonomy.tsv",
  metadata = "data/small/permanova-dispersion/metadata.tsv",
  source_summary = "data/small/permanova-dispersion/source-summary.json"
)
expected_sha256 <- c(
  otutab = "de52dc7457f15f07776c216561803afd92782af7d6723ec93deb8dc0e932b900",
  taxonomy = "599ed73e3d1f59bc2867782c58113100def184c3ebd2cab36286b2fc96dbaee0",
  metadata = "5ebf2fd60582e1a0ef80fe20d251594bec3607df11d8bb74d0981236e03f982c",
  source_summary = "0b7f7c5a6db0c3d9c09d9ae147ff926cb1ce8cbbed804290c074c088100e076d"
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

otutab_frame <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])
source_summary <- jsonlite::read_json(
  input_paths[["source_summary"]],
  simplifyVector = TRUE
)

otutab <- matrix(
  suppressWarnings(as.numeric(as.matrix(otutab_frame))),
  nrow = nrow(otutab_frame),
  ncol = ncol(otutab_frame),
  dimnames = dimnames(otutab_frame)
)
feature_ids <- rownames(otutab)
library_ids <- colnames(otutab)
taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[library_ids, , drop = FALSE]

stopifnot(
  nrow(otutab) == 16825L,
  ncol(otutab) == 56L,
  sum(otutab) == 98022,
  setequal(feature_ids, rownames(taxonomy)),
  setequal(library_ids, rownames(metadata)),
  all(is.finite(otutab)),
  all(otutab >= 0),
  all(otutab == floor(otutab)),
  identical(
    source_summary$source_sha256,
    "47b25a0a033b794fca1db30ae6fdb68c55a44e6c4cd0d1be76788d1d77b51f36"
  )
)

dim(otutab)
sum(otutab)
otutab[1:5, 1:5]
taxonomy[1:5, , drop = FALSE]
head(metadata)

treatment_levels <- c(
  "Ambient + Unclipped", "Ambient + Clipped",
  "Warmed + Unclipped", "Warmed + Clipped"
)
metadata$Block <- factor(metadata$Block, levels = paste("Block", 1:6))
metadata$Warming <- factor(metadata$Warming, levels = c("Ambient", "Warmed"))
metadata$Clipping <- factor(metadata$Clipping, levels = c("Unclipped", "Clipped"))
metadata$Treatment <- factor(metadata$Treatment, levels = treatment_levels)

biological_ids <- unique(as.character(metadata$BiologicalSampleID))
biological_meta <- do.call(rbind, lapply(biological_ids, function(id) {
  rows <- metadata$BiologicalSampleID == id
  fields <- c("Block", "MainPlot", "Warming", "Clipping", "Treatment")
  stopifnot(all(vapply(metadata[rows, fields, drop = FALSE], function(x) {
    length(unique(as.character(x))) == 1L
  }, logical(1))))
  data.frame(
    BiologicalSampleID = id,
    Block = as.character(metadata$Block[rows][[1L]]),
    MainPlot = as.character(metadata$MainPlot[rows][[1L]]),
    Warming = as.character(metadata$Warming[rows][[1L]]),
    Clipping = as.character(metadata$Clipping[rows][[1L]]),
    Treatment = as.character(metadata$Treatment[rows][[1L]]),
    TechnicalLibraries = sum(rows),
    AggregateReads = sum(otutab[, rows, drop = FALSE]),
    stringsAsFactors = FALSE
  )
}))

biological_meta$Block <- factor(
  biological_meta$Block, levels = paste("Block", 1:6)
)
biological_meta$Warming <- factor(
  biological_meta$Warming, levels = c("Ambient", "Warmed")
)
biological_meta$Clipping <- factor(
  biological_meta$Clipping, levels = c("Unclipped", "Clipped")
)
biological_meta$Treatment <- factor(
  biological_meta$Treatment, levels = treatment_levels
)
biological_meta <- biological_meta[
  order(biological_meta$Block, biological_meta$Warming, biological_meta$Clipping),
  , drop = FALSE
]
rownames(biological_meta) <- biological_meta$BiologicalSampleID
biological_ids <- rownames(biological_meta)

aggregate_counts <- vapply(biological_ids, function(id) {
  rowSums(otutab[, metadata$BiologicalSampleID == id, drop = FALSE])
}, numeric(nrow(otutab)))
rownames(aggregate_counts) <- feature_ids
colnames(aggregate_counts) <- biological_ids

stopifnot(
  ncol(aggregate_counts) == 24L,
  sum(aggregate_counts) == sum(otutab),
  all(with(biological_meta, table(Block, Warming, Clipping)) == 1L),
  identical(as.integer(range(biological_meta$TechnicalLibraries)), c(1L, 3L))
)

design_audit <- data.frame(
  BiologicalSampleID = biological_meta$BiologicalSampleID,
  Block = as.character(biological_meta$Block),
  MainPlot = biological_meta$MainPlot,
  Warming = as.character(biological_meta$Warming),
  Clipping = as.character(biological_meta$Clipping),
  Treatment = as.character(biological_meta$Treatment),
  TechnicalLibraries = biological_meta$TechnicalLibraries,
  AggregateReads = biological_meta$AggregateReads,
  stringsAsFactors = FALSE
)
design_audit

presence_by_treatment <- vapply(treatment_levels, function(group) {
  rowSums(otutab[, metadata$Treatment == group, drop = FALSE] > 0)
}, numeric(nrow(otutab)))
colnames(presence_by_treatment) <- treatment_levels
primary_keep <- apply(presence_by_treatment >= 4L, 1L, any)

primary_counts <- aggregate_counts[primary_keep, , drop = FALSE]
primary_community <- t(primary_counts)
primary_relative <- sweep(
  primary_community, 1L, rowSums(primary_community), "/"
)
primary_distance <- vegan::vegdist(primary_relative, method = "bray")
distance_matrix <- as.matrix(primary_distance)

stopifnot(
  sum(primary_keep) == 1824L,
  max(abs(rowSums(primary_relative) - 1)) < 1e-12,
  max(abs(distance_matrix - t(distance_matrix))) < 1e-12,
  max(abs(diag(distance_matrix))) < 1e-12,
  identical(rownames(distance_matrix), biological_ids)
)

data.frame(
  Branch = c(
    "Primary: detected in >=4 tagged libraries in any treatment",
    "Sensitivity: no prevalence filter"
  ),
  InputFeatures = nrow(otutab),
  RetainedFeatures = c(sum(primary_keep), nrow(otutab)),
  UsedForPrimaryInference = c(TRUE, FALSE)
)

make_splitplot_permutations <- function(meta, nset, seed) {
  block_levels <- levels(meta$Block)
  total_possible <- 8^length(block_levels)
  stopifnot(nset <= total_possible - 1L)
  set.seed(seed)
  codes <- sample.int(total_possible - 1L, nset, replace = FALSE)
  permutations <- matrix(
    rep(seq_len(nrow(meta)), each = nset),
    nrow = nset,
    byrow = FALSE
  )

  for (row_i in seq_len(nset)) {
    code <- codes[[row_i]]
    for (block_i in seq_along(block_levels)) {
      state <- code %% 8L
      code <- code %/% 8L
      destination <- which(meta$Block == block_levels[[block_i]])
      ambient <- destination[meta$Warming[destination] == "Ambient"]
      warmed <- destination[meta$Warming[destination] == "Warmed"]
      ambient <- ambient[order(meta$Clipping[ambient])]
      warmed <- warmed[order(meta$Clipping[warmed])]
      if (bitwAnd(state, 2L) != 0L) ambient <- rev(ambient)
      if (bitwAnd(state, 4L) != 0L) warmed <- rev(warmed)
      source <- if (bitwAnd(state, 1L) != 0L) {
        c(warmed, ambient)
      } else {
        c(ambient, warmed)
      }
      permutations[row_i, destination] <- source
    }
  }
  list(matrix = permutations, codes = codes, possible = total_possible)
}

is_legal_splitplot <- function(permutation, meta) {
  for (block_level in levels(meta$Block)) {
    destination <- which(meta$Block == block_level)
    source <- permutation[destination]
    if (!all(meta$Block[source] == block_level)) return(FALSE)

    destination_ambient <- destination[meta$Warming[destination] == "Ambient"]
    destination_warmed <- destination[meta$Warming[destination] == "Warmed"]
    source_plot_a <- unique(as.character(
      meta$MainPlot[permutation[destination_ambient]]
    ))
    source_plot_w <- unique(as.character(
      meta$MainPlot[permutation[destination_warmed]]
    ))
    if (length(source_plot_a) != 1L || length(source_plot_w) != 1L) {
      return(FALSE)
    }
    if (identical(source_plot_a, source_plot_w)) return(FALSE)
    if (!setequal(
      meta$Clipping[permutation[destination_ambient]],
      levels(meta$Clipping)
    )) return(FALSE)
    if (!setequal(
      meta$Clipping[permutation[destination_warmed]],
      levels(meta$Clipping)
    )) return(FALSE)
  }
  TRUE
}

primary_seed <- 20260723L
primary_nperm <- 9999L
splitplot <- make_splitplot_permutations(
  biological_meta, primary_nperm, primary_seed
)
permutation_matrix <- splitplot$matrix

control_design <- permute::how(
  blocks = biological_meta$Block,
  plots = permute::Plots(
    strata = factor(biological_meta$MainPlot),
    type = "free"
  ),
  within = permute::Within(type = "free"),
  nperm = primary_nperm
)
possible_from_permute <- permute::numPerms(
  nrow(biological_meta), control = control_design
)
identity_rows <- sum(apply(permutation_matrix, 1L, function(x) {
  identical(as.integer(x), seq_len(nrow(biological_meta)))
}))
legal_rows <- vapply(seq_len(nrow(permutation_matrix)), function(i) {
  is_legal_splitplot(permutation_matrix[i, ], biological_meta)
}, logical(1))

permutation_contract <- data.frame(
  Seed = primary_seed,
  PossibleIncludingIdentity = splitplot$possible,
  PossibleFromPermute = possible_from_permute,
  PermutationsUsed = nrow(permutation_matrix),
  UniqueRows = nrow(unique(permutation_matrix)),
  IdentityRows = identity_rows,
  LegalRows = sum(legal_rows),
  MinimumP = 1 / (primary_nperm + 1),
  MatrixSHA256 = digest(
    permutation_matrix, algo = "sha256", serialize = TRUE
  )
)
stopifnot(
  splitplot$possible == 262144,
  possible_from_permute == 262144,
  nrow(unique(permutation_matrix)) == 9999L,
  identity_rows == 0L,
  all(legal_rows)
)
permutation_contract

set.seed(primary_seed)
omnibus_fit <- vegan::adonis2(
  primary_distance ~ Block + Treatment,
  data = biological_meta,
  permutations = permutation_matrix,
  by = "terms",
  parallel = 1
)
omnibus_ss <- omnibus_fit["Treatment", "SumOfSqs"]
omnibus_residual_ss <- omnibus_fit["Residual", "SumOfSqs"]
permanova_omnibus <- data.frame(
  Term = "Treatment (4 levels)",
  Df = omnibus_fit["Treatment", "Df"],
  SumOfSquares = omnibus_ss,
  ConditioningBlockR2 = omnibus_fit["Block", "R2"],
  R2Total = omnibus_fit["Treatment", "R2"],
  R2Partial = omnibus_ss / (omnibus_ss + omnibus_residual_ss),
  PseudoF = omnibus_fit["Treatment", "F"],
  PValue = omnibus_fit["Treatment", "Pr(>F)"],
  Permutations = primary_nperm,
  MinimumP = 1 / (primary_nperm + 1)
)

set.seed(primary_seed)
factorial_fit <- vegan::adonis2(
  primary_distance ~ Block + Warming * Clipping,
  data = biological_meta,
  permutations = permutation_matrix,
  by = "terms",
  parallel = 1
)
factorial_terms <- c("Warming", "Clipping", "Warming:Clipping")
factorial_residual_ss <- factorial_fit["Residual", "SumOfSqs"]
permanova_factorial <- data.frame(
  Term = factorial_terms,
  Df = factorial_fit[factorial_terms, "Df"],
  SumOfSquares = factorial_fit[factorial_terms, "SumOfSqs"],
  R2Total = factorial_fit[factorial_terms, "R2"],
  R2Partial = factorial_fit[factorial_terms, "SumOfSqs"] /
    (factorial_fit[factorial_terms, "SumOfSqs"] + factorial_residual_ss),
  PseudoF = factorial_fit[factorial_terms, "F"],
  PValue = factorial_fit[factorial_terms, "Pr(>F)"],
  Permutations = primary_nperm
)
permanova_factorial$PAdjustedHolm <- p.adjust(
  permanova_factorial$PValue, method = "holm"
)
permanova_factorial$RejectHolm05 <-
  permanova_factorial$PAdjustedHolm < 0.05

permanova_omnibus
permanova_factorial

make_paired_permutations <- function(meta) {
  block_levels <- levels(meta$Block)
  permutations <- matrix(
    rep(seq_len(nrow(meta)), each = 2^length(block_levels) - 1L),
    nrow = 2^length(block_levels) - 1L,
    byrow = FALSE
  )
  for (code in seq_len(2^length(block_levels) - 1L)) {
    for (block_i in seq_along(block_levels)) {
      if (bitwAnd(code, bitwShiftL(1L, block_i - 1L)) != 0L) {
        destination <- which(meta$Block == block_levels[[block_i]])
        permutations[code, destination] <- rev(destination)
      }
    }
  }
  permutations
}

planned_contrasts <- list(
  list(
    label = "Warming | Unclipped",
    group_a = "Ambient + Unclipped",
    group_b = "Warmed + Unclipped"
  ),
  list(
    label = "Warming | Clipped",
    group_a = "Ambient + Clipped",
    group_b = "Warmed + Clipped"
  ),
  list(
    label = "Clipping | Ambient",
    group_a = "Ambient + Unclipped",
    group_b = "Ambient + Clipped"
  ),
  list(
    label = "Clipping | Warmed",
    group_a = "Warmed + Unclipped",
    group_b = "Warmed + Clipped"
  )
)

pairwise_rows <- vector("list", length(planned_contrasts))
for (contrast_i in seq_along(planned_contrasts)) {
  contrast <- planned_contrasts[[contrast_i]]
  keep_samples <- biological_meta$Treatment %in%
    c(contrast$group_a, contrast$group_b)
  pair_meta <- biological_meta[keep_samples, , drop = FALSE]
  pair_meta$ContrastGroup <- factor(
    as.character(pair_meta$Treatment),
    levels = c(contrast$group_a, contrast$group_b)
  )
  pair_meta <- pair_meta[
    order(pair_meta$Block, pair_meta$ContrastGroup),
    , drop = FALSE
  ]
  pair_distance <- as.dist(
    distance_matrix[rownames(pair_meta), rownames(pair_meta), drop = FALSE]
  )
  pair_permutations <- make_paired_permutations(pair_meta)
  set.seed(primary_seed + contrast_i)
  pair_fit <- vegan::adonis2(
    pair_distance ~ Block + ContrastGroup,
    data = pair_meta,
    permutations = pair_permutations,
    by = "terms",
    parallel = 1
  )
  pair_ss <- pair_fit["ContrastGroup", "SumOfSqs"]
  pair_residual_ss <- pair_fit["Residual", "SumOfSqs"]
  pairwise_rows[[contrast_i]] <- data.frame(
    Contrast = contrast$label,
    GroupA = contrast$group_a,
    GroupB = contrast$group_b,
    Samples = nrow(pair_meta),
    Blocks = nlevels(droplevels(pair_meta$Block)),
    R2Total = pair_fit["ContrastGroup", "R2"],
    R2Partial = pair_ss / (pair_ss + pair_residual_ss),
    PseudoF = pair_fit["ContrastGroup", "F"],
    PValue = pair_fit["ContrastGroup", "Pr(>F)"],
    Permutations = nrow(pair_permutations),
    MinimumP = 1 / (nrow(pair_permutations) + 1)
  )
}
pairwise_permanova <- do.call(rbind, pairwise_rows)
pairwise_permanova$PAdjustedHolm <- p.adjust(
  pairwise_permanova$PValue, method = "holm"
)
pairwise_permanova$RejectHolm05 <-
  pairwise_permanova$PAdjustedHolm < 0.05
pairwise_permanova

dispersion_model <- vegan::betadisper(
  primary_distance,
  group = biological_meta$Treatment,
  type = "median",
  bias.adjust = TRUE
)
set.seed(primary_seed)
dispersion_test <- vegan::permutest(
  dispersion_model,
  permutations = permutation_matrix,
  pairwise = FALSE,
  parallel = 1
)
dispersion_global <- data.frame(
  Grouping = "Treatment (4 levels)",
  Center = "Spatial median",
  BiasAdjustment = TRUE,
  FValue = dispersion_test$tab["Groups", "F"],
  PValue = dispersion_test$tab["Groups", "Pr(>F)"],
  Permutations = dispersion_test$tab["Groups", "N.Perm"],
  MinimumP = 1 / (primary_nperm + 1)
)

dispersion_distances <- data.frame(
  BiologicalSampleID = names(dispersion_model$distances),
  Treatment = as.character(biological_meta[
    names(dispersion_model$distances), "Treatment"
  ]),
  DistanceToSpatialMedian = as.numeric(dispersion_model$distances)
)

set.seed(primary_seed)
anosim_fit <- vegan::anosim(
  primary_distance,
  grouping = biological_meta$Treatment,
  permutations = permutation_matrix,
  parallel = 1
)
set.seed(primary_seed)
mrpp_fit <- vegan::mrpp(
  primary_distance,
  grouping = biological_meta$Treatment,
  permutations = permutation_matrix,
  weight.type = 1,
  parallel = 1
)
alternative_tests <- data.frame(
  Method = c("PERMANOVA", "ANOSIM", "MRPP"),
  Statistic = c("Pseudo-F", "R", "A"),
  StatisticValue = c(
    permanova_omnibus$PseudoF,
    unname(anosim_fit$statistic),
    unname(mrpp_fit$A)
  ),
  PValue = c(
    permanova_omnibus$PValue,
    anosim_fit$signif,
    mrpp_fit$Pvalue
  ),
  Role = c("Primary", "Directional sensitivity", "Directional sensitivity")
)

dispersion_global
alternative_tests

make_free_permutations <- function(n, nset, seed) {
  set.seed(seed)
  out <- t(replicate(nset, sample.int(n), simplify = "matrix"))
  stopifnot(nrow(unique(out)) == nset)
  out
}

run_omnibus <- function(distance, meta, permutations, include_block = TRUE) {
  if (include_block) {
    fit <- vegan::adonis2(
      distance ~ Block + Treatment,
      data = meta, permutations = permutations,
      by = "terms", parallel = 1
    )
  } else {
    fit <- vegan::adonis2(
      distance ~ Treatment,
      data = meta, permutations = permutations,
      by = "terms", parallel = 1
    )
  }
  ss <- fit["Treatment", "SumOfSqs"]
  residual_ss <- fit["Residual", "SumOfSqs"]
  c(
    R2Total = fit["Treatment", "R2"],
    R2Partial = ss / (ss + residual_ss),
    PseudoF = fit["Treatment", "F"],
    PValue = fit["Treatment", "Pr(>F)"]
  )
}

unfiltered_relative <- sweep(
  t(aggregate_counts), 1L, rowSums(t(aggregate_counts)), "/"
)
unfiltered_distance <- vegan::vegdist(
  unfiltered_relative, method = "bray"
)
unfiltered_result <- run_omnibus(
  unfiltered_distance, biological_meta, permutation_matrix, TRUE
)

free_biological_permutations <- make_free_permutations(
  nrow(biological_meta), primary_nperm, primary_seed + 1L
)
free_biological_result <- run_omnibus(
  primary_distance, biological_meta, free_biological_permutations, TRUE
)

library_relative <- sweep(
  t(otutab[primary_keep, , drop = FALSE]),
  1L,
  colSums(otutab[primary_keep, , drop = FALSE]),
  "/"
)
library_distance <- vegan::vegdist(library_relative, method = "bray")
free_library_permutations <- make_free_permutations(
  nrow(metadata), primary_nperm, primary_seed + 2L
)
technical_result <- run_omnibus(
  library_distance, metadata, free_library_permutations, FALSE
)

sensitivity_analysis <- data.frame(
  Branch = c(
    "Paper filter + split-plot restriction",
    "No prevalence filter + split-plot restriction",
    "Paper filter + unrestricted biological units",
    "Paper filter + unrestricted technical libraries"
  ),
  Features = c(sum(primary_keep), nrow(otutab), sum(primary_keep), sum(primary_keep)),
  Units = c(24L, 24L, 24L, 56L),
  ValidForPrimaryInference = c(TRUE, TRUE, FALSE, FALSE),
  R2Total = c(
    permanova_omnibus$R2Total,
    unfiltered_result[["R2Total"]],
    free_biological_result[["R2Total"]],
    technical_result[["R2Total"]]
  ),
  R2Partial = c(
    permanova_omnibus$R2Partial,
    unfiltered_result[["R2Partial"]],
    free_biological_result[["R2Partial"]],
    technical_result[["R2Partial"]]
  ),
  PseudoF = c(
    permanova_omnibus$PseudoF,
    unfiltered_result[["PseudoF"]],
    free_biological_result[["PseudoF"]],
    technical_result[["PseudoF"]]
  ),
  PValue = c(
    permanova_omnibus$PValue,
    unfiltered_result[["PValue"]],
    free_biological_result[["PValue"]],
    technical_result[["PValue"]]
  )
)
sensitivity_analysis

treatment_palette <- c(
  "Ambient + Unclipped" = pal_pub[["blue"]],
  "Ambient + Clipped" = pal_pub[["sky"]],
  "Warmed + Unclipped" = pal_pub[["vermillion"]],
  "Warmed + Clipped" = pal_pub[["orange"]]
)
warming_palette <- c(
  "Ambient" = pal_pub[["blue"]],
  "Warmed" = pal_pub[["vermillion"]]
)
clipping_shapes <- c("Unclipped" = 21, "Clipped" = 24)

design_plot_data <- design_audit
design_plot_data$Block <- factor(
  design_plot_data$Block, levels = rev(paste("Block", 1:6))
)
design_plot_data$Treatment <- factor(
  design_plot_data$Treatment, levels = treatment_levels
)
design_plot_data$Warming <- factor(
  design_plot_data$Warming, levels = c("Ambient", "Warmed")
)
design_plot_data$CellLabel <- paste0(
  design_plot_data$TechnicalLibraries,
  ifelse(design_plot_data$TechnicalLibraries == 1L, " library", " libraries"),
  "\n1 analysis unit"
)

design_plot <- ggplot(
  design_plot_data, aes(Treatment, Block, fill = Warming)
) +
  geom_tile(
    colour = "white", linewidth = 1.1, width = 0.96, height = 0.92
  ) +
  geom_text(
    aes(label = CellLabel), colour = "white", size = 2.65,
    lineheight = 0.9, family = font_pub
  ) +
  scale_fill_manual(values = warming_palette, drop = FALSE) +
  scale_x_discrete(labels = c(
    "Ambient + Unclipped" = "Ambient\nUnclipped",
    "Ambient + Clipped" = "Ambient\nClipped",
    "Warmed + Unclipped" = "Warmed\nUnclipped",
    "Warmed + Clipped" = "Warmed\nClipped"
  )) +
  labs(
    title = "Technical libraries collapse to 24 biological units",
    subtitle = paste(
      "One row is one block; whole plots may swap within blocks,",
      "and clipping may swap only within each main plot"
    ),
    x = "Field treatment", y = NULL, fill = "Whole-plot factor",
    caption = paste(
      "Labels show 1–3 PCR/tag libraries aggregated before inference.",
      "The resulting split-plot space contains 8^6 = 262,144 assignments."
    )
  ) +
  theme_pub(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8.3, lineheight = 0.9),
    axis.text.y = element_text(size = 8.3),
    legend.position = "top"
  )

save_pub(
  design_plot, "figures/23-exchangeability-design",
  width = 183, height = 104
)
design_plot

site_vectors <- as.data.frame(
  dispersion_model$vectors[, seq_len(2L), drop = FALSE]
)
colnames(site_vectors) <- c("PCoA1", "PCoA2")
site_vectors$BiologicalSampleID <- rownames(site_vectors)
site_vectors$Treatment <- as.character(biological_meta[
  site_vectors$BiologicalSampleID, "Treatment"
])
site_vectors$Clipping <- as.character(biological_meta[
  site_vectors$BiologicalSampleID, "Clipping"
])

center_vectors <- as.data.frame(
  dispersion_model$centroids[, seq_len(2L), drop = FALSE]
)
colnames(center_vectors) <- c("Center1", "Center2")
center_vectors$Treatment <- rownames(center_vectors)
site_vectors$Center1 <- center_vectors$Center1[
  match(site_vectors$Treatment, center_vectors$Treatment)
]
site_vectors$Center2 <- center_vectors$Center2[
  match(site_vectors$Treatment, center_vectors$Treatment)
]
positive_eigenvalues <- dispersion_model$eig[dispersion_model$eig > 0]
axis_percent <- 100 * dispersion_model$eig[seq_len(2L)] /
  sum(positive_eigenvalues)

pcoa_plot <- ggplot(site_vectors, aes(PCoA1, PCoA2)) +
  geom_segment(
    aes(xend = Center1, yend = Center2, colour = Treatment),
    linewidth = 0.42, alpha = 0.55, show.legend = FALSE
  ) +
  geom_point(
    aes(fill = Treatment, shape = Clipping),
    colour = "#1A1A1A", size = 3.0, stroke = 0.55
  ) +
  geom_point(
    data = center_vectors,
    aes(Center1, Center2, colour = Treatment),
    inherit.aes = FALSE, shape = 4, size = 4.2, stroke = 1.2,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = treatment_palette, drop = FALSE) +
  scale_colour_manual(values = treatment_palette, drop = FALSE) +
  scale_shape_manual(values = clipping_shapes, drop = FALSE) +
  labs(
    title = "Location and spread must be read together",
    subtitle = sprintf(
      "Blocked PERMANOVA: R² = %.3f, P = %s; PERMDISP: P = %s",
      permanova_omnibus$R2Total,
      format_p(permanova_omnibus$PValue),
      format_p(dispersion_global$PValue)
    ),
    x = sprintf("PCoA 1 (%.1f%%)", axis_percent[[1L]]),
    y = sprintf("PCoA 2 (%.1f%%)", axis_percent[[2L]]),
    fill = "Treatment", shape = "Clipping",
    caption = paste(
      "Crosses are treatment spatial medians; segments show the two-dimensional projection.",
      "PERMDISP uses full-space distances to spatial medians.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  theme(legend.position = "top", legend.box = "vertical") +
  guides(
    fill = guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(shape = 21, size = 3.0, colour = "#1A1A1A")
    ),
    shape = guide_legend(nrow = 1)
  )

save_pub(
  pcoa_plot, "figures/23-pcoa-centroid-dispersion",
  width = 183, height = 120
)
pcoa_plot

pairwise_plot_data <- pairwise_permanova
pairwise_plot_data$Contrast <- factor(
  pairwise_plot_data$Contrast,
  levels = rev(vapply(planned_contrasts, `[[`, character(1), "label"))
)
pairwise_plot_data$Decision <- ifelse(
  pairwise_plot_data$RejectHolm05,
  "Holm-adjusted P < 0.05",
  "Holm-adjusted P >= 0.05"
)
pairwise_plot_data$Label <- sprintf(
  "Partial R² = %.3f  |  Holm P %s",
  pairwise_plot_data$R2Partial,
  format_p(pairwise_plot_data$PAdjustedHolm)
)
pairwise_decision_palette <- c(
  "Holm-adjusted P < 0.05" = pal_pub[["vermillion"]],
  "Holm-adjusted P >= 0.05" = pal_pub[["grey"]]
)

pairwise_plot <- ggplot(
  pairwise_plot_data, aes(R2Partial, Contrast, colour = Decision)
) +
  geom_segment(
    aes(x = 0, xend = R2Partial, yend = Contrast),
    linewidth = 0.7, show.legend = FALSE
  ) +
  geom_point(size = 3.2) +
  geom_text(
    aes(label = Label), hjust = -0.08, colour = "#333333",
    size = 2.65, family = font_pub
  ) +
  scale_colour_manual(values = pairwise_decision_palette) +
  scale_x_continuous(
    limits = c(0, max(pairwise_plot_data$R2Partial) * 1.75),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    title = "Planned simple effects use exact paired permutations",
    subtitle = "Each contrast enumerates 63 non-identity swaps across six blocks",
    x = "Partial variance explained", y = NULL,
    colour = "Multiplicity decision",
    caption = paste(
      "The four prespecified P values form one Holm family.",
      "With 64 assignments including the observed order, the minimum P is 0.015625.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  theme(panel.grid.major.y = element_blank(), legend.position = "top")

save_pub(
  pairwise_plot, "figures/23-pairwise-permanova",
  width = 183, height = 96
)
pairwise_plot

dispersion_plot_data <- dispersion_distances
dispersion_plot_data$Treatment <- factor(
  dispersion_plot_data$Treatment, levels = treatment_levels
)

dispersion_plot <- ggplot(
  dispersion_plot_data,
  aes(Treatment, DistanceToSpatialMedian, fill = Treatment)
) +
  geom_boxplot(
    width = 0.58, outlier.shape = NA, alpha = 0.35,
    colour = "#333333", linewidth = 0.45
  ) +
  geom_point(
    shape = 21, colour = "#1A1A1A", size = 2.4, stroke = 0.45,
    position = position_jitter(width = 0.09, height = 0, seed = primary_seed)
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 4, size = 3.4,
    stroke = 1.0, colour = "#1A1A1A"
  ) +
  scale_fill_manual(values = treatment_palette, guide = "none") +
  scale_x_discrete(labels = c(
    "Ambient + Unclipped" = "Ambient\nUnclipped",
    "Ambient + Clipped" = "Ambient\nClipped",
    "Warmed + Unclipped" = "Warmed\nUnclipped",
    "Warmed + Clipped" = "Warmed\nClipped"
  )) +
  labs(
    title = "PERMDISP tests spread, not group location",
    subtitle = sprintf(
      "Spatial median, bias-adjusted; F = %.3f, restricted P = %s",
      dispersion_global$FValue, format_p(dispersion_global$PValue)
    ),
    x = "Treatment", y = "Distance to spatial median",
    caption = paste(
      "Points are 24 biological soil samples; crosses are group means of the distances.",
      "A non-significant result is not proof that dispersions are identical.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 8.2, lineheight = 0.9)
  )

save_pub(
  dispersion_plot, "figures/23-dispersion-by-treatment",
  width = 183, height = 104
)
dispersion_plot

audit_summary_path <- paste0(
  "results/23-permanova-dispersion/",
  "permanova-dispersion-summary.json"
)
if (file.exists(audit_summary_path)) {
  audit_summary <- jsonlite::read_json(
    audit_summary_path, simplifyVector = TRUE
  )
  stopifnot(
    audit_summary$checks_total == 214L,
    audit_summary$checks_passed == 214L,
    audit_summary$checks_failed == 0L,
    audit_summary$dataset$features == 16825L,
    audit_summary$dataset$technical_libraries == 56L,
    audit_summary$dataset$biological_samples == 24L,
    audit_summary$standardization$features_retained == 1824L,
    audit_summary$permutations$possible == 262144L,
    audit_summary$permutations$primary == 9999L,
    audit_summary$analysis$factorial_rejections_holm == 0L,
    audit_summary$analysis$pairwise_rejections_holm == 0L,
    audit_summary$graphics$primary_figures == 4L,
    audit_summary$graphics$format_files == 16L
  )
  unlist(audit_summary[c("checks_total", "checks_passed", "checks_failed")])
}
