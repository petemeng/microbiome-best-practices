# 批次效应：来源、检测与校正策略
# Run sequentially in a new working directory.
# Required packages: MMUPHin, decontam, ggplot2, scales, vegan.

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
library(vegan)
library(ggplot2)

stopifnot(
  identical(
    as.character(utils::packageVersion("decontam")),
    "1.24.0"
  ),
  identical(
    as.character(utils::packageVersion("vegan")),
    "2.6.6.1"
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
  all(c(
    "PlateNumber",
    "Subject",
    "Habitat",
    "quant_reading",
    "IsNegative"
  ) %in% colnames(metadata)),
  all(metadata$quant_reading > 0)
)

dim(otutab)
otutab[1:5, 1:5]
taxonomy[1:5, , drop = FALSE]
head(metadata)

seqtab <- t(otutab)
contam_combined <- decontam::isContaminant(
  seqtab,
  method = "combined",
  conc = metadata$quant_reading,
  neg = metadata$IsNegative,
  threshold = 0.10
)

contaminant_ids <- rownames(contam_combined)[
  contam_combined$contaminant
]
biological_ids <- rownames(metadata)[
  !metadata$IsNegative
]

otutab_bio <- otutab[
  setdiff(rownames(otutab), contaminant_ids),
  biological_ids,
  drop = FALSE
]
taxonomy_bio <- taxonomy[
  rownames(otutab_bio),
  ,
  drop = FALSE
]
metadata_bio <- metadata[
  biological_ids,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(otutab_bio),
    rownames(taxonomy_bio)
  ),
  identical(
    colnames(otutab_bio),
    rownames(metadata_bio)
  )
)

data.frame(
  ContaminantsRemoved = length(contaminant_ids),
  BiologicalSamples = ncol(otutab_bio),
  RetainedASVs = nrow(otutab_bio)
)

habitat_labels <- c(
  "Cheek_Left" = "Cheek left",
  "Cheek_Right" = "Cheek right",
  "Floor" = "Mouth floor",
  "Lip_Lower" = "Lower lip",
  "Lip_Upper" = "Upper lip",
  "Roof" = "Mouth roof",
  "Tongue" = "Tongue"
)
habitat_levels <- unname(habitat_labels)

metadata_bio$Habitat <- factor(
  unname(habitat_labels[metadata_bio$Habitat]),
  levels = habitat_levels
)
metadata_bio$Subject <- factor(
  paste("Subject", metadata_bio$Subject)
)
metadata_bio$Plate <- factor(
  paste("Plate", metadata_bio$PlateNumber),
  levels = paste("Plate", 1:6)
)

plate_subject_table <- with(
  metadata_bio,
  table(Plate, Subject)
)
plate_habitat_table <- with(
  metadata_bio,
  table(Plate, Habitat)
)

plate_subject_table
plate_habitat_table

stopifnot(
  all(rowSums(plate_subject_table > 0) == 1),
  all(colSums(plate_subject_table > 0) == 2),
  all(plate_habitat_table > 0)
)

plate_subject_map <- unique(
  metadata_bio[c("Plate", "Subject")]
)
plate_subject_map <- plate_subject_map[
  order(plate_subject_map$Plate),
  ,
  drop = FALSE
]

design_counts <- as.data.frame(
  plate_habitat_table,
  responseName = "Samples"
)
design_counts$Subject <- plate_subject_map$Subject[
  match(design_counts$Plate, plate_subject_map$Plate)
]
design_counts$PlateSubject <- factor(
  paste(
    design_counts$Plate,
    design_counts$Subject,
    sep = " · "
  ),
  levels = rev(
    paste(
      plate_subject_map$Plate,
      plate_subject_map$Subject,
      sep = " · "
    )
  )
)
head(design_counts)

relative_abundance <- sweep(
  otutab_bio,
  2,
  colSums(otutab_bio),
  "/"
)
bray_distance <- vegan::vegdist(
  t(relative_abundance),
  method = "bray"
)

feature_prevalence <- rowMeans(otutab_bio > 0)
clr_feature_ids <- names(feature_prevalence)[
  feature_prevalence >= 0.05
]
clr_counts <- otutab_bio[
  clr_feature_ids,
  ,
  drop = FALSE
]

clr_matrix <- t(
  apply(
    t(clr_counts) + 0.5,
    1,
    function(x) {
      log(x) - mean(log(x))
    }
  )
)
rownames(clr_matrix) <- colnames(clr_counts)
colnames(clr_matrix) <- rownames(clr_counts)

stopifnot(
  nrow(clr_matrix) == nrow(metadata_bio),
  max(abs(rowSums(clr_matrix))) < 1e-10
)

aitchison_distance <- stats::dist(clr_matrix)
dim(clr_matrix)

run_permanova <- function(distance_object) {
  set.seed(20260719)
  vegan::adonis2(
    distance_object ~ Habitat + Subject + Plate,
    data = metadata_bio,
    permutations = 999,
    by = "terms",
    strata = metadata_bio$Subject
  )
}

extract_permanova <- function(model, metric) {
  term_names <- c("Habitat", "Subject", "Plate")
  data.frame(
    Metric = metric,
    Test = "PERMANOVA",
    Term = term_names,
    Df = model[term_names, "Df"],
    R2 = model[term_names, "R2"],
    F = model[term_names, "F"],
    P = model[term_names, "Pr(>F)"],
    row.names = NULL
  )
}

bray_permanova <- run_permanova(bray_distance)
aitchison_permanova <- run_permanova(
  aitchison_distance
)

community_terms <- rbind(
  extract_permanova(
    bray_permanova,
    "Bray-Curtis"
  ),
  extract_permanova(
    aitchison_permanova,
    "Aitchison"
  )
)
community_terms

run_dispersion <- function(distance_object, metric) {
  dispersion <- vegan::betadisper(
    distance_object,
    metadata_bio$Plate,
    bias.adjust = TRUE
  )
  set.seed(20260719)
  test <- vegan::permutest(
    dispersion,
    permutations = 999
  )
  data.frame(
    Metric = metric,
    Test = "PERMDISP",
    Term = "Plate dispersion",
    Df = unname(test$tab[1, "Df"]),
    R2 = NA_real_,
    F = unname(test$tab[1, "F"]),
    P = unname(test$tab[1, "Pr(>F)"])
  )
}

dispersion_results <- rbind(
  run_dispersion(
    bray_distance,
    "Bray-Curtis"
  ),
  run_dispersion(
    aitchison_distance,
    "Aitchison"
  )
)
dispersion_results

community_diagnostics <- rbind(
  community_terms,
  dispersion_results
)

fit_feature_model <- function(y) {
  fit <- stats::lm(
    y ~ Habitat + Subject + Plate,
    data = metadata_bio
  )
  tab <- stats::anova(fit)
  total_ss <- sum(tab[, "Sum Sq"])
  c(
    HabitatP = tab["Habitat", "Pr(>F)"],
    HabitatR2 = tab["Habitat", "Sum Sq"] / total_ss,
    SubjectP = tab["Subject", "Pr(>F)"],
    SubjectR2 = tab["Subject", "Sum Sq"] / total_ss,
    PlateP = tab["Plate", "Pr(>F)"],
    PlateR2 = tab["Plate", "Sum Sq"] / total_ss
  )
}

feature_model_matrix <- t(
  apply(
    clr_matrix,
    2,
    fit_feature_model
  )
)
feature_model_table <- as.data.frame(
  feature_model_matrix
)
feature_model_table <- feature_model_table[
  colnames(clr_matrix),
  ,
  drop = FALSE
]
feature_model_table$HabitatQ <- stats::p.adjust(
  feature_model_table$HabitatP,
  method = "BH"
)
feature_model_table$PlateQ <- stats::p.adjust(
  feature_model_table$PlateP,
  method = "BH"
)

feature_diagnostics <- data.frame(
  FeatureID = colnames(clr_matrix),
  taxonomy_bio[
    colnames(clr_matrix),
    ,
    drop = FALSE
  ],
  Prevalence = feature_prevalence[
    colnames(clr_matrix)
  ],
  TotalReads = rowSums(
    otutab_bio[
      colnames(clr_matrix),
      ,
      drop = FALSE
    ]
  ),
  feature_model_table,
  check.names = FALSE,
  row.names = NULL
)

feature_summary <- data.frame(
  TestedASVs = nrow(feature_diagnostics),
  HabitatFDR05 = sum(
    feature_diagnostics$HabitatQ < 0.05
  ),
  PlateFDR05 = sum(
    feature_diagnostics$PlateQ < 0.05
  ),
  MedianHabitatR2 = stats::median(
    feature_diagnostics$HabitatR2
  ),
  MedianPlateR2 = stats::median(
    feature_diagnostics$PlateR2
  )
)
feature_summary

leave_one_plate_out <- do.call(
  rbind,
  lapply(
    levels(metadata_bio$Plate),
    function(dropped_plate) {
      keep_samples <- (
        metadata_bio$Plate != dropped_plate
      )
      metadata_subset <- droplevels(
        metadata_bio[
          keep_samples,
          ,
          drop = FALSE
        ]
      )
      model <- vegan::adonis2(
        stats::dist(
          clr_matrix[
            keep_samples,
            ,
            drop = FALSE
          ]
        ) ~ Habitat + Subject + Plate,
        data = metadata_subset,
        permutations = 0,
        by = "terms"
      )
      data.frame(
        DroppedPlate = dropped_plate,
        Samples = sum(keep_samples),
        HabitatR2 = model["Habitat", "R2"],
        SubjectR2 = model["Subject", "R2"],
        PlateR2 = model["Plate", "R2"]
      )
    }
  )
)
leave_one_plate_out

x_biology <- stats::model.matrix(
  ~ Habitat + Subject,
  data = metadata_bio
)

subject_levels <- levels(metadata_bio$Subject)
plate_contrasts <- sapply(
  subject_levels,
  function(subject_level) {
    subject_plates <- sort(
      unique(
        as.character(
          metadata_bio$Plate[
            metadata_bio$Subject == subject_level
          ]
        )
      )
    )
    stopifnot(length(subject_plates) == 2)
    ifelse(
      metadata_bio$Subject == subject_level &
        metadata_bio$Plate == subject_plates[1],
      -0.5,
      ifelse(
        metadata_bio$Subject == subject_level &
          metadata_bio$Plate == subject_plates[2],
        0.5,
        0
      )
    )
  }
)
colnames(plate_contrasts) <- paste0(
  "Within_",
  gsub(" ", "_", subject_levels)
)

x_full <- cbind(
  x_biology,
  plate_contrasts
)
stopifnot(qr(x_full)$rank == ncol(x_full))

all_coefficients <- qr.coef(
  qr(x_full),
  clr_matrix
)
batch_rows <- (
  ncol(x_biology) + 1
):ncol(x_full)
batch_coefficients <- all_coefficients[
  batch_rows,
  ,
  drop = FALSE
]

clr_adjusted_diagnostic <- (
  clr_matrix -
    plate_contrasts %*% batch_coefficients
)

adjustment_coefficients <- data.frame(
  FeatureID = colnames(clr_matrix),
  t(batch_coefficients),
  check.names = FALSE,
  row.names = NULL
)
head(adjustment_coefficients)

adjusted_aitchison <- stats::dist(
  clr_adjusted_diagnostic
)
adjusted_permanova <- run_permanova(
  adjusted_aitchison
)
adjusted_terms <- extract_permanova(
  adjusted_permanova,
  "Aitchison after diagnostic adjustment"
)

adjustment_diagnostic <- rbind(
  transform(
    community_terms[
      community_terms$Metric == "Aitchison",
      ,
      drop = FALSE
    ],
    Stage = "Before adjustment"
  ),
  transform(
    adjusted_terms,
    Stage = "After adjustment"
  )
)
adjustment_diagnostic

residualize_biology <- function(x) {
  stats::lm.fit(
    x = x_biology,
    y = x
  )$residuals
}

residual_before <- residualize_biology(
  clr_matrix
)
residual_after <- residualize_biology(
  clr_adjusted_diagnostic
)

pca_reference <- stats::prcomp(
  residual_before,
  center = TRUE,
  scale. = FALSE
)
scores_before <- pca_reference$x[, 1:2]
scores_after <- (
  scale(
    residual_after,
    center = pca_reference$center,
    scale = FALSE
  ) %*%
    pca_reference$rotation[, 1:2]
)

pca_variance <- (
  pca_reference$sdev^2 /
    sum(pca_reference$sdev^2)
)

ordination_scores <- rbind(
  data.frame(
    SampleID = rownames(metadata_bio),
    Stage = "Before adjustment",
    PC1 = scores_before[, 1],
    PC2 = scores_before[, 2],
    Plate = metadata_bio$Plate,
    Subject = metadata_bio$Subject
  ),
  data.frame(
    SampleID = rownames(metadata_bio),
    Stage = "After adjustment",
    PC1 = scores_after[, 1],
    PC2 = scores_after[, 2],
    Plate = metadata_bio$Plate,
    Subject = metadata_bio$Subject
  )
)
ordination_scores$Stage <- factor(
  ordination_scores$Stage,
  levels = c(
    "Before adjustment",
    "After adjustment"
  )
)

ordination_centroids <- aggregate(
  cbind(PC1, PC2) ~ Stage + Plate,
  data = ordination_scores,
  FUN = mean
)
head(ordination_scores)

result_dir <- "results/05-batch"
dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_result <- function(x, filename) {
  utils::write.table(
    x,
    file.path(result_dir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

write_result(
  design_counts,
  "batch-design-counts.tsv"
)
write_result(
  community_diagnostics,
  "community-batch-diagnostics.tsv"
)
write_result(
  feature_diagnostics,
  "feature-batch-diagnostics.tsv"
)
write_result(
  leave_one_plate_out,
  "leave-one-plate-out.tsv"
)
write_result(
  adjustment_diagnostic,
  "adjustment-diagnostic.tsv"
)
write_result(
  adjustment_coefficients,
  "plate-adjustment-coefficients.tsv"
)
write_result(
  ordination_scores,
  "ordination-scores.tsv"
)

p_batch_design <- ggplot(
  design_counts,
  aes(
    Habitat,
    PlateSubject,
    fill = Samples
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = Samples),
    color = "black",
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_gradient(
    low = "#EAF4FA",
    high = "#0072B2",
    limits = range(design_counts$Samples),
    name = "Samples"
  ) +
  labs(
    title = "Every plate covers all seven oral habitats",
    subtitle = paste0(
      "Plate is nested within subject · ",
      nrow(metadata_bio),
      " biological samples · 6 plates"
    ),
    x = "Oral habitat",
    y = "Technical plate · biological subject"
  ) +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = element_text(
      color = "grey35"
    ),
    axis.text.x = element_text(
      angle = 32,
      hjust = 1
    ),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )

p_batch_design
save_pub(
  p_batch_design,
  "figures/05-batch-design-audit",
  width = 180,
  height = 115
)

community_plot_data <- community_terms
community_plot_data$Term <- factor(
  community_plot_data$Term,
  levels = c("Plate", "Subject", "Habitat")
)

p_community_batch <- ggplot(
  community_plot_data,
  aes(
    R2,
    Term,
    fill = Metric
  )
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.64,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = scales::percent(
        R2,
        accuracy = 0.1
      )
    ),
    position = position_dodge(width = 0.72),
    hjust = -0.12,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Bray-Curtis" = "#0072B2",
      "Aitchison" = "#D55E00"
    ),
    name = "Distance"
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.37),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Plate is detectable but explains little community variation",
    subtitle = "Sequential PERMANOVA · permutations blocked within subject",
    x = "Variance explained (R²)",
    y = NULL
  ) +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = element_text(
      color = "grey35"
    ),
    axis.title.x = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

p_community_batch
save_pub(
  p_community_batch,
  "figures/05-community-batch-effects",
  width = 180,
  height = 115
)

feature_diagnostics$PlateStatus <- factor(
  ifelse(
    feature_diagnostics$PlateQ < 0.05,
    "Plate FDR < 0.05",
    "Not plate-associated"
  ),
  levels = c(
    "Not plate-associated",
    "Plate FDR < 0.05"
  )
)

label_feature_ids <- head(
  feature_diagnostics$FeatureID[
    order(feature_diagnostics$PlateQ)
  ],
  3
)
feature_labels <- feature_diagnostics[
  feature_diagnostics$FeatureID %in%
    label_feature_ids,
  ,
  drop = FALSE
]
feature_labels$LabelHjust <- ifelse(
  feature_labels$PlateR2 > 0.15,
  1.08,
  -0.08
)

p_feature_batch <- ggplot(
  feature_diagnostics,
  aes(
    PlateR2,
    HabitatR2,
    color = PlateStatus
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    linewidth = 0.5,
    color = "grey55"
  ) +
  geom_point(
    size = 2.2,
    alpha = 0.72
  ) +
  geom_text(
    data = feature_labels,
    aes(
      label = FeatureID,
      hjust = LabelHjust
    ),
    nudge_y = 0.012,
    color = "black",
    size = 3.2,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Not plate-associated" = "#B8B8B8",
      "Plate FDR < 0.05" = "#D55E00"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.21),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.45),
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_fixed(ratio = 0.21 / 0.45) +
  labs(
    title = "Batch effects concentrate in a subset of ASVs",
    subtitle = paste0(
      feature_summary$PlateFDR05,
      " of ",
      feature_summary$TestedASVs,
      " ASVs associated with plate after BH correction"
    ),
    x = "Plate variance explained",
    y = "Habitat variance explained"
  ) +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = element_text(
      color = "grey35"
    ),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_feature_batch
save_pub(
  p_feature_batch,
  "figures/05-feature-batch-effects",
  width = 165,
  height = 135
)

plate_colors <- setNames(
  pal_pub[1:6],
  levels(metadata_bio$Plate)
)

p_adjustment <- ggplot(
  ordination_scores,
  aes(
    PC1,
    PC2,
    color = Plate
  )
) +
  geom_point(
    size = 1.25,
    alpha = 0.45
  ) +
  stat_ellipse(
    level = 0.68,
    linewidth = 0.55,
    alpha = 0.75
  ) +
  geom_point(
    data = ordination_centroids,
    shape = 4,
    size = 3.2,
    stroke = 1.1,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ Stage,
    nrow = 1
  ) +
  scale_color_manual(
    values = plate_colors
  ) +
  labs(
    title = "Design-aware residualization removes within-subject plate means",
    subtitle = "Diagnostic only: habitat and subject were preserved in the design matrix",
    x = paste0(
      "Residual PC1 (",
      scales::percent(
        pca_variance[1],
        accuracy = 0.1
      ),
      ")"
    ),
    y = paste0(
      "Residual PC2 (",
      scales::percent(
        pca_variance[2],
        accuracy = 0.1
      ),
      ")"
    ),
    color = "Plate"
  ) +
  theme_pub(base_size = 10.5) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = element_text(
      color = "grey35"
    ),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

p_adjustment
save_pub(
  p_adjustment,
  "figures/05-adjustment-diagnostic",
  width = 200,
  height = 115
)
