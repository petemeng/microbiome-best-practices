#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(permute)
  library(vegan)
})

set.seed(20260726)

root <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root, "data", "small", "multi-kingdom-duran")
out_dir <- file.path(root, "results", "49-multi-kingdom")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_outputs <- c(
  "cross-kingdom-associations.tsv",
  "feature-selection.tsv",
  "ordination-scores.tsv",
  "pcoa-variance.tsv",
  "permanova.tsv",
  "procrustes-arrows.tsv",
  "procrustes-summary.tsv",
  "summary.json"
)
unexpected_outputs <- setdiff(
  list.files(out_dir, all.files = FALSE),
  expected_outputs
)
if (length(unexpected_outputs) > 0L) {
  stop(
    "Refusing to overwrite unexpected files in result directory: ",
    paste(unexpected_outputs, collapse = ", ")
  )
}

metadata <- read.delim(
  file.path(data_dir, "metadata.tsv"),
  row.names = 1,
  check.names = FALSE
)
metadata$Soil <- factor(metadata$Soil, levels = c("GE", "PU", "SD"))
metadata$Compartment <- factor(
  metadata$Compartment,
  levels = c("Soil", "Rhizosphere", "Root")
)

view_files <- c(
  Bacteria = "otutab.tsv",
  Fungi = "fungi-otutab.tsv",
  Oomycete = "oomycete-otutab.tsv"
)
taxonomy_files <- c(
  Bacteria = "taxonomy.tsv",
  Fungi = "fungi-taxonomy.tsv",
  Oomycete = "oomycete-taxonomy.tsv"
)

counts <- lapply(view_files, function(filename) {
  value <- as.matrix(read.delim(
    file.path(data_dir, filename),
    row.names = 1,
    check.names = FALSE
  ))
  storage.mode(value) <- "numeric"
  stopifnot(identical(colnames(value), rownames(metadata)), all(value >= 0))
  value
})

taxonomy <- Map(function(filename, count_matrix) {
  value <- read.delim(
    file.path(data_dir, filename),
    row.names = 1,
    check.names = FALSE
  )
  stopifnot(identical(rownames(value), rownames(count_matrix)))
  value
}, taxonomy_files, counts)

hellinger <- lapply(counts, function(count_matrix) {
  vegan::decostand(t(count_matrix), method = "hellinger")
})
bray <- lapply(hellinger, function(sample_by_feature) {
  vegan::vegdist(sample_by_feature, method = "bray")
})

ordination <- lapply(names(bray), function(view_name) {
  fit <- stats::cmdscale(bray[[view_name]], k = 5, eig = TRUE, add = TRUE)
  scores <- as.data.frame(fit$points)
  colnames(scores) <- paste0("Axis", seq_len(ncol(scores)))
  scores$SampleID <- rownames(scores)
  scores$Kingdom <- view_name
  positive_eigenvalues <- fit$eig[fit$eig > 0]
  variance <- data.frame(
    Kingdom = view_name,
    Axis = paste0("Axis", seq_len(5)),
    Percent = 100 * fit$eig[seq_len(5)] / sum(positive_eigenvalues)
  )
  list(scores = scores, variance = variance)
})
names(ordination) <- names(bray)

ordination_scores <- do.call(rbind, lapply(ordination, function(item) {
  item$scores[, c("SampleID", "Kingdom", paste0("Axis", 1:5))]
}))
ordination_scores <- cbind(
  ordination_scores,
  metadata[ordination_scores$SampleID, c("Soil", "Compartment"), drop = FALSE]
)
pcoa_variance <- do.call(rbind, lapply(ordination, `[[`, "variance"))
write.table(
  ordination_scores,
  file.path(out_dir, "ordination-scores.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  pcoa_variance,
  file.path(out_dir, "pcoa-variance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

set.seed(20260726)
permanova_rows <- lapply(names(bray), function(view_name) {
  fit <- vegan::adonis2(
    bray[[view_name]] ~ Soil * Compartment,
    data = metadata,
    permutations = 999,
    by = "margin"
  )
  fit_table <- as.data.frame(fit)
  fit_table$Term <- rownames(fit_table)
  fit_table$Kingdom <- view_name
  rownames(fit_table) <- NULL
  fit_table[, c("Kingdom", "Term", "Df", "R2", "F", "Pr(>F)")]
})
permanova <- do.call(rbind, permanova_rows)
colnames(permanova)[colnames(permanova) == "Pr(>F)"] <- "PValue"
permanova$QValue <- ave(
  permanova$PValue,
  permanova$Kingdom,
  FUN = function(value) stats::p.adjust(value, method = "BH")
)
write.table(
  permanova,
  file.path(out_dir, "permanova.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

residualize_matrix <- function(value, design) {
  stats::lm.fit(design, value)$residuals
}
niche_design <- stats::model.matrix(~ Soil * Compartment, data = metadata)
ordination_matrices <- lapply(ordination, function(item) {
  as.matrix(item$scores[rownames(metadata), paste0("Axis", 1:5)])
})
ordination_adjusted <- lapply(
  ordination_matrices,
  residualize_matrix,
  design = niche_design
)

view_pairs <- combn(names(ordination_matrices), 2, simplify = FALSE)
procrustes_rows <- list()
arrow_rows <- list()
for (pair in view_pairs) {
  pair_name <- paste(pair, collapse = " vs ")
  for (model_name in c("Observed", "Niche-adjusted")) {
    x <- if (model_name == "Observed") {
      ordination_matrices[[pair[1]]]
    } else {
      ordination_adjusted[[pair[1]]]
    }
    y <- if (model_name == "Observed") {
      ordination_matrices[[pair[2]]]
    } else {
      ordination_adjusted[[pair[2]]]
    }
    permutation_control <- if (model_name == "Observed") {
      permute::how(nperm = 999)
    } else {
      permute::how(
        nperm = 999,
        blocks = interaction(metadata$Soil, metadata$Compartment, drop = TRUE)
      )
    }
    set.seed(20260726 + length(procrustes_rows))
    test <- vegan::protest(
      x,
      y,
      permutations = permutation_control,
      symmetric = TRUE
    )
    procrustes_rows[[length(procrustes_rows) + 1L]] <- data.frame(
      Pair = pair_name,
      Model = model_name,
      Correlation = unname(test$t0),
      PValue = unname(test$signif)
    )

    if (identical(pair, c("Bacteria", "Fungi"))) {
      fit <- vegan::procrustes(x, y, symmetric = TRUE)
      arrow_rows[[length(arrow_rows) + 1L]] <- data.frame(
        SampleID = rownames(metadata),
        Model = model_name,
        StartAxis1 = fit$X[, 1],
        StartAxis2 = fit$X[, 2],
        EndAxis1 = fit$Yrot[, 1],
        EndAxis2 = fit$Yrot[, 2],
        Soil = metadata$Soil,
        Compartment = metadata$Compartment
      )
    }
  }
}
procrustes_summary <- do.call(rbind, procrustes_rows)
procrustes_summary$QValue <- ave(
  procrustes_summary$PValue,
  procrustes_summary$Model,
  FUN = function(value) stats::p.adjust(value, method = "BH")
)
procrustes_arrows <- do.call(rbind, arrow_rows)
write.table(
  procrustes_summary,
  file.path(out_dir, "procrustes-summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  procrustes_arrows,
  file.path(out_dir, "procrustes-arrows.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

clr_feature_matrix <- function(count_matrix) {
  log_counts <- log(t(count_matrix) + 0.5)
  log_counts - rowMeans(log_counts)
}

selected_ids <- list()
clr_selected <- list()
feature_rows <- list()
for (view_name in names(counts)) {
  prevalence <- rowMeans(counts[[view_name]] > 0)
  retained <- names(prevalence)[prevalence >= 0.20]
  clr_all <- clr_feature_matrix(counts[[view_name]][retained, , drop = FALSE])
  feature_variance <- apply(clr_all, 2, stats::var)
  selected <- names(sort(feature_variance, decreasing = TRUE))[seq_len(20)]
  selected_ids[[view_name]] <- selected
  clr_selected[[view_name]] <- clr_all[, selected, drop = FALSE]
  feature_rows[[view_name]] <- data.frame(
    Kingdom = view_name,
    FeatureID = selected,
    DisplayName = taxonomy[[view_name]][selected, "DisplayName"],
    Prevalence = unname(prevalence[selected]),
    CLRVariance = unname(feature_variance[selected])
  )
}
feature_selection <- do.call(rbind, feature_rows)
feature_selection$Label <- paste0(
  feature_selection$DisplayName,
  " (",
  feature_selection$FeatureID,
  ")"
)
write.table(
  feature_selection,
  file.path(out_dir, "feature-selection.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

correlation_table <- function(x, y, view_1, view_2) {
  result <- expand.grid(
    Feature1 = colnames(x),
    Feature2 = colnames(y),
    stringsAsFactors = FALSE
  )
  raw_test <- lapply(seq_len(nrow(result)), function(index) {
    stats::cor.test(
      x[, result$Feature1[index]],
      y[, result$Feature2[index]],
      method = "spearman",
      exact = FALSE
    )
  })
  adjusted_x <- residualize_matrix(x, niche_design)
  adjusted_y <- residualize_matrix(y, niche_design)
  adjusted_test <- lapply(seq_len(nrow(result)), function(index) {
    stats::cor.test(
      adjusted_x[, result$Feature1[index]],
      adjusted_y[, result$Feature2[index]],
      method = "spearman",
      exact = FALSE
    )
  })
  result$Pair <- paste(c(view_1, view_2), collapse = " vs ")
  result$Kingdom1 <- view_1
  result$Kingdom2 <- view_2
  result$Label1 <- feature_selection$Label[
    match(
      paste(view_1, result$Feature1),
      paste(feature_selection$Kingdom, feature_selection$FeatureID)
    )
  ]
  result$Label2 <- feature_selection$Label[
    match(
      paste(view_2, result$Feature2),
      paste(feature_selection$Kingdom, feature_selection$FeatureID)
    )
  ]
  result$RawRho <- vapply(raw_test, `[[`, numeric(1), "estimate")
  result$RawPValue <- vapply(raw_test, `[[`, numeric(1), "p.value")
  result$AdjustedRho <- vapply(adjusted_test, `[[`, numeric(1), "estimate")
  result$AdjustedPValue <- vapply(adjusted_test, `[[`, numeric(1), "p.value")
  result
}

association_rows <- lapply(view_pairs, function(pair) {
  correlation_table(
    clr_selected[[pair[1]]],
    clr_selected[[pair[2]]],
    pair[1],
    pair[2]
  )
})
associations <- do.call(rbind, association_rows)
associations$RawQValue <- stats::p.adjust(associations$RawPValue, method = "BH")
associations$AdjustedQValue <- stats::p.adjust(
  associations$AdjustedPValue,
  method = "BH"
)
associations$SignChanged <- sign(associations$RawRho) != sign(associations$AdjustedRho)

candidate_order <- order(
  associations$AdjustedQValue,
  -abs(associations$AdjustedRho)
)
bootstrap_candidates <- candidate_order[seq_len(min(30, nrow(associations)))]
associations$BootstrapLower <- NA_real_
associations$BootstrapUpper <- NA_real_
associations$SignStability <- NA_real_

strata <- interaction(metadata$Soil, metadata$Compartment, drop = TRUE)
strata_indices <- split(seq_len(nrow(metadata)), strata)
set.seed(20260726)
for (association_index in bootstrap_candidates) {
  feature_1 <- associations$Feature1[association_index]
  feature_2 <- associations$Feature2[association_index]
  view_1 <- associations$Kingdom1[association_index]
  view_2 <- associations$Kingdom2[association_index]
  bootstrap_rho <- replicate(500, {
    sampled <- unlist(lapply(strata_indices, function(index) {
      sample(index, length(index), replace = TRUE)
    }), use.names = FALSE)
    bootstrap_metadata <- metadata[sampled, , drop = FALSE]
    bootstrap_design <- stats::model.matrix(
      ~ Soil * Compartment,
      data = bootstrap_metadata
    )
    residual_1 <- as.numeric(residualize_matrix(
      clr_selected[[view_1]][sampled, feature_1, drop = FALSE],
      bootstrap_design
    ))
    residual_2 <- as.numeric(residualize_matrix(
      clr_selected[[view_2]][sampled, feature_2, drop = FALSE],
      bootstrap_design
    ))
    suppressWarnings(stats::cor(residual_1, residual_2, method = "spearman"))
  })
  bootstrap_rho <- bootstrap_rho[is.finite(bootstrap_rho)]
  associations$BootstrapLower[association_index] <- stats::quantile(
    bootstrap_rho,
    0.025,
    names = FALSE
  )
  associations$BootstrapUpper[association_index] <- stats::quantile(
    bootstrap_rho,
    0.975,
    names = FALSE
  )
  associations$SignStability[association_index] <- mean(
    sign(bootstrap_rho) == sign(associations$AdjustedRho[association_index])
  )
}

associations <- associations[, c(
  "Pair", "Kingdom1", "Feature1", "Label1", "Kingdom2", "Feature2",
  "Label2", "RawRho", "RawPValue", "RawQValue", "AdjustedRho",
  "AdjustedPValue", "AdjustedQValue", "SignChanged", "BootstrapLower",
  "BootstrapUpper", "SignStability"
)]
associations <- associations[order(
  associations$AdjustedQValue,
  -abs(associations$AdjustedRho)
), ]
write.table(
  associations,
  file.path(out_dir, "cross-kingdom-associations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary <- list(
  status = "passed",
  seed = 20260726,
  samples = nrow(metadata),
  balanced_design_cells = nlevels(interaction(metadata$Soil, metadata$Compartment)),
  replicates_per_cell = unname(unique(table(metadata$Soil, metadata$Compartment))),
  input_features = as.list(vapply(counts, nrow, integer(1))),
  selected_features = as.list(vapply(selected_ids, length, integer(1))),
  procrustes_observed_q_lt_0_05 = sum(
    procrustes_summary$Model == "Observed" & procrustes_summary$QValue < 0.05
  ),
  procrustes_adjusted_q_lt_0_05 = sum(
    procrustes_summary$Model == "Niche-adjusted" &
      procrustes_summary$QValue < 0.05
  ),
  tested_cross_kingdom_pairs = nrow(associations),
  raw_cross_kingdom_q_lt_0_05 = sum(associations$RawQValue < 0.05),
  adjusted_cross_kingdom_q_lt_0_05 = sum(associations$AdjustedQValue < 0.05),
  sign_changed_after_adjustment = sum(associations$SignChanged),
  bootstrapped_candidates = sum(is.finite(associations$SignStability)),
  stable_candidates_q_lt_0_05 = sum(
    associations$AdjustedQValue < 0.05 &
      associations$SignStability >= 0.90,
    na.rm = TRUE
  ),
  versions = list(
    R = as.character(getRversion()),
    vegan = as.character(packageVersion("vegan")),
    permute = as.character(packageVersion("permute"))
  )
)
write_json(
  summary,
  file.path(out_dir, "summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
cat(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
