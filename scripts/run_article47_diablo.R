#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mixOmics)
  library(pROC)
  library(jsonlite)
})

set.seed(20260726)

root <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root, "data", "small", "paired-ibd-multiomics")
out_dir <- file.path(root, "results", "47-spls-diablo")
expected_outputs <- c(
  "diablo-selected-loadings.tsv",
  "diablo-selection-stability.tsv",
  "diablo-train-scores.tsv",
  "split-ledger.tsv",
  "spls-correlations.tsv",
  "summary.json",
  "test-auc.tsv",
  "test-confusion.tsv",
  "test-predictions.tsv",
  "tuned-keepx.tsv"
)
existing_outputs <- if (dir.exists(out_dir)) {
  list.files(out_dir, all.files = FALSE)
} else {
  character()
}
unexpected_outputs <- setdiff(existing_outputs, expected_outputs)
if (length(unexpected_outputs) > 0L) {
  stop(
    "Refusing to overwrite unexpected files in result directory: ",
    paste(unexpected_outputs, collapse = ", ")
  )
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

otutab <- read.delim(
  file.path(data_dir, "otutab.tsv"),
  row.names = 1,
  check.names = FALSE
)
taxonomy <- read.delim(
  file.path(data_dir, "taxonomy.tsv"),
  row.names = 1,
  check.names = FALSE
)
metadata <- read.delim(
  file.path(data_dir, "metadata.tsv"),
  row.names = 1,
  check.names = FALSE
)
metabolites <- read.delim(
  file.path(data_dir, "metabolites.tsv"),
  row.names = 1,
  check.names = FALSE
)
metabolite_annotation <- read.delim(
  file.path(data_dir, "metabolite-annotation.tsv"),
  row.names = 1,
  check.names = FALSE
)

otutab <- as.matrix(otutab)
metabolites <- as.matrix(metabolites)
storage.mode(otutab) <- "numeric"
storage.mode(metabolites) <- "numeric"
stopifnot(
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata)),
  identical(colnames(metabolites), rownames(metadata)),
  all(otutab >= 0),
  all(metabolites >= 0)
)

metadata$StudyGroup <- factor(
  metadata$StudyGroup,
  levels = c("Control", "CD", "UC")
)

set.seed(20260726)
train_ids <- unlist(lapply(levels(metadata$StudyGroup), function(group_name) {
  ids <- rownames(metadata)[metadata$StudyGroup == group_name]
  sample(ids, size = floor(0.70 * length(ids)), replace = FALSE)
}), use.names = FALSE)
test_ids <- setdiff(rownames(metadata), train_ids)
train_ids <- rownames(metadata)[rownames(metadata) %in% train_ids]
test_ids <- rownames(metadata)[rownames(metadata) %in% test_ids]

train_group <- droplevels(metadata[train_ids, "StudyGroup"])
test_group <- factor(metadata[test_ids, "StudyGroup"], levels = levels(train_group))

min_train_samples <- ceiling(0.10 * length(train_ids))
prevalent_microbes <- rownames(otutab)[
  rowSums(otutab[, train_ids, drop = FALSE] > 0) >= min_train_samples
]

train_log_counts_all <- log(t(otutab[prevalent_microbes, train_ids, drop = FALSE]) + 0.5)
train_clr_all <- train_log_counts_all - rowMeans(train_log_counts_all)
microbe_variance <- apply(train_clr_all, 2, stats::var)
microbe_features <- names(sort(microbe_variance, decreasing = TRUE))[1:60]

train_log_metabolite_all <- log1p(t(metabolites[, train_ids, drop = FALSE]))
metabolite_variance <- apply(train_log_metabolite_all, 2, stats::var)
metabolite_features <- names(sort(metabolite_variance, decreasing = TRUE))[1:80]

clr_transform <- function(count_matrix) {
  log_matrix <- log(t(count_matrix) + 0.5)
  log_matrix - rowMeans(log_matrix)
}

train_microbe_raw <- clr_transform(
  otutab[microbe_features, train_ids, drop = FALSE]
)
test_microbe_raw <- clr_transform(
  otutab[microbe_features, test_ids, drop = FALSE]
)
train_metabolite_raw <- log1p(t(
  metabolites[metabolite_features, train_ids, drop = FALSE]
))
test_metabolite_raw <- log1p(t(
  metabolites[metabolite_features, test_ids, drop = FALSE]
))

scale_from_train <- function(train_matrix, test_matrix) {
  center <- colMeans(train_matrix)
  spread <- apply(train_matrix, 2, stats::sd)
  spread[!is.finite(spread) | spread == 0] <- 1
  list(
    train = sweep(sweep(train_matrix, 2, center, "-"), 2, spread, "/"),
    test = sweep(sweep(test_matrix, 2, center, "-"), 2, spread, "/"),
    center = center,
    scale = spread
  )
}

microbe_scaled <- scale_from_train(train_microbe_raw, test_microbe_raw)
metabolite_scaled <- scale_from_train(
  train_metabolite_raw,
  test_metabolite_raw
)

train_blocks <- list(
  Microbiome = microbe_scaled$train,
  Metabolome = metabolite_scaled$train
)
test_blocks <- list(
  Microbiome = microbe_scaled$test,
  Metabolome = metabolite_scaled$test
)

spls_fit <- mixOmics::spls(
  train_blocks$Microbiome,
  train_blocks$Metabolome,
  ncomp = 2,
  keepX = c(15, 15),
  keepY = c(20, 20),
  mode = "regression",
  scale = FALSE
)

spls_train_cor <- vapply(1:2, function(component) {
  stats::cor(
    spls_fit$variates$X[, component],
    spls_fit$variates$Y[, component],
    method = "pearson"
  )
}, numeric(1))
spls_test_x <- test_blocks$Microbiome %*% spls_fit$loadings$X
spls_test_y <- test_blocks$Metabolome %*% spls_fit$loadings$Y
spls_test_cor <- vapply(1:2, function(component) {
  stats::cor(
    spls_test_x[, component],
    spls_test_y[, component],
    method = "pearson"
  )
}, numeric(1))
spls_correlations <- data.frame(
  Component = paste0("Component", 1:2),
  TrainCorrelation = spls_train_cor,
  TestCorrelation = spls_test_cor
)
write.table(
  spls_correlations,
  file.path(out_dir, "spls-correlations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

design <- matrix(
  0.10,
  nrow = 2,
  ncol = 2,
  dimnames = list(names(train_blocks), names(train_blocks))
)
diag(design) <- 0

set.seed(20260726)
tuned <- mixOmics::tune.block.splsda(
  X = train_blocks,
  Y = train_group,
  ncomp = 2,
  test.keepX = list(
    Microbiome = c(5, 10, 15),
    Metabolome = c(5, 10, 15)
  ),
  validation = "Mfold",
  folds = 5,
  dist = "centroids.dist",
  measure = "BER",
  nrepeat = 3,
  design = design,
  scale = FALSE,
  progressBar = FALSE
)

keep_x <- tuned$choice.keepX
diablo_fit <- mixOmics::block.splsda(
  X = train_blocks,
  Y = train_group,
  ncomp = 2,
  keepX = keep_x,
  design = design,
  scale = FALSE
)

set.seed(20260726)
cv_perf <- mixOmics::perf(
  diablo_fit,
  validation = "Mfold",
  folds = 5,
  nrepeat = 5,
  dist = "centroids.dist",
  progressBar = FALSE
)

prediction <- predict(
  diablo_fit,
  newdata = test_blocks,
  dist = "centroids.dist"
)
predicted_class <- factor(
  prediction$WeightedVote$centroids.dist[, 2],
  levels = levels(train_group)
)
confusion <- table(Truth = test_group, Predicted = predicted_class)
accuracy <- sum(diag(confusion)) / sum(confusion)
recall <- diag(confusion) / rowSums(confusion)
balanced_accuracy <- mean(recall)
test_ber <- 1 - balanced_accuracy

weighted_scores <- prediction$WeightedPredict[, , 2, drop = FALSE][, , 1]
auc_rows <- lapply(levels(train_group), function(group_name) {
  roc_fit <- pROC::roc(
    response = as.integer(test_group == group_name),
    predictor = weighted_scores[, group_name],
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  data.frame(
    Class = group_name,
    AUC = as.numeric(pROC::auc(roc_fit)),
    CILower = as.numeric(pROC::ci.auc(roc_fit)[1]),
    CIUpper = as.numeric(pROC::ci.auc(roc_fit)[3])
  )
})
auc_table <- do.call(rbind, auc_rows)
macro_auc <- mean(auc_table$AUC)

prediction_table <- data.frame(
  SampleID = test_ids,
  Truth = test_group,
  Predicted = predicted_class,
  weighted_scores,
  check.names = FALSE
)
write.table(
  prediction_table,
  file.path(out_dir, "test-predictions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  as.data.frame.matrix(confusion),
  file.path(out_dir, "test-confusion.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)
write.table(
  auc_table,
  file.path(out_dir, "test-auc.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

score_rows <- lapply(names(train_blocks), function(block_name) {
  data.frame(
    SampleID = train_ids,
    StudyGroup = train_group,
    Block = block_name,
    Component1 = diablo_fit$variates[[block_name]][, 1],
    Component2 = diablo_fit$variates[[block_name]][, 2]
  )
})
diablo_scores <- do.call(rbind, score_rows)
write.table(
  diablo_scores,
  file.path(out_dir, "diablo-train-scores.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

feature_labels <- c(
  setNames(taxonomy[microbe_features, "DisplayName"], microbe_features),
  setNames(
    metabolite_annotation[metabolite_features, "DisplayName"],
    metabolite_features
  )
)
loading_rows <- list()
for (block_name in names(train_blocks)) {
  for (component in 1:2) {
    loading_vector <- diablo_fit$loadings[[block_name]][, component]
    selected <- loading_vector[loading_vector != 0]
    loading_rows[[length(loading_rows) + 1L]] <- data.frame(
      Block = block_name,
      Component = paste0("Component", component),
      FeatureID = names(selected),
      DisplayName = unname(feature_labels[names(selected)]),
      Loading = as.numeric(selected)
    )
  }
}
loading_table <- do.call(rbind, loading_rows)
write.table(
  loading_table,
  file.path(out_dir, "diablo-selected-loadings.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

stability_rows <- list()
for (repeat_name in names(cv_perf$features$stable)) {
  repeat_features <- cv_perf$features$stable[[repeat_name]]
  for (block_name in names(repeat_features)) {
    for (component_name in names(repeat_features[[block_name]])) {
      selection_table <- repeat_features[[block_name]][[component_name]]
      stability_rows[[length(stability_rows) + 1L]] <- data.frame(
        Repeat = repeat_name,
        Block = block_name,
        Component = component_name,
        FeatureID = names(selection_table),
        FoldSelectionRate = as.numeric(selection_table)
      )
    }
  }
}
stability_long <- do.call(rbind, stability_rows)
stability_summary <- aggregate(
  FoldSelectionRate ~ Block + Component + FeatureID,
  data = stability_long,
  FUN = mean
)
stability_summary$DisplayName <- unname(
  feature_labels[stability_summary$FeatureID]
)
write.table(
  stability_summary,
  file.path(out_dir, "diablo-selection-stability.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

training_ledger <- rbind(
  data.frame(Split = "Train", as.list(table(train_group)), check.names = FALSE),
  data.frame(Split = "Test", as.list(table(test_group)), check.names = FALSE)
)
write.table(
  training_ledger,
  file.path(out_dir, "split-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

tuned_ledger <- do.call(rbind, lapply(names(keep_x), function(block_name) {
  data.frame(
    Block = block_name,
    Component = paste0("Component", seq_along(keep_x[[block_name]])),
    KeepX = as.integer(keep_x[[block_name]])
  )
}))
write.table(
  tuned_ledger,
  file.path(out_dir, "tuned-keepx.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cv_weighted_ber <- cv_perf$WeightedVote.error.rate[["centroids.dist"]][
  "Overall.BER",
  "comp2"
]

summary <- list(
  status = "passed",
  seed = 20260726,
  samples = nrow(metadata),
  train_samples = length(train_ids),
  test_samples = length(test_ids),
  train_group_counts = as.list(table(train_group)),
  test_group_counts = as.list(table(test_group)),
  retained_microbe_features = length(microbe_features),
  retained_metabolite_features = length(metabolite_features),
  spls_train_correlations = as.list(unname(spls_train_cor)),
  spls_test_correlations = as.list(unname(spls_test_cor)),
  tune_grid = list(Microbiome = c(5, 10, 15), Metabolome = c(5, 10, 15)),
  tuned_keepX = keep_x,
  repeated_cv_weighted_vote_ber_component2 = unname(cv_weighted_ber),
  test_accuracy = unname(accuracy),
  test_balanced_accuracy = unname(balanced_accuracy),
  test_ber = unname(test_ber),
  test_macro_auc = unname(macro_auc),
  stable_selected_features_rate_ge_0_8 = sum(
    stability_summary$FoldSelectionRate >= 0.80
  ),
  versions = list(
    R = as.character(getRversion()),
    mixOmics = as.character(packageVersion("mixOmics")),
    pROC = as.character(packageVersion("pROC"))
  )
)
write_json(
  summary,
  file.path(out_dir, "summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
cat(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
