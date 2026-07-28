#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(FEAST)
  library(jsonlite)
})

set.seed(20260726)

root <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root, "data", "small", "source-tracking-feast")
out_dir <- file.path(root, "results", "50-source-tracking")
native_dir <- file.path(out_dir, "feast-native")
dir.create(native_dir, recursive = TRUE, showWarnings = FALSE)

expected_outputs <- c(
  "feast-loo-predictions.tsv",
  "feast-main-class.tsv",
  "feast-main-sample.tsv",
  "feast-missing-source.tsv",
  "feast-native",
  "source-tracker-loo-predictions.tsv",
  "source-tracker-main.tsv",
  "summary-feast.json",
  "summary-source-tracker.json"
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

otutab <- as.matrix(read.delim(
  file.path(data_dir, "otutab.tsv"),
  row.names = 1,
  check.names = FALSE
))
storage.mode(otutab) <- "numeric"
metadata <- read.delim(
  file.path(data_dir, "metadata.tsv"),
  row.names = 1,
  check.names = FALSE
)
stopifnot(
  identical(colnames(otutab), rownames(metadata)),
  all(otutab >= 0),
  all(otutab == floor(otutab)),
  sum(metadata$SourceSink == "Sink") == 1
)

count_by_sample <- t(otutab)
source_ids <- rownames(metadata)[metadata$SourceSink == "Source"]
sink_id <- rownames(metadata)[metadata$SourceSink == "Sink"]
source_classes <- unique(metadata[source_ids, "SourceClass"])

aggregate_contributions <- function(sample_contributions, sample_metadata) {
  classes <- c(unique(sample_metadata$SourceClass), "Unknown")
  values <- setNames(numeric(length(classes)), classes)
  for (source_name in names(sample_contributions)) {
    if (source_name == "Unknown") {
      values["Unknown"] <- values["Unknown"] + sample_contributions[source_name]
    } else {
      sample_id <- sub("_.*$", "", source_name)
      source_class <- sample_metadata[sample_id, "SourceClass"]
      values[source_class] <- values[source_class] + sample_contributions[source_name]
    }
  }
  values
}

run_feast <- function(sample_ids, sink_sample, run_label, seed) {
  run_dir <- file.path(native_dir, run_label)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_metadata <- metadata[sample_ids, , drop = FALSE]
  run_metadata$SourceSink <- ifelse(
    rownames(run_metadata) == sink_sample,
    "Sink",
    "Source"
  )
  feast_metadata <- data.frame(
    Env = run_metadata$SourceClass,
    SourceSink = run_metadata$SourceSink,
    id = 1,
    row.names = rownames(run_metadata),
    check.names = FALSE
  )
  output_file <- file.path(
    run_dir,
    paste0(run_label, "_source_contributions_matrix.txt")
  )
  set.seed(seed)
  previous_directory <- getwd()
  on.exit(setwd(previous_directory), add = TRUE)
  invisible(FEAST::FEAST(
    C = count_by_sample[sample_ids, , drop = FALSE],
    metadata = feast_metadata,
    EM_iterations = 1000,
    COVERAGE = 5000,
    different_sources_flag = 0,
    dir_path = run_dir,
    outfile = run_label
  ))
  setwd(previous_directory)
  result <- read.delim(
    output_file,
    row.names = 1,
    check.names = FALSE
  )
  stopifnot(nrow(result) == 1, abs(sum(result[1, ]) - 1) < 1e-6)
  unlist(result[1, ], use.names = TRUE)
}

main_samples <- c(sink_id, source_ids)
main_sample_contribution <- run_feast(
  main_samples,
  sink_id,
  "main",
  20260726
)
main_class_contribution <- aggregate_contributions(
  main_sample_contribution,
  metadata[source_ids, , drop = FALSE]
)

main_sample_table <- data.frame(
  SinkSampleID = sink_id,
  Source = names(main_sample_contribution),
  Contribution = as.numeric(main_sample_contribution)
)
main_sample_table$SourceSampleID <- ifelse(
  main_sample_table$Source == "Unknown",
  NA_character_,
  sub("_.*$", "", main_sample_table$Source)
)
main_sample_table$SourceClass <- ifelse(
  main_sample_table$Source == "Unknown",
  "Unknown",
  metadata[main_sample_table$SourceSampleID, "SourceClass"]
)
main_class_table <- data.frame(
  SinkSampleID = sink_id,
  SourceClass = names(main_class_contribution),
  Contribution = as.numeric(main_class_contribution)
)
write.table(
  main_sample_table,
  file.path(out_dir, "feast-main-sample.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  main_class_table,
  file.path(out_dir, "feast-main-class.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

class_sizes <- table(metadata[source_ids, "SourceClass"])
eligible_loo_ids <- source_ids[
  class_sizes[metadata[source_ids, "SourceClass"]] >= 2
]
loo_rows <- lapply(seq_along(eligible_loo_ids), function(index) {
  held_out <- eligible_loo_ids[index]
  training_sources <- setdiff(source_ids, held_out)
  sample_contribution <- run_feast(
    c(held_out, training_sources),
    held_out,
    paste0("loo-", held_out),
    20260726 + index
  )
  class_contribution <- aggregate_contributions(
    sample_contribution,
    metadata[training_sources, , drop = FALSE]
  )
  all_classes <- c(source_classes, "Unknown")
  aligned <- setNames(numeric(length(all_classes)), all_classes)
  aligned[names(class_contribution)] <- class_contribution
  known_values <- aligned[source_classes]
  predicted_known <- names(which.max(known_values))
  predicted_with_unknown <- names(which.max(aligned))
  data.frame(
    SampleID = held_out,
    TrueClass = metadata[held_out, "SourceClass"],
    PredictedKnownClass = predicted_known,
    PredictedIncludingUnknown = predicted_with_unknown,
    CorrectKnownClass = predicted_known == metadata[held_out, "SourceClass"],
    InfantGut = unname(aligned["Infant gut"]),
    AdultGut = unname(aligned["Adult gut"]),
    AdultSkin = unname(aligned["Adult skin"]),
    Soil = unname(aligned["Soil"]),
    Unknown = unname(aligned["Unknown"]),
    check.names = FALSE
  )
})
loo_predictions <- do.call(rbind, loo_rows)
write.table(
  loo_predictions,
  file.path(out_dir, "feast-loo-predictions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

missing_rows <- list(data.frame(
  RemovedClass = "None",
  SourceClass = names(main_class_contribution),
  Contribution = as.numeric(main_class_contribution)
))
for (class_index in seq_along(source_classes)) {
  removed_class <- source_classes[class_index]
  retained_sources <- source_ids[
    metadata[source_ids, "SourceClass"] != removed_class
  ]
  sample_contribution <- run_feast(
    c(sink_id, retained_sources),
    sink_id,
    paste0("remove-", gsub(" ", "-", tolower(removed_class))),
    20260800 + class_index
  )
  class_contribution <- aggregate_contributions(
    sample_contribution,
    metadata[retained_sources, , drop = FALSE]
  )
  missing_rows[[length(missing_rows) + 1L]] <- data.frame(
    RemovedClass = removed_class,
    SourceClass = names(class_contribution),
    Contribution = as.numeric(class_contribution)
  )
}
missing_source <- do.call(rbind, missing_rows)
write.table(
  missing_source,
  file.path(out_dir, "feast-missing-source.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

loo_accuracy <- mean(loo_predictions$CorrectKnownClass)
unknown_by_scenario <- missing_source[
  missing_source$SourceClass == "Unknown",
  c("RemovedClass", "Contribution")
]
summary <- list(
  status = "passed",
  seed = 20260726,
  samples = nrow(metadata),
  sink_samples = length(sink_id),
  source_samples = length(source_ids),
  source_classes = as.list(source_classes),
  taxa = nrow(otutab),
  coverage = 5000,
  em_iterations = 1000,
  main_largest_source_class = names(which.max(
    main_class_contribution[names(main_class_contribution) != "Unknown"]
  )),
  main_unknown = unname(main_class_contribution["Unknown"]),
  leave_one_out_samples = nrow(loo_predictions),
  leave_one_out_known_class_accuracy = unname(loo_accuracy),
  missing_source_unknown = as.list(setNames(
    unknown_by_scenario$Contribution,
    unknown_by_scenario$RemovedClass
  )),
  versions = list(
    R = as.character(getRversion()),
    FEAST = as.character(packageVersion("FEAST"))
  )
)
write_json(
  summary,
  file.path(out_dir, "summary-feast.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
cat(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
