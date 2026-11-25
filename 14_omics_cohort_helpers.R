# Shared cohort utilities for the T1 DIABLO mediation analyses.

d14c_need <- function(x, columns, label) {
  absent <- setdiff(columns, names(x))
  if (length(absent)) stop(label, " is missing: ", paste(absent, collapse = ", "))
}

d14c_unique_id <- function(x, label) {
  x$id <- as.character(x$id)
  if (anyNA(x$id) || any(!nzchar(x$id)) || anyDuplicated(x$id))
    stop(label, " has missing or duplicated participant IDs")
  x
}

# Load the already processed T1 16S data. Metabolomics data are deliberately
# not merged here, so a bacteria-only mediation model retains the full 16S
# cohort instead of being reduced to the paired multi-omics cohort.
d14c_load_full_t1_16s <- function(outcomes, baseline, pregnancy) {
  files <- c(
    microbe = "data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData",
    population = "data/01_thsbc_cleaned_data/total_pop_5212.RData",
    covariates = "data/01_thsbc_cleaned_data/covariates_5212.RData",
    outcomes = "data/01_thsbc_cleaned_data/outcomes_5212.RData"
  )
  if (any(!file.exists(files)))
    stop("Missing full-cohort input: ", paste(files[!file.exists(files)], collapse = ", "))
  e <- new.env(parent = emptyenv())
  for (file in files) load(file, envir = e)
  required_objects <- c("bact_fung_clrz", "total_pop", "cov_tab", "otc")
  absent <- required_objects[!vapply(required_objects, exists, logical(1),
    envir = e, inherits = FALSE)]
  if (length(absent)) stop("Missing RData object(s): ", paste(absent, collapse = ", "))

  pop <- get("total_pop", envir = e)
  pop$id <- as.character(pop$id)
  pop <- pop[as.character(pop$pbmi_cut3) %in% c("Normal", "OWO"), , drop = FALSE]
  pop$OWO <- as.integer(as.character(pop$pbmi_cut3) == "OWO")
  pop <- d14c_unique_id(pop[, c("id", "OWO"), drop = FALSE], "population")

  micro <- get("bact_fung_clrz", envir = e)
  micro$id <- as.character(micro$id)
  micro <- micro[as.character(micro$period) == "T1", , drop = FALSE]
  micro$period <- NULL
  if ("bact_batch" %in% names(micro)) names(micro)[names(micro) == "bact_batch"] <- "microbe_batch"
  micro <- d14c_unique_id(micro, "T1 16S data")

  cov <- get("cov_tab", envir = e)
  cov$id <- as.character(cov$id)
  cov <- cov[as.character(cov$period) == "T1", , drop = FALSE]
  covariate_names <- unique(c("id", baseline, pregnancy, "wk"))
  d14c_need(cov, covariate_names, "T1 covariates")
  cov <- d14c_unique_id(cov[, covariate_names, drop = FALSE], "T1 covariates")

  outcome_data <- get("otc", envir = e)
  outcome_data$id <- as.character(outcome_data$id)
  d14c_need(outcome_data, c("id", outcomes), "outcomes")
  outcome_data <- d14c_unique_id(
    outcome_data[, c("id", outcomes), drop = FALSE], "outcomes")

  data <- merge(pop, micro, by = "id", all = FALSE, sort = FALSE)
  data <- merge(data, cov, by = "id", all.x = TRUE, sort = FALSE)
  data <- merge(data, outcome_data, by = "id", all.x = TRUE, sort = FALSE)
  data <- data[order(data$id), , drop = FALSE]
  rownames(data) <- NULL
  if ("microbe_batch" %in% names(data)) data$microbe_batch <- factor(data$microbe_batch)
  for (outcome in outcomes)
    data[[outcome]] <- suppressWarnings(as.integer(as.character(data[[outcome]])))
  list(
    data = data,
    files = unname(files),
    input_md5 = tools::md5sum(files),
    source = "processed_T1_16S_full_cohort"
  )
}

# Apply the frozen 14.1 score transformation to another cohort. The loading,
# input centering/scaling and score centering/scaling all come from the paired
# DIABLO refit and are not re-estimated in the larger 16S cohort.
d14c_project_score <- function(data, loadings, block = "Bacteria") {
  z <- loadings[loadings$block == block, , drop = FALSE]
  d14c_need(z, c("feature_id", "loading", "input_center", "input_scale",
    "score_center", "score_scale"), "module_loadings")
  if (!nrow(z) || anyDuplicated(z$feature_id))
    stop("Missing or duplicated ", block, " loadings")
  d14c_need(data, z$feature_id, paste(block, "projection data"))
  numeric_columns <- vapply(data[z$feature_id], is.numeric, logical(1))
  if (!all(numeric_columns))
    stop("Non-numeric score input(s): ", paste(z$feature_id[!numeric_columns], collapse = ", "))
  if (any(!is.finite(z$input_scale)) || any(z$input_scale <= 0) ||
      length(unique(z$score_center)) != 1L || length(unique(z$score_scale)) != 1L ||
      !is.finite(z$score_scale[1]) || z$score_scale[1] <= 0)
    stop("Invalid frozen score scaling parameters")
  x <- as.matrix(data[z$feature_id])
  x <- sweep(x, 2L, z$input_center, "-")
  x <- sweep(x, 2L, z$input_scale, "/")
  raw <- as.numeric(x %*% z$loading)
  (raw - z$score_center[1]) / z$score_scale[1]
}
