# Shared implementation for 14.2 T1 score-level mediation.
#
# Scientific scope
#   1. OWO -> bacterial DIABLO score -> outcome
#      OWO -> metabolite DIABLO score -> outcome
#   2. Bacterial score -> metabolite score -> outcome, and the reverse direction
#      (OWO is adjusted in both regression models)


if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  invisible(suppressWarnings(try(
    Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE
  )))
}

m142_bind <- function(xs) {
  xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
  if (!length(xs)) return(data.frame())
  columns <- unique(unlist(lapply(xs, names)))
  xs <- lapply(xs, function(x) {
    for (column in setdiff(columns, names(x))) x[[column]] <- NA
    x[columns]
  })
  ans <- do.call(rbind, xs)
  rownames(ans) <- NULL
  ans
}

m142_hash <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 2L)
  unname(tools::md5sum(path))
}

m142_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  saveRDS(x, temporary, compress = TRUE)
  if (file.exists(path) && !file.remove(path)) stop("Cannot replace: ", path)
  if (!file.rename(temporary, path)) stop("Cannot save: ", path)
}

m142_need <- function(x, columns, label) {
  absent <- setdiff(columns, names(x))
  if (length(absent)) stop(label, " is missing: ", paste(absent, collapse = ", "))
}

m142_binary <- function(x, label) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) x <- as.integer(x)
  x <- suppressWarnings(as.numeric(x))
  if (anyNA(x) || !all(x %in% 0:1) || length(unique(x)) != 2L)
    stop(label, " must contain both explicitly coded values 0 and 1")
  x
}

m142_config <- function(outcomes = c("GDM", "HDP"), workers = 8L,
                        nsims = 5000L) {
  outcomes <- unique(as.character(outcomes))
  if (!length(outcomes) || any(!outcomes %in% c("GDM", "HDP")))
    stop("outcomes must be GDM, HDP, or both")
  outcome_tag <- paste(outcomes, collapse = "_")
  list(
    version = "14.2_score_module_mediation_v3.0_mixed_cohorts",
    outcomes = outcomes,
    input = "result/PART4/DIABLO/14.1_DIABLO_T1_results.rds",
    output = file.path("result/PART4/DIABLO",
      "14.2_score_module_mediation_mixed_cohorts", outcome_tag),
    cohort_helper_file = "clean_code/14_omics_cohort_helpers.R",
    helper_file = "clean_code/14_module_mediation_helpers.R",
    entry_file = "clean_code/1_analysis/4_outcome_associations_mediation.R",
    baseline = c("age", "parity", "smk", "drk", "edu",
      "medical_insurance", "income"),
    pregnancy = c("physic_activity", "fd_grain", "fd_meat", "fd_vege",
      "fd_fruit", "fd_egg", "fd_milkpro", "antibiotic"),
    technical = c("wk", "metabolite_batch", "microbe_batch"),
    bacteria_technical = c("wk", "microbe_batch"),
    workers = as.integer(workers), nsims = as.integer(nsims),
    seed = 20261022L, q_cutoff = 0.10,
    module_status = "stable_crossomic_module",
    score_type = "full_sample_consensus_refit_not_OOF",
    fdr_denominator = "planned", weighting = "none", fpc = "none",
    adjustment = "unweighted_extended"
  )
}

# Read frozen 14.1 module scores and build the paired complete-case cohort per
# outcome. No paired score is recalculated, standardized again, or multiplied
# by the sign of its loading.
m142_prepare_paired_outcome <- function(discovery, outcome, cfg) {
  z <- discovery$outcomes[[outcome]]
  if (is.null(z)) stop("14.1 has no outcome result for ", outcome)
  m142_need(z$module_summary, c("module_id", "status"), "module_summary")
  m142_need(z$module_scores,
    c("id", "module_id", "OWO", "y", "B_score", "M_score", "score_type"),
    "module_scores")
  m142_need(z$clinical_data, c("id", "OWO", outcome), "clinical_data")

  modules <- z$module_summary[
    z$module_summary$status == cfg$module_status, , drop = FALSE
  ]
  if (!nrow(modules)) stop(outcome, " has no stable cross-omics module")
  if (anyDuplicated(modules$module_id)) stop("Duplicated stable module IDs")

  scores <- z$module_scores[z$module_scores$module_id %in% modules$module_id,
    , drop = FALSE]
  if (!nrow(scores) || !all(scores$score_type == cfg$score_type))
    stop("Unexpected or missing frozen score type for ", outcome)
  scores$id <- as.character(scores$id)
  for (module_id in modules$module_id) {
    current <- scores[scores$module_id == module_id, , drop = FALSE]
    if (!nrow(current) || anyNA(current$id) || anyDuplicated(current$id))
      stop("Missing or duplicated score IDs for ", module_id)
  }

  common_ids <- Reduce(intersect, lapply(modules$module_id, function(module_id)
    scores$id[scores$module_id == module_id]))
  clinical <- z$clinical_data
  clinical$id <- as.character(clinical$id)
  if (anyNA(clinical$id) || anyDuplicated(clinical$id))
    stop("Missing or duplicated clinical IDs for ", outcome)
  common_ids <- clinical$id[clinical$id %in% common_ids]
  if (!length(common_ids)) stop("No participants shared by all modules")
  d <- clinical[match(common_ids, clinical$id), , drop = FALSE]
  d$A <- m142_binary(d$OWO, "OWO")
  d$Y <- m142_binary(d[[outcome]], outcome)

  score_map <- data.frame(
    module_id = modules$module_id,
    B_alias = sprintf("B%03d", seq_len(nrow(modules))),
    M_alias = sprintf("M%03d", seq_len(nrow(modules))),
    stringsAsFactors = FALSE
  )
  scaling <- list()
  for (i in seq_len(nrow(score_map))) {
    current <- scores[scores$module_id == score_map$module_id[i], , drop = FALSE]
    index <- match(d$id, current$id)
    if (anyNA(index)) stop("Score alignment failed for ", score_map$module_id[i])
    if (!all(d$A == m142_binary(current$OWO[index], "score OWO")) ||
        !all(d$Y == m142_binary(current$y[index], "score outcome")))
      stop("Frozen scores and clinical labels disagree for ", score_map$module_id[i])
    d[[score_map$B_alias[i]]] <- current$B_score[index]
    d[[score_map$M_alias[i]]] <- current$M_score[index]
    scaling[[length(scaling) + 1L]] <- data.frame(
      outcome = outcome, module_id = score_map$module_id[i], block = "Bacteria",
      alias = score_map$B_alias[i], mean = mean(current$B_score[index], na.rm = TRUE),
      sd = sd(current$B_score[index], na.rm = TRUE),
      minimum = min(current$B_score[index], na.rm = TRUE),
      maximum = max(current$B_score[index], na.rm = TRUE),
      score_type = cfg$score_type, restandardized = FALSE,
      stringsAsFactors = FALSE
    )
    scaling[[length(scaling) + 1L]] <- data.frame(
      outcome = outcome, module_id = score_map$module_id[i], block = "Metabolites",
      alias = score_map$M_alias[i], mean = mean(current$M_score[index], na.rm = TRUE),
      sd = sd(current$M_score[index], na.rm = TRUE),
      minimum = min(current$M_score[index], na.rm = TRUE),
      maximum = max(current$M_score[index], na.rm = TRUE),
      score_type = cfg$score_type, restandardized = FALSE,
      stringsAsFactors = FALSE
    )
  }

  score_aliases <- c(score_map$B_alias, score_map$M_alias)
  covariates <- c(cfg$baseline, cfg$pregnancy, cfg$technical)
  m142_need(d, c(covariates, score_aliases), paste(outcome, "analysis cohort"))
  required <- c(covariates, score_aliases)
  complete <- complete.cases(d[required])
  for (column in score_aliases)
    complete <- complete & is.finite(d[[column]])
  excluded <- data.frame(
    outcome = rep(outcome, sum(!complete)), id = d$id[!complete],
    reason = rep("incomplete_extended_adjustment_or_frozen_score", sum(!complete)),
    stringsAsFactors = FALSE
  )
  d <- droplevels(d[complete, , drop = FALSE])
  d$A <- m142_binary(d$A, "complete-case OWO")
  d$Y <- m142_binary(d$Y, paste("complete-case", outcome))
  for (column in score_aliases) {
    if (!is.numeric(d[[column]]) || any(!is.finite(d[[column]])) ||
        !is.finite(sd(d[[column]])) || sd(d[[column]]) <= 0)
      stop("Non-numeric or constant frozen score: ", column)
  }

  # Expand factors once. Every model for an outcome then uses the same numeric
  # nuisance columns and exactly the same participants.
  nonconstant <- covariates[vapply(d[covariates], function(x)
    length(unique(x[!is.na(x)])) > 1L, logical(1))]
  constant_covariates <- setdiff(covariates, nonconstant)
  dropped <- data.frame(
    outcome = rep(outcome, length(constant_covariates)),
    column = constant_covariates,
    reason = rep("constant", length(constant_covariates)),
    stringsAsFactors = FALSE
  )
  matrix_covariates <- if (length(nonconstant))
    model.matrix(reformulate(nonconstant), d)[, -1, drop = FALSE] else
    matrix(nrow = nrow(d), ncol = 0L)
  base <- cbind(intercept = 1, A = d$A)
  covariate_map <- data.frame(
    outcome = character(), original_column = character(), alias = character(),
    stringsAsFactors = FALSE
  )
  for (j in seq_len(ncol(matrix_covariates))) {
    if (qr(cbind(base, matrix_covariates[, j]))$rank <= ncol(base)) {
      dropped <- rbind(dropped, data.frame(
        outcome = outcome, column = colnames(matrix_covariates)[j],
        reason = "linearly_redundant", stringsAsFactors = FALSE
      ))
    } else {
      alias <- sprintf("C%03d", nrow(covariate_map) + 1L)
      d[[alias]] <- matrix_covariates[, j]
      base <- cbind(base, matrix_covariates[, j])
      covariate_map <- rbind(covariate_map, data.frame(
        outcome = outcome, original_column = colnames(matrix_covariates)[j],
        alias = alias, stringsAsFactors = FALSE
      ))
    }
  }

  module_metadata <- modules
  module_metadata$outcome <- outcome
  keep <- c("id", "A", "Y", score_aliases, covariate_map$alias)
  list(
    outcome = outcome, data = d[keep], modules = module_metadata,
    score_map = score_map, score_scaling = m142_bind(scaling),
    covars = covariate_map$alias, covariate_map = covariate_map,
    dropped_covariates = dropped, excluded = excluded,
    n = nrow(d), events = sum(d$Y), n_OWO = sum(d$A)
  )
}

m142_encode_full_bacteria_cohort <- function(discovery, outcome, cfg,
                                              full_16s, paired) {
  z <- discovery$outcomes[[outcome]]
  d <- full_16s$data
  d$id <- as.character(d$id)
  d <- d[!is.na(d[[outcome]]) & d[[outcome]] %in% 0:1, , drop = FALSE]
  d$A <- m142_binary(d$OWO, paste(outcome, "full-16S OWO"))
  d$Y <- m142_binary(d[[outcome]], paste(outcome, "full-16S outcome"))

  audits <- list()
  for (i in seq_len(nrow(paired$score_map))) {
    module_id <- paired$score_map$module_id[i]
    alias <- paired$score_map$B_alias[i]
    loading <- z$module_loadings[
      z$module_loadings$module_id == module_id &
        z$module_loadings$block == "Bacteria", , drop = FALSE]
    d[[alias]] <- d14c_project_score(d, loading, "Bacteria")

    stored <- z$module_scores[z$module_scores$module_id == module_id,
      c("id", "B_score"), drop = FALSE]
    stored$id <- as.character(stored$id)
    overlap <- intersect(d$id, stored$id)
    projected <- d[[alias]][match(overlap, d$id)]
    reference <- stored$B_score[match(overlap, stored$id)]
    finite <- is.finite(projected) & is.finite(reference)
    max_error <- if (any(finite)) max(abs(projected[finite] - reference[finite])) else NA_real_
    audits[[length(audits) + 1L]] <- data.frame(
      outcome = outcome, module_id = module_id, block = "Bacteria",
      validation_overlap_n = sum(finite), max_abs_score_error = max_error,
      tolerance = 1e-8, passed = is.finite(max_error) && max_error < 1e-8,
      projection = "frozen_14.1_paired_DIABLO_parameters",
      stringsAsFactors = FALSE
    )
    if (!is.finite(max_error) || max_error >= 1e-8)
      stop("Frozen bacterial score reconstruction failed for ", outcome,
        " / ", module_id, "; max error=", format(max_error, digits = 6))
  }

  score_aliases <- paired$score_map$B_alias
  covariates <- c(cfg$baseline, cfg$pregnancy, cfg$bacteria_technical)
  m142_need(d, c(covariates, score_aliases), paste(outcome, "full T1 16S cohort"))
  required <- c(covariates, score_aliases)
  complete <- complete.cases(d[required])
  for (column in score_aliases) complete <- complete & is.finite(d[[column]])
  excluded <- data.frame(
    outcome = rep(outcome, sum(!complete)), id = d$id[!complete],
    cohort = rep("full_T1_16S", sum(!complete)),
    reason = rep("incomplete_bacteria_adjustment_or_projected_score", sum(!complete)),
    stringsAsFactors = FALSE
  )
  d <- droplevels(d[complete, , drop = FALSE])
  d$A <- m142_binary(d$A, paste(outcome, "complete full-16S OWO"))
  d$Y <- m142_binary(d$Y, paste(outcome, "complete full-16S outcome"))

  nonconstant <- covariates[vapply(d[covariates], function(x)
    length(unique(x[!is.na(x)])) > 1L, logical(1))]
  constant <- setdiff(covariates, nonconstant)
  dropped <- data.frame(
    outcome = rep(outcome, length(constant)), cohort = rep("full_T1_16S", length(constant)),
    column = constant, reason = rep("constant", length(constant)),
    stringsAsFactors = FALSE
  )
  mm <- if (length(nonconstant))
    model.matrix(reformulate(nonconstant), d)[, -1, drop = FALSE] else
    matrix(nrow = nrow(d), ncol = 0L)
  base <- cbind(intercept = 1, A = d$A)
  covariate_map <- data.frame(
    outcome = character(), cohort = character(), original_column = character(),
    alias = character(), stringsAsFactors = FALSE)
  for (j in seq_len(ncol(mm))) {
    if (qr(cbind(base, mm[, j]))$rank <= ncol(base)) {
      dropped <- rbind(dropped, data.frame(
        outcome = outcome, cohort = "full_T1_16S", column = colnames(mm)[j],
        reason = "linearly_redundant", stringsAsFactors = FALSE))
    } else {
      alias <- sprintf("C%03d", nrow(covariate_map) + 1L)
      d[[alias]] <- mm[, j]
      base <- cbind(base, mm[, j])
      covariate_map <- rbind(covariate_map, data.frame(
        outcome = outcome, cohort = "full_T1_16S",
        original_column = colnames(mm)[j], alias = alias,
        stringsAsFactors = FALSE))
    }
  }
  keep <- c("id", "A", "Y", score_aliases, covariate_map$alias)
  score_scaling <- m142_bind(lapply(seq_len(nrow(paired$score_map)), function(i) {
    alias <- paired$score_map$B_alias[i]
    data.frame(
      outcome = outcome, module_id = paired$score_map$module_id[i],
      block = "Bacteria", alias = alias, mean = mean(d[[alias]]),
      sd = sd(d[[alias]]), minimum = min(d[[alias]]), maximum = max(d[[alias]]),
      score_type = "frozen_14.1_score_projected_to_full_T1_16S",
      restandardized = FALSE, cohort = "full_T1_16S",
      stringsAsFactors = FALSE)
  }))
  list(
    data = d[keep], covars = covariate_map$alias,
    covariate_map = covariate_map, dropped_covariates = dropped,
    excluded = excluded, n = nrow(d), events = sum(d$Y), n_OWO = sum(d$A),
    cohort_key = "bacteria_full", analysis_cohort = "full_T1_16S",
    cohort_source = full_16s$source,
    actual_covariates = covariates,
    score_scaling = score_scaling,
    score_reconstruction_audit = m142_bind(audits)
  )
}

m142_prepare_outcome <- function(discovery, outcome, cfg, full_16s) {
  paired <- m142_prepare_paired_outcome(discovery, outcome, cfg)
  paired$analysis_cohort <- "paired_T1_16S_metabolomics"
  paired$cohort_source <- "14.1_DIABLO_paired_clinical_data"
  paired$actual_covariates <- c(cfg$baseline, cfg$pregnancy, cfg$technical)
  paired$cohort_key <- "paired"
  if (nrow(paired$score_scaling)) paired$score_scaling$cohort <- paired$analysis_cohort
  if (nrow(paired$excluded)) paired$excluded$cohort <- paired$analysis_cohort
  if (nrow(paired$covariate_map)) paired$covariate_map$cohort <- paired$analysis_cohort
  if (nrow(paired$dropped_covariates))
    paired$dropped_covariates$cohort <- paired$analysis_cohort
  bacteria_full <- m142_encode_full_bacteria_cohort(
    discovery, outcome, cfg, full_16s, paired)
  paired$cohorts <- list(paired = paired, bacteria_full = bacteria_full)
  paired$cohort_summary <- m142_bind(lapply(paired$cohorts, function(x)
    data.frame(outcome = outcome, cohort_key = x$cohort_key,
      analysis_cohort = x$analysis_cohort, cohort_source = x$cohort_source,
      n = x$n, events = x$events, n_OWO = x$n_OWO,
      covariates = paste(x$actual_covariates, collapse = ";"),
      stringsAsFactors = FALSE)))
  paired
}

m142_make_tasks <- function(obj, cfg) {
  tasks <- list()
  add <- function(module_id, analysis_type, direction, method, exposure,
                  mediator, path, cohort_key) {
    key <- list(outcome = obj$outcome, module_id = module_id,
      analysis_type = analysis_type, direction = direction,
      exposure = exposure, mediator = mediator)
    task_id <- paste0("m142_", substr(m142_hash(key), 1L, 14L))
    tasks[[task_id]] <<- list(
      task_id = task_id, outcome = obj$outcome, period = "T1",
      module_id = module_id, analysis_type = analysis_type,
      direction = direction, method = method, exposure = exposure,
      mediator = mediator, path = path, cohort_key = cohort_key,
      seed = as.integer((cfg$seed + strtoi(substr(m142_hash(task_id), 1L, 7L),
        base = 16L)) %% 2147483646L)
    )
  }
  for (i in seq_len(nrow(obj$score_map))) {
    current <- obj$score_map[i, ]
    y <- obj$outcome
    add(current$module_id, "single_score", "OWO_to_B_score", "mediate_QB",
      "A", current$B_alias, paste0("OWO -> B_score -> ", y), "bacteria_full")
    add(current$module_id, "single_score", "OWO_to_M_score", "mediate_QB",
      "A", current$M_alias, paste0("OWO -> M_score -> ", y), "paired")
    add(current$module_id, "bidirectional_score", "B_score_to_M_score",
      "mediate_QB", current$B_alias, current$M_alias,
      paste0("B_score -> M_score -> ", y), "paired")
    add(current$module_id, "bidirectional_score", "M_score_to_B_score",
      "mediate_QB", current$M_alias, current$B_alias,
      paste0("M_score -> B_score -> ", y), "paired")
  }
  tasks
}

m142_validate_model <- function(fit, binary = FALSE) {
  if (fit$rank < length(coef(fit)) || any(!is.finite(coef(fit))))
    stop("Rank-deficient regression")
  if (binary && !isTRUE(fit$converged))
    stop("Logistic regression did not converge")
  if (any(!is.finite(vcov(fit))))
    stop("Non-finite coefficient covariance")
  invisible(fit)
}

m142_mediate_effect <- function(fit, key, estimate_key = key) {
  interval <- fit[[paste0(key, ".ci")]]
  if (length(interval) != 2L) interval <- c(NA_real_, NA_real_)
  c(estimate = as.numeric(fit[[estimate_key]]), low = as.numeric(interval[1]),
    high = as.numeric(interval[2]), p = as.numeric(fit[[paste0(key, ".p")]]))
}

m142_empty_row <- function(task, obj, cfg) {
  data.frame(
    outcome = task$outcome, period = task$period, module_id = task$module_id,
    module_status = cfg$module_status, analysis_type = task$analysis_type,
    task_id = task$task_id, direction = task$direction, path = task$path,
    method = task$method, exposure = task$exposure, mediator = task$mediator,
    analysis_cohort = obj$analysis_cohort,
    cohort_source = obj$cohort_source,
    n = obj$n, events = obj$events, n_OWO = obj$n_OWO,
    exposure_contrast = if (task$exposure == "A") "OWO: 0 vs 1" else
      "frozen DIABLO score: 0 vs +1 score unit",
    nsims = cfg$nsims,
    nsims_success = NA_integer_,
    exposure_to_mediator_beta = NA_real_,
    exposure_to_mediator_low = NA_real_,
    exposure_to_mediator_high = NA_real_,
    exposure_to_mediator_p = NA_real_,
    mediator_to_outcome_logOR = NA_real_,
    mediator_to_outcome_OR = NA_real_,
    mediator_to_outcome_low = NA_real_,
    mediator_to_outcome_high = NA_real_,
    mediator_to_outcome_p = NA_real_,
    TE = NA_real_, TE_low = NA_real_, TE_high = NA_real_, TE_p = NA_real_,
    DE = NA_real_, DE_low = NA_real_, DE_high = NA_real_, DE_p = NA_real_,
    IE = NA_real_, IE_low = NA_real_, IE_high = NA_real_, IE_p = NA_real_,
    ADE = NA_real_, ADE_low = NA_real_, ADE_high = NA_real_, ADE_p = NA_real_,
    ACME = NA_real_, ACME_low = NA_real_, ACME_high = NA_real_, ACME_p = NA_real_,
    PM = NA_real_, PM_low = NA_real_, PM_high = NA_real_, PM_p = NA_real_,
    IE_q = NA_real_, ACME_q = NA_real_, PM_q = NA_real_,
    IE_q10 = FALSE, PM_positive_q10 = FALSE, PM_interpretable = FALSE,
    effect_scale = "", adjustment = cfg$adjustment,
    OWO_role = if (task$analysis_type == "bidirectional_score")
      "adjustment_covariate_in_both_models" else "exposure",
    OWO_adjusted = task$analysis_type == "bidirectional_score",
    weighting_used = FALSE, FPC_used = FALSE,
    covariates = paste(obj$actual_covariates, collapse = ";"),
    interpretation_flags = "", status = "not_run", warnings = "", error = "",
    seed = task$seed, stringsAsFactors = FALSE
  )
}

m142_base_flags <- function(task) {
  unique(c(
    "conditional_on_outcome_supervised_DIABLO_selection",
    "full_sample_frozen_scores_not_out_of_fold",
    if (task$cohort_key == "bacteria_full") c(
      "unweighted_full_T1_16S_cohort",
      "frozen_B_score_projected_from_paired_DIABLO") else
      "unweighted_paired_16S_metabolomics_cohort",
    "pregnancy_covariates_may_be_post_exposure",
    "T1_temporality_requires_biospecimen_before_outcome_onset",
    if (task$analysis_type == "bidirectional_score")
      "same_timepoint_score_order_not_identified"
  ))
}

m142_fit_mediate <- function(task, obj, cfg) {
  d <- obj$data
  covars <- c(if (task$exposure != "A") "A", obj$covars)
  mediator_formula <- reformulate(c(task$exposure, covars), task$mediator)
  outcome_formula <- reformulate(c(task$exposure, task$mediator, covars), "Y")
  mediator_model <- lm(mediator_formula, data = d, na.action = na.fail)
  outcome_model <- glm(outcome_formula, data = d, family = binomial(),
    na.action = na.fail, control = glm.control(maxit = 100L))
  mediator_model$call$formula <- mediator_formula
  outcome_model$call$formula <- outcome_formula
  m142_validate_model(mediator_model)
  m142_validate_model(outcome_model, binary = TRUE)

  fit <- suppressMessages(mediation::mediate(
    mediator_model, outcome_model, treat = task$exposure,
    mediator = task$mediator, control.value = 0, treat.value = 1,
    boot = FALSE, sims = cfg$nsims, robustSE = FALSE, long = TRUE
  ))
  effects <- list(
    TE = m142_mediate_effect(fit, "tau", "tau.coef"),
    ADE = m142_mediate_effect(fit, "z.avg", "z.avg"),
    ACME = m142_mediate_effect(fit, "d.avg", "d.avg"),
    PM = m142_mediate_effect(fit, "n.avg", "n.avg")
  )
  effects$DE <- effects$ADE
  effects$IE <- effects$ACME
  valid <- is.finite(fit$tau.sims) & is.finite(fit$z.avg.sims) &
    is.finite(fit$d.avg.sims)
  if (sum(valid) != cfg$nsims) stop("Incomplete quasi-Bayesian simulations")

  a <- coef(summary(mediator_model))[task$exposure, ]
  b <- coef(summary(outcome_model))[task$mediator, ]
  paths <- data.frame(
    path = c("exposure_to_mediator", "mediator_to_outcome"),
    estimate = c(a[1], b[1]), se = c(a[2], b[2]),
    low = c(a[1] - qnorm(.975) * a[2], b[1] - qnorm(.975) * b[2]),
    high = c(a[1] + qnorm(.975) * a[2], b[1] + qnorm(.975) * b[2]),
    p = c(a[4], b[4]),
    scale = c("linear_beta", "log_odds"), stringsAsFactors = FALSE
  )
  flags <- m142_base_flags(task)
  if (any(abs(coef(outcome_model)) > 20) ||
      any(fitted(outcome_model) < 1e-8 | fitted(outcome_model) > 1 - 1e-8))
    flags <- c(flags, "possible_separation_or_extreme_prediction")
  if (min(sum(d$Y == 1), sum(d$Y == 0)) / outcome_model$rank < 10)
    flags <- c(flags, "sparse_events_per_parameter")
  list(
    effects = effects, fit = fit, mediator_model = mediator_model,
    outcome_model = outcome_model, paths = paths,
    nsims_success = sum(valid), flags = unique(flags),
    effect_scale = "marginal risk difference; PM=mediate average proportion"
  )
}

m142_add_effect <- function(row, effects, name) {
  if (is.null(effects[[name]])) return(row)
  current <- effects[[name]]
  row[[name]] <- unname(current["estimate"])
  row[[paste0(name, "_low")]] <- unname(current["low"])
  row[[paste0(name, "_high")]] <- unname(current["high"])
  row[[paste0(name, "_p")]] <- unname(current["p"])
  row
}

m142_fit_task <- function(task, obj, cfg) {
  analysis_obj <- obj$cohorts[[task$cohort_key]]
  if (is.null(analysis_obj)) stop("Unknown analysis cohort: ", task$cohort_key)
  row <- m142_empty_row(task, analysis_obj, cfg)
  warnings <- character()
  result <- NULL
  tryCatch(withCallingHandlers({
    set.seed(task$seed)
    fitted <- m142_fit_mediate(task, analysis_obj, cfg)
    for (effect in c("TE", "DE", "IE", "ADE", "ACME", "PM"))
      row <- m142_add_effect(row, fitted$effects, effect)
    a <- fitted$paths[fitted$paths$path == "exposure_to_mediator", ]
    b <- fitted$paths[fitted$paths$path == "mediator_to_outcome", ]
    row$exposure_to_mediator_beta <- a$estimate
    row$exposure_to_mediator_low <- a$low
    row$exposure_to_mediator_high <- a$high
    row$exposure_to_mediator_p <- a$p
    row$mediator_to_outcome_logOR <- b$estimate
    row$mediator_to_outcome_OR <- exp(b$estimate)
    row$mediator_to_outcome_low <- exp(b$low)
    row$mediator_to_outcome_high <- exp(b$high)
    row$mediator_to_outcome_p <- b$p
    row$nsims_success <- fitted$nsims_success
    flags <- fitted$flags
    if (is.finite(row$TE_low) && is.finite(row$TE_high) &&
        row$TE_low <= 0 && row$TE_high >= 0)
      flags <- c(flags, "TE_CI_crosses_zero_PM_unstable")
    if (is.finite(row$PM) && (row$PM < 0 || row$PM > 1))
      flags <- c(flags, "PM_outside_0_1")
    if (is.finite(row$PM_low) && is.finite(row$PM_high) &&
        (row$PM_low < 0 || row$PM_high > 1))
      flags <- c(flags, "PM_CI_outside_0_1")
    if (is.finite(row$IE) && is.finite(row$DE) && row$IE * row$DE < 0)
      flags <- c(flags, "opposing_direct_indirect_effects")
    if (is.finite(row$IE) && is.finite(row$TE) && row$IE * row$TE < 0)
      flags <- c(flags, "IE_TE_direction_discordant")
    row$PM_interpretable <- is.finite(row$PM) && row$PM >= 0 && row$PM <= 1 &&
      is.finite(row$TE_low) && is.finite(row$TE_high) &&
      !(row$TE_low <= 0 && row$TE_high >= 0) &&
      is.finite(row$IE) && row$IE * row$TE > 0
    row$effect_scale <- fitted$effect_scale
    row$interpretation_flags <- paste(unique(flags), collapse = ";")
    row$status <- "ok"
    result <- list(row = row, fitted = fitted)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    row$status <<- "failed"
    row$error <<- conditionMessage(e)
    row$interpretation_flags <<- paste(c(m142_base_flags(task), "model_failed"),
      collapse = ";")
    NULL
  })
  row$warnings <- paste(unique(warnings), collapse = " | ")
  if (length(warnings))
    row$interpretation_flags <- paste(row$interpretation_flags,
      "model_warnings_present", sep = ";")
  if (is.null(result)) result <- list(row = row, fitted = NULL)
  result$row <- row
  result
}

m142_run_task <- function(task, obj, cfg) {
  do.call(RNGkind, as.list(c("L'Ecuyer-CMRG", "Inversion", "Rejection")))
  m142_fit_task(task, obj, cfg)
}

m142_adjust_fdr <- function(rows, cfg) {
  if (!nrow(rows)) return(list(results = rows, families = data.frame()))
  keys <- c("outcome", "analysis_type")
  rows$IE_q <- rows$ACME_q <- rows$PM_q <- NA_real_
  families <- list()
  groups <- interaction(rows[keys], drop = TRUE, lex.order = TRUE)
  for (index in split(seq_len(nrow(rows)), groups)) {
    for (effect in c("IE", "PM")) {
      p <- rows[[paste0(effect, "_p")]][index]
      good <- rows$status[index] == "ok" & is.finite(p) & p >= 0 & p <= 1
      denominator <- if (cfg$fdr_denominator == "planned") length(index) else sum(good)
      if (any(good)) rows[[paste0(effect, "_q")]][index[good]] <-
        p.adjust(p[good], method = "BH", n = denominator)
      families[[length(families) + 1L]] <- cbind(
        rows[index[1], keys, drop = FALSE],
        data.frame(effect = effect, planned = length(index),
          evaluable = sum(good), fdr_n = denominator, method = "BH",
          stringsAsFactors = FALSE)
      )
    }
  }
  mediate_rows <- rows$method == "mediate_QB"
  rows$ACME_q[mediate_rows] <- rows$IE_q[mediate_rows]
  rows$IE_q10 <- rows$status == "ok" & is.finite(rows$IE_q) &
    rows$IE_q < cfg$q_cutoff
  # Main-figure selection uses ACME q, not PM q.
  rows$PM_positive_q10 <- rows$status == "ok" & is.finite(rows$PM) &
    rows$PM > 0 & is.finite(rows$ACME_q) & rows$ACME_q < cfg$q_cutoff
  list(results = rows, families = m142_bind(families))
}

m142_ci <- function(estimate, low, high, multiplier = 1) {
  ifelse(is.finite(estimate), sprintf("%.3f [%.3f, %.3f]",
    multiplier * estimate, multiplier * low, multiplier * high), NA_character_)
}

m142_compact <- function(rows) {
  if (!nrow(rows)) return(rows)
  data.frame(
    outcome = rows$outcome, period = rows$period, module = rows$module_id,
    analysis = rows$analysis_type, direction = rows$direction,
    path = rows$path, method = rows$method,
    analysis_cohort = rows$analysis_cohort,
    N_events = paste(rows$n, rows$events, sep = "/"),
    contrast = rows$exposure_contrast,
    PM_percent_CI95 = m142_ci(rows$PM, rows$PM_low, rows$PM_high, 100),
    PM_P = rows$PM_p, PM_q = rows$PM_q,
    PM_positive_q10 = rows$PM_positive_q10,
    IE_CI95 = m142_ci(rows$IE, rows$IE_low, rows$IE_high),
    IE_P = rows$IE_p, IE_q = rows$IE_q,
    ACME_CI95 = m142_ci(rows$ACME, rows$ACME_low, rows$ACME_high),
    ACME_P = rows$ACME_p, ACME_q = rows$ACME_q,
    TE_CI95 = m142_ci(rows$TE, rows$TE_low, rows$TE_high), TE_P = rows$TE_p,
    ADE_CI95 = m142_ci(rows$ADE, rows$ADE_low, rows$ADE_high), ADE_P = rows$ADE_p,
    DE_CI95 = m142_ci(rows$DE, rows$DE_low, rows$DE_high), DE_P = rows$DE_p,
    PM_interpretable = rows$PM_interpretable,
    effect_scale = rows$effect_scale,
    interpretation_flags = rows$interpretation_flags,
    status = rows$status, error = rows$error,
    stringsAsFactors = FALSE
  )
}

m142_export <- function(rows, families, objects, tasks, cfg, results) {
  compact <- m142_compact(rows)
  group_keys <- c("outcome", "analysis_type")
  groups <- unique(rows[group_keys])
  counts <- m142_bind(lapply(seq_len(nrow(groups)), function(i) {
    keep <- rows$outcome == groups$outcome[i] &
      rows$analysis_type == groups$analysis_type[i]
    current <- rows[keep, , drop = FALSE]
    cbind(groups[i, , drop = FALSE], data.frame(
      planned = nrow(current), successful = sum(current$status == "ok"),
      PM_positive_q10 = sum(current$PM_positive_q10, na.rm = TRUE),
      IE_q10 = sum(current$IE_q10, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }))
  by_object <- function(field) m142_bind(lapply(objects, function(x) {
    current <- x[[field]]
    if (is.null(current) || !is.data.frame(current) || !nrow(current))
      return(data.frame())
    if (!"outcome" %in% names(current)) current$outcome <- x$outcome
    current
  }))
  by_cohort <- function(field) m142_bind(lapply(objects, function(x)
    m142_bind(lapply(x$cohorts, function(z) {
      current <- z[[field]]
      if (is.null(current) || !is.data.frame(current) || !nrow(current))
        return(data.frame())
      if (!"outcome" %in% names(current)) current$outcome <- x$outcome
      current
    }))))
  mediation_paths <- m142_bind(lapply(results, function(result) {
    if (is.null(result$fitted) || is.null(result$fitted$paths)) return(data.frame())
    current <- result$fitted$paths
    current$outcome <- result$row$outcome
    current$module_id <- result$row$module_id
    current$analysis_type <- result$row$analysis_type
    current$direction <- result$row$direction
    current
  }))
  tables <- list(
    summary = compact,
    PM_candidates = compact[rows$PM_positive_q10, , drop = FALSE],
    all_results = rows, FDR_families = families, counts = counts,
    mediation_paths = mediation_paths,
    score_scaling = by_cohort("score_scaling"),
    cohort_summary = by_object("cohort_summary"),
    score_reconstruction_audit = m142_bind(lapply(objects, function(x)
      x$cohorts$bacteria_full$score_reconstruction_audit)),
    module_summary = by_object("modules"),
    covariate_map = by_cohort("covariate_map"),
    dropped_covariates = by_cohort("dropped_covariates"),
    excluded_participants = by_cohort("excluded"),
    model_failures = rows[rows$status != "ok", , drop = FALSE]
  )
  dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) utils::write.csv(
    tables[[name]], file.path(cfg$output, paste0(name, ".csv")),
    row.names = FALSE, fileEncoding = "UTF-8", na = ""
  )
  writexl::write_xlsx(tables, file.path(cfg$output, "mediation_summary.xlsx"))
  m142_save_rds(list(
    config = cfg, tables = tables,
    prepared_inputs = objects, tasks = tasks, results = results,
    input_md5 = tools::md5sum(cfg$input),
    code_md5 = tools::md5sum(unique(c(cfg$cohort_helper_file,
      cfg$helper_file, cfg$entry_file))),
    full_16S_input_md5 = tools::md5sum(unique(unlist(lapply(objects,
      function(x) x$full_16S_files)))),
    session_info = sessionInfo(), completed_at = Sys.time()
  ), file.path(cfg$output, "mediation_results.rds"))
  writeLines(c(
    paste("Version:", cfg$version),
    paste("Outcomes:", paste(cfg$outcomes, collapse = ", ")),
    paste("Planned / successful:", nrow(rows), "/", sum(rows$status == "ok")),
    paste("PM>0 and ACME q<0.10:", sum(rows$PM_positive_q10, na.rm = TRUE)),
    paste("mediate simulations per applicable model:", cfg$nsims),
    "Blocks: single score (Figure 5c); bidirectional scores (Figure 5d)",
    paste("Output:", normalizePath(cfg$output, winslash = "/", mustWork = FALSE))
  ), file.path(cfg$output, "run_summary.txt"), useBytes = TRUE)
  invisible(tables)
}

m142_main <- function(outcomes = c("GDM", "HDP"), workers = 8L,
                      nsims = 5000L) {
  cfg <- m142_config(outcomes, workers, nsims)
  if (!file.exists(cfg$input)) stop("Missing DIABLO result: ", cfg$input)
  for (file in unique(c(cfg$cohort_helper_file, cfg$helper_file, cfg$entry_file)))
    if (!file.exists(file)) stop("Missing analysis code: ", file)
  packages <- c("mediation", "writexl")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  if (!is.finite(cfg$workers) || cfg$workers < 1L ||
      cfg$workers != floor(cfg$workers))
    stop("workers must be a positive integer")
  if (!is.finite(cfg$nsims) || cfg$nsims < 100L || cfg$nsims != floor(cfg$nsims))
    stop("nsims must be an integer >= 100")

  source(cfg$cohort_helper_file, local = environment(m142_main))
  discovery <- readRDS(cfg$input)
  full_16s <- d14c_load_full_t1_16s(cfg$outcomes, cfg$baseline, cfg$pregnancy)
  objects <- setNames(lapply(cfg$outcomes, function(outcome)
    m142_prepare_outcome(discovery, outcome, cfg, full_16s)), cfg$outcomes)
  for (outcome in names(objects)) objects[[outcome]]$full_16S_files <- full_16s$files
  tasks <- setNames(lapply(objects, m142_make_tasks, cfg = cfg), cfg$outcomes)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1")
  total_tasks <- sum(lengths(tasks))
  n_workers <- min(cfg$workers, max(1L, total_tasks))
  cluster <- NULL
  if (n_workers > 1L) {
    cluster <- parallel::makePSOCKcluster(n_workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(paths, wd) {
      .libPaths(paths)
      setwd(wd)
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1")
      NULL
    }, .libPaths(), getwd())
    parallel::clusterExport(cluster,
      ls(environment(m142_main), pattern = "^m142_"),
      envir = environment(m142_main))
  }

  results <- list()
  for (outcome in cfg$outcomes) {
    message(outcome, ": ", length(tasks[[outcome]]),
      " planned score-level models; workers=", n_workers)
    fitted <- if (is.null(cluster)) {
      lapply(tasks[[outcome]], m142_run_task, obj = objects[[outcome]],
        cfg = cfg)
    } else {
      parallel::parLapplyLB(cluster, tasks[[outcome]], m142_run_task,
        obj = objects[[outcome]], cfg = cfg)
    }
    results <- c(results, fitted)
  }
  rows <- m142_bind(lapply(results, `[[`, "row"))
  adjusted <- m142_adjust_fdr(rows, cfg)
  rows <- adjusted$results
  stopifnot(nrow(rows) == total_tasks, !anyDuplicated(rows$task_id))
  tables <- m142_export(rows, adjusted$families, objects, tasks, cfg, results)
  invisible(list(config = cfg, tables = tables, results = results,
    prepared_inputs = objects, tasks = tasks))
}
