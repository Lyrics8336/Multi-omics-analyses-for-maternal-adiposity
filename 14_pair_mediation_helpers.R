# Shared implementation for 14.3 T1 member-level mediation.

if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  invisible(suppressWarnings(try(
    Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE
  )))
}

m14_bind <- function(xs) {
  xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
  if (!length(xs)) return(data.frame())
  cols <- unique(unlist(lapply(xs, names)))
  xs <- lapply(xs, function(x) {
    for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
    x[cols]
  })
  ans <- do.call(rbind, xs)
  rownames(ans) <- NULL
  ans
}

m14_hash <- function(x) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, version = 2L)
  unname(tools::md5sum(f))
}

m14_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp, compress = TRUE)
  if (file.exists(path) && !file.remove(path)) stop("Cannot replace: ", path)
  if (!file.rename(tmp, path)) stop("Cannot save: ", path)
}

m14_need <- function(x, columns, label) {
  absent <- setdiff(columns, names(x))
  if (length(absent)) stop(label, " is missing: ", paste(absent, collapse = ", "))
}

m14_binary <- function(x, label) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) x <- as.integer(x)
  x <- suppressWarnings(as.numeric(x))
  if (anyNA(x) || !all(x %in% 0:1) || length(unique(x)) != 2L)
    stop(label, " must contain both explicitly coded values 0 and 1")
  x
}

m14_config <- function(analysis = c("pair", "single"), outcomes = c("GDM", "HDP"),
                       workers = 8L, nsims = 5000L) {
  analysis <- match.arg(analysis)
  outcomes <- unique(as.character(outcomes))
  if (!length(outcomes) || any(!outcomes %in% c("GDM", "HDP")))
    stop("outcomes must be GDM, HDP, or both")
  list(
    version = "14.3_v3.0_mixed_cohorts", analysis = analysis, outcomes = outcomes,
    input = "result/PART4/DIABLO/14.1_DIABLO_T1_results.rds",
    output = file.path("result/PART4/DIABLO", if (analysis == "pair")
      "14.3_bidir_feature_mediation_mixed_cohorts" else
      "14.3_OWO_feature_mediation_mixed_cohorts",
      paste(outcomes, collapse = "_")),
    cohort_helper_file = "clean_code/14_omics_cohort_helpers.R",
    entry_file = "clean_code/1_analysis/4_outcome_associations_mediation.R",
    helper_file = "clean_code/14_pair_mediation_helpers.R",
    baseline = c("age", "parity", "smk", "drk", "edu",
      "medical_insurance", "income"),
    pregnancy = c("physic_activity", "fd_grain", "fd_meat", "fd_vege",
      "fd_fruit", "fd_egg", "fd_milkpro", "antibiotic"),
    technical = c("wk", "metabolite_batch", "microbe_batch"),
    bacteria_technical = c("wk", "microbe_batch"),
    nsims = as.integer(nsims), workers = as.integer(workers),
    seed = 20261021L, q_cutoff = 0.10,
    fdr_denominator = "planned",
    weighting = "none", fpc = "none",
    module_status = "stable_crossomic_module",
    effect_definition = "mediate average effects on marginal risk-difference scale"
  )
}

# Read only the frozen 14.1 objects. No DIABLO model is refitted.
m14_prepare_paired_outcome <- function(discovery, outcome, cfg) {
  z <- discovery$outcomes[[outcome]]
  if (is.null(z)) stop("14.1 has no outcome result for ", outcome)
  m14_need(z$module_summary, c("module_id", "status"), "module_summary")
  m14_need(z$stable_module_features,
    c("module_id", "block", "feature_id", "feature_name", "stable_feature"),
    "stable_module_features")
  m14_need(z$module_scores, c("id", "module_id"), "module_scores")
  m14_need(z$clinical_data, c("id", "OWO", outcome), "clinical_data")

  modules <- z$module_summary[
    z$module_summary$status == cfg$module_status, , drop = FALSE
  ]
  if (!nrow(modules)) stop(outcome, " has no stable cross-omics module")
  members <- z$stable_module_features[
    z$stable_module_features$module_id %in% modules$module_id &
      z$stable_module_features$stable_feature %in% TRUE, , drop = FALSE
  ]
  if (!nrow(members) || !all(members$block %in% c("Bacteria", "Metabolites")))
    stop("Invalid stable members for ", outcome)
  for (module in modules$module_id) {
    if (!all(c("Bacteria", "Metabolites") %in%
             members$block[members$module_id == module]))
      stop("Module lacks one omics block: ", module)
  }
  if (anyDuplicated(members[c("module_id", "block", "feature_id")]))
    stop("Duplicated stable member keys for ", outcome)

  features <- unique(members[c("block", "feature_id", "feature_name")])
  if (anyDuplicated(features[c("block", "feature_id")]))
    stop("Conflicting stable feature annotations for ", outcome)
  features <- features[order(features$block, features$feature_id), , drop = FALSE]
  features$alias <- sprintf("F%03d", seq_len(nrow(features)))

  ids <- unique(as.character(z$module_scores$id[
    z$module_scores$module_id %in% modules$module_id
  ]))
  d <- z$clinical_data
  d$id <- as.character(d$id)
  if (anyNA(d$id) || anyDuplicated(d$id)) stop("Invalid clinical IDs for ", outcome)
  d <- d[match(ids, d$id), , drop = FALSE]
  if (anyNA(d$id)) stop("A module-score participant is absent from clinical_data")

  covariates <- c(cfg$baseline, cfg$pregnancy, cfg$technical)
  needed <- unique(c("id", "OWO", outcome, covariates, features$feature_id))
  m14_need(d, needed, paste(outcome, "analysis cohort"))
  d <- d[needed]
  complete <- complete.cases(d)
  for (nm in features$feature_id) complete <- complete & is.finite(d[[nm]])
  excluded_ids <- d$id[!complete]
  excluded <- data.frame(
    outcome = rep(outcome, length(excluded_ids)), id = excluded_ids,
    reason = rep("incomplete_extended_adjustment_or_stable_feature",
                 length(excluded_ids)),
    stringsAsFactors = FALSE
  )
  d <- droplevels(d[complete, , drop = FALSE])
  d$A <- m14_binary(d$OWO, "OWO")
  d$Y <- m14_binary(d[[outcome]], outcome)

  scalers <- vector("list", nrow(features))
  for (i in seq_len(nrow(features))) {
    raw <- d[[features$feature_id[i]]]
    if (!is.numeric(raw) || any(!is.finite(raw)) || !is.finite(sd(raw)) || sd(raw) <= 0)
      stop("Non-numeric or constant stable feature: ", features$feature_id[i])
    center <- mean(raw)
    scale <- sd(raw)
    d[[features$alias[i]]] <- (raw - center) / scale
    scalers[[i]] <- data.frame(
      outcome = outcome, block = features$block[i],
      feature_id = features$feature_id[i], feature_name = features$feature_name[i],
      alias = features$alias[i], center = center, scale = scale,
      n = nrow(d), stringsAsFactors = FALSE
    )
  }

  nonconstant <- covariates[vapply(d[covariates], function(x)
    length(unique(x[!is.na(x)])) > 1L, logical(1))]
  constant_covariates <- setdiff(covariates, nonconstant)
  dropped <- data.frame(
    outcome = rep(outcome, length(constant_covariates)),
    column = constant_covariates,
    reason = rep("constant", length(constant_covariates)),
    stringsAsFactors = FALSE
  )
  mm <- if (length(nonconstant))
    model.matrix(reformulate(nonconstant), d)[, -1, drop = FALSE] else
    matrix(nrow = nrow(d), ncol = 0L)
  base <- cbind(intercept = 1, A = d$A)
  covariate_map <- data.frame(
    outcome = character(), original_column = character(), alias = character(),
    stringsAsFactors = FALSE
  )
  for (j in seq_len(ncol(mm))) {
    if (qr(cbind(base, mm[, j]))$rank <= ncol(base)) {
      dropped <- rbind(dropped, data.frame(
        outcome = outcome, column = colnames(mm)[j],
        reason = "linearly_redundant", stringsAsFactors = FALSE
      ))
    } else {
      alias <- sprintf("C%03d", nrow(covariate_map) + 1L)
      d[[alias]] <- mm[, j]
      base <- cbind(base, mm[, j])
      covariate_map <- rbind(covariate_map, data.frame(
        outcome = outcome, original_column = colnames(mm)[j], alias = alias,
        stringsAsFactors = FALSE
      ))
    }
  }

  members$alias <- features$alias[
    match(paste(members$block, members$feature_id),
          paste(features$block, features$feature_id))
  ]
  keep <- c("id", "A", "Y", features$alias, covariate_map$alias)
  list(
    outcome = outcome, data = d[keep], modules = modules, members = members,
    features = features, scalers = m14_bind(scalers), covars = covariate_map$alias,
    covariate_map = covariate_map, dropped_covariates = dropped,
    excluded = excluded, n = nrow(d), events = sum(d$Y), n_OWO = sum(d$A)
  )
}

m14_prepare_full_bacteria_cohort <- function(outcome, cfg, full_16s, paired) {
  bacteria <- paired$features[paired$features$block == "Bacteria", , drop = FALSE]
  if (!nrow(bacteria)) stop("No stable bacterial features for ", outcome)
  d <- full_16s$data
  d$id <- as.character(d$id)
  d <- d[!is.na(d[[outcome]]) & d[[outcome]] %in% 0:1, , drop = FALSE]
  covariates <- c(cfg$baseline, cfg$pregnancy, cfg$bacteria_technical)
  m14_need(d, c("id", "OWO", outcome, covariates, bacteria$feature_id),
    paste(outcome, "full T1 16S cohort"))
  required <- c(covariates, bacteria$feature_id)
  complete <- complete.cases(d[required])
  for (feature in bacteria$feature_id)
    complete <- complete & is.finite(d[[feature]])
  excluded <- data.frame(
    outcome = rep(outcome, sum(!complete)), id = d$id[!complete],
    cohort = rep("full_T1_16S", sum(!complete)),
    reason = rep("incomplete_bacteria_adjustment_or_stable_feature", sum(!complete)),
    stringsAsFactors = FALSE)
  d <- droplevels(d[complete, , drop = FALSE])
  d$A <- m14_binary(d$OWO, paste(outcome, "full-16S OWO"))
  d$Y <- m14_binary(d[[outcome]], paste(outcome, "full-16S outcome"))

  scalers <- list()
  for (i in seq_len(nrow(bacteria))) {
    raw <- d[[bacteria$feature_id[i]]]
    if (!is.numeric(raw) || !is.finite(sd(raw)) || sd(raw) <= 0)
      stop("Non-numeric or constant full-cohort bacterial feature: ",
        bacteria$feature_id[i])
    center <- mean(raw)
    scale <- sd(raw)
    d[[bacteria$alias[i]]] <- (raw - center) / scale
    scalers[[length(scalers) + 1L]] <- data.frame(
      outcome = outcome, cohort = "full_T1_16S", block = "Bacteria",
      feature_id = bacteria$feature_id[i], feature_name = bacteria$feature_name[i],
      alias = bacteria$alias[i], center = center, scale = scale, n = nrow(d),
      stringsAsFactors = FALSE)
  }

  nonconstant <- covariates[vapply(d[covariates], function(x)
    length(unique(x[!is.na(x)])) > 1L, logical(1))]
  constant <- setdiff(covariates, nonconstant)
  dropped <- data.frame(
    outcome = rep(outcome, length(constant)), cohort = rep("full_T1_16S", length(constant)),
    column = constant, reason = rep("constant", length(constant)),
    stringsAsFactors = FALSE)
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
  keep <- c("id", "A", "Y", bacteria$alias, covariate_map$alias)
  list(
    data = d[keep], covars = covariate_map$alias,
    covariate_map = covariate_map, dropped_covariates = dropped,
    excluded = excluded, scalers = m14_bind(scalers),
    n = nrow(d), events = sum(d$Y), n_OWO = sum(d$A),
    cohort_key = "bacteria_full", analysis_cohort = "full_T1_16S",
    cohort_source = full_16s$source,
    actual_covariates = covariates)
}

m14_prepare_outcome <- function(discovery, outcome, cfg, full_16s) {
  paired <- m14_prepare_paired_outcome(discovery, outcome, cfg)
  paired$cohort_key <- "paired"
  paired$analysis_cohort <- "paired_T1_16S_metabolomics"
  paired$cohort_source <- "14.1_DIABLO_paired_clinical_data"
  paired$actual_covariates <- c(cfg$baseline, cfg$pregnancy, cfg$technical)
  if (nrow(paired$scalers)) paired$scalers$cohort <- paired$analysis_cohort
  if (nrow(paired$excluded)) paired$excluded$cohort <- paired$analysis_cohort
  if (nrow(paired$covariate_map)) paired$covariate_map$cohort <- paired$analysis_cohort
  if (nrow(paired$dropped_covariates))
    paired$dropped_covariates$cohort <- paired$analysis_cohort

  cohorts <- list(paired = paired)
  if (cfg$analysis == "single")
    cohorts$bacteria_full <- m14_prepare_full_bacteria_cohort(
      outcome, cfg, full_16s, paired)
  paired$cohorts <- cohorts
  paired$cohort_summary <- m14_bind(lapply(cohorts, function(x)
    data.frame(outcome = outcome, cohort_key = x$cohort_key,
      analysis_cohort = x$analysis_cohort, cohort_source = x$cohort_source,
      n = x$n, events = x$events, n_OWO = x$n_OWO,
      covariates = paste(x$actual_covariates, collapse = ";"),
      stringsAsFactors = FALSE)))
  if (cfg$analysis == "single" && !is.null(cohorts$bacteria_full)) {
    paired$scalers <- m14_bind(list(
      cohorts$bacteria_full$scalers,
      paired$scalers[paired$scalers$block == "Metabolites", , drop = FALSE]))
  }
  paired
}

m14_modules_for <- function(members, aliases) {
  z <- members[members$alias %in% aliases, c("alias", "module_id"), drop = FALSE]
  modules <- Reduce(intersect, lapply(aliases, function(a) z$module_id[z$alias == a]))
  paste(sort(unique(modules)), collapse = ";")
}

m14_feature_name <- function(obj, alias) {
  obj$features$feature_name[match(alias, obj$features$alias)]
}

m14_feature_id <- function(obj, alias) {
  obj$features$feature_id[match(alias, obj$features$alias)]
}

m14_make_tasks <- function(obj, analysis, cfg) {
  tasks <- list()
  add <- function(meta) {
    if (is.null(meta$cohort_key)) meta$cohort_key <- if (analysis == "pair")
      "paired" else if (meta$mediator_type == "Bacteria")
      "bacteria_full" else "paired"
    scientific_key <- meta[setdiff(names(meta), c("task_id", "seed"))]
    id <- paste0(if (analysis == "pair") "pair_" else "single_",
                 substr(m14_hash(scientific_key), 1, 12))
    if (id %in% names(tasks)) return(invisible(NULL))
    meta$task_id <- id
    meta$seed <- as.integer((cfg$seed + strtoi(substr(m14_hash(id), 1, 7), 16L)) %%
                              2147483646L)
    tasks[[id]] <<- meta
  }

  if (analysis == "pair") {
    for (module in obj$modules$module_id) {
      z <- obj$members[obj$members$module_id == module, , drop = FALSE]
      bacteria <- sort(unique(z$alias[z$block == "Bacteria"]))
      metabolites <- sort(unique(z$alias[z$block == "Metabolites"]))
      for (b in bacteria) for (m in metabolites) {
        common_modules <- m14_modules_for(obj$members, c(b, m))
        pair_id <- paste(obj$outcome, m14_feature_id(obj, b),
                         m14_feature_id(obj, m), sep = "__")
        add(list(
          analysis_type = "feature_bidirectional", outcome = obj$outcome,
          period = "T1", module_ids = common_modules, pair_id = pair_id,
          direction = "B_to_M", exposure = b, mediator = m,
          exposure_type = "Bacteria", mediator_type = "Metabolites",
          exposure_id = m14_feature_id(obj, b), exposure_name = m14_feature_name(obj, b),
          mediator_id = m14_feature_id(obj, m), mediator_name = m14_feature_name(obj, m),
          bacteria_id = m14_feature_id(obj, b), bacteria_name = m14_feature_name(obj, b),
          metabolite_id = m14_feature_id(obj, m), metabolite_name = m14_feature_name(obj, m),
          path = paste0(m14_feature_name(obj, b), " -> ", m14_feature_name(obj, m),
                        " -> ", obj$outcome)
        ))
        add(list(
          analysis_type = "feature_bidirectional", outcome = obj$outcome,
          period = "T1", module_ids = common_modules, pair_id = pair_id,
          direction = "M_to_B", exposure = m, mediator = b,
          exposure_type = "Metabolites", mediator_type = "Bacteria",
          exposure_id = m14_feature_id(obj, m), exposure_name = m14_feature_name(obj, m),
          mediator_id = m14_feature_id(obj, b), mediator_name = m14_feature_name(obj, b),
          bacteria_id = m14_feature_id(obj, b), bacteria_name = m14_feature_name(obj, b),
          metabolite_id = m14_feature_id(obj, m), metabolite_name = m14_feature_name(obj, m),
          path = paste0(m14_feature_name(obj, m), " -> ", m14_feature_name(obj, b),
                        " -> ", obj$outcome)
        ))
      }
    }
  } else {
    for (i in seq_len(nrow(obj$features))) {
      f <- obj$features[i, ]
      modules <- paste(sort(unique(obj$members$module_id[
        obj$members$alias == f$alias])), collapse = ";")
      add(list(
        analysis_type = "OWO_feature_single", outcome = obj$outcome,
        period = "T1", module_ids = modules, pair_id = NA_character_,
        direction = "OWO_to_feature", exposure = "A", mediator = f$alias,
        exposure_type = "OWO", mediator_type = f$block,
        exposure_id = "OWO", exposure_name = "OWO",
        mediator_id = f$feature_id, mediator_name = f$feature_name,
        bacteria_id = if (f$block == "Bacteria") f$feature_id else NA_character_,
        bacteria_name = if (f$block == "Bacteria") f$feature_name else NA_character_,
        metabolite_id = if (f$block == "Metabolites") f$feature_id else NA_character_,
        metabolite_name = if (f$block == "Metabolites") f$feature_name else NA_character_,
        path = paste0("OWO -> ", f$feature_name, " -> ", obj$outcome)
      ))
    }
  }
  tasks
}

m14_validate_model <- function(fit, binary = FALSE) {
  if (fit$rank < length(coef(fit)) || any(!is.finite(coef(fit))))
    stop("Rank-deficient regression")
  if (binary && !isTRUE(fit$converged)) stop("Logistic regression did not converge")
  if (any(!is.finite(vcov(fit)))) stop("Non-finite coefficient covariance")
  invisible(fit)
}

m14_effect <- function(fit, key, estimate_key = key) {
  ci <- fit[[paste0(key, ".ci")]]
  if (length(ci) != 2L) ci <- c(NA_real_, NA_real_)
  data.frame(
    estimate = as.numeric(fit[[estimate_key]]), low = as.numeric(ci[1]),
    high = as.numeric(ci[2]), p = as.numeric(fit[[paste0(key, ".p")]]),
    stringsAsFactors = FALSE
  )
}

m14_empty_row <- function(task, obj, cfg) {
  data.frame(
    outcome = task$outcome, period = task$period,
    analysis_type = task$analysis_type, module_ids = task$module_ids,
    task_id = task$task_id, pair_id = task$pair_id,
    direction = task$direction, path = task$path,
    exposure_type = task$exposure_type, exposure_id = task$exposure_id,
    exposure_name = task$exposure_name, mediator_type = task$mediator_type,
    mediator_id = task$mediator_id, mediator_name = task$mediator_name,
    bacteria_id = task$bacteria_id, bacteria_name = task$bacteria_name,
    metabolite_id = task$metabolite_id, metabolite_name = task$metabolite_name,
    analysis_cohort = obj$analysis_cohort,
    cohort_source = obj$cohort_source,
    n = obj$n, events = obj$events, n_OWO = obj$n_OWO,
    exposure_contrast = if (task$exposure == "A") "OWO: 0 vs 1" else
      "standardized feature: 0 vs +1 SD",
    nsims = cfg$nsims, nsims_success = NA_integer_,
    exposure_to_mediator_beta = NA_real_, exposure_to_mediator_low = NA_real_,
    exposure_to_mediator_high = NA_real_, exposure_to_mediator_p = NA_real_,
    mediator_to_outcome_logOR = NA_real_, mediator_to_outcome_OR = NA_real_,
    mediator_to_outcome_low = NA_real_, mediator_to_outcome_high = NA_real_,
    mediator_to_outcome_p = NA_real_,
    TE = NA_real_, TE_low = NA_real_, TE_high = NA_real_, TE_p = NA_real_,
    ADE = NA_real_, ADE_low = NA_real_, ADE_high = NA_real_, ADE_p = NA_real_,
    ACME = NA_real_, ACME_low = NA_real_, ACME_high = NA_real_, ACME_p = NA_real_,
    PM = NA_real_, PM_low = NA_real_, PM_high = NA_real_, PM_p = NA_real_,
    ACME_q = NA_real_, PM_q = NA_real_, ACME_q10 = FALSE,
    PM_positive_q10 = FALSE, PM_interpretable = FALSE,
    effect_scale = "marginal risk difference; PM=proportion",
    adjustment = "unweighted_extended",
    OWO_role = if (task$exposure == "A") "exposure" else "adjustment_covariate",
    OWO_adjusted = task$exposure != "A",
    FPC_used = FALSE, weighting_used = FALSE,
    covariates = paste(obj$actual_covariates, collapse = ";"),
    interpretation_flags = "", status = "not_run", warnings = "", error = "",
    seed = task$seed, stringsAsFactors = FALSE
  )
}

m14_fit_task <- function(task, obj, cfg) {
  analysis_obj <- obj$cohorts[[task$cohort_key]]
  if (is.null(analysis_obj)) stop("Unknown analysis cohort: ", task$cohort_key)
  row <- m14_empty_row(task, analysis_obj, cfg)
  warnings <- character()
  result <- NULL
  tryCatch(withCallingHandlers({
    d <- analysis_obj$data
    covars <- c(if (task$exposure != "A") "A", analysis_obj$covars)
    fm <- reformulate(c(task$exposure, covars), response = task$mediator)
    fy <- reformulate(c(task$exposure, task$mediator, covars), response = "Y")
    mediator_model <- lm(fm, data = d, na.action = na.fail)
    outcome_model <- glm(fy, data = d, family = binomial(), na.action = na.fail,
      control = glm.control(maxit = 100L))
    # mediate() updates stored calls; keep literal formulas in those calls.
    mediator_model$call$formula <- fm
    outcome_model$call$formula <- fy
    m14_validate_model(mediator_model)
    m14_validate_model(outcome_model, binary = TRUE)

    set.seed(task$seed)
    fit <- suppressMessages(mediation::mediate(
      mediator_model, outcome_model, treat = task$exposure,
      mediator = task$mediator, control.value = 0, treat.value = 1,
      boot = FALSE, sims = cfg$nsims, robustSE = FALSE, long = TRUE
    ))

    a <- coef(summary(mediator_model))[task$exposure, ]
    b <- coef(summary(outcome_model))[task$mediator, ]
    row$exposure_to_mediator_beta <- a[1]
    row$exposure_to_mediator_low <- a[1] - qnorm(.975) * a[2]
    row$exposure_to_mediator_high <- a[1] + qnorm(.975) * a[2]
    row$exposure_to_mediator_p <- a[4]
    row$mediator_to_outcome_logOR <- b[1]
    row$mediator_to_outcome_OR <- exp(b[1])
    row$mediator_to_outcome_low <- exp(b[1] - qnorm(.975) * b[2])
    row$mediator_to_outcome_high <- exp(b[1] + qnorm(.975) * b[2])
    row$mediator_to_outcome_p <- b[4]

    effect_specs <- list(
      TE = c("tau", "tau.coef"), ADE = c("z.avg", "z.avg"),
      ACME = c("d.avg", "d.avg"), PM = c("n.avg", "n.avg")
    )
    for (nm in names(effect_specs)) {
      e <- m14_effect(fit, effect_specs[[nm]][1], effect_specs[[nm]][2])
      row[[nm]] <- e$estimate
      row[[paste0(nm, "_low")]] <- e$low
      row[[paste0(nm, "_high")]] <- e$high
      row[[paste0(nm, "_p")]] <- e$p
    }
    valid <- is.finite(fit$tau.sims) & is.finite(fit$z.avg.sims) &
      is.finite(fit$d.avg.sims)
    row$nsims_success <- sum(valid)
    if (row$nsims_success != cfg$nsims) stop("Incomplete quasi-Bayesian simulations")

    flags <- c(
      "conditional_on_outcome_supervised_DIABLO_selection",
      if (task$cohort_key == "bacteria_full")
        "full_T1_16S_cohort" else
        "paired_16S_metabolomics_cohort",
      "pregnancy_covariates_may_be_post_exposure",
      "T1_temporality_requires_outcome_after_biospecimen"
    )
    if (task$analysis_type == "feature_bidirectional")
      flags <- c(flags, "same_timepoint_direction_not_identified")
    if (any(abs(coef(outcome_model)) > 20) ||
        any(fitted(outcome_model) < 1e-8 | fitted(outcome_model) > 1 - 1e-8))
      flags <- c(flags, "possible_separation_or_extreme_prediction")
    if (min(sum(d$Y == 1), sum(d$Y == 0)) / outcome_model$rank < 10)
      flags <- c(flags, "sparse_events_per_parameter")
    if (is.finite(row$TE_low) && is.finite(row$TE_high) &&
        row$TE_low <= 0 && row$TE_high >= 0)
      flags <- c(flags, "TE_CI_crosses_zero_PM_unstable")
    if (is.finite(row$PM) && (row$PM < 0 || row$PM > 1))
      flags <- c(flags, "PM_outside_0_1")
    if (is.finite(row$PM_low) && is.finite(row$PM_high) &&
        (row$PM_low < 0 || row$PM_high > 1))
      flags <- c(flags, "PM_CI_outside_0_1")
    if (is.finite(row$ACME) && is.finite(row$ADE) && row$ACME * row$ADE < 0)
      flags <- c(flags, "opposing_direct_indirect_effects")
    if (is.finite(row$ACME) && is.finite(row$TE) && row$ACME * row$TE < 0)
      flags <- c(flags, "ACME_TE_direction_discordant")
    row$PM_interpretable <- is.finite(row$PM) && row$PM >= 0 && row$PM <= 1 &&
      is.finite(row$TE_low) && is.finite(row$TE_high) &&
      !(row$TE_low <= 0 && row$TE_high >= 0) && row$ACME * row$TE > 0
    row$interpretation_flags <- paste(unique(flags), collapse = ";")
    row$status <- "ok"
    result <- list(
      row = row, fit = fit, mediator_model = mediator_model,
      outcome_model = outcome_model,
      formulas = list(mediator = fm, outcome = fy)
    )
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    row$status <<- "failed"
    row$error <<- conditionMessage(e)
    NULL
  })
  row$warnings <- paste(unique(warnings), collapse = " | ")
  if (is.null(result)) result <- list(row = row, fit = NULL,
    mediator_model = NULL, outcome_model = NULL, formulas = NULL)
  result$row <- row
  result
}

m14_run_task <- function(task, obj, cfg) {
  do.call(RNGkind, as.list(c("L'Ecuyer-CMRG", "Inversion", "Rejection")))
  m14_fit_task(task, obj, cfg)
}

m14_adjust_fdr <- function(rows, analysis, q_cutoff = 0.10) {
  if (!nrow(rows)) return(list(results = rows, families = data.frame()))
  keys <- if (analysis == "pair") c("outcome", "direction") else
    c("outcome", "mediator_type")
  rows$ACME_q <- rows$PM_q <- NA_real_
  families <- list()
  group <- interaction(rows[keys], drop = TRUE, lex.order = TRUE)
  for (idx in split(seq_len(nrow(rows)), group)) {
    for (effect in c("ACME", "PM")) {
      p <- rows[[paste0(effect, "_p")]][idx]
      good <- rows$status[idx] == "ok" & is.finite(p) & p >= 0 & p <= 1
      # p.adjust(n=planned) treats failed planned tests as p=1 while leaving
      # their displayed q values missing.
      if (any(good)) rows[[paste0(effect, "_q")]][idx[good]] <-
        p.adjust(p[good], method = "BH", n = length(idx))
      families[[length(families) + 1L]] <- cbind(
        rows[idx[1], keys, drop = FALSE],
        data.frame(effect = effect, planned = length(idx), evaluable = sum(good),
          fdr_n = length(idx), method = "BH", stringsAsFactors = FALSE)
      )
    }
  }
  rows$ACME_q10 <- rows$status == "ok" & is.finite(rows$ACME_q) &
    rows$ACME_q < q_cutoff
  # Main-figure selection uses ACME q, not PM q.
  rows$PM_positive_q10 <- rows$status == "ok" & is.finite(rows$PM) &
    rows$PM > 0 & is.finite(rows$ACME_q) & rows$ACME_q < q_cutoff
  list(results = rows, families = m14_bind(families))
}

m14_ci <- function(est, low, high, multiplier = 1) {
  ifelse(is.finite(est),
    sprintf("%.3f [%.3f, %.3f]", multiplier * est,
            multiplier * low, multiplier * high), NA_character_)
}

m14_compact <- function(rows) {
  if (!nrow(rows)) return(rows)
  data.frame(
    outcome = rows$outcome, period = rows$period, module = rows$module_ids,
    direction = rows$direction, path = rows$path,
    analysis_cohort = rows$analysis_cohort,
    bacteria = rows$bacteria_name, metabolite = rows$metabolite_name,
    mediator_omics = rows$mediator_type, feature = rows$mediator_name,
    N_events = paste(rows$n, rows$events, sep = "/"),
    contrast = rows$exposure_contrast,
    PM_percent_CI95 = m14_ci(rows$PM, rows$PM_low, rows$PM_high, 100),
    PM_P = rows$PM_p, PM_q = rows$PM_q,
    PM_positive_q10 = rows$PM_positive_q10,
    ACME_percentage_points_CI95 = m14_ci(rows$ACME, rows$ACME_low,
                                         rows$ACME_high, 100),
    ACME_P = rows$ACME_p, ACME_q = rows$ACME_q,
    TE_percentage_points_CI95 = m14_ci(rows$TE, rows$TE_low, rows$TE_high, 100),
    TE_P = rows$TE_p,
    ADE_percentage_points_CI95 = m14_ci(rows$ADE, rows$ADE_low, rows$ADE_high, 100),
    ADE_P = rows$ADE_p,
    PM_interpretable = rows$PM_interpretable,
    effect_scale = rows$effect_scale,
    interpretation_flags = rows$interpretation_flags,
    status = rows$status, error = rows$error,
    stringsAsFactors = FALSE
  )
}

m14_export <- function(rows, families, objects, tasks, cfg, results) {
  folder <- cfg$output
  compact <- m14_compact(rows)
  counts_keys <- if (cfg$analysis == "pair") c("outcome", "direction") else
    c("outcome", "mediator_type")
  groups <- unique(rows[counts_keys])
  counts <- m14_bind(lapply(seq_len(nrow(groups)), function(i) {
    keep <- rep(TRUE, nrow(rows))
    for (nm in counts_keys) keep <- keep & rows[[nm]] == groups[[nm]][i]
    z <- rows[keep, , drop = FALSE]
    cbind(groups[i, , drop = FALSE], data.frame(
      planned = nrow(z), successful = sum(z$status == "ok"),
      PM_positive_q10 = sum(z$PM_positive_q10, na.rm = TRUE),
      ACME_q10 = sum(z$ACME_q10, na.rm = TRUE), stringsAsFactors = FALSE
    ))
  }))
  by_object <- function(field) m14_bind(lapply(objects, function(x) {
    z <- x[[field]]
    if (is.null(z) || !nrow(z)) return(data.frame())
    if (!"outcome" %in% names(z)) z$outcome <- x$outcome
    z
  }))
  by_cohort <- function(field) m14_bind(lapply(objects, function(x)
    m14_bind(lapply(x$cohorts, function(z) {
      value <- z[[field]]
      if (is.null(value) || !is.data.frame(value) || !nrow(value)) return(data.frame())
      if (!"outcome" %in% names(value)) value$outcome <- x$outcome
      value
    }))))
  tables <- list(
    summary = compact,
    PM_candidates = compact[rows$PM_positive_q10, , drop = FALSE],
    all_results = rows, FDR_families = families, counts = counts,
    stable_features = by_object("members"), feature_scaling = by_object("scalers"),
    cohort_summary = by_object("cohort_summary"),
    covariate_map = by_cohort("covariate_map"),
    dropped_covariates = by_cohort("dropped_covariates"),
    excluded_participants = by_cohort("excluded"),
    model_failures = rows[rows$status != "ok", , drop = FALSE]
  )
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(tables)) utils::write.csv(
    tables[[nm]], file.path(folder, paste0(nm, ".csv")), row.names = FALSE,
    fileEncoding = "UTF-8", na = ""
  )
  writexl::write_xlsx(tables, file.path(folder, "mediation_summary.xlsx"))
  m14_save_rds(list(
    config = cfg, tables = tables,
    prepared_inputs = objects, tasks = tasks, results = results,
    input_md5 = tools::md5sum(cfg$input), code_md5 = tools::md5sum(
      unique(c(cfg$cohort_helper_file, cfg$helper_file, cfg$entry_file))),
    full_16S_input_md5 = tools::md5sum(unique(unlist(lapply(objects,
      function(x) x$full_16S_files)))),
    session_info = sessionInfo(), completed_at = Sys.time()
  ), file.path(folder, "mediation_results.rds"))
  writeLines(c(
    paste("Version:", cfg$version),
    paste("Analysis:", cfg$analysis),
    paste("Outcomes:", paste(cfg$outcomes, collapse = ", ")),
    paste("Planned / successful:", nrow(rows), "/", sum(rows$status == "ok")),
    paste("PM>0 and ACME q<0.10:", sum(rows$PM_positive_q10, na.rm = TRUE)),
    paste("Direct mediation::mediate simulations per model:", cfg$nsims),
    paste("Output:", normalizePath(folder, winslash = "/", mustWork = FALSE))
  ), file.path(folder, "run_summary.txt"), useBytes = TRUE)
  invisible(tables)
}

m14_main <- function(analysis = c("pair", "single"), outcomes = c("GDM", "HDP"),
                     workers = 8L, nsims = 5000L) {
  analysis <- match.arg(analysis)
  cfg <- m14_config(analysis, outcomes, workers, nsims)
  if (!file.exists(cfg$input)) stop("Missing DIABLO result: ", cfg$input)
  for (file in unique(c(cfg$cohort_helper_file, cfg$helper_file, cfg$entry_file)))
    if (!file.exists(file)) stop("Missing analysis code: ", file)
  packages <- c("mediation", "writexl")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  if (!is.finite(cfg$workers) || cfg$workers < 1L || cfg$workers != floor(cfg$workers))
    stop("workers must be a positive integer")
  if (!is.finite(cfg$nsims) || cfg$nsims < 100L || cfg$nsims != floor(cfg$nsims))
    stop("nsims must be an integer >= 100")

  source(cfg$cohort_helper_file, local = environment(m14_main))
  discovery <- readRDS(cfg$input)
  full_16s <- d14c_load_full_t1_16s(cfg$outcomes, cfg$baseline, cfg$pregnancy)
  objects <- setNames(lapply(cfg$outcomes, function(y)
    m14_prepare_outcome(discovery, y, cfg, full_16s)), cfg$outcomes)
  for (outcome in names(objects)) objects[[outcome]]$full_16S_files <- full_16s$files
  tasks <- setNames(lapply(objects, m14_make_tasks, analysis = analysis, cfg = cfg),
                    cfg$outcomes)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  total_tasks <- sum(lengths(tasks))
  n_workers <- min(cfg$workers, max(1L, total_tasks))
  cluster <- NULL
  if (n_workers > 1L) {
    cluster <- parallel::makePSOCKcluster(n_workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(paths, wd) {
      .libPaths(paths)
      setwd(wd)
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
      NULL
    }, .libPaths(), getwd())
    parallel::clusterExport(cluster,
      ls(environment(m14_main), pattern = "^m14_"), envir = environment(m14_main))
  }

  results <- list()
  for (outcome in cfg$outcomes) {
    message(outcome, ": ", length(tasks[[outcome]]),
            " planned models; workers=", n_workers)
    if (!length(tasks[[outcome]])) next
    fitted <- if (is.null(cluster)) {
      lapply(tasks[[outcome]], m14_run_task, obj = objects[[outcome]], cfg = cfg)
    } else {
      parallel::parLapplyLB(cluster, tasks[[outcome]], m14_run_task,
        obj = objects[[outcome]], cfg = cfg)
    }
    results <- c(results, fitted)
  }
  rows <- m14_bind(lapply(results, `[[`, "row"))
  adjusted <- m14_adjust_fdr(rows, analysis, cfg$q_cutoff)
  rows <- adjusted$results
  stopifnot(nrow(rows) == total_tasks, !anyDuplicated(rows$task_id))
  m14_export(rows, adjusted$families, objects, tasks, cfg, results)
  message("Saved: ", normalizePath(cfg$output, winslash = "/"))
  invisible(rows)
}

