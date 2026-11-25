########################################################################
# Function for analysis 01_multiomics_profile
########################################################################


#' Integrated PERMANOVA (Adonis2) Engine for Multi-Omics Data
#'
#' @param data Input clinical and omics dataframe.
#' @param feature_cols Character vector of omics feature names.
#' @param target_cols Character vector of clinical phenotype names.
#' @param period_term Character string for longitudinal timepoint (e.g., "T1", "TA").
#' @param label_type Character string for omics modality label.
#' @param label_var Character string for phenotype label.
#' @param dist_method Distance metric ("euclidean" for metab; "bray" for micro).
#' @param apply_tss Logical; whether to apply Total Sum Scaling (relative abundance).
#' @param id_col Character string for subject ID column. Default is "id".
#' @param period_col Character string for timepoint column. Default is "period".
#' @param all_period_label Character string indicating unfiltered timepoint. Default is "TA".
#' @param by_mode adonis2 test mode: NULL (overall), "terms" (sequential), or "margin".
#' @param compute_adjusted_r2 Logical; if TRUE, fit the same distance-based model
#'   with vegan::dbrda() and report vegan::RsquareAdj() as R2_adj.
#' @param n_perms Number of permutations for significance testing. Default is 999.
#' @param n_cores Number of parallel processing threads. Default is 8.
#' @param seed Random seed for computational reproducibility. Default is 10.
func_permanova_integrated <- function(data, 
                                      feature_cols, 
                                      target_cols, 
                                      period_term, 
                                      label_type, 
                                      label_var,
                                      dist_method = "bray",
                                      apply_tss = FALSE,
                                      id_col = "id",
                                      period_col = "period",
                                      all_period_label = "TA",
                                      by_mode = NULL,
                                      compute_adjusted_r2 = FALSE,
                                      n_perms = 999, 
                                      n_cores = 8, 
                                      seed = 10) {
  
  # 1. Dynamic filtering of complete clinical cases
  temp_df <- data %>%
    filter(if_all(all_of(target_cols), ~ !is.na(.x)))
  
  # 2. Longitudinal timepoint subsetting
  if (period_term != all_period_label) {
    temp_df <- temp_df %>% 
      filter(.data[[period_col]] == period_term)
  }
  
  # 3. Statistical rigor check: Abort if sample size is insufficient
  if (nrow(temp_df) < 5) {
    warning(paste0("Insufficient sample size (N = ", nrow(temp_df), ") for: ", label_var, " [Period: ", period_term, "]"))
    return(NULL)
  }
  
  # 4. Extract feature matrix
  abund_mat <- temp_df %>%
    mutate(id_period = str_c(.data[[id_col]], "_", .data[[period_col]])) %>%
    column_to_rownames("id_period") %>%
    select(all_of(feature_cols))
  
  # 5. Optional: Total Sum Scaling (TSS) for compositional microbiome data
  if (apply_tss) {
    abund_mat <- vegan::decostand(abund_mat, method = "total", MARGIN = 1)
  }
  
  # 6. Defensive programming: Remove zero-sum rows to prevent NaN in Bray-Curtis
  if (dist_method == "bray") {
    non_zero_indices <- rowSums(abund_mat) > 0
    if (sum(!non_zero_indices) > 0) {
      abund_mat <- abund_mat[non_zero_indices, , drop = FALSE]
      temp_df <- temp_df[non_zero_indices, , drop = FALSE]
    }
  }
  
  # 7. Compute dissimilarity matrix
  dist_obj <- vegan::vegdist(abund_mat, method = dist_method)
  dist_mat <- as.matrix(dist_obj)
  
  # 8. Construct dynamic statistical formula
  formula_str <- paste0("dist_mat ~ ", str_flatten(target_cols, collapse = " + "))
  
  # 9. Execute permutation multivariate analysis of variance
  set.seed(seed)
  test_res <- vegan::adonis2(
    as.formula(formula_str),
    data = temp_df,
    by = by_mode, 
    permutations = n_perms,
    parallel = n_cores
  )
  
  # 10. Compile standardized output summary
  res <- data.frame(
    type = label_type,
    variable = label_var,
    period = period_term,
    sample_size = nrow(temp_df),
    F_value = test_res$F[1],
    R2 = test_res$R2[1],
    P_value = test_res$`Pr(>F)`[1],
    stringsAsFactors = FALSE
  )

  # 11. Adjust the group-level explained variance for model complexity.
  # The permutation test above remains based on adonis2; dbrda is used only
  # to obtain the adjusted R-squared for the same distance and predictors.
  if (compute_adjusted_r2) {
    dbrda_formula_str <- paste0(
      "dist_obj ~ ",
      str_flatten(target_cols, collapse = " + ")
    )
    dbrda_fit <- vegan::dbrda(
      as.formula(dbrda_formula_str),
      data = temp_df
    )
    res$R2_adj <- unname(vegan::RsquareAdj(dbrda_fit)$adj.r.squared)
  }
  
  return(res)
}

########################################################################
# Functions for differential multi-omics analysis
########################################################################

#' Fit one inverse-probability-weighted pBMI-omics-feature model
#'
#' @param data Input data frame.
#' @param feature Standardized omics-feature column.
#' @param exposure pBMI exposure column.
#' @param covars Adjustment variables.
#' @param model Either "cross_sectional" or "mixed".
#' @param period Timepoint for a cross-sectional model.
#' @param exposure_label Human-readable exposure label.
#' @param target_population_size Size of the original target population
#'   represented by the analysis exposure.
#'
#' @return A tibble containing the exposure coefficient(s).
func_omics_regression <- function(data,
                                  feature,
                                  exposure,
                                  covars,
                                  model = c("cross_sectional", "mixed"),
                                  period = NULL,
                                  exposure_label = exposure,
                                  id_col = "id",
                                  period_col = "period",
                                  probability_col = "prob",
                                  target_population_size = 5212) {
  model <- match.arg(model)
  required_cols <- unique(c(
    feature, exposure, covars, id_col, period_col, probability_col
  ))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Variables absent from input data: ", paste(missing_cols, collapse = ", "))
  }
  if (length(target_population_size) != 1 ||
      is.na(target_population_size) ||
      target_population_size <= 0) {
    stop("`target_population_size` must be one positive number.")
  }

  analysis_data <- data
  if (model == "cross_sectional") {
    if (is.null(period)) stop("`period` is required for a cross-sectional model.")
    analysis_data <- analysis_data %>%
      dplyr::filter(.data[[period_col]] == .env$period)
  }
  analysis_data <- analysis_data %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(required_cols), ~ !is.na(.x)))

  empty_result <- function(status_text) {
    tibble::tibble(
      feature = feature,
      exposure = exposure_label,
      term = NA_character_,
      model = model,
      period = ifelse(model == "mixed", "Pooled", period),
      n_obs = nrow(analysis_data),
      n_id = dplyr::n_distinct(analysis_data[[id_col]]),
      target_population_size = target_population_size,
      beta = NA_real_,
      se = NA_real_,
      lower_ci = NA_real_,
      upper_ci = NA_real_,
      p_value = NA_real_,
      status = status_text
    )
  }
  if (nrow(analysis_data) < 20 ||
      dplyr::n_distinct(analysis_data[[exposure]]) < 2) {
    return(empty_result("insufficient_data"))
  }

  formula_text <- paste(
    feature, "~", paste(c(exposure, covars), collapse = " + ")
  )
  tryCatch({
    if (model == "cross_sectional") {
      analysis_data$.sampling_probability <- analysis_data[[probability_col]]
      analysis_data$.target_population_fpc <- target_population_size
      survey_design <- survey::svydesign(
        id = ~1,
        prob = ~.sampling_probability,
        fpc = ~.target_population_fpc,
        data = analysis_data
      )
      fit <- survey::svyglm(
        stats::as.formula(formula_text),
        design = survey_design,
        family = stats::gaussian()
      )
    } else {
      analysis_data <- analysis_data %>%
        dplyr::group_by(.data[[period_col]]) %>%
        dplyr::mutate(
          .analysis_weight = (1 / .data[[probability_col]]) *
            target_population_size /
            sum(1 / .data[[probability_col]])
        ) %>%
        dplyr::ungroup()
      fit <- lmerTest::lmer(
        stats::as.formula(paste0(formula_text, " + (1|", id_col, ")")),
        data = analysis_data,
        weights = .analysis_weight,
        REML = FALSE
      )
    }

    coefficient_table <- summary(fit)$coefficients %>%
      as.data.frame() %>%
      tibble::rownames_to_column("term")
    p_column <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), names(coefficient_table))[1]

    coefficient_table %>%
      dplyr::filter(stringr::str_starts(.data$term, exposure)) %>%
      dplyr::transmute(
        feature = feature,
        exposure = exposure_label,
        term = .data$term,
        model = model,
        period = ifelse(model == "mixed", "Pooled", period),
        n_obs = stats::nobs(fit),
        n_id = dplyr::n_distinct(analysis_data[[id_col]]),
        target_population_size = target_population_size,
        beta = .data$Estimate,
        se = .data$`Std. Error`,
        lower_ci = .data$Estimate - stats::qnorm(0.975) * .data$`Std. Error`,
        upper_ci = .data$Estimate + stats::qnorm(0.975) * .data$`Std. Error`,
        p_value = .data[[p_column]],
        status = "ok"
      )
  }, error = function(error_condition) {
    empty_result(paste0("error: ", conditionMessage(error_condition)))
  })
}




#' Compare the relative explanatory contribution of pBMI and gestational week
#' in weighted longitudinal mixed-effects models.
#'
#' The four primary models use the same complete-case data and calibrated
#' period-specific inverse-probability weights. Marginal Nakagawa R2 is used
#' because the requested comparison concerns fixed-effect explanatory power.
#' Weighted ICCs are reported from both M_wk and M_both.
func_omics_variance_contribution <- function(
    data, feature, pbmi_var = "pbmi_z", wk_var = "wk", covars,
    id_col = "id", period_col = "period", probability_col = "prob",
    target_population_size = 5212) {
  required <- unique(c(feature, pbmi_var, wk_var, covars, id_col,
                       period_col, probability_col))
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing variables: ", paste(missing, collapse = ", "))
  }
  if (wk_var %in% covars || pbmi_var %in% covars) {
    stop("`covars` must not contain pbmi_var or wk_var.")
  }

  d <- data %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(required), ~ !is.na(.x)))
  empty <- tibble::tibble(
    feature = feature,
    n_obs = nrow(d),
    n_id = dplyr::n_distinct(d[[id_col]]),
    target_population_size = target_population_size,
    r2_m0 = NA_real_, r2_m_pbmi = NA_real_, r2_m_wk = NA_real_,
    r2_m_both = NA_real_, delta_r2_pbmi = NA_real_,
    delta_r2_wk = NA_real_, shapley_pbmi = NA_real_,
    shapley_wk = NA_real_, r2_m_interaction = NA_real_,
    delta_r2_interaction = NA_real_, p_delta_r2_pbmi_lrt = NA_real_,
    p_delta_r2_wk_lrt = NA_real_, p_interaction_lrt = NA_real_,
    interaction_term = NA_character_, interaction_beta = NA_real_, 
    interaction_se = NA_real_, p_interaction_coef = NA_real_,
    icc_m_wk = NA_real_, icc_m_both = NA_real_,
    singular_m_wk = NA, singular_m_both = NA,
    status = "insufficient_data"
  )
  if (nrow(d) < 100 || dplyr::n_distinct(d[[id_col]]) < 30) return(empty)

  tryCatch({
    d <- d %>%
      dplyr::group_by(.data[[period_col]]) %>%
      dplyr::mutate(
        .analysis_weight = (1 / .data[[probability_col]]) *
          target_population_size / sum(1 / .data[[probability_col]])
      ) %>%
      dplyr::ungroup()
    # ICC is estimated with the normalized design weight used for variance
    # decomposition: wt_B = wt * sum(wt) / sum(wt^2), where wt = 1 / prob.
    raw_weight <- 1 / d[[probability_col]]
    d$.icc_weight <- raw_weight * sum(raw_weight) / sum(raw_weight^2)
    covar_rhs <- if (length(covars) == 0) "1" else paste(covars, collapse = " + ")
    random_rhs <- paste0(" + (1|", id_col, ")")
    fit <- function(rhs) lmerTest::lmer(
      stats::as.formula(paste0(feature, " ~ ", rhs, random_rhs)),
      data = d, weights = .analysis_weight, REML = FALSE
    )
    fit_m0 <- fit(covar_rhs)
    fit_pbmi <- fit(paste(pbmi_var, covar_rhs, sep = " + "))
    fit_wk <- fit(paste(wk_var, covar_rhs, sep = " + "))
    fit_both <- fit(paste(pbmi_var, wk_var, covar_rhs, sep = " + "))
    fit_interaction <- fit(paste0(pbmi_var, " * ", wk_var, " + ", covar_rhs))
    fit_wk_icc <- lmerTest::lmer(
      stats::as.formula(
        paste0(feature, " ~ ", wk_var, " + ", covar_rhs, random_rhs)
      ),
      data = d, weights = .icc_weight, REML = FALSE
    )
    fit_both_icc <- lmerTest::lmer(
      stats::as.formula(paste0(feature, " ~ ", pbmi_var, " + ", wk_var, " + ", covar_rhs, random_rhs)),
      data = d, weights = .icc_weight, REML = FALSE
    )

    marginal_r2 <- function(model) {
      tryCatch(
        as.numeric(performance::r2_nakagawa(model)[["R2_marginal"]]),
        error = function(e) NA_real_
      )
    }
    r2_m0 <- marginal_r2(fit_m0)
    r2_pbmi <- marginal_r2(fit_pbmi)
    r2_wk <- marginal_r2(fit_wk)
    r2_both <- marginal_r2(fit_both)
    r2_interaction <- marginal_r2(fit_interaction)
    # All nested comparisons use maximum-likelihood fits. The comparisons
    # test the unique contribution of pBMI, wk, and their interaction.
    pbmi_test <- stats::anova(fit_wk, fit_both)
    wk_test <- stats::anova(fit_pbmi, fit_both)
    interaction_test <- stats::anova(fit_both, fit_interaction)
    interaction_coef <- summary(fit_interaction)$coefficients %>%
      as.data.frame() %>%
      tibble::rownames_to_column("term") %>%
      dplyr::filter(stringr::str_detect(.data$term, ":"))
    interaction_term <- if (nrow(interaction_coef) == 0) NA_character_ else
      interaction_coef$term[[1]]
    interaction_beta <- if (nrow(interaction_coef) == 0) NA_real_ else
      interaction_coef$Estimate[[1]]
    interaction_se <- if (nrow(interaction_coef) == 0) NA_real_ else
      interaction_coef$`Std. Error`[[1]]
    interaction_df <- if (nrow(interaction_coef) == 0) NA_real_ else
      interaction_coef$df[[1]]
    interaction_p <- if (nrow(interaction_coef) == 0) NA_real_ else
      interaction_coef$`Pr(>|t|)`[[1]]
    p_delta_r2_pbmi_lrt <- pbmi_test$`Pr(>Chisq)`[2]
    p_delta_r2_wk_lrt <- wk_test$`Pr(>Chisq)`[2]
    p_interaction_lrt <- interaction_test$`Pr(>Chisq)`[2]

    calculate_icc <- function(model) {
      vc <- as.data.frame(lme4::VarCorr(model))
      subject_variance <- vc$vcov[vc$grp == id_col][1]
      subject_variance /
        (subject_variance + stats::sigma(model)^2)
    }
    icc_m_wk <- calculate_icc(fit_wk_icc)
    icc_m_both <- calculate_icc(fit_both_icc)
    tibble::tibble(
      feature = feature,
      n_obs = stats::nobs(fit_both),
      n_id = dplyr::n_distinct(d[[id_col]]),
      target_population_size = target_population_size,
      r2_m0 = r2_m0,
      r2_m_pbmi = r2_pbmi,
      r2_m_wk = r2_wk,
      r2_m_both = r2_both,
      delta_r2_pbmi = r2_both - r2_wk,
      delta_r2_wk = r2_both - r2_pbmi,
      shapley_pbmi = 0.5 * ((r2_pbmi - r2_m0) + (r2_both - r2_wk)),
      shapley_wk = 0.5 * ((r2_wk - r2_m0) + (r2_both - r2_pbmi)),
      r2_m_interaction = r2_interaction,
      delta_r2_interaction = r2_interaction - r2_both,
      p_delta_r2_pbmi_lrt = p_delta_r2_pbmi_lrt,
      p_delta_r2_wk_lrt = p_delta_r2_wk_lrt,
      p_interaction_lrt = p_interaction_lrt,
      interaction_term = interaction_term,
      interaction_beta = interaction_beta,
      interaction_se = interaction_se,
      p_interaction_coef = interaction_p,
      icc_m_wk = icc_m_wk,
      icc_m_both = icc_m_both,
      singular_m_wk = lme4::isSingular(fit_wk_icc, tol = 1e-4),
      singular_m_both = lme4::isSingular(fit_both_icc, tol = 1e-4),
      status = "ok"
    )
  }, error = function(e) {
    empty %>% dplyr::mutate(status = paste0("error: ", conditionMessage(e)))
  })
}


#' Run pBMI-feature models across omics features and timepoints
func_run_omics_models <- function(data,
                                  feature_cols,
                                  exposure,
                                  covars,
                                  periods = c("T1", "T2", "T3"),
                                  exposure_label = exposure,
                                  target_population_size = 5212,
                                  period_covars = NULL) {
  if (is.null(period_covars)) {
    period_covars <- stats::setNames(
      rep(list(covars), length(periods)),
      periods
    )
  }
  missing_period_covars <- setdiff(periods, names(period_covars))
  if (length(missing_period_covars) > 0) {
    stop(
      "`period_covars` is missing: ",
      paste(missing_period_covars, collapse = ", ")
    )
  }

  cross_results <- tidyr::expand_grid(
    feature = feature_cols,
    period = periods
  ) %>%
    purrr::pmap_dfr(function(feature, period) {
      func_omics_regression(
        data = data,
        feature = feature,
        exposure = exposure,
        covars = period_covars[[period]],
        model = "cross_sectional",
        period = period,
        exposure_label = exposure_label,
        target_population_size = target_population_size
      )
    })

  mixed_results <- purrr::map_dfr(
    feature_cols,
    ~ func_omics_regression(
      data = data,
      feature = .x,
      exposure = exposure,
      covars = covars,
      model = "mixed",
      exposure_label = exposure_label,
      target_population_size = target_population_size
    )
  )

  dplyr::bind_rows(cross_results, mixed_results) %>%
    dplyr::group_by(.data$exposure, .data$model, .data$period, .data$term) %>%
    dplyr::mutate(
      p_fdr = stats::p.adjust(.data$p_value, method = "BH"),
      p_bonferroni = stats::p.adjust(.data$p_value, method = "bonferroni")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$exposure, .data$period, .data$p_value)
}




#' Append feature annotations to model results
#'
#' @param results Model-result data frame containing `feature`.
#' @param info_data Feature annotation data frame.
#' @param info_key Feature identifier column in `info_data`.
#' @param omics_type Label identifying the omics layer.
#'
#' @return Results left-joined to a deduplicated annotation table.
func_append_feature_info <- function(results,
                                     info_data,
                                     info_key,
                                     omics_type) {
  if (!"feature" %in% names(results)) {
    stop("`results` must contain a `feature` column.")
  }
  if (!info_key %in% names(info_data)) {
    stop("Annotation key absent from `info_data`: ", info_key)
  }

  info_unique <- info_data %>%
    dplyr::distinct(.data[[info_key]], .keep_all = TRUE)

  results %>%
    dplyr::mutate(omics_type = omics_type, .after = "feature") %>%
    dplyr::left_join(
      info_unique,
      by = stats::setNames(info_key, "feature")
    )
}


########################################################################
# Functions for unweighted external-cohort multi-omics analysis
########################################################################

#' Fit one unweighted external-cohort omics association model
func_omics_regression_unweighted <- function(data,
                                             feature,
                                             exposure,
                                             covars,
                                             model = c("cross_sectional", "mixed"),
                                             period = NULL,
                                             exposure_label = exposure,
                                             id_col = "id",
                                             period_col = "period") {
  model <- match.arg(model)
  required_cols <- unique(c(feature, exposure, covars, id_col, period_col))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Variables absent from input data: ", paste(missing_cols, collapse = ", "))
  }

  analysis_data <- data
  if (model == "cross_sectional") {
    if (is.null(period)) stop("`period` is required for a cross-sectional model.")
    analysis_data <- analysis_data %>%
      dplyr::filter(.data[[period_col]] == .env$period)
  }
  analysis_data <- analysis_data %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(required_cols), ~ !is.na(.x)))

  empty_result <- function(status_text) {
    tibble::tibble(
      feature = feature,
      exposure = exposure_label,
      term = NA_character_,
      model = model,
      period = ifelse(model == "mixed", "Pooled", period),
      n_obs = nrow(analysis_data),
      n_id = dplyr::n_distinct(analysis_data[[id_col]]),
      beta = NA_real_,
      se = NA_real_,
      lower_ci = NA_real_,
      upper_ci = NA_real_,
      p_value = NA_real_,
      status = status_text
    )
  }
  if (nrow(analysis_data) < 20 ||
      dplyr::n_distinct(analysis_data[[exposure]]) < 2) {
    return(empty_result("insufficient_data"))
  }

  formula_text <- paste(
    feature, "~", paste(c(exposure, covars), collapse = " + ")
  )
  tryCatch({
    if (model == "cross_sectional") {
      fit <- stats::lm(stats::as.formula(formula_text), data = analysis_data)
    } else {
      fit <- lmerTest::lmer(
        stats::as.formula(paste0(formula_text, " + (1|", id_col, ")")),
        data = analysis_data,
        REML = FALSE
      )
    }

    coefficient_table <- summary(fit)$coefficients %>%
      as.data.frame() %>%
      tibble::rownames_to_column("term")
    p_column <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), names(coefficient_table))[1]

    coefficient_table %>%
      dplyr::filter(stringr::str_starts(.data$term, exposure)) %>%
      dplyr::transmute(
        feature = feature,
        exposure = exposure_label,
        term = .data$term,
        model = model,
        period = ifelse(model == "mixed", "Pooled", period),
        n_obs = stats::nobs(fit),
        n_id = dplyr::n_distinct(analysis_data[[id_col]]),
        beta = .data$Estimate,
        se = .data$`Std. Error`,
        lower_ci = .data$Estimate - stats::qnorm(0.975) * .data$`Std. Error`,
        upper_ci = .data$Estimate + stats::qnorm(0.975) * .data$`Std. Error`,
        p_value = .data[[p_column]],
        status = "ok"
      )
  }, error = function(error_condition) {
    empty_result(paste0("error: ", conditionMessage(error_condition)))
  })
}


#' Run unweighted external-cohort models across features and timepoints
func_run_omics_models_unweighted <- function(data,
                                             feature_cols,
                                             exposure,
                                             covars,
                                             periods = c("T1", "T2", "T3"),
                                             exposure_label = exposure) {
  cross_results <- tidyr::expand_grid(
    feature = feature_cols,
    period = periods
  ) %>%
    purrr::pmap_dfr(function(feature, period) {
      func_omics_regression_unweighted(
        data = data,
        feature = feature,
        exposure = exposure,
        covars = covars,
        model = "cross_sectional",
        period = period,
        exposure_label = exposure_label
      )
    })

  mixed_results <- purrr::map_dfr(
    feature_cols,
    ~ func_omics_regression_unweighted(
      data = data,
      feature = .x,
      exposure = exposure,
      covars = covars,
      model = "mixed",
      exposure_label = exposure_label
    )
  )

  dplyr::bind_rows(cross_results, mixed_results) %>%
    dplyr::group_by(.data$exposure, .data$model, .data$period, .data$term) %>%
    dplyr::mutate(
      p_fdr = stats::p.adjust(.data$p_value, method = "BH"),
      p_bonferroni = stats::p.adjust(.data$p_value, method = "bonferroni")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$exposure, .data$period, .data$p_value)
}


