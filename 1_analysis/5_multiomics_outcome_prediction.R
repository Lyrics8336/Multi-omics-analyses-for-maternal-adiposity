##############################################################################################
# Analyses for Figure 5
# Outcome prediction using OWO-associated multi-omics features
# See plot_fig5.R for the corresponding figure-generation code.
###############################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(writexl)
  library(caret)
  library(glmnet)
  library(pROC)
  library(PRROC)
  library(doParallel)
})

## 1. Load and assemble THSBC analysis data ----------------------------------
load("data/01_thsbc_cleaned_data/metabolites_1292_991.RData")
load("data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData")
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")
load("data/01_thsbc_cleaned_data/outcomes_5212.RData")
load("data/01_thsbc_cleaned_data/wt_omics.rdata")

DF_metab_micro <- total_pop %>%
  inner_join(otc, by = "id") %>%
  inner_join(metab_tab_st, by = "id") %>%
  inner_join(bact_fung_clrz, by = c("id", "period")) %>%
  inner_join(cov_tab, by = c("id", "period")) %>%
  inner_join(wt_metab, by = c("id", "period"))

DF_metab_micro_complete <- total_pop %>%
  inner_join(otc, by = "id") %>%
  inner_join(metab_tab_st, by = "id") %>%
  inner_join(bact_fung_clrz, by = c("id", "period")) %>%
  inner_join(cov_tab, by = c("id", "period")) %>%
  inner_join(wt_inter747, by = c("id", "period"))

## 2. Prespecified outcomes, features and incremental models -----------------
binary_outcomes <- c("GDM", "HDP", "PTB", "MAC", "LGA", "owo_m36")
offspring_outcomes <- c("PTB", "MAC", "LGA", "owo_m36")

load("result/PART2/1_diff_metabolites_sig_sets_thsbc.rdata")
full_metab <- unique(sig_metab_sets$sigid_metab_owo)

load("result/PART3/1_diff_16s_its_sig_sets_thsbc.rdata")
full_bact <- unique(sig_16s_sets$sigid_16s_owo)
full_fung <- unique(sig_its_sets$sigid_its_owo)
full_micro <- unique(c(full_bact, full_fung))
full_multiomics <- unique(c(full_metab, full_micro))

basic_common <- c(
  "age", "parity", "smk", "drk", "physic_activity",
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
  "edu", "medical_insurance", "income", "wk", "antibiotic"
)

# B0 is the clinically screened baseline. Every model for the same
# outcome/window uses this identical selected clinical set. pBMI is fixed
# (not screened) in B1 and all omics-augmented models, so the incremental
# predictive value of each omics layer is assessed beyond B1. In elastic-net
# models clinical variables have penalty.factor = 0; only omics are penalized.
model_definitions <- list(
  B0 = list(include_pbmi = FALSE, omics = character(0)),
  B1 = list(include_pbmi = TRUE, omics = character(0)),
  B1_Microbiome = list(include_pbmi = TRUE, omics = full_micro),
  B1_Metabolome = list(include_pbmi = TRUE, omics = full_metab),
  B1_Multiomics = list(include_pbmi = TRUE, omics = full_multiomics)
)

prediction_dir <- "result/PART4/prediction/3.6.1elastic_net/"
dir.create(prediction_dir, recursive = TRUE, showWarnings = FALSE)

## 3. Utilities ---------------------------------------------------------------

as_binary_outcome <- function(x) {
  x_chr <- as.character(x)
  ifelse(x_chr %in% c("1", "TRUE", "yes", "Yes"), 1L,
         ifelse(x_chr %in% c("0", "FALSE", "no", "No"), 0L, NA_integer_))
}

safe_first <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) x[[1]] else NA
}

safe_auc <- function(obs, prob) {
  if (length(unique(obs)) < 2L) return(NA_real_)
  tryCatch(as.numeric(pROC::auc(pROC::roc(
    obs, prob, levels = c("x0", "x1"), direction = "<", quiet = TRUE
  ))), error = function(e) NA_real_)
}

safe_pr_auc <- function(obs, prob) {
  if (length(unique(obs)) < 2L) return(NA_real_)
  tryCatch(as.numeric(PRROC::pr.curve(
    scores.class0 = prob[obs == "x1"], scores.class1 = prob[obs == "x0"],
    curve = FALSE
  )$auc.integral), error = function(e) NA_real_)
}

prediction_metrics <- function(obs, prob, threshold = 0.5) {
  pred <- factor(ifelse(prob >= threshold, "x1", "x0"), levels = c("x0", "x1"))
  cm <- caret::confusionMatrix(pred, obs, positive = "x1")
  precision <- unname(cm$byClass["Pos Pred Value"])
  recall <- unname(cm$byClass["Sensitivity"])
  f1 <- if (is.na(precision) || is.na(recall) || precision + recall == 0) NA_real_ else {
    2 * precision * recall / (precision + recall)
  }
  c(
    ROC = safe_auc(obs, prob), PR_AUC = safe_pr_auc(obs, prob),
    Accuracy = unname(cm$overall["Accuracy"]),
    Balanced_Accuracy = mean(unname(cm$byClass[c("Sensitivity", "Specificity")]), na.rm = TRUE),
    Sensitivity = unname(cm$byClass["Sensitivity"]),
    Specificity = unname(cm$byClass["Specificity"]), Precision = precision, F1 = f1
  )
}

make_stratified_folds <- function(y, k = 5L, repeats = 10L, seed = 111L) {
  if (length(levels(factor(y))) != 2L || min(table(y)) < k) {
    stop("Each outcome class must contain at least k observations.")
  }
  set.seed(seed)
  folds <- vector("list", repeats)
  for (repeat_id in seq_len(repeats)) {
    folds[[repeat_id]] <- caret::createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  }
  names(folds) <- paste0("Repeat", seq_len(repeats))
  folds
}

screen_basic <- function(data, outcome = "group", basic_candidates, p_cut = 0.25) {
  map_dfr(intersect(unique(basic_candidates), names(data)), function(covariate) {
    fit <- tryCatch(glm(stats::as.formula(paste0(outcome, " ~ `", covariate, "`")),
                        data = data, family = binomial()), error = function(e) NULL)
    p_value <- if (is.null(fit)) NA_real_ else {
      tryCatch(drop1(fit, test = "Chisq")[2, "Pr(>Chi)"], error = function(e) NA_real_)
    }
    tibble(clinical_feature = covariate, p_value = p_value,
           selected = !is.na(p_value) && p_value < p_cut)
  })
}

make_inner_foldid <- function(y, k = 5L, seed = 111L) {
  y <- factor(y, levels = c("x0", "x1"))
  if (length(unique(y)) != 2L || min(table(y)) < k) {
    stop("Each outer-training class must contain at least inner_k observations.")
  }
  set.seed(seed)
  folds <- caret::createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  foldid <- integer(length(y))
  for (i in seq_along(folds)) foldid[folds[[i]]] <- i
  foldid
}

screen_training_clinical <- function(
    data, train_idx, basic_candidates, forced_clinical = character(0),
    method = c("univariate_p", "lasso"), p_cut = 0.25,
    lasso_k = 5L, lasso_nlambda = 100L, seed = 111L) {
  method <- match.arg(method)
  train_data <- data[train_idx, , drop = FALSE]
  candidates <- intersect(unique(basic_candidates), names(train_data))

  if (method == "univariate_p") {
    screen <- screen_basic(train_data, "group", candidates, p_cut) %>%
      mutate(screen_method = method, screen_lambda = NA_real_,
             screen_coefficient = NA_real_)
  } else if (!length(candidates)) {
    screen <- tibble(
      clinical_feature = character(), p_value = double(), selected = logical(),
      screen_method = character(), screen_lambda = double(),
      screen_coefficient = double()
    )
  } else {
    clinical_frame <- train_data %>% select(all_of(candidates))
    clinical_matrix <- stats::model.matrix(~ ., data = clinical_frame)
    term_labels <- attr(stats::terms(~ ., data = clinical_frame), "term.labels")
    non_intercept <- colnames(clinical_matrix) != "(Intercept)"
    source_variable <- term_labels[attr(clinical_matrix, "assign")[non_intercept]]
    clinical_x <- clinical_matrix[, non_intercept, drop = FALSE]
    varying <- apply(clinical_x, 2, function(z) {
      z <- z[is.finite(z)]
      length(z) > 1L && is.finite(stats::sd(z)) && stats::sd(z) > 0
    })
    clinical_x <- clinical_x[, varying, drop = FALSE]
    source_variable <- source_variable[varying]
    y <- as.integer(train_data$group == "x1")
    if (ncol(clinical_x) < 2L) {
      # glmnet requires at least two columns. With zero/one varying clinical
      # column, keep the column only if a univariate likelihood-ratio test passes.
      fallback <- screen_basic(train_data, "group", candidates, p_cut)
      source_beta <- setNames(rep(NA_real_, length(candidates)), candidates)
      selected_sources <- fallback$clinical_feature[fallback$selected]
      chosen_lambda <- NA_real_
    } else {
      nfolds <- min(as.integer(lasso_k), min(table(y)))
      if (nfolds < 3L) stop("Too few events for clinical LASSO screening.")
      foldid <- make_inner_foldid(train_data$group, nfolds, seed)
      cv_fit <- glmnet::cv.glmnet(
        x = clinical_x, y = y, family = "binomial", alpha = 1,
        foldid = foldid, type.measure = "auc", standardize = TRUE,
        nlambda = lasso_nlambda, maxit = 1e6,
        type.logistic = "modified.Newton"
      )
      chosen_lambda <- cv_fit$lambda.1se
      beta <- as.matrix(coef(cv_fit, s = chosen_lambda))[-1, 1]
      beta_by_source <- split(unname(beta), source_variable)
      source_beta <- vapply(candidates, function(v) {
        b <- beta_by_source[[v]]
        if (is.null(b) || !length(b)) 0 else b[[which.max(abs(b))]]
      }, numeric(1))
      selected_sources <- names(source_beta)[source_beta != 0]
    }
    screen <- tibble(
      clinical_feature = candidates, p_value = NA_real_,
      selected = candidates %in% selected_sources,
      screen_method = method, screen_lambda = chosen_lambda,
      screen_coefficient = unname(source_beta[candidates])
    )
  }

  screen <- screen %>% mutate(forced = FALSE, used_in_all_models = .data$selected)
  if (length(forced_clinical)) {
    screen <- screen %>% filter(!.data$clinical_feature %in% forced_clinical)
    screen <- bind_rows(screen, tibble(
      clinical_feature = forced_clinical, p_value = NA_real_, selected = TRUE,
      screen_method = "forced", screen_lambda = NA_real_,
      screen_coefficient = NA_real_, forced = TRUE, used_in_all_models = TRUE
    ))
  }
  list(table = screen, keep = screen$clinical_feature[screen$used_in_all_models])
}

## 4. Period aggregation ------------------------------------------------------
aggregate_period_data <- function(data, period = "T1", outcome,
                                  feature_names = character(0), basic_names = character(0)) {
  feature_names <- intersect(unique(feature_names), names(data))
  basic_names <- intersect(unique(basic_names), names(data))
  required <- unique(c("id", "period", outcome, feature_names, basic_names))
  missing_variables <- setdiff(required, names(data))
  if (length(missing_variables)) stop("Missing variables: ", paste(missing_variables, collapse = ", "))
  d <- data %>% select(all_of(required)) %>% filter(.data$period %in% .env$period)
  if (length(period) == 1L) {
    return(d %>% filter(.data$period == .env$period[[1]]) %>%
      select(-any_of("period")) %>% distinct(.data$id, .keep_all = TRUE))
  }
  d %>%
    mutate(.period_order = match(.data$period, .env$period)) %>%
    arrange(.data$id, .data$.period_order) %>%
    group_by(.data$id) %>%
    summarise(across(all_of(feature_names), ~ mean(.x, na.rm = TRUE)),
              !!outcome := safe_first(.data[[outcome]]),
              across(all_of(basic_names), safe_first), .groups = "drop")
}

## 5. Nested elastic-net model ------------------------------------------------
# `model.matrix()` permits categorical clinical covariates. Its assign vector
# maps each resulting column back to the originating clinical/omics variable.
make_design_matrices <- function(train_data, test_data, predictors, clinical_features) {
  if (!length(predictors)) {
    return(list(
      x_train = matrix(numeric(0), nrow = nrow(train_data), ncol = 0L),
      x_test = matrix(numeric(0), nrow = nrow(test_data), ncol = 0L),
      penalty_factor = numeric(0), source_variable = character(0)
    ))
  }
  formula_x <- stats::reformulate(predictors, intercept = FALSE)
  x_train <- stats::model.matrix(formula_x, data = train_data)
  x_test <- stats::model.matrix(formula_x, data = test_data,
                                xlev = stats::.getXlevels(stats::terms(formula_x), train_data))
  if (!setequal(colnames(x_train), colnames(x_test))) {
    stop("Outer-training and test design matrices have incompatible columns.")
  }
  x_test <- x_test[, colnames(x_train), drop = FALSE]
  term_labels <- attr(stats::terms(formula_x), "term.labels")
  source_variable <- term_labels[attr(x_train, "assign")]
  penalty_factor <- ifelse(source_variable %in% clinical_features, 0, 1)
  list(x_train = x_train, x_test = x_test, penalty_factor = penalty_factor,
       source_variable = source_variable)
}

fit_unpenalized_logistic <- function(x_train, y_train, x_test) {
  train_with_intercept <- cbind(`(Intercept)` = 1, x_train)
  test_with_intercept <- cbind(`(Intercept)` = 1, x_test)
  fit <- suppressWarnings(glm.fit(
    x = train_with_intercept, y = y_train, family = binomial()
  ))
  beta <- fit$coefficients
  beta[!is.finite(beta)] <- 0
  names(beta) <- colnames(train_with_intercept)
  list(
    probability = as.numeric(plogis(test_with_intercept %*% beta)),
    coefficient = beta
  )
}

run_elastic_net_model <- function(
    analysis, outcome, window_label, model_name, configuration,
    outer_folds, outer_screens, population = "available_case", inner_k = 5L,
    alpha_grid = c(0.1, 0.25, 0.5, 0.75, 1),
    lambda_choice = c("lambda.1se", "lambda.min"), nlambda = 100L,
    seed = 111L, allow_parallel = FALSE, verbose = TRUE) {
  lambda_choice <- match.arg(lambda_choice)
  if (!length(alpha_grid) || any(!is.finite(alpha_grid)) ||
      any(alpha_grid <= 0 | alpha_grid > 1)) {
    stop("alpha_grid must contain values in (0, 1].")
  }
  alpha_grid <- sort(unique(alpha_grid))
  omics_features <- configuration$omics
  all_predictions <- list()
  all_tuning <- list()
  all_selection <- list()
  all_coefficients <- list()
  row_id <- 0L

  for (repeat_id in seq_along(outer_folds)) {
    for (outer_fold in seq_along(outer_folds[[repeat_id]])) {
      test_idx <- outer_folds[[repeat_id]][[outer_fold]]
      train_idx <- setdiff(seq_len(nrow(analysis)), test_idx)
      clinical_features <- outer_screens[[repeat_id]][[outer_fold]]$keep
      if (isTRUE(configuration$include_pbmi)) {
        clinical_features <- unique(c(clinical_features, "pbmi"))
      }
      predictors <- unique(c(clinical_features, omics_features))
      train_dat <- analysis[train_idx, c("id", "group", predictors), drop = FALSE]
      test_dat <- analysis[test_idx, c("id", "group", predictors), drop = FALSE]
      design <- make_design_matrices(
        train_dat %>% select(-id, -group), test_dat %>% select(-id, -group),
        predictors, clinical_features
      )
      keep_col <- if (ncol(design$x_train)) {
        apply(design$x_train, 2, function(z) {
          z <- z[is.finite(z)]
          length(z) > 1L && is.finite(stats::sd(z)) && stats::sd(z) > 0
        })
      } else logical(0)
      x_train <- design$x_train[, keep_col, drop = FALSE]
      x_test <- design$x_test[, keep_col, drop = FALSE]
      penalty_factor <- design$penalty_factor[keep_col]
      source_variable <- design$source_variable[keep_col]
      y_train <- as.integer(train_dat$group == "x1")

      if (verbose) message(sprintf(
        "%s | %s | %s | repeat %d/%d, fold %d/%d | train=%d, test=%d",
        outcome, window_label, model_name, repeat_id, length(outer_folds),
        outer_fold, length(outer_folds[[repeat_id]]), nrow(train_dat), nrow(test_dat)
      ))

      if (!any(penalty_factor > 0)) {
        final <- fit_unpenalized_logistic(x_train, y_train, x_test)
        probability <- final$probability
        fold_tuning <- tibble(
          alpha = NA_real_, lambda = NA_real_, inner_auc = NA_real_,
          lambda_rule = "unpenalized_clinical", n_selected_omics = 0L,
          selected_candidate = TRUE
        )
        omics_beta <- setNames(rep(0, length(omics_features)), omics_features)
      } else {
        foldid <- make_inner_foldid(
          train_dat$group, inner_k,
          seed + repeat_id * 10000L + outer_fold * 100L
        )
        cv_fits <- map(alpha_grid, function(alpha) {
          stabilized <- FALSE
          fit <- tryCatch(
            glmnet::cv.glmnet(
              x = x_train, y = y_train, family = "binomial", alpha = alpha,
              foldid = foldid, type.measure = "auc", standardize = TRUE,
              penalty.factor = penalty_factor, parallel = allow_parallel,
              nlambda = nlambda, lambda.min.ratio = 1e-3, maxit = 1e6,
              type.logistic = "modified.Newton", lower.limits = -10,
              upper.limits = 10
            ),
            error = function(e) {
              stabilized <<- TRUE
              if (verbose) message(
                "  alpha=", alpha,
                " required weak regularization of the clinical block: ",
                conditionMessage(e)
              )
              fallback_penalty <- penalty_factor
              fallback_penalty[fallback_penalty == 0] <- 1e-4
              glmnet::cv.glmnet(
                x = x_train, y = y_train, family = "binomial", alpha = alpha,
                foldid = foldid, type.measure = "auc", standardize = TRUE,
                penalty.factor = fallback_penalty, parallel = allow_parallel,
                nlambda = max(30L, as.integer(nlambda / 2L)),
                lambda.min.ratio = 1e-2, maxit = 2e6,
                type.logistic = "modified.Newton", lower.limits = -10,
                upper.limits = 10
              )
            }
          )
          list(fit = fit, stabilized = stabilized)
        })
        alpha_results <- map2_dfr(alpha_grid, cv_fits, function(alpha, cv_result) {
          cv_fit <- cv_result$fit
          lambda <- cv_fit[[lambda_choice]]
          lambda_index <- which.min(abs(cv_fit$lambda - lambda))
          beta <- as.matrix(coef(cv_fit, s = lambda))[-1, 1]
          selected_sources <- unique(source_variable[
            penalty_factor > 0 & is.finite(beta) & beta != 0
          ])
          tibble(
            alpha = alpha, lambda = lambda,
            inner_auc = cv_fit$cvm[[lambda_index]],
            lambda_rule = lambda_choice,
            n_selected_omics = length(selected_sources),
            cv_stabilized = cv_result$stabilized
          )
        })
        best_index <- alpha_results %>%
          mutate(.candidate_index = row_number()) %>%
          arrange(
            desc(.data$inner_auc), .data$n_selected_omics,
            desc(.data$alpha), desc(.data$lambda)
          ) %>%
          slice(1L) %>%
          pull(.data$.candidate_index)
        tuning_row <- alpha_results[best_index, , drop = FALSE]
        fold_tuning <- alpha_results %>%
          mutate(selected_candidate = row_number() == best_index)
        final_penalty <- penalty_factor
        if (isTRUE(tuning_row$cv_stabilized[[1]])) {
          final_penalty[final_penalty == 0] <- 1e-4
        }
        final_fit <- glmnet::glmnet(
          x = x_train, y = y_train, family = "binomial",
          alpha = tuning_row$alpha[[1]], lambda = tuning_row$lambda[[1]],
          standardize = TRUE, penalty.factor = final_penalty,
          maxit = 1e6, type.logistic = "modified.Newton",
          lower.limits = -10, upper.limits = 10
        )
        probability <- as.numeric(predict(
          final_fit, newx = x_test, s = tuning_row$lambda[[1]], type = "response"
        ))
        beta <- as.matrix(coef(final_fit, s = tuning_row$lambda[[1]]))[-1, 1]
        beta_by_source <- split(unname(beta[penalty_factor > 0]),
                                source_variable[penalty_factor > 0])
        omics_beta <- setNames(rep(0, length(omics_features)), omics_features)
        for (feature in intersect(names(beta_by_source), omics_features)) {
          omics_beta[[feature]] <- sum(beta_by_source[[feature]], na.rm = TRUE)
        }
        if (verbose) message(sprintf(
          "  selected alpha=%.2f, %s=%.6g, omics=%d, inner AUROC=%.4f",
          tuning_row$alpha[[1]], lambda_choice, tuning_row$lambda[[1]],
          sum(omics_beta != 0), tuning_row$inner_auc[[1]]
        ))
      }

      row_id <- row_id + 1L
      all_predictions[[row_id]] <- tibble(
        id = test_dat$id, outcome = outcome, window = window_label,
        model = model_name, population = population, repeat_id = repeat_id,
        outer_fold = outer_fold, observed = test_dat$group,
        probability = probability
      )
      all_tuning[[row_id]] <- fold_tuning %>% mutate(
        outcome = outcome, window = window_label, model = model_name,
        population = population, repeat_id = repeat_id,
        outer_fold = outer_fold, n_outer_train = nrow(train_dat),
        n_outer_test = nrow(test_dat), n_design_columns = ncol(x_train)
      )
      if (length(omics_features)) {
        selection_row <- tibble(
          feature = omics_features, coefficient = unname(omics_beta),
          selected = unname(omics_beta != 0), outcome = outcome,
          window = window_label, model = model_name, population = population,
          repeat_id = repeat_id, outer_fold = outer_fold
        )
        all_selection[[row_id]] <- selection_row
        all_coefficients[[row_id]] <- selection_row %>% filter(.data$selected)
      }
    }
  }
  predictions <- bind_rows(all_predictions)
  if (!nrow(predictions)) stop("No valid outer-fold predictions were generated.")
  prediction_counts <- predictions %>% count(.data$repeat_id, .data$id)
  if (any(prediction_counts$n != 1L)) {
    stop("Each participant must have exactly one out-of-fold prediction per repeat.")
  }
  metrics <- predictions %>%
    group_by(.data$repeat_id) %>%
    group_modify(~ as_tibble(as.list(prediction_metrics(
      .x$observed, .x$probability
    )))) %>%
    ungroup() %>%
    mutate(
      outcome = outcome, window = window_label, model = model_name,
      population = population, sample_size = nrow(analysis),
      events = sum(analysis$group == "x1")
    )
  list(
    predictions = predictions, metrics = metrics,
    tuning = bind_rows(all_tuning), selection = bind_rows(all_selection),
    coefficients = bind_rows(all_coefficients)
  )
}

auc_contrasts <- tribble(
  ~comparison, ~new_model, ~reference_model, ~interpretation,
  "M1-M0", "B1", "B0", "pBMI gain",
  "M2-M1", "B1_Microbiome", "B1", "microbiome gain",
  "M3-M1", "B1_Metabolome", "B1", "metabolome gain",
  "M4-M1", "B1_Multiomics", "B1", "complete multiomics gain",
  "M4-M3", "B1_Multiomics", "B1_Metabolome",
  "microbiome gain beyond metabolome"
)

make_paired_auc_differences <- function(metrics) {
  id_columns <- c(
    "population", "outcome", "window", "repeat_id", "sample_size", "events"
  )
  wide <- metrics %>%
    select(all_of(c(id_columns, "model", "ROC"))) %>%
    pivot_wider(names_from = "model", values_from = "ROC")
  missing_models <- setdiff(
    unique(c(auc_contrasts$new_model, auc_contrasts$reference_model)), names(wide)
  )
  if (length(missing_models)) {
    stop("Cannot calculate prespecified AUC contrasts; missing models: ",
         paste(missing_models, collapse = ", "))
  }
  map_dfr(seq_len(nrow(auc_contrasts)), function(i) {
    new_name <- auc_contrasts$new_model[[i]]
    reference_name <- auc_contrasts$reference_model[[i]]
    wide %>% transmute(
      across(all_of(id_columns)), comparison = auc_contrasts$comparison[[i]],
      interpretation = auc_contrasts$interpretation[[i]],
      new_model = new_name, reference_model = reference_name,
      AUC_new = .data[[new_name]], AUC_reference = .data[[reference_name]],
      delta_AUC = .data[[new_name]] - .data[[reference_name]]
    )
  })
}

## 6. Paired model set and feature-stability summaries -----------------------
run_model_set <- function(data, outcomes = binary_outcomes,
                          period_windows = list(T1 = "T1", T2 = "T2", T3 = "T3"),
                          population = "available_case", outer_k = 5L, outer_repeats = 10L,
                          inner_k = 5L, clinical_p_cut = 0.25,
                          clinical_screen_method = c("univariate_p", "lasso"),
                          clinical_lasso_k = 5L, clinical_lasso_nlambda = 100L,
                          alpha_grid = c(0.1, 0.25, 0.5, 0.75, 1),
                          lambda_choice = c("lambda.1se", "lambda.min"),
                          nlambda = 100L, seed = 111L,
                          n_cores = 8L, verbose = TRUE) {
  clinical_screen_method <- match.arg(clinical_screen_method)
  lambda_choice <- match.arg(lambda_choice)
  model_results <- list(); fold_registry <- list(); clinical_screening <- list()
  all_omics <- unique(unlist(map(model_definitions, "omics"), use.names = FALSE))
  n_available_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(n_available_cores)) n_available_cores <- 1L
  n_workers <- max(1L, min(as.integer(n_cores), as.integer(n_available_cores)))
  cluster <- NULL
  if (n_workers > 1L) {
    cluster <- parallel::makePSOCKcluster(n_workers)
    doParallel::registerDoParallel(cluster)
    on.exit({ parallel::stopCluster(cluster); foreach::registerDoSEQ() }, add = TRUE)
    if (verbose) message("Registered ", n_workers, " workers for glmnet inner CV.")
  } else {
    foreach::registerDoSEQ()
  }
  for (outcome_index in seq_along(outcomes)) {
    outcome <- outcomes[[outcome_index]]
    if (!outcome %in% names(data)) { warning("Skipping absent outcome: ", outcome); next }
    forced_clinical <- if (outcome %in% offspring_outcomes) "baby_sex" else character(0)
    for (window_index in seq_along(period_windows)) {
      window_label <- names(period_windows)[[window_index]]; period <- period_windows[[window_index]]
      if (is.null(window_label) || !nzchar(window_label)) window_label <- paste(period, collapse = "+")
      analysis <- aggregate_period_data(data, period, outcome, all_omics,
        unique(c(basic_common, forced_clinical, "pbmi"))) %>%
        mutate(group = factor(as_binary_outcome(.data[[outcome]]), levels = c(0, 1), labels = c("x0", "x1"))) %>%
        select(-all_of(outcome))
      required <- unique(c(basic_common, forced_clinical, "pbmi", all_omics))
      if (length(setdiff(required, names(analysis)))) stop("Missing predictors after aggregation: ",
        paste(setdiff(required, names(analysis)), collapse = ", "))
      analysis <- analysis %>% filter(!is.na(.data$group)) %>%
        filter(if_all(all_of(required), ~ !is.na(.x))) %>% arrange(.data$id)
      if (min(table(analysis$group)) < outer_k) { warning("Skipping ", outcome, " / ", window_label, ": too few events."); next }
      fold_seed <- seed + outcome_index * 1000L + window_index * 10L
      outer_folds <- make_stratified_folds(analysis$group, outer_k, outer_repeats, fold_seed)
      fold_registry[[paste(outcome, window_label, sep = "__")]] <- map2_dfr(
        outer_folds, seq_along(outer_folds), function(repeat_folds, repeat_id) {
          map2_dfr(repeat_folds, seq_along(repeat_folds), function(test_idx, outer_fold) {
            tibble(id = analysis$id[test_idx], repeat_id = repeat_id,
                   outer_fold = outer_fold)
          })
        }
      )
      outer_screens <- map(seq_along(outer_folds), function(repeat_id) {
        map(seq_along(outer_folds[[repeat_id]]), function(outer_fold) {
          test_idx <- outer_folds[[repeat_id]][[outer_fold]]
          train_idx <- setdiff(seq_len(nrow(analysis)), test_idx)
          screened <- screen_training_clinical(
            data = analysis, train_idx = train_idx,
            basic_candidates = basic_common, forced_clinical = forced_clinical,
            method = clinical_screen_method, p_cut = clinical_p_cut,
            lasso_k = clinical_lasso_k,
            lasso_nlambda = clinical_lasso_nlambda,
            seed = fold_seed + repeat_id * 1000L + outer_fold * 10L
          )
          clinical_screening[[paste(
            outcome, window_label, repeat_id, outer_fold, sep = "__"
          )]] <<- screened$table %>% mutate(
            outcome = outcome, window = window_label, population = population,
            repeat_id = repeat_id, outer_fold = outer_fold,
            n_outer_train = length(train_idx)
          )
          screened
        })
      })
      if (verbose) {
        selected_counts <- unlist(map(outer_screens, ~ map_int(.x, ~ length(.x$keep))))
        message(
          population, " | ", outcome, " | ", window_label,
          " | N=", nrow(analysis), " events=", sum(analysis$group == "x1"),
          " | clinical screen=", clinical_screen_method,
          " | selected per outer fold median=", median(selected_counts),
          " [", min(selected_counts), ", ", max(selected_counts), "]"
        )
      }
      for (model_name in names(model_definitions)) {
        configuration <- model_definitions[[model_name]]
        key <- paste(population, outcome, window_label, model_name, sep = "__")
        model_results[[key]] <- run_elastic_net_model(
          analysis = analysis, outcome = outcome, window_label = window_label,
          model_name = model_name, configuration = configuration,
          outer_folds = outer_folds, outer_screens = outer_screens,
          population = population, inner_k = inner_k, alpha_grid = alpha_grid,
          lambda_choice = lambda_choice, nlambda = nlambda, seed = fold_seed,
          allow_parallel = !is.null(cluster), verbose = verbose
        )
      }
    }
  }
  metrics <- map_dfr(model_results, "metrics"); predictions <- map_dfr(model_results, "predictions")
  tuning <- map_dfr(model_results, "tuning"); selection <- map_dfr(model_results, "selection")
  coefficients <- map_dfr(model_results, "coefficients"); clinical_screening <- bind_rows(clinical_screening)
  auc_differences <- make_paired_auc_differences(metrics)
  auc_difference_summary <- auc_differences %>%
    group_by(.data$population, .data$outcome, .data$window, .data$sample_size,
      .data$events, .data$comparison, .data$interpretation,
      .data$new_model, .data$reference_model) %>%
    summarise(
      mean_delta_AUC = mean(.data$delta_AUC, na.rm = TRUE),
      sd_delta_AUC = sd(.data$delta_AUC, na.rm = TRUE),
      repeats = sum(is.finite(.data$delta_AUC)),
      positive_repeats = sum(.data$delta_AUC > 0, na.rm = TRUE),
      .groups = "drop"
    )
  metric_names <- c("ROC", "PR_AUC", "Accuracy", "Balanced_Accuracy", "Sensitivity", "Specificity", "Precision", "F1")
  metric_summary <- metrics %>% pivot_longer(all_of(metric_names), names_to = "metric", values_to = "value") %>%
    group_by(.data$population, .data$outcome, .data$window, .data$model, .data$sample_size, .data$events, .data$metric) %>%
    summarise(estimate = mean(.data$value, na.rm = TRUE), sd_across_repeats = sd(.data$value, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(mean_sd = sprintf("%.3f (%.3f)", .data$estimate, .data$sd_across_repeats))
  coefficient_summary <- selection %>%
    filter(.data$selected) %>%
    group_by(.data$population, .data$outcome, .data$window, .data$model, .data$feature) %>%
    summarise(
      mean_coefficient_when_selected = mean(.data$coefficient, na.rm = TRUE),
      median_coefficient_when_selected = median(.data$coefficient, na.rm = TRUE),
      mean_abs_coefficient_when_selected = mean(abs(.data$coefficient), na.rm = TRUE),
      .groups = "drop"
    )
  feature_stability <- selection %>%
    group_by(.data$population, .data$outcome, .data$window, .data$model, .data$feature) %>%
    summarise(
      selection_frequency = mean(.data$selected),
      selected_outer_models = sum(.data$selected), total_outer_models = n(),
      positive_selection_frequency = mean(.data$coefficient > 0),
      negative_selection_frequency = mean(.data$coefficient < 0),
      .groups = "drop"
    ) %>%
    left_join(coefficient_summary,
      by = c("population", "outcome", "window", "model", "feature")) %>%
    arrange(.data$outcome, .data$window, .data$model, desc(.data$selection_frequency))
  list(models = model_results, metrics = metrics, metric_summary = metric_summary,
       auc_differences = auc_differences, auc_difference_summary = auc_difference_summary,
       predictions = predictions,
       tuning = tuning, selection = selection, feature_stability = feature_stability, coefficients = coefficients,
       clinical_screening = clinical_screening, outer_folds = fold_registry,
       settings = list(outer_k = outer_k, outer_repeats = outer_repeats, inner_k = inner_k,
         clinical_p_cut = clinical_p_cut,
         clinical_screen_method = clinical_screen_method,
         clinical_lasso_k = clinical_lasso_k,
         clinical_lasso_nlambda = clinical_lasso_nlambda,
         alpha_grid = alpha_grid, lambda_choice = lambda_choice,
         nlambda = nlambda, n_cores = n_workers, seed = seed))
}

save_model_set <- function(results, file_stem) {
  saveRDS(results, paste0(file_stem, ".rds"))
  writexl::write_xlsx(list(metric_summary = results$metric_summary, repeat_metrics = results$metrics,
    auc_difference_summary = results$auc_difference_summary,
    repeat_auc_differences = results$auc_differences,
    outer_predictions = results$predictions, tuning = results$tuning,
    clinical_screening = results$clinical_screening, feature_stability = results$feature_stability,
    elastic_net_coefficients = results$coefficients), paste0(file_stem, ".xlsx"))
}



## 7. Formal analysis commands -----------------------------------------------
thsbc_results_GDM_HDP <- run_model_set(
  data = DF_metab_micro, outcomes = c('GDM', 'HDP'), 
  period_windows = list(T1 = "T1", T2 = "T2", Mean_T1_T2 = c("T1", "T2")),
  population = "available_case", outer_k = 5, outer_repeats = 10, inner_k = 5,
  clinical_p_cut = 0.25, clinical_screen_method = "univariate_p",
  alpha_grid = c(0.1, 0.25, 0.5, 0.75, 1),
  lambda_choice = "lambda.1se", nlambda = 100,
  seed = 311, n_cores = 12, verbose = TRUE
)
save_model_set(thsbc_results_GDM_HDP, file.path(prediction_dir, "3.6.1_elastic_net_model_specific_lambda_GDM_HDP"))

thsbc_results_MAC_LGA_PTB_owo <- run_model_set(
  data = DF_metab_micro, outcomes = c('MAC', 'LGA', 'PTB', 'owo_m36'), 
  period_windows = list(T1 = "T1", T2 = "T2", T3 = "T3", Mean_T1_T2_T3 = c("T1", "T2", "T3")),
  population = "available_case", outer_k = 5, outer_repeats = 10, inner_k = 5,
  clinical_p_cut = 0.25, clinical_screen_method = "univariate_p",
  alpha_grid = c(0.1, 0.25, 0.5, 0.75, 1),
  lambda_choice = "lambda.1se", nlambda = 100,
  seed = 311, n_cores = 12, verbose = TRUE
)
save_model_set(thsbc_results_MAC_LGA_PTB_owo, file.path(prediction_dir, "3.6.1_elastic_net_model_specific_lambda_MAC_LGA_PTB_owo_m36"))

