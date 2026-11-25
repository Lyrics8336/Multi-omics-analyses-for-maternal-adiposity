##############################################################################################
# Analyses for Figure 6
# Development and internal evaluation of multi-omics-inferred pBMI
# Models:
#   Basic: clinical covariates.
#   GenoBMI: clinical covariates and BMI polygenic risk score (PRS).
#   MetBMI: clinical covariates and 123 OWO-associated metabolites.
#   MicroBMI: clinical covariates, 74 bacterial genera and 6 fungal genera.
#   Multi_omics_pBMI: clinical covariates, PRS and all selected omics features.
# See plot_fig6.R for the corresponding figure-generation code.
###############################################################################################

rm(list = ls())

library(tidyverse)
library(glmnet)
library(writexl)

## 1. Load and assemble analysis data ----------------------------------------

load("data/01_thsbc_cleaned_data/metabolites_1292_991.RData")
load("data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData")
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")
load("data/01_thsbc_cleaned_data/outcomes_5212.RData")

DF_total <- total_pop %>%
  inner_join(otc, by = "id") %>%
  inner_join(metab_tab, by = "id") %>%
  inner_join(bact_fung_tab, by = c("id", "period")) %>%
  inner_join(cov_tab, by = c("id", "period")) %>%
  filter(period == "T1", !is.na(PRS))

## pBMI-related omics features
load("result/PART2/1_diff_metabolites_sig_sets_thsbc.rdata")
rep_metab123 = sig_metab_sets$sigid_metab_owo
load("result/PART3/1_diff_16s_its_sig_sets_thsbc.rdata")
rep_bact74_fung6 = c(sig_16s_sets$sigid_16s_owo, sig_its_sets$sigid_its_owo)

## basic covariates
Covars_common <- c(
  "age", "parity", "smk", "drk", "physic_activity",
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
  "edu", "medical_insurance", "income"
)

## 2. Configuration and input validation ------------------------------------

outcome_var <- "pbmi"
id_var <- "id"
prs_var <- "PRS"
outer_v <- 10L
inner_v <- 10L
random_seed <- 10L
bootstrap_reps <- 2000L
micro_pseudocount <- 1
output_dir <- "result/PART5/pbmi_model"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rep_metab123 <- unique(rep_metab123)
rep_bact74_fung6 <- unique(rep_bact74_fung6)
rep_bact74 <- intersect(rep_bact74_fung6, bact_uni)
rep_fung6 <- intersect(rep_bact74_fung6, fung_uni)

required_vars <- unique(c(
  id_var, outcome_var, prs_var, Covars_common,
  metab_uni, bact_uni, fung_uni,
  rep_metab123, rep_bact74_fung6
))
missing_vars <- setdiff(required_vars, names(DF_total))

if (length(missing_vars) > 0L) {
  stop("Variables absent from DF_total: ", paste(missing_vars, collapse = ", "))
}
if (anyDuplicated(DF_total[[id_var]]) > 0L) {
  stop("DF_total contains duplicated participant IDs.")
}
if (length(rep_metab123) != 123L ||
    length(rep_bact74) != 74L ||
    length(rep_fung6) != 6L) {
  stop(
    "Unexpected feature counts: metabolites = ", length(rep_metab123),
    ", bacteria = ", length(rep_bact74),
    ", fungi = ", length(rep_fung6), "."
  )
}

model_complete_vars <- unique(c(
  id_var, outcome_var, prs_var, Covars_common,
  rep_metab123, rep_bact74_fung6
))
if (anyNA(DF_total[model_complete_vars])) {
  stop(
    "Missing values are present in variables required by the models. ",
    "Define an imputation strategy or a complete-case population before fitting."
  )
}
if (any(DF_total[[outcome_var]] <= 0)) {
  stop(outcome_var, " must be positive because it is modelled on the log scale.")
}

## 3. Transform the complete T1 multi-omics data ----------------------------

z_trans <- function(x) {
  (x - mean(x)) / sd(x)
}

get_z_parameters <- function(data, feature_type, transformation, pseudocount = NA_real_) {
  tibble(
    feature = names(data),
    feature_type = feature_type,
    transformation = transformation,
    center = map_dbl(data, mean),
    scale = map_dbl(data, sd),
    pseudocount = pseudocount
  )
}

# Metabolites: log transformation followed by Z-standardization.
if (any(as.matrix(DF_total[metab_uni]) <= 0)) {
  stop("All metabolite values must be positive for the requested log transformation.")
}
metab_log <- DF_total %>%
  select(all_of(metab_uni)) %>%
  mutate(across(everything(), log))

# Bacteria and fungi: CLR transformation performed separately, then Z-standardized.
bact_clr <- DF_total %>%
  select(all_of(bact_uni)) %>%
  vegan::decostand(method = "clr", MARGIN = 1, pseudocount = 1) %>%
  as_tibble()

fung_clr <- DF_total %>%
  select(all_of(fung_uni)) %>%
  vegan::decostand(method = "clr", MARGIN = 1, pseudocount = micro_pseudocount) %>%
  as_tibble()

# Z-standardize the transformed omics data.
metab_logz <- metab_log %>% mutate(across(everything(), z_trans))
bact_clrz <- bact_clr %>% mutate(across(everything(), z_trans))
fung_clrz <- fung_clr %>% mutate(across(everything(), z_trans))

# Remove the raw omics columns and bind the transformed columns explicitly.
DF_model <- DF_total %>%
  select(-all_of(metab_uni), -matches("bact_g"), -matches("fung_g")) %>%
  bind_cols(metab_logz, bact_clrz, fung_clrz)

preprocessing_parameters <- bind_rows(
  get_z_parameters(metab_log, "metabolite", "log_z"),
  get_z_parameters(bact_clr, "bacterium", "clr_z", micro_pseudocount),
  get_z_parameters(fung_clr, "fungus", "clr_z", micro_pseudocount)
)

if (any(!is.finite(preprocessing_parameters$scale) |
        preprocessing_parameters$scale == 0)) {
  stop("At least one omics feature has zero or undefined SD after transformation.")
}

## 4. Create one set of BMI-stratified outer folds --------------------------

make_stratified_folds <- function(y, v = 10L, seed = 10L) {
  if (length(y) < v) {
    stop("The number of observations must be at least the number of folds.")
  }

  breaks <- unique(quantile(
    y, probs = seq(0, 1, length.out = min(11L, length(y) + 1L)),
    na.rm = TRUE, type = 8
  ))
  strata <- cut(y, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  set.seed(seed)
  fold_id <- integer(length(y))
  for (stratum in sort(unique(strata))) {
    stratum_rows <- sample(which(strata == stratum))
    fold_labels <- rep(seq_len(v), length.out = length(stratum_rows))
    fold_id[stratum_rows] <- sample(fold_labels)
  }
  fold_id
}

outer_fold <- make_stratified_folds(
  DF_model[[outcome_var]], v = outer_v, seed = random_seed
)

folds <- tibble(
  id = DF_model[[id_var]],
  outer_fold = outer_fold,
  observed_pbmi = DF_model[[outcome_var]]
)

## 5. Construct model matrices ----------------------------------------------

make_design_matrix <- function(data, variables) {
  design_formula <- reformulate(variables, response = NULL)
  model.matrix(design_formula, data = data)[, -1, drop = FALSE]
}

x_basic <- make_design_matrix(DF_model, Covars_common)
x_geno <- cbind(x_basic, PRS = DF_model[[prs_var]])
x_metab <- cbind(x_basic, as.matrix(DF_model[rep_metab123]))
x_micro <- cbind(x_basic, as.matrix(DF_model[rep_bact74_fung6]))
x_multi <- cbind(
  x_basic,
  PRS = DF_model[[prs_var]],
  as.matrix(DF_model[rep_metab123]),
  as.matrix(DF_model[rep_bact74_fung6])
)

if (anyDuplicated(colnames(x_multi)) > 0L) {
  stop("Duplicated columns were generated in the multi-omics design matrix.")
}

## 6. Fit models and generate out-of-fold predictions -----------------------

fit_ols_oof <- function(x, y_log, fold_id, model_name) {
  predictions <- rep(NA_real_, length(y_log))
  coefficient_rows <- vector("list", max(fold_id))
  fit_objects <- vector("list", max(fold_id))

  for (fold in seq_len(max(fold_id))) {
    test_index <- which(fold_id == fold)
    train_index <- which(fold_id != fold)
    x_train <- cbind(`(Intercept)` = 1, x[train_index, , drop = FALSE])
    x_test <- cbind(`(Intercept)` = 1, x[test_index, , drop = FALSE])

    fit <- lm.fit(x = x_train, y = y_log[train_index])
    fold_prediction <- drop(x_test %*% fit$coefficients)
    predictions[test_index] <- fold_prediction

    coefficient_rows[[fold]] <- tibble(
      model = model_name,
      outer_fold = fold,
      term = names(fit$coefficients),
      coefficient = unname(fit$coefficients),
      selected = !is.na(coefficient) & coefficient != 0
    )
    fit_objects[[fold]] <- list(coefficients = fit$coefficients, rank = fit$rank)
  }

  list(
    prediction_log = predictions,
    coefficients = bind_rows(coefficient_rows),
    lambdas = tibble(),
    fits = fit_objects
  )
}

fit_lasso_oof <- function( x, y_log, fold_id, unpenalized_columns, model_name, inner_folds = 10L, seed = 10L) {
  predictions <- rep(NA_real_, length(y_log))
  coefficient_rows <- vector("list", max(fold_id))
  lambda_rows <- vector("list", max(fold_id))
  fit_objects <- vector("list", max(fold_id))
  penalty_factor <- ifelse(colnames(x) %in% unpenalized_columns, 0, 1)

  for (fold in seq_len(max(fold_id))) {
    test_index <- which(fold_id == fold)
    train_index <- which(fold_id != fold)
    inner_fold_id <- make_stratified_folds(
      y_log[train_index], v = inner_folds, seed = seed + fold
    )

    set.seed(seed + fold)
    fit <- cv.glmnet(
      x = x[train_index, , drop = FALSE],
      y = y_log[train_index],
      family = "gaussian",
      alpha = 1,
      nlambda = 200,
      lambda.min.ratio = 0.001,
      standardize = TRUE,
      intercept = TRUE,
      penalty.factor = penalty_factor,
      foldid = inner_fold_id,
      type.measure = "mse"
    )

    predictions[test_index] <- drop(predict(
      fit,
      newx = x[test_index, , drop = FALSE],
      s = "lambda.min"
    ))

    fold_coefficients <- as.matrix(coef(fit, s = "lambda.min"))[, 1]
    coefficient_rows[[fold]] <- tibble(
      model = model_name,
      outer_fold = fold,
      term = names(fold_coefficients),
      coefficient = unname(fold_coefficients),
      selected = coefficient != 0
    )
    lambda_rows[[fold]] <- tibble(
      model = model_name,
      outer_fold = fold,
      lambda_min = fit$lambda.min,
      lambda_1se = fit$lambda.1se,
      n_nonzero_lambda_min = sum(fold_coefficients[-1] != 0),
      n_unpenalized = sum(penalty_factor == 0)
    )
    fit_objects[[fold]] <- fit
  }

  list(
    prediction_log = predictions,
    coefficients = bind_rows(coefficient_rows),
    lambdas = bind_rows(lambda_rows),
    fits = fit_objects
  )
}

y_log <- log(DF_model[[outcome_var]])
basic_columns <- colnames(x_basic)

model_results <- list(
  Basic = fit_ols_oof(x_basic, y_log, outer_fold, "Basic"),
  GenoBMI = fit_ols_oof(x_geno, y_log, outer_fold, "GenoBMI"),
  MetBMI = fit_lasso_oof(
    x_metab, y_log, outer_fold, basic_columns, "MetBMI",
    inner_folds = inner_v, seed = random_seed
  ),
  MicroBMI = fit_lasso_oof(
    x_micro, y_log, outer_fold, basic_columns, "MicroBMI",
    inner_folds = inner_v, seed = random_seed
  ),
  Multi_omics_pBMI = fit_lasso_oof(
    x_multi, y_log, outer_fold, basic_columns, "Multi_omics_pBMI",
    inner_folds = inner_v, seed = random_seed
  )
)

## 7. Combine OOF predictions and calculate performance ---------------------

predictions <- folds
for (model_name in names(model_results)) {
  predictions[[paste0(model_name, "_pred_log")]] <-
    model_results[[model_name]]$prediction_log
  predictions[[paste0(model_name, "_pred_pbmi")]] <-
    exp(model_results[[model_name]]$prediction_log)
}

calculate_metrics <- function(observed, predicted) {
  residual <- observed - predicted
  tibble(
    R2 = cor(observed, predicted)^2,
    R2_predictive = 1 - sum(residual^2) / sum((observed - mean(observed))^2),
    RMSE = sqrt(mean(residual^2)),
    MAE = mean(abs(residual))
  )
}

bootstrap_metric_ci <- function(observed, predicted, reps = 2000L, seed = 10L, conf_level = 0.95) {
  set.seed(seed)
  n <- length(observed)
  boot_metrics <- map_dfr(seq_len(reps), function(i) {
    index <- sample.int(n, size = n, replace = TRUE)
    calculate_metrics(observed[index], predicted[index])
  })
  alpha <- (1 - conf_level) / 2

  boot_metrics %>%
    pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
    group_by(metric) %>%
    summarise(
      ci_lower = quantile(value, alpha, na.rm = TRUE),
      ci_upper = quantile(value, 1 - alpha, na.rm = TRUE),
      .groups = "drop"
    )
}

overall_performance <- map_dfr(seq_along(model_results), function(model_index) {
  model_name <- names(model_results)[model_index]
  predicted <- predictions[[paste0(model_name, "_pred_pbmi")]]
  point_estimates <- calculate_metrics(predictions$observed_pbmi, predicted) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "estimate")
  confidence_intervals <- bootstrap_metric_ci(
    predictions$observed_pbmi,
    predicted,
    reps = bootstrap_reps,
    seed = random_seed + model_index
  )

  left_join(point_estimates, confidence_intervals, by = "metric") %>%
    mutate(
      model = model_name,
      n = nrow(predictions),
      confidence_level = 0.95,
      bootstrap_reps = bootstrap_reps,
      .before = 1
    )
})

fold_performance <- map_dfr(names(model_results), function(model_name) {
  map_dfr(seq_len(outer_v), function(fold) {
    fold_data <- predictions %>% filter(outer_fold == fold)
    calculate_metrics(
      fold_data$observed_pbmi,
      fold_data[[paste0(model_name, "_pred_pbmi")]]
    ) %>% mutate(model = model_name, outer_fold = fold, n = nrow(fold_data), .before = 1)
  })
})

coefficients <- map_dfr(model_results, "coefficients")
lambdas <- map_dfr(model_results, "lambdas")

analysis_summary <- tibble(
  item = c(
    "analysis_population", "outcome", "outcome_model_scale",
    "outer_folds", "inner_folds_lasso", "lambda_rule",
    "metabolites", "bacteria", "fungi", "microbiome_pseudocount",
    "basic_design_columns", "bootstrap_reps", "random_seed"
  ),
  value = as.character(c(
    nrow(DF_model), outcome_var, "log",
    outer_v, inner_v, "lambda.min",
    length(rep_metab123), length(rep_bact74), length(rep_fung6),
    micro_pseudocount, length(basic_columns), bootstrap_reps, random_seed
  ))
)

## 8. Save reproducible model outputs ---------------------------------------

write_xlsx(
  list(
    analysis_summary = analysis_summary,
    overall_performance = overall_performance,
    fold_performance = fold_performance,
    predictions = predictions,
    folds = folds,
    lambdas = lambdas,
    coefficients = coefficients,
    preprocessing = preprocessing_parameters
  ),
  path = file.path(output_dir, "pbmi_model_results.xlsx")
)

write_csv(predictions, file.path(output_dir, "oof_predictions.csv"))
write_csv(folds, file.path(output_dir, "outer_folds.csv"))
write_csv(overall_performance, file.path(output_dir, "overall_performance.csv"))
write_csv(fold_performance, file.path(output_dir, "fold_performance.csv"))
write_csv(lambdas, file.path(output_dir, "lasso_lambdas.csv"))
write_csv(coefficients, file.path(output_dir, "model_coefficients.csv"))
write_csv(
  preprocessing_parameters,
  file.path(output_dir, "preprocessing_parameters.csv")
)

saveRDS(
  list(
    analysis_summary = analysis_summary,
    feature_sets = list(
      covariates = Covars_common,
      metabolites = rep_metab123,
      bacteria = rep_bact74,
      fungi = rep_fung6
    ),
    folds = folds,
    preprocessing_parameters = preprocessing_parameters,
    predictions = predictions,
    performance = list(
      overall = overall_performance,
      by_fold = fold_performance
    ),
    coefficients = coefficients,
    lambdas = lambdas,
    fitted_models = map(model_results, "fits"),
    session_info = sessionInfo()
  ),
  file = file.path(output_dir, "pbmi_model_objects.rds")
)

print(overall_performance)