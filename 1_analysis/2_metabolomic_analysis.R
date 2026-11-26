##############################################################################################
# Analyses for Figure 3
# Characterization of serum metabolomic signatures associated with maternal pre-pregnancy adiposity.
# Metabolite pathway enrichment for Figure 3b was performed using the
# MetaboAnalyst 6.0 web platform. KEGG pathway over-representation analysis
# used the hypergeometric test against the default human KEGG pathway library.
# See plot_fig3.R for the corresponding figure-generation code.
###############################################################################################

# Differential metabolite analyses in THSBC -----------------------------------------
rm(list = ls())

library(tidyverse)
library(purrr)
library(data.table)
library(readxl)
library(writexl)
library(vegan)
library(lmerTest)

source('clean_code/01_help_func.R')

## load data
load('data/01_thsbc_cleaned_data/metabolites_1292_991.RData')
load('data/01_thsbc_cleaned_data/total_pop_5212.RData')
load('data/01_thsbc_cleaned_data/covariates_5212.RData')
load('data/01_thsbc_cleaned_data/wt_omics.rdata')

DF_metab <- total_pop %>%
  inner_join(metab_tab_st, by = "id") %>%
  inner_join(cov_tab, by = c("id", "period")) %>%
  inner_join(wt_metab, by = c("id", "period"))


# Regression of pBMI and metabolites 
Covars <- c("age", "parity", "smk", "drk",  "physic_activity", 
            "fd_grain", "fd_meat", "fd_vege" ,"fd_fruit", "fd_egg", "fd_milkpro", 
            "edu", "medical_insurance", "income",
            "batch", "wk", "antibiotic")

periods <- c("T1", "T2", "T3")

# Standardize pBMI across unique participants rather than repeated observations.
# Coefficients therefore represent the metabolite SD difference per 1-SD pBMI.
pbmi_mean <- total_pop %>%
  distinct(id, pbmi) %>%
  summarise(value = mean(pbmi, na.rm = TRUE)) %>%
  pull(value)

pbmi_sd <- total_pop %>%
  distinct(id, pbmi) %>%
  summarise(value = sd(pbmi, na.rm = TRUE)) %>%
  pull(value)

factor_covars <- c(
  "parity", "smk", "drk", "fd_grain", "fd_meat", "fd_vege",
  "fd_fruit", "fd_egg", "fd_milkpro", "edu", "medical_insurance",
  "income", "batch", "antibiotic"
)

DF_metab <- DF_metab %>%
  mutate(
    pbmi_z = (pbmi - pbmi_mean) / pbmi_sd,
    across(all_of(factor_covars), as.factor)
  )

# Define pairwise pBMI contrasts using the appropriate categorical variable
comparison_config <- tribble(
  ~exposure,                 ~reference,   ~comparison,   ~label,
  "pbmi_owo_normal",         "Normal",     "OWO",         "OWO vs normal weight",
  "pbmi_overweight_normal",  "Normal",     "Overweight",  "Overweight vs normal weight",
  "pbmi_obesity_normal",     "Normal",     "Obese",       "Obesity vs normal weight",
  "pbmi_underweight_normal", "Normal",     "Underweight", "Underweight vs normal weight",
  "pbmi_obesity_overweight", "Overweight", "Obese",       "Obesity vs overweight"
) %>%
  mutate(
    target_population_size = pmap_int(
      list(exposure, reference, comparison),
      function(exposure, reference, comparison) {
        if (exposure == "pbmi_owo_normal") {
          sum(total_pop$pbmi_cut3 %in% c(reference, comparison), na.rm = TRUE)
        } else {
          sum(total_pop$pbmi_cut4 %in% c(reference, comparison), na.rm = TRUE)
        }
      }
    )
  )

DF_metab$pbmi_owo_normal <- factor(
  if_else(
    DF_metab$pbmi_cut3 %in% c("Normal", "OWO"),
    as.character(DF_metab$pbmi_cut3),
    NA_character_
  ),
  levels = c("Normal", "OWO")
)

DF_metab <- comparison_config %>%
  filter(exposure != "pbmi_owo_normal") %>%
  pmap(function(exposure, reference, comparison, label, target_population_size) {
    factor(
      if_else(
        DF_metab$pbmi_cut4 %in% c(reference, comparison),
        as.character(DF_metab$pbmi_cut4),
        NA_character_
      ),
      levels = c(reference, comparison)
    )
  }) %>%
  set_names(comparison_config$exposure[-1]) %>%
  bind_cols(DF_metab, .)


# 1 Linear association per 1-SD higher pBMI 
message("Running continuous pBMI-metabolite association models...")

res_pbmi_continuous <- func_run_omics_models(
  data = DF_metab,
  feature_cols = metab_uni,
  exposure = "pbmi_z",
  covars = Covars,
  periods = periods,
  exposure_label = "pBMI per 1 SD",
  target_population_size = nrow(total_pop)
) %>%
  func_append_feature_info(
    info_data = info_metab,
    info_key = "Index",
    omics_type = "Metabolites"
  )

# 2 Categorical pBMI contrasts 
message("Running categorical pBMI-metabolite association models...")

res_pbmi_categorical <- pmap_dfr(
  comparison_config,
  function(exposure, reference, comparison, label, target_population_size) {
    func_run_omics_models(
      data = DF_metab,
      feature_cols = metab_uni,
      exposure = exposure,
      covars = Covars,
      periods = periods,
      exposure_label = label,
      target_population_size = target_population_size
    )
  }
) %>%
  func_append_feature_info(
    info_data = info_metab,
    info_key = "Index",
    omics_type = "Metabolites"
  )

# 3 Save analysis-ready results
dir.create("result/PART2", recursive = TRUE, showWarnings = FALSE)

save(
  res_pbmi_continuous,
  res_pbmi_categorical,
  comparison_config,
  Covars,
  pbmi_mean,
  pbmi_sd,
  file = "result/PART2/1_diff_metabolites.rdata"
)
message("Differential metabolite analyses completed.")


# Regression of other anthropometric traits and metabolites 
secondary_continuous_config <- tribble(
  ~exposure_raw, ~exposure_z, ~label,
  "wc_t1",       "wc_t1_z",   "Waist circumference per 1 SD",
  "whr_t1",      "whr_t1_z",  "Waist-to-hip ratio per 1 SD",
  "bfr_t1",      "bfr_t1_z",  "Body fat rate per 1 SD"
) %>%
  mutate(
    target_population_size = map_int(
      exposure_raw,
      ~ sum(!is.na(total_pop[[.x]]))
    )
  )

secondary_binary_config <- tribble(
  ~exposure,    ~reference, ~comparison,          ~label,
  "wc_t1_cut",  "Normal",   "Abdominal_obesity", "WC abdominal obesity vs normal",
  "whr_t1_cut", "Normal",   "Abdominal_obesity", "WHR abdominal obesity vs normal",
  "bfr_t1_cut", "Normal",   "Obesity",           "BFR obesity vs normal"
) %>%
  mutate(
    target_population_size = map_int(
      exposure,
      ~ sum(!is.na(total_pop[[.x]]))
    )
  )

for (config_index in seq_len(nrow(secondary_continuous_config))) {
  exposure_raw <- secondary_continuous_config$exposure_raw[[config_index]]
  exposure_z <- secondary_continuous_config$exposure_z[[config_index]]
  exposure_mean <- mean(total_pop[[exposure_raw]], na.rm = TRUE)
  exposure_sd <- sd(total_pop[[exposure_raw]], na.rm = TRUE)
  DF_metab[[exposure_z]] <-
    (DF_metab[[exposure_raw]] - exposure_mean) / exposure_sd
}

DF_metab <- DF_metab %>%
  mutate(
    wc_t1_cut = factor(
      wc_t1_cut,
      levels = c("Normal", "Abdominal_obesity")
    ),
    whr_t1_cut = factor(
      whr_t1_cut,
      levels = c("Normal", "Abdominal_obesity")
    ),
    bfr_t1_cut = factor(
      bfr_t1_cut,
      levels = c("Normal", "Obesity")
    )
  )

message("Running continuous secondary anthropometric-metabolite models...")

res_secondary_metab_continuous <- pmap_dfr(
  secondary_continuous_config,
  function(exposure_raw, exposure_z, label, target_population_size) {
    func_run_omics_models(
      data = DF_metab,
      feature_cols = metab_uni,
      exposure = exposure_z,
      covars = Covars,
      periods = periods,
      exposure_label = label,
      target_population_size = target_population_size
    )
  }
) %>%
  func_append_feature_info(
    info_data = info_metab,
    info_key = "Index",
    omics_type = "Metabolites"
  )

message("Running binary secondary anthropometric-metabolite models...")

res_secondary_metab_binary <- pmap_dfr(
  secondary_binary_config,
  function(exposure, reference, comparison, label, target_population_size) {
    func_run_omics_models(
      data = DF_metab,
      feature_cols = metab_uni,
      exposure = exposure,
      covars = Covars,
      periods = periods,
      exposure_label = label,
      target_population_size = target_population_size
    )
  }
) %>%
  func_append_feature_info(
    info_data = info_metab,
    info_key = "Index",
    omics_type = "Metabolites"
  )

write_xlsx(
  list(
    continuous = res_secondary_metab_continuous,
    binary = res_secondary_metab_binary
  ),
  path = "result/PART2/1-2_diff_metabolites_secdwt.xlsx"
)


# Differential metabolite analyses in HBC-----------------------------------------
## load data
load('data/02_hbc_cleaned_data/metabolites_393_158.RData')
load('data/02_hbc_cleaned_data/total_pop_hbc_393.Rdata')
load('data/02_hbc_cleaned_data/covariates_393.RData')

DF_metab <- total_pop_hbc %>%
  inner_join(metab_tab_hbc_st, by = "id") %>%
  inner_join(cov_tab, by = c("id", "period")) 

# Regression of pBMI and metabolites 

Covars <- c("age", "parity", "smk", "drk",  "physic_activity", 
            "fd_grain", "fd_meat", "fd_vege" ,"fd_fruit", "fd_egg", "fd_milkpro", 
            "edu", "medical_insurance", "income",
            "wk", "antibiotic")

periods <- c("T1", "T2", "T3")
factor_covars <- c(
  "parity", "smk", "drk", "fd_grain", "fd_meat", "fd_vege",
  "fd_fruit", "fd_egg", "fd_milkpro", "edu", "medical_insurance",
  "income", "antibiotic"
)

continuous_config <- tribble(
  ~exposure_raw, ~exposure_z, ~label,
  "pbmi",        "pbmi_z",    "pBMI per 1 SD",
  "wc_t1",       "wc_t1_z",   "Waist circumference per 1 SD",
  "whr_t1",      "whr_t1_z",  "Waist-to-hip ratio per 1 SD"
)

for (config_index in seq_len(nrow(continuous_config))) {
  exposure_raw <- continuous_config$exposure_raw[[config_index]]
  exposure_z <- continuous_config$exposure_z[[config_index]]
  exposure_mean <- mean(total_pop_hbc[[exposure_raw]], na.rm = TRUE)
  exposure_sd <- sd(total_pop_hbc[[exposure_raw]], na.rm = TRUE)
  DF_metab[[exposure_z]] <-
    (DF_metab[[exposure_raw]] - exposure_mean) / exposure_sd
}

pbmi_comparison_config <- tribble(
  ~exposure,                 ~reference,   ~comparison,   ~label,
  "pbmi_owo_normal",         "Normal",     "OWO",         "OWO vs normal weight",
  "pbmi_overweight_normal",  "Normal",     "Overweight",  "Overweight vs normal weight",
  "pbmi_obesity_normal",     "Normal",     "Obese",       "Obesity vs normal weight",
  "pbmi_underweight_normal", "Normal",     "Underweight", "Underweight vs normal weight",
  "pbmi_obesity_overweight", "Overweight", "Obese",       "Obesity vs overweight"
)

DF_metab$pbmi_owo_normal <- factor(
  if_else(
    DF_metab$pbmi_cut3 %in% c("Normal", "OWO"),
    as.character(DF_metab$pbmi_cut3),
    NA_character_
  ),
  levels = c("Normal", "OWO")
)

DF_metab <- pbmi_comparison_config %>%
  filter(exposure != "pbmi_owo_normal") %>%
  pmap(function(exposure, reference, comparison, label) {
    factor(
      if_else(
        DF_metab$pbmi_cut4 %in% c(reference, comparison),
        as.character(DF_metab$pbmi_cut4),
        NA_character_
      ),
      levels = c(reference, comparison)
    )
  }) %>%
  set_names(pbmi_comparison_config$exposure[-1]) %>%
  bind_cols(DF_metab, .)

binary_anthropometric_config <- tribble(
  ~exposure,    ~label,
  "wc_t1_cut",  "WC abdominal obesity vs normal",
  "whr_t1_cut", "WHR abdominal obesity vs normal"
)

DF_metab <- DF_metab %>%
  mutate(
    across(all_of(factor_covars), as.factor),
    wc_t1_cut = factor(
      wc_t1_cut,
      levels = c("Normal", "Abdominal_obesity")
    ),
    whr_t1_cut = factor(
      whr_t1_cut,
      levels = c("Normal", "Abdominal_obesity")
    )
  )

message("Running unweighted HBC continuous exposure-metabolite models...")

res_hbc_metab_continuous <- pmap_dfr(
  continuous_config,
  function(exposure_raw, exposure_z, label) {
    func_run_omics_models_unweighted(
      data = DF_metab,
      feature_cols = metab_uni,
      exposure = exposure_z,
      covars = Covars,
      periods = periods,
      exposure_label = label
    )
  }
) %>%
  func_append_feature_info(info_metab_hbc, "Index", "Metabolites_HBC")

message("Running unweighted HBC categorical pBMI-metabolite models...")

res_hbc_metab_pbmi_categorical <- pmap_dfr(
  pbmi_comparison_config,
  function(exposure, reference, comparison, label) {
    func_run_omics_models_unweighted(
      data = DF_metab,
      feature_cols = metab_uni,
      exposure = exposure,
      covars = Covars,
      periods = periods,
      exposure_label = label
    )
  }
) %>%
  func_append_feature_info(info_metab_hbc, "Index", "Metabolites_HBC")

message("Running unweighted HBC binary anthropometric-metabolite models...")

res_hbc_metab_anthropometric_binary <- pmap_dfr(
  binary_anthropometric_config,
  function(exposure, label) {
    func_run_omics_models_unweighted(
      data = DF_metab,
      feature_cols = metab_uni,
      exposure = exposure,
      covars = Covars,
      periods = periods,
      exposure_label = label
    )
  }
) %>%
  func_append_feature_info(info_metab_hbc, "Index", "Metabolites_HBC")

dir.create("result/PART2", recursive = TRUE, showWarnings = FALSE)

save(
  res_hbc_metab_continuous,
  res_hbc_metab_pbmi_categorical,
  res_hbc_metab_anthropometric_binary,
  continuous_config,
  pbmi_comparison_config,
  binary_anthropometric_config,
  file = "result/PART2/2_diff_metabolites_HBC.rdata"
)

# OPLS-DA of longitudinal metabolites associated with binary adiposity traits ----------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lmerTest)
  library(ropls)
  library(pROC)
  library(PRROC)
  library(writexl)
})

# Configuration 
SEED <- 111L
N_REPEATS <- as.integer(Sys.getenv("OPLS_N_REPEATS", "1000"))
TRAIN_PROP <- 0.75
INNER_FOLDS <- 5L
N_PERM_FINAL <- as.integer(Sys.getenv("OPLS_N_PERM", "1000"))
VIP_THRESHOLD <- 1
VIP_FREQUENCY_THRESHOLD <- 0.70
PREDICTIVE_COMPONENTS <- 1L
ORTHOGONAL_COMPONENTS <- 1L
USE_IPW_FOR_RESIDUALISATION <- TRUE

OUT_DIR <- "result/PART2/4.1_OPLS_DA_longitudinal"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (is.na(N_REPEATS) || N_REPEATS < 2L) {
  stop("OPLS_N_REPEATS must be an integer >= 2.")
}

# Data 
load("data/01_thsbc_cleaned_data/metabolites_1292_991.RData")
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")
load("data/01_thsbc_cleaned_data/wt_omics.rdata")
# Precomputed THSBC metabolite sets used as inputs for downstream analyses.
# For each adiposity contrast, metabolites were selected if the
# Benjamini-Hochberg-adjusted P value was < 0.05 in at least one
# timepoint-specific weighted linear model (T1, T2 or T3) or the pooled
# weighted linear mixed-effects model. Each set is the deduplicated union
# of significant metabolites across these four analyses.
# sigid_metab_owo: OWO versus normal weight (123 metabolites).
# sigid_metab_wc_cut: WC-defined abdominal obesity versus normal WC.
# sigid_metab_bfr_cut: BFR-defined obesity versus normal BFR.
load("result/PART2/1_diff_metabolites_sig_sets_thsbc.rdata")

exposure_config <- tribble(
  ~analysis_key, ~exposure, ~reference, ~comparison, ~metabolite_set,
  "owo",     "pbmi_cut3",  "Normal", "OWO",               "sigid_metab_owo",
  "wc_cut",  "wc_t1_cut", "Normal", "Abdominal_obesity", "sigid_metab_wc_cut",
  "bfr_cut", "bfr_t1_cut","Normal", "Obesity",           "sigid_metab_bfr_cut"
) %>%
  mutate(
    candidate_metabolites = map(
      metabolite_set,
      ~ unique(sig_metab_sets[[.x]])
    )
  )

all_candidate_metabolites <- exposure_config$candidate_metabolites %>%
  unlist(use.names = FALSE) %>%
  unique()
all_exposures <- exposure_config$exposure

covariates <- c(
  "age", "parity", "smk", "drk", "physic_activity",
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
  "edu", "medical_insurance", "income", "batch", "wk", "antibiotic"
)

required_columns <- c(
  "id", "period", all_exposures, "prob", covariates,
  all_candidate_metabolites
)

analysis_long <- total_pop %>%
  select(id, all_of(all_exposures)) %>%
  inner_join(metab_tab_st, by = "id") %>%
  inner_join(cov_tab, by = c("id", "period")) %>%
  inner_join(wt_metab, by = c("id", "period")) %>%
  mutate(
    period = factor(period, levels = c("T1", "T2", "T3"))
  )

missing_columns <- setdiff(required_columns, names(analysis_long))
if (length(missing_columns) > 0L) {
  stop("Required columns are missing: ", paste(missing_columns, collapse = ", "))
}
if (length(exposure_config$candidate_metabolites[[1L]]) != 123L) {
  warning("Expected 123 candidate metabolites, but found ",
          length(exposure_config$candidate_metabolites[[1L]]), ".")
}

complete_metabolomics_ids <- analysis_long %>%
  distinct(id, period) %>%
  count(id, name = "n_periods") %>%
  filter(n_periods == 3L) %>%
  pull(id)

complete_multiomics_ids <- wt_inter747 %>%
  distinct(id, period) %>%
  count(id, name = "n_periods") %>%
  filter(n_periods == 3L) %>%
  pull(id)

# Residualisation 
make_adjustment_formula <- function(outcome, covariate_names) {
  reformulate(
    termlabels = c(covariate_names, "(1 | id)"),
    response = outcome
  )
}

residualise_one_metabolite <- function(data, metabolite, covariate_names, use_ipw = TRUE) {
  model_data <- data %>%
    select(id, period, prob, all_of(covariate_names), all_of(metabolite)) %>%
    drop_na() %>%
    group_by(period) %>%
    # Normalize inverse-probability weights within each period.
    # Multiplication by a common target population size is omitted because it only
    # rescales all weights by the same constant and does not affect model estimates.
    mutate(
      analysis_weight = (1 / prob) / sum(1 / prob)
    ) %>%
    ungroup()

  if (nrow(model_data) == 0L) {
    stop("No complete observations for residualising ", metabolite, ".")
  }

  fit <- if (use_ipw) {
    lmerTest::lmer(
      formula = make_adjustment_formula(metabolite, covariate_names),
      data = model_data,
      weights = analysis_weight,
      REML = FALSE,
      na.action = na.exclude,
      control = lme4::lmerControl(
        check.nobs.vs.nRE = "ignore",
        check.nobs.vs.nlev = "ignore"
      )
    )
  } else {
    lmerTest::lmer(
      formula = make_adjustment_formula(metabolite, covariate_names),
      data = model_data,
      REML = FALSE,
      na.action = na.exclude,
      control = lme4::lmerControl(
        check.nobs.vs.nRE = "ignore",
        check.nobs.vs.nlev = "ignore"
      )
    )
  }

  # Remove fixed covariate effects while retaining participant-specific random
  # intercept variation for subsequent participant-level longitudinal averaging.
  fixed_effect_prediction <- predict(
    fit,
    newdata = model_data,
    re.form = NA
  )
  adjusted_value <-
    model_data[[metabolite]] - as.numeric(fixed_effect_prediction)

  tibble(
    id = model_data$id,
    period = model_data$period,
    metabolite = metabolite,
    adjusted_residual = adjusted_value
  )
}

build_participant_matrix <- function(long_data, metabolite_names, covariate_names, exposure_names, use_ipw = TRUE) {
  residual_long <- map_dfr(
    metabolite_names,
    ~ residualise_one_metabolite(
      data = long_data,
      metabolite = .x,
      covariate_names = covariate_names,
      use_ipw = use_ipw
    )
  )

  residual_wide <- residual_long %>%
    group_by(id, metabolite) %>%
    summarise(
      adjusted_residual_mean = mean(adjusted_residual, na.rm = TRUE),
      n_observed_periods = n_distinct(period),
      .groups = "drop"
    ) %>%
    select(id, metabolite, adjusted_residual_mean) %>%
    pivot_wider(
      names_from = metabolite,
      values_from = adjusted_residual_mean
    )

  participant_info <- long_data %>%
    distinct(id, period, across(all_of(exposure_names))) %>%
    group_by(id, across(all_of(exposure_names))) %>%
    summarise(n_metabolomics_periods = n_distinct(period), .groups = "drop")

  participant_info %>%
    inner_join(residual_wide, by = "id") %>%
    drop_na(all_of(metabolite_names)) %>%
    arrange(id)
}

# Repeated train/test OPLS-DA 
stratified_train_ids <- function(group, train_prop) {
  split(seq_along(group), group) %>%
    map(~ sample(.x, size = floor(length(.x) * train_prop), replace = FALSE)) %>%
    unlist(use.names = FALSE)
}

classification_metrics <- function(observed, predicted_score, threshold = 0.5) {
  observed <- as.integer(observed)
  predicted_class <- if_else(predicted_score >= threshold, 1L, 0L)

  tp <- sum(predicted_class == 1L & observed == 1L)
  tn <- sum(predicted_class == 0L & observed == 0L)
  fp <- sum(predicted_class == 1L & observed == 0L)
  fn <- sum(predicted_class == 0L & observed == 1L)

  sensitivity_value <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
  specificity_value <- if ((tn + fp) > 0L) tn / (tn + fp) else NA_real_

  roc_auc_value <- tryCatch(
    as.numeric(
      pROC::auc(
        response = observed,
        predictor = predicted_score,
        levels = c(0, 1),
        direction = "<",
        quiet = TRUE
      )
    ),
    error = function(e) NA_real_
  )

  pr_auc_value <- tryCatch(
    PRROC::pr.curve(
      scores.class0 = predicted_score[observed == 1L],
      scores.class1 = predicted_score[observed == 0L],
      curve = FALSE
    )$auc.integral,
    error = function(e) NA_real_
  )

  tibble(
    balanced_accuracy = mean(c(sensitivity_value, specificity_value),
                             na.rm = TRUE),
    sensitivity = sensitivity_value,
    specificity = specificity_value,
    roc_auc = roc_auc_value,
    pr_auc = pr_auc_value
  )
}

extract_model_statistics <- function(opls_fit) {
  # Read cumulative statistics for the whole model from the summary table.
  # modelDF also contains component-specific values, which are not the
  # overall R2Y or Q2 of the fitted model.
  summary_df <- as.data.frame(ropls::getSummaryDF(opls_fit))
  required_columns <- c("R2X(cum)", "R2Y(cum)", "Q2(cum)")
  if (nrow(summary_df) != 1L ||
      !all(required_columns %in% names(summary_df))) {
    stop("Expected one OPLS summary row containing R2X(cum), R2Y(cum) and Q2(cum).")
  }

  tibble(
    r2x = as.numeric(summary_df[["R2X(cum)"]][1L]),
    r2y = as.numeric(summary_df[["R2Y(cum)"]][1L]),
    q2 = as.numeric(summary_df[["Q2(cum)"]][1L])
  )
}

fit_one_split <- function(participant_data, metabolite_names, group_var, comparison, split_index, seed) {
  set.seed(seed + split_index)
  group_numeric <- if_else(participant_data[[group_var]] == comparison, 1L, 0L)
  train_index <- stratified_train_ids(group_numeric, TRAIN_PROP)
  test_index <- setdiff(seq_len(nrow(participant_data)), train_index)

  x_train <- as.matrix(participant_data[train_index, metabolite_names])
  x_test <- as.matrix(participant_data[test_index, metabolite_names])
  y_train <- group_numeric[train_index]
  y_test <- group_numeric[test_index]

  fit <- ropls::opls(
    x = x_train,
    y = y_train,
    predI = PREDICTIVE_COMPONENTS,
    orthoI = ORTHOGONAL_COMPONENTS,
    crossvalI = INNER_FOLDS,
    permI = 0,
    scaleC = "standard",
    fig.pdfC = "none",
    info.txtC = "none"
  )

  predicted_score <- as.numeric(predict(fit, x_test))

  performance <- bind_cols(
    tibble(
      iteration = split_index,
      n_train = length(train_index),
      n_test = length(test_index),
      n_train_reference = sum(y_train == 0L),
      n_train_comparison = sum(y_train == 1L),
      n_test_reference = sum(y_test == 0L),
      n_test_comparison = sum(y_test == 1L)
    ),
    classification_metrics(y_test, predicted_score),
    extract_model_statistics(fit)
  )

  vip <- tibble(
    iteration = split_index,
    metabolite = names(fit@vipVn),
    vip = as.numeric(fit@vipVn)
  )

  list(performance = performance, vip = vip)
}

summarise_performance <- function(performance) {
  metric_names <- c(
    "balanced_accuracy", "sensitivity", "specificity", "roc_auc", "pr_auc",
    "q2", "r2x", "r2y"
  )

  performance %>%
    summarise(
      across(
        all_of(metric_names),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          median = ~ median(.x, na.rm = TRUE),
          q025 = ~ quantile(.x, 0.025, na.rm = TRUE),
          q975 = ~ quantile(.x, 0.975, na.rm = TRUE)
        ),
        .names = "{.col}__{.fn}"
      )
    ) %>%
    pivot_longer(
      everything(),
      names_to = c("metric", ".value"),
      names_sep = "__"
    )
}

summarise_vip <- function(vip_long) {
  vip_long %>%
    group_by(metabolite) %>%
    summarise(
      median_vip = median(vip, na.rm = TRUE),
      mean_vip = mean(vip, na.rm = TRUE),
      vip_q025 = quantile(vip, 0.025, na.rm = TRUE),
      vip_q975 = quantile(vip, 0.975, na.rm = TRUE),
      vip_gt1_frequency = mean(vip > VIP_THRESHOLD, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      important_metabolite =
        median_vip > VIP_THRESHOLD &
        vip_gt1_frequency >= VIP_FREQUENCY_THRESHOLD
    ) %>%
    arrange(desc(important_metabolite), desc(median_vip),
            desc(vip_gt1_frequency))
}

run_opls_analysis <- function(participant_data, metabolite_names, analysis_name,
                              group_var, reference, comparison, output_dir,
                              analysis_seed) {
  message("Running ", analysis_name, ": ", nrow(participant_data),
          " participants; ", sum(participant_data[[group_var]] == comparison),
          " ", comparison, " and ",
          sum(participant_data[[group_var]] == reference), " ", reference, ".")

  if (any(table(participant_data[[group_var]]) < 10L)) {
    stop("Too few participants in at least one exposure group for ", analysis_name)
  }

  split_results <- map(
    seq_len(N_REPEATS),
    ~ fit_one_split(
      participant_data, metabolite_names, group_var, comparison, .x,
      analysis_seed
    )
  )

  performance <- map_dfr(split_results, "performance")
  vip_long <- map_dfr(split_results, "vip")
  performance_summary <- summarise_performance(performance)
  vip_summary <- summarise_vip(vip_long)

  x_full <- as.matrix(participant_data[, metabolite_names])
  y_full <- if_else(participant_data[[group_var]] == comparison, 1L, 0L)

  set.seed(analysis_seed)
  final_fit <- ropls::opls(
    x = x_full,
    y = y_full,
    predI = PREDICTIVE_COMPONENTS,
    orthoI = ORTHOGONAL_COMPONENTS,
    crossvalI = INNER_FOLDS,
    permI = N_PERM_FINAL,
    scaleC = "standard",
    fig.pdfC = "none",
    info.txtC = "none"
  )

  final_vip <- tibble(
    metabolite = names(final_fit@vipVn),
    final_full_data_vip = as.numeric(final_fit@vipVn)
  )
  vip_summary <- vip_summary %>%
    left_join(final_vip, by = "metabolite")

  sample_flow <- participant_data %>%
    count(.data[[group_var]], n_metabolomics_periods,
          name = "n_participants") %>%
    rename(exposure_group = all_of(group_var)) %>%
    arrange(exposure_group, n_metabolomics_periods)

  result <- list(
    analysis_name = analysis_name,
    exposure = group_var,
    reference = reference,
    comparison = comparison,
    parameters = list(
      seed = analysis_seed,
      n_repeats = N_REPEATS,
      train_proportion = TRAIN_PROP,
      internal_cv_folds = INNER_FOLDS,
      n_final_permutations = N_PERM_FINAL,
      vip_threshold = VIP_THRESHOLD,
      vip_frequency_threshold = VIP_FREQUENCY_THRESHOLD,
      predictive_components = PREDICTIVE_COMPONENTS,
      orthogonal_components = ORTHOGONAL_COMPONENTS,
      residualisation_covariates = covariates,
      gestational_week_function = "linear",
      ipw_residualisation = USE_IPW_FOR_RESIDUALISATION
    ),
    sample_flow = sample_flow,
    participant_data = participant_data,
    performance = performance,
    performance_summary = performance_summary,
    vip_long = vip_long,
    vip_summary = vip_summary,
    important_metabolites = vip_summary %>%
      filter(important_metabolite) %>%
      pull(metabolite),
    final_model_statistics = extract_model_statistics(final_fit),
    final_model = final_fit
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, file.path(output_dir, paste0(analysis_name, ".rds")))

  write_xlsx(
    list(
      sample_flow = sample_flow,
      performance_summary = performance_summary,
      performance_repeats = performance,
      vip_summary = vip_summary,
      important_metabolites = vip_summary %>% filter(important_metabolite)
    ),
    file.path(output_dir, paste0(analysis_name, ".xlsx"))
  )

  result
}

participant_matrix <- build_participant_matrix(
  long_data = analysis_long,
  metabolite_names = all_candidate_metabolites,
  covariate_names = covariates,
  exposure_names = all_exposures,
  use_ipw = USE_IPW_FOR_RESIDUALISATION
)

run_one_exposure <- function(config_row, config_index) {
  analysis_key <- config_row$analysis_key[[1L]]
  group_var <- config_row$exposure[[1L]]
  reference <- config_row$reference[[1L]]
  comparison <- config_row$comparison[[1L]]
  metabolite_names <- config_row$candidate_metabolites[[1L]]
  exposure_output_dir <- file.path(OUT_DIR, analysis_key)

  if (length(metabolite_names) == 0L) {
    stop("No candidate metabolites were supplied for ", analysis_key, ".")
  }

  exposure_matrix <- participant_matrix %>%
    filter(.data[[group_var]] %in% c(reference, comparison)) %>%
    mutate(
      "{group_var}" := factor(
        .data[[group_var]], levels = c(reference, comparison)
      )
    ) %>%
    drop_na(all_of(metabolite_names))

  analysis_sets <- list(
    primary_available_case = exposure_matrix,
    sensitivity_complete_metabolomics_991 = exposure_matrix %>%
      filter(id %in% complete_metabolomics_ids),
    sensitivity_complete_multiomics_747 = exposure_matrix %>%
      filter(id %in% complete_multiomics_ids)
  )

  analysis_set_summary <- imap_dfr(
    analysis_sets,
    ~ .x %>%
      count(.data[[group_var]], name = "n_participants") %>%
      rename(exposure_group = all_of(group_var)) %>%
      mutate(
        exposure = group_var,
        analysis_key = analysis_key,
        analysis_population = .y,
        .before = 1
      )
  )
  dir.create(exposure_output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(
    analysis_set_summary,
    file.path(exposure_output_dir, "analysis_set_summary.csv")
  )

  exposure_results <- imap(
    analysis_sets,
    ~ run_opls_analysis(
      participant_data = .x,
      metabolite_names = metabolite_names,
      analysis_name = .y,
      group_var = group_var,
      reference = reference,
      comparison = comparison,
      output_dir = exposure_output_dir,
      analysis_seed = SEED + config_index * 10000L
    )
  )

  important_overlap <- list(
    primary = exposure_results$primary_available_case$important_metabolites,
    complete_metabolomics_991 =
      exposure_results$sensitivity_complete_metabolomics_991$important_metabolites,
    complete_multiomics_747 =
      exposure_results$sensitivity_complete_multiomics_747$important_metabolites
  )

  important_overlap_table <- tibble(metabolite = metabolite_names) %>%
    mutate(
      primary = metabolite %in% important_overlap$primary,
      complete_metabolomics_991 =
        metabolite %in% important_overlap$complete_metabolomics_991,
      complete_multiomics_747 =
        metabolite %in% important_overlap$complete_multiomics_747,
      n_analyses_selected = rowSums(across(where(is.logical)))
    ) %>%
    arrange(desc(n_analyses_selected), metabolite)

  write_csv(
    important_overlap_table,
    file.path(exposure_output_dir, "important_metabolite_overlap.csv")
  )

  list(
    config = config_row %>% select(-candidate_metabolites),
    analysis_set_summary = analysis_set_summary,
    results = exposure_results,
    important_overlap_table = important_overlap_table
  )
}

all_results <- map2(
  split(exposure_config, seq_len(nrow(exposure_config))),
  seq_len(nrow(exposure_config)),
  run_one_exposure
)
names(all_results) <- exposure_config$analysis_key

analysis_set_summary <- map_dfr(all_results, "analysis_set_summary")
important_overlap_table <- imap_dfr(
  all_results,
  ~ .x$important_overlap_table %>% mutate(analysis_key = .y, .before = 1)
)

write_csv(analysis_set_summary, file.path(OUT_DIR, "analysis_set_summary.csv"))
write_csv(
  important_overlap_table,
  file.path(OUT_DIR, "important_metabolite_overlap_all_exposures.csv")
)
save(
  all_results,
  analysis_set_summary,
  important_overlap_table,
  file = file.path(OUT_DIR, "4.1_OPLS_DA_longitudinal_results.RData")
)
message("OPLS-DA analyses completed. Results saved to: ", OUT_DIR)

vip_sets <- list(
  owo_vip  = all_results$owo$results$primary_available_case$important_metabolites,
  wc_vip  = all_results$wc_cut$results$primary_available_case$important_metabolites,
  bfr_vip  = all_results$bfr_cut$results$primary_available_case$important_metabolites
)
save(vip_sets, file = file.path(OUT_DIR, "4.1_OPLS_DA_vip_sets.RData"))


# Temporal clustering of OWO-metabolite association coefficients ----------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(latrend)
  library(ggplot2)
  library(ggforce)
  library(ggrepel)
})

## 1. Settings 

random_seed <- 111L
expected_periods <- c("T1", "T2", "T3")
candidate_k <- 2:10
number_of_redrawings <- 50L

result_dir <- "result/PART2/6_metabolites_trend/"
figure_dir <- "figure/PART2/6_metabolites_trend"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

## 2. Load and validate the coefficient data 
load("result/PART2/1_diff_metabolites.rdata")
load("result/PART2/1_diff_metabolites_sig_sets_thsbc.rdata")
load("result/PART2/4.1_OPLS_DA_longitudinal/4.1_OPLS_DA_vip_sets.RData")
load("data/project_col_set.Rdata")

sig_metabolites <- unique(sig_metab_sets$sigid_metab_owo)
owo_vip_metabolites <- unique(vip_sets$owo_vip)

required_columns <- c(
  "feature", "exposure", "period", "beta", "se", "p_value", "p_fdr",
  "status"
)

coefficient_data <- res_pbmi_categorical %>%
  filter(
    exposure == "OWO vs normal weight",
    period %in% expected_periods,
    feature %in% sig_metabolites,
    status == "ok"
  ) %>%
  mutate(
    period = factor(period, levels = expected_periods),
    time = match(as.character(period), expected_periods)
  ) %>%
  arrange(feature, time)

trend_data <- coefficient_data %>%
  transmute(
    ID = feature, time, period,
    value = beta, se,
    lower_ci, upper_ci,
    p_value, p_fdr
  )


## 3. Longitudinal K-means clustering 
options(latrend.id = "ID", latrend.time = "time")
selected_k <- 3L
kml_method <- lcMethodKML(
  response = "value",
  nClusters = selected_k,
  nbRedrawing = number_of_redrawings,
  distanceName = "euclidean",
  scale = FALSE
)
kml_methods <- lcMethods(kml_method, nClusters = candidate_k)

set.seed(random_seed)
coefficient_models <- latrendBatch(kml_methods, data = trend_data)

validity_metric_names <- c("ASW", "BIC", "WMAE", "WRSS")
validity_metric_matrix <- metric(coefficient_models, validity_metric_names)
model_specifications <- as.data.frame(coefficient_models)

cluster_size_by_model <- map_int(
  coefficient_models,
  ~ min(table(trajectoryAssignments(.x)))
)

cluster_validity_metrics <- as_tibble(validity_metric_matrix) %>%
  mutate(
    n_clusters = as.integer(model_specifications$nClusters),
    minimum_cluster_size = cluster_size_by_model,
    minimum_cluster_proportion = minimum_cluster_size /
      length(sig_metabolites),
    .before = 1
  ) %>%
  arrange(n_clusters) %>%
  mutate(
    relative_wmae_reduction = c(NA_real_, -diff(WMAE) / head(WMAE, -1)),
    relative_wrss_reduction = c(NA_real_, -diff(WRSS) / head(WRSS, -1))
  )

# Retain three clusters based on fit metrics, cluster sizes, and interpretability
best_fit <- subset(
  coefficient_models,
  nClusters == selected_k,
  drop = TRUE
)

raw_cluster_assignments <- tibble(
  ID = as.character(ids(best_fit)),
  kml_cluster = as.character(trajectoryAssignments(best_fit))
)

if (nrow(raw_cluster_assignments) != length(sig_metabolites)) {
  stop("The selected model did not assign every included metabolite.")
}

trend_with_raw_clusters <- trend_data %>%
  left_join(raw_cluster_assignments, by = "ID")

raw_cluster_centroids <- trend_with_raw_clusters %>%
  group_by(kml_cluster, period, time) %>%
  summarise(mean_beta = mean(value), .groups = "drop")

# KML cluster letters are arbitrary. Reorder clusters by their T3 mean beta so
# that cluster numbering is reproducible and biologically easier to interpret.
cluster_order <- raw_cluster_centroids %>%
  filter(period == "T3") %>%
  arrange(desc(mean_beta)) %>%
  transmute(
    kml_cluster,
    cluster_number = row_number(),
    cluster = paste0("Cluster ", cluster_number)
  )

classify_coefficient_pattern <- function(beta_t1, beta_t2, beta_t3) {
  betas <- c(beta_t1, beta_t2, beta_t3)

  if (all(betas > 0)) {
    return("Persistently positive")
  }
  if (all(betas < 0)) {
    return("Persistently negative")
  }
  if (beta_t1 > 0 && beta_t3 <= 0) {
    return("Positive-to-null/reversed")
  }
  if (beta_t1 < 0 && beta_t3 >= 0) {
    return("Negative-to-null/reversed")
  }
  "Non-monotonic or mixed"
}

cluster_centroids <- raw_cluster_centroids %>%
  left_join(cluster_order, by = "kml_cluster") %>%
  arrange(cluster_number, time)

cluster_summary <- cluster_centroids %>%
  select(cluster_number, cluster, period, mean_beta) %>%
  pivot_wider(names_from = period, values_from = mean_beta) %>%
  left_join(
    raw_cluster_assignments %>%
      count(kml_cluster, name = "n_metabolites") %>%
      left_join(cluster_order, by = "kml_cluster") %>%
      select(cluster_number, n_metabolites),
    by = "cluster_number"
  ) %>%
  rowwise() %>%
  mutate(
    beta_range = max(c(T1, T2, T3)) - min(c(T1, T2, T3)),
    beta_change_t3_vs_t1 = T3 - T1,
    pattern = classify_coefficient_pattern(T1, T2, T3)
  ) %>%
  ungroup() %>%
  arrange(cluster_number)

annotation_columns <- intersect(
  c("feature", "sugg_analysis", "sugg_cmpd_name", "Compounds",
    "class_shrt", "class1", "class2", "HMDB"),
  names(coefficient_data)
)

metabolite_annotations <- coefficient_data %>%
  select(all_of(annotation_columns)) %>%
  distinct(feature, .keep_all = TRUE)

cluster_assignments <- raw_cluster_assignments %>%
  left_join(cluster_order, by = "kml_cluster") %>%
  rename(feature = ID) %>%
  left_join(metabolite_annotations, by = "feature") %>%
  arrange(cluster_number, feature)

trend_clustered <- trend_data %>%
  rename(feature = ID, beta = value) %>%
  left_join(
    cluster_assignments %>%
      select(feature, cluster_number, cluster),
    by = "feature"
  ) %>%
  mutate(cluster = factor(cluster, levels = paste0("Cluster ", 1:selected_k)))

## 4. Feature-level descriptions of coefficient change 
metabolite_coefficient_patterns <- trend_clustered %>%
  select(feature, cluster_number, cluster, period, beta, se, p_value, p_fdr) %>%
  pivot_wider(
    names_from = period,
    values_from = c(beta, se, p_value, p_fdr),
    names_glue = "{.value}_{period}"
  ) %>%
  rowwise() %>%
  mutate(
    beta_range = max(c(beta_T1, beta_T2, beta_T3)) -
      min(c(beta_T1, beta_T2, beta_T3)),
    beta_change_t3_vs_t1 = beta_T3 - beta_T1,
    same_direction_all_periods =
      all(c(beta_T1, beta_T2, beta_T3) > 0) |
      all(c(beta_T1, beta_T2, beta_T3) < 0),
    coefficient_pattern = classify_coefficient_pattern(
      beta_T1, beta_T2, beta_T3
    )
  ) %>%
  ungroup() %>%
  left_join(metabolite_annotations, by = "feature") %>%
  arrange(cluster_number, desc(beta_range))

cluster_significance_summary <- trend_clustered %>%
  group_by(cluster_number, cluster, period) %>%
  summarise(
    n_metabolites = n_distinct(feature),
    n_fdr_significant = sum(p_fdr < 0.05, na.rm = TRUE),
    proportion_fdr_significant = n_fdr_significant / n_metabolites,
    .groups = "drop"
  )

## 5. Export results 
analysis_set_summary <- tibble(
  n_sig_metabolites_input = length(sig_metabolites),
  coefficient_definition = "Adjusted beta for OWO versus normal weight",
  periods = paste(expected_periods, collapse = ", "),
  candidate_k = paste(range(candidate_k), collapse = "-"),
  selected_k = selected_k,
  number_of_redrawings = number_of_redrawings,
  random_seed = random_seed,
  distance = "Euclidean distance on unscaled standardized beta coefficients",
  selection_rationale = paste(
    "k = 3 retained as an interpretable refinement of k = 2 with balanced",
    "cluster sizes; k >= 4 produced a small cluster at k = 4"
  )
)

write_csv(
  analysis_set_summary,
  file.path(result_dir, "6.2 analysis_set_summary.csv")
)
write_csv(
  coefficient_data,
  file.path(result_dir, "6.2 coefficient_input_data.csv")
)
write_csv(
  cluster_validity_metrics,
  file.path(result_dir, "6.2 cluster_validity_metrics.csv")
)
write_csv(
  cluster_assignments,
  file.path(result_dir, "6.2 metabolite_cluster_assignments.csv")
)
write_csv(
  cluster_centroids,
  file.path(result_dir, "6.2 cluster_centroids.csv")
)
write_csv(
  cluster_summary,
  file.path(result_dir, "6.2 cluster_summary.csv")
)
write_csv(
  metabolite_coefficient_patterns,
  file.path(result_dir, "6.2 metabolite_coefficient_patterns.csv")
)
write_csv(
  cluster_significance_summary,
  file.path(result_dir, "6.2 cluster_significance_summary.csv")
)

save(
  sig_metabolites,
  coefficient_data,
  trend_data,
  coefficient_models,
  best_fit,
  cluster_validity_metrics,
  cluster_assignments,
  cluster_centroids,
  cluster_summary,
  metabolite_coefficient_patterns,
  cluster_significance_summary,
  analysis_set_summary,
  file = file.path(result_dir, "6.2_metabolites_coeffi_trend.RData")
)



# Relative explanatory contribution of pBMI and gestational week ----------------------
rm(list = ls())

library(tidyverse)
library(purrr)
library(data.table)
library(readxl)
library(writexl)
library(lmerTest)
source("clean_code/01_help_func.R")

## Load data 
load("data/01_thsbc_cleaned_data/metabolites_1292_991.RData")
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")
load("data/01_thsbc_cleaned_data/wt_omics.rdata")

DF_metab <- total_pop %>%
  inner_join(metab_tab_st, by = "id") %>%
  inner_join(cov_tab, by = c("id", "period")) %>%
  inner_join(wt_metab, by = c("id", "period"))

Covars <- c("age", "parity", "smk", "drk", "physic_activity",
            "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
            "edu", "medical_insurance", "income", "batch", "antibiotic")

DF_metab <- DF_metab %>%
  mutate(
    pbmi_z = (pbmi - mean(total_pop$pbmi, na.rm = TRUE)) /
      sd(total_pop$pbmi, na.rm = TRUE),
    across(all_of(c("parity", "smk", "drk", "fd_grain", "fd_meat",
                    "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro", "edu",
                    "medical_insurance", "income", "batch", "antibiotic")),
           as.factor),
    pbmi_owo = factor(
      if_else(pbmi_cut3 %in% c("Normal", "OWO"),
              as.character(pbmi_cut3), NA_character_),
      levels = c("Normal", "OWO")
    )
  )

target_owo_normal <- sum(
  total_pop$pbmi_cut3 %in% c("Normal", "OWO"), na.rm = TRUE
)

## Four nested models plus the interaction model 
variance_results <- purrr::map_dfr(
  metab_uni,
  ~ func_omics_variance_contribution(
    data = DF_metab,
    feature = .x,
    pbmi_var = "pbmi_owo",
    wk_var = "wk",
    covars = Covars,
    target_population_size = target_owo_normal
  )
) %>%
  mutate(
    p_delta_r2_pbmi_lrt_fdr = p.adjust(p_delta_r2_pbmi_lrt, method = "fdr"),
    p_delta_r2_wk_lrt_fdr = p.adjust(p_delta_r2_wk_lrt, method = "fdr"),
    p_interaction_fdr = p.adjust(p_interaction_lrt, method = "fdr"),
    p_interaction_coef_fdr = p.adjust(p_interaction_coef, method = "fdr")
  ) %>%
  left_join(info_metab, by = c("feature" = "Index")) %>%
  arrange(p_interaction_lrt)

## Period-specific survey-weighted models 
# All five models for a feature-period combination are fitted to exactly the
# same complete-case sample:
# M0            = covariates only
# M_pBMI        = M0 + pBMI
# M_wk          = M0 + linear gestational week
# M_both        = M0 + pBMI + linear gestational week
# M_interaction = M0 + pBMI * linear gestational week
func_omics_period_variance_contribution <- function(
    data,
    feature,
    period,
    pbmi_var = "pbmi_owo",
    wk_var = "wk",
    covars,
    id_col = "id",
    period_col = "period",
    probability_col = "prob",
    target_population_size = 5212,
    minimum_observations = 20) {
  required_variables <- unique(c(
    feature, pbmi_var, wk_var, covars,
    id_col, period_col, probability_col
  ))
  missing_variables <- setdiff(required_variables, names(data))

  if (length(missing_variables) > 0L) {
    stop(
      "Missing variables for ", feature, "/", period, ": ",
      paste(missing_variables, collapse = ", "), "."
    )
  }
  if (wk_var %in% covars || pbmi_var %in% covars) {
    stop("`covars` must not contain `pbmi_var` or `wk_var`.")
  }

  analysis_data <- data %>%
    filter(.data[[period_col]] == .env$period) %>%
    filter(if_all(all_of(required_variables), ~ !is.na(.x)))

  empty_result <- tibble(
    feature = feature,
    period = period,
    model = "cross_sectional",
    n_obs = nrow(analysis_data),
    n_id = n_distinct(analysis_data[[id_col]]),
    target_population_size = target_population_size,
    r2_m0 = NA_real_,
    r2_m_pbmi = NA_real_,
    r2_m_wk = NA_real_,
    r2_m_both = NA_real_,
    delta_r2_pbmi = NA_real_,
    delta_r2_wk = NA_real_,
    shapley_pbmi = NA_real_,
    shapley_wk = NA_real_,
    r2_m_interaction = NA_real_,
    delta_r2_interaction = NA_real_,
    p_delta_r2_pbmi_wald = NA_real_,
    p_delta_r2_wk_wald = NA_real_,
    p_interaction_wald = NA_real_,
    interaction_term = NA_character_,
    interaction_beta = NA_real_,
    interaction_se = NA_real_,
    p_interaction_coef = NA_real_,
    test_method = "survey-weighted Wald F test",
    status = "insufficient_data"
  )

  if (nrow(analysis_data) < minimum_observations ||
      n_distinct(analysis_data[[pbmi_var]]) < 2L ||
      n_distinct(analysis_data[[wk_var]]) < 2L) {
    return(empty_result)
  }

  tryCatch({
    if (any(!is.finite(analysis_data[[feature]])) ||
        any(!is.finite(analysis_data[[wk_var]])) ||
        any(!is.finite(analysis_data[[probability_col]])) ||
        any(analysis_data[[probability_col]] <= 0) ||
        any(analysis_data[[probability_col]] > 1)) {
      stop("Non-finite feature/week values or invalid inclusion probabilities.")
    }

    analysis_data$.sampling_probability <-
      analysis_data[[probability_col]]
    analysis_data$.target_population_fpc <- target_population_size

    survey_design <- survey::svydesign(
      id = ~1,
      prob = ~.sampling_probability,
      fpc = ~.target_population_fpc,
      data = analysis_data
    )

    covariate_rhs <- if (length(covars) == 0L) {
      "1"
    } else {
      paste(covars, collapse = " + ")
    }

    fit_weighted_model <- function(fixed_rhs) {
      survey::svyglm(
        formula = as.formula(paste0(feature, " ~ ", fixed_rhs)),
        design = survey_design,
        family = gaussian()
      )
    }

    fit_m0 <- fit_weighted_model(covariate_rhs)
    fit_m_pbmi <- fit_weighted_model(
      paste(pbmi_var, covariate_rhs, sep = " + ")
    )
    fit_m_wk <- fit_weighted_model(
      paste(wk_var, covariate_rhs, sep = " + ")
    )
    fit_m_both <- fit_weighted_model(
      paste(pbmi_var, wk_var, covariate_rhs, sep = " + ")
    )
    fit_m_interaction <- fit_weighted_model(
      paste0(pbmi_var, " * ", wk_var, " + ", covariate_rhs)
    )

    sampling_weights <- 1 / analysis_data[[probability_col]]
    response <- analysis_data[[feature]]
    response_mean <- weighted.mean(response, sampling_weights)
    total_ss <- sum(sampling_weights * (response - response_mean)^2)

    weighted_r2 <- function(model) {
      if (!is.finite(total_ss) || total_ss <= 0) {
        return(NA_real_)
      }
      residual_ss <- sum(
        sampling_weights * residuals(model, type = "response")^2
      )
      1 - residual_ss / total_ss
    }

    extract_wald_p <- function(model, term) {
      test_result <- survey::regTermTest(
        model,
        as.formula(paste("~", term)),
        method = "Wald"
      )
      as.numeric(test_result$p)
    }

    r2_m0 <- weighted_r2(fit_m0)
    r2_m_pbmi <- weighted_r2(fit_m_pbmi)
    r2_m_wk <- weighted_r2(fit_m_wk)
    r2_m_both <- weighted_r2(fit_m_both)
    r2_m_interaction <- weighted_r2(fit_m_interaction)

    interaction_coefficient <- summary(fit_m_interaction)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      filter(str_detect(term, fixed(":")))

    interaction_term <- if (nrow(interaction_coefficient) == 0L) {
      NA_character_
    } else {
      interaction_coefficient$term[[1L]]
    }
    interaction_beta <- if (nrow(interaction_coefficient) == 0L) {
      NA_real_
    } else {
      interaction_coefficient$Estimate[[1L]]
    }
    interaction_se <- if (nrow(interaction_coefficient) == 0L) {
      NA_real_
    } else {
      interaction_coefficient$`Std. Error`[[1L]]
    }
    interaction_p_value <- if (nrow(interaction_coefficient) == 0L) {
      NA_real_
    } else {
      interaction_coefficient$`Pr(>|t|)`[[1L]]
    }

    tibble(
      feature = feature,
      period = period,
      model = "cross_sectional",
      n_obs = nrow(analysis_data),
      n_id = n_distinct(analysis_data[[id_col]]),
      target_population_size = target_population_size,
      r2_m0 = r2_m0,
      r2_m_pbmi = r2_m_pbmi,
      r2_m_wk = r2_m_wk,
      r2_m_both = r2_m_both,
      delta_r2_pbmi = r2_m_both - r2_m_wk,
      delta_r2_wk = r2_m_both - r2_m_pbmi,
      shapley_pbmi = 0.5 * (
        (r2_m_pbmi - r2_m0) + (r2_m_both - r2_m_wk)
      ),
      shapley_wk = 0.5 * (
        (r2_m_wk - r2_m0) + (r2_m_both - r2_m_pbmi)
      ),
      r2_m_interaction = r2_m_interaction,
      delta_r2_interaction = r2_m_interaction - r2_m_both,
      p_delta_r2_pbmi_wald = extract_wald_p(fit_m_both, pbmi_var),
      p_delta_r2_wk_wald = extract_wald_p(fit_m_both, wk_var),
      p_interaction_wald = extract_wald_p(
        fit_m_interaction,
        paste0(pbmi_var, ":", wk_var)
      ),
      interaction_term = interaction_term,
      interaction_beta = interaction_beta,
      interaction_se = interaction_se,
      p_interaction_coef = interaction_p_value,
      test_method = "survey-weighted Wald F test",
      status = "ok"
    )
  }, error = function(e) {
    empty_result %>%
      mutate(status = paste0("error: ", conditionMessage(e)))
  })
}

periods <- c("T1", "T2", "T3")

period_variance_results <- tidyr::crossing(
  feature = metab_uni,
  period = periods
) %>%
  pmap_dfr(
    ~ func_omics_period_variance_contribution(
      data = DF_metab,
      feature = ..1,
      period = ..2,
      pbmi_var = "pbmi_owo",
      wk_var = "wk",
      covars = Covars,
      target_population_size = target_owo_normal
    )
  ) %>%
  group_by(period) %>%
  mutate(
    p_delta_r2_pbmi_wald_fdr = p.adjust(
      p_delta_r2_pbmi_wald, method = "fdr"
    ),
    p_delta_r2_wk_wald_fdr = p.adjust(
      p_delta_r2_wk_wald, method = "fdr"
    ),
    p_interaction_wald_fdr = p.adjust(
      p_interaction_wald, method = "fdr"
    ),
    p_interaction_coef_fdr = p.adjust(
      p_interaction_coef, method = "fdr"
    )
  ) %>%
  ungroup() %>%
  left_join(info_metab, by = c("feature" = "Index")) %>%
  arrange(period, p_interaction_wald, feature)

period_analysis_summary <- period_variance_results %>%
  group_by(period) %>%
  summarise(
    n_features = n(),
    n_successful = sum(status == "ok"),
    n_failed = sum(status != "ok"),
    n_delta_r2_pbmi_gt_wk = sum(
      delta_r2_pbmi > delta_r2_wk, na.rm = TRUE
    ),
    n_pbmi_fdr = sum(p_delta_r2_pbmi_wald_fdr < 0.05, na.rm = TRUE),
    n_wk_fdr = sum(p_delta_r2_wk_wald_fdr < 0.05, na.rm = TRUE),
    n_interaction_fdr = sum(
      p_interaction_wald_fdr < 0.05, na.rm = TRUE
    ),
    median_delta_r2_pbmi = median(delta_r2_pbmi, na.rm = TRUE),
    median_delta_r2_wk = median(delta_r2_wk, na.rm = TRUE),
    .groups = "drop"
  )

mixed_analysis_summary <- variance_results %>%
  summarise(
    n_features = n(),
    n_successful = sum(status == "ok"),
    n_failed = sum(status != "ok"),
    n_singular_m_wk = sum(singular_m_wk %in% TRUE, na.rm = TRUE),
    n_singular_m_both = sum(singular_m_both %in% TRUE, na.rm = TRUE),
    median_icc_m_wk = median(icc_m_wk, na.rm = TRUE),
    median_icc_m_both = median(icc_m_both, na.rm = TRUE)
  )

dir.create("result/PART2", recursive = TRUE, showWarnings = FALSE)
save(
  variance_results,
  period_variance_results,
  mixed_analysis_summary,
  period_analysis_summary,
  file = "result/PART2/3.3_variance_icc.rdata"
)
message("Relative pBMI and gestational-week variance analysis completed.")

# Clinical relevance in THSBC----------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(survey)
  library(writexl)
})

# Input 
load("data/01_thsbc_cleaned_data/metabolites_1292_991.RData")
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")
load("data/01_thsbc_cleaned_data/wt_omics.rdata")
load("data/01_thsbc_cleaned_data/clinical_traits_t1t2_4562_4874.RData")
load("result/PART2/4.1_OPLS_DA_longitudinal/4.1_OPLS_DA_vip_sets.RData")

owo_metabolites <- unique(vip_sets$owo_vip)
target_population_size <- 5212L

clini_var_t1 <- c(
  "SBP", "DBP", "ALB", "TP", "a_HBDH", "A/G Ratio", "CHE", 
  "GGT", "ALT", "AST", "IBIL", "ALP", "GLB", "LAC", "LDH", "ADA",
  "DBIL", "TBIL", "Cr", "Ccr", "BUN", "UA", "CysC", "FPG",
  "FT3", "FT4", "Tg", "TPOAb", "TSH"
)
clini_var_t2 <- c("OGTT_0h", "OGTT_1h", "OGTT_2h")

covariates <- c(
  "age", "parity", "smk", "drk", "physic_activity",
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
  "edu", "medical_insurance", "income", "batch", "wk", "antibiotic"
)

# Build analysis datasets 
make_metabolite_period_data <- function(metabolite_period) {
  metab_tab_st %>%
    select(id, period, batch, all_of(owo_metabolites)) %>%
    inner_join(cov_tab,by = c("id", "period")) %>%
    inner_join(wt_metab,by = c("id", "period")) %>%
    filter(period == metabolite_period) %>%
    mutate(metabolite_period = metabolite_period) %>%
    select(-period)
}

add_metabolite_pc1_score <- function(data, metabolite_names,
                                     metabolite_period) {
  complete_index <- complete.cases(data[, metabolite_names, drop = FALSE])

  metabolite_matrix <- as.matrix(
    data[complete_index, metabolite_names, drop = FALSE]
  )

  # Metabolites have already been log-transformed and standardized upstream.
  pc_fit <- prcomp(
    metabolite_matrix,
    center = FALSE,
    scale. = FALSE
  )

  score_data <- tibble(
    id = data$id[complete_index],
    metabolite_score = as.numeric(pc_fit$x[, 1])
  ) %>%
    left_join(
      total_pop %>% select(id, pbmi_cut3),
      by = "id"
    )

  group_means <- score_data %>%
    filter(pbmi_cut3 %in% c("Normal", "OWO")) %>%
    group_by(pbmi_cut3) %>%
    summarise(mean_score = mean(metabolite_score, na.rm = TRUE),
              .groups = "drop")

  normal_mean <- group_means %>%
    filter(pbmi_cut3 == "Normal") %>%
    pull(mean_score)
  owo_mean <- group_means %>%
    filter(pbmi_cut3 == "OWO") %>%
    pull(mean_score)

  score_sign <- if (
    length(normal_mean) == 1L && length(owo_mean) == 1L &&
      is.finite(normal_mean) && is.finite(owo_mean) && owo_mean < normal_mean
  ) -1 else 1

  score_data <- score_data %>%
    transmute(
      id,
      metabolite_score = as.numeric(
        scale(score_sign * metabolite_score)
      )
    )

  pc1_loadings <- tibble(
    metabolite = rownames(pc_fit$rotation),
    loading = score_sign * pc_fit$rotation[, 1],
    metabolite_period = metabolite_period
  )

  pc1_summary <- tibble(
    metabolite_period = metabolite_period,
    n_participants = sum(complete_index),
    n_metabolites = length(metabolite_names),
    variance_explained =
      pc_fit$sdev[1]^2 / sum(pc_fit$sdev^2),
    score_standardization = "Mean 0 and SD 1 within metabolite period",
    direction_rule = "Mean PC1 in OWO >= mean PC1 in Normal",
    sign_multiplier = score_sign
  )

  list(
    data = data %>% left_join(score_data, by = "id"),
    loadings = pc1_loadings,
    summary = pc1_summary
  )
}

metabolites_t1_pc <- make_metabolite_period_data("T1") %>%
  add_metabolite_pc1_score(owo_metabolites, "T1")
metabolites_t2_pc <- make_metabolite_period_data("T2") %>%
  add_metabolite_pc1_score(owo_metabolites, "T2")

metabolites_t1 <- metabolites_t1_pc$data
metabolites_t2 <- metabolites_t2_pc$data

pc1_loadings <- bind_rows(
  metabolites_t1_pc$loadings,
  metabolites_t2_pc$loadings
)
pc1_summary <- bind_rows(
  metabolites_t1_pc$summary,
  metabolites_t2_pc$summary
)

analysis_t1_metab_t1_clinical <- metabolites_t1 %>%
  inner_join(clini_traits_T1, by = "id")

analysis_t1_metab_t2_clinical <- metabolites_t1 %>%
  inner_join(clini_traits_T2, by = "id")

analysis_t2_metab_t2_clinical <- metabolites_t2 %>%
  inner_join(clini_traits_T2, by = "id")

# Weighted linear regression
fit_one_association <- function(data, clinical_trait, metabolite,
                                clinical_period, metabolite_period,
                                target_n = target_population_size) {
  required_columns <- unique(c(
    "id", "prob", clinical_trait, metabolite, covariates
  ))

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Variables absent from analysis data: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  model_data <- data %>%
    select(all_of(required_columns)) %>%
    filter(if_all(everything(), ~ !is.na(.x))) %>%
    mutate(
      .clinical_outcome = .data[[clinical_trait]],
      .metabolite_exposure = .data[[metabolite]]
    )

  empty_result <- function(status_text) {
    tibble(
      clinical_trait = clinical_trait,
      metabolite = metabolite,
      clinical_period = clinical_period,
      metabolite_period = metabolite_period,
      n_obs = nrow(model_data),
      beta = NA_real_,
      se = NA_real_,
      lower_ci = NA_real_,
      upper_ci = NA_real_,
      p_value = NA_real_,
      status = status_text
    )
  }

  if (nrow(model_data) < 20L) {
    return(empty_result("insufficient_data"))
  }
  if (n_distinct(model_data[[clinical_trait]]) < 2L ||
      n_distinct(model_data[[metabolite]]) < 2L) {
    return(empty_result("no_variation"))
  }
  if (any(model_data$prob <= 0 | model_data$prob > 1)) {
    return(empty_result("invalid_sampling_probability"))
  }

  model_formula <- reformulate(
    termlabels = c(".metabolite_exposure", covariates),
    response = ".clinical_outcome"
  )

  tryCatch({
    model_data$.sampling_probability <- model_data$prob
    model_data$.target_population_fpc <- target_n

    survey_design <- survey::svydesign(
      id = ~1,
      prob = ~.sampling_probability,
      fpc = ~.target_population_fpc,
      data = model_data
    )

    fit <- survey::svyglm(
      stats::as.formula(model_formula),
      design = survey_design,
      family = stats::gaussian()
    )

    coefficient_table <- summary(fit)$coefficients
    coefficient_row <- match(
      ".metabolite_exposure",
      rownames(coefficient_table)
    )

    if (is.na(coefficient_row)) {
      return(empty_result("metabolite_coefficient_not_estimable"))
    }

    estimate <- coefficient_table[coefficient_row, "Estimate"]
    standard_error <- coefficient_table[coefficient_row, "Std. Error"]
    p_column <- grep("Pr\\(", colnames(coefficient_table), value = TRUE)[1L]

    tibble(
      clinical_trait = clinical_trait,
      metabolite = metabolite,
      clinical_period = clinical_period,
      metabolite_period = metabolite_period,
      n_obs = nrow(model_data),
      beta = unname(estimate),
      se = unname(standard_error),
      lower_ci = unname(estimate - qnorm(0.975) * standard_error),
      upper_ci = unname(estimate + qnorm(0.975) * standard_error),
      p_value = unname(coefficient_table[coefficient_row, p_column]),
      status = "ok"
    )
  }, error = function(error_condition) {
    empty_result(paste0("error: ", conditionMessage(error_condition)))
  })
}

run_association_grid <- function(data, clinical_traits, clinical_period,
                                 metabolite_period,
                                 metabolite_features = owo_metabolites) {
  tidyr::crossing(
    clinical_trait = clinical_traits,
    metabolite = metabolite_features
  ) %>%
    mutate(
      result = map2(
        clinical_trait,
        metabolite,
        ~ fit_one_association(
          data = data,
          clinical_trait = .x,
          metabolite = .y,
          clinical_period = clinical_period,
          metabolite_period = metabolite_period
        )
      )
    ) %>%
    select(result) %>%
    unnest(result) %>%
    group_by(clinical_period, metabolite_period, clinical_trait) %>%
    mutate(
      p_fdr = p.adjust(p_value, method = "fdr"),
      p_bonferroni = p.adjust(p_value, method = "bonferroni")
    ) %>%
    ungroup() %>%
    arrange(clinical_period, metabolite_period, clinical_trait, p_value)
}

# Analyses 
message("Running T1 metabolites -> T1 clinical traits...")
results_t1_metab_t1_clinical <- run_association_grid(
  data = analysis_t1_metab_t1_clinical,
  clinical_traits = clini_var_t1,
  clinical_period = "T1",
  metabolite_period = "T1"
)

message("Running T1 metabolites -> T2 OGTT traits...")
results_t1_metab_t2_clinical <- run_association_grid(
  data = analysis_t1_metab_t2_clinical,
  clinical_traits = clini_var_t2,
  clinical_period = "T2",
  metabolite_period = "T1"
)

message("Running T1 metabolite PC1 score -> T1 clinical traits...")
results_t1_score_t1_clinical <- run_association_grid(
  data = analysis_t1_metab_t1_clinical,
  clinical_traits = clini_var_t1,
  clinical_period = "T1",
  metabolite_period = "T1",
  metabolite_features = "metabolite_score"
)

message("Running T1 metabolite PC1 score -> T2 OGTT traits...")
results_t1_score_t2_clinical <- run_association_grid(
  data = analysis_t1_metab_t2_clinical,
  clinical_traits = clini_var_t2,
  clinical_period = "T2",
  metabolite_period = "T1",
  metabolite_features = "metabolite_score"
)

score_association_results <- bind_rows(
  results_t1_score_t1_clinical,
  results_t1_score_t2_clinical
) %>%
  mutate(feature_type = "PC1 score", .before = metabolite)

individual_metabolite_results <- bind_rows(
  results_t1_metab_t1_clinical,
  results_t1_metab_t2_clinical
) %>%
  mutate(feature_type = "Individual metabolite", .before = metabolite)

association_results <- bind_rows(
  individual_metabolite_results,
  score_association_results
)

analysis_sample_summary <- bind_rows(
  tibble(
    clinical_period = "T1",
    metabolite_period = "T1",
    n_joined = n_distinct(analysis_t1_metab_t1_clinical$id)
  ),
  tibble(
    clinical_period = "T2",
    metabolite_period = "T1",
    n_joined = n_distinct(analysis_t1_metab_t2_clinical$id)
  )
)

status_summary <- association_results %>%
  count(clinical_period, metabolite_period, status, name = "n_models")

# Output
dir.create("result/PART2", recursive = TRUE, showWarnings = FALSE)
save(
  association_results,
  results_t1_metab_t1_clinical,
  results_t1_metab_t2_clinical,
  results_t1_score_t1_clinical,
  results_t1_score_t2_clinical,
  score_association_results,
  pc1_summary,
  pc1_loadings,
  analysis_sample_summary,
  status_summary,
  file = "result/PART2/5.1_association_metab_clinical_traits_thsbc.RData"
)

# Clinical relevance in HBC --------------------------------------------
rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(writexl)
})

# Input 
load("data/02_hbc_cleaned_data/metabolites_393_158.RData")
load("data/02_hbc_cleaned_data/total_pop_hbc_393.Rdata")
load("data/02_hbc_cleaned_data/covariates_393.RData")
load("data/02_hbc_cleaned_data/clinical_traits_t1t2_393.RData")
load("result/PART2/4.1_OPLS_DA_longitudinal/4.1_OPLS_DA_vip_sets.RData")

owo_metabolites <- info_metab_hbc$Index[
  which(info_metab_hbc$Index_xh %in% unique(vip_sets$owo_vip))
] %>%
  unique()

stopifnot(
  length(owo_metabolites) == 27L,
  all(owo_metabolites %in% names(metab_tab_hbc_st))
)

metabolite_annotation <- info_metab_hbc %>%
  filter(Index %in% owo_metabolites) %>%
  transmute(
    metabolite = Index,
    thsbc_metabolite = Index_xh,
    metabolite_name = sugg_cmpd_name
  ) %>%
  distinct(metabolite, .keep_all = TRUE)

clini_var_t1_hbc <- c(
  "SBP", "DBP", "FPG",
  "ALP", "ALT", "GGT", "AST", "DBIL", "IBIL", "A/G Ratio", "TBIL",
  "GLB", "ALB", "TP", "LDH", "ADA", "Cr", "UA", "CysC",
  "FT4", "TPOAb", "TSH"
)
clini_var_t2_hbc <- c("OGTT_0h", "OGTT_1h", "OGTT_2h")

covariates <- c(
  "age", "parity", "smk", "drk", "physic_activity",
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
  "edu", "medical_insurance", "income", "wk", "antibiotic"
)

factor_covariates <- c(
  "parity", "smk", "drk", "fd_grain", "fd_meat", "fd_vege",
  "fd_fruit", "fd_egg", "fd_milkpro", "edu", "medical_insurance",
  "income", "antibiotic"
)

# Build analysis datasets 
make_metabolite_period_data <- function(metabolite_period) {
  metab_tab_hbc_st %>%
    select(id, period, wk, all_of(owo_metabolites)) %>%
    inner_join(cov_tab, by = c("id", "period")) %>%
    filter(period == metabolite_period) %>%
    mutate(
      across(any_of(factor_covariates), as.factor),
      metabolite_period = metabolite_period
    ) %>%
    select(-period)
}

add_metabolite_pc1_score <- function(data, metabolite_names,
                                     metabolite_period) {
  complete_index <- complete.cases(data[, metabolite_names, drop = FALSE])
  metabolite_matrix <- as.matrix(
    data[complete_index, metabolite_names, drop = FALSE]
  )

  pc_fit <- prcomp(
    metabolite_matrix,
    center = FALSE,
    scale. = FALSE
  )

  score_data <- tibble(
    id = data$id[complete_index],
    metabolite_score = as.numeric(pc_fit$x[, 1])
  ) %>%
    left_join(
      total_pop_hbc %>% select(id, pbmi_cut3),
      by = "id"
    )

  group_means <- score_data %>%
    filter(pbmi_cut3 %in% c("Normal", "OWO")) %>%
    group_by(pbmi_cut3) %>%
    summarise(
      mean_score = mean(metabolite_score, na.rm = TRUE),
      .groups = "drop"
    )

  normal_mean <- group_means %>%
    filter(pbmi_cut3 == "Normal") %>%
    pull(mean_score)
  owo_mean <- group_means %>%
    filter(pbmi_cut3 == "OWO") %>%
    pull(mean_score)

  score_sign <- if (
    length(normal_mean) == 1L && length(owo_mean) == 1L &&
      is.finite(normal_mean) && is.finite(owo_mean) && owo_mean < normal_mean
  ) -1 else 1

  score_data <- score_data %>%
    transmute(
      id,
      metabolite_score = as.numeric(scale(score_sign * metabolite_score))
    )

  pc1_loadings <- tibble(
    metabolite = rownames(pc_fit$rotation),
    loading = score_sign * pc_fit$rotation[, 1],
    metabolite_period = metabolite_period
  ) %>%
    left_join(metabolite_annotation, by = "metabolite")

  pc1_summary <- tibble(
    metabolite_period = metabolite_period,
    n_participants = sum(complete_index),
    n_metabolites = length(metabolite_names),
    variance_explained = pc_fit$sdev[1]^2 / sum(pc_fit$sdev^2),
    score_standardization = "Mean 0 and SD 1 within metabolite period",
    direction_rule = "Mean PC1 in OWO >= mean PC1 in Normal",
    sign_multiplier = score_sign
  )

  list(
    data = data %>% left_join(score_data, by = "id"),
    loadings = pc1_loadings,
    summary = pc1_summary
  )
}

metabolites_t1_pc <- make_metabolite_period_data("T1") %>%
  add_metabolite_pc1_score(owo_metabolites, "T1")
metabolites_t2_pc <- make_metabolite_period_data("T2") %>%
  add_metabolite_pc1_score(owo_metabolites, "T2")

metabolites_t1 <- metabolites_t1_pc$data
metabolites_t2 <- metabolites_t2_pc$data

pc1_loadings <- bind_rows(
  metabolites_t1_pc$loadings,
  metabolites_t2_pc$loadings
)
pc1_summary <- bind_rows(
  metabolites_t1_pc$summary,
  metabolites_t2_pc$summary
)

analysis_t1_metab_t1_clinical <- metabolites_t1 %>%
  inner_join(clini_traits_T1 %>% select(-period), by = "id")

analysis_t1_metab_t2_clinical <- metabolites_t1 %>%
  inner_join(clini_traits_T2 %>% select(-period), by = "id")

# Ordinary linear regression 
fit_one_association <- function(data, clinical_trait, metabolite,
                                clinical_period, metabolite_period) {
  required_columns <- unique(c("id", clinical_trait, metabolite, covariates))

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Variables absent from analysis data: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  model_data <- data %>%
    select(all_of(required_columns)) %>%
    filter(if_all(everything(), ~ !is.na(.x))) %>%
    mutate(
      .clinical_outcome = .data[[clinical_trait]],
      .metabolite_exposure = .data[[metabolite]]
    )

  empty_result <- function(status_text) {
    tibble(
      clinical_trait = clinical_trait,
      metabolite = metabolite,
      clinical_period = clinical_period,
      metabolite_period = metabolite_period,
      n_obs = nrow(model_data),
      beta = NA_real_,
      se = NA_real_,
      lower_ci = NA_real_,
      upper_ci = NA_real_,
      p_value = NA_real_,
      status = status_text
    )
  }

  if (nrow(model_data) < 20L) {
    return(empty_result("insufficient_data"))
  }
  if (n_distinct(model_data[[clinical_trait]]) < 2L ||
      n_distinct(model_data[[metabolite]]) < 2L) {
    return(empty_result("no_variation"))
  }

  estimable_covariates <- covariates[
    map_lgl(
      covariates,
      ~ n_distinct(model_data[[.x]], na.rm = TRUE) > 1L
    )
  ]

  model_formula <- reformulate(
    termlabels = c(".metabolite_exposure", estimable_covariates),
    response = ".clinical_outcome"
  )

  tryCatch({
    fit <- lm(stats::as.formula(model_formula), data = model_data)
    coefficient_table <- summary(fit)$coefficients
    coefficient_row <- match(
      ".metabolite_exposure",
      rownames(coefficient_table)
    )

    if (is.na(coefficient_row)) {
      return(empty_result("metabolite_coefficient_not_estimable"))
    }

    estimate <- coefficient_table[coefficient_row, "Estimate"]
    standard_error <- coefficient_table[coefficient_row, "Std. Error"]
    critical_value <- qt(0.975, df = fit$df.residual)

    tibble(
      clinical_trait = clinical_trait,
      metabolite = metabolite,
      clinical_period = clinical_period,
      metabolite_period = metabolite_period,
      n_obs = nrow(model_data),
      beta = unname(estimate),
      se = unname(standard_error),
      lower_ci = unname(estimate - critical_value * standard_error),
      upper_ci = unname(estimate + critical_value * standard_error),
      p_value = unname(coefficient_table[coefficient_row, "Pr(>|t|)"]),
      status = "ok"
    )
  }, error = function(error_condition) {
    empty_result(paste0("error: ", conditionMessage(error_condition)))
  })
}

run_association_grid <- function(data, clinical_traits, clinical_period,
                                 metabolite_period,
                                 metabolite_features = owo_metabolites) {
  tidyr::crossing(
    clinical_trait = clinical_traits,
    metabolite = metabolite_features
  ) %>%
    mutate(
      result = map2(
        clinical_trait,
        metabolite,
        ~ fit_one_association(
          data = data,
          clinical_trait = .x,
          metabolite = .y,
          clinical_period = clinical_period,
          metabolite_period = metabolite_period
        )
      )
    ) %>%
    select(result) %>%
    unnest(result) %>%
    group_by(clinical_period, metabolite_period, clinical_trait) %>%
    mutate(
      p_fdr = p.adjust(p_value, method = "fdr"),
      p_bonferroni = p.adjust(p_value, method = "bonferroni")
    ) %>%
    ungroup() %>%
    arrange(clinical_period, metabolite_period, clinical_trait, p_value)
}

# Analyses 
message("Running HBC T1 metabolites -> T1 clinical traits...")
results_t1_metab_t1_clinical <- run_association_grid(
  data = analysis_t1_metab_t1_clinical,
  clinical_traits = clini_var_t1_hbc,
  clinical_period = "T1",
  metabolite_period = "T1"
)

message("Running HBC T1 metabolites -> T2 OGTT traits...")
results_t1_metab_t2_clinical <- run_association_grid(
  data = analysis_t1_metab_t2_clinical,
  clinical_traits = clini_var_t2_hbc,
  clinical_period = "T2",
  metabolite_period = "T1"
)

message("Running HBC T1 metabolite PC1 score -> T1 clinical traits...")
results_t1_score_t1_clinical <- run_association_grid(
  data = analysis_t1_metab_t1_clinical,
  clinical_traits = clini_var_t1_hbc,
  clinical_period = "T1",
  metabolite_period = "T1",
  metabolite_features = "metabolite_score"
)

message("Running HBC T1 metabolite PC1 score -> T2 OGTT traits...")
results_t1_score_t2_clinical <- run_association_grid(
  data = analysis_t1_metab_t2_clinical,
  clinical_traits = clini_var_t2_hbc,
  clinical_period = "T2",
  metabolite_period = "T1",
  metabolite_features = "metabolite_score"
)

individual_metabolite_results <- bind_rows(
  results_t1_metab_t1_clinical,
  results_t1_metab_t2_clinical
) %>%
  left_join(metabolite_annotation, by = "metabolite") %>%
  mutate(feature_type = "Individual metabolite", .before = metabolite) %>%
  relocate(thsbc_metabolite, metabolite_name, .after = metabolite)

score_association_results <- bind_rows(
  results_t1_score_t1_clinical,
  results_t1_score_t2_clinical
) %>%
  mutate(
    feature_type = "PC1 score",
    thsbc_metabolite = NA_character_,
    metabolite_name = "PC1 of 27 overlapping OWO-related metabolites",
    .before = metabolite
  ) %>%
  relocate(thsbc_metabolite, metabolite_name, .after = metabolite)

association_results <- bind_rows(
  individual_metabolite_results,
  score_association_results
)

analysis_sample_summary <- bind_rows(
  tibble(
    clinical_period = "T1",
    metabolite_period = "T1",
    n_joined = n_distinct(analysis_t1_metab_t1_clinical$id)
  ),
  tibble(
    clinical_period = "T2",
    metabolite_period = "T1",
    n_joined = n_distinct(analysis_t1_metab_t2_clinical$id)
  )
)

feature_summary <- tibble(
  feature_set = "HBC-detectable THSBC OWO-related metabolites",
  n_features = length(owo_metabolites),
  source = "info_metab_hbc$Index_xh overlapping vip_sets$owo_vip",
  weighting = "None (ordinary linear regression)"
)

status_summary <- association_results %>%
  count(
    feature_type, clinical_period, metabolite_period, status,
    name = "n_models"
  )

# Output 
dir.create("result/PART2", recursive = TRUE, showWarnings = FALSE)
save(
  association_results,
  individual_metabolite_results,
  score_association_results,
  results_t1_metab_t1_clinical,
  results_t1_metab_t2_clinical,
  results_t1_score_t1_clinical,
  results_t1_score_t2_clinical,
  owo_metabolites,
  metabolite_annotation,
  pc1_summary,
  pc1_loadings,
  feature_summary,
  analysis_sample_summary,
  status_summary,
  file = "result/PART2/5.2_association_metab_clinical_traits_hbc.RData"
)

# Distance- and model-based variance estimation of the determinants ----------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(writexl)
  library(vegan)
  library(glmnet)
  library(lmerTest)
})

# 0. Load and harmonize data 
load("data/01_thsbc_cleaned_data/metabolites_1292_991.RData")
load("data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData")
load("data/01_thsbc_cleaned_data/gut_microbe_MetaG_1250_473.RData")
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")
load("data/01_thsbc_cleaned_data/wt_omics.rdata")

DF_metab_16SITS <- total_pop %>%
  filter(!is.na(PRS)) %>%
  inner_join(metab_tab_st, by = "id") %>%
  inner_join(bact_fung_clrz, by = c("id", "period")) %>%
  inner_join(cov_tab, by = c("id", "period")) %>%
  inner_join(wt_metab, by = c("id", "period"))

pc_file <- "tongji_5598_PC.eigenvec"
if (!file.exists(pc_file)) stop("Genetic PC file not found: ", pc_file)

pc_data <- fread(pc_file, header = FALSE) %>%
  transmute(
    id = str_pad(as.character(V1), width = 7, side = "left", pad = "0"),
    across(V3:V12)
  )
pc_data <- pc_data[, 1:4]
pc_vars <- paste0("PC", 1:3)
names(pc_data) <- c("id", pc_vars)

load("result/PART2/1_diff_metabolites_sig_sets_thsbc.rdata")
load("result/PART2/4.1_OPLS_DA_longitudinal/4.1_OPLS_DA_vip_sets.RData")
owo_vip_metabolites <- unique(vip_sets$owo_vip)

load("result/PART3/1_diff_16s_its_sig_sets_thsbc.rdata")
load("result/PART2/metab_gene/rec_vip40_1e5_clumped.Rdata")

demo_var <- c("age", "parity", "pbmi")
life_var <- c("smk", "drk", "physic_activity")
diet_var <- c( "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro")
social_var <- c("edu", "income", "medical_insurance")
micro_var <- unique(c(bact_uni, fung_uni))
gene_var <- unique(snp_df$var2)

period_levels <- c("T1", "T2", "T3")
screen_p_cutoff <- 0.05
nperm_screen <- 9999L
nperm_final <- 9999L
n_cores <- 12L

lasso_alpha <- 1
lasso_nfolds <- 5L
lasso_lambda_rule <- "lambda.min"
target_population_size <- 5212

# The same base covariates are used throughout the dbRDA and individual-
# metabolite analyses; microbiome/genetic models additionally control for
# antibiotic exposure/genetic ancestry PCs, respectively.
base_covars <- c("batch", "wk")
micro_covars <- c(base_covars, "antibiotic")
gene_covars <- c(base_covars, pc_vars)

output_dir <- "result/PART2/7_determinant_metabolites"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_data <- DF_metab_16SITS %>%
  inner_join(recode_snps, by = "id") %>%
  inner_join(pc_data, by = "id") %>%
  filter(period %in% period_levels) %>%
  mutate(
    period = factor(period, levels = period_levels),
    batch = factor(batch)
  ) %>%
  arrange(id, period)

required_vars <- unique(c(
  "id", "period", "prob", owo_vip_metabolites,
  demo_var, life_var, diet_var, social_var, micro_var, gene_var,
  base_covars, micro_covars, gene_covars
))
missing_vars <- setdiff(required_vars, names(analysis_data))
if (length(missing_vars) > 0L) {
  stop("Variables missing from analysis_data: ", paste(missing_vars, collapse = ", "))
}
if (anyDuplicated(analysis_data[c("id", "period")])) {
  stop("analysis_data contains duplicated id-period records")
}
if (length(owo_vip_metabolites) != 40L) {
  warning("Expected 40 OWO-related metabolites; found ", length(owo_vip_metabolites))
}

sample_summary <- bind_rows(
  analysis_data %>%
    count(period, name = "n_observations") %>%
    transmute(summary = as.character(period), n = n_observations),
  tibble(
    summary = c("total observations", "unique participants"),
    n = c(nrow(analysis_data), n_distinct(analysis_data$id))
  )
)
print(sample_summary)
readr::write_csv(sample_summary, file.path(output_dir, "7-00_sample_summary.csv"))

dbrda_required_vars <- unique(c(
  "id", "period", owo_vip_metabolites,
  demo_var, life_var, diet_var, social_var, micro_var, gene_var,
  base_covars, micro_covars, gene_covars
))
if (anyNA(analysis_data[dbrda_required_vars])) {
  stop("Missing values detected in variables required for the dbRDA analysis")
}

period_data <- setNames(
  lapply(period_levels, function(period_value) {
    analysis_data %>% filter(period == period_value) %>% droplevels()
  }),
  period_levels
)


# 1. Matrix-level dbRDA helper functions 
quote_var <- function(x) paste0("`", gsub("`", "", x), "`")

make_dbrda_formula <- function(predictors, conditions = character()) {
  rhs <- quote_var(predictors)
  if (length(conditions) > 0L) {
    rhs <- c(
      rhs,
      paste0(
        "Condition(", paste(quote_var(conditions), collapse = " + "), ")"
      )
    )
  }
  stats::as.formula(
    paste("Y ~", paste(rhs, collapse = " + ")),
    env = parent.frame()
  )
}

extract_dbrda_result <- function(fit, test, period, group, predictors, n) {
  r2_stats <- vegan::RsquareAdj(fit)
  tibble(
    period = period,
    group = group,
    n = n,
    n_predictors = length(predictors),
    predictors = paste(predictors, collapse = ";"),
    df = unname(test$Df[1]),
    r2 = unname(r2_stats[["r.squared"]]),
    adjusted_r2 = unname(r2_stats[["adj.r.squared"]]),
    pseudo_F = unname(test$F[1]),
    pvalue = unname(test$`Pr(>F)`[1])
  )
}

screen_one_matrix_predictor <- function(
    data, predictor, conditions, group, predictor_index,
    nperm = nperm_screen, seed = 1000L) {
  period_value <- as.character(unique(data$period))
  tryCatch({
    Y <- as.matrix(data[owo_vip_metabolites])
    fit <- vegan::dbrda(
      make_dbrda_formula(predictor, conditions),
      data = data,
      distance = "euclidean"
    )
    set.seed(seed + predictor_index)
    test <- anova(
      fit,
      permutations = nperm,
      parallel = n_cores,
      by = NULL
    )
    r2_stats <- vegan::RsquareAdj(fit)

    tibble(
      period = period_value,
      group = group,
      predictor = predictor,
      n = nrow(data),
      df = unname(test$Df[1]),
      r2 = unname(r2_stats[["r.squared"]]),
      adjusted_r2 = unname(r2_stats[["adj.r.squared"]]),
      pseudo_F = unname(test$F[1]),
      pvalue = unname(test$`Pr(>F)`[1]),
      status = "ok",
      error = NA_character_
    )
  }, error = function(e) {
    tibble(
      period = period_value,
      group = group,
      predictor = predictor,
      n = nrow(data),
      df = NA_real_,
      r2 = NA_real_,
      adjusted_r2 = NA_real_,
      pseudo_F = NA_real_,
      pvalue = NA_real_,
      status = "error",
      error = conditionMessage(e)
    )
  })
}

screen_matrix_predictors <- function(data, predictors, conditions, group, period_seed) {
  purrr::map2_dfr(
    predictors,
    seq_along(predictors),
    ~ screen_one_matrix_predictor(
      data = data,
      predictor = .x,
      conditions = conditions,
      group = group,
      predictor_index = .y,
      seed = period_seed
    )
  ) %>%
    mutate(
      qvalue_bh = p.adjust(pvalue, method = "BH"),
      selected_nominal_p05 = status == "ok" & pvalue < screen_p_cutoff
    ) %>%
    arrange(pvalue)
}

fit_matrix_group <- function(data, predictors, conditions, group,
                             nperm = nperm_final, seed = 10L) {
  period_value <- as.character(unique(data$period))
  if (length(predictors) == 0L) {
    return(list(
      result = tibble(
        period = period_value,
        group = group,
        n = nrow(data),
        n_predictors = 0L,
        predictors = NA_character_,
        df = NA_real_,
        r2 = NA_real_,
        adjusted_r2 = NA_real_,
        pseudo_F = NA_real_,
        pvalue = NA_real_
      ),
      model = NULL
    ))
  }

  Y <- as.matrix(data[owo_vip_metabolites])
  fit <- vegan::dbrda(
    make_dbrda_formula(predictors, conditions),
    data = data,
    distance = "euclidean"
  )
  set.seed(seed)
  test <- anova(
    fit,
    permutations = nperm,
    parallel = n_cores,
    by = NULL
  )

  list(
    result = extract_dbrda_result(
      fit, test, period_value, group, predictors, nrow(data)
    ),
    model = fit
  )
}


# 2. Period-specific microbiota/SNP dbRDA screening 
# Screen microbial taxa and SNPs independently within each timepoint
micro_screen_all_periods <- purrr::map_dfr(
  period_levels,
  function(period_value) {
    screen_matrix_predictors(
      data = period_data[[period_value]],
      predictors = micro_var,
      conditions = micro_covars,
      group = "microbiota",
      period_seed = 1000L + match(period_value, period_levels) * 10000L
    )
  }
)

gene_screen_all_periods <- purrr::map_dfr(
  period_levels,
  function(period_value) {
    screen_matrix_predictors(
      data = period_data[[period_value]],
      predictors = gene_var,
      conditions = gene_covars,
      group = "genetics",
      period_seed = 2000L + match(period_value, period_levels) * 10000L
    )
  }
)

micro_selected_by_period <- setNames(
  lapply(period_levels, function(period_value) {
    micro_screen_all_periods %>%
      filter(period == period_value, selected_nominal_p05) %>%
      pull(predictor)
  }),
  period_levels
)

gene_selected_by_period <- setNames(
  lapply(period_levels, function(period_value) {
    gene_screen_all_periods %>%
      filter(period == period_value, selected_nominal_p05) %>%
      pull(predictor)
  }),
  period_levels
)

matrix_screen_summary <- bind_rows(
  purrr::map_dfr(period_levels, function(period_value) {
    tibble(
      period = period_value,
      group = "microbiota",
      n_candidates = length(micro_var),
      n_selected_nominal_p05 = length(micro_selected_by_period[[period_value]]),
      selected_predictors = paste(
        micro_selected_by_period[[period_value]], collapse = ";"
      )
    )
  }),
  purrr::map_dfr(period_levels, function(period_value) {
    tibble(
      period = period_value,
      group = "genetics",
      n_candidates = length(gene_var),
      n_selected_nominal_p05 = length(gene_selected_by_period[[period_value]]),
      selected_predictors = paste(
        gene_selected_by_period[[period_value]], collapse = ";"
      )
    )
  })
) %>%
  arrange(period, group)

matrix_predictor_period_patterns <- bind_rows(
  micro_screen_all_periods,
  gene_screen_all_periods
) %>%
  group_by(group, predictor) %>%
  summarise(
    n_periods_tested = n_distinct(period),
    n_periods_selected = sum(selected_nominal_p05, na.rm = TRUE),
    selected_periods = paste(period[selected_nominal_p05 %in% TRUE], collapse = ";"),
    selection_pattern = if_else(
      n_periods_selected == 0L, "not selected", selected_periods
    ),
    period_specific = n_periods_selected == 1L,
    period_specific_period = if_else(
      period_specific, selected_periods, NA_character_
    ),
    .groups = "drop"
  ) %>%
  arrange(group, desc(period_specific), period_specific_period, predictor)

readr::write_csv(
  micro_screen_all_periods,
  file.path(output_dir, "7-01_microbiota_single_predictor_dbrda.csv")
)
readr::write_csv(
  gene_screen_all_periods,
  file.path(output_dir, "7-02_genetics_single_predictor_dbrda.csv")
)
readr::write_csv(
  matrix_screen_summary,
  file.path(output_dir, "7-03_matrix_selected_summary.csv")
)
readr::write_csv(
  matrix_predictor_period_patterns,
  file.path(output_dir, "7-04_matrix_predictor_period_patterns.csv")
)

save(
  micro_screen_all_periods,
  gene_screen_all_periods,
  micro_selected_by_period,
  gene_selected_by_period,
  matrix_screen_summary,
  matrix_predictor_period_patterns,
  file = file.path(output_dir, "7-05_matrix_screening_and_selected_sets.RData")
)


# 3. Final period-specific dbRDA models 
matrix_results_list <- list()
matrix_models <- list()
matrix_group_specs_by_period <- list()

for (period_value in period_levels) {
  matrix_group_specs <- list(
    demographics = list(vars = demo_var, covars = base_covars),
    lifestyle = list(vars = life_var, covars = base_covars),
    diet = list(vars = diet_var, covars = base_covars),
    social = list(vars = social_var, covars = base_covars),
    microbiota = list(
      vars = micro_selected_by_period[[period_value]],
      covars = micro_covars
    ),
    genetics = list(
      vars = gene_selected_by_period[[period_value]],
      covars = gene_covars
    )
  )
  matrix_group_specs_by_period[[period_value]] <- matrix_group_specs

  for (group_name in names(matrix_group_specs)) {
    spec <- matrix_group_specs[[group_name]]
    fit_object <- fit_matrix_group(
      data = period_data[[period_value]],
      predictors = spec$vars,
      conditions = spec$covars,
      group = group_name,
      nperm = nperm_final,
      seed = 3000L + match(period_value, period_levels) * 100L +
        match(group_name, names(matrix_group_specs))
    )
    result_name <- paste(period_value, group_name, sep = "_")
    matrix_results_list[[result_name]] <- fit_object$result
    matrix_models[[result_name]] <- fit_object$model
  }
}

matrix_dbrda_results <- bind_rows(matrix_results_list) %>%
  group_by(group) %>%
  mutate(qvalue_bh_across_periods = p.adjust(pvalue, method = "BH")) %>%
  ungroup()

print(matrix_dbrda_results)
readr::write_csv(
  matrix_dbrda_results,
  file.path(output_dir, "7-06_matrix_dbrda_period_specific.csv")
)


# 4. Weighted longitudinal mixed-model helper functions
make_mixed_formula <- function(conditions, predictors = character()) {
  fixed_terms <- c(conditions, predictors)
  stats::as.formula(
    paste0(
      ".outcome ~ ", paste(quote_var(fixed_terms), collapse = " + "),
      " + (1 | id)"
    ),
    env = parent.frame()
  )
}

prepare_model_data <- function(data, outcome, predictors, conditions) {
  required <- unique(c("id", "period", "prob", outcome, predictors, conditions))
  model_data <- data %>%
    select(all_of(required)) %>%
    filter(if_all(everything(), ~ !is.na(.x))) %>%
    filter(is.finite(prob), prob > 0, prob <= 1) %>%
    mutate(
      period = factor(period, levels = period_levels),
      .outcome = .data[[outcome]]
    ) %>%
    group_by(period) %>%
    mutate(
      .analysis_weight = (1 / prob) *
        target_population_size / sum(1 / prob)
    ) %>%
    ungroup() %>%
    droplevels()

  if (nrow(model_data) < 30L || n_distinct(model_data$id) < 20L) {
    stop("Insufficient complete observations")
  }
  if (length(unique(model_data$.outcome)) <= 1L) stop("Outcome is constant")
  model_data
}

weighted_variance <- function(x, weights) {
  weighted_mean <- stats::weighted.mean(x, weights)
  sum(weights * (x - weighted_mean)^2) / sum(weights)
}

mixed_r2 <- function(model, weights) {
  fixed_prediction <- stats::predict(model, re.form = NA)
  fixed_variance <- weighted_variance(fixed_prediction, weights)
  variance_components <- as.data.frame(lme4::VarCorr(model))
  random_variance <- sum(
    variance_components$vcov[variance_components$grp != "Residual"]
  )
  residual_variance <- stats::weighted.mean(
    stats::sigma(model)^2 / weights,
    weights
  )
  total_variance <- fixed_variance + random_variance + residual_variance

  c(
    marginal = fixed_variance / total_variance,
    conditional = (fixed_variance + random_variance) / total_variance
  )
}

adjust_marginal_r2 <- function(r2, n, p) {
  if (!is.finite(r2) || n <= p + 1L) return(NA_real_)
  1 - (1 - r2) * (n - 1) / (n - p - 1)
}

make_id_foldid <- function(ids, nfolds, seed) {
  unique_ids <- unique(as.character(ids))
  set.seed(seed)
  id_folds <- sample(rep(seq_len(nfolds), length.out = length(unique_ids)))
  id_folds[match(as.character(ids), unique_ids)]
}

empty_mixed_result <- function(outcome, group, n_candidates, status, error) {
  tibble(
    period = "T1-T3",
    metabolite = outcome,
    group = group,
    n_observations = NA_integer_,
    n_participants = NA_integer_,
    n_candidate_predictors = n_candidates,
    n_selected_predictors = NA_integer_,
    selected_predictors = NA_character_,
    base_marginal_r2 = NA_real_,
    full_marginal_r2 = NA_real_,
    delta_marginal_r2 = NA_real_,
    base_adjusted_marginal_r2 = NA_real_,
    full_adjusted_marginal_r2 = NA_real_,
    delta_adjusted_marginal_r2 = NA_real_,
    delta_adjusted_marginal_r2_nonnegative = NA_real_,
    base_conditional_r2 = NA_real_,
    full_conditional_r2 = NA_real_,
    likelihood_ratio = NA_real_,
    group_pvalue = NA_real_,
    singular_base = NA,
    singular_full = NA,
    status = status,
    error = error
  )
}

fit_weighted_mixed_group <- function(data, outcome, predictors, conditions,
                                     group,
                                     n_candidates = length(predictors)) {
  tryCatch({
    model_data <- prepare_model_data(data, outcome, predictors, conditions)

    base_model <- lmerTest::lmer(
      make_mixed_formula(conditions),
      data = model_data,
      weights = .analysis_weight,
      REML = FALSE
    )
    base_r2 <- mixed_r2(base_model, model_data$.analysis_weight)

    if (length(predictors) > 0L) {
      full_model <- lmerTest::lmer(
        make_mixed_formula(conditions, predictors),
        data = model_data,
        weights = .analysis_weight,
        REML = FALSE
      )
      full_r2 <- mixed_r2(full_model, model_data$.analysis_weight)
      comparison <- stats::anova(base_model, full_model)
      likelihood_ratio <- unname(comparison$Chisq[2])
      group_pvalue <- unname(comparison$`Pr(>Chisq)`[2])
      singular_full <- lme4::isSingular(full_model, tol = 1e-4)
      p_full <- length(lme4::fixef(full_model)) - 1L
    } else {
      full_model <- base_model
      full_r2 <- base_r2
      likelihood_ratio <- NA_real_
      group_pvalue <- NA_real_
      singular_full <- lme4::isSingular(base_model, tol = 1e-4)
      p_full <- length(lme4::fixef(base_model)) - 1L
    }

    p_base <- length(lme4::fixef(base_model)) - 1L
    base_adjusted <- adjust_marginal_r2(
      unname(base_r2["marginal"]), nrow(model_data), p_base
    )
    full_adjusted <- adjust_marginal_r2(
      unname(full_r2["marginal"]), nrow(model_data), p_full
    )
    delta_adjusted <- full_adjusted - base_adjusted

    tibble(
      period = "T1-T3",
      metabolite = outcome,
      group = group,
      n_observations = nrow(model_data),
      n_participants = n_distinct(model_data$id),
      n_candidate_predictors = n_candidates,
      n_selected_predictors = length(predictors),
      selected_predictors = if (
        length(predictors) == 0L
      ) NA_character_ else paste(predictors, collapse = ";"),
      base_marginal_r2 = unname(base_r2["marginal"]),
      full_marginal_r2 = unname(full_r2["marginal"]),
      delta_marginal_r2 = unname(
        full_r2["marginal"] - base_r2["marginal"]
      ),
      base_adjusted_marginal_r2 = base_adjusted,
      full_adjusted_marginal_r2 = full_adjusted,
      delta_adjusted_marginal_r2 = delta_adjusted,
      delta_adjusted_marginal_r2_nonnegative = pmax(delta_adjusted, 0),
      base_conditional_r2 = unname(base_r2["conditional"]),
      full_conditional_r2 = unname(full_r2["conditional"]),
      likelihood_ratio = likelihood_ratio,
      group_pvalue = group_pvalue,
      singular_base = lme4::isSingular(base_model, tol = 1e-4),
      singular_full = singular_full,
      status = "ok",
      error = NA_character_
    )
  }, error = function(e) {
    empty_mixed_result(
      outcome, group, n_candidates, "model_error", conditionMessage(e)
    )
  })
}


# 5. Low-dimensional groups: direct weighted mixed models 
mixed_lowdim_specs <- list(
  demographics = list(vars = demo_var, covars = base_covars),
  lifestyle = list(vars = life_var, covars = base_covars),
  diet = list(vars = diet_var, covars = base_covars),
  social = list(vars = social_var, covars = base_covars)
)

mixed_lowdim_results <- purrr::imap_dfr(
  mixed_lowdim_specs,
  function(spec, group_name) {
    purrr::map_dfr(owo_vip_metabolites, function(outcome) {
      fit_weighted_mixed_group(
        data = analysis_data,
        outcome = outcome,
        predictors = spec$vars,
        conditions = spec$covars,
        group = group_name
      )
    })
  }
) %>%
  group_by(group) %>%
  mutate(group_qvalue_bh_40_metabolites = p.adjust(group_pvalue, method = "BH")) %>%
  ungroup()

readr::write_csv(
  mixed_lowdim_results,
  file.path(output_dir, "7-07_mixed_lowdim_results.csv")
)


# 6. ID-folded LASSO for microbiota and genetics 
mixed_lasso_specs <- list(
  microbiota = list(vars = micro_var, covars = micro_covars),
  genetics = list(vars = gene_var, covars = gene_covars)
)

select_lasso_predictors <- function(data, outcome, candidate_predictors,
                                    conditions, group, fold_seed) {
  tryCatch({
    selection_data <- prepare_model_data(
      data, outcome, candidate_predictors, conditions
    )
    variable_predictors <- candidate_predictors[
      vapply(
        selection_data[candidate_predictors],
        function(x) length(unique(x)) > 1L,
        logical(1)
      )
    ]
    if (length(variable_predictors) < 1L) {
      stop("No non-constant candidate predictors")
    }

    x_condition <- model.matrix(
      stats::reformulate(conditions),
      data = selection_data
    )[, -1, drop = FALSE]
    x_candidate <- data.matrix(selection_data[variable_predictors])
    x <- cbind(x_condition, x_candidate)
    penalty_factor <- c(
      rep(0, ncol(x_condition)),
      rep(1, ncol(x_candidate))
    )
    foldid <- make_id_foldid(selection_data$id, lasso_nfolds, fold_seed)

    set.seed(fold_seed)
    cvfit <- glmnet::cv.glmnet(
      x = x,
      y = selection_data$.outcome,
      family = "gaussian",
      alpha = lasso_alpha,
      nfolds = lasso_nfolds,
      foldid = foldid,
      penalty.factor = penalty_factor,
      standardize = TRUE,
      type.measure = "mse"
    )

    coefficient_vector <- as.matrix(
      coef(cvfit, s = lasso_lambda_rule)
    )[, 1]
    selected_features <- intersect(
      names(coefficient_vector)[coefficient_vector != 0],
      variable_predictors
    )
    selected_coefficients <- coefficient_vector[selected_features]

    tibble(
      metabolite = outcome,
      group = group,
      n_observations = nrow(selection_data),
      n_participants = n_distinct(selection_data$id),
      n_candidate_predictors = length(variable_predictors),
      n_selected_predictors = length(selected_features),
      selected_predictors = if (
        length(selected_features) == 0L
      ) NA_character_ else paste(selected_features, collapse = ";"),
      selected_coefficients = if (
        length(selected_coefficients) == 0L
      ) NA_character_ else paste(
        paste0(
          names(selected_coefficients), "=",
          signif(as.numeric(selected_coefficients), 6)
        ),
        collapse = ";"
      ),
      alpha = lasso_alpha,
      nfolds = lasso_nfolds,
      lambda_rule = lasso_lambda_rule,
      lambda_selected = cvfit$lambda.min,
      lambda_min = cvfit$lambda.min,
      lambda_1se = cvfit$lambda.1se,
      cvm_min = min(cvfit$cvm),
      selection_status = "ok",
      selection_error = NA_character_
    )
  }, error = function(e) {
    tibble(
      metabolite = outcome,
      group = group,
      n_observations = NA_integer_,
      n_participants = NA_integer_,
      n_candidate_predictors = length(candidate_predictors),
      n_selected_predictors = NA_integer_,
      selected_predictors = NA_character_,
      selected_coefficients = NA_character_,
      alpha = lasso_alpha,
      nfolds = lasso_nfolds,
      lambda_rule = lasso_lambda_rule,
      lambda_selected = NA_real_,
      lambda_min = NA_real_,
      lambda_1se = NA_real_,
      cvm_min = NA_real_,
      selection_status = "selection_error",
      selection_error = conditionMessage(e)
    )
  })
}

mixed_lasso_selection <- purrr::imap_dfr(
  mixed_lasso_specs,
  function(spec, group_name) {
    group_seed <- ifelse(group_name == "microbiota", 100L, 200L)
    purrr::map2_dfr(
      owo_vip_metabolites,
      seq_along(owo_vip_metabolites),
      function(outcome, outcome_index) {
        select_lasso_predictors(
          data = analysis_data,
          outcome = outcome,
          candidate_predictors = spec$vars,
          conditions = spec$covars,
          group = group_name,
          fold_seed = 9000L + group_seed + outcome_index
        )
      }
    )
  }
)

if (anyDuplicated(mixed_lasso_selection[c("metabolite", "group")])) {
  stop("Duplicated metabolite-group rows in mixed_lasso_selection")
}

readr::write_csv(
  mixed_lasso_selection,
  file.path(output_dir, "7-08_lasso_selection_results.csv")
)


# 7. Refit selected high-dimensional groups with weighted mixed models 
split_selected_predictors <- function(x) {
  if (length(x) == 0L || is.na(x) || trimws(x) == "") return(character())
  strsplit(x, ";", fixed = TRUE)[[1]]
}

mixed_lasso_results <- purrr::pmap_dfr(
  mixed_lasso_selection,
  function(metabolite, group, n_observations, n_participants,
           n_candidate_predictors, n_selected_predictors,
           selected_predictors, selected_coefficients, alpha, nfolds,
           lambda_rule, lambda_selected, lambda_min, lambda_1se, cvm_min,
           selection_status, selection_error) {
    if (selection_status != "ok") {
      return(
        empty_mixed_result(
          metabolite, group, n_candidate_predictors,
          "selection_error", selection_error
        ) %>%
          mutate(
            selected_coefficients = selected_coefficients,
            alpha = alpha,
            nfolds = nfolds,
            lambda_rule = lambda_rule,
            lambda_selected = lambda_selected,
            lambda_min = lambda_min,
            lambda_1se = lambda_1se,
            cvm_min = cvm_min
          )
      )
    }

    spec <- mixed_lasso_specs[[group]]
    selected_features <- split_selected_predictors(selected_predictors)
    unknown_selected <- setdiff(selected_features, spec$vars)
    if (length(unknown_selected) > 0L) {
      stop(
        "Selected predictors absent from candidate set for ", metabolite,
        "/", group, ": ", paste(unknown_selected, collapse = ", ")
      )
    }

    # Recreate the same complete-case analysis sample used at selection.
    selection_data <- prepare_model_data(
      analysis_data, metabolite, spec$vars, spec$covars
    )
    fit_result <- fit_weighted_mixed_group(
      data = selection_data,
      outcome = metabolite,
      predictors = selected_features,
      conditions = spec$covars,
      group = group,
      n_candidates = n_candidate_predictors
    )

    fit_result %>%
      mutate(
        n_observations = nrow(selection_data),
        n_participants = n_distinct(selection_data$id),
        n_selected_predictors = length(selected_features),
        selected_predictors = if (
          length(selected_features) == 0L
        ) NA_character_ else paste(selected_features, collapse = ";"),
        selected_coefficients = selected_coefficients,
        alpha = alpha,
        nfolds = nfolds,
        lambda_rule = lambda_rule,
        lambda_selected = lambda_selected,
        lambda_min = lambda_min,
        lambda_1se = lambda_1se,
        cvm_min = cvm_min
      )
  }
)

readr::write_csv(
  mixed_lasso_results,
  file.path(output_dir, "7-09_mixed_lasso_results.csv")
)

lasso_selection_summary <- mixed_lasso_results %>%
  group_by(group) %>%
  summarise(
    n_metabolites = n(),
    n_success = sum(status == "ok"),
    n_zero_selected = sum(status == "ok" & n_selected_predictors == 0L),
    median_selected = median(n_selected_predictors[status == "ok"], na.rm = TRUE),
    min_selected = min(n_selected_predictors[status == "ok"], na.rm = TRUE),
    max_selected = max(n_selected_predictors[status == "ok"], na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  lasso_selection_summary,
  file.path(output_dir, "7-10_lasso_selection_summary.csv")
)


# 8. Combined individual-metabolite variance output 
individual_metabolite_variance <- bind_rows(
  mixed_lowdim_results %>%
    transmute(
      period, metabolite, group,
      n_observations, n_participants,
      n_candidate_predictors, n_selected_predictors,
      selected_predictors,
      base_marginal_r2, full_marginal_r2, delta_marginal_r2,
      base_adjusted_marginal_r2, full_adjusted_marginal_r2,
      delta_adjusted_marginal_r2,
      delta_adjusted_marginal_r2_nonnegative,
      base_conditional_r2, full_conditional_r2,
      group_pvalue, group_qvalue_bh_40_metabolites,
      alpha = NA_real_,
      method = "weighted linear mixed-effects model",
      singular_base, singular_full, status, error
    ),
  mixed_lasso_results %>%
    transmute(
      period, metabolite, group,
      n_observations, n_participants,
      n_candidate_predictors, n_selected_predictors,
      selected_predictors,
      base_marginal_r2, full_marginal_r2, delta_marginal_r2,
      base_adjusted_marginal_r2, full_adjusted_marginal_r2,
      delta_adjusted_marginal_r2,
      delta_adjusted_marginal_r2_nonnegative,
      base_conditional_r2, full_conditional_r2,
      group_pvalue,
      group_qvalue_bh_40_metabolites = NA_real_,
      alpha,
      method = paste0(
        "ID-folded LASSO ", lasso_lambda_rule,
        " + weighted mixed model"
      ),
      singular_base, singular_full, status, error
    )
) %>%
  arrange(metabolite, group)

readr::write_csv(
  individual_metabolite_variance,
  file.path(output_dir, "7-11_individual_metabolite_marginal_r2.csv")
)

writexl::write_xlsx(
  list(
    sample_summary = sample_summary,
    matrix_screen_summary = matrix_screen_summary,
    matrix_dbrda = matrix_dbrda_results,
    microbiota_dbrda_screen = micro_screen_all_periods,
    genetics_dbrda_screen = gene_screen_all_periods,
    mixed_lowdim = mixed_lowdim_results,
    lasso_selection = mixed_lasso_selection,
    mixed_lasso = mixed_lasso_results,
    lasso_selection_summary = lasso_selection_summary,
    individual_marginal_r2 = individual_metabolite_variance
  ),
  path = file.path(output_dir, "7_determinant_metabolites_results.xlsx")
)

save(
  analysis_data,
  sample_summary,
  period_data,
  matrix_group_specs_by_period,
  micro_screen_all_periods,
  gene_screen_all_periods,
  micro_selected_by_period,
  gene_selected_by_period,
  matrix_screen_summary,
  matrix_predictor_period_patterns,
  matrix_dbrda_results,
  matrix_models,
  mixed_lowdim_specs,
  mixed_lasso_specs,
  mixed_lowdim_results,
  mixed_lasso_selection,
  mixed_lasso_results,
  lasso_selection_summary,
  individual_metabolite_variance,
  nperm_screen,
  nperm_final,
  screen_p_cutoff,
  lasso_alpha,
  lasso_nfolds,
  lasso_lambda_rule,
  target_population_size,
  file = file.path(output_dir, "7_determinant_metabolites_results.RData")
)

message(
  "Determinant analysis completed. Results saved to: ", output_dir
)

