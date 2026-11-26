##############################################################################################
# Analyses for Figure 5
# Associations between prepregnancy BMI (pBMI) and study outcomes
# DIABLO score- and member-level mediation for Figures 5c,d
# See plot_fig5.R for the corresponding figure-generation code.
###############################################################################################

# Associations between pBMI and maternal and offspring outcomes ------------------------------------

rm(list = ls())

library(tidyverse)
library(purrr)
library(readxl)
library(writexl)

## 0. Configuration and helper functions 
binary_otc <- c("GDM", "HDP", "PE", "ICP", "LBW", "MAC", "PTB", "SGA", "LGA", "owo_m36")
contin_otc <- c("Delivery_week", "Birth_weight", "Birth_length", "zbmi_m36")
pregnancy_complications <- c("GDM", "HDP", "PE", "ICP")
early_gwg_outcomes <- c(pregnancy_complications, "PTB", "Delivery_week")

Covars <- c(
  "age", "parity", "smk", "drk", "physic_activity",
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
  "edu", "medical_insurance", "income"
)

pbmi_exposures <- c(
  "pbmi_z",
  "pbmi_owo_normal",
  "pbmi_overweight_normal",
  "pbmi_obesity_normal",
  "pbmi_underweight_normal",
  "pbmi_obesity_overweight"
)

z_trans <- function(x) {
  if (all(is.na(x)) || is.na(sd(x, na.rm = TRUE)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

model_covars <- function(outcome, include_sex = TRUE, total_gwg_in_model = FALSE) {
  out <- Covars
  if (include_sex && !outcome %in% pregnancy_complications) out <- c(out, "baby_sex")
  if (total_gwg_in_model) out <- c(out, "time_last")
  unique(out)
}

empty_result <- function(analysis, outcome, exposure, model_type, status, n = NA_integer_) {
  tibble(
    analysis = analysis, outcome = outcome, exposure = exposure,
    model_type = model_type, term = NA_character_, n = as.integer(n),
    beta = NA_real_, se = NA_real_, lci = NA_real_, uci = NA_real_,
    effect = NA_real_, effect_lci = NA_real_, effect_uci = NA_real_,
    p = NA_real_, status = status
  )
}

# Complete-case glm/lm with a consistent output. For logistic models, beta is
# the log odds ratio and effect is the odds ratio; for linear models they match.
fit_terms <- function(data, outcome, exposure, covars, model_type, analysis, term_selector = c("exposure", "interaction")) {
  term_selector <- match.arg(term_selector)
  required <- unique(c(outcome, exposure, covars))
  missing_vars <- setdiff(required, names(data))
  if (length(missing_vars) > 0) {
    return(empty_result(
      analysis, outcome, exposure, model_type,
      paste0("missing_variables: ", paste(missing_vars, collapse = ", "))
    ))
  }

  d <- data %>% filter(if_all(all_of(required), ~ !is.na(.x)))
  exposure_has_variation <- all(map_lgl(exposure, ~ n_distinct(d[[.x]]) >= 2))
  if (nrow(d) < 20 || n_distinct(d[[outcome]]) < 2 || !exposure_has_variation) {
    return(empty_result(analysis, outcome, exposure, model_type,
                        "insufficient_data", nrow(d)))
  }

  rhs <- if (term_selector == "interaction") {
    paste0(exposure[1], " * ", exposure[2],
           if (length(covars) > 0) paste0(" + ", paste(covars, collapse = " + ")) else "")
  } else {
    paste(c(exposure, covars), collapse = " + ")
  }

  tryCatch({
    formula <- as.formula(paste(outcome, "~", rhs))
    fit <- if (model_type == "logistic") {
      glm(formula, family = binomial(link = "logit"), data = d)
    } else {
      lm(formula, data = d)
    }
    coef_tab <- summary(fit)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("term")
    p_col <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), names(coef_tab))[1]
    keep <- if (term_selector == "interaction") {
      str_detect(coef_tab$term, fixed(":"))
    } else {
      str_starts(coef_tab$term, fixed(exposure[1]))
    }
    coef_tab <- coef_tab[keep, , drop = FALSE]
    if (nrow(coef_tab) == 0) {
      return(empty_result(analysis, outcome, paste(exposure, collapse = " * "),
                          model_type, "term_not_estimable", nobs(fit)))
    }
    coef_tab %>%
      transmute(
        analysis = analysis,
        outcome = outcome,
        exposure = paste(exposure, collapse = " * "),
        model_type = model_type,
        term = term,
        n = as.integer(nobs(fit)),
        beta = Estimate,
        se = `Std. Error`,
        lci = Estimate - qnorm(0.975) * `Std. Error`,
        uci = Estimate + qnorm(0.975) * `Std. Error`,
        effect = if_else(model_type == "logistic", exp(Estimate), Estimate),
        effect_lci = if_else(model_type == "logistic", exp(lci), lci),
        effect_uci = if_else(model_type == "logistic", exp(uci), uci),
        p = .data[[p_col]],
        status = "ok"
      )
  }, error = function(e) {
    empty_result(analysis, outcome, paste(exposure, collapse = " * "),
                 model_type, paste0("error: ", conditionMessage(e)), nrow(d))
  })
}

outcome_model_type <- function(outcome) {
  if_else(outcome %in% binary_otc, "logistic", "linear")
}

run_pbmi_outcomes <- function(data, analysis, include_sex = TRUE) {
  crossing(outcome = c(binary_otc, contin_otc), exposure = pbmi_exposures) %>%
    pmap_dfr(function(outcome, exposure) {
      fit_terms(data, outcome, exposure,
                model_covars(outcome, include_sex = include_sex),
                outcome_model_type(outcome), analysis)
    })
}

## 1. Load and prepare participant-level data
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")
load("data/01_thsbc_cleaned_data/outcomes_5212.RData")
load("expo_cov_otc_dat_6943.Rdata")

weight_gain_dat <- expo_cov_otc_dat %>%
  select(id, rate_fst) %>%
  distinct(id, .keep_all = TRUE)

DF_totalpop <- total_pop %>%
  inner_join(cov_tab %>% filter(period == "T1") %>% distinct(id, .keep_all = TRUE), by = "id") %>%
  inner_join(otc, by = "id") %>%
  left_join(weight_gain_dat, by = "id")

# Each binary pBMI contrast excludes all women outside its two requested groups.
DF_totalpop <- DF_totalpop %>%
  mutate(
    pbmi_z = z_trans(pbmi),
    Delivery_week = z_trans(Delivery_week),
    Birth_weight = z_trans(Birth_weight),
    Birth_length = z_trans(Birth_length),
    pbmi_owo_normal = factor(
      if_else(pbmi_cut3 %in% c("Normal", "OWO"), as.character(pbmi_cut3), NA_character_),
      levels = c("Normal", "OWO")
    ),
    pbmi_overweight_normal = factor(
      if_else(pbmi_cut4 %in% c("Normal", "Overweight"), as.character(pbmi_cut4), NA_character_),
      levels = c("Normal", "Overweight")
    ),
    pbmi_obesity_normal = factor(
      if_else(pbmi_cut4 %in% c("Normal", "Obese"), as.character(pbmi_cut4), NA_character_),
      levels = c("Normal", "Obese")
    ),
    pbmi_underweight_normal = factor(
      if_else(pbmi_cut4 %in% c("Normal", "Underweight"), as.character(pbmi_cut4), NA_character_),
      levels = c("Normal", "Underweight")
    ),
    pbmi_obesity_overweight = factor(
      if_else(pbmi_cut4 %in% c("Overweight", "Obese"), as.character(pbmi_cut4), NA_character_),
      levels = c("Overweight", "Obese")
    ),
    WG_cut = factor(WG_cut, levels = c("Adequate", "Low", "Excessive")),
    baby_sex = as.factor(baby_sex),
    across(all_of(Covars), ~ if (is.character(.x)) as.factor(.x) else .x)
  )

## Analysis 1: pBMI and study outcomes 

res_a1_pbmi_outcomes <- run_pbmi_outcomes(
  DF_totalpop, "A1_pBMI_outcomes", include_sex = TRUE
)

## Analysis 1 sensitivity: additional adjustment for weight gain 

# To isolate the impact of adding the weight-gain covariate, the base and
# additionally adjusted models below use exactly the same complete-case sample.
# GDM/HDP/PE/ICP, PTB, and delivery week are adjusted for early-pregnancy
# rate_fst. Other outcomes are adjusted for total weight_gain and time_last,
# the gestational week of the last prenatal weight measurement.
a1_sensitivity_pairs <- crossing(
  outcome = c(binary_otc, contin_otc), exposure = pbmi_exposures
) %>%
  pmap(function(outcome, exposure) {
    gwg_covar <- if (outcome %in% early_gwg_outcomes) "rate_fst" else "weight_gain"
    base_covars <- model_covars(outcome, include_sex = TRUE)
    added_covars <- unique(c(
      base_covars, gwg_covar,
      if (identical(gwg_covar, "weight_gain")) "time_last" else character(0)
    ))
    same_sample <- DF_totalpop %>%
      filter(if_all(all_of(unique(c(outcome, exposure, added_covars))), ~ !is.na(.x)))

    list(
      base = fit_terms(
        same_sample, outcome, exposure, base_covars,
        outcome_model_type(outcome), "A1S_base_same_sample"
      ) %>% mutate(added_weight_gain_covar = gwg_covar, .after = exposure),
      adjusted = fit_terms(
        same_sample, outcome, exposure, added_covars,
        outcome_model_type(outcome), "A1S_weight_gain_adjusted"
      ) %>% mutate(added_weight_gain_covar = gwg_covar, .after = exposure)
    )
  })

res_a1s_base_same_sample <- map_dfr(a1_sensitivity_pairs, "base")
res_a1s_weight_gain_adjusted <- map_dfr(a1_sensitivity_pairs, "adjusted")

## Analysis 2: pBMI and gestational weight-gain indicators 
res_a2_pbmi_weight_gain <- crossing(
  outcome = c("weight_gain", "rate_fst"), exposure = pbmi_exposures
) %>%
  pmap_dfr(function(outcome, exposure) {
    extra_covars <- if (outcome == "weight_gain") c(Covars, "time_last") else Covars
    fit_terms(DF_totalpop, outcome, exposure, extra_covars, "linear",
              "A2_pBMI_weight_gain")
  })

## Analysis 3: multiplicative pBMI x weight-gain interaction 
# Pregnancy complications, PTB, and delivery week use early-pregnancy gain
# rate. Other outcomes use total gestational weight gain, with time_last added
# to account for the length of the weight-gain observation window.
res_a3_pbmi_gwg_interaction <- crossing(
  outcome = c(binary_otc, contin_otc), pbmi_exposure = pbmi_exposures
) %>%
  pmap_dfr(function(outcome, pbmi_exposure) {
    gwg_exposure <- if (outcome %in% early_gwg_outcomes) "rate_fst" else "weight_gain"
    covars <- model_covars(
      outcome, include_sex = TRUE,
      total_gwg_in_model = identical(gwg_exposure, "weight_gain")
    )
    fit_terms(DF_totalpop, outcome, c(pbmi_exposure, gwg_exposure), covars,
              outcome_model_type(outcome), "A3_pBMI_x_weight_gain", "interaction")
  })

## Analysis 4: multiplicative pBMI x newborn-sex interaction 
# baby_sex is already included through the interaction expression and is not
# duplicated in the covariate vector.
res_a4_pbmi_sex_interaction <- crossing(
  outcome = c(binary_otc, contin_otc), pbmi_exposure = pbmi_exposures
) %>%
  pmap_dfr(function(outcome, pbmi_exposure) {
    fit_terms(DF_totalpop, outcome, c(pbmi_exposure, "baby_sex"), Covars,
              outcome_model_type(outcome), "A4_pBMI_x_newborn_sex", "interaction")
  })

## Analysis 5: pBMI-outcome associations within newborn-sex strata
res_a5_sex_subgroup <- levels(droplevels(DF_totalpop$baby_sex)) %>%
  map_dfr(function(group) {
    run_pbmi_outcomes(
      DF_totalpop %>% filter(baby_sex == group),
      paste0("A5_newborn_sex_subgroup_", group), include_sex = FALSE
    ) %>% mutate(subgroup_variable = "baby_sex", subgroup = group, .after = analysis)
  })

## Save results 
dir.create("result/PART4", recursive = TRUE, showWarnings = FALSE)

result_list <- list(
  A1_pbmi_outcomes = res_a1_pbmi_outcomes,
  A1S_base_same_sample = res_a1s_base_same_sample,
  A1S_weight_gain_adjusted = res_a1s_weight_gain_adjusted,
  A2_pbmi_weight_gain = res_a2_pbmi_weight_gain,
  A3_pbmi_gwg_interaction = res_a3_pbmi_gwg_interaction,
  A4_pbmi_sex_interaction = res_a4_pbmi_sex_interaction,
  A5_sex_subgroups = res_a5_sex_subgroup
)

save(result_list, DF_totalpop,
     file = "result/PART4/1_trad_pbmi_outcomes.RData")

message("Completed: result/PART4/1_trad_pbmi_outcomes.RData")





# Score-level mediation ------------------------------------
# T1 score-level mediation for frozen GDM/HDP DIABLO modules (Figures 5c,d).
#
# Part 1 (Figure 5c): OWO -> bacterial/metabolite score -> outcome.
# Part 2 (Figure 5d): bacterial score -> metabolite score -> outcome, and reverse.
#
# OWO -> B score uses the full processed T1 16S cohort. Analyses containing
# both B and M scores, plus OWO -> M score, use the paired T1 multi-omics
# cohort. Frozen 14.1 B-score parameters are projected without refitting.


if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  invisible(suppressWarnings(try(
    Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE
  )))
}

source("clean_code/14_module_mediation_helpers.R",
  local = environment())

run_14_2_module_mediation <- function(outcomes = c("GDM", "HDP"),
                                      workers = 8L, nsims = 5000L) {
  message("Running score-level module mediation blocks:")
  message("  1/2 OWO -> B score (full T1 16S); OWO -> M score (paired cohort)")
  message("  2/2 B score -> M score -> outcome, and reverse")
  m142_main(outcomes = outcomes, workers = workers,
    nsims = nsims)
}

# Member-level mediation ------------------------------------
# T1 member-level mediation for GDM/HDP DIABLO modules (Figures 5c,d).
#
# Part 1 (Figure 5d): bacterium -> metabolite -> outcome, and the reverse.
# Part 2 (Figure 5c): OWO -> one stable bacterium/metabolite -> outcome.
# Both outcomes retain all planned member tests, including nonsignificant paths.
#
# Bacteria-metabolite pairs use the paired T1 multi-omics cohort. Within the
# OWO-single-feature block, bacterial mediators use the full processed T1 16S
# cohort and metabolite mediators retain the paired cohort. Both parts use
# direct mediation::mediate() quasi-Bayesian simulation defined in 14_pair_mediation_helpers.R.
# Results are saved in separate pair and OWO-feature result folders.


if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  invisible(suppressWarnings(try(
    Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE
  )))
}

source("clean_code/14_pair_mediation_helpers.R", local = environment())

run_14_3_member_mediation <- function(outcomes = c("GDM", "HDP"), workers = 8L,
                                      nsims = 5000L) {
  message("Part 1/2: bidirectional bacterium-metabolite mediation")
  pair_results <- m14_main(
    analysis = "pair", outcomes = outcomes, workers = workers,
    nsims = nsims
  )
  invisible(gc())

  message("Part 2/2: OWO-single-feature mediation (bacteria=full 16S; metabolites=paired)")
  single_results <- m14_main(
    analysis = "single", outcomes = outcomes, workers = workers,
    nsims = nsims
  )

  invisible(list(pair = pair_results, OWO_single_feature = single_results))
}

# Shared entry point: run score-level and member-level analyses in sequence.
parse_mediation_arguments <- function(args = commandArgs(trailingOnly = TRUE)) {
  get_arg <- function(prefix, default) {
    z <- args[startsWith(args, prefix)]
    if (!length(z)) default else sub(prefix, "", z[1], fixed = TRUE)
  }
  list(
    outcomes = strsplit(get_arg("--outcomes=", "GDM,HDP"), ",", fixed = TRUE)[[1]],
    workers = as.integer(get_arg("--workers=", "8")),
    nsims = as.integer(get_arg("--nsims=", "5000"))
  )
}

run_options <- parse_mediation_arguments()
do.call(run_14_2_module_mediation, run_options)
do.call(run_14_3_member_mediation, run_options)
