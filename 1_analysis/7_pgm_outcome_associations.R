##############################################################################################
# Analyses for Figure 6
# Phenotype-genotype mismatch (PGM) and outcome associations
# See plot_fig6.R for the corresponding figure-generation code.
###############################################################################################

rm(list = ls())

library(tidyverse)
library(data.table)
library(writexl)

output_dir <- "result/PART5/PGM"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outcome_vars <- c("GDM", "HDP", "PE", "ICP", "LBW", "MAC", "PTB", "SGA", "LGA", "owo_m36")
outcome_labels <- c(
  GDM = "GDM", HDP = "HDP", PE = "PE", ICP = "ICP", LBW = "LBW",
  MAC = "MAC", PTB = "PTB", SGA = "SGA", LGA = "LGA",
  owo_m36 = "Child OWO\nm36"
)
Covars_common <- c(
  "age", "parity", "smk", "drk", "physic_activity",
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
  "edu", "medical_insurance", "income"
)
pc_vars <- paste0("PC", 1:10)

## 1. Assemble measured-pBMI and multi-omics-pBMI populations ---------------

load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")
load("data/01_thsbc_cleaned_data/outcomes_5212.RData")

pc_file <- "data/genetics/genetic_pcs.eigenvec"
if (!file.exists(pc_file)) stop("Genetic PC file not found: ", pc_file)

pc_data <- fread(pc_file, header = FALSE) %>%
  transmute(
    id = str_pad(as.character(V1), width = 7, side = "left", pad = "0"),
    across(V3:V12)
  )
names(pc_data) <- c("id", pc_vars)

predictions <- readRDS("result/PART5/pbmi_model/pbmi_model_objects.rds")$predictions %>%
  select(id, omics_pbmi = Multi_omics_pBMI_pred_pbmi)

base_population <- total_pop %>%
  select(id, pbmi, PRS) %>%
  inner_join(cov_tab %>% filter(period == "T1"), by = "id") %>%
  inner_join(otc %>% select(id, all_of(outcome_vars)), by = "id") %>%
  inner_join(pc_data, by = "id")

required_measured <- c("id", "pbmi", "PRS", Covars_common, pc_vars)
measured_population <- base_population %>%
  filter(if_all(all_of(required_measured), ~ !is.na(.x)))

omics_population <- base_population %>%
  inner_join(predictions, by = "id") %>%
  filter(if_all(all_of(c(required_measured, "omics_pbmi")), ~ !is.na(.x)))

if (anyDuplicated(measured_population$id) || anyDuplicated(omics_population$id)) {
  stop("Duplicated participant IDs were generated.")
}

## 2. Construct mismatch scores from covariate-adjusted percentile ranks ----

percentile_rank <- function(x) {
  if (length(x) < 2L) stop("At least two observations are required.")
  (rank(x, ties.method = "average") - 1) / (length(x) - 1)
}

derive_pgm <- function(data, phenotype, pgm_type) {
  phenotype_formula <- reformulate(Covars_common, response = phenotype)
  prs_formula <- reformulate(pc_vars, response = "PRS")

  phenotype_model <- lm(phenotype_formula, data = data, na.action = na.fail)
  prs_model <- lm(prs_formula, data = data, na.action = na.fail)

  phenotype_percentile <- percentile_rank(residuals(phenotype_model))
  genetic_percentile <- percentile_rank(residuals(prs_model))
  pgm_score <- 100 * (phenotype_percentile - genetic_percentile)
  cutoffs <- quantile(pgm_score, probs = c(0.25, 0.75), names = FALSE)

  result <- data %>%
    mutate(
      phenotype_residual = residuals(phenotype_model),
      prs_residual = residuals(prs_model),
      phenotype_percentile = phenotype_percentile,
      genetic_percentile = genetic_percentile,
      pgm_score = pgm_score,
      pgm_group = case_when(
        pgm_score < cutoffs[1] ~ "Discordantly low",
        pgm_score >= cutoffs[2] ~ "Discordantly high",
        TRUE ~ "Concordant"
      ),
      pgm_group = factor(
        pgm_group,
        levels = c("Concordant", "Discordantly high", "Discordantly low")
      ),
      pgm_type = pgm_type
    )

  list(
    data = result,
    cutoffs = tibble(
      pgm_type = pgm_type,
      lower_quartile = cutoffs[1],
      upper_quartile = cutoffs[2]
    ),
    phenotype_model = phenotype_model,
    prs_model = prs_model
  )
}

measured_pgm <- derive_pgm(
  measured_population, "pbmi", "Measured pBMI vs genetic BMI"
)
omics_pgm <- derive_pgm(
  omics_population, "omics_pbmi", "Multi-omics-inferred pBMI vs genetic BMI"
)

pgm_data <- bind_rows(measured_pgm$data, omics_pgm$data)
pgm_cutoffs <- bind_rows(measured_pgm$cutoffs, omics_pgm$cutoffs)

## 3. Calculate occurrence and Chi-squared test / Fisher exact tests ------------------

to_binary <- function(x) {
  if (is.factor(x) || is.character(x)) x <- as.character(x)
  as.integer(x)
}

calculate_occurrence <- function(data) {
  data %>%
    select(pgm_type, pgm_group, all_of(outcome_vars)) %>%
    pivot_longer(all_of(outcome_vars), names_to = "outcome", values_to = "value") %>%
    mutate(value = to_binary(value)) %>%
    group_by(pgm_type, pgm_group, outcome) %>%
    summarise(
      n_available = sum(!is.na(value)),
      n_event = sum(value == 1, na.rm = TRUE),
      occurrence_percent = 100 * n_event / n_available,
      .groups = "drop"
    ) %>%
    mutate(outcome_label = unname(outcome_labels[outcome]))
}

pairwise_fisher <- function(data, comparison_group) {
  map_dfr(outcome_vars, function(outcome) {
    test_data <- data %>%
      filter(pgm_group %in% c("Concordant", comparison_group)) %>%
      transmute(
        pgm_group = droplevels(pgm_group),
        outcome_value = to_binary(.data[[outcome]])
      ) %>%
      drop_na()

    contingency <- table(test_data$pgm_group, test_data$outcome_value)
    valid_table <- nrow(contingency) == 2L && ncol(contingency) == 2L
    if (valid_table) {
      chi_check <- suppressWarnings(chisq.test(contingency))
      use_chisq <- all(chi_check$expected >= 5)
      effect_test <- fisher.test(contingency)
      test <- if (use_chisq) {
        chisq.test(contingency)
      } else {
        fisher.test(contingency)
      }
      test_method <- if (use_chisq) "Pearson's Chi-squared test" else "Fisher's exact test"
    } else {
      test <- NULL
      effect_test <- NULL
      test_method <- NA_character_
    }

    tibble(
      outcome = outcome,
      comparison = paste(comparison_group, "vs Concordant"),
      comparison_group = comparison_group,
      test_method = test_method,
      statistic = if (is.null(test) || is.null(test$statistic)) {
        NA_real_
      } else {
        unname(test$statistic)
      },
      p_value = if (is.null(test)) NA_real_ else test$p.value,
      odds_ratio = if (is.null(effect_test)) {
        NA_real_
      } else {
        unname(effect_test$estimate)
      },
      conf_low = if (is.null(effect_test)) {
        NA_real_
      } else {
        effect_test$conf.int[1]
      },
      conf_high = if (is.null(effect_test)) {
        NA_real_
      } else {
        effect_test$conf.int[2]
      }
    )
  })
}

comparison_results <- pgm_data %>%
  group_split(pgm_type) %>%
  map_dfr(function(cohort_data) {
    pgm_type_value <- as.character(cohort_data$pgm_type[1])
    bind_rows(
      pairwise_fisher(cohort_data, "Discordantly high"),
      pairwise_fisher(cohort_data, "Discordantly low")
    ) %>%
      mutate(pgm_type = pgm_type_value, .before = 1)
  }) %>%
  group_by(pgm_type, comparison_group) %>%
  mutate(p_adjust_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    outcome_label = unname(outcome_labels[outcome])
  )

outcome_occurrence <- calculate_occurrence(pgm_data) %>%
  left_join(
    comparison_results %>%
      select(pgm_type, outcome, comparison_group, p_value, p_adjust_bh, significance),
    by = c("pgm_type", "outcome", "pgm_group" = "comparison_group")
  )

group_summary <- pgm_data %>%
  count(pgm_type, pgm_group, name = "n")

## 4. Prespecified main, subgroup, and offspring-sex interaction analyses ----

subgroup_outcomes <- c(
  "GDM", "HDP", "PE", "ICP", "LBW", "MAC", "PTB", "SGA", "LGA", "owo_m36"
)

# Use one PRS threshold from the measured-pBMI population for both PGM datasets.
prs_median <- median(measured_population$PRS, na.rm = TRUE)

analysis_data <- pgm_data %>%
  mutate(
    pgm_per10 = pgm_score / 10,
    pbmi_group = case_when(
      pbmi < 18.5 ~ "Underweight",
      pbmi < 24 ~ "Normal",
      TRUE ~ "OWO"
    ),
    pbmi_group = factor(pbmi_group, levels = c("Normal", "Underweight", "OWO")),
    pbmi_owo_group = factor(
      if_else(pbmi >= 24, "OWO", "non-OWO"),
      levels = c("non-OWO", "OWO")
    ),
    prs_group = factor(
      if_else(PRS >= prs_median, "HPRS", "LPRS"),
      levels = c("LPRS", "HPRS")
    ),
    baby_sex = factor(baby_sex),
    pgm_group = factor(
      pgm_group,
      levels = c("Concordant", "Discordantly high", "Discordantly low")
    )
  )

for (outcome_name in subgroup_outcomes) {
  analysis_data[[outcome_name]] <- to_binary(analysis_data[[outcome_name]])
}

varying_terms <- function(data, terms) {
  terms[map_lgl(terms, function(term) {
    n_distinct(data[[term]][!is.na(data[[term]])]) > 1L
  })]
}

fit_logistic <- function(data, outcome, exposure, covariates, analysis, subgroup = "All") {
  exposure_variables <- all.vars(as.formula(paste("~", exposure)))
  required_variables <- unique(c(outcome, exposure_variables, covariates))

  model_data <- data %>%
    select(all_of(required_variables)) %>%
    drop_na()

  n_total <- nrow(model_data)
  n_event <- sum(model_data[[outcome]] == 1L)

  invalid_data <- n_total == 0L ||
    n_distinct(model_data[[outcome]]) < 2L ||
    any(map_lgl(exposure_variables, function(variable) {
      n_distinct(model_data[[variable]]) < 2L
    }))

  if (invalid_data) {
    return(tibble(
      outcome = outcome, analysis = analysis, subgroup = subgroup,
      term = NA_character_, term_type = NA_character_,
      beta = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      p_value = NA_real_, n = n_total, events = n_event,
      status = "Outcome or exposure has no variation"
    ))
  }

  usable_covariates <- varying_terms(model_data, covariates)
  model_formula <- as.formula(paste(
    outcome, "~", paste(c(exposure, usable_covariates), collapse = " + ")
  ))

  model <- tryCatch(
    glm(model_formula, family = binomial(), data = model_data),
    error = function(error) error
  )

  if (inherits(model, "error")) {
    return(tibble(
      outcome = outcome, analysis = analysis, subgroup = subgroup,
      term = NA_character_, term_type = NA_character_,
      beta = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      p_value = NA_real_, n = n_total, events = n_event,
      status = conditionMessage(model)
    ))
  }

  coefficient_table <- summary(model)$coefficients
  keep <- str_detect(rownames(coefficient_table), "pgm_per10|pgm_group")
  coefficient_table <- coefficient_table[keep, , drop = FALSE]

  tibble(
    outcome = outcome,
    analysis = analysis,
    subgroup = subgroup,
    term = rownames(coefficient_table),
    term_type = if_else(str_detect(term, ":"), "interaction", "main effect"),
    beta = coefficient_table[, "Estimate"],
    OR = exp(beta),
    CI_low = exp(beta - 1.96 * coefficient_table[, "Std. Error"]),
    CI_high = exp(beta + 1.96 * coefficient_table[, "Std. Error"]),
    p_value = coefficient_table[, "Pr(>|z|)"],
    n = n_total,
    events = n_event,
    status = if_else(model$converged, "OK", "Model did not converge")
  )
}

run_outcome_models <- function(data, pgm_label) {
  map_dfr(subgroup_outcomes, function(outcome) {
    bind_rows(
      # Primary associations: retain adjustment for measured pBMI group.
      fit_logistic(
        data, outcome, "pgm_per10", c(Covars_common, "pbmi_group"),
        "Continuous PGM (per 10 units), adjusted for pBMI group"
      ),
      fit_logistic(
        data, outcome, "pgm_group", c(Covars_common, "pbmi_group"),
        "Categorical PGM, adjusted for pBMI group"
      ),

      # Subgroup associations
      map_dfr(levels(data$pbmi_owo_group), function(group) {
        fit_logistic(
          filter(data, pbmi_owo_group == group),
          outcome, "pgm_group", Covars_common,
          "pBMI subgroup", group
        )
      }),
      map_dfr(levels(data$prs_group), function(group) {
        fit_logistic(
          filter(data, prs_group == group),
          outcome, "pgm_group", Covars_common,
          "PRS subgroup", group
        )
      }),

      # Exploratory effect modification by offspring sex. Continuous PGM gives
      # one interaction coefficient and preserves information.
      fit_logistic(
        data, outcome, "pgm_per10 * baby_sex",
        c(Covars_common, "pbmi_group"),
        "Continuous PGM-by-offspring-sex interaction"
      ),
      map_dfr(levels(data$baby_sex), function(group) {
        fit_logistic(
          filter(data, baby_sex == group),
          outcome, "pgm_group", c(Covars_common, "pbmi_group"),
          "Offspring-sex subgroup", paste0("baby_sex=", group)
        )
      })
    ) %>%
      mutate(pgm_type = pgm_label, .before = 1)
  })
}

measured_subgroup_results <- analysis_data %>%
  filter(pgm_type == "Measured pBMI vs genetic BMI") %>%
  run_outcome_models("Measured pBMI vs genetic BMI")

omics_subgroup_results <- analysis_data %>%
  filter(pgm_type == "Multi-omics-inferred pBMI vs genetic BMI") %>%
  run_outcome_models("Multi-omics-inferred pBMI vs genetic BMI")

add_significance <- function(data) {
  data %>%
    mutate(
      significance = case_when(
        is.na(p_value) ~ NA_character_,
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

measured_subgroup_results <- add_significance(measured_subgroup_results)
omics_subgroup_results <- add_significance(omics_subgroup_results)
subgroup_interaction_results <- bind_rows(
  measured_subgroup_results, omics_subgroup_results
)

write_model_workbook <- function(results, path) {
  write_xlsx(
    list(
      all_results = results,
      continuous_PGM = filter(
        results,
        analysis == "Continuous PGM (per 10 units), adjusted for pBMI group"
      ),
      categorical_PGM = filter(
        results,
        analysis == "Categorical PGM, adjusted for pBMI group"
      ),
      pBMI_subgroups = filter(results, analysis == "pBMI subgroup"),
      PRS_subgroups = filter(results, analysis == "PRS subgroup"),
      sex_interaction = filter(
        results,
        analysis == "Continuous PGM-by-offspring-sex interaction"
      ),
      sex_subgroups = filter(
        results,
        analysis == "Offspring-sex subgroup"
      )
    ),
    path
  )
}

write_model_workbook(
  measured_subgroup_results,
  file.path(output_dir, "measured_pbmi_pgm_subgroup_interaction.xlsx")
)
write_model_workbook(
  omics_subgroup_results,
  file.path(output_dir, "omics_pbmi_pgm_subgroup_interaction.xlsx")
)

saveRDS(
  list(
    prs_median = prs_median,
    measured = measured_subgroup_results,
    multi_omics = omics_subgroup_results
  ),
  file.path(output_dir, "PGM_subgroup_interaction_analysis.rds")
)

## 5. Save analysis results --------------------------------------------------

write_xlsx(
  list(
    group_summary = group_summary,
    cutoffs = pgm_cutoffs,
    occurrence = outcome_occurrence,
    fisher_comparisons = comparison_results,
    participant_groups = pgm_data %>%
      select(id, pgm_type, pgm_score, pgm_group,
             phenotype_percentile, genetic_percentile)
  ),
  file.path(output_dir, "PGM_outcome_analysis.xlsx")
)

write_csv(outcome_occurrence, file.path(output_dir, "PGM_outcome_occurrence.csv"))
write_csv(comparison_results, file.path(output_dir, "PGM_fisher_comparisons.csv"))

saveRDS(
  list(
    group_summary = group_summary,
    cutoffs = pgm_cutoffs,
    outcome_occurrence = outcome_occurrence,
    comparison_results = comparison_results,
    participant_data = pgm_data,
    models = list(
      measured_phenotype = measured_pgm$phenotype_model,
      measured_prs = measured_pgm$prs_model,
      omics_phenotype = omics_pgm$phenotype_model,
      omics_prs = omics_pgm$prs_model
    )
  ),
  file.path(output_dir, "PGM_outcome_analysis.rds")
)

print(group_summary)
message("PGM outcome analysis completed. Results saved to ", output_dir, ".")
