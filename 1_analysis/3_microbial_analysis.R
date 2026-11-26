##############################################################################################
# Analyses for Figure 4
# Characterization of gut microbial signatures associated with maternal pre-pregnancy adiposity.
# See plot_fig4.R for the corresponding figure-generation code.
###############################################################################################

# Differential microbiota analyses in THSBC -----------------------------------------
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
load('data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData')
load('data/01_thsbc_cleaned_data/gut_microbe_MetaG_1250_473.RData')
load('data/01_thsbc_cleaned_data/total_pop_5212.RData')
load('data/01_thsbc_cleaned_data/covariates_5212.RData')
load('data/01_thsbc_cleaned_data/wt_omics.rdata')

DF_16SITS <- total_pop %>%
  inner_join(bact_fung_clrz, by = "id") %>%
  inner_join(cov_tab, by = c("id", "period")) %>%
  inner_join(wt_16SITS, by = c("id", "period"))

DF_MetaG <- total_pop %>%
  inner_join(mgs_clrz, by = "id") %>%
  inner_join(cov_tab, by = c("id", "period")) %>%
  inner_join(wt_MetaG, by = c("id", "period"))


# Associations between pBMI and gut microbial features
Covars_common <- c(
  "age", "parity", "smk", "drk", "physic_activity",
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
  "edu", "medical_insurance", "income"
)

periods <- c("T1", "T2", "T3")

pbmi_unique <- total_pop %>% distinct(id, pbmi)
pbmi_mean <- mean(pbmi_unique$pbmi, na.rm = TRUE)
pbmi_sd <- sd(pbmi_unique$pbmi, na.rm = TRUE)

factor_covars <- c(
  "parity", "smk", "drk", "fd_grain", "fd_meat", "fd_vege",
  "fd_fruit", "fd_egg", "fd_milkpro", "edu", "medical_insurance",
  "income", "antibiotic"
)

prepare_pbmi_exposures <- function(data) {
  data <- data %>%
    mutate(
      pbmi_z = (pbmi - pbmi_mean) / pbmi_sd,
      across(all_of(factor_covars), as.factor)
    )

  data$pbmi_owo_normal <- factor(
    if_else(
      data$pbmi_cut3 %in% c("Normal", "OWO"),
      as.character(data$pbmi_cut3),
      NA_character_
    ),
    levels = c("Normal", "OWO")
  )

  contrast_definitions <- tribble(
    ~exposure,                 ~reference,   ~comparison,
    "pbmi_overweight_normal",  "Normal",     "Overweight",
    "pbmi_obesity_normal",     "Normal",     "Obese",
    "pbmi_underweight_normal", "Normal",     "Underweight",
    "pbmi_obesity_overweight", "Overweight", "Obese"
  )

  contrast_columns <- pmap(
    contrast_definitions,
    function(exposure, reference, comparison) {
      factor(
        if_else(
          data$pbmi_cut4 %in% c(reference, comparison),
          as.character(data$pbmi_cut4),
          NA_character_
        ),
        levels = c(reference, comparison)
      )
    }
  ) %>%
    set_names(contrast_definitions$exposure)

  bind_cols(data, contrast_columns)
}

DF_16SITS <- prepare_pbmi_exposures(DF_16SITS) %>%
  mutate(
    bact_batch = as.factor(bact_batch),
    fung_batch = as.factor(fung_batch)
  )

DF_MetaG <- prepare_pbmi_exposures(DF_MetaG) %>%
  mutate(mgs_batch = as.factor(mgs_batch))

comparison_config <- tribble(
  ~exposure,                 ~label,
  "pbmi_owo_normal",         "OWO vs normal weight",
  "pbmi_overweight_normal",  "Overweight vs normal weight",
  "pbmi_obesity_normal",     "Obesity vs normal weight",
  "pbmi_underweight_normal", "Underweight vs normal weight",
  "pbmi_obesity_overweight", "Obesity vs overweight"
) %>%
  mutate(
    target_population_size = map_int(exposure, ~ switch(
      .x,
      pbmi_owo_normal = sum(
        total_pop$pbmi_cut3 %in% c("Normal", "OWO"), na.rm = TRUE
      ),
      pbmi_overweight_normal = sum(
        total_pop$pbmi_cut4 %in% c("Normal", "Overweight"), na.rm = TRUE
      ),
      pbmi_obesity_normal = sum(
        total_pop$pbmi_cut4 %in% c("Normal", "Obese"), na.rm = TRUE
      ),
      pbmi_underweight_normal = sum(
        total_pop$pbmi_cut4 %in% c("Normal", "Underweight"), na.rm = TRUE
      ),
      pbmi_obesity_overweight = sum(
        total_pop$pbmi_cut4 %in% c("Overweight", "Obese"), na.rm = TRUE
      )
    ))
  )

omics_config <- list(
  `16S_g` = list(
    data = DF_16SITS,
    features = bact_uni,
    batch = "bact_batch",
    info = info_bact_fung,
    info_key = "ID"
  ),
  ITS_g = list(
    data = DF_16SITS,
    features = fung_uni,
    batch = "fung_batch",
    info = info_bact_fung,
    info_key = "ID"
  ),
  MetaG_g = list(
    data = DF_MetaG,
    features = mgs_g_uni,
    batch = "mgs_batch",
    info = info_mgs,
    info_key = "ID"
  ),
  MetaG_s = list(
    data = DF_MetaG,
    features = mgs_s_uni,
    batch = "mgs_batch",
    info = info_mgs,
    info_key = "ID"
  )
)

run_microbiome_layer <- function(omics_name, config) {
  message("Running differential analysis for ", omics_name, "...")

  covars <- c(Covars_common, config$batch, "wk", "antibiotic")
  period_covars <- set_names(
    map(periods, function(period) {
      if (omics_name %in% c("MetaG_g", "MetaG_s") &&
          period %in% c("T2", "T3")) {
        setdiff(covars, "mgs_batch")
      } else {
        covars
      }
    }),
    periods
  )

  continuous <- func_run_omics_models(
    data = config$data,
    feature_cols = config$features,
    exposure = "pbmi_z",
    covars = covars,
    periods = periods,
    exposure_label = "pBMI per 1 SD",
    target_population_size = nrow(total_pop),
    period_covars = period_covars
  ) %>%
    func_append_feature_info(config$info, config$info_key, omics_name)

  categorical <- pmap_dfr(
    comparison_config,
    function(exposure, label, target_population_size) {
      func_run_omics_models(
      data = config$data,
      feature_cols = config$features,
      exposure = exposure,
      covars = covars,
      periods = periods,
      exposure_label = label,
      target_population_size = target_population_size,
      period_covars = period_covars
    )
    }
  ) %>%
    func_append_feature_info(config$info, config$info_key, omics_name)

  list(
    continuous = continuous,
    categorical = categorical
  )
}

microbiome_results <- imap(
  omics_config,
  ~ run_microbiome_layer(omics_name = .y, config = .x)
)

res_pbmi_microbe_continuous <- map_dfr(microbiome_results, "continuous")
res_pbmi_microbe_categorical <- map_dfr(microbiome_results, "categorical")

dir.create("result/PART3", recursive = TRUE, showWarnings = FALSE)

save(
  res_pbmi_microbe_continuous,
  res_pbmi_microbe_categorical,
  microbiome_results,
  comparison_config,
  file = "result/PART3/1_diff_microbe.rdata"
)

message("Differential microbiome analyses completed.")


# Associations between other anthropometric traits and gut microbial features
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

prepare_secondary_anthropometrics <- function(data) {
  for (config_index in seq_len(nrow(secondary_continuous_config))) {
    exposure_raw <- secondary_continuous_config$exposure_raw[[config_index]]
    exposure_z <- secondary_continuous_config$exposure_z[[config_index]]
    exposure_mean <- mean(total_pop[[exposure_raw]], na.rm = TRUE)
    exposure_sd <- sd(total_pop[[exposure_raw]], na.rm = TRUE)
    data[[exposure_z]] <- (data[[exposure_raw]] - exposure_mean) / exposure_sd
  }

  data %>%
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
}

DF_16SITS <- prepare_secondary_anthropometrics(DF_16SITS)
DF_MetaG <- prepare_secondary_anthropometrics(DF_MetaG)
omics_config$`16S_g`$data <- DF_16SITS
omics_config$ITS_g$data <- DF_16SITS
omics_config$MetaG_g$data <- DF_MetaG
omics_config$MetaG_s$data <- DF_MetaG

run_secondary_anthropometric_layer <- function(omics_name, config) {
  message("Running secondary anthropometric analysis for ", omics_name, "...")

  covars <- c(Covars_common, config$batch, "wk", "antibiotic")
  period_covars <- set_names(
    map(periods, function(period) {
      if (omics_name %in% c("MetaG_g", "MetaG_s") &&
          period %in% c("T2", "T3")) {
        setdiff(covars, "mgs_batch")
      } else {
        covars
      }
    }),
    periods
  )

  continuous <- pmap_dfr(
    secondary_continuous_config,
    function(exposure_raw, exposure_z, label, target_population_size) {
      func_run_omics_models(
        data = config$data,
        feature_cols = config$features,
        exposure = exposure_z,
        covars = covars,
        periods = periods,
        exposure_label = label,
        target_population_size = target_population_size,
        period_covars = period_covars
      )
    }
  ) %>%
    func_append_feature_info(config$info, config$info_key, omics_name)

  binary <- pmap_dfr(
    secondary_binary_config,
    function(exposure, reference, comparison, label, target_population_size) {
      func_run_omics_models(
        data = config$data,
        feature_cols = config$features,
        exposure = exposure,
        covars = covars,
        periods = periods,
        exposure_label = label,
        target_population_size = target_population_size,
        period_covars = period_covars
      )
    }
  ) %>%
    func_append_feature_info(config$info, config$info_key, omics_name)

  list(continuous = continuous, binary = binary)
}

secondary_microbiome_results <- imap(
  omics_config,
  ~ run_secondary_anthropometric_layer(omics_name = .y, config = .x)
)

res_secondary_microbe_continuous <- map_dfr(
  secondary_microbiome_results,
  "continuous"
)
res_secondary_microbe_binary <- map_dfr(
  secondary_microbiome_results,
  "binary"
)

save(
  res_secondary_microbe_binary,
  res_secondary_microbe_continuous,
  file = "result/PART3/1-2_diff_microbe_secdwt.rdata"
)

message("Secondary anthropometric-microbiome analyses completed.")





# Differential microbiota analyses in HBC -----------------------------------------
## load data
load('data/02_hbc_cleaned_data/gut_microbe_fung_393_158.RData')
load('data/02_hbc_cleaned_data/gut_microbe_MetaG_393_158.RData')
load('data/02_hbc_cleaned_data/total_pop_hbc_393.Rdata')
load('data/02_hbc_cleaned_data/covariates_393.RData')

DF_ITS <- total_pop_hbc %>%
  inner_join(fung_clrz_hbc, by = "id") %>%
  inner_join(cov_tab, by = c("id", "period")) 

DF_MetaG <- total_pop_hbc %>%
  inner_join(mgs_clrz_hbc, by = "id") %>%
  inner_join(cov_tab, by = c("id", "period")) 


# Regression of pBMI and other anthropometric traits  

Covars <- c(
  "age", "parity", "smk", "drk", "physic_activity",
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
  "edu", "medical_insurance", "income", "wk", "antibiotic"
)
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

pbmi_comparison_config <- tribble(
  ~exposure,                 ~reference,   ~comparison,   ~label,
  "pbmi_owo_normal",         "Normal",     "OWO",         "OWO vs normal weight",
  "pbmi_overweight_normal",  "Normal",     "Overweight",  "Overweight vs normal weight",
  "pbmi_obesity_normal",     "Normal",     "Obese",       "Obesity vs normal weight",
  "pbmi_underweight_normal", "Normal",     "Underweight", "Underweight vs normal weight",
  "pbmi_obesity_overweight", "Overweight", "Obese",       "Obesity vs overweight"
)

binary_anthropometric_config <- tribble(
  ~exposure,    ~label,
  "wc_t1_cut",  "WC abdominal obesity vs normal",
  "whr_t1_cut", "WHR abdominal obesity vs normal"
)

prepare_hbc_exposures <- function(data) {
  for (config_index in seq_len(nrow(continuous_config))) {
    exposure_raw <- continuous_config$exposure_raw[[config_index]]
    exposure_z <- continuous_config$exposure_z[[config_index]]
    exposure_mean <- mean(total_pop_hbc[[exposure_raw]], na.rm = TRUE)
    exposure_sd <- sd(total_pop_hbc[[exposure_raw]], na.rm = TRUE)
    data[[exposure_z]] <- (data[[exposure_raw]] - exposure_mean) / exposure_sd
  }

  data$pbmi_owo_normal <- factor(
    if_else(
      data$pbmi_cut3 %in% c("Normal", "OWO"),
      as.character(data$pbmi_cut3),
      NA_character_
    ),
    levels = c("Normal", "OWO")
  )

  contrast_columns <- pbmi_comparison_config %>%
    filter(exposure != "pbmi_owo_normal") %>%
    pmap(function(exposure, reference, comparison, label) {
      factor(
        if_else(
          data$pbmi_cut4 %in% c(reference, comparison),
          as.character(data$pbmi_cut4),
          NA_character_
        ),
        levels = c(reference, comparison)
      )
    }) %>%
    set_names(pbmi_comparison_config$exposure[-1])

  bind_cols(data, contrast_columns) %>%
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
}

DF_ITS <- prepare_hbc_exposures(DF_ITS)
DF_MetaG <- prepare_hbc_exposures(DF_MetaG)

info_fung_hbc_keyed <- info_fung_its_hbc %>%
  mutate(feature = paste0("fung_", ID))

info_mgs_hbc_keyed <- info_mgs %>%
  mutate(
    feature = case_when(
      level == "Genus" ~ paste0("mgs_", ID),
      level == "Species" ~ paste0("mgs_", ID),
      TRUE ~ NA_character_
    )
  )

omics_config_hbc <- list(
  ITS_g = list(
    data = DF_ITS,
    features = fung_uni,
    info = info_fung_hbc_keyed,
    info_key = "feature"
  ),
  MetaG_g = list(
    data = DF_MetaG,
    features = mgs_g_uni,
    info = info_mgs_hbc_keyed,
    info_key = "feature"
  ),
  MetaG_s = list(
    data = DF_MetaG,
    features = mgs_s_uni,
    info = info_mgs_hbc_keyed,
    info_key = "feature"
  )
)

run_hbc_microbiome_layer <- function(omics_name, config) {
  message("Running unweighted HBC analysis for ", omics_name, "...")

  continuous <- pmap_dfr(
    continuous_config,
    function(exposure_raw, exposure_z, label) {
      func_run_omics_models_unweighted(
        data = config$data,
        feature_cols = config$features,
        exposure = exposure_z,
        covars = Covars,
        periods = periods,
        exposure_label = label
      )
    }
  ) %>%
    func_append_feature_info(config$info, config$info_key, omics_name)

  pbmi_categorical <- pmap_dfr(
    pbmi_comparison_config,
    function(exposure, reference, comparison, label) {
      func_run_omics_models_unweighted(
        data = config$data,
        feature_cols = config$features,
        exposure = exposure,
        covars = Covars,
        periods = periods,
        exposure_label = label
      )
    }
  ) %>%
    func_append_feature_info(config$info, config$info_key, omics_name)

  anthropometric_binary <- pmap_dfr(
    binary_anthropometric_config,
    function(exposure, label) {
      func_run_omics_models_unweighted(
        data = config$data,
        feature_cols = config$features,
        exposure = exposure,
        covars = Covars,
        periods = periods,
        exposure_label = label
      )
    }
  ) %>%
    func_append_feature_info(config$info, config$info_key, omics_name)

  list(
    continuous = continuous,
    pbmi_categorical = pbmi_categorical,
    anthropometric_binary = anthropometric_binary
  )
}

hbc_microbiome_results <- imap(
  omics_config_hbc,
  ~ run_hbc_microbiome_layer(omics_name = .y, config = .x)
)

res_hbc_microbe_continuous <- map_dfr(hbc_microbiome_results, "continuous")
res_hbc_microbe_pbmi_categorical <- map_dfr(
  hbc_microbiome_results, "pbmi_categorical"
)
res_hbc_microbe_anthropometric_binary <- map_dfr(
  hbc_microbiome_results, "anthropometric_binary"
)
dir.create("result/PART3", recursive = TRUE, showWarnings = FALSE)

save(
  res_hbc_microbe_continuous,
  res_hbc_microbe_pbmi_categorical,
  res_hbc_microbe_anthropometric_binary,
  hbc_microbiome_results,
  continuous_config,
  pbmi_comparison_config,
  binary_anthropometric_config,
  file = "result/PART3/2_diff_microb_HBC.rdata"
)

message("HBC differential microbiome analyses completed.")


# Clinical relevance in THSBC----------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(survey)
  library(writexl)
})

# Input 
load("data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData")
load('data/01_thsbc_cleaned_data/gut_microbe_MetaG_1250_473.RData')
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")
load("data/01_thsbc_cleaned_data/wt_omics.rdata")
load("data/01_thsbc_cleaned_data/clinical_traits_t1t2_4562_4874.RData")
load("result/PART3/thsbc_hbc_rep_genus_species_owo_sig_sets.rdata")
# Precomputed THSBC discovery sets for OWO versus normal weight.
# Bacterial and fungal genera were selected separately at
# Benjamini-Hochberg-adjusted P < 0.05 in at least one timepoint-specific
# weighted linear model (T1, T2 or T3) or the pooled weighted
# linear mixed-effects model.
# The deduplicated discovery sets contain 74 bacterial genera
# (sig_16s_sets$sigid_16s_owo) and 6 fungal genera
# (sig_its_sets$sigid_its_owo).
# The 23 metagenomically supported bacterial genera used below are
# supplied separately through owo_rep_id$rep_genus.
load("result/PART3/1_diff_16s_its_sig_sets_thsbc.rdata")

owo_16s <- unique(owo_rep_id$rep_genus)
owo_its <- unique(sig_its_sets$sigid_its_owo)
owo_mgs_g <- unique(owo_rep_id$rep_genus_mgs)
owo_mgs_s <- unique(owo_rep_id$rep_species)
target_16s_its <- unique(c(owo_16s, owo_its))
# Use the stored CLR-transformed abundances and require all selected taxa to be available
target_mgs_g <- owo_mgs_g
target_mgs_s <- owo_mgs_s
target_metag <- unique(c(target_mgs_g, target_mgs_s))
target_microbiota <- unique(c(target_16s_its, target_metag))
target_population_size <- 5212L
missing_metag_features <- setdiff(c(owo_mgs_g, owo_mgs_s), names(mgs_clrz))

feature_availability <- tibble(
  microbiome_type = c(
    rep("MetaG_g", length(owo_mgs_g)),
    rep("MetaG_s", length(owo_mgs_s))
  ),
  microbe = c(owo_mgs_g, owo_mgs_s)
) %>%
  mutate(
    available_in_saved_clr = microbe %in% names(mgs_clrz),
    exclusion_reason = if_else(
      available_in_saved_clr,
      NA_character_,
      "not_present_in_saved_mgs_clrz"
    )
  ) %>%
  left_join(
    info_mgs %>% select(microbe = ID, microbe_name = name),
    by = "microbe"
  )

stopifnot(
  length(owo_16s) == 23L,
  length(owo_its) == 6L,
  length(owo_mgs_g) == 23L,
  length(owo_mgs_s) == 37L,
  !anyDuplicated(owo_16s),
  !anyDuplicated(owo_mgs_g),
  !anyDuplicated(owo_mgs_s),
  all(target_16s_its %in% names(bact_fung_clrz)),
  all(c(owo_mgs_g, owo_mgs_s) %in% info_mgs$ID),
  all(target_mgs_g %in% names(mgs_clrz)),
  all(target_mgs_s %in% names(mgs_clrz))
)

if (length(missing_metag_features) > 0L) {
  message(
    "Excluded from clinical-trait models because the feature is absent from ",
    "the saved THSBC CLR matrix: ",
    paste(missing_metag_features, collapse = ", ")
  )
}

microbiota_annotation <- bind_rows(
  info_bact_fung %>%
    filter(ID %in% target_16s_its) %>%
    transmute(
      microbiome_type = if_else(ID %in% owo_16s, "16S_g", "ITS_g"),
      microbe = ID,
      microbe_name = name
    ),
  info_mgs %>%
    filter(ID %in% target_metag) %>%
    transmute(
      microbiome_type = if_else(ID %in% target_mgs_g, "MetaG_g", "MetaG_s"),
      microbe = ID,
      microbe_name = name
    )
) %>%
  distinct(microbiome_type, microbe, .keep_all = TRUE)

clini_var_t1 <- c(
  "SBP", "DBP", "ALB", "TP", "a_HBDH", "A/G Ratio", "CHE", 
  "GGT", "ALT", "AST", "IBIL", "ALP", "GLB", "LAC", "LDH", "ADA",
  "DBIL", "TBIL", "Cr", "Ccr", "BUN", "UA", "CysC", "FPG",
  "FT3", "FT4", "Tg", "TPOAb", "TSH"
)
clini_var_t2 <- c("OGTT_0h", "OGTT_1h", "OGTT_2h")

covariates_common <- c(
  "age", "parity", "smk", "drk", "physic_activity",
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro",
  "edu", "medical_insurance", "income", "wk", "antibiotic"
)

factor_covariates <- c(
  "parity", "smk", "drk", "fd_grain", "fd_meat", "fd_vege",
  "fd_fruit", "fd_egg", "fd_milkpro", "edu", "medical_insurance",
  "income", "antibiotic", "bact_batch", "fung_batch", "mgs_batch"
)

# Build analysis datasets 
make_microbiome_period_data <- function(microbiome_period) {
  bact_fung_clrz %>%
    select(
      id, period, bact_batch, fung_batch,
      all_of(target_16s_its)
    ) %>%
    inner_join(cov_tab, by = c("id", "period")) %>%
    inner_join(wt_16SITS, by = c("id", "period")) %>%
    filter(period == microbiome_period) %>%
    mutate(
      across(any_of(factor_covariates), as.factor),
      microbiome_period = microbiome_period
    ) %>%
    select(-period)
}

make_metag_period_data <- function(microbiome_period) {
  mgs_clrz %>%
    select(id, period, mgs_batch, all_of(target_metag)) %>%
    inner_join(cov_tab, by = c("id", "period")) %>%
    inner_join(wt_MetaG, by = c("id", "period")) %>%
    filter(period == microbiome_period) %>%
    mutate(
      across(any_of(factor_covariates), as.factor),
      microbiome_period = microbiome_period
    ) %>%
    select(-period)
}

microbiome_t1 <- make_microbiome_period_data("T1")
microbiome_t2 <- make_microbiome_period_data("T2")
metag_t1 <- make_metag_period_data("T1")
metag_t2 <- make_metag_period_data("T2")

analysis_t1_micro_t1_clinical <- microbiome_t1 %>%
  inner_join(clini_traits_T1, by = "id")

analysis_t1_micro_t2_clinical <- microbiome_t1 %>%
  inner_join(clini_traits_T2, by = "id")

analysis_t1_metag_t1_clinical <- metag_t1 %>%
  inner_join(clini_traits_T1, by = "id")

analysis_t1_metag_t2_clinical <- metag_t1 %>%
  inner_join(clini_traits_T2, by = "id")

# Weighted linear regression

fit_one_association <- function(data, clinical_trait, microbe, clinical_period, microbiome_period, target_n = target_population_size) {
  microbiome_type <- case_when(
    microbe %in% owo_16s ~ "16S_g",
    microbe %in% owo_its ~ "ITS_g",
    microbe %in% owo_mgs_g ~ "MetaG_g",
    microbe %in% owo_mgs_s ~ "MetaG_s",
    TRUE ~ NA_character_
  )
  if (is.na(microbiome_type)) {
    stop("Unknown microbiome feature: ", microbe)
  }
  batch_variable <- switch(
    microbiome_type,
    `16S_g` = "bact_batch",
    ITS_g = "fung_batch",
    MetaG_g = "mgs_batch",
    MetaG_s = "mgs_batch"
  )
  model_covariates <- c(covariates_common, batch_variable)
  required_columns <- unique(c(
    "id", "prob", clinical_trait, microbe, model_covariates
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
      .microbe_exposure = .data[[microbe]]
    )

  empty_result <- function(status_text) {
    tibble(
      clinical_trait = clinical_trait,
      microbiome_type = microbiome_type,
      microbe = microbe,
      clinical_period = clinical_period,
      microbiome_period = microbiome_period,
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
      n_distinct(model_data[[microbe]]) < 2L) {
    return(empty_result("no_variation"))
  }
  if (any(model_data$prob <= 0 | model_data$prob > 1)) {
    return(empty_result("invalid_sampling_probability"))
  }

  # A batch/covariate represented by only one level in a period cannot be
  # estimated (e.g. all T2 MetaG samples are from mgs_batch 2).
  estimable_covariates <- model_covariates[
    map_lgl(
      model_covariates,
      ~ n_distinct(model_data[[.x]], na.rm = TRUE) > 1L
    )
  ]

  model_formula <- reformulate(
    termlabels = c(".microbe_exposure", estimable_covariates),
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
      ".microbe_exposure",
      rownames(coefficient_table)
    )

    if (is.na(coefficient_row)) {
      return(empty_result("microbe_coefficient_not_estimable"))
    }

    estimate <- coefficient_table[coefficient_row, "Estimate"]
    standard_error <- coefficient_table[coefficient_row, "Std. Error"]
    p_column <- grep("Pr\\(", colnames(coefficient_table), value = TRUE)[1L]

    tibble(
      clinical_trait = clinical_trait,
      microbiome_type = microbiome_type,
      microbe = microbe,
      clinical_period = clinical_period,
      microbiome_period = microbiome_period,
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

run_association_grid <- function(data, clinical_traits, clinical_period, microbiome_period, microbiota_features) {
  tidyr::crossing(
    clinical_trait = clinical_traits,
    microbe = microbiota_features
  ) %>%
    mutate(
      result = map2(
        clinical_trait,
        microbe,
        ~ fit_one_association(
          data = data,
          clinical_trait = .x,
          microbe = .y,
          clinical_period = clinical_period,
          microbiome_period = microbiome_period
        )
      )
    ) %>%
    select(result) %>%
    unnest(result) %>%
    group_by(
      clinical_period, microbiome_period, microbiome_type, clinical_trait
    ) %>%
    mutate(
      p_fdr = p.adjust(p_value, method = "fdr"),
      p_bonferroni = p.adjust(p_value, method = "bonferroni")
    ) %>%
    ungroup() %>%
    left_join(
      microbiota_annotation,
      by = c("microbiome_type", "microbe")
    ) %>%
    relocate(microbe_name, .after = microbe) %>%
    arrange(
      clinical_period, microbiome_period, microbiome_type,
      clinical_trait, p_value
    )
}

# Analyses 

message("Running T1 16S/ITS -> T1 clinical traits...")
results_t1_16sits_t1_clinical <- run_association_grid(
  data = analysis_t1_micro_t1_clinical,
  clinical_traits = clini_var_t1,
  clinical_period = "T1",
  microbiome_period = "T1",
  microbiota_features = target_16s_its
)

message("Running T1 16S/ITS -> T2 OGTT traits...")
results_t1_16sits_t2_clinical <- run_association_grid(
  data = analysis_t1_micro_t2_clinical,
  clinical_traits = clini_var_t2,
  clinical_period = "T2",
  microbiome_period = "T1",
  microbiota_features = target_16s_its
)

message("Running T1 MetaG genus/species -> T1 clinical traits...")
results_t1_metag_t1_clinical <- run_association_grid(
  data = analysis_t1_metag_t1_clinical,
  clinical_traits = clini_var_t1,
  clinical_period = "T1",
  microbiome_period = "T1",
  microbiota_features = target_metag
)

message("Running T1 MetaG genus/species -> T2 OGTT traits...")
results_t1_metag_t2_clinical <- run_association_grid(
  data = analysis_t1_metag_t2_clinical,
  clinical_traits = clini_var_t2,
  clinical_period = "T2",
  microbiome_period = "T1",
  microbiota_features = target_metag
)

results_t1_micro_t1_clinical <- bind_rows(
  results_t1_16sits_t1_clinical,
  results_t1_metag_t1_clinical
)
results_t1_micro_t2_clinical <- bind_rows(
  results_t1_16sits_t2_clinical,
  results_t1_metag_t2_clinical
)

association_results <- bind_rows(
  results_t1_micro_t1_clinical,
  results_t1_micro_t2_clinical
)

results_16s <- association_results %>%
  filter(microbiome_type == "16S_g")

results_its <- association_results %>%
  filter(microbiome_type == "ITS_g")

results_metag_g <- association_results %>%
  filter(microbiome_type == "MetaG_g")

results_metag_s <- association_results %>%
  filter(microbiome_type == "MetaG_s")

analysis_sample_summary <- bind_rows(
  tibble(
    microbiome_platform = "16S/ITS",
    clinical_period = "T1",
    microbiome_period = "T1",
    n_joined = n_distinct(analysis_t1_micro_t1_clinical$id)
  ),
  tibble(
    microbiome_platform = "16S/ITS",
    clinical_period = "T2",
    microbiome_period = "T1",
    n_joined = n_distinct(analysis_t1_micro_t2_clinical$id)
  ),
  tibble(
    microbiome_platform = "MetaG",
    clinical_period = "T1",
    microbiome_period = "T1",
    n_joined = n_distinct(analysis_t1_metag_t1_clinical$id)
  ),
  tibble(
    microbiome_platform = "MetaG",
    clinical_period = "T2",
    microbiome_period = "T1",
    n_joined = n_distinct(analysis_t1_metag_t2_clinical$id)
  )
)

feature_summary <- tibble(
  microbiome_type = c("16S_g", "ITS_g", "MetaG_g", "MetaG_s"),
  n_features = c(
    length(owo_16s), length(owo_its), length(owo_mgs_g), length(owo_mgs_s)
  ),
  feature_source = c(
    "owo_rep_id$rep_genus",
    "sig_its_sets$sigid_its_owo",
    "owo_rep_id$rep_genus_mgs",
    "owo_rep_id$rep_species"
  ),
  n_available = c(
    sum(owo_16s %in% names(bact_fung_clrz)),
    sum(owo_its %in% names(bact_fung_clrz)),
    sum(owo_mgs_g %in% names(mgs_clrz)),
    sum(owo_mgs_s %in% names(mgs_clrz))
  )
)

status_summary <- association_results %>%
  count(
    clinical_period, microbiome_period, microbiome_type, status,
    name = "n_models"
  )

# Output 
dir.create("result/PART3", recursive = TRUE, showWarnings = FALSE)

save(
  association_results,
  results_16s,
  results_its,
  results_metag_g,
  results_metag_s,
  results_t1_micro_t1_clinical,
  results_t1_micro_t2_clinical,
  results_t1_16sits_t1_clinical,
  results_t1_16sits_t2_clinical,
  results_t1_metag_t1_clinical,
  results_t1_metag_t2_clinical,
  owo_16s,
  owo_its,
  owo_mgs_g,
  owo_mgs_s,
  target_mgs_g,
  target_mgs_s,
  missing_metag_features,
  feature_availability,
  microbiota_annotation,
  feature_summary,
  analysis_sample_summary,
  status_summary,
  file = "result/PART3/3.1_association_micro_clinical_traits_thsbc.RData"
)
message("Clinical-trait association analyses completed.")

# Clinical relevance in HBC----------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(writexl)
})

# Input 
load("data/02_hbc_cleaned_data/gut_microbe_fung_393_158.RData")
load("data/02_hbc_cleaned_data/gut_microbe_MetaG_393_158.RData")
load("data/02_hbc_cleaned_data/total_pop_hbc_393.Rdata")
load("data/02_hbc_cleaned_data/covariates_393.RData")
load("data/02_hbc_cleaned_data/clinical_traits_t1t2_393.RData")
load("result/PART3/thsbc_hbc_rep_genus_species_owo_sig_sets.rdata")
load("result/PART3/1_diff_16s_its_sig_sets_thsbc.rdata")

owo_its <- info_fung_its_hbc %>%
  filter(ID_thsbc %in% str_remove(unique(sig_its_sets$sigid_its_owo), "fung_")) %>%
  transmute(feature = str_c("fung_", ID)) %>%
  pull(feature) %>%
  unique()

# Use shared MetaG feature IDs after verifying cross-cohort taxonomy concordance
owo_mgs_g <- unique(owo_rep_id$rep_genus_mgs)
owo_mgs_s <- unique(owo_rep_id$rep_species)
target_mgs_g <- owo_mgs_g
target_mgs_s <- owo_mgs_s
target_metag <- unique(c(target_mgs_g, target_mgs_s))
target_microbiota <- unique(c(owo_its, target_metag))
missing_metag_features <- setdiff(c(owo_mgs_g, owo_mgs_s), names(mgs_clrz_hbc))

feature_availability <- bind_rows(
  tibble(microbiome_type = "MetaG_g", microbe = owo_mgs_g),
  tibble(microbiome_type = "MetaG_s", microbe = owo_mgs_s)
) %>%
  mutate(
    available_in_saved_clr = microbe %in% names(mgs_clrz_hbc),
    exclusion_reason = if_else(
      available_in_saved_clr,
      NA_character_,
      "not_present_in_saved_mgs_clrz_hbc"
    )
  ) %>%
  left_join(
    info_mgs %>%
      transmute(
        microbe = str_c("mgs_", ID),
        thsbc_microbe = microbe,
        microbe_name = name
      ),
    by = "microbe"
  )

stopifnot(
  length(owo_its) == 6L,
  length(owo_mgs_g) == 23L,
  length(owo_mgs_s) == 37L,
  !anyDuplicated(owo_mgs_g),
  !anyDuplicated(owo_mgs_s),
  all(owo_its %in% names(fung_clrz_hbc)),
  all(target_mgs_g %in% names(mgs_clrz_hbc)),
  all(target_mgs_s %in% names(mgs_clrz_hbc))
)

if (length(missing_metag_features) > 0L) {
  message(
    "Excluded from clinical-trait models because the feature is absent from ",
    "the saved HBC CLR matrix: ",
    paste(missing_metag_features, collapse = ", ")
  )
}

microbiota_annotation <- bind_rows(
  info_fung_its_hbc %>%
    mutate(microbe = str_c("fung_", ID)) %>%
    filter(microbe %in% owo_its) %>%
    transmute(
      microbiome_type = "ITS_g",
      microbe,
      thsbc_microbe = str_c("fung_", ID_thsbc),
      microbe_name = name
    ),
  info_mgs %>%
    mutate(microbe = str_c("mgs_", ID)) %>%
    filter(microbe %in% target_metag) %>%
    transmute(
      microbiome_type = if_else(microbe %in% target_mgs_g, "MetaG_g", "MetaG_s"),
      microbe,
      thsbc_microbe = microbe,
      microbe_name = name
    )
) %>%
  distinct(microbiome_type, microbe, .keep_all = TRUE)

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
make_its_period_data <- function(microbiome_period) {
  fung_clrz_hbc %>%
    select(id, period, wk, all_of(owo_its)) %>%
    inner_join(cov_tab, by = c("id", "period")) %>%
    filter(period == microbiome_period) %>%
    mutate(
      across(any_of(factor_covariates), as.factor),
      microbiome_period = microbiome_period
    ) %>%
    select(-period)
}

make_metag_period_data <- function(microbiome_period) {
  mgs_clrz_hbc %>%
    select(id, period, wk, all_of(target_metag)) %>%
    inner_join(cov_tab, by = c("id", "period")) %>%
    filter(period == microbiome_period) %>%
    mutate(
      across(any_of(factor_covariates), as.factor),
      microbiome_period = microbiome_period
    ) %>%
    select(-period)
}

its_t1 <- make_its_period_data("T1")
its_t2 <- make_its_period_data("T2")
metag_t1 <- make_metag_period_data("T1")
metag_t2 <- make_metag_period_data("T2")

analysis_t1_its_t1_clinical <- its_t1 %>%
  inner_join(clini_traits_T1 %>% select(-period), by = "id")
analysis_t1_its_t2_clinical <- its_t1 %>%
  inner_join(clini_traits_T2 %>% select(-period), by = "id")

analysis_t1_metag_t1_clinical <- metag_t1 %>%
  inner_join(clini_traits_T1 %>% select(-period), by = "id")
analysis_t1_metag_t2_clinical <- metag_t1 %>%
  inner_join(clini_traits_T2 %>% select(-period), by = "id")

# Ordinary linear regression 
fit_one_association <- function(data, clinical_trait, microbe,
                                clinical_period, microbiome_period) {
  microbiome_type <- case_when(
    microbe %in% owo_its ~ "ITS_g",
    microbe %in% owo_mgs_g ~ "MetaG_g",
    microbe %in% owo_mgs_s ~ "MetaG_s",
    TRUE ~ NA_character_
  )
  if (is.na(microbiome_type)) {
    stop("Unknown microbiome feature: ", microbe)
  }

  required_columns <- unique(c("id", clinical_trait, microbe, covariates))
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
      .microbe_exposure = .data[[microbe]]
    )

  empty_result <- function(status_text) {
    tibble(
      clinical_trait = clinical_trait,
      microbiome_type = microbiome_type,
      microbe = microbe,
      clinical_period = clinical_period,
      microbiome_period = microbiome_period,
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
      n_distinct(model_data[[microbe]]) < 2L) {
    return(empty_result("no_variation"))
  }

  estimable_covariates <- covariates[
    map_lgl(
      covariates,
      ~ n_distinct(model_data[[.x]], na.rm = TRUE) > 1L
    )
  ]
  model_formula <- reformulate(
    termlabels = c(".microbe_exposure", estimable_covariates),
    response = ".clinical_outcome"
  )

  tryCatch({
    fit <- lm(stats::as.formula(model_formula), data = model_data)
    coefficient_table <- summary(fit)$coefficients
    coefficient_row <- match(".microbe_exposure", rownames(coefficient_table))

    if (is.na(coefficient_row)) {
      return(empty_result("microbe_coefficient_not_estimable"))
    }

    estimate <- coefficient_table[coefficient_row, "Estimate"]
    standard_error <- coefficient_table[coefficient_row, "Std. Error"]
    critical_value <- qt(0.975, df = fit$df.residual)

    tibble(
      clinical_trait = clinical_trait,
      microbiome_type = microbiome_type,
      microbe = microbe,
      clinical_period = clinical_period,
      microbiome_period = microbiome_period,
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
                                 microbiome_period, microbiota_features) {
  tidyr::crossing(
    clinical_trait = clinical_traits,
    microbe = microbiota_features
  ) %>%
    mutate(
      result = map2(
        clinical_trait,
        microbe,
        ~ fit_one_association(
          data = data,
          clinical_trait = .x,
          microbe = .y,
          clinical_period = clinical_period,
          microbiome_period = microbiome_period
        )
      )
    ) %>%
    select(result) %>%
    unnest(result) %>%
    group_by(
      clinical_period, microbiome_period, microbiome_type, clinical_trait
    ) %>%
    mutate(
      p_fdr = p.adjust(p_value, method = "fdr"),
      p_bonferroni = p.adjust(p_value, method = "bonferroni")
    ) %>%
    ungroup() %>%
    left_join(
      microbiota_annotation,
      by = c("microbiome_type", "microbe")
    ) %>%
    relocate(thsbc_microbe, microbe_name, .after = microbe) %>%
    arrange(
      clinical_period, microbiome_period, microbiome_type,
      clinical_trait, p_value
    )
}

# Analyses 
message("Running HBC T1 ITS -> T1 clinical traits...")
results_t1_its_t1_clinical <- run_association_grid(
  analysis_t1_its_t1_clinical, clini_var_t1_hbc, "T1", "T1", owo_its
)
message("Running HBC T1 ITS -> T2 OGTT traits...")
results_t1_its_t2_clinical <- run_association_grid(
  analysis_t1_its_t2_clinical, clini_var_t2_hbc, "T2", "T1", owo_its
)

message("Running HBC T1 MetaG genus/species -> T1 clinical traits...")
results_t1_metag_t1_clinical <- run_association_grid(
  analysis_t1_metag_t1_clinical, clini_var_t1_hbc, "T1", "T1", target_metag
)
message("Running HBC T1 MetaG genus/species -> T2 OGTT traits...")
results_t1_metag_t2_clinical <- run_association_grid(
  analysis_t1_metag_t2_clinical, clini_var_t2_hbc, "T2", "T1", target_metag
)

results_t1_micro_t1_clinical <- bind_rows(
  results_t1_its_t1_clinical,
  results_t1_metag_t1_clinical
)
results_t1_micro_t2_clinical <- bind_rows(
  results_t1_its_t2_clinical,
  results_t1_metag_t2_clinical
)

association_results <- bind_rows(
  results_t1_micro_t1_clinical,
  results_t1_micro_t2_clinical
)

results_its <- association_results %>%
  filter(microbiome_type == "ITS_g")
results_metag_g <- association_results %>%
  filter(microbiome_type == "MetaG_g")
results_metag_s <- association_results %>%
  filter(microbiome_type == "MetaG_s")

feature_summary <- tibble(
  microbiome_type = c("ITS_g", "MetaG_g", "MetaG_s"),
  n_features = c(length(owo_its), length(owo_mgs_g), length(owo_mgs_s)),
  feature_source = c(
    "HBC features mapped from sig_its_sets$sigid_its_owo",
    "Shared MetaG IDs from owo_rep_id$rep_genus_mgs",
    "Shared MetaG IDs from owo_rep_id$rep_species"
  ),
  n_available = c(
    sum(owo_its %in% names(fung_clrz_hbc)),
    sum(owo_mgs_g %in% names(mgs_clrz_hbc)),
    sum(owo_mgs_s %in% names(mgs_clrz_hbc))
  ),
  weighting = "None (ordinary linear regression)"
)

analysis_sample_summary <- bind_rows(
  tibble(
    microbiome_platform = "ITS",
    clinical_period = "T1",
    microbiome_period = "T1",
    n_joined = n_distinct(analysis_t1_its_t1_clinical$id)
  ),
  tibble(
    microbiome_platform = "ITS",
    clinical_period = "T2",
    microbiome_period = "T1",
    n_joined = n_distinct(analysis_t1_its_t2_clinical$id)
  ),
  tibble(
    microbiome_platform = "MetaG",
    clinical_period = "T1",
    microbiome_period = "T1",
    n_joined = n_distinct(analysis_t1_metag_t1_clinical$id)
  ),
  tibble(
    microbiome_platform = "MetaG",
    clinical_period = "T2",
    microbiome_period = "T1",
    n_joined = n_distinct(analysis_t1_metag_t2_clinical$id)
  )
)

status_summary <- association_results %>%
  count(
    clinical_period, microbiome_period, microbiome_type, status,
    name = "n_models"
  )

# Output 
dir.create("result/PART3", recursive = TRUE, showWarnings = FALSE)

save(
  association_results,
  results_its,
  results_metag_g,
  results_metag_s,
  results_t1_micro_t1_clinical,
  results_t1_micro_t2_clinical,
  results_t1_its_t1_clinical,
  results_t1_its_t2_clinical,
  results_t1_metag_t1_clinical,
  results_t1_metag_t2_clinical,
  owo_its,
  owo_mgs_g,
  owo_mgs_s,
  target_mgs_g,
  target_mgs_s,
  missing_metag_features,
  feature_availability,
  microbiota_annotation,
  feature_summary,
  analysis_sample_summary,
  status_summary,
  file = "result/PART3/3.2_association_micro_clinical_traits_hbc.RData"
)

print(feature_summary)
print(analysis_sample_summary)
print(status_summary)
message("HBC microbiota-clinical trait association analyses completed.")
