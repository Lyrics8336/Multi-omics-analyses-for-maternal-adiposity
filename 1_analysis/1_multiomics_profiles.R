##############################################################################################
# Analyses for Figure 2
# Profiling maternal serum metabolomic and gut microbial composition and temporal variation
# PCA and PERMANOVA for Figure 2b are implemented in plot_fig2.R
###############################################################################################

# Contributions of maternal factors to variation in metabolome and microbiome ----------------------
rm(list = ls())

library(tidyverse)
library(purrr)
library(data.table)
library(readxl)
library(writexl)
library(vegan)
library(pairwiseAdonis)
source('clean_code/01_help_func.R')

load('data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData')
load('data/01_thsbc_cleaned_data/gut_microbe_MetaG_1250_473.RData')
load('data/01_thsbc_cleaned_data/metabolites_1292_991.RData')
load('data/01_thsbc_cleaned_data/total_pop_5212.RData')
load('data/01_thsbc_cleaned_data/covariates_5212.RData')

DF_metab <- total_pop %>% inner_join(metab_tab_st) %>% inner_join(cov_tab)
DF_16SITS <- total_pop %>% inner_join(bact_fung_tab) %>% inner_join(cov_tab)
DF_metaG <- total_pop %>% inner_join(mgs_tab) %>% inner_join(cov_tab)
table(DF_metab$period); table(DF_16SITS$period); table(DF_metaG$period)

phenotype <- c(
  "pbmi", "wc_t1", "whr_t1", "bfr_t1",
  "pbmi_cut3", "wc_t1_cut", "whr_t1_cut", "bfr_t1_cut",
  "age", "parity", "smk", "drk", "edu", 
  "fd_grain", "fd_meat", "fd_vege", "fd_fruit",
  "fd_egg", "fd_milkpro", "physic_activity",
  "medical_insurance", "income", 'wk', 'antibiotic'
)

phenotype_group <- list(
  anthropometric = c("pbmi", "wc_t1", "whr_t1", "bfr_t1"),
  demography     = c("age", "parity"),
  dietary        = c("fd_grain", "fd_meat", "fd_vege", "fd_fruit", "fd_egg", "fd_milkpro"),
  behavioral     = c("smk", "drk", "physic_activity"),
  social_determ  = c("edu", "medical_insurance", "income"),
  all            = c("pbmi", "wc_t1", "whr_t1", "bfr_t1", 
                     "age", "parity", "smk", "drk", "edu", 
                     "fd_grain", "fd_meat", "fd_vege", "fd_fruit",
                     "fd_egg", "fd_milkpro", "physic_activity",
                     "medical_insurance", "income")
)

timepoints <- c("T1", "T2", "T3", "TA")

omics_config <- list(
  metabolites = list(
    data        = DF_metab,
    features    = metab_uni,
    label_type  = "Metabolites",
    dist_method = "euclidean",   # Euclidean distance for continuous metabolomic concentrations
    apply_tss   = FALSE          # No TSS scaling for metabolomics
  ),
  `16S_g` = list(
    data        = DF_16SITS,
    features    = bact_uni,
    label_type  = "16S_g",
    dist_method = "bray",        # Bray-Curtis dissimilarity for compositional microbial data
    apply_tss   = TRUE           # Enable Total Sum Scaling (relative abundance normalization)
  ),
  ITS_g = list(
    data        = DF_16SITS,
    features    = fung_uni,
    label_type  = "ITS_g",
    dist_method = "bray",
    apply_tss   = TRUE
  ),
  MetaG_g = list(
    data        = DF_metaG,
    features    = mgs_g_uni,
    label_type  = "MetaG_g",
    dist_method = "bray",
    apply_tss   = TRUE           # Shotgun metagenomics: Genus level
  ),
  MetaG_s = list(
    data        = DF_metaG,
    features    = mgs_s_uni,
    label_type  = "MetaG_s",
    dist_method = "bray",
    apply_tss   = TRUE           # Shotgun metagenomics: Species level
  )
)

message("Starting batch analysis for single clinical phenotypes across all omics layers...")

grid_single <- expand_grid(
  omics_name  = names(omics_config),
  target_cols = phenotype,
  period_term = timepoints
)

Permanova_Single_Phenotypes <- purrr::pmap_dfr(
  .l = grid_single,
  .f = function(omics_name, target_cols, period_term) {
    cfg <- omics_config[[omics_name]]
    func_permanova_integrated(
      data         = cfg$data,
      feature_cols = cfg$features,
      target_cols  = target_cols,
      period_term  = period_term,
      label_type   = cfg$label_type,
      label_var    = target_cols,
      dist_method  = cfg$dist_method,
      apply_tss    = cfg$apply_tss,
      by_mode      = NULL
    )
  },
  .progress = TRUE
) %>% 
  arrange(type, period)

message("Single phenotype analysis complete. Starting phenotypic cluster analysis...")

grid_group <- expand_grid(
  omics_name  = names(omics_config),
  group_name  = names(phenotype_group),
  period_term = timepoints
)

Permanova_Phenotype_Clusters <- purrr::pmap_dfr(
  .l = grid_group,
  .f = function(omics_name, group_name, period_term) {
    cfg <- omics_config[[omics_name]]
    func_permanova_integrated(
      data         = cfg$data,
      feature_cols = cfg$features,
      target_cols  = phenotype_group[[group_name]], 
      period_term  = period_term,
      label_type   = cfg$label_type,
      label_var    = group_name,
      dist_method  = cfg$dist_method,
      apply_tss    = cfg$apply_tss,
      by_mode      = NULL,
      compute_adjusted_r2 = TRUE
    )
  },
  .progress = TRUE
) %>% 
  arrange(type, period)

save(Permanova_Single_Phenotypes, Permanova_Phenotype_Clusters, 
  file = "result/PART1/1-1 Permanova_Single_Phenotypes_and_Phenotype_Clusters_avaliable_pop.rdata")

# Alpha diversity, F/B log-ratio, beta diversity of the THSBC and HBC gut microbiome ----------------------
rm(list = ls())

library(dplyr)
library(vegan)
library(tidyr)
library(purrr)
library(RNOmni)
library(readxl)

# Functions 
# Calculate alpha diversity using all taxa at the specified taxonomic level.
calculate_alpha_diversity <- function(data, taxon_prefix, data_type = c("count", "relative_abundance")) {
  data_type <- match.arg(data_type)
  key_vars <- c("id", "period")
  taxon_pattern <- paste0("^", taxon_prefix, "\\d+$")
  taxon_vars <- grep(taxon_pattern, names(data), value = TRUE)

  missing_keys <- setdiff(key_vars, names(data))
  if (length(missing_keys) > 0L) {
    stop("Missing key column(s): ", paste(missing_keys, collapse = ", "))
  }
  if (anyDuplicated(data[key_vars]) > 0L) {
    stop("Each id-period combination must uniquely identify one sample.")
  }
  if (length(taxon_vars) == 0L) {
    stop("No taxon columns matched: ", taxon_pattern)
  }

  taxon <- data[taxon_vars]
  non_numeric <- taxon_vars[!vapply(taxon, is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop("Non-numeric taxon column(s): ", paste(non_numeric, collapse = ", "))
  }

  taxon <- as.matrix(taxon)
  if (anyNA(taxon)) {
    stop("Taxon abundance data contain missing values.")
  }
  if (any(!is.finite(taxon))) {
    stop("Taxon abundance data contain non-finite values.")
  }
  if (any(taxon < 0)) {
    stop("Taxon abundance data contain negative values.")
  }
  sample_total <- rowSums(taxon)
  if (any(sample_total == 0)) {
    stop("At least one sample has zero total abundance.")
  }
  if (
    data_type == "count" &&
      any(abs(taxon - round(taxon)) > sqrt(.Machine$double.eps))
  ) {
    stop("Count-based alpha diversity requires integer taxon abundances.")
  }

  # Shannon and Simpson indices depend on within-sample proportions. Explicit
  # normalization makes count data, proportions (sum = 1), and percentages
  # (sum = 100) follow exactly the same calculation.
  relative_abundance <- vegan::decostand(
    taxon,
    method = "total",
    MARGIN = 1
  )
  richness <- rowSums(taxon > 0)
  shannon <- vegan::diversity(relative_abundance,index = "shannon",base = exp(1))
  simpson <- vegan::diversity(relative_abundance, index = "simpson")
  invsimpson <- vegan::diversity(relative_abundance,index = "invsimpson")

  if (data_type == "count") {
    chao1 <- unname(vegan::estimateR(taxon)["S.chao1", ])
    pielou <- ifelse(richness > 1, shannon / log(richness), NA_real_)

    result <- data.frame(
      data[key_vars],
      richness, chao1, shannon, simpson, invsimpson, pielou, check.names = FALSE
    )
  } else {
    result <- data.frame(
      data[key_vars],
      richness, shannon, simpson, invsimpson, check.names = FALSE
    )
  }

  rownames(result) <- NULL
  result
}

# Calculate Bray-Curtis PCoA from all unfiltered taxa. All periods are
# included in one ordination so that their coordinates are comparable.
calculate_beta_pcoa <- function(data, taxon_prefix) {
  key_vars <- c("id", "period")
  taxon_pattern <- paste0("^", taxon_prefix, "\\d+$")
  taxon_vars <- grep(taxon_pattern, names(data), value = TRUE)

  missing_keys <- setdiff(key_vars, names(data))
  if (length(missing_keys) > 0L) {
    stop("Missing key column(s): ", paste(missing_keys, collapse = ", "))
  }
  if (anyDuplicated(data[key_vars]) > 0L) {
    stop("Each id-period combination must uniquely identify one sample.")
  }
  if (length(taxon_vars) == 0L) {
    stop("No taxon columns matched: ", taxon_pattern)
  }

  taxon <- data[taxon_vars]
  non_numeric <- taxon_vars[!vapply(taxon, is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop("Non-numeric taxon column(s): ", paste(non_numeric, collapse = ", "))
  }

  taxon <- as.matrix(taxon)
  if (anyNA(taxon)) {
    stop("Taxon abundance data contain missing values.")
  }
  if (any(!is.finite(taxon))) {
    stop("Taxon abundance data contain non-finite values.")
  }
  if (any(taxon < 0)) {
    stop("Taxon abundance data contain negative values.")
  }

  sample_total <- rowSums(taxon)
  if (any(sample_total == 0)) {
    stop("At least one sample has zero total abundance.")
  }

  relative_abundance <- vegan::decostand(
    taxon,
    method = "total",
    MARGIN = 1
  )
  rownames(relative_abundance) <- paste(data$id, data$period, sep = "_")
  bray_curtis <- vegan::vegdist(relative_abundance, method = "bray")

  pcoa <- stats::cmdscale(
    bray_curtis,
    k = 2,
    eig = TRUE
  )
  variance_explained <- 100 * pcoa$eig[1:2] / sum(pcoa$eig)

  scores <- data.frame(
    data[key_vars],
    PCoA1 = unname(pcoa$points[, 1]),
    PCoA2 = unname(pcoa$points[, 2]),
    check.names = FALSE
  )

  list(
    scores = scores,
    variance_explained = unname(variance_explained),
    bray_curtis = bray_curtis
  )
}

# Prepare and compare alpha-diversity and F/B-ratio outcomes 
group_variables <- c( "pbmi_cut3", "wc_t1_cut", "whr_t1_cut", "bfr_t1_cut", "pbmi_cut4" )

rank_inverse_normal <- function(x) {
  non_missing <- is.finite(x)
  transformed <- rep(NA_real_, length(x))
  if (any(non_missing)) {
    transformed[non_missing] <- RNOmni::RankNorm(
      x[non_missing],
      k = 3 / 8
    )
  }
  transformed
}

get_requested_comparisons <- function(group_variable, observed_groups) {
  requested <- switch(
    group_variable,
    pbmi_cut3 = list(
      c("Underweight", "Normal"),
      c("Underweight", "OWO"),
      c("Normal", "OWO")
    ),
    wc_t1_cut = list(c("Normal", "Abdominal_obesity")),
    whr_t1_cut = list(c("Normal", "Abdominal_obesity")),
    bfr_t1_cut = list(c("Normal", "Obesity")),
    pbmi_cut4 = list(
      c("Normal", "Overweight"),
      c("Normal", "Obese"),
      c("Overweight", "Obese")
    ),
    list()
  )

  keep(requested, ~ all(.x %in% observed_groups))
}

run_one_group_comparison <- function(data, group_variable, strata) {
  if (!group_variable %in% names(data) || all(is.na(data[[group_variable]]))) {
    return(list(overall = tibble(), pairwise = tibble(), trend = tibble()))
  }

  nested_data <- data %>%
    filter( is.finite(value), !is.na(.data[[group_variable]]) ) %>%
    mutate(group = droplevels(factor(.data[[group_variable]]))) %>%
    group_by(across(all_of(strata))) %>%
    nest()

  if (group_variable %in% c("pbmi_cut3", "pbmi_cut4")) {
    overall <- nested_data %>%
      mutate(result = map(data, function(df) {
        if (n_distinct(df$group) < 2L) {
          return(tibble(
            n = nrow(df), n_group = n_distinct(df$group),
            statistic = NA_real_, p_value = NA_real_
          ))
        }
        fit <- kruskal.test(value ~ group, data = df)
        tibble(
          n = nrow(df), n_group = n_distinct(df$group),
          statistic = unname(fit$statistic), p_value = fit$p.value
        )
      })) %>%
      select(-data) %>%
      unnest(result) %>%
      ungroup() %>%
      mutate(group_variable = group_variable, .before = 1)
  } else {
    overall <- tibble()
  }

  pairwise <- nested_data %>%
    mutate(result = map(data, function(df) {
      observed_groups <- unique(as.character(df$group))
      comparisons <- get_requested_comparisons(
        group_variable,
        observed_groups
      )
      if (length(comparisons) == 0L) {
        return(tibble())
      }

      map_dfr(comparisons, function(groups) {
        group_1_values <- df$value[df$group == groups[1]]
        group_2_values <- df$value[df$group == groups[2]]
        fit <- wilcox.test(
          group_1_values,
          group_2_values,
          exact = FALSE
        )
        tibble(
          group_1 = groups[1],
          group_2 = groups[2],
          comparison = paste(groups, collapse = " vs. "),
          n_1 = length(group_1_values),
          median_1 = median(group_1_values),
          n_2 = length(group_2_values),
          median_2 = median(group_2_values),
          p_value = fit$p.value
        )
      }) %>%
        mutate(
          p_adjusted = p.adjust(p_value, method = "fdr"),
          significance = case_when(
            p_adjusted < 0.001 ~ "***", p_adjusted < 0.01 ~ "**",
            p_adjusted < 0.05 ~ "*", TRUE ~ "ns"
          )
        )
    })) %>%
    select(-data) %>%
    unnest(result) %>%
    ungroup() %>%
    mutate(group_variable = group_variable, .before = 1)

  if (group_variable == "pbmi_cut3") {
    trend <- nested_data %>%
      mutate(result = map(data, function(df) {
        trend_data <- df %>%
          mutate(
            group_score = match(
              as.character(group),
              c("Underweight", "Normal", "OWO")
            )
          )
        fit <- lm(value ~ group_score, data = trend_data)
        coefficient <- summary(fit)$coefficients["group_score", ]
        tibble(
          n = nrow(trend_data),
          beta_per_category = unname(coefficient["Estimate"]),
          standard_error = unname(coefficient["Std. Error"]),
          statistic = unname(coefficient["t value"]),
          p_trend = unname(coefficient["Pr(>|t|)"])
        )
      })) %>%
      select(-data) %>%
      unnest(result) %>%
      ungroup() %>%
      mutate(group_variable = group_variable, .before = 1)
  } else {
    trend <- tibble()
  }

  list(overall = overall, pairwise = pairwise, trend = trend)
}

run_all_group_comparisons <- function(data) {
  strata <- c("cohort", "period", "data_type", "omics", "metric")
  results <- map(
    group_variables,
    ~ run_one_group_comparison(data, .x, strata)
  )

  list(
    overall = map_dfr(results, "overall"),
    pairwise = map_dfr(results, "pairwise"),
    trend = map_dfr(results, "trend")
  )
}

prepare_microbiome_comparison_data <- function( alpha_objects, fb_objects, population, cohort_name) {
  available_group_variables <- intersect( group_variables, names(population) )
  group_metadata <- population %>%
    select(id, all_of(available_group_variables))
  missing_group_variables <- setdiff( group_variables, available_group_variables )
  for (variable in missing_group_variables) {
    group_metadata[[variable]] <- NA_character_
  }

  alpha_data <- imap_dfr(alpha_objects, function(data, omics_name) {
    data %>%
      inner_join(group_metadata, by = "id") %>%
      select(
        id, period, all_of(group_variables),
        richness, shannon
      ) %>%
      pivot_longer(
        cols = c(richness, shannon), names_to = "metric", values_to = "value"
      ) %>%
      mutate(
        cohort = cohort_name, data_type = "alpha_diversity", omics = omics_name, .before = 1
      )
  })

  fb_data <- imap_dfr(fb_objects, function(data, omics_name) {
    data %>%
      inner_join(group_metadata, by = "id") %>%
      transmute(
        cohort = cohort_name,
        data_type = "fb_ratio",
        omics = omics_name,
        id, period,
        across(all_of(group_variables)),
        metric = "fb_ratio",
        value = fb_ratio
      )
  })

  raw_data <- bind_rows(alpha_data, fb_data) %>%
    mutate(period = factor(period, levels = c("T1", "T2", "T3")))

  rank_data <- raw_data %>%
    group_by(cohort, period, data_type, omics, metric) %>%
    mutate(
      value_raw = value,
      value = rank_inverse_normal(value)
    ) %>%
    ungroup()

  list(
    raw_data = raw_data,
    rank_data = rank_data,
    raw_results = run_all_group_comparisons(raw_data),
    rank_results = run_all_group_comparisons(rank_data),
    group_variable_availability = tibble(
      cohort = cohort_name,
      group_variable = group_variables,
      available = group_variables %in% available_group_variables
    )
  )
}


# THSBC cohort 
load("data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData")
load("data/01_thsbc_cleaned_data/gut_microbe_MetaG_1250_473.RData")
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")
load("data/01_thsbc_cleaned_data/wt_omics.rdata")

DF_16SITS <- total_pop %>%
  inner_join(bact_fung_tab, by = "id") 

DF_MetaG <- total_pop %>%
  inner_join(mgs_tab, by = "id") 

## alpha diversity
alpha_16s_thsbc <- calculate_alpha_diversity(data = DF_16SITS, taxon_prefix = "bact_g", data_type = "count")
alpha_its_thsbc <- calculate_alpha_diversity(data = DF_16SITS, taxon_prefix = "fung_g", data_type = "count")
alpha_mgs_g_thsbc <- calculate_alpha_diversity(data = DF_MetaG, taxon_prefix = "mgs_g", data_type = "relative_abundance")
alpha_mgs_s_thsbc <- calculate_alpha_diversity(data = DF_MetaG, taxon_prefix = "mgs_s", data_type = "relative_abundance")
save(alpha_16s_thsbc, alpha_its_thsbc, alpha_mgs_g_thsbc, alpha_mgs_s_thsbc,
  file = "result/PART1/3-1_alpha_diversity_thsbc.RData")

## beta diversity
beta_16s_fit <- calculate_beta_pcoa(DF_16SITS, "bact_g")
beta_its_fit <- calculate_beta_pcoa(DF_16SITS, "fung_g")
beta_mgs_g_fit <- calculate_beta_pcoa(DF_MetaG, "mgs_g")
beta_mgs_s_fit <- calculate_beta_pcoa(DF_MetaG, "mgs_s")

beta_16s_thsbc <- beta_16s_fit$scores
beta_its_thsbc <- beta_its_fit$scores
beta_mgs_g_thsbc <- beta_mgs_g_fit$scores
beta_mgs_s_thsbc <- beta_mgs_s_fit$scores

bc_16s_thsbc <- beta_16s_fit$bray_curtis
bc_its_thsbc <- beta_its_fit$bray_curtis
bc_mgs_g_thsbc <- beta_mgs_g_fit$bray_curtis
bc_mgs_s_thsbc <- beta_mgs_s_fit$bray_curtis

beta_pcoa_variance_thsbc <- data.frame(
  data_type = c("16S_genus", "ITS_genus", "MetaG_genus", "MetaG_species"),
  PCoA1_percent = c(
    beta_16s_fit$variance_explained[1],
    beta_its_fit$variance_explained[1],
    beta_mgs_g_fit$variance_explained[1],
    beta_mgs_s_fit$variance_explained[1]
  ),
  PCoA2_percent = c(
    beta_16s_fit$variance_explained[2],
    beta_its_fit$variance_explained[2],
    beta_mgs_g_fit$variance_explained[2],
    beta_mgs_s_fit$variance_explained[2]
  )
)

save( beta_16s_thsbc, beta_its_thsbc, beta_mgs_g_thsbc, beta_mgs_s_thsbc,
  bc_16s_thsbc, bc_its_thsbc, bc_mgs_g_thsbc, bc_mgs_s_thsbc, beta_pcoa_variance_thsbc,
  file = "result/PART1/3-3_beta_diversity_thsbc.RData"
)

## F/B ratio
bact_tab_p <- read_xlsx('16S_Batch123.xlsx') %>% 
  distinct(id, period, .keep_all = T) %>% 
  mutate(period = str_replace(period, 'V', 'T')) %>% 
  filter(id %in% total_pop$id) %>% 
  select(id, period, bact_batch = batch, matches('^p\\d+'))

abund_bact_p <- bact_tab_p %>% 
  dplyr::select(matches('^p\\d+')) %>% 
  vegan::decostand('clr', 1, pseudocount = 1) %>% 
  mutate(fb_ratio = p1-p3)
summary(abund_bact_p$fb_ratio)

fb_16s_thsbc = bact_tab_p %>% select(id, period) %>% 
  mutate(fb_ratio = abund_bact_p$fb_ratio,
         fb_ratio = case_when(
           fb_ratio <  quantile(fb_ratio, 0.02) ~ NA,
           fb_ratio >  quantile(fb_ratio, 0.98) ~ NA,
           TRUE ~ fb_ratio
         )) 

mgs_tab_p = read_xlsx('mgs_taxonomy_abundance.xlsx') %>% 
  distinct(id, period, .keep_all = T) %>% 
  select(-sample_id) %>% 
  filter(id %in% total_pop$id) %>%
  select(id, period, mgs_batch = batch, matches('^p\\d+'))

abund_mgs_p <- mgs_tab_p %>% dplyr::select(matches('^p\\d+')) %>% 
  vegan::decostand('clr', 1, pseudocount = 1e-5) %>% 
  mutate(fb_ratio = p1-p4)
summary(abund_mgs_p$fb_ratio)

fb_mgs_thsbc = mgs_tab_p %>% select(id, period) %>% 
  mutate(fb_ratio = abund_mgs_p$fb_ratio,
         fb_ratio = case_when(
           fb_ratio <  quantile(fb_ratio, 0.02) ~ NA,
           fb_ratio >  quantile(fb_ratio, 0.98) ~ NA,
           TRUE ~ fb_ratio
         )) 

save( fb_16s_thsbc, fb_mgs_thsbc, file = "result/PART1/3-5_fb_ratio_thsbc.RData")

## group comparison
group_comparison_thsbc <- prepare_microbiome_comparison_data(
  alpha_objects = list(
    `16S_g` = alpha_16s_thsbc,
    `ITS_g` = alpha_its_thsbc,
    `MetaG_g` = alpha_mgs_g_thsbc,
    `MetaG_s` = alpha_mgs_s_thsbc
  ),
  fb_objects = list(
    `16S_g` = fb_16s_thsbc,
    `MetaG_g` = fb_mgs_thsbc
  ),
  population = total_pop,
  cohort_name = "THSBC"
)

micro_group_data_raw_thsbc <- group_comparison_thsbc$raw_data
micro_group_data_rank_thsbc <- group_comparison_thsbc$rank_data
micro_group_overall_raw_thsbc <- group_comparison_thsbc$raw_results$overall
micro_group_pairwise_raw_thsbc <- group_comparison_thsbc$raw_results$pairwise
micro_group_trend_raw_thsbc <- group_comparison_thsbc$raw_results$trend
micro_group_overall_rank_thsbc <- group_comparison_thsbc$rank_results$overall
micro_group_pairwise_rank_thsbc <- group_comparison_thsbc$rank_results$pairwise
micro_group_trend_rank_thsbc <- group_comparison_thsbc$rank_results$trend
group_variable_availability_thsbc <- group_comparison_thsbc$group_variable_availability

save(
  micro_group_data_raw_thsbc, micro_group_data_rank_thsbc,
  micro_group_overall_raw_thsbc, micro_group_pairwise_raw_thsbc,
  micro_group_trend_raw_thsbc,
  micro_group_overall_rank_thsbc, micro_group_pairwise_rank_thsbc,
  micro_group_trend_rank_thsbc,
  group_variable_availability_thsbc,
  file = "result/PART1/3-7_microbiome_group_comparisons_thsbc.RData"
)






# HBC cohort 

rm(
  bact_fung_tab,
  mgs_tab,
  total_pop,
  cov_tab,
  wt_16SITS,
  wt_MetaG,
  DF_16SITS,
  DF_MetaG
)

load("data/02_hbc_cleaned_data/gut_microbe_fung_393_158.RData")
load("data/02_hbc_cleaned_data/gut_microbe_MetaG_393_158.RData")
load("data/02_hbc_cleaned_data/total_pop_hbc_393.Rdata")
load("data/02_hbc_cleaned_data/covariates_393.RData")

DF_ITS <- total_pop_hbc %>%
  inner_join(fung_tab_hbc, by = "id")

DF_MetaG <- total_pop_hbc %>%
  inner_join(mgs_tab_hbc, by = "id") 

alpha_fung_hbc <- calculate_alpha_diversity(data = DF_ITS, taxon_prefix = "fung_g", data_type = "count")
alpha_mgs_g_hbc <- calculate_alpha_diversity(data = DF_MetaG, taxon_prefix = "mgs_g", data_type = "relative_abundance")
alpha_mgs_s_hbc <- calculate_alpha_diversity(data = DF_MetaG, taxon_prefix = "mgs_s", data_type = "relative_abundance")

save(alpha_fung_hbc, alpha_mgs_g_hbc, alpha_mgs_s_hbc,
  file = "result/PART1/3-2_alpha_diversity_hbc.RData"
)

beta_its_fit <- calculate_beta_pcoa(DF_ITS, "fung_g")
beta_mgs_g_fit <- calculate_beta_pcoa(DF_MetaG, "mgs_g")
beta_mgs_s_fit <- calculate_beta_pcoa(DF_MetaG, "mgs_s")

beta_its_hbc <- beta_its_fit$scores
beta_mgs_g_hbc <- beta_mgs_g_fit$scores
beta_mgs_s_hbc <- beta_mgs_s_fit$scores

bc_its_hbc <- beta_its_fit$bray_curtis
bc_mgs_g_hbc <- beta_mgs_g_fit$bray_curtis
bc_mgs_s_hbc <- beta_mgs_s_fit$bray_curtis

beta_pcoa_variance_hbc <- data.frame(
  data_type = c("ITS_genus", "MetaG_genus", "MetaG_species"),
  PCoA1_percent = c(
    beta_its_fit$variance_explained[1],
    beta_mgs_g_fit$variance_explained[1],
    beta_mgs_s_fit$variance_explained[1]
  ),
  PCoA2_percent = c(
    beta_its_fit$variance_explained[2],
    beta_mgs_g_fit$variance_explained[2],
    beta_mgs_s_fit$variance_explained[2]
  )
)

save(
  beta_its_hbc, beta_mgs_g_hbc, beta_mgs_s_hbc,
  bc_its_hbc, bc_mgs_g_hbc, bc_mgs_s_hbc,
  beta_pcoa_variance_hbc,
  file = "result/PART1/3-4_beta_diversity_hbc.RData"
)


# F/B ratio
wk_facel <- read_xlsx('HBC_fecal_wk_cleaned_1042.xlsx') %>% 
  select(-sample_id) %>% filter(period %in% c('T1', 'T3', 'T4')) %>% 
  mutate(period = case_match(period, 'T3'~'T2', 'T4'~'T3', .default = period)) 

mgs_tab_p_hbc <- read_xlsx('mgs_taxonomy_abundance_fecal_1042_HBC.xlsx') %>% 
  select(-c(sample_id, sample_type)) %>% 
  distinct(id, period, .keep_all = T) %>% 
  mutate(period = case_match(period, 'A'~'T1', 'C'~'T2', 'D'~'T3', .default = period)) %>% 
  filter(period %in% c('T1', 'T2', 'T3')) %>% 
  left_join(wk_facel) %>% 
  relocate(wk, .after = period) %>% 
  filter(id %in% total_pop_hbc$id) %>% 
  select(id, period, wk, matches('^p\\d+')) 

abund_mgs_p_hbc <- mgs_tab_p_hbc %>% dplyr::select(matches('^p\\d+')) %>% 
  vegan::decostand('clr', 1, pseudocount = 1e-5) %>% 
  mutate(fb_ratio = p1-p4)

fb_mgs_hbc = mgs_tab_p_hbc %>% select(id, period) %>% 
  mutate(fb_ratio = abund_mgs_p_hbc$fb_ratio,
         fb_ratio = case_when(
           fb_ratio <  quantile(fb_ratio, 0.02) ~ NA,
           fb_ratio >  quantile(fb_ratio, 0.98) ~ NA,
           TRUE ~ fb_ratio
         )) 

save( fb_mgs_hbc, file = "result/PART1/3-6_fb_ratio_hbc.RData")

## group comparison

group_comparison_hbc <- prepare_microbiome_comparison_data(
  alpha_objects = list(
    `ITS_g` = alpha_fung_hbc,
    `MetaG_g` = alpha_mgs_g_hbc,
    `MetaG_s` = alpha_mgs_s_hbc
  ),
  fb_objects = list(
    `MetaG_g` = fb_mgs_hbc
  ),
  population = total_pop_hbc,
  cohort_name = "HBC"
)

micro_group_data_raw_hbc <- group_comparison_hbc$raw_data
micro_group_data_rank_hbc <- group_comparison_hbc$rank_data
micro_group_overall_raw_hbc <- group_comparison_hbc$raw_results$overall
micro_group_pairwise_raw_hbc <- group_comparison_hbc$raw_results$pairwise
micro_group_trend_raw_hbc <- group_comparison_hbc$raw_results$trend
micro_group_overall_rank_hbc <- group_comparison_hbc$rank_results$overall
micro_group_pairwise_rank_hbc <- group_comparison_hbc$rank_results$pairwise
micro_group_trend_rank_hbc <- group_comparison_hbc$rank_results$trend
group_variable_availability_hbc <-
  group_comparison_hbc$group_variable_availability

save(
  micro_group_data_raw_hbc,
  micro_group_data_rank_hbc,
  micro_group_overall_raw_hbc,
  micro_group_pairwise_raw_hbc,
  micro_group_trend_raw_hbc,
  micro_group_overall_rank_hbc,
  micro_group_pairwise_rank_hbc,
  micro_group_trend_rank_hbc,
  group_variable_availability_hbc,
  file = "result/PART1/3-8_microbiome_group_comparisons_hbc.RData"
)

# Multi-omics distance analysis ----------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(writexl)
  library(vegan)
})

periods <- c("T1", "T2", "T3")
period_pairs <- list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3"))
n_perm <- as.integer(Sys.getenv("OMICS_N_PERM", "10000"))
max_between_per_anchor <- as.integer(Sys.getenv("OMICS_MAX_BETWEEN_PER_ANCHOR", "25"))
dir.create("result/PART1", recursive = TRUE, showWarnings = FALSE)
dir.create("figure/PART1", recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
}

as_id_period <- function(x) {
  x %>% mutate(id = as.character(id), period = as.character(period))
}

attach_metadata <- function(feature_data, features, population, covariates) {
  stop_if_missing(feature_data, c("id", "period", features), "Feature table")
  base_names <- intersect(c("id", "period", "wk", "batch", "bact_batch",
                            "fung_batch", "mgs_batch"), names(feature_data))
  base <- feature_data %>%
    as_id_period() %>%
    select(all_of(base_names), all_of(features)) %>%
    distinct(id, period, .keep_all = TRUE)

  pop <- population %>%
    mutate(id = as.character(id)) %>%
    select(id, any_of(c("pbmi_cut4", "pbmi_cut3", "pbmi_cut2",
                        "pbmi_subOverNor", "pbmi_subObseNor", "pbmi_subOverObse"))) %>%
    distinct(id, .keep_all = TRUE)
  out <- base %>% left_join(pop, by = "id")

  cov <- covariates %>% as_id_period()
  cov_names <- setdiff(names(cov), names(out))
  cov_names <- setdiff(cov_names, c("id", "period"))
  out <- out %>% left_join(cov %>% select(id, period, all_of(cov_names)),
                           by = c("id", "period"))

  if ("pbmi_cut3" %in% names(out) || "pbmi_cut4" %in% names(out)) {
    pbmi3_value <- if ("pbmi_cut3" %in% names(out)) as.character(out$pbmi_cut3) else rep(NA_character_, nrow(out))
    pbmi4_value <- if ("pbmi_cut4" %in% names(out)) as.character(out$pbmi_cut4) else rep(NA_character_, nrow(out))
    pbmi_value <- dplyr::coalesce(pbmi3_value, pbmi4_value)
    out <- out %>% mutate(
      # Retain all pBMI categories; restrict subgroup comparisons to Normal and OWO
      pBMI_group = case_when(
        pbmi_value %in% c("Normal", "normal") ~ "Normal",
        pbmi_value %in% c("OWO", "owo", "Overweight/obesity") ~ "OWO",
        TRUE ~ pbmi_value
      )
    )
  } else {
    out$pBMI_group <- "All"
  }
  out
}

cross_euclidean <- function(x, y) {
  d2 <- outer(rowSums(x * x), rowSums(y * y), "+") - 2 * tcrossprod(x, y)
  sqrt(pmax(d2, 0))
}

cross_distance <- function(x, y, metric) {
  if (metric == "euclidean") return(cross_euclidean(x, y))
  if (metric == "bray_curtis") {
    # Select the corresponding raw taxon features and convert each sample to
    # relative abundance before calculating Bray-Curtis distance.
    if (any(!is.finite(x)) || any(!is.finite(y)) ||
        any(x < 0, na.rm = TRUE) || any(y < 0, na.rm = TRUE) ||
        any(rowSums(x) <= 0) || any(rowSums(y) <= 0)) {
      stop("Bray-Curtis requires finite, non-negative raw abundances with positive sample totals.")
    }
    x <- vegan::decostand(x, method = "total", MARGIN = 1)
    y <- vegan::decostand(y, method = "total", MARGIN = 1)
    n_x <- nrow(x)
    d <- as.matrix(vegan::vegdist(rbind(x, y), method = "bray"))
    return(d[seq_len(n_x), n_x + seq_len(nrow(y)), drop = FALSE])
  }
  stop("Unknown distance metric: ", metric)
}

distance_components <- function(dm) {
  if (!is.matrix(dm) || nrow(dm) != ncol(dm)) {
    stop("Identity permutation requires a square cross-period distance matrix.")
  }
  within <- diag(dm)
  between <- dm[row(dm) != col(dm)]
  list(within = within[is.finite(within)],
       between = between[is.finite(between)])
}

mean_distance_difference <- function(dm) {
  x <- distance_components(dm)
  if (!length(x$within) || !length(x$between)) return(NA_real_)
  mean(x$within) - mean(x$between)
}

# Compute the permuted mean distance difference without copying the full matrix
mean_difference_for_match <- function(dm, perm, total_sum = NULL, total_n = NULL) {
  selected <- dm[cbind(seq_len(nrow(dm)), perm)]
  selected <- selected[is.finite(selected)]
  if (is.null(total_sum) || is.null(total_n)) {
    finite <- is.finite(dm)
    total_sum <- sum(dm[finite])
    total_n <- sum(finite)
  }
  n_between <- total_n - length(selected)
  if (!length(selected) || n_between < 1) return(NA_real_)
  mean(selected) - (total_sum - sum(selected)) / n_between
}

permute_identity <- function(n, strata = NULL) {
  if (is.null(strata)) return(sample.int(n))
  if (length(strata) != n || anyNA(strata)) {
    stop("Permutation strata must be complete and match the distance matrix.")
  }
  perm <- seq_len(n)
  for (index in split(seq_len(n), strata, drop = TRUE)) {
    perm[index] <- sample(index, length(index), replace = FALSE)
  }
  perm
}

# Adapted from the distance-comparison framework of Sinha et al. (2026).
# Participant identities at the second timepoint (matrix columns) are
# permuted, preserving the dependence among distances involving a participant.
permutation_identity_test <- function(dm, n_perm = 10000, seed = 1, strata = NULL) {
  observed_delta <- mean_distance_difference(dm)
  if (!is.finite(observed_delta)) {
    return(list(observed_delta = NA_real_,
                null_delta = rep(NA_real_, n_perm), p_value = NA_real_))
  }

  finite <- is.finite(dm)
  total_sum <- sum(dm[finite])
  total_n <- sum(finite)
  set.seed(seed)
  null_delta <- replicate(n_perm, {
    perm <- permute_identity(ncol(dm), strata)
    mean_difference_for_match(dm, perm, total_sum, total_n)
  })
  p_value <- (1 + sum(null_delta <= observed_delta, na.rm = TRUE)) /
    (1 + sum(is.finite(null_delta)))
  list(observed_delta = observed_delta, null_delta = null_delta,
       p_value = p_value)
}

complete_ids_for_periods <- function(data, features) {
  data %>%
    filter(period %in% periods) %>%
    filter(if_all(all_of(features), ~ is.finite(.x))) %>%
    count(id, period, name = "n_samples") %>%
    filter(n_samples == 1) %>%
    group_by(id) %>%
    summarise(has_all_periods = all(periods %in% period), .groups = "drop") %>%
    filter(has_all_periods) %>%
    pull(id)
}

period_matrix <- function(data, features, period_label, ids) {
  x <- data %>%
    filter(id %in% ids, .data$period == period_label) %>%
    arrange(match(id, ids))
  if (nrow(x) != length(ids) || anyDuplicated(x$id)) {
    stop("The period table is not one-sample-per-participant.")
  }
  list(ids = x$id, matrix = as.matrix(x[, features]))
}

run_pair_distance <- function(data, features, cohort, omics, metric,
                              distance_version, group_scope = "All",
                              n_perm = 10000, max_between = 25, seed = 1) {
  ids <- complete_ids_for_periods(data, features)
  if (group_scope != "All") {
    group_ids <- data %>% filter(pBMI_group == group_scope) %>% distinct(id) %>% pull(id)
    ids <- intersect(ids, group_ids)
  }
  if (length(ids) < 3) {
    return(list(summary = tibble(), pairs = tibble(), null = tibble()))
  }

  output <- map2(period_pairs, seq_along(period_pairs), function(pp, pair_index) {
    d1 <- period_matrix(data, features, pp[1], ids)
    d2 <- period_matrix(data, features, pp[2], ids)
    stopifnot(identical(d1$ids, d2$ids))
    dm <- cross_distance(d1$matrix, d2$matrix, metric)
    n <- length(ids)
    within <- diag(dm)
    off_diag <- row(dm) != col(dm)
    between_all <- dm[off_diag]
    distance_test <- permutation_identity_test(
      dm, n_perm = n_perm, seed = seed + pair_index
    )

    pair_label <- paste(pp, collapse = "-")
    valid_pairs <- map_dfr(seq_len(n), function(i) {
      candidates <- setdiff(seq_len(n), i)
      set.seed(seed + 10000 + pair_index * 100 + i)
      chosen <- sample(candidates, min(max_between, length(candidates)))
      bind_rows(
        tibble(cohort = cohort, omics = omics, distance_version = distance_version,
               metric = metric, group_scope = group_scope, period_pair = pair_label,
               pair_class = "within", anchor_id = ids[i], comparator_id = ids[i],
               distance = within[i]),
        tibble(cohort = cohort, omics = omics, distance_version = distance_version,
               metric = metric, group_scope = group_scope, period_pair = pair_label,
               pair_class = "between", anchor_id = ids[i], comparator_id = ids[chosen],
               distance = dm[i, chosen])
      )
    })

    summary <- tibble(
      cohort = cohort, omics = omics, distance_version = distance_version,
      metric = metric, group_scope = group_scope, period_pair = pair_label,
      n_individuals = n, n_features = length(features),
      n_between_distances = length(between_all), n_permutations = n_perm,
      within_mean = mean(within),
      within_q025 = unname(quantile(within, .025)), within_median = median(within),
      within_q975 = unname(quantile(within, .975)),
      between_mean = mean(between_all),
      between_q025 = unname(quantile(between_all, .025)),
      between_median = median(between_all),
      between_q975 = unname(quantile(between_all, .975)),
      mean_difference = distance_test$observed_delta,
      test_statistic = distance_test$observed_delta,
      permutation_unit = "participant_identity_at_second_timepoint",
      alternative = "within_mean_less_than_between_mean",
      permutation_null_delta_median = median(distance_test$null_delta, na.rm = TRUE),
      p_distance = distance_test$p_value
    )
    null <- tibble(
      cohort = cohort, omics = omics, distance_version = distance_version,
      metric = metric, group_scope = group_scope, period_pair = pair_label,
      permutation = seq_len(n_perm),
      permuted_mean_difference = distance_test$null_delta
    )
    list(summary = summary, pairs = valid_pairs, null = null,
         within = tibble(cohort = cohort, omics = omics,
                         distance_version = distance_version, metric = metric,
                         group_scope = group_scope, period_pair = pair_label,
                         id = ids, within_distance = within,
                         pBMI_group = data$pBMI_group[match(ids, data$id)]))
  })

  list(summary = map_dfr(output, "summary"), pairs = map_dfr(output, "pairs"),
       null = map_dfr(output, "null"), within = map_dfr(output, "within"))
}

permuted_pbmi_contrasts <- function(dm, group_labels, n_perm, seed, batch_size = 250) {
  if (any(!is.finite(dm))) {
    stop("The pBMI-label permutation requires a complete finite distance matrix.")
  }
  n <- length(group_labels)
  n_normal <- sum(group_labels == "Normal")
  n_owo <- sum(group_labels == "OWO")
  if (nrow(dm) != n || ncol(dm) != n || min(n_normal, n_owo) < 2) {
    stop("Invalid distance matrix or pBMI group sizes for label permutation.")
  }

  diagonal <- diag(dm)
  total_sum <- sum(dm)
  total_diagonal <- sum(diagonal)
  matrix_row_sums <- rowSums(dm)
  matrix_column_sums <- colSums(dm)
  base_membership <- c(rep(1, n_normal), rep(0, n_owo))
  null_contrast <- numeric(n_perm)

  set.seed(seed)
  for (batch_start in seq.int(1L, n_perm, by = batch_size)) {
    batch_index <- batch_start:min(batch_start + batch_size - 1L, n_perm)
    membership <- vapply(batch_index, function(i) {
      sample(base_membership, n, replace = FALSE)
    }, numeric(n))
    if (is.null(dim(membership))) membership <- matrix(membership, ncol = 1)

    normal_all_sum <- colSums(membership * (dm %*% membership))
    normal_diagonal_sum <- as.numeric(crossprod(diagonal, membership))
    normal_row_sum <- as.numeric(crossprod(matrix_row_sums, membership))
    normal_column_sum <- as.numeric(crossprod(matrix_column_sums, membership))
    owo_all_sum <- total_sum - normal_row_sum - normal_column_sum + normal_all_sum
    owo_diagonal_sum <- total_diagonal - normal_diagonal_sum

    delta_normal <- normal_diagonal_sum / n_normal -
      (normal_all_sum - normal_diagonal_sum) / (n_normal * (n_normal - 1))
    delta_owo <- owo_diagonal_sum / n_owo -
      (owo_all_sum - owo_diagonal_sum) / (n_owo * (n_owo - 1))
    null_contrast[batch_index] <- delta_normal - delta_owo
  }
  null_contrast
}

run_pbmi_specificity_tests <- function(data, features, cohort, omics, metric,
                                       distance_version, n_perm = 10000,
                                       seed = 1) {
  ids <- complete_ids_for_periods(data, features)
  group_map <- data %>%
    filter(id %in% ids, pBMI_group %in% c("Normal", "OWO")) %>%
    distinct(id, pBMI_group) %>%
    count(id, name = "n_group_labels") %>%
    filter(n_group_labels == 1) %>%
    select(id) %>%
    left_join(data %>%
                filter(pBMI_group %in% c("Normal", "OWO")) %>%
                distinct(id, pBMI_group), by = "id")
  ids <- ids[ids %in% group_map$id]
  group_labels <- group_map$pBMI_group[match(ids, group_map$id)]
  if (length(ids) < 4 || any(table(factor(group_labels,
                                          levels = c("Normal", "OWO"))) < 2)) {
    return(tibble())
  }

  map2_dfr(period_pairs, seq_along(period_pairs), function(pp, pair_index) {
    d1 <- period_matrix(data, features, pp[1], ids)
    d2 <- period_matrix(data, features, pp[2], ids)
    stopifnot(identical(d1$ids, d2$ids))
    dm <- cross_distance(d1$matrix, d2$matrix, metric)

    group_delta <- function(labels, group) {
      index <- which(labels == group)
      mean_distance_difference(dm[index, index, drop = FALSE])
    }
    delta_normal <- group_delta(group_labels, "Normal")
    delta_owo <- group_delta(group_labels, "OWO")
    observed_contrast <- delta_normal - delta_owo

    null_contrast <- permuted_pbmi_contrasts(
      dm, group_labels, n_perm = n_perm, seed = seed + pair_index
    )
    # Use a two-sided permutation test for the pBMI-group contrast
    p_value <- (1 + sum(abs(null_contrast) >= abs(observed_contrast),
                       na.rm = TRUE)) /
      (1 + sum(is.finite(null_contrast)))

    normal_component <- distance_components(
      dm[group_labels == "Normal", group_labels == "Normal", drop = FALSE]
    )
    owo_component <- distance_components(
      dm[group_labels == "OWO", group_labels == "OWO", drop = FALSE]
    )
    tibble(
      cohort = cohort, omics = omics, distance_version = distance_version,
      metric = metric, period_pair = paste(pp, collapse = "-"),
      n_normal = sum(group_labels == "Normal"),
      n_owo = sum(group_labels == "OWO"),
      normal_within_mean = mean(normal_component$within),
      normal_between_mean = mean(normal_component$between),
      delta_normal = delta_normal,
      owo_within_mean = mean(owo_component$within),
      owo_between_mean = mean(owo_component$between),
      delta_owo = delta_owo,
      delta_normal_minus_owo = observed_contrast,
      n_permutations = n_perm,
      permutation_unit = "participant_pBMI_label",
      alternative = "two_sided_pBMI_effect",
      permutation_null_contrast_median = median(null_contrast, na.rm = TRUE),
      p_pbmi_specificity = p_value
    )
  })
}

run_omics_analysis <- function(data, features, cohort, omics, seed = 1) {
  metric <- if (omics == "Metabolite") "euclidean" else "bray_curtis"
  all_results <- list()
  for (g in c("All", "Normal", "OWO")) {
    all_results[[g]] <- run_pair_distance(
      data, features, cohort, omics, metric, "primary", group_scope = g,
      n_perm = n_perm, max_between = max_between_per_anchor,
      seed = seed + match(g, c("All", "Normal", "OWO"))
    )
  }
  list(
    results = all_results,
    group_tests = run_pbmi_specificity_tests(
      data, features, cohort, omics, metric, "primary",
      n_perm = n_perm, seed = seed + 1000
    )
  )
}

# THSBC 
load("data/01_thsbc_cleaned_data/metabolites_1292_991.RData")
load("data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData")
load("data/01_thsbc_cleaned_data/gut_microbe_MetaG_1250_473.RData")
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")

thsbc_metab_features <- metab_uni
thsbc_16s_features <- bact_uni
thsbc_its_features <- fung_uni
thsbc_mgs_g_features <- mgs_g_uni
thsbc_mgs_s_features <- mgs_s_uni
thsbc_metab <- attach_metadata(metab_tab_st, thsbc_metab_features, total_pop, cov_tab)
thsbc_16s_raw <- attach_metadata(bact_fung_tab, thsbc_16s_features, total_pop, cov_tab)
thsbc_its_raw <- attach_metadata(bact_fung_tab, thsbc_its_features, total_pop, cov_tab)
thsbc_mgs_g_raw <- attach_metadata(mgs_tab, thsbc_mgs_g_features, total_pop, cov_tab)
thsbc_mgs_s_raw <- attach_metadata(mgs_tab, thsbc_mgs_s_features, total_pop, cov_tab)

thsbc_specs <- list(
  Metabolite = list(primary_data = thsbc_metab, primary_features = thsbc_metab_features),
  `16S_g` = list(primary_data = thsbc_16s_raw, primary_features = thsbc_16s_features),
  ITS_g = list(primary_data = thsbc_its_raw, primary_features = thsbc_its_features),
  MetaG_g = list(primary_data = thsbc_mgs_g_raw, primary_features = thsbc_mgs_g_features),
  MetaG_s = list(primary_data = thsbc_mgs_s_raw, primary_features = thsbc_mgs_s_features)
)
thsbc_out <- imap(thsbc_specs, ~ run_omics_analysis(
  .x$primary_data, .x$primary_features, "THSBC", .y,
  seed = 100 + match(.y, names(thsbc_specs))
))

# HBC 
load("data/02_hbc_cleaned_data/metabolites_393_158.RData")
load("data/02_hbc_cleaned_data/gut_microbe_fung_393_158.RData")
load("data/02_hbc_cleaned_data/gut_microbe_MetaG_393_158.RData")
load("data/02_hbc_cleaned_data/total_pop_hbc_393.Rdata")
load("data/02_hbc_cleaned_data/covariates_393.RData")

hbc_metab_features <- metab_uni
hbc_its_features <- fung_uni
hbc_mgs_g_features <- mgs_g_uni
hbc_mgs_s_features <- mgs_s_uni
hbc_metab <- attach_metadata(metab_tab_hbc_st, hbc_metab_features, total_pop_hbc, cov_tab)
hbc_its_raw <- attach_metadata(fung_tab_hbc, hbc_its_features, total_pop_hbc, cov_tab)
hbc_mgs_g_raw <- attach_metadata(mgs_tab_hbc, hbc_mgs_g_features, total_pop_hbc, cov_tab)
hbc_mgs_s_raw <- attach_metadata(mgs_tab_hbc, hbc_mgs_s_features, total_pop_hbc, cov_tab)

hbc_specs <- list(
  Metabolite = list(primary_data = hbc_metab, primary_features = hbc_metab_features),
  ITS_g = list(primary_data = hbc_its_raw, primary_features = hbc_its_features),
  MetaG_g = list(primary_data = hbc_mgs_g_raw, primary_features = hbc_mgs_g_features),
  MetaG_s = list(primary_data = hbc_mgs_s_raw, primary_features = hbc_mgs_s_features)
)
hbc_out <- imap(hbc_specs, ~ run_omics_analysis(
  .x$primary_data, .x$primary_features, "HBC", .y,
  seed = 500 + match(.y, names(hbc_specs))
))

flatten_out <- function(x) {
  rs <- map(x, "results") %>% flatten()
  list(
    summary = map_dfr(rs, "summary"), pairs = map_dfr(rs, "pairs"),
    null = map_dfr(rs, "null"), within = map_dfr(rs, "within"),
    group_tests = map_dfr(x, "group_tests")
  )
}
thsbc_results <- flatten_out(thsbc_out)
hbc_results <- flatten_out(hbc_out)

add_fdr <- function(x) {
  x %>% group_by(group_scope) %>% mutate(
    p_distance_fdr = p.adjust(p_distance, method = "fdr")
  ) %>% ungroup()
}
thsbc_results$summary <- add_fdr(thsbc_results$summary)
hbc_results$summary <- add_fdr(hbc_results$summary)
adjust_group_fdr <- function(x) {
  if (!nrow(x)) return(x)
  x %>% mutate(p_fdr = p.adjust(p_pbmi_specificity, method = "fdr"))
}

thsbc_results$group_tests <- adjust_group_fdr(thsbc_results$group_tests)
hbc_results$group_tests <- adjust_group_fdr(hbc_results$group_tests)

save(thsbc_results, file = "result/PART1/5.1-1_thsbc_omics_distance_comparison.RData")
save(hbc_results, file = "result/PART1/5.1-1_hbc_omics_distance_comparison.RData")

# External pregnancy metabolomics: target timepoints and gestational intervals
z_score <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
external <- readxl::read_xlsx("data/03_external_data/mmc2.xlsx")
external_features <- paste0("Metabolite", seq_len(264))
external_info <- tibble(Index = external_features, compound = names(external)[6:269])
names(external)[6:269] <- external_features
external <- external %>%
  transmute(id = as.character(`Subject ID`), wk = `Gestational age (GA)/weeks`,
            bw = `Birth GA/weeks`, Cohort, Labor_onset = `Labor onset`,
            across(all_of(external_features))) %>%
  filter(wk < bw) %>% arrange(Cohort, id, wk) %>%
  group_by(Cohort) %>% mutate(across(all_of(external_features), z_score)) %>% ungroup() %>%
  mutate(period = NA_character_, pBMI_group = "All")

thsbc_target_weeks <- thsbc_metab %>%
  filter(period %in% periods) %>% group_by(period) %>%
  summarise(target_wk = median(wk, na.rm = TRUE), .groups = "drop")
external_three_period <- tidyr::crossing(external %>% select(-period), thsbc_target_weeks) %>%
  mutate(wk_difference = abs(wk - target_wk)) %>%
  group_by(id, period) %>% slice_min(wk_difference, n = 1, with_ties = FALSE) %>%
  ungroup() %>% filter(wk_difference <= 2)
external_complete_ids <- complete_ids_for_periods(external_three_period, external_features)
external_three_period <- external_three_period %>% filter(id %in% external_complete_ids)
external_three <- run_pair_distance(external_three_period, external_features, "Liang_external",
                                     "Metabolite", "euclidean", "three_period_windows",
                                     n_perm = n_perm, max_between = max_between_per_anchor, seed = 700)

external_lag_breaks <- c(0.5, 1.5, 4, 8, 13, Inf)
external_lag_labels <- c("1_week", "2_4_weeks", "5_8_weeks", "9_13_weeks", "14plus_weeks")

run_external_high_frequency <- function(data, features, n_perm = 10000,
                                        ga_caliper = 1.5, max_pairs_per_subject = 50,
                                        seed = 800) {
  data <- data %>% arrange(Cohort, id, wk) %>% mutate(sample_index = row_number())
  ids <- unique(data$id)
  info <- data %>% distinct(id, Cohort) %>% arrange(match(id, ids))
  dm <- cross_distance(as.matrix(data[, features]), as.matrix(data[, features]), "euclidean")
  pair_data <- map_dfr(ids, function(subject_id) {
    s <- data %>% filter(id == subject_id) %>% arrange(wk)
    ij <- which(upper.tri(matrix(FALSE, nrow(s), nrow(s))), arr.ind = TRUE)
    tibble(id = subject_id, Cohort = s$Cohort[1],
           start_sample = s$sample_index[ij[, 1]], end_sample = s$sample_index[ij[, 2]],
           start_wk = s$wk[ij[, 1]], end_wk = s$wk[ij[, 2]],
           lag_weeks = s$wk[ij[, 2]] - s$wk[ij[, 1]]) %>%
      mutate(lag_group = cut(lag_weeks, external_lag_breaks, external_lag_labels,
                             right = TRUE))
  }) %>% filter(!is.na(lag_group))

  lag_out <- map2(external_lag_labels, seq_along(external_lag_labels), function(lag_label, k) {
    # Use a lag-specific seed and retain the selected pairs in the returned
    # object so the high-frequency analysis is exactly reproducible.
    set.seed(seed + 1000 + k)
    lp <- pair_data %>% filter(lag_group == lag_label) %>% group_by(id) %>%
      group_modify(~ slice_sample(.x, n = min(nrow(.x), max_pairs_per_subject))) %>%
      ungroup()
    subject_distance <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
    for (i in seq_along(ids)) {
      ap <- lp %>% filter(id == ids[i])
      if (!nrow(ap)) next
      for (j in seq_along(ids)) {
        if (info$Cohort[i] != info$Cohort[j]) next
        donor <- data %>% filter(id == ids[j])
        nearest_start <- vapply(ap$start_wk, function(w) {
          q <- abs(donor$wk - w); if (min(q) <= ga_caliper) donor$sample_index[which.min(q)] else NA_integer_
        }, integer(1))
        nearest_end <- vapply(ap$end_wk, function(w) {
          q <- abs(donor$wk - w); if (min(q) <= ga_caliper) donor$sample_index[which.min(q)] else NA_integer_
        }, integer(1))
        ok <- !is.na(nearest_start) & !is.na(nearest_end) & nearest_start != nearest_end
        if (any(ok)) subject_distance[i, j] <- median(dm[cbind(ap$start_sample[ok], nearest_end[ok])])
      }
    }
    cohort_size <- table(info$Cohort)
    valid_cohorts <- names(cohort_size)[cohort_size >= 2]
    eligible <- which(
      info$Cohort %in% valid_cohorts & is.finite(diag(subject_distance)) &
        rowSums(is.finite(subject_distance) & row(subject_distance) != col(subject_distance)) > 0
    )
    if (length(eligible) < 2) return(tibble())
    subject_distance <- subject_distance[eligible, eligible, drop = FALSE]
    ids <- ids[eligible]
    info <- info[eligible, , drop = FALSE]
    lp_used <- lp %>% filter(id %in% ids)
    within <- diag(subject_distance)
    between <- vapply(seq_along(ids), function(i) {
      ind <- which(info$Cohort == info$Cohort[i])
      median(subject_distance[i, setdiff(ind, i)], na.rm = TRUE)
    }, numeric(1))
    off_diag <- row(subject_distance) != col(subject_distance)
    between_all <- subject_distance[off_diag & is.finite(subject_distance)]
    distance_test <- permutation_identity_test(
      subject_distance, n_perm = n_perm, seed = seed + k,
      strata = info$Cohort
    )
    list(
      summary = tibble(analysis = "high_frequency", lag_group = lag_label,
        n_individuals = length(ids), n_within_pairs = nrow(lp_used),
        n_between_distances = length(between_all), n_permutations = n_perm,
        within_mean = mean(within),
        within_q025 = quantile(within, .025), within_median = median(within),
        within_q975 = quantile(within, .975),
        between_mean = mean(between_all),
        between_q025 = quantile(between_all, .025), between_median = median(between_all),
        between_q975 = quantile(between_all, .975),
        mean_difference = distance_test$observed_delta,
        test_statistic = distance_test$observed_delta,
        permutation_unit = "participant_identity_within_cohort",
        alternative = "within_mean_less_than_between_mean",
        sampling_seed = seed + 1000 + k,
        permutation_null_delta_median = median(distance_test$null_delta, na.rm = TRUE),
        p_distance = distance_test$p_value),
      subject = tibble(lag_group = lag_label, id = ids, Cohort = info$Cohort,
                       within_distance = within, between_reference = between),
      selected_pairs = lp_used,
      null = tibble(lag_group = lag_label, permutation = seq_len(n_perm),
                    permuted_mean_difference = distance_test$null_delta)
    )
  })
  list(summary = map_dfr(lag_out, "summary") %>% mutate(
          p_distance_fdr = p.adjust(p_distance, "fdr")),
       subject = map_dfr(lag_out, "subject"),
       selected_pairs = map_dfr(lag_out, "selected_pairs"),
       null = map_dfr(lag_out, "null"))
}

external_high_frequency <- run_external_high_frequency(external, external_features, n_perm = n_perm)
save(external_info, thsbc_target_weeks, external_three, external_high_frequency,
     file = "result/PART1/5.1-1_external_liang_distance_comparison.RData")