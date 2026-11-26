####################################################
# Visualization and statistical testing for Figure 2
#####################################################

# Figure 2a Contributions of maternal factors to variation in metabolome and microbiome ------------------------------
rm(list = ls())
library(tidyverse)
library(purrr)
library(data.table)
library(readxl)
library(writexl)
library(ggplot2)
library(ggh4x)
library(patchwork)

load("result/PART1/1-1 Permanova_Single_Phenotypes_and_Phenotype_Clusters_avaliable_pop.rdata")

df = Permanova_Single_Phenotypes %>% 
  filter(period != 'TA') %>% 
  filter(! variable %in% c("pbmi_cut3", "wc_t1_cut", "whr_t1_cut", "bfr_t1_cut")) %>% 
  mutate(var_name = case_when(
    variable == 'pbmi' ~ 'pBMI', 
    variable == 'wc_t1' ~ 'Waist (T1)',
    variable == 'whr_t1' ~ 'WHR (T1)',
    variable == 'bfr_t1' ~ 'BFR (T1)', 
    variable == 'age' ~ 'Age', 
    variable == 'parity' ~ 'Parity',
    variable == 'smk' ~ 'Smoking status', 
    variable == 'drk' ~ 'Alcohol intake', 
    variable == 'edu' ~ 'Education level', 
    variable == 'fd_grain' ~ 'Grain intake',
    variable == 'fd_meat' ~ 'Meat intake', 
    variable == 'fd_vege' ~ 'Vegetable intake', 
    variable == 'fd_fruit' ~ 'Fruit intake', 
    variable == 'fd_egg' ~ 'Egg intake',
    variable == 'fd_milkpro' ~ 'Milk intake',
    variable == 'physic_activity' ~ 'Physical activity',
    variable == 'medical_insurance' ~ 'Medical insurance',
    variable == 'income' ~ 'Annual income',
    variable == 'wk' ~ 'Gestational week',
    variable == 'antibiotic' ~ 'Antibiotic use',
    TRUE~ NA
  ), 
  text = case_when(P_value < 0.001 ~ '***',
                  P_value < 0.01 ~ '**',
                  P_value < 0.05 ~ '*', 
                  P_value < 0.1 ~ '+', .default = ''))
df$group = factor(df$type, levels = c('Metabolites','16S_g', 'ITS_g', 'MetaG_g', 'MetaG_s'),
                  labels = c('Metabolome','16S_g', 'ITS_g', 'MetaG_g', 'MetaG_s'))
var_order = c('pBMI', 'Waist (T1)', 'WHR (T1)', 'BFR (T1)', 
            'Age', 'Parity', 
            'Smoking status', 'Alcohol intake', 'Physical activity', 
            'Grain intake', 'Meat intake', 'Vegetable intake', 'Fruit intake', 'Egg intake', 'Milk intake', 
            'Education level', 'Medical insurance', 'Annual income', 'Antibiotic use', 'Gestational week')
df$var_name = factor(df$var_name, levels=rev(var_order))

# period-specific
range(df$R2*100)
col_group = colorspace::lighten(c("#7EABCA","#FFB2B9", "#95C6C6","#ebdcb2", "#d4cfdf"), amount = 0.3)
plot_upper <- ggplot(df, aes(x=period, y=var_name))+
  facet_wrap2(~group, nrow = 1, 
              strip = strip_themed(background_x = elem_list_rect(fill = col_group, size=rep(0.5,2))))+
  geom_point(aes(color=R2*100, size=R2*100), shape = 16)+
  scale_size_continuous(range = c(0.1,6.5))+
  scale_color_gradient(low = "#f7f7f7", high = "#bf3e3e", breaks=c(0.3,0.6,0.9))+
  geom_text(data=df %>% filter(text!='+'), aes(label=text), vjust = 0.7, hjust = 0.5)+
  geom_text(data=df %>% filter(text=='+'), aes(label=text), hjust = 0.5, size = 2.9)+
  labs(color = 'Variance explained (%)')+
  guides(size = 'none', color = guide_colorbar(
    barheight = unit(4, "mm"),
    barwidth = unit(23, "mm"),
    ticks.linewidth = 0.2
  ))+
  theme_bw(base_line_size = 0.3)+
  theme(axis.text = element_text(color='black', size=9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 2),
        legend.direction = 'horizontal',
        legend.title.position = 'top',
        legend.background = element_rect(fill='transparent'),
        legend.text = element_text(size=6.5, vjust = 3),
        legend.title = element_text(size=7.5),
        legend.key.size = unit(3,'mm'),
        axis.title = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size = 9, color = 'white',face = 'bold'))

df = Permanova_Phenotype_Clusters %>% 
  filter(period != 'TA', variable != 'all') %>% 
  mutate(var_name = case_when(
    variable == 'anthropometric' ~ 'Anthropometric', 
    variable == 'demography' ~ 'Demographic',
    variable == 'behavioral' ~ 'Behavioral',  
    variable == 'dietary' ~ 'Dietary',
    variable == 'social_determ' ~ 'Socioeconomic',
    TRUE~ NA
  ), R2_adj = ifelse(R2_adj < 0, 0, R2_adj),
  text = case_when(P_value < 0.001 ~ '***',
                  P_value < 0.01 ~ '**',
                  P_value < 0.05 ~ '*', 
                  P_value < 0.1 ~ '+', .default = ''))
df$group = factor(df$type, levels = c('Metabolites','16S_g', 'ITS_g', 'MetaG_g', 'MetaG_s'),
                  labels = c('Metabolome','16S_g', 'ITS_g', 'MetaG_g', 'MetaG_s'))
var_order = c('Anthropometric', 'Demographic', 'Behavioral', 'Dietary', 'Socioeconomic')
df$var_name = factor(df$var_name, levels=rev(var_order))

# period-specific
range(df$R2_adj)
plot_lower <- ggplot(df, aes(x=period, y=var_name))+
  facet_wrap2(~group, nrow = 1)+
  geom_point(aes(color=R2_adj*100, size=R2_adj*100))+
  scale_size_continuous(range = c(0.1,6.5))+
  scale_color_gradient(low = "#f7f7f7", high = "#7b3b58", breaks=c(0.3,0.6,0.9))+
  geom_text(data=df %>% filter(text!='+'), aes(label=text), vjust = 0.7, hjust = 0.5)+
  geom_text(data=df %>% filter(text=='+'), aes(label=text), hjust = 0.5, size = 2.9)+
  labs(color = 'Variance explained by group (%)')+
  guides(size = 'none', color = guide_colorbar(
    barheight = unit(4, "mm"),
    barwidth = unit(23, "mm"),
    ticks.linewidth = 0.2
  )) +
  theme_bw(base_line_size = 0.3)+
  theme(axis.text = element_text(color='black', size=9),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 2),
        legend.direction = 'horizontal',
        legend.title.position = 'top',
        legend.background = element_rect(fill='transparent'),
        legend.text = element_text(size=6.5, vjust = 3),
        legend.title = element_text(size=7.5),
        legend.key.size = unit(3,'mm'),
        axis.title = element_blank(),
        plot.background = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank())

plot_combined <- plot_upper / plot_lower +
  plot_layout(heights = c(18, 5), guides = 'collect') &
  theme(legend.position = 'bottom')
ggsave('figure/PART1/fig2a_variance_explained_metab_micro_combined.pdf', plot_combined, width = 5.5, height = 6)


# Figure 2b serum metabolomic composition ------------------------------

rm(list = ls())

library(tidyverse)
library(vegan)
library(permute)
library(ggpubr)
library(patchwork)

set.seed(1111)
n_permutations <- 999
n_cores <- 12

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", paste0("=", formatC(p, digits = 3, format = "f")))
}

# Permute complete subject profiles between participants while retaining the
# original order of periods within each subject.
make_subject_permutations <- function(metadata, nperm = 999, seed = 1111) {
  rows_by_id <- split(seq_len(nrow(metadata)), metadata$id)
  period_patterns <- map_chr(
    rows_by_id,
    ~ paste(as.character(metadata$period[.x]), collapse = "|")
  )
  if (length(unique(lengths(rows_by_id))) != 1 ||
      length(unique(period_patterns)) != 1) {
    stop("Whole-subject permutations require the same periods for every id.")
  }

  set.seed(seed)
  subject_orders <- permute::shuffleSet(length(rows_by_id), nset = nperm)
  t(apply(subject_orders, 1, function(subject_order) {
    unlist(rows_by_id[subject_order], use.names = FALSE)
  }))
}

run_joint_permanova <- function(distance, metadata, nperm = 999, cores = 1) {
  # Period: exchange timepoints only within the same participant.
  ctrl_period <- permute::how( nperm = nperm, 
    plots = permute::Plots( strata = as.factor(metadata$id), type = "none" ),
    within = permute::Within(type = "free")
  )

  set.seed(1111)
  period_fit <- vegan::adonis2(
    distance ~ period + pbmi_cut3, data = metadata,
    permutations = ctrl_period, by = "margin",
    parallel = cores
  )

  # pBMI: exchange all repeated observations from a participant as one block;
  # matching period positions are retained, so the pBMI test controls period.
  subject_permutations <- make_subject_permutations(metadata, nperm, 1111)
  set.seed(1111)
  pbmi_fit <- vegan::adonis2(
    distance ~ period + pbmi_cut3, data = metadata,
    permutations = subject_permutations, by = "margin",
    parallel = cores
  )

  results <- bind_rows(
    tibble( variable = "period", F_value = period_fit["period", "F"],
      R2 = period_fit["period", "R2"], p_value = period_fit["period", "Pr(>F)"],
      permutation = "within participant"
    ),
    tibble( variable = "pbmi_cut3", F_value = pbmi_fit["pbmi_cut3", "F"],
      R2 = pbmi_fit["pbmi_cut3", "R2"], p_value = pbmi_fit["pbmi_cut3", "Pr(>F)"],
      permutation = "whole participant"
    )
  ) %>%
    mutate(
      n = nrow(metadata),
      n_id = n_distinct(metadata$id),
      model = "distance ~ period + pbmi_cut3; marginal tests"
    )

  list(
    results = results,
    period_fit = period_fit,
    pbmi_fit = pbmi_fit,
    subject_permutations = subject_permutations
  )
}
# Data preparation 
load("data/01_thsbc_cleaned_data/metabolites_1292_991.RData")
load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/covariates_5212.RData")

DF_metab <- total_pop %>%
  inner_join(metab_tab_st, by = "id") %>%
  inner_join(cov_tab, by = c("id", "period"))

metab_features <- intersect(metab_uni, names(DF_metab))

analysis_data <- DF_metab %>%
  filter( period %in% c("T1", "T2", "T3"), !is.na(pbmi_cut3),
    if_all(all_of(metab_features), ~ is.finite(.x))
  ) %>%
  group_by(id) %>%
  filter(n_distinct(period) == 3, n() == 3) %>%
  ungroup() %>%
  mutate(
    period = factor(period, levels = c("T1", "T2", "T3")),
    pbmi_cut3 = factor(pbmi_cut3, levels = c("Underweight", "Normal", "OWO"))
  ) %>%
  arrange(id, period)

if (anyDuplicated(analysis_data[c("id", "period")])) {
  stop("The analysis data contain duplicated id-period rows.")
}

metab_matrix <- analysis_data %>%
  select(all_of(metab_features)) %>%
  as.matrix()

# Metabolites are already log-z transformed; no additional scaling is applied.
pca_fit <- prcomp(metab_matrix, center = FALSE, scale. = FALSE)
variance_explained <- 100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2)
pca_scores <- analysis_data %>%
  select(id, period, pbmi_cut3) %>%
  bind_cols(as_tibble(pca_fit$x[, 1:2, drop = FALSE]))

euclidean_distance <- vegan::vegdist(metab_matrix, method = "euclidean")
joint_permanova <- run_joint_permanova(
  euclidean_distance, analysis_data, n_permutations, n_cores
)
permanova_results <- joint_permanova$results

dir.create("result/PART1", recursive = TRUE, showWarnings = FALSE)
save(
  pca_fit, pca_scores, variance_explained, euclidean_distance,
  permanova_results,
  file = "result/PART1/fig2b_serum_metabolite_composition.rdata"
)

disper_pbmi <- vegan::betadisper(euclidean_distance, analysis_data$pbmi_cut3)
anova(disper_pbmi)

disper_period <- vegan::betadisper(euclidean_distance, analysis_data$period)
anova(disper_period)

# Plotting 
group_colours <- c( Underweight = "#c0b9e6", Normal = "#b0d2e8", OWO = "#f7c8ae" )
period_shapes <- c(T1 = 16, T2 = 17, T3 = 15)

pbmi_result <- permanova_results %>% filter(variable == "pbmi_cut3")
period_result <- permanova_results %>% filter(variable == "period")
permanova_annotation <- paste0(
  "R²=", sprintf("%.3f", pbmi_result$R2),
  ", P", format_p(pbmi_result$p_value), " for pBMI group\n",
  "R²=", sprintf("%.3f", period_result$R2),
  ", P", format_p(period_result$p_value), " for timepoints"
)

p_pca <- ggplot(pca_scores, aes(PC1, PC2, colour = pbmi_cut3)) +
  geom_point(aes(shape = period), size = 1.5, alpha = 0.8) +
  stat_ellipse(
    aes(fill = pbmi_cut3), geom = "polygon", level = 0.95,
    linetype = 2, linewidth= 0.8, alpha = 0.2, show.legend = FALSE
  ) +
  annotate(
    "text", x = -21, y = 16, label = permanova_annotation,
    hjust = 0, vjust = 0.8, size = 6, lineheight = 1.05
  ) +
  scale_x_continuous(limits = c(-21, 20))+
  scale_y_continuous(limits = c(-15, 16))+
  scale_colour_manual(values = group_colours, guide = "none") +
  scale_fill_manual(values = group_colours, guide = "none") +
  scale_shape_manual(values = period_shapes, name = 'Timepoint') +
  labs( x = sprintf("PC1 (%.1f%%)", variance_explained[1]), y = sprintf("PC2 (%.1f%%)", variance_explained[2]), shape = NULL ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 11) +
  theme(
    panel.background = element_blank(), 
    plot.background = element_blank(), 
    panel.grid = element_blank(),
    legend.position = c(0.92, 0.1),
    legend.title = element_text(size = 12),
    legend.background = element_blank(), 
    legend.key.size = unit(14, "pt"), 
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 13, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    plot.margin = margin(3,3,3,3)
  ) 

library(ggpubr)
p_pc1 <- ggplot( pca_scores, aes(pbmi_cut3, PC1, fill = pbmi_cut3, colour = pbmi_cut3)) +
  facet_wrap(. ~ period, nrow = 1, ncol = 3) +
  geom_boxplot( linewidth = 0.8, width = 0.3, alpha = 0.25, outlier.shape = NA ) +
  stat_compare_means(aes(group = pbmi_cut3), method = "wilcox.test", 
                     label.x.npc = 0.5,label.y.npc = 0.96,
                     ref.group = 'Normal',
                     label = "p.signif", show.legend = F, 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                     size = 5.5) +
  scale_colour_manual(values = group_colours, guide = "none") +
  scale_fill_manual(values = group_colours, guide = "none") +
  labs(x = NULL, y = "PC1") +
  theme_bw(base_size = 11) +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8),
        panel.background = element_blank(),
        axis.title = element_text(color = "black",size = 13),
        axis.text = element_text(color = "black",size = 12),
        axis.text.x = element_text(color = "black",angle = 30, size = 12, hjust = 1),
        strip.text = element_text(color = "black",size = 13),
        legend.title = element_blank(),
        legend.key.size = unit(12, "pt"), 
        legend.text = element_text(size = 10),
        legend.position = 'none')

# Fig. 2b layout: PCA above and the three PC1 boxplots below.
fig2b <- p_pca / p_pc1 + patchwork::plot_layout(heights = c(2, 1))

dir.create("figure/PART1", recursive = TRUE, showWarnings = FALSE)
ggsave(
  "figure/PART1/fig2b_serum_metabolite_composition.pdf",
  fig2b, width = 6.2, height = 9
)

# Figure 2c-d: alpha diversity and F/B log-ratio ridge plots ------------------------------
rm(list = ls())

library(tidyverse)
library(ggridges)
library(ggpubr)
library(patchwork)

pbmi_levels <- c("Underweight", "Normal", "OWO")
pbmi_colours <- c( Underweight = "#c0b9e6", Normal = "#b0d2e8", OWO = "#f7c8ae")

# 1. Statistical helper functions 
p_to_symbol <- function(p) {
  case_when( is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns" )
}

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", paste0("=", formatC(p, digits = 3, format = "f")))
}

# 2. Load precomputed alpha-diversity and F/B comparison results 
load("result/PART1/3-7_microbiome_group_comparisons_thsbc.RData")
load("result/PART1/3-8_microbiome_group_comparisons_hbc.RData")

micro_group_data_raw <- bind_rows(
  micro_group_data_raw_thsbc,
  micro_group_data_raw_hbc
)
micro_group_data_rank <- bind_rows(
  micro_group_data_rank_thsbc,
  micro_group_data_rank_hbc
)
micro_group_overall_raw <- bind_rows(
  micro_group_overall_raw_thsbc,
  micro_group_overall_raw_hbc
)
micro_group_pairwise_raw <- bind_rows(
  micro_group_pairwise_raw_thsbc,
  micro_group_pairwise_raw_hbc
)
micro_group_trend_raw <- bind_rows(
  micro_group_trend_raw_thsbc,
  micro_group_trend_raw_hbc
)
micro_group_overall_rank <- bind_rows(
  micro_group_overall_rank_thsbc,
  micro_group_overall_rank_hbc
)
micro_group_pairwise_rank <- bind_rows(
  micro_group_pairwise_rank_thsbc,
  micro_group_pairwise_rank_hbc
)
micro_group_trend_rank <- bind_rows(
  micro_group_trend_rank_thsbc,
  micro_group_trend_rank_hbc
)

standardize_plot_data <- function(data) {
  data %>%
    mutate(
      cohort = factor(cohort, levels = c("THSBC", "HBC")),
      period = factor(period, levels = c("T1", "T2", "T3")),
      pbmi_cut3 = factor(pbmi_cut3, levels = pbmi_levels),
      omics = factor(
        omics,
        levels = c("16S_g", "ITS_g", "MetaG_g", "MetaG_s")
      ),
      metric = factor(metric, levels = c("richness", "shannon", "fb_ratio"))
    )
}

alpha_data <- micro_group_data_raw %>%
  filter(data_type == "alpha_diversity") %>%
  standardize_plot_data()
alpha_rank_data <- micro_group_data_rank %>%
  filter(data_type == "alpha_diversity") %>%
  standardize_plot_data()
fb_data <- micro_group_data_raw %>%
  filter(data_type == "fb_ratio") %>%
  standardize_plot_data()
fb_rank_data <- micro_group_data_rank %>%
  filter(data_type == "fb_ratio") %>%
  standardize_plot_data()

alpha_group_overall <- micro_group_overall_raw %>%
  filter(group_variable == "pbmi_cut3", data_type == "alpha_diversity")
alpha_group_pairwise <- micro_group_pairwise_raw %>%
  filter(group_variable == "pbmi_cut3", data_type == "alpha_diversity")
alpha_trend <- micro_group_trend_raw %>%
  filter(group_variable == "pbmi_cut3", data_type == "alpha_diversity")
alpha_rank_group_overall <- micro_group_overall_rank %>%
  filter(group_variable == "pbmi_cut3", data_type == "alpha_diversity")
alpha_rank_group_pairwise <- micro_group_pairwise_rank %>%
  filter(group_variable == "pbmi_cut3", data_type == "alpha_diversity")
alpha_rank_trend <- micro_group_trend_rank %>%
  filter(group_variable == "pbmi_cut3", data_type == "alpha_diversity")

fb_group_overall <- micro_group_overall_raw %>%
  filter(group_variable == "pbmi_cut3", data_type == "fb_ratio")
fb_group_pairwise <- micro_group_pairwise_raw %>%
  filter(group_variable == "pbmi_cut3", data_type == "fb_ratio")
fb_trend <- micro_group_trend_raw %>%
  filter(group_variable == "pbmi_cut3", data_type == "fb_ratio")
fb_rank_group_overall <- micro_group_overall_rank %>%
  filter(group_variable == "pbmi_cut3", data_type == "fb_ratio")
fb_rank_group_pairwise <- micro_group_pairwise_rank %>%
  filter(group_variable == "pbmi_cut3", data_type == "fb_ratio")
fb_rank_trend <- micro_group_trend_rank %>%
  filter(group_variable == "pbmi_cut3", data_type == "fb_ratio")

# 3. Plotting helper functions 
make_annotation_data <- function(plot_data, overall, trend, facet_vars) {
  positions <- plot_data %>%
    group_by(across(all_of(facet_vars))) %>%
    summarise( x = median(value, na.rm = TRUE), .groups = "drop" )

  positions %>%
    left_join(overall, by = facet_vars) %>%
    left_join(trend, by = facet_vars) %>%
    mutate( 
      y = 3.35,
      label = paste0("Ptrend", format_p(p_trend))
    )
}

make_star_data <- function(plot_data, pairwise, facet_vars) {
  positions <- plot_data %>%
    group_by(across(all_of(facet_vars))) %>%
    summarise(x = quantile(value, 0.98, na.rm = TRUE), .groups = "drop")

  pairwise %>%
    filter(
      significance != "ns",
      (group_1 == "Underweight" & group_2 == "Normal") |
        (group_1 == "Normal" & group_2 == "OWO")
    ) %>%
    mutate(
      target_group = if_else(group_1 == "Normal", group_2, group_1),
      y = if_else(target_group == "Underweight", 1, 3)
    ) %>%
    left_join(positions, by = facet_vars)
}

plot_ridges <- function(plot_data, overall, pairwise, trend, facet_var, x_label, show_facet_labels = TRUE) {
  facet_vars <- facet_var
  annotation_data <- make_annotation_data(
    plot_data, overall, trend, facet_vars
  )
  star_data <- make_star_data(plot_data, pairwise, facet_vars)

  ggplot(plot_data, aes(value, pbmi_cut3, fill = pbmi_cut3)) +
    ggridges::geom_density_ridges(
      colour = "white", alpha = 0.9, rel_min_height = 0.02, 
      quantile_lines = TRUE, quantiles = 2, scale = 1.2
    ) +
    geom_text( data = annotation_data,
      aes(x = Inf, y = y, label = label),
      inherit.aes = FALSE, hjust = 1.2, vjust = -0.5, size = 4.5
    ) +
    geom_text( data = star_data, aes(x = x, y = y, label = significance),
      inherit.aes = FALSE, hjust = 0, size = 6 ) +
    scale_fill_manual(values = pbmi_colours) +
    labs(x = x_label, y = NULL) +
    coord_cartesian(clip = "off") +
    ggpubr::theme_classic2(base_size = 10) +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(color = "black", size = 11),
      axis.title.x = element_text(size = 13, colour = "black"),
      strip.placement = "outside",
      strip.text.y.left = if (show_facet_labels) {
        element_text(angle = 90, size = 13, colour = "black")
      } else {
        element_blank()
      },
      strip.background = element_blank(),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    )
}

# 4. Fig. 2c: THSBC and HBC alpha diversity at T1 
alpha_t1_rank_data <- alpha_rank_data %>%
  filter(period == "T1", !is.na(value)) %>%
  droplevels()
alpha_t1_pairwise <- alpha_rank_group_pairwise %>%
  filter(period == "T1")
alpha_t1_trend <- alpha_rank_trend %>%
  filter(period == "T1")

alpha_t1_summary <- alpha_t1_rank_data %>%
  group_by(cohort, metric, omics, pbmi_cut3) %>%
  summarise( n = n(), mean = mean(value),
    standard_error = sd(value) / sqrt(n),
    .groups = "drop"
  )

# Mark significant Underweight/OWO comparisons against Normal above the
# corresponding box. Pairwise P values are FDR adjusted within each stratum.
alpha_t1_pair_significance <- alpha_t1_pairwise %>%
  filter( p_adjusted < 0.05,
    (group_1 == "Underweight" & group_2 == "Normal") | (group_1 == "Normal" & group_2 == "OWO")
  ) %>%
  mutate(
    pbmi_cut3 = factor( if_else(group_1 == "Normal", group_2, group_1), levels = pbmi_levels ),
    label = p_to_symbol(p_adjusted)
  ) %>%
  ungroup() %>%
  select(cohort, metric, omics, pbmi_cut3, label)

# Include blank rows for all three groups so position_dodge retains the same
# horizontal slots even when only one comparison is significant.
alpha_t1_pair_labels <- alpha_t1_summary %>%
  left_join(
    alpha_t1_pair_significance,
    by = c("cohort", "metric", "omics", "pbmi_cut3")
  ) %>%
  mutate(
    label = replace_na(label, ""),
    label_y = pmax(mean + standard_error, 0) + 0.10
  )

# Display the trend P value above the Normal box for each omics layer.
alpha_t1_trend_text <- alpha_t1_trend %>%
  mutate(
    pbmi_cut3 = factor("Normal", levels = pbmi_levels),
    label = paste0("Ptrend", format_p(p_trend))
  ) %>%
  ungroup() %>%
  select(cohort, metric, omics, pbmi_cut3, label)

alpha_t1_trend_labels <- alpha_t1_summary %>%
  left_join(
    alpha_t1_trend_text,
    by = c("cohort", "metric", "omics", "pbmi_cut3")
  ) %>%
  mutate(
    label = replace_na(label, ""),
    label_y = pmax(mean + standard_error, 0) + 0.10
  )

position_pbmi <- position_dodge(width = 0.62)
fig2c <- ggplot( alpha_t1_rank_data, aes(x = omics, y = value, fill = pbmi_cut3) ) +
  geom_boxplot(
    position = position_pbmi, width = 0.5,
    outlier.shape = NA, alpha = 0.82,
    colour = "black", linewidth = 0.30
  ) +
  geom_text(
    data = alpha_t1_pair_labels,
    aes(x = omics, y = Inf, label = label, group = pbmi_cut3),
    position = position_pbmi, inherit.aes = FALSE, size = 6, vjust = 2.2,
  ) +
  geom_text(
    data = alpha_t1_trend_labels,
    aes(x = omics, y = -Inf, label = label, group = pbmi_cut3),
    position = position_pbmi, inherit.aes = FALSE, vjust = -1, hjust = 0.5, size = 3
  ) +
  facet_grid(
    rows = vars(cohort),
    cols = vars(metric),
    as.table = FALSE,
    labeller = labeller( metric = c(richness = "Richness (T1)", shannon = "Shannon (T1)") )
  ) +
  scale_fill_manual(values = pbmi_colours, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)), limits = c(-3, 3)) +
  labs(x = NULL, y = "Rank-normalized alpha diversity") +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 10) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 13, color = "black"),
    strip.background = element_rect(colour = "black", linewidth = 0.6),
    axis.text = element_text(color = "black", size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = c(0.08, 0.1), 
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_blank()
  )
fig2c

# 5. Fig. 2d: THSBC F/B log-ratio across three periods 
fb_plot_data <- fb_data %>%
  filter(cohort == "THSBC", !is.na(value)) %>%
  droplevels()
fb_overall_thsbc <- fb_group_overall %>% filter(cohort == "THSBC")
fb_pairwise_thsbc <- fb_group_pairwise %>% filter(cohort == "THSBC")
fb_trend_thsbc <- fb_trend %>% filter(cohort == "THSBC")

p_fb_16s <- plot_ridges(
  plot_data = fb_plot_data %>% filter(omics == "16S_g"),
  overall = fb_overall_thsbc %>% filter(omics == "16S_g"),
  pairwise = fb_pairwise_thsbc %>% filter(omics == "16S_g"),
  trend = fb_trend_thsbc %>% filter(omics == "16S_g"),
  facet_var = "period",
  x_label = "F/B log-ratio (16S rRNA)",
  show_facet_labels = TRUE
) +
  facet_grid(period ~ ., switch = "y", scales = "free_x") +
  scale_x_continuous(expand = c(0.02,0), limits = c(0, 8))

p_fb_metag <- plot_ridges(
  plot_data = fb_plot_data %>% filter(omics == "MetaG_g"),
  overall = fb_overall_thsbc %>% filter(omics == "MetaG_g"),
  pairwise = fb_pairwise_thsbc %>% filter(omics == "MetaG_g"),
  trend = fb_trend_thsbc %>% filter(omics == "MetaG_g"),
  facet_var = "period",
  x_label = "F/B log-ratio (Metagenome)",
  show_facet_labels = FALSE
) +
  facet_grid(period ~ ., switch = "y", scales = "free_x") +
  scale_x_continuous(expand = c(0.02,0), limits = c(-0.5, 11))

fig2d <- p_fb_16s + p_fb_metag + patchwork::plot_layout(ncol = 2)

dir.create("figure/PART1", recursive = TRUE, showWarnings = FALSE)
ggsave("figure/PART1/fig2c_alpha_diversity_boxplot.pdf", fig2c, width = 8, height = 4.5)
ggsave("figure/PART1/fig2d_fb_ratio_ridgeplot.pdf", fig2d, width = 6, height = 4.5)


# Figure 2e: multi-omics distance analysis ------------------------------
rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

result_dir <- "result/PART1"
figure_dir <- "figure/PART1"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

load_object <- function(path, object_name) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists(object_name, envir = e, inherits = FALSE)) {
    stop("Object '", object_name, "' was not found in ", path)
  }
  get(object_name, envir = e, inherits = FALSE)
}

thsbc_results <- load_object(
  file.path(result_dir, "5.1-1_thsbc_omics_distance_comparison.RData"),
  "thsbc_results"
)
hbc_results <- load_object(
  file.path(result_dir, "5.1-1_hbc_omics_distance_comparison.RData"),
  "hbc_results"
)
external_objects <- new.env(parent = emptyenv())
load(file.path(result_dir, "5.1-1_external_liang_distance_comparison.RData"),
     envir = external_objects)
external_three <- get("external_three", envir = external_objects)
external_high_frequency <- get("external_high_frequency", envir = external_objects)

period_levels <- c("T1-T2", "T2-T3", "T1-T3")
period_labels <- c("T1-T2", "T2-T3", "T1-T3")
omics_levels <- c("Metabolite", "16S_g", "ITS_g", "MetaG_g", "MetaG_s")
omics_labels <- c(Metabolite = "Metabolome", `16S_g` = "16S genus",
                  ITS_g = "Fungal genus", MetaG_g = "Metagenomic genus",
                  MetaG_s = "Metagenomic species")
omics_display_levels <- unname(omics_labels)
omics_factor <- function(x) {
  factor(unname(omics_labels[as.character(x)]), levels = omics_display_levels)
}
omics_facet_labels <- c(
  Metabolome = "Metabolome", `16S genus` = "16S-g",
  `Fungal genus` = "ITS-g", `Metagenomic genus` = "MetaG-g",
  `Metagenomic species` = "MetaG-s"
)
pair_levels <- c("within", "between")
pair_labels <- c(within = "Within individual", between = "Between individuals")
pair_colours <- c("Within individual" = "#E6956F", "Between individuals" = "#788FCE")
pbmi_colours <- c(Normal = "#b0d2e8", OWO = "#f7d4bc")
interval_levels <- c("1_week", "2_4_weeks", "5_8_weeks", "9_13_weeks", "14plus_weeks")
interval_labels <- c("1", "2-4", "5-8", "9-13", ">=14")

p_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

base_theme <- theme_bw(base_size = 9) + theme(
  plot.background = element_blank(),
  panel.background = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text = element_text(size = 12, color = "black"),
  strip.background = element_rect(colour = "black"),
  axis.text = element_text(color = "black", size = 11),
  axis.title = element_text(size = 12),
  legend.title = element_blank(),
  legend.text = element_text(size = 9),
  legend.key = element_blank(),
  legend.key.size = unit(0.4, "cm"),
  legend.background = element_blank()
)

thin_plot_data <- function(x, max_n = 3000, seed = 111) {
  set.seed(seed)
  x %>% group_by(omics, period_pair, pair_class) %>%
    group_modify(~ slice_sample(.x, n = min(nrow(.x), max_n))) %>% ungroup()
}

format_pairs <- function(x, cohort_label) {
  x %>%
    filter(group_scope == "All", distance_version == "primary") %>%
    thin_plot_data() %>%
    mutate(
      cohort = cohort_label,
      omics = omics_factor(omics),
      period_pair = factor(period_pair, levels = period_levels, labels = period_labels),
      pair_class = factor(pair_class, levels = pair_levels, labels = pair_labels)
    ) %>% droplevels()
}

format_pvalues <- function(x, cohort_label) {
  x %>%
    filter(group_scope == "All", distance_version == "primary") %>%
    mutate(
      cohort = cohort_label,
      omics = omics_factor(omics),
      period_pair = factor(period_pair, levels = period_levels, labels = period_labels),
      p_plot = coalesce(p_distance_fdr, p_distance),
      stars = p_stars(p_plot)
    )
}

star_positions <- function(data, pdat, group_vars = c("omics", "period_pair")) {
  data <- data %>% mutate(across(all_of(group_vars), as.character))
  pdat <- pdat %>% mutate(across(all_of(group_vars), as.character))
  ymax <- data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(y = max(distance, na.rm = TRUE), .groups = "drop")
  pdat %>%
    left_join(ymax, by = group_vars) %>%
    filter(stars != "", is.finite(y)) %>%
    mutate(y = y * 0.95)
}

make_omics_panel <- function(pair_data, p_data, omics_value,
                             show_y_title = FALSE) {
  pair_data <- pair_data %>% filter(as.character(omics) == omics_value)
  p_data <- p_data %>% filter(as.character(omics) == omics_value)
  ann <- star_positions(pair_data, p_data)
  ann$omics <- factor(ann$omics, levels = omics_display_levels)
  ann$period_pair <- factor(ann$period_pair, levels = period_labels)
  p <- ggplot(pair_data, aes(period_pair, distance, fill = pair_class, colour = pair_class)) +
    facet_wrap(. ~ unname(omics_facet_labels[omics_value])) +
    geom_violin(trim = TRUE, scale = "width", linewidth = 0.25, alpha = 0.25,
                position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.12,
                 linewidth = 0.25,
                 position = position_dodge(width = 0.8)) +
    geom_text(data = ann, aes(period_pair, y, label = stars), inherit.aes = FALSE,
              size = 5, vjust = 0, colour = "black") +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = pair_colours, drop = FALSE) +
    scale_colour_manual(values = pair_colours, drop = FALSE) +
    labs(x = NULL, y = "Distance") +
    coord_cartesian(clip = "off") +
    base_theme +
    theme(legend.position = "none",
          plot.title = element_text(size = 13, hjust = 0.5))
  if (omics_value != "Metabolome") {
    p <- p + scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, by = 0.25),
                                expand = expansion(mult = c(0.01, 0.04)))
  }
  if (!show_y_title) {
    p <- p + theme(axis.title.y = element_blank())
  }
  p
}

# Figure 2e: THSBC and HBC
thsbc_pairs <- format_pairs(thsbc_results$pairs, "THSBC")
hbc_pairs <- format_pairs(hbc_results$pairs, "HBC")
thsbc_p <- format_pvalues(thsbc_results$summary, "THSBC")
hbc_p <- format_pvalues(hbc_results$summary, "HBC")

thsbc_omics_panels <- lapply(omics_display_levels, function(z) {
  make_omics_panel(thsbc_pairs, thsbc_p, z,
                   show_y_title = z == "Metabolome")
})
hbc_omics_present <- omics_display_levels[
  omics_display_levels %in% unique(as.character(hbc_pairs$omics))
]
hbc_omics_panels <- lapply(hbc_omics_present, function(z) {
  make_omics_panel(hbc_pairs, hbc_p, z,
                   show_y_title = z == "Metabolome")
})

# Use the first panel in each row to carry the cohort label, while retaining
# the corresponding omics label as its subtitle.
thsbc_omics_panels[[1]] <- thsbc_omics_panels[[1]] +
  labs(title = "THSBC", tag = 'e') +
  theme(plot.title = element_text(size = 13, hjust = 0), 
          plot.tag.position = c(0.05, 0.97),
          plot.tag = element_text(size = 22, colour = "black"))
hbc_omics_panels[[1]] <- hbc_omics_panels[[1]] +
  labs(title = "HBC") +
  theme(plot.title = element_text(size = 13, hjust = 0))

# Fig. 2e, panel c: external cohort, two prespecified schemes 
external_three_pairs <- external_three$pairs %>%
  filter(group_scope == "All") %>%
  transmute(analysis_scheme = "Three target-period windows",
            x_label = factor(period_pair, levels = period_levels, labels = period_labels),
            pair_class = factor(pair_class, levels = pair_levels, labels = pair_labels),
            distance)
external_three_p <- external_three$summary %>%
  filter(group_scope == "All") %>%
  mutate(p_plot = p.adjust(p_distance, method = "fdr"),
         stars = p_stars(p_plot),
         analysis_scheme = "Three target-period windows",
         x_label = factor(period_pair, levels = period_levels, labels = period_labels))

external_hf_pairs <- external_high_frequency$subject %>%
  select(lag_group, within_distance, between_reference) %>%
  pivot_longer(c(within_distance, between_reference),
               names_to = "pair_class", values_to = "distance") %>%
  mutate(
    analysis_scheme = "All high-frequency pairs by lag",
    pair_class = recode(pair_class, within_distance = "within", between_reference = "between"),
    pair_class = factor(pair_class, levels = pair_levels, labels = pair_labels),
    x_label = factor(lag_group, levels = interval_levels, labels = interval_labels)
  )
external_hf_p <- external_high_frequency$summary %>%
  mutate(analysis_scheme = "All high-frequency pairs by lag",
         x_label = factor(lag_group, levels = interval_levels, labels = interval_labels),
         stars = p_stars(p_distance_fdr))

external_pairs <- bind_rows(external_three_pairs, external_hf_pairs)
external_three_ann <- external_three_p %>%
  select(x_label, stars) %>%
  left_join(external_three_pairs %>% group_by(x_label) %>%
              summarise(y = max(distance, na.rm = TRUE), .groups = "drop"),
            by = "x_label") %>%
  filter(stars != "", is.finite(y)) %>% mutate(y = y * 0.97)
external_hf_ann <- external_hf_p %>%
  select(x_label, stars) %>%
  left_join(external_hf_pairs %>% group_by(x_label) %>%
              summarise(y = max(distance, na.rm = TRUE), .groups = "drop"),
            by = "x_label") %>%
  filter(stars != "", is.finite(y)) %>% mutate(y = y * 0.97)

external_plot_theme <- theme(legend.position = "none",
                             axis.title.y = element_blank(),
                             plot.title = element_text(size = 13))
fig2e_external_three <- ggplot(external_three_pairs, aes(x_label, distance, fill = pair_class, colour = pair_class)) +
   facet_wrap(. ~ "Three target timepoints") +
  geom_violin(trim = TRUE, scale = "width", linewidth = 0.25, alpha = 0.25,
              position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.12, linewidth = 0.25,
               position = position_dodge(width = 0.8)) +
  geom_text(data = external_three_ann, aes(x_label, y, label = stars), inherit.aes = FALSE,
            size = 5, vjust = 0, colour = "black") +
  scale_fill_manual(values = pair_colours, drop = FALSE) +
  scale_colour_manual(values = pair_colours, drop = FALSE) +
  scale_y_continuous(limits = c(12, 33), breaks = c(15, 20, 25, 30),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Liang et al.",  x = NULL, y = "Distance") +
  coord_cartesian(clip = "off") + base_theme + external_plot_theme

fig2e_external_hf <- ggplot(external_hf_pairs, aes(x_label, distance, fill = pair_class, colour = pair_class)) +
  facet_wrap(. ~ "High-frequency pairs") + 
  geom_violin(trim = TRUE, scale = "width", linewidth = 0.25, alpha = 0.25,
              position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.12, linewidth = 0.25,
               position = position_dodge(width = 0.8)) +
  geom_text(data = external_hf_ann, aes(x_label, y, label = stars), inherit.aes = FALSE,
            size = 5, vjust = 0, colour = "black") +
  scale_fill_manual(values = pair_colours, drop = FALSE) +
  scale_colour_manual(values = pair_colours, drop = FALSE) +
  labs(title = "Liang et al.", x = "Gestational-week interval", y = "Distance") +
  coord_cartesian(clip = "off") + base_theme + external_plot_theme +
  theme(legend.position = c(0.86, 0.12), legend.direction = "vertical")

# Use six equal-width columns. The first row contains the five THSBC panels
# and the three-timepoint Liang panel. In the second row, the four HBC panels
# each occupy one column and the high-frequency Liang panel spans columns 5-6.
fig2e_items <- c(thsbc_omics_panels, list(fig2e_external_three),
                 hbc_omics_panels, list(fig2e_external_hf))
fig2e_design <- "
ABCDEF
GHIJKK
"
fig2e <- wrap_plots(fig2e_items, design = fig2e_design, guides = "keep")

# Save plotting data and figures 
write.csv(thsbc_pairs, file.path(result_dir, "fig2e_thsbc_within_between_distance_data.csv"), row.names = FALSE)
write.csv(hbc_pairs, file.path(result_dir, "fig2e_hbc_within_between_distance_data.csv"), row.names = FALSE)
write.csv(external_pairs, file.path(result_dir, "fig2e_external_distance_data.csv"), row.names = FALSE)
if (identical(Sys.getenv("SKIP_DISTANCE_FIGURES"), "1")) {
  quit(save = "no", status = 0)
}

message("saving Fig. 2e PDF")
ggsave(file.path(figure_dir, "fig2e_v2_multi-omics_distance_combined.pdf"),
       fig2e, width = 15, height = 6, device = grDevices::cairo_pdf)