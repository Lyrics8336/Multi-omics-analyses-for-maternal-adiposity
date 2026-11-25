####################################################
# Visualization and statistical testing for Figure 6
#####################################################

## Figure 6a: model performance ---------------------------------------------
rm(list = ls())

library(tidyverse)
library(patchwork)
library(cowplot)
library(scales)

model_dir <- "result/PART5/pbmi_model"
figure_dir <- "figure/PART5"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

predictions <- read_csv(
  file.path(model_dir, "oof_predictions.csv"), show_col_types = FALSE
)
fold_performance <- read_csv(
  file.path(model_dir, "fold_performance.csv"), show_col_types = FALSE
)

model_levels <- c(
  "Basic", "GenoBMI", "MicroBMI", "MetBMI", "Multi_omics_pBMI"
)
model_labels <- c(
  Basic = "Basic",
  GenoBMI = "Genomics",
  MicroBMI = "Microbiome",
  MetBMI = "Metabolomics",
  Multi_omics_pBMI = "Multi-omics"
)
model_colors <- c(
  Basic = "#0072B5", 
  GenoBMI = "#BC3C29",
  MicroBMI = "#E18727",
  MetBMI = "#20854E",
  Multi_omics_pBMI = "#7876B1"
)
performance_colors <- c(
  Basic = "#0072B5", 
  GenoBMI = "#BC3C29",
  MicroBMI = "#E18727",
  MetBMI = "#20854E",
  Multi_omics_pBMI = "#7876B1"
)

performance_plot_data <- fold_performance %>%
  mutate(model = factor(model, levels = model_levels))

performance_summary <- performance_plot_data %>%
  group_by(model) %>%
  summarise(
    mean_R2 = mean(R2),
    se_R2 = sd(R2) / sqrt(n()),
    lower_ci = mean_R2 - 1.96 * se_R2,
    upper_ci = mean_R2 + 1.96 * se_R2,
    .groups = "drop"
  )

p1 <- ggplot(performance_summary, aes(model, mean_R2, fill = model)) +
  geom_col(color = "black", width = 0.62) +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci), width = 0.18, linewidth = 0.45
  ) +
  geom_jitter(
    data = performance_plot_data,
    aes(y = R2),
    width = 0.08, height = 0, size = 1.6, alpha = 0.45,
    inherit.aes = TRUE
  ) +
  scale_fill_manual(values = performance_colors) +
  scale_x_discrete(labels = model_labels) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 0.6),breaks=seq(0, 0.6,0.1))+
  labs(x = NULL, y = expression("Out-of-sample " * R^2), tag = 'a') +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12),
    legend.position = "none", panel.grid = element_blank(),
    plot.tag = element_text(color = "black", size = 16),
    plot.tag.position = c(0.03, 0.98),
    plot.background = element_blank(), panel.background = element_blank()
  )

make_prediction_plot <- function(model, color, xlab = NULL, ylab = NULL, title = NULL) {
  prediction_var <- paste0(model, "_pred_pbmi")
  plot_data <- predictions %>%
    transmute(
      measured_pbmi = observed_pbmi,
      inferred_pbmi = .data[[prediction_var]]
    )
  correlation <- cor.test(plot_data$measured_pbmi, plot_data$inferred_pbmi)

  ggplot(plot_data, aes(inferred_pbmi, measured_pbmi)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.45) +
    geom_point(color = color, alpha = 0.2, size = 2, shape = 21, fill = "white", stroke = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = color, fill = alpha(color, 0.2)) +
    annotate(
      "text", x = 16, y = 34,
      label = paste0(
        "Pearson's r = ", round(unname(correlation$estimate), 3),
        "\nP = ", scientific(correlation$p.value, digits = 2)
      ),
      hjust = 0, vjust = 1, size = 4
    ) +
    scale_x_continuous(limits = c(16, 28), breaks=seq(16, 28, 4))+
    scale_y_continuous(limits = c(15, 34), breaks=seq(15, 35, 5))+
    labs( x = xlab, y = ylab, title = title ) +
    theme_classic(base_size = 10) +
    theme(plot.background = element_blank(),  
          panel.background = element_blank(),  
          legend.background = element_blank(),
          panel.grid = element_blank(),
          axis.title.x = element_text(size = 13, hjust = 0.5, vjust = 15),
          axis.title.y = element_text(size = 12, hjust = 0.5),
          plot.title = element_text(size = 14, hjust = 0.5),
          axis.text = element_text(size = 11, colour = 'black'),
          legend.position = 'none')
}

p2 <- make_prediction_plot(
  "Basic", model_colors[["Basic"]],
  ylab = expression("Observed pBMI (kg/m"^2 * ")"), title = "Basic"
)
p3 <- make_prediction_plot(
  "GenoBMI", model_colors[["GenoBMI"]], title = "Genomics"
)
p4 <- make_prediction_plot(
  "MicroBMI", model_colors[["MicroBMI"]],
  xlab = expression("Inferred pBMI (kg/m"^2 * ")"), title = "Microbiome"
)
p5 <- make_prediction_plot(
  "MetBMI", model_colors[["MetBMI"]], title = "Metabolomics"
)
p6 <- make_prediction_plot(
  "Multi_omics_pBMI", model_colors[["Multi_omics_pBMI"]], title = "Multi-omics"
)

fig6a <- (p1 + p2 + p3 + p4 + p5 + p6) +
  plot_layout(ncol = 6) &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(
  file.path(figure_dir, "fig6a_pbmi_model_performance.pdf"),
  fig6a, width = 13.5, height = 4, bg = "transparent"
)

## Figure 6b: PGM outcome radar plots -----------------
rm(list = ls())

library(tidyverse)
library(ggradar)
library(patchwork)

analysis_file <- "result/PART5/PGM/PGM_outcome_analysis.rds"
figure_dir <- "figure/PART5"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

pgm_results <- readRDS(analysis_file)
occurrence <- pgm_results$outcome_occurrence
group_summary <- pgm_results$group_summary

outcome_order <- c(
  "GDM", "HDP", "PE", "ICP", "LBW", "SGA", "MAC", "LGA", "PTB", "owo_m36"
)
outcome_labels <- c(
  GDM = "GDM", HDP = "HDP", PE = "PE", ICP = "ICP", LBW = "LBW",
  SGA = "SGA", MAC = "MAC", LGA = "LGA", PTB = "PTB",
  owo_m36 = "Offspring OWO\nm36"
)
group_order <- c("Concordant", "Discordantly high", "Discordantly low")

prepare_ggradar_data <- function(pgm_type_value) {
  sample_sizes <- group_summary %>%
    filter(pgm_type == pgm_type_value) %>%
    mutate(pgm_group = factor(pgm_group, levels = group_order)) %>%
    arrange(pgm_group)

  legend_labels <- setNames(
    paste0(as.character(sample_sizes$pgm_group), " (n=", sample_sizes$n, ")"),
    as.character(sample_sizes$pgm_group)
  )

  radar_data <- occurrence %>%
    filter(pgm_type == pgm_type_value) %>%
    select(pgm_group, outcome, occurrence_percent) %>%
    mutate(
      pgm_group = factor(pgm_group, levels = group_order),
      outcome = factor(outcome, levels = outcome_order)
    ) %>%
    arrange(pgm_group, outcome) %>%
    pivot_wider(
      names_from = outcome,
      values_from = occurrence_percent
    ) %>%
    arrange(pgm_group) %>%
    mutate(
      group = unname(legend_labels[as.character(pgm_group)]),
      .before = 1
    ) %>%
    select(group, all_of(outcome_order))

  names(radar_data)[-1] <- unname(outcome_labels[outcome_order])
  list(data = radar_data, legend_labels = legend_labels)
}

make_ggradar_plot <- function(prepared_data, pgm_type_value, title,
                              grid_mid, grid_max, colours) {
  star_data <- occurrence %>%
    filter(
      pgm_type == pgm_type_value,
      pgm_group != "Concordant",
      significance != ""
    ) %>%
    mutate(
      axis_index = match(outcome, outcome_order),
      angle = 2 * pi * (axis_index - 1) / length(outcome_order) +
        if_else(pgm_group == "Discordantly high", -0.018, 0.018),
      radius = pmin(occurrence_percent + grid_max * 0.055, grid_max * 1.06),
      x = radius * sin(angle),
      y = radius * cos(angle),
      group = unname(prepared_data$legend_labels[as.character(pgm_group)])
    )

  radar_plot <- ggradar(
    prepared_data$data,
    background.circle.colour = "transparent",
    values.radar = c("0", paste0(grid_mid, "%"), paste0(grid_max, "%")),
    grid.min = 0,
    grid.mid = grid_mid,
    grid.max = grid_max,
    group.line.width = 0.8,
    group.point.size = 2,
    group.colours = colours,
    fill = TRUE,
    fill.alpha = 0.2,
    plot.title = title,
    grid.label.size = 4,
    gridline.label.offset = 0.1,
    axis.label.size = 5,
    axis.label.offset = 1.05,
    axis.line.colour = "grey",
    grid.line.width = 0.4,
    legend.text.size = 10,
    legend.position = c(0.95, 0.2),
    gridline.min.linetype = "solid",
    gridline.max.linetype = "solid",
    gridline.max.colour = "black"
  )

  radar_plot +
    geom_text(
      data = star_data,
      aes(x = x, y = y, label = significance, colour = group),
      inherit.aes = FALSE,
      size = 8,
      show.legend = FALSE
    ) +
    theme(
      plot.background = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5, vjust = -5),
      plot.title.position = "plot",
      plot.margin = margin(0, 0, 0, 0)
    )
}

measured_radar_data <- prepare_ggradar_data("Measured pBMI vs genetic BMI")
omics_radar_data <- prepare_ggradar_data(
  "Multi-omics-inferred pBMI vs genetic BMI"
)

p_measured <- make_ggradar_plot(
  measured_radar_data,
  pgm_type_value = "Measured pBMI vs genetic BMI",
  title = "PGM derived from observed pBMI",
  grid_mid = 10,
  grid_max = 20,
  colours = c("grey", "#F5A3A3", "#A4D7F3")
)

p_omics <- make_ggradar_plot(
  omics_radar_data,
  pgm_type_value = "Multi-omics-inferred pBMI vs genetic BMI",
  title = "PGM derived from multi-omics-inferred pBMI",
  grid_mid = 15,
  grid_max = 30,
  colours = c("grey", "#FCBB71", "#8FAADC")
)

fig6b <- p_measured / p_omics + plot_layout(heights = c(1, 1)) &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(
  file.path(figure_dir, "fig6b_PGM_outcome_radar.pdf"),
  fig6b, width = 8, height = 10, bg = "transparent"
)



## Figure 6c: health outcomes by observed and inferred pBMI -----------------
rm(list = ls())

library(tidyverse)
library(patchwork)
library(cowplot)
library(scales)

model_dir <- "result/PART5/pbmi_model"
figure_dir <- "figure/PART5"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

predictions <- read_csv(
  file.path(model_dir, "oof_predictions.csv"), show_col_types = FALSE
)
fold_performance <- read_csv(
  file.path(model_dir, "fold_performance.csv"), show_col_types = FALSE
)

load("data/01_thsbc_cleaned_data/total_pop_5212.RData")
load("data/01_thsbc_cleaned_data/outcomes_5212.RData")

outcome_vars <- c("GDM", "HDP", "PE", "ICP", "LBW", "MAC", "PTB", "SGA", "LGA", "owo_m36")
outcome_labels <- c(
  GDM = "GDM", HDP = "HDP", PE = "PE", ICP = "ICP",
  LBW = "LBW", MAC = "MAC", PTB = "PTB", SGA = "SGA", LGA = "LGA",
  owo_m36 = "Offspring OWO\nm36"
)
outcome_levels <- rev(unname(outcome_labels))

a_comb <- predictions %>%
  select(id, inferred_pbmi = Multi_omics_pBMI_pred_pbmi) %>%
  inner_join(total_pop %>% select(id, measured_group = pbmi_cut3), by = "id") %>%
  inner_join(otc %>% select(id, all_of(outcome_vars)), by = "id") %>%
  mutate(
    measured_owo = measured_group == "OWO",
    inferred_owo = inferred_pbmi >= 24,
    popclass = case_when(
      !measured_owo & !inferred_owo ~ "N_N",
      !measured_owo & inferred_owo ~ "N_O",
      measured_owo & !inferred_owo ~ "O_N",
      measured_owo & inferred_owo ~ "O_O"
    ),
    popclass = factor(popclass, levels = c("N_N", "N_O", "O_N", "O_O"))
  )

population_counts <- a_comb %>% count(popclass, name = "n")

outcome_plot_data <- a_comb %>%
  select(popclass, all_of(outcome_vars)) %>%
  pivot_longer(-popclass, names_to = "outcome_var", values_to = "outcome_value") %>%
  mutate(event = as.character(outcome_value) == "1") %>%
  group_by(popclass, outcome_var) %>%
  summarise(
    n_available = sum(!is.na(event)),
    n_event = sum(event, na.rm = TRUE),
    Frequency = mean(event, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    outcome = recode(outcome_var, !!!outcome_labels),
    outcome = factor(outcome, levels = outcome_levels)
  )

compare_outcomes_by_group <- function(data, classes, measured_stratum) {
  map_dfr(outcome_vars, function(outcome) {
    test_data <- data %>%
      filter(popclass %in% classes, !is.na(.data[[outcome]])) %>%
      mutate(
        popclass = droplevels(popclass),
        event = factor(as.character(.data[[outcome]]), levels = c("0", "1"))
      )
    contingency_table <- table(test_data$popclass, test_data$event)
    valid_table <- all(dim(contingency_table) == c(2, 2))

    if (valid_table) {
      chi_check <- suppressWarnings(chisq.test(contingency_table))
      use_chisq <- all(chi_check$expected >= 5)
      test_result <- if (use_chisq) {
        chisq.test(contingency_table)
      } else {
        fisher.test(contingency_table)
      }
      test_method <- if (use_chisq) {
        "Pearson's Chi-squared test"
      } else {
        "Fisher's exact test"
      }
      statistic <- if (use_chisq) unname(test_result$statistic) else NA_real_
      p_value <- test_result$p.value
    } else {
      test_method <- NA_character_
      statistic <- NA_real_
      p_value <- NA_real_
    }

    tibble(
      measured_stratum = measured_stratum,
      comparison = paste(classes, collapse = " vs "),
      outcome_var = outcome,
      outcome = unname(outcome_labels[[outcome]]),
      test_method = test_method,
      statistic = statistic,
      p_value = p_value,
      significance = case_when(
        is.na(p_value) ~ NA_character_,
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  })
}

comparison_results <- bind_rows(
  compare_outcomes_by_group(a_comb, c("N_N", "N_O"), "Non-OWO"),
  compare_outcomes_by_group(a_comb, c("O_O", "O_N"), "OWO")
) %>%
  mutate(outcome = factor(outcome, levels = outcome_levels))

make_significance_data <- function(plot_data, test_data, side) {
  plot_data %>%
    group_by(outcome) %>%
    summarise(max_frequency = max(Frequency), .groups = "drop") %>%
    left_join(test_data, by = "outcome") %>%
    filter(significance != "ns", !is.na(significance)) %>%
    mutate(
      y = as.numeric(outcome),
      x = if (side == "right") max_frequency + 3 else -(max_frequency + 3),
      tick_end = if (side == "right") x - 1.5 else x + 1.5,
      label_x = if (side == "right") x + 1.2 else x - 1.2
    )
}

right_data <- outcome_plot_data %>%
  filter(popclass %in% c("N_N", "N_O")) %>%
  mutate(popclass = factor(popclass, levels = c("N_N", "N_O")))
left_data <- outcome_plot_data %>%
  filter(popclass %in% c("O_O", "O_N")) %>%
  mutate(popclass = factor(popclass, levels = c("O_O", "O_N")))

right_significance <- make_significance_data(
  right_data, comparison_results %>% filter(measured_stratum == "Non-OWO"), "right"
)
left_significance <- make_significance_data(
  left_data, comparison_results %>% filter(measured_stratum == "OWO"), "left"
)

rightplot <- ggplot(right_data, aes(Frequency, outcome, fill = popclass)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.76, color = "white") +
  geom_segment(
    data = right_significance,
    aes(x = x, xend = x, y = y - 0.20, yend = y + 0.20),
    inherit.aes = FALSE, linewidth = 0.45
  ) +
  geom_segment(
    data = right_significance,
    aes(x = x, xend = tick_end, y = y - 0.20, yend = y - 0.20),
    inherit.aes = FALSE, linewidth = 0.45
  ) +
  geom_segment(
    data = right_significance,
    aes(x = x, xend = tick_end, y = y + 0.20, yend = y + 0.20),
    inherit.aes = FALSE, linewidth = 0.45
  ) +
  geom_text(
    data = right_significance,
    aes(x = label_x, y = y, label = significance),
    inherit.aes = FALSE, hjust = 0, size = 6, vjust = 0.7
  ) +
  scale_fill_manual(
    values = c(N_N = "#A1C6E7", N_O = "#B9D9EB"),
    labels = c(
      N_N = "Observed non-OWO & inferred non-OWO",
      N_O = "Observed non-OWO & inferred OWO"
    ),
    name = NULL
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 60), breaks = seq(0, 60, 15)) +
  scale_y_discrete(expand = expansion(add = 0.55)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 10, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(r = 8, b = 24),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

leftplot <- ggplot(left_data, aes(-Frequency, outcome, fill = popclass)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.76, color = "white") +
  geom_segment(
    data = left_significance,
    aes(x = x, xend = x, y = y - 0.20, yend = y + 0.20),
    inherit.aes = FALSE, linewidth = 0.45
  ) +
  geom_segment(
    data = left_significance,
    aes(x = x, xend = tick_end, y = y - 0.20, yend = y - 0.20),
    inherit.aes = FALSE, linewidth = 0.45
  ) +
  geom_segment(
    data = left_significance,
    aes(x = x, xend = tick_end, y = y + 0.20, yend = y + 0.20),
    inherit.aes = FALSE, linewidth = 0.45
  ) +
  geom_text(
    data = left_significance,
    aes(x = label_x, y = y, label = significance),
    inherit.aes = FALSE, hjust = 1, size = 6, vjust = 0.7
  ) +
  scale_fill_manual(
    values = c(O_O = "#F4A261", O_N = "#F8C3B1"),
    labels = c(
      O_O = "Observed OWO & inferred OWO",
      O_N = "Observed OWO & inferred non-OWO"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    expand = c(0, 0), limits = c(-60, 0), breaks = seq(-60, 0, 15),
    labels = abs
  ) +
  scale_y_discrete(position = "right", expand = expansion(add = 0.55)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(l = 8, b = 24),
    plot.background = element_blank()
  )

group_colors <- c(
  O_N = "#F8C3B1", O_O = "#F4A261",
  N_O = "#B9D9EB", N_N = "#A1C6E7"
)
group_labels <- c(
  O_N = "Observed OWO & inferred non-OWO",
  O_O = "Observed OWO & inferred OWO",
  N_O = "Observed non-OWO & inferred OWO",
  N_N = "Observed non-OWO & inferred non-OWO"
  
)

legend_plot <- ggplot(
  tibble(
    popclass = factor(names(group_colors), levels = names(group_colors)),
    x = seq_along(group_colors), y = 1
  ),
  aes(x, y, fill = popclass)
) +
  geom_col() +
  scale_fill_manual(values = group_colors, labels = group_labels, name = NULL) +
  guides(fill = guide_legend(ncol = 2, byrow = FALSE)) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.key.size = unit(11, "pt"),
    legend.spacing.x = unit(5, "pt"),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )
outcome_legend <- get_legend(legend_plot)

combined_plot <- ggdraw() +
  draw_plot(leftplot, x = 0, y = 0, width = 0.42, height = 0.90) +
  draw_plot(rightplot, x = 0.42, y = 0, width = 0.54, height = 0.90) +
  draw_plot(outcome_legend, x = 0.18, y = 0.89, width = 0.72, height = 0.10) +
  draw_label(
    "OWO", x = 0.1, y = 0.47, angle = 90,
    hjust = 0.5, vjust = 0.5, size = 16
  ) +
  draw_label(
    "Non-OWO", x = 0.9, y = 0.47, angle = 270,
    hjust = 0.5, vjust = 0.5, size = 16
  ) +
  draw_label(
    "Incidence of adverse outcome (%)",
    x = 0.5, y = 0.03, hjust = 0.5, vjust = 0.5, size = 13
  )

fig6c <- combined_plot
ggsave(
  file.path(figure_dir, "fig6c_pbmi_outcome_comparison.pdf"),
  fig6c, width = 8, height = 7, bg = "transparent"
)

