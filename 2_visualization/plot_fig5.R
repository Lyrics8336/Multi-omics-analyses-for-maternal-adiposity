####################################################
# Visualization and statistical testing for Figure 5
#####################################################

# Fig 5a: associations between pBMI and pregnancy/offspring outcomes ---------------------------------------------
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(patchwork)

load("result/PART4/1_trad_pbmi_outcomes.RData")

pregnancy_complications <- c("GDM", "HDP", "PE", "ICP")
birth_outcomes <- c(
  "LBW", "MAC", "PTB", "SGA", "LGA",
  "Delivery_week", "Birth_weight", "Birth_length"
)
child_outcomes <- c("owo_m36", "zbmi_m36")

outcome_labels <- c(
  GDM = "GDM", HDP = "HDP", PE = "PE", ICP = "ICP",
  LBW = "LBW", MAC = "MAC", PTB = "PTB", SGA = "SGA", LGA = "LGA",
  Delivery_week = "Delivery\nweek", Birth_weight = "Birth\nweight",
  Birth_length = "Birth\nlength", owo_m36 = "Offspring \nOWO(m36)   ",
  zbmi_m36 = "Offspring \nBMI(m36)   "
)


pdata_fig5a <- result_list$A1_pbmi_outcomes %>%
  filter(exposure == "pbmi_owo_normal") %>%
  mutate(
    pv_sig = case_when(
      p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"
    ),
    outcome_type = if_else(model_type == "logistic", "Binary", "Continuous")
  )

pdata_fig5a_preg <- pdata_fig5a %>%
  filter(outcome %in% pregnancy_complications) %>%
  mutate(outcome = factor(outcome, levels = pregnancy_complications))

pdata_fig5a_birth <- pdata_fig5a %>%
  filter(outcome %in% birth_outcomes) %>%
  mutate(outcome = factor(outcome, levels = birth_outcomes))

pdata_fig5a_child <- pdata_fig5a %>%
  filter(outcome %in% child_outcomes) %>%
  mutate(outcome = factor(outcome, levels = child_outcomes))

common_theme <- theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank(),
    panel.grid.major = element_line(linetype = 2, linewidth = 0.35),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 13, color = "black"),
    strip.background = element_rect(fill = "#faecb1"),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  )

## Horizontal forest plots 
# Continuous-outcome coefficients are plotted on the primary coordinate as
# 1 + beta_scale * beta. Hence OR = 1 and beta = 0 share one reference line.
beta_scale <- 5
to_primary <- function(beta) 1 + beta_scale * beta
to_beta <- function(primary) (primary - 1) / beta_scale
beta_breaks <- seq(-0.2, 0.6, 0.2)

prepare_dual_axis_data <- function(data) {
  data %>%
    mutate(
      plot_value = if_else(model_type == "logistic", effect, to_primary(beta)),
      plot_lci = if_else(model_type == "logistic", effect_lci, to_primary(lci)),
      plot_uci = if_else(model_type == "logistic", effect_uci, to_primary(uci))
    )
}

plot_horizontal_outcomes <- function(
    data, panel_title,
    show_bottom_x_axis = FALSE, show_top_x_axis = FALSE,
    use_expand = TRUE) 
    {
  plot_data <- prepare_dual_axis_data(data) %>%
    mutate(panel = panel_title)
  outcome_order <- rev(levels(droplevels(plot_data$outcome)))

  ggplot(
    plot_data,
    aes(
      x = plot_value, y = outcome,
      xmin = plot_lci, xmax = plot_uci,
      color = outcome_type
    )
  ) +
    facet_wrap(~ panel, strip.position = "left") +
    # One reference line represents both OR = 1 and beta = 0.
    geom_vline(xintercept = 1, color = "black", linewidth = 0.8, linetype = 3) +
    geom_pointrange(orientation = "y", linewidth = 0.7, size = 1) +
    geom_text(
      data = subset(plot_data, pv_sig != "ns"),
      aes(x = plot_value, label = pv_sig),
      nudge_y = 0.15,
      size = 5, color = "black", show.legend = FALSE
    ) +
    scale_color_manual(values = c(Binary = "#ff7f50", Continuous = "#ff7f50")) +
    scale_y_discrete(
      limits = outcome_order,
      labels = outcome_labels,
      expand = if (use_expand) {
        expansion(add = c(0.15, 0.45))
      } else {
        expansion(mult = c(0.5, 0.5))
      }
    ) +
    scale_x_continuous(
      name = if (show_bottom_x_axis) "Beta" else NULL,
      limits = c(0, 4.3),
      breaks = to_primary(beta_breaks),
      labels = function(x) scales::number(to_beta(x), accuracy = 0.1),
      expand = expansion(mult = c(0, 0.02)),
      sec.axis = sec_axis(
        ~ .,
        name = if (show_top_x_axis) "Odds ratio" else NULL,
        breaks = 0:4
      )
    ) +
    labs(y = NULL) +
    common_theme +
    theme(
      axis.title.x.bottom = if (show_bottom_x_axis) {
        element_text(color = "black")
      } else {
        element_blank()
      },
      axis.text.x.bottom = if (show_bottom_x_axis) {
        element_text(color = "black")
      } else {
        element_blank()
      },
      axis.ticks.x.bottom = if (show_bottom_x_axis) {
        element_line(color = "black")
      } else {
        element_blank()
      },
      axis.title.x.top = if (show_top_x_axis) {
        element_text(color = "black")
      } else {
        element_blank()
      },
      axis.text.x.top = if (show_top_x_axis) {
        element_text(color = "black")
      } else {
        element_blank()
      },
      axis.ticks.x.top = if (show_top_x_axis) {
        element_line(color = "black")
      } else {
        element_blank()
      }
    )
}

## First column: pregnancy complications above child outcomes 
p1 <- plot_horizontal_outcomes(
  pdata_fig5a_preg, "Pregnancy complications",
  show_bottom_x_axis = FALSE, show_top_x_axis = TRUE,
  use_expand = TRUE
) +
  labs(tag = "a") +
  theme(
    plot.tag = element_text(size = 20),
    plot.tag.position = c(0.05, 0.98)
  )

p3 <- plot_horizontal_outcomes(
  pdata_fig5a_child, "Child outcomes",
  show_bottom_x_axis = TRUE, show_top_x_axis = FALSE,
  use_expand = FALSE
)

first_column <- p1 / p3 +
  plot_layout(heights = c(length(pregnancy_complications), length(child_outcomes)))

## Second column: birth outcomes 
p2 <- plot_horizontal_outcomes(
  pdata_fig5a_birth, "Birth outcomes",
  show_bottom_x_axis = TRUE, show_top_x_axis = TRUE,
  use_expand = TRUE
)

## Combine 
fig5a <- (first_column | p2) +
  plot_layout(widths = c(
    length(pregnancy_complications) + length(child_outcomes),
    length(birth_outcomes)
  ))


# Fig 5b: elastic net prediction of GDM and HDP ------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(pROC)
  library(patchwork)
})

prediction_dir <- "result/PART4/prediction/3.6.1elastic_net/"

prediction_file <- file.path(
  prediction_dir, "3.6.1_elastic_net_model_specific_lambda_GDM_HDP.rds"
)
if (!file.exists(prediction_file)) {
  stop("Missing prediction result: ", prediction_file)
}

prediction_model_levels <- c(
  "B0", "B1", "B1_Microbiome", "B1_Metabolome", "B1_Multiomics"
)
prediction_model_codes <- c(
  B0 = "M0", B1 = "M1", B1_Microbiome = "M2",
  B1_Metabolome = "M3", B1_Multiomics = "M4"
)
prediction_model_labels <- c(
  B0 = "M0: Clinical",
  B1 = "M1: Clinical + pBMI",
  B1_Microbiome = "M2: M1 + microbiome",
  B1_Metabolome = "M3: M1 + metabolome",
  B1_Multiomics = "M4: M1 + multiomics"
)
prediction_model_colours <- c(
  B0 = "#0072B5", B1 = "#E18727", B1_Microbiome = "#BC3C29",
  B1_Metabolome = "#20854E", B1_Multiomics = "#7876B1"
)
prediction_window_labels <- c(
  T1 = "T1", T2 = "T2", Mean_T1_T2 = "Pooled"
)
prediction_outcomes <- c("GDM", "HDP")

prediction_results <- readRDS(prediction_file)
prediction_metrics <- prediction_results$metrics
prediction_values <- prediction_results$predictions

required_metric_columns <- c(
  "population", "outcome", "window", "model", "repeat_id", "ROC"
)
missing_metric_columns <- setdiff(required_metric_columns, names(prediction_metrics))
if (length(missing_metric_columns)) {
  stop(
    basename(prediction_file), " lacks metric columns: ",
    paste(missing_metric_columns, collapse = ", ")
  )
}

required_prediction_columns <- c(
  "id", "population", "outcome", "window", "model", "repeat_id",
  "outer_fold", "observed", "probability"
)
missing_prediction_columns <- setdiff(
  required_prediction_columns, names(prediction_values)
)
if (length(missing_prediction_columns)) {
  stop(
    basename(prediction_file), " lacks prediction columns: ",
    paste(missing_prediction_columns, collapse = ", ")
  )
}

if (!setequal(unique(prediction_metrics$outcome), prediction_outcomes)) {
  stop("The prediction result must contain GDM and HDP only.")
}
if (!identical(unique(prediction_metrics$population), "available_case")) {
  stop("Expected available-case prediction results only.")
}
if (any(!is.finite(prediction_metrics$ROC))) {
  stop("The prediction metrics contain missing or non-finite AUROC values.")
}
if (any(!prediction_metrics$model %in% prediction_model_levels) ||
    any(!prediction_metrics$window %in% names(prediction_window_labels))) {
  stop("Unexpected model or period label in the prediction results.")
}
if (anyDuplicated(prediction_metrics[
  c("population", "outcome", "window", "model", "repeat_id")
])) {
  stop("Duplicate AUROC result for an outcome/period/model/repeat.")
}

## Mean AUROC across the 10 outer-CV repeats 
prediction_auc_summary <- prediction_metrics %>%
  group_by(.data$population, .data$outcome, .data$window, .data$model) %>%
  summarise(
    mean_auc = mean(.data$ROC),
    sd_auc = sd(.data$ROC),
    n_repeats = n_distinct(.data$repeat_id),
    .groups = "drop"
  )
if (any(prediction_auc_summary$n_repeats != 10L)) {
  stop("Expected 10 repeat-level AUROC values in every GDM/HDP tile.")
}

expected_prediction_grid <- tidyr::expand_grid(
  outcome = prediction_outcomes,
  window = names(prediction_window_labels),
  model = prediction_model_levels
)
actual_prediction_grid <- prediction_auc_summary %>%
  distinct(.data$outcome, .data$window, .data$model)
if (nrow(anti_join(
  expected_prediction_grid, actual_prediction_grid,
  by = c("outcome", "window", "model")
)) > 0L ||
    nrow(anti_join(
      actual_prediction_grid, expected_prediction_grid,
      by = c("outcome", "window", "model")
    )) > 0L) {
  stop("The prediction result lacks the expected GDM/HDP period/model grid.")
}

## Pooled-window ROC curves from out-of-fold predictions -
pooled_predictions <- prediction_values %>%
  filter(
    .data$outcome %in% prediction_outcomes,
    .data$window == "Mean_T1_T2",
    .data$population == "available_case"
  )

expected_roc_grid <- tidyr::expand_grid(
  outcome = prediction_outcomes,
  model = prediction_model_levels,
  repeat_id = seq_len(10L)
)
actual_roc_grid <- pooled_predictions %>%
  distinct(.data$outcome, .data$model, .data$repeat_id)
if (nrow(anti_join(
  expected_roc_grid, actual_roc_grid,
  by = c("outcome", "model", "repeat_id")
)) > 0L ||
    nrow(anti_join(
      actual_roc_grid, expected_roc_grid,
      by = c("outcome", "model", "repeat_id")
    )) > 0L) {
  stop("Pooled GDM/HDP predictions lack the expected model/repeat grid.")
}
if (anyDuplicated(pooled_predictions[
  c("outcome", "model", "repeat_id", "id")
])) {
  stop("A participant has multiple pooled predictions within one repeat.")
}
if (any(!is.finite(pooled_predictions$probability)) ||
    any(!pooled_predictions$observed %in% c("x0", "x1"))) {
  stop("Invalid pooled GDM/HDP outcomes or predicted probabilities.")
}

pooled_group_sizes <- pooled_predictions %>%
  group_by(.data$outcome, .data$model, .data$repeat_id) %>%
  summarise(
    n = n(),
    n_classes = n_distinct(.data$observed),
    .groups = "drop"
  ) %>%
  left_join(
    prediction_metrics %>%
      filter(
        .data$window == "Mean_T1_T2",
        .data$outcome %in% prediction_outcomes
      ) %>%
      select(outcome, model, repeat_id, sample_size),
    by = c("outcome", "model", "repeat_id")
  )
if (anyNA(pooled_group_sizes$sample_size) ||
    any(pooled_group_sizes$n_classes != 2L) ||
    any(pooled_group_sizes$n != pooled_group_sizes$sample_size)) {
  stop("A pooled ROC repeat does not cover its full analysis cohort.")
}

fpr_grid <- seq(0, 1, length.out = 201L)
roc_curve_repeats <- pooled_predictions %>%
  group_by(.data$outcome, .data$model, .data$repeat_id) %>%
  group_modify(~ {
    roc_object <- pROC::roc(
      response = .x$observed,
      predictor = .x$probability,
      levels = c("x0", "x1"),
      direction = "<",
      quiet = TRUE
    )
    sensitivity <- as.numeric(pROC::coords(
      roc_object,
      x = 1 - fpr_grid,
      input = "specificity",
      ret = "sensitivity",
      transpose = FALSE
    )$sensitivity)
    sensitivity[c(1L, length(sensitivity))] <- c(0, 1)

    tibble(
      fpr = fpr_grid,
      sensitivity = sensitivity,
      calculated_auc = as.numeric(pROC::auc(roc_object))
    )
  }) %>%
  ungroup()

roc_auc_check <- roc_curve_repeats %>%
  distinct(
    .data$outcome, .data$model, .data$repeat_id, .data$calculated_auc
  ) %>%
  left_join(
    prediction_metrics %>%
      filter(
        .data$window == "Mean_T1_T2",
        .data$outcome %in% prediction_outcomes
      ) %>%
      select(outcome, model, repeat_id, ROC),
    by = c("outcome", "model", "repeat_id")
  )
if (anyNA(roc_auc_check$ROC) ||
    any(abs(roc_auc_check$calculated_auc - roc_auc_check$ROC) > 1e-6)) {
  stop("Recalculated pooled ROC AUC does not match the saved metrics.")
}

roc_curve_mean <- roc_curve_repeats %>%
  group_by(.data$outcome, .data$model, .data$fpr) %>%
  summarise(
    mean_sensitivity = mean(.data$sensitivity),
    .groups = "drop"
  ) %>%
  mutate(
    outcome = factor(.data$outcome, levels = prediction_outcomes),
    model = factor(.data$model, levels = prediction_model_levels)
  )

make_roc_plot <- function(outcome_name, show_legend = FALSE, tag = NULL) {
  plot_data <- roc_curve_mean %>%
    filter(.data$outcome == outcome_name)

  plot <- ggplot(
    plot_data,
    aes(
      x = .data$fpr,
      y = .data$mean_sensitivity,
      colour = .data$model
    )
  ) +
    geom_abline(
      intercept = 0, slope = 1, linetype = "dashed",
      colour = "grey55", linewidth = 0.45
    ) +
    geom_line(linewidth = 0.85, alpha = 0.95) +
    facet_wrap(~ outcome, labeller = labeller(outcome = outcome_labels)) +
    scale_colour_manual(
      values = prediction_model_colours,
      breaks = prediction_model_levels,
      labels = unname(prediction_model_labels[prediction_model_levels]),
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      expand = expansion(mult = 0)
    ) +
    guides(colour = guide_legend(
      ncol = 1, byrow = TRUE,
      override.aes = list(linewidth = 1.1)
    )) +
    labs(x = "1 - Specificity", y = "Sensitivity", tag = tag) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_line(linetype = "dotted", colour = "grey93"),
      axis.text = element_text(size = 11, colour = "black"),
      axis.title = element_text(size = 12, colour = "black"),
      strip.background = element_rect(fill = "#faecb1"),
      strip.text = element_text(size = 13, colour = "black"),
      plot.tag = element_text(size = 20),
      plot.tag.position = c(0.02, 0.99)
    )

  if (show_legend) {
    plot +
      theme(
        legend.position = "inside",
        legend.background = element_blank(),
        legend.position.inside = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.margin = margin(3, 5, 3, 4),
        legend.text = element_text(size = 9),
        legend.key.height = grid::unit(0.8, "lines"),
        legend.key.width = grid::unit(1.2, "lines")
      )
  } else {
    plot + theme(legend.position = "none")
  }
}

## Mean AUROC (SD) heatmaps 
auc_plot_data <- prediction_auc_summary %>%
  mutate(
    period = unname(prediction_window_labels[.data$window]),
    tile_label = sprintf("%.3f\n(%.3f)", .data$mean_auc, .data$sd_auc)
  )

auc_fill_limits <- range(auc_plot_data$mean_auc, finite = TRUE)
auc_fill_breaks <- c(0.5, 0.6, 0.7)
auc_label_cutoff <- auc_fill_limits[[1]] + 0.72 * diff(auc_fill_limits)

auc_plot_data <- auc_plot_data %>%
  mutate(label_colour = if_else(
    .data$mean_auc >= auc_label_cutoff, "white", "#18232b"
  ))

make_auc_heatmap <- function(outcome_name, show_legend = FALSE) {
  plot_data <- auc_plot_data %>%
    filter(.data$outcome == outcome_name) %>%
    mutate(
      outcome = factor(.data$outcome, levels = prediction_outcomes),
      period = factor(.data$period, levels = c("T1", "T2", "Pooled"))
    )

  plot <- ggplot(
    plot_data,
    aes(x = .data$period, y = .data$model, fill = .data$mean_auc)
  ) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(
      aes(label = .data$tile_label, colour = .data$label_colour),
      size = 4, lineheight = 0.95, show.legend = FALSE
    ) +
    facet_wrap(~ outcome, labeller = labeller(outcome = outcome_labels)) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = 0)) +
    scale_y_discrete(
      limits = rev(prediction_model_levels),
      drop = FALSE,
      labels = function(x) unname(prediction_model_codes[x]),
      expand = expansion(add = 0)
    ) +
    scale_fill_gradientn(
      colours = c(
        "#F2F7FA", "#D4E4ED", "#AFC9D9", "#7DA8C3",
        "#4A7B9D", "#28506B", "#0E2F40"
      ),
      limits = auc_fill_limits,
      breaks = auc_fill_breaks,
      labels = scales::label_number(accuracy = 0.1),
      oob = scales::squish,
      name = "Mean AUC"
    ) +
    scale_colour_identity() +
    labs(x = "Timepoint", y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "#faecb1"),
      strip.text = element_text(size = 13, colour = "black"),
      axis.text = element_text(colour = "black", size = 11),
      axis.ticks = element_blank(),
      legend.position = if (show_legend) "bottom" else "none"
    )

  plot
}

## Stack GDM above HDP within the ROC and heatmap columns 
roc_gdm <- make_roc_plot("GDM", show_legend = FALSE, tag = "b")
heatmap_gdm <- make_auc_heatmap("GDM", show_legend = FALSE)

roc_hdp <- make_roc_plot("HDP", show_legend = TRUE)
heatmap_hdp <- make_auc_heatmap("HDP", show_legend = TRUE)

roc_column <- (roc_gdm / roc_hdp) +
  plot_layout(
    heights = c(1, 1.12),
    axes = "collect_x",
    axis_titles = "collect"
  )

heatmap_column <- (heatmap_gdm / heatmap_hdp) +
  plot_layout(
    heights = c(1, 1.12),
    axes = "collect_x"
  )

fig5b <- (roc_column | heatmap_column) +
  plot_layout(widths = c(1.35, 1))

## Combine panels a and b 
fig5ab <- wrap_plots(
  patchwork::free(fig5a, side = "t"), plot_spacer(), fig5b,
  nrow = 1,
  widths = c(1.1, 0.05, 1)
)

## Save the combined a/b figure 
figure_dir <- "figure/PART4"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(figure_dir, "fig5ab_association_prediction.pdf"),
  plot = fig5ab,
  width = 14, height = 6.8, units = "in"
)


# Fig 5c: single mediation ---------
rm(list = ls())
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"),
                                silent = TRUE)))
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# Configuration 
fig5c_cfg <- utils::modifyList(list(
  root = NULL,
  input_dir = "result/PART4/DIABLO/14.3_OWO_feature_mediation_mixed_cohorts/GDM_HDP",
  score_input_dir = "result/PART4/DIABLO/14.2_score_module_mediation_mixed_cohorts/GDM_HDP",
  figure_dir = "figure/PART4/DIABLO",
  outcomes = c("GDM", "HDP"), period = "T1", q_cutoff = 0.10,
  name_inside_threshold_percent = 15,
  name_size = 3.5, number_size = 3.5, label_offset = 0.65,
  equal_facet_width = FALSE,
  y_upper = 50, panel_tag = "c",
  colors = c(Metabolites = "#7EABCA", "16S_g" = "#FFB2B9")
), getOption("fig5c.diablo.single", list()))


# Data preparation helpers 
f5c_find_root <- function(root = NULL) {
  if (!is.null(root)) return(normalizePath(root, winslash = "/", mustWork = TRUE))
  starts <- getwd()
  for (fr in sys.frames()) if (!is.null(fr$ofile)) starts <- c(starts, dirname(fr$ofile))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) starts <- c(starts, dirname(sub("^--file=", "", arg[1L])))
  for (s in starts) {
    p <- normalizePath(s, winslash = "/", mustWork = FALSE)
    for (i in seq_len(10L)) {
      if (file.exists(file.path(p, "clean_code/14_pair_mediation_helpers.R"))) return(p)
      up <- dirname(p)
      if (identical(up, p)) break
      p <- up
    }
  }
  stop("Project root not found; set root in options(fig5c.diablo.single = list(root = ...)).")
}

f5c_check_fdr <- function(dat, families, level) {
  # Audit the stored q values using every planned test. Do not pool score and
  # feature tests, and do not replace the saved q values after filtering.
  groups <- if (level == "score") "outcome" else c("outcome", "mediator_type")
  keys <- unique(dat[groups])
  for (i in seq_len(nrow(keys))) {
    keep <- dat$outcome == keys$outcome[i]
    fkeep <- families$outcome == keys$outcome[i]
    if (level == "feature") {
      keep <- keep & dat$mediator_type == keys$mediator_type[i]
      fkeep <- fkeep & families$mediator_type == keys$mediator_type[i] & families$effect == "ACME"
    } else {
      # 14.2 labels the single-score indirect-effect family as IE; ACME is the
      # same mediate d.avg estimate and shares this q value.
      fkeep <- fkeep & families$analysis_type == "single_score" & families$effect == "IE"
    }
    z <- dat[keep, , drop = FALSE]
    fam <- families[fkeep, , drop = FALSE]
    if (nrow(fam) != 1L || fam$method != "BH" || fam$planned != nrow(z) || fam$fdr_n != nrow(z))
      stop("Unexpected ", level, " FDR family for ", keys$outcome[i])
    good <- z$status == "ok" & is.finite(z$ACME_p) & z$ACME_p >= 0 & z$ACME_p <= 1
    if (any(good)) {
      expected <- stats::p.adjust(z$ACME_p[good], "BH", n = fam$fdr_n)
      if (any(!is.finite(z$ACME_q[good])) || max(abs(expected - z$ACME_q[good])) > 1e-10)
        stop("Saved ", level, " ACME q values failed the complete-family audit.")
    }
  }
  invisible(TRUE)
}

f5c_prepare_data <- function(cfg) {
  fp <- file.path(cfg$root, cfg$input_dir)
  sp <- file.path(cfg$root, cfg$score_input_dir)
  files <- c(feature_results = file.path(fp, "all_results.csv"),
    stable_features = file.path(fp, "stable_features.csv"),
    feature_families = file.path(fp, "FDR_families.csv"),
    score_results = file.path(sp, "all_results.csv"),
    score_families = file.path(sp, "FDR_families.csv"))
  if (any(!file.exists(files))) stop("Missing inputs: ", paste(files[!file.exists(files)], collapse = ", "))
  features <- read_csv(files["feature_results"], show_col_types = FALSE) |>
    filter(outcome %in% cfg$outcomes, period == cfg$period,
      analysis_type == "OWO_feature_single", direction == "OWO_to_feature",
      exposure_id == "OWO", mediator_type %in% c("Metabolites", "Bacteria"))
  stable <- read_csv(files["stable_features"], show_col_types = FALSE) |>
    filter(stable_feature, module_status == "stable_crossomic_module")
  feature_families <- read_csv(files["feature_families"], show_col_types = FALSE)
  scores <- read_csv(files["score_results"], show_col_types = FALSE) |>
    filter(outcome %in% cfg$outcomes, period == cfg$period,
      analysis_type == "single_score",
      direction %in% c("OWO_to_M_score", "OWO_to_B_score")) |>
    mutate(mediator_type = if_else(direction == "OWO_to_M_score", "Metabolites", "Bacteria"),
      mediator_id = paste0(module_id, if_else(mediator_type == "Metabolites", "_M_score", "_B_score")),
      mediator_name = if_else(mediator_type == "Metabolites", "Metabolite score", "Bacterial score"))
  score_families <- read_csv(files["score_families"], show_col_types = FALSE)
  if (!nrow(features) || !nrow(scores)) stop("Both single-feature and single-score results are required.")
  if (anyDuplicated(features[c("outcome", "mediator_type", "mediator_id")]) ||
      anyDuplicated(scores[c("outcome", "direction")])) stop("Duplicate mediation models.")
  key <- function(y, b, id) paste(y, b, id, sep = "|")
  if (any(!key(features$outcome, features$mediator_type, features$mediator_id) %in%
          key(stable$outcome, stable$block, stable$feature_id)))
    stop("A single feature is outside the saved stable 14.1 modules.")
  if (any(scores$module_status != "stable_crossomic_module") ||
      any(!paste(scores$outcome, scores$module_id) %in% paste(stable$outcome, stable$module_id)))
    stop("Score and feature results do not refer to the same stable 14.1 modules.")
  f5c_check_fdr(features, feature_families, "feature")
  f5c_check_fdr(scores, score_families, "score")
  features <- features |>
    mutate(source_level = "Feature", is_score = FALSE, module_id = module_ids,
      FDR_scope = paste(outcome, mediator_type, "single_features", sep = "|"))
  scores <- scores |>
    mutate(source_level = "Score", is_score = TRUE,
      FDR_scope = paste(outcome, "single_scores_B_and_M", sep = "|"))
  all <- bind_rows(features, scores) |>
    mutate(selected_for_plot = status == "ok" & is.finite(PM) & PM > 0 &
             is.finite(ACME_q) & ACME_q < cfg$q_cutoff,
           selected_for_plot = coalesce(selected_for_plot, FALSE))
  z <- all |>
    filter(selected_for_plot) |>
    mutate(outcome = factor(outcome, levels = cfg$outcomes),
      mediator_type = factor(mediator_type, levels = c("Metabolites", "Bacteria")),
      omics = factor(if_else(mediator_type == "Metabolites", "Metabolites", "16S_g"),
                     levels = c("Metabolites", "16S_g")),
      PM_percent = 100 * PM,
      feature_label = trimws(gsub("[[:space:]\u00a0]+", " ", mediator_name)),
      feature_label = case_when(
        mediator_id == "MEDP1913" ~ "NAEPL",
        mediator_id == "MEDN0066" ~ "PAGln",
        TRUE ~ feature_label
      )) |>
    # Fixed score first, followed by descending individual-feature PM.
    arrange(outcome, mediator_type, desc(is_score), desc(PM_percent), feature_label, mediator_id) |>
    group_by(outcome) |>
    mutate(plot_position = row_number()) |>
    ungroup() |>
    mutate(PM_label = sprintf("%.2f", PM_percent),
      # Long names stay outside; short-bar values must remain above zero.
      name_inside = PM_percent > cfg$name_inside_threshold_percent &
        nchar(feature_label) <= PM_percent,
      name_y = if_else(name_inside, PM_percent - cfg$label_offset,
        pmax(PM_percent + cfg$label_offset, 4 + cfg$label_offset)),
      name_hjust = if_else(name_inside, 1, 0),
      name_color = if_else(name_inside, "white", "black"),
      name_fontface = if_else(!is_score & mediator_type == "Bacteria", "italic", "plain"),
      PM_y = if_else(name_inside, PM_percent + cfg$label_offset,
        pmax(PM_percent - cfg$label_offset, 4)),
      PM_hjust = if_else(name_inside, 0, 1),
      PM_color = if_else(name_inside | PM_y > PM_percent, "black", "white"))
  if (!nrow(z)) stop("No positive PM with ACME q < ", cfg$q_cutoff)
  if (any(z$PM_percent >= cfg$y_upper)) stop("A plotted PM reaches the y upper limit; increase y_upper.")
  if (anyNA(z$feature_label) || any(!nzchar(z$feature_label))) stop("Missing labels.")
  list(plot_data = z, all_results = all, feature_FDR_families = feature_families,
       score_FDR_families = score_families, input_files = files)
}


# Prepare the plotting data 
fig5c_cfg$root <- f5c_find_root(fig5c_cfg$root)
if (fig5c_cfg$q_cutoff <= 0 || fig5c_cfg$q_cutoff >= 1) stop("Invalid q cutoff.")
fig5c_input <- f5c_prepare_data(fig5c_cfg)
fig5c_plot_data <- fig5c_input$plot_data
fig5c_counts <- fig5c_plot_data |>
  count(outcome, mediator_type, source_level, name = "n_bars")


# Plot
fig5c <- ggplot(
  fig5c_plot_data,
  aes(x = plot_position, y = PM_percent, fill = omics)
) +
  geom_col(width = 0.7, color = NA) +
  geom_text(
    aes(y = name_y, label = feature_label, hjust = name_hjust,
        colour = name_color, fontface = name_fontface),
    angle = 90, vjust = 0.5, size = fig5c_cfg$name_size,
    show.legend = FALSE
  ) +
  geom_text(
    aes(y = PM_y, label = PM_label, hjust = PM_hjust, colour = PM_color),
    angle = 90, vjust = 0.5, size = fig5c_cfg$number_size,
    show.legend = FALSE
  ) +
  facet_grid(
    cols = vars(outcome), scales = "free_x",
    space = if (isTRUE(fig5c_cfg$equal_facet_width)) "fixed" else "free",
    drop = FALSE
  ) +
  scale_fill_manual(values = fig5c_cfg$colors, drop = FALSE) +
  scale_colour_identity() +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 0.40)) +
  scale_y_continuous(
    breaks = seq(0, fig5c_cfg$y_upper, 5), expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, fig5c_cfg$y_upper), clip = "off") +
  labs(x = NULL, y = "Proportion mediated (%)", fill = NULL, tag = fig5c_cfg$panel_tag) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank(),
    axis.title.y = element_text(size = 11, color = "black", margin = margin(r = 7)),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.4),
    axis.ticks.y = element_line(linewidth = 0.35),
    strip.background = element_rect(fill = "#faecb1"),
    strip.text = element_text(size = 11.5, color = "black"),
    panel.spacing.x = grid::unit(4, "mm"),
    legend.position = "none",
    plot.tag.position = c(0.005, 0.995),
    plot.tag = element_text(size = 18)
  ) +
  guides(fill = "none")


# Export 
figure_dir <- file.path(fig5c_cfg$root, fig5c_cfg$figure_dir)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
figure_stem <- file.path(figure_dir, "fig5c_DIABLO_single_mediation")

ggsave(
  paste0(figure_stem, ".pdf"), plot = fig5c,
  width = 7, height = 4.5, units = "in",
  device = grDevices::cairo_pdf, bg = "transparent"
)

# Fig 5d: pair mediation ---------
rm(list = ls())
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"),
                                silent = TRUE)))
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(circlize)
})


# Configuration 
fig5d_cfg <- utils::modifyList(list(
  root = NULL,
  input_dir = "result/PART4/DIABLO/14.3_bidir_feature_mediation_mixed_cohorts/GDM_HDP",
  figure_dir = "figure/PART4/DIABLO",
  outcome = "GDM", period = "T1", q_cutoff = 0.10,
  colors = c(Bact = "#FFB2B9", Met = "#7EABCA"),
  start_degree = 270, big_gap = c(forw = 45, back = 75), small_gap = 5,
  arrow_length = 0.09, label_size = 0.88, label_radius = 1.8,
  label_max_length = 28
), getOption("fig5d.diablo.pair", list()))


# Data preparation 
f5d_find_root <- function(root = NULL) {
  if (!is.null(root)) return(normalizePath(root, winslash = "/", mustWork = TRUE))
  starts <- getwd()
  for (fr in sys.frames()) if (!is.null(fr$ofile)) starts <- c(starts, dirname(fr$ofile))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) starts <- c(starts, dirname(sub("^--file=", "", arg[1L])))
  for (s in starts) {
    p <- normalizePath(s, winslash = "/", mustWork = FALSE)
    for (i in seq_len(10L)) {
      if (file.exists(file.path(p, fig5d_cfg$input_dir, "all_results.csv"))) return(p)
      up <- dirname(p)
      if (identical(up, p)) break
      p <- up
    }
  }
  stop("Project root not found; set root in options(fig5d.diablo.pair = list(root = ...)).")
}

f5d_check_fdr <- function(dat, families, outcome) {
  for (direction in c("B_to_M", "M_to_B")) {
    z <- dat |>
      filter(.data$outcome == .env$outcome, .data$direction == .env$direction)
    fam <- families |>
      filter(
        .data$outcome == .env$outcome, .data$direction == .env$direction,
        .data$effect == "ACME"
      )
    if (nrow(fam) != 1L || fam$method != "BH" || fam$planned != nrow(z) || fam$fdr_n != nrow(z)) {
      stop("Unexpected ACME FDR family for ", outcome, " / ", direction, ".")
    }
    good <- z$status == "ok" & is.finite(z$ACME_p) & z$ACME_p >= 0 & z$ACME_p <= 1
    expected <- p.adjust(z$ACME_p[good], method = "BH", n = fam$fdr_n)
    if (any(!is.finite(z$ACME_q[good])) || max(abs(expected - z$ACME_q[good])) > 1e-10) {
      stop("Saved ACME q values failed the complete-family audit for ", direction, ".")
    }
  }
  invisible(TRUE)
}

f5d_prepare_data <- function(cfg) {
  input_files <- c(
    results = file.path(cfg$root, cfg$input_dir, "all_results.csv"),
    families = file.path(cfg$root, cfg$input_dir, "FDR_families.csv")
  )
  if (any(!file.exists(input_files))) {
    stop("Missing inputs: ", paste(input_files[!file.exists(input_files)], collapse = ", "))
  }
  all_results <- read_csv(input_files["results"], show_col_types = FALSE) |>
    filter(
      .data$outcome == cfg$outcome, .data$period == cfg$period,
      .data$analysis_type == "feature_bidirectional",
      .data$direction %in% c("B_to_M", "M_to_B")
    )
  families <- read_csv(input_files["families"], show_col_types = FALSE)
  if (nrow(all_results) == 0L || anyDuplicated(all_results[c("pair_id", "direction")])) {
    stop("Missing or duplicated GDM bidirectional mediation models.")
  }
  f5d_check_fdr(all_results, families, cfg$outcome)

  all_results <- all_results |>
    mutate(
      selected_for_plot = status == "ok" & is.finite(PM) & PM > 0 &
        is.finite(ACME_q) & ACME_q < cfg$q_cutoff,
      selected_for_plot = coalesce(selected_for_plot, FALSE)
    )
  selected <- all_results |>
    filter(selected_for_plot) |>
    mutate(
      type = "Bact",
      name = str_squish(str_replace_all(bacteria_name, fixed(intToUtf8(160)), " ")),
      Compounds = str_squish(str_replace_all(metabolite_name, fixed(intToUtf8(160)), " ")),
      Compounds = case_when(
        metabolite_id == "MEDP1913" ~ "NAEPL",
        metabolite_id == "MEDN1162" ~ "p-tolyl sulfate",
        TRUE ~ Compounds
      ),
      prop.mediated = PM,
      pv_adj_FDR = ACME_q,
      direction_label = recode(
        direction,
        B_to_M = "Bacteria -> metabolite -> GDM",
        M_to_B = "Metabolite -> bacteria -> GDM"
      )
    ) |>
    arrange(direction, type, name, Compounds)
  if (!all(c("B_to_M", "M_to_B") %in% selected$direction)) {
    stop("At least one selected result is required in each mediation direction.")
  }
  list(selected = selected, all_results = all_results, families = families,
       input_files = input_files)
}

Prepare_Micros <- function(input, dir) {
  data <- input |> arrange(type, name, Compounds)
  value <- "prop.mediated"
  if (dir == "forw") {
    name_x <- "name"
    name_y <- "Compounds"
    sector_name <- c(data[[name_x]], data[[name_y]])
    sec_source <- c(data$type, rep("Met", length(data[[name_y]])))
  } else if (dir == "back") {
    name_x <- "Compounds"
    name_y <- "name"
    sector_name <- c(data[[name_x]], data[[name_y]])
    sec_source <- c(rep("Met", length(data[[name_x]])), data$type)
  } else {
    stop("dir must be 'forw' or 'back'.")
  }

  dat <- data |>
    mutate(trans = ifelse(pv_adj_FDR < fig5d_cfg$q_cutoff, 0.1, 0.8)) |>
    select(from = all_of(name_x), to = all_of(name_y), value = all_of(value), trans)
  dat$from <- as.character(dat$from)
  dat$to <- as.character(dat$to)

  gaps <- unique(data.frame(name = sector_name, source = sec_source)) |>
    mutate(color = case_match(
      source,
      "Bact" ~ unname(fig5d_cfg$colors["Bact"]),
      "Met" ~ unname(fig5d_cfg$colors["Met"])
    ))
  gaps$gap <- 1
  source_end <- c(which(gaps$source[-nrow(gaps)] != gaps$source[-1]), nrow(gaps))
  gaps$gap[source_end] <- 10
  list(dat = dat, gaps = gaps)
}

split_text <- function(text, max_length = fig5d_cfg$label_max_length) {
  if (nchar(text) <= max_length) return(text)
  parts <- strsplit(text, "")[[1]]
  lines <- character()
  current_line <- ""
  for (part in parts) {
    if (nchar(paste0(current_line, part)) + 1 > max_length) {
      lines <- c(lines, current_line)
      current_line <- part
    } else {
      current_line <- paste0(current_line, part)
    }
  }
  paste(c(lines, current_line), collapse = "\n")
}

Draw_chrod <- function(dat, gaps, dir, panel_title, start_deg = 270,
                       big_gap = 15, small_gap = 2) {
  gaps <- unique(gaps)
  grid.col <- setNames(gaps$color, gaps$name)
  if (dir == "forw") {
    arrow_color <- dat |> left_join(gaps, by = c("from" = "name")) |> pull(color)
  } else if (dir == "back") {
    arrow_color <- dat |> left_join(gaps, by = c("from" = "name")) |> pull(color)
  } else {
    stop("dir must be 'forw' or 'back'.")
  }

  circos.clear()
  circos.par(
    start.degree = start_deg,
    canvas.xlim = c(-1.55, 1.55), canvas.ylim = c(-1.55, 1.55),
    points.overflow.warning = FALSE
  )
  chordDiagram(
    dat[c("from", "to", "value")], big.gap = big_gap, small.gap = small_gap,
    grid.col = grid.col, col = arrow_color, transparency = dat$trans,
    directional = 1, direction.type = c("arrows", "diffHeight"),
    diffHeight = -0.02, annotationTrack = "grid",
    annotationTrackHeight = c(0.08, 0.1),
    link.arr.type = "big.arrow", link.arr.length = fig5d_cfg$arrow_length,
    link.largest.ontop = TRUE,
    preAllocateTracks = list(track.height = 0.1, track.margin = c(0.01, 0))
  )
  circos.trackPlotRegion(
    track.index = 2, bg.border = NA,
    panel.fun = function(x, y) {
      sector.index <- sapply(CELL_META$sector.index, split_text)
      circos.text(
        CELL_META$xcenter, fig5d_cfg$label_radius, sector.index,
        facing = "bending.outside", niceFacing = TRUE,
        cex = fig5d_cfg$label_size, col = "black", adj = c(0.5, 0.5)
      )
    }
  )
  title(main = panel_title, cex.main = 1.05, line = 0.3, col.main = "black")
  invisible(NULL)
}


# Prepare 
fig5d_cfg$root <- f5d_find_root(fig5d_cfg$root)
if (fig5d_cfg$q_cutoff <= 0 || fig5d_cfg$q_cutoff >= 1) stop("Invalid q cutoff.")
fig5d_input <- f5d_prepare_data(fig5d_cfg)
fig5d_selected <- fig5d_input$selected

forw_GDM <- fig5d_selected |> filter(direction == "B_to_M")
back_GDM <- fig5d_selected |> filter(direction == "M_to_B")
forw_DATA <- Prepare_Micros(forw_GDM, "forw")
back_DATA <- Prepare_Micros(back_GDM, "back")
fig5d_counts <- fig5d_selected |> count(direction, direction_label, name = "n_links")


# Plot and export 
figure_dir <- file.path(fig5d_cfg$root, fig5d_cfg$figure_dir)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

forw_stem <- file.path(figure_dir, "fig5d_DIABLO_pair_mediation_B_to_M")
grDevices::cairo_pdf(paste0(forw_stem, ".pdf"), width = 6.2, height = 6.2, bg = "transparent")
par(mar = c(0.2, 0.2, 2.2, 0.2))
Draw_chrod(forw_DATA$dat, forw_DATA$gaps, "forw", "Bacteria -> metabolite -> GDM",
           fig5d_cfg$start_degree, fig5d_cfg$big_gap["forw"], fig5d_cfg$small_gap)
dev.off()

back_stem <- file.path(figure_dir, "fig5d_DIABLO_pair_mediation_M_to_B")
grDevices::cairo_pdf(paste0(back_stem, ".pdf"), width = 6.2, height = 6.2, bg = "transparent")
par(mar = c(0.2, 0.2, 2.2, 0.2))
Draw_chrod(back_DATA$dat, back_DATA$gaps, "back", "Metabolite -> bacteria -> GDM",
           fig5d_cfg$start_degree, fig5d_cfg$big_gap["back"], fig5d_cfg$small_gap)
dev.off()

combined_stem <- file.path(figure_dir, "fig5d_DIABLO_pair_mediation")
grDevices::cairo_pdf(paste0(combined_stem, ".pdf"), width = 12, height = 6.2,
                     bg = "transparent", onefile = TRUE)
par(mfrow = c(1, 2), mar = c(0.2, 0.2, 2.2, 0.2), oma = c(0, 0, 1.8, 0))
Draw_chrod(forw_DATA$dat, forw_DATA$gaps, "forw", "Bacteria -> metabolite -> GDM",
           fig5d_cfg$start_degree, fig5d_cfg$big_gap["forw"], fig5d_cfg$small_gap)
Draw_chrod(back_DATA$dat, back_DATA$gaps, "back", "Metabolite -> bacteria -> GDM",
           fig5d_cfg$start_degree, fig5d_cfg$big_gap["back"], fig5d_cfg$small_gap)
mtext("d", side = 3, outer = TRUE, adj = 0.01, line = 0.1, cex = 1.5)
dev.off()
circos.clear()

