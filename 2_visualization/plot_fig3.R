####################################################
# Visualization and statistical testing for Figure 3
#####################################################


# Figure 3a: cluster-ordered circular heatmap of OWO-associated metabolites---------------
rm(list = ls())

library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(patchwork)
library(grid)

load("result/PART2/1_diff_metabolites.rdata")
load("result/PART2/4.1_OPLS_DA_longitudinal/4.1_OPLS_DA_longitudinal_results.RData")
load("data/01_thsbc_cleaned_data/metabolites_1292_991.RData")
load("data/project_col_set.Rdata")

cluster_assignment_file <- file.path(
  "result/PART2/6_metabolites_trend",
  "6.2 metabolite_cluster_assignments.csv"
)
cluster_assignments <- read_csv(cluster_assignment_file, show_col_types = FALSE) %>%
  select(feature, cluster_number, cluster) %>%
  distinct(feature, .keep_all = TRUE)

figure_dir <- "figure/PART2"
result_dir <- "result/PART2"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

period_order <- c("Pooled", "T3", "T2", "T1")
period_levels_bar <- c("Pooled", "T3", "T2", "T1")
association_exposure <- "OWO vs normal weight"
fdr_cutoff <- 0.05

## 1. Extract OWO-versus-normal association results 
owo_results <- res_pbmi_categorical %>%
  filter(
    exposure == association_exposure, period %in% period_order, status == "ok"
  ) %>%
  select(
    feature, period, beta, p_value, p_fdr,
    sugg_cmpd_name, Compounds, class_shrt
  ) %>%
  mutate(period = factor(period, levels = period_order))

# Include every metabolite significant in at least one of T1, T2, T3 or Pooled.
significant_metabolites <- owo_results %>%
  filter(!is.na(p_fdr), p_fdr < fdr_cutoff) %>%
  distinct(feature) %>%
  pull(feature)

metabolite_info <- info_metab %>%
  filter(Index %in% significant_metabolites) %>%
  transmute(
    feature = Index,
    # Use curated compound names for display labels
    metabolite_name = sugg_cmpd_name,
    class_shrt
  ) %>%
  distinct(feature, .keep_all = TRUE)

if (nrow(metabolite_info) != length(significant_metabolites)) {
  stop("Not all significant metabolites have a unique entry in info_metab.")
}
if (any(is.na(metabolite_info$metabolite_name) |
        str_trim(metabolite_info$metabolite_name) == "")) {
  stop(
    "sugg_cmpd_name is missing for: ",
    paste(
      metabolite_info$feature[
        is.na(metabolite_info$metabolite_name) |
          str_trim(metabolite_info$metabolite_name) == ""
      ],
      collapse = ", "
    )
  )
}

# Use the supplied class-colour dictionary. 
metabolite_info <- metabolite_info %>%
  mutate(
    metab_class = if_else(
      !is.na(class_shrt) & class_shrt %in% names(metab_class_col),
      class_shrt, "Others"
    )
  )

association_wide <- owo_results %>%
  filter(feature %in% significant_metabolites) %>%
  select(feature, period, beta, p_fdr) %>%
  pivot_wider(
    names_from = period,
    values_from = c(beta, p_fdr),
    names_glue = "{.value}_{period}"
  ) %>%
  left_join(metabolite_info, by = "feature")

expected_columns <- c(
  paste0("beta_", period_order), paste0("p_fdr_", period_order)
)
if (!all(expected_columns %in% names(association_wide))) {
  stop("One or more T1/T2/T3/Pooled association columns are absent.")
}

## 2. Join median VIP from the longitudinal OPLS-DA analysis 
vip_summary <- all_results$owo$results$primary_available_case$vip_summary %>%
  select(feature = metabolite, median_vip)

plot_data <- association_wide %>%
  left_join(vip_summary, by = "feature")

if (anyNA(plot_data$median_vip)) {
  stop(
    "Median VIP is unavailable for: ",
    paste(plot_data$feature[is.na(plot_data$median_vip)], collapse = ", ")
  )
}

## 3. Apply the coefficient-trajectory clusters and order by median VIP -
# Cluster membership is taken directly from 6.2_metabolites_coeffi_trend.R.
# Sectors are displayed as Cluster 1 -> Cluster 3; within each sector,
# metabolites are ordered by decreasing median VIP (name resolves exact ties).
plot_data <- plot_data %>%
  left_join(cluster_assignments, by = "feature")

missing_cluster <- plot_data %>%
  filter(is.na(cluster_number) | is.na(cluster)) %>%
  pull(feature)
if (length(missing_cluster) > 0L) {
  stop(
    "Coefficient-trajectory cluster is unavailable for: ",
    paste(missing_cluster, collapse = ", ")
  )
}

unexpected_cluster <- setdiff(unique(plot_data$cluster_number), 1:3)
if (length(unexpected_cluster) > 0L) {
  stop("Unexpected cluster number(s): ", paste(unexpected_cluster, collapse = ", "))
}

cluster_levels <- paste0("Cluster ", 1:3)
cluster_colours <- c(
  "Cluster 1" = "#D8745A",
  "Cluster 2" = "#FCC575",
  "Cluster 3" = "#A8C5A0"
)

plot_data <- plot_data %>%
  mutate(cluster = factor(cluster, levels = cluster_levels)) %>%
  arrange(cluster, desc(median_vip), metabolite_name) %>%
  mutate(order = row_number())

cluster_split <- factor(plot_data$cluster, levels = cluster_levels)

beta_matrix <- plot_data %>%
  select(all_of(paste0("beta_", period_order))) %>%
  as.matrix()
rownames(beta_matrix) <- plot_data$metabolite_name

fdr_matrix <- plot_data %>%
  select(all_of(paste0("p_fdr_", period_order))) %>%
  as.matrix()
rownames(fdr_matrix) <- plot_data$metabolite_name

finite_beta <- beta_matrix[is.finite(beta_matrix)]
beta_limit <- max(abs(finite_beta))
beta_limit
range(beta_matrix)
beta_col_fun <- colorRamp2(c(min(beta_matrix), 0, max(beta_matrix)), c("#3D4F86", "#F7F7F7", "#A51C30"))

class_colours <- metab_class_col
if (!"Others" %in% names(class_colours)) class_colours["Others"] <- "#D9D9D9"
active_classes <- names(class_colours)[names(class_colours) %in% unique(plot_data$metab_class)]

# Use one stack order for every period and direction. ggplot2 stacks factor
# levels in reverse order, so ascending class size puts the largest class at
# the bottom of every stacked bar. Alphabetical order resolves tied counts.
bar_stack_levels <- metabolite_info %>%
  count(metab_class, name = "n_metabolites") %>%
  arrange(n_metabolites, metab_class) %>%
  pull(metab_class)

## 4. Stacked counts of up- and downregulated significant metabolites 
bar_data <- owo_results %>%
  filter(
    feature %in% significant_metabolites,
    !is.na(p_fdr), p_fdr < fdr_cutoff,
    !is.na(beta), beta != 0
  ) %>%
  left_join(
    metabolite_info %>% select(feature, metab_class),
    by = "feature"
  ) %>%
  mutate(
    direction = if_else(beta > 0, "Up", "Down"),
    period = factor(as.character(period), levels = period_levels_bar),
    metab_class = factor(metab_class, levels = bar_stack_levels)
  ) %>%
  count(direction, period, metab_class, name = "n") %>%
  complete(
    direction = c("Down", "Up"),
    period = factor(period_levels_bar, levels = period_levels_bar),
    metab_class = factor(bar_stack_levels, levels = bar_stack_levels),
    fill = list(n = 0L)
  )

make_count_plot <- function(direction_value, show_x_axis) {
  count_data <- bar_data %>%
    filter(direction == direction_value) %>%
    group_by(period) %>%
    mutate(total = sum(n)) %>%
    ungroup()
  ymax <- max(count_data$total, na.rm = TRUE)

  ggplot(count_data, aes(period, n, fill = metab_class)) +
    geom_col(width = 0.8) +
    geom_text(
      data = distinct(count_data, period, total),
      aes(x = period, y = total, label = total),
      inherit.aes = FALSE, vjust = -0.25, size = 2.2
    ) +
    annotate(
      "text", x = 2.5, y = ymax * 1.22,
      label = direction_value, size = 2.5
    ) +
    scale_fill_manual(values = class_colours, drop = FALSE) +
    scale_x_discrete(expand = c(0.08, 0.08))+
    scale_y_continuous(
      limits = c(0, max(1, ymax * 1.32)), expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 8) +
    theme(
      plot.background = element_blank(),  
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text = element_text(colour = "black", size = 6),
      axis.text.x = if (show_x_axis) {
        element_text(colour = "black", size = 6)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
      plot.margin = margin(1, 2, 1, 2)
    ) 
}

down_plot <- make_count_plot("Down", show_x_axis = FALSE)
up_plot <- make_count_plot("Up", show_x_axis = TRUE)
bar_plot <- down_plot / up_plot +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

## 5. Circular heatmap, VIP and metabolite-class tracks 
draw_fig3a <- function() {
  circos.clear()
  circos.par(
    track.height = 0.20,
    start.degree = 90,
    clock.wise = TRUE,
    gap.after = c(2, 2, 88),
    track.margin = c(0.004, 0.004),
    cell.padding = c(0, 0.7, 0, 0.7),
    points.overflow.warning = FALSE,
    canvas.xlim = c(-0.85, 0.85), 
    canvas.ylim = c(-0.85, 0.85)
  )

  circos.heatmap.initialize(
    beta_matrix,
    split = cluster_split,
    cluster = FALSE
  )

  # Heatmap is ordered from outside to inside as Pooled, T3, T2 and T1.
  circos.heatmap(
    beta_matrix,
    col = beta_col_fun,
    bg.border = "white",
    cell.border = "white",
    cluster = FALSE,
    rownames.side = "outside",
    rownames.cex = 0.55,
    rownames.font = 1,
    track.height = 0.205
  )

  # Mark model-specific significance using the model-specific FDR values.
  heatmap_track_index <- get.current.track.index()
  for (cluster_name in cluster_levels) {
    sector_rows <- which(as.character(cluster_split) == cluster_name)
    sector_fdr <- fdr_matrix[sector_rows, , drop = FALSE]
    for (column_index in seq_len(ncol(sector_fdr))) {
      significant_index <- which(
        !is.na(sector_fdr[, column_index]) &
          sector_fdr[, column_index] < fdr_cutoff
      )
      if (length(significant_index) > 0L) {
        circos.points(
          significant_index - 0.5,
          ncol(sector_fdr) - column_index + 0.5,
          sector.index = cluster_name,
          track.index = heatmap_track_index,
          pch = "*", col = "black", cex = 0.8
        )
      }
    }
  }

  # Metabolite class is positioned directly between the association heatmap
  # and the median-VIP point-line track.
  circos.heatmap(
    matrix(plot_data$metab_class, ncol = 1),
    col = class_colours,
    bg.border = "white",
    cell.border = "white",
    cluster = FALSE,
    track.height = 0.025
  )

  # Display median VIP values across repeated training splits
  vip_values <- plot_data$median_vip
  vip_upper <- max(vip_values, na.rm = TRUE) * 1.08
  circos.track(
    ylim = c(0, vip_upper),
    bg.border = "black",
    track.height = 0.1,
    panel.fun = function(x, y) {
      track_vip <- vip_values[CELL_META$subset]
      track_x <- seq_along(track_vip) - 0.5
      if (CELL_META$sector.index == cluster_levels[1]) {
        circos.yaxis(
          side = "left",
          at = unique(c(0, 1, round(max(vip_values), 1))),
          labels.cex = 0.4,
          labels.niceFacing = FALSE,
          tick.length = convert_x(0.7, "mm")
        )
      }
      circos.segments(
        CELL_META$xlim[1], 1, CELL_META$xlim[2], 1,
        col = "black", lty = 3, lwd = 0.8
      )
      circos.lines(
        track_x, track_vip, type = "h",
        col = "#B98FAA", lwd = 1.3
      )
      circos.points(
        track_x, track_vip,
        pch = 16, cex = 0.6, col = "#B98FAA"
      )
    }
  )

  # A coloured annotation track identifies the three temporal coefficient
  # clusters. Labels follow the curvature of their respective sectors.
  circos.track(
    ylim = c(0, 1),
    bg.col = unname(cluster_colours[cluster_levels]),
    bg.border = "white",
    track.height = 0.035,
    panel.fun = function(x, y) {
      circos.text(
        CELL_META$xcenter, 0.5, CELL_META$sector.index,
        facing = "bending.inside", niceFacing = TRUE,
        cex = 0.55, col = "white", font = 2
      )
    }
  )

  # Central annotations are placed in the open sector.
  text(
    0, 0.02,
    paste0(length(significant_metabolites), " serum\nmetabolites"),
    cex = 1.22, adj = 0.5
  )
  text(0, -0.1, expression("* " * italic(P)[FDR] < 0.05), cex = 0.92)

  circos.clear()

  # Insert the vertically stacked Down/Up bar plots into the upper-left gap.
  print(
    bar_plot,
    newpage = FALSE,
    vp = viewport(
      x = unit(0.25, "npc"), y = unit(0.70, "npc"),
      width = unit(0.21, "npc"), height = unit(0.38, "npc")
    )
  )
  grid.text(
    "Number of significant metabolites",
    x = unit(0.135, "npc"), y = unit(0.7, "npc"), rot = 90,
    gp = gpar(fontsize = 9)
  )
  grid.text(
    "a",
    x = unit(0.14, "npc"), y = unit(0.87, "npc"),
    just = c("right", "bottom"),
    gp = gpar(fontsize = 18)
  )
  grid.text(
    period_order,
    x = unit(rep(0.493, length(period_order)), "npc"),
    y = unit(c(0.817, 0.787, 0.757, 0.727), "npc"),
    just = "right", gp = gpar(fontsize = 6.5)
  )

  # Heatmap and metabolite-class legends are shown together at bottom-left.
  beta_legend <- Legend(
    title = "Beta", col_fun = beta_col_fun,
    direction = "vertical",
    at = c(-0.3, 0, 0.5),
    labels = sprintf("%.1f", c(-0.3, 0, 0.5)),
    legend_height = unit(16, "mm"),
    grid_width = unit(2.5, "mm"),
    title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6.2)
  )
  class_legend <- Legend(
    title = "Class", at = active_classes,
    legend_gp = gpar(fill = class_colours[active_classes], col = NA),
    ncol = 2,
    title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6.2),
    grid_width = unit(2.5, "mm"), grid_height = unit(2.5, "mm"),
    gap = unit(1.2, "mm")
  )
  packed_legend <- packLegend(
    class_legend, beta_legend,
    direction = "horizontal", gap = unit(3, "mm")
  )
  draw(
    packed_legend,
    x = unit(0.16, "npc"), y = unit(0.1, "npc"),
    just = c("center", "center")
  )
}

cairo_pdf(
  file.path(figure_dir, "fig3av2_metab_circle_heatmap.pdf"),
  width = 7.5, height = 7.5
)
draw_fig3a()
dev.off()


# Figure 3b: pathway enrichment analysis ---------------
rm(list = ls())

library(tidyverse)
library(readxl)
library(ggsankey)
library(patchwork)
library(cowplot)
library(cols4all)

# Input 

input_dir <- "result/PART2/metab_pathway/Download"
output_dir <- "figure/PART2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

enrichment_raw <- read.csv(
  file.path(input_dir, "enrichment_pathora.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(enrichment_raw)[1] <- "pathway"

# Only pathways with a nominal/raw P value below 0.05 are displayed.
pathway_enrichment <- enrichment_raw %>%
  transmute(
    pathway = pathway,
    total = as.numeric(Total),
    expected = as.numeric(Expected),
    hits = as.numeric(Hits),
    raw_p = as.numeric(Raw_p),
    p_log10 = -log10(raw_p),
    impact = as.numeric(Impact),
    hit_ratio = hits / total
  ) %>%
  filter(is.finite(raw_p), raw_p < 0.05, hits > 0) %>%
  arrange(raw_p, pathway)

if (nrow(pathway_enrichment) == 0L) {
  stop("No pathway has a raw P value below 0.05.")
}

sankey_links <- read_xlsx(file.path(input_dir, "sankeydt.xlsx")) %>%
  transmute(
    pathway = as.character(pathway),
    metab = as.character(metab) %>%
      str_replace_all(fixed("<U+00A0>"), " ") %>%
      str_replace_all(fixed("<U+03B1>"), "alpha"),
    count = as.numeric(count)
  ) %>%
  filter(pathway %in% pathway_enrichment$pathway, is.finite(count), count > 0)

missing_pathways <- setdiff(pathway_enrichment$pathway, sankey_links$pathway)
if (length(missing_pathways) > 0L) {
  stop(
    "Significant pathways missing from sankeydt.xlsx: ",
    paste(missing_pathways, collapse = ", ")
  )
}

# Sankey panel 
# Higher-frequency metabolites are placed toward the top. Significant pathways
# are ordered from the smallest raw P value at the top to the largest at bottom.
metabolite_order <- sankey_links %>%
  group_by(metab) %>%
  summarise(link_count = sum(count), .groups = "drop") %>%
  arrange(link_count, metab) %>%
  pull(metab)

pathway_order <- pathway_enrichment$pathway
sankey_long <- sankey_links %>%
  ggsankey::make_long(metab, pathway)

node_levels <- c(metabolite_order, rev(pathway_order))
sankey_long <- sankey_long %>%
  mutate(node = factor(node, levels = node_levels))

node_palette <- cols4all::c4a("rainbow_wh_rd", length(node_levels))
names(node_palette) <- node_levels

p_sankey <- ggplot(
  sankey_long,
  aes(
    x = x, next_x = next_x,
    node = node, next_node = next_node,
    fill = node, label = node
  )
) +
  geom_sankey(
    space = 0.1,
    flow.alpha = 0.55,
    flow.fill = "grey82",
    flow.color = "grey75",
    flow.linetype = 1,
    flow.size = 0.2,
    smooth = 8,
    width = 0.08
  ) +
  geom_sankey_text(
    space = 0.1,
    size = 3.45,
    color = "black",
    hjust = 1
  ) +
  scale_fill_manual(values = node_palette, drop = FALSE) +
  scale_x_discrete(expand = expansion(add = c(0.75, 0.12))) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 0, 27, 12)
  )

# Extract the exact vertical centres of the right-hand pathway nodes generated
# by geom_sankey. The bubble panel therefore aligns with the Sankey layout even
# if the number of significant pathways or hits changes.
sankey_build <- ggplot_build(p_sankey)
sankey_node_data <- sankey_build$data[[2]]

pathway_layout <- sankey_node_data %>%
  filter(x == max(x)) %>%
  transmute(
    pathway = as.character(node),
    pathway_y = (ymin + ymax) / 2
  ) %>%
  distinct(pathway, .keep_all = TRUE) %>%
  inner_join(pathway_enrichment, by = "pathway")

if (nrow(pathway_layout) != nrow(pathway_enrichment)) {
  stop("Failed to recover a Sankey node position for every significant pathway.")
}

sankey_y_limits <- range(
  c(sankey_node_data$ymin, sankey_node_data$ymax),
  finite = TRUE
)

p_sankey <- p_sankey +
  scale_y_continuous(limits = sankey_y_limits, expand = c(0, 0))

# Bubble panel 
x_lower <- max(0, floor(min(pathway_layout$p_log10) * 2) / 2 - 0.15)
x_upper <- ceiling(max(pathway_layout$p_log10) * 2) / 2 + 0.35
x_breaks <- pretty(c(x_lower, x_upper), n = 5)
x_breaks <- x_breaks[x_breaks >= x_lower & x_breaks <= x_upper]

p_bubble <- ggplot(
  pathway_layout,
  aes(x = p_log10, y = pathway_y, size = hits, color = hit_ratio)
) +
  geom_point(alpha = 0.95) +
  scale_x_continuous(
    limits = c(x_lower, x_upper),
    breaks = x_breaks,
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = sankey_y_limits, expand = c(0, 0)) +
  scale_size_continuous(
    name = "Count",
    breaks = sort(unique(pathway_layout$hits)),
    range = c(2.8, 9)
  ) +
  scale_color_distiller(
    name = "Hit ratio",
    palette = "Reds",
    direction = 1
  ) +
  labs(
    x = expression(-log[10](italic(P)~"value")),
    y = NULL
  ) +
  guides(
    color = guide_colorbar(
      order = 1,
      title.position = "top",
      barheight = unit(3.1, "cm"),
      barwidth = unit(0.35, "cm")
    ),
    size = guide_legend(order = 2, title.position = "top")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10.5, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(8, 8, 27, 0)
  )

# Combined figure 
function_enrichment_core <-
  p_sankey + p_bubble +
  plot_layout(nrow = 1, widths = c(2.55, 1.0))

function_enrichment <- cowplot::ggdraw(function_enrichment_core) +
  cowplot::draw_label(
    "b",
    x = 0.085, y = 0.985,
    hjust = 0, vjust = 1,
    fontface = "plain",
    size = 22
  )

# Output 
ggsave(
  file.path(output_dir, "fig3b_function_enrichment.pdf"),
  function_enrichment,
  width = 9, height = 5.1,
  bg = "white"
)

# Figure 3c: cross-cohort consistency ---------------
rm(list = ls())
library(tidyverse)

# Input 
load("result/PART2/4.1_OPLS_DA_longitudinal/4.1_OPLS_DA_vip_sets.RData")
load("result/PART2/1_diff_metabolites_sig_sets_thsbc.rdata")
load("result/PART2/1_diff_metabolites.rdata")
load("result/PART2/2_diff_metabolites_HBC.rdata")

owo_vip_metabolites <- vip_sets$owo_vip
sig_metabolites <- sig_metab_sets$sigid_metab_owo
expected_periods <- c("Pooled", "T1", "T2", "T3")

# Prepare harmonized THSBC and HBC association results 
thsbc_results <- res_pbmi_categorical %>%
  filter(
    exposure == "OWO vs normal weight",
    period %in% expected_periods,
    status == "ok"
  ) %>%
  transmute(
    Index_xh = feature,
    period,
    beta_thsbc = beta,
    p_thsbc = p_value,
    p_fdr_thsbc = p_fdr
  )

hbc_results <- res_hbc_metab_pbmi_categorical %>%
  filter(
    exposure == "OWO vs normal weight",
    period %in% expected_periods,
    status == "ok",
    !is.na(Index_xh)
  ) %>%
  transmute(
    feature_hbc = feature,
    Index_xh,
    Compounds,
    period,
    beta_hbc = beta,
    p_hbc = p_value,
    p_fdr_hbc = p_fdr
  )

duplicate_thsbc <- thsbc_results %>%
  count(Index_xh, period) %>%
  filter(n > 1L)
duplicate_hbc <- hbc_results %>%
  count(Index_xh, period) %>%
  filter(n > 1L)

if (nrow(duplicate_thsbc) > 0L) {
  stop("Duplicated THSBC feature-period keys were detected.")
}
if (nrow(duplicate_hbc) > 0L) {
  stop("Duplicated HBC Index_xh-period keys were detected.")
}

shared_metabolites <- intersect(
  owo_vip_metabolites,
  unique(hbc_results$Index_xh)
)

message(
  "OWO OPLS-DA metabolites shared by THSBC and HBC: ",
  length(shared_metabolites), "/", length(owo_vip_metabolites)
)

replication_data <- hbc_results %>%
  filter(Index_xh %in% shared_metabolites) %>%
  inner_join(
    thsbc_results %>% filter(Index_xh %in% shared_metabolites),
    by = c("Index_xh", "period")
  ) %>%
  mutate(
    period = factor(period, levels = expected_periods),
    direction = case_when(
      beta_thsbc * beta_hbc <= 0 ~ "Null",
      beta_hbc > 0 ~ "Positive",
      beta_hbc < 0 ~ "Negative",
      TRUE ~ "Null"
    ),
    color_group = case_when(
      direction == "Null" ~ "null",
      direction == "Positive" & period == "Pooled" ~ "pos_pooled",
      direction == "Positive" & period == "T1" ~ "pos_t1",
      direction == "Positive" & period == "T2" ~ "pos_t2",
      direction == "Positive" & period == "T3" ~ "pos_t3",
      direction == "Negative" & period == "Pooled" ~ "neg_pooled",
      direction == "Negative" & period == "T1" ~ "neg_t1",
      direction == "Negative" & period == "T2" ~ "neg_t2",
      direction == "Negative" & period == "T3" ~ "neg_t3"
    ),
    replicated_fdr =
      p_fdr_thsbc < 0.05 & p_fdr_hbc < 0.05 & direction != "Null",
    point_stroke = if_else(replicated_fdr, 1.2, 0.35)
  )

period_counts <- replication_data %>% count(period, name = "n_pairs")
print(period_counts)

if (n_distinct(replication_data$Index_xh) != length(shared_metabolites)) {
  stop("The joined data do not contain every shared metabolite.")
}
if (any(period_counts$n_pairs != length(shared_metabolites))) {
  stop("At least one period does not contain exactly one pair per shared metabolite.")
}

# Correlation across all shared metabolites and all four model periods 
correlation_test <- cor.test(
  replication_data$beta_thsbc,
  replication_data$beta_hbc,
  method = "spearman",
  exact = FALSE
)

correlation_r <- unname(correlation_test$estimate)
correlation_p <- correlation_test$p.value

# Period-specific correlations use the same shared-metabolite set and the same
# Spearman method as the pooled four-period analysis above.
period_correlation_tests <- set_names(
  map(expected_periods, function(period_name) {
    period_data <- replication_data %>%
      filter(period == period_name) %>%
      select(beta_thsbc, beta_hbc) %>%
      filter(if_all(everything(), is.finite))

    if (nrow(period_data) < 3L) {
      stop("Fewer than three complete effect-size pairs were available for ", period_name, ".")
    }

    cor.test(
      period_data$beta_thsbc,
      period_data$beta_hbc,
      method = "spearman",
      exact = FALSE
    )
  }),
  expected_periods
)

period_correlation_results <- imap_dfr(
  period_correlation_tests,
  function(test, period_name) {
    tibble(
      period = period_name,
      n_pairs = sum(
        replication_data$period == period_name &
          is.finite(replication_data$beta_thsbc) &
          is.finite(replication_data$beta_hbc)
      ),
      spearman_r = unname(test$estimate),
      p_value = test$p.value,
      p_value_text = if_else(
        test$p.value < 0.001,
        "<0.001",
        sprintf("%.3f", test$p.value)
      )
    )
  }
) %>%
  mutate(period = factor(period, levels = expected_periods)) %>%
  arrange(period) %>%
  mutate(period = as.character(period))

format_period_correlation_label <- function(results, heading) {
  period_lines <- results %>%
    mutate(
      p_display = if_else(
        p_value < 0.001,
        "P < 0.001",
        sprintf("P = %.3f", p_value)
      ),
      line = sprintf("%s: r = %.3f, %s;", period, spearman_r, p_display)
    ) %>%
    pull(line)

  paste(c(heading, period_lines), collapse = "\n")
}

correlation_label <- format_period_correlation_label(
  period_correlation_results,
  "All shared key metabolites"
)

# Plot 
point_colors <- c(
  "pos_pooled" = "#6D011F",
  "pos_t1" = "#B72230",
  "pos_t2" = "#DC6D57",
  "pos_t3" = "#F6B293",
  "neg_pooled" = "#104680",
  "neg_t1" = "#317CB7",
  "neg_t2" = "#6DADD1",
  "neg_t3" = "#B6D7E8",
  "null" = "#A9A9A9"
)

legend_breaks <- c(
  "neg_pooled", "pos_pooled",
  "neg_t1", "pos_t1",
  "neg_t2", "pos_t2",
  "neg_t3", "pos_t3",
  "null"
)
legend_labels <- c(
  "Pooled", "Pooled", "T1", "T1", "T2", "T2", "T3", "T3", "Null"
)

make_replication_plot <- function(plot_data, annotation_label, point_size = 3.5) {
  ggplot(
    plot_data,
    aes(x = beta_thsbc, y = beta_hbc, fill = color_group)
  ) +
    geom_smooth(
      method = "lm", formula = y ~ x,
      color = "black", fill = "grey80", alpha = 0.45,
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.5
    ) +
    geom_vline(
      xintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.5
    ) +
    geom_point(
      shape = 21, size = point_size, color = "black",
      stroke = plot_data$point_stroke
    ) +
    scale_fill_manual(
      values = point_colors,
      breaks = legend_breaks,
      labels = legend_labels,
      name = "Neg.     Pos.",
      drop = FALSE
    ) +
    annotate(
      "text", x = -Inf, y = Inf,
      hjust = -0.04, vjust = 1.1,
      label = annotation_label, size = 3.8, lineheight = 1.15
    ) +
    labs(
      x = "Effect size (THSBC discovery)",
      y = "Effect size (HBC validation)"
    ) +
    guides(
      fill = guide_legend(
        ncol = 2, byrow = TRUE,
        override.aes = list(shape = 21, size = 3.2, color = "black")
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.background = element_blank(),
      legend.position = c(0.82, 0.13),
      legend.title = element_text(size = 9, hjust = 0.5),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10, color = "black")
    )
}

fig3c_base <- make_replication_plot(replication_data, correlation_label)
fig3c <- cowplot::ggdraw(fig3c_base) +
  cowplot::draw_label(
    "c",
    x = 0.012, y = 0.992,
    hjust = 0, vjust = 1,
    size = 20
  )

# Output 
dir.create("figure/PART2", recursive = TRUE, showWarnings = FALSE)
save(
  replication_data,
  shared_metabolites,
  correlation_test,
  period_correlation_tests,
  period_correlation_results,
  fig3c,
  file = "result/PART2/fig3c_metabolites_replication.RData"
)
ggsave(
  "figure/PART2/fig3c_metabolites_replication.pdf",
  plot = fig3c, width = 4.8, height = 4.2
)


# Figure 3d: temporal clustering ---------------
# Run this plotting section in the same R session immediately after the
# temporal clustering analysis in 2_metabolomic_analysis.R.
# Required clustering objects are used directly from the analysis environment.

metric_plot_data <- cluster_validity_metrics %>%
  pivot_longer(
    cols = all_of(validity_metric_names),
    names_to = "metric",
    values_to = "value"
  )

validity_plot <- ggplot(
  metric_plot_data,
  aes(x = n_clusters, y = value)
) +
  geom_line(linewidth = 0.55, colour = "#365F73") +
  geom_point(aes(size = minimum_cluster_size), colour = "#365F73") +
  geom_vline(
    xintercept = selected_k,
    linetype = 2,
    linewidth = 0.5,
    colour = "#D6604D"
  ) +
  facet_wrap(~metric, scales = "free_y", nrow = 2) +
  scale_x_continuous(breaks = candidate_k) +
  scale_size_continuous(range = c(1.5, 3.5)) +
  labs(
    x = "Number of clusters",
    y = NULL,
    size = "Minimum\ncluster size",
    tag = 'a'
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    plot.background = element_blank(), 
    legend.background = element_blank(), 
    legend.position = c(0.9, 0.25),
    axis.text = element_text(size = 11, color = 'black'),
    axis.title = element_text(size = 12, color = 'black'),
    legend.text = element_text(size = 10, color = 'black'),
    legend.title = element_text(size = 11, color = 'black'),
    strip.text = element_text(size = 12, color = 'black'),
    plot.tag = element_text(size = 16, color = 'black'),
    plot.tag.position = c(0.01, 0.99)
  )

ggsave(
  file.path(figure_dir, "6.2 coefficient_cluster_validity_metrics.pdf"),
  validity_plot,
  width = 6,
  height = 4.6
)

cluster_labels <- cluster_summary %>%
  transmute(
    cluster,
    label = sprintf(
      "%s (n = %d)\n%s",
      cluster, n_metabolites, pattern
    )
  ) %>%
  deframe()

cluster_colours <- setNames(
  c("#D8745A", "#FCC575", "#A8C5A0")[seq_len(selected_k)],
  paste0("Cluster ", seq_len(selected_k))
)

# Summarise the representation of the 40 OWO VIP metabolites in each
# coefficient-trajectory cluster. The Fisher test compares the distribution
# of VIP versus non-VIP metabolites across all selected clusters.
vip_cluster_annotations <- cluster_assignments %>%
  mutate(is_owo_vip = feature %in% owo_vip_metabolites) %>%
  group_by(cluster_number, cluster) %>%
  summarise(
    n_cluster = n_distinct(feature),
    n_owo_vip = n_distinct(feature[is_owo_vip]),
    .groups = "drop"
  ) %>%
  mutate(
    vip_proportion = n_owo_vip / n_cluster,
    label = sprintf(
      "Key metabolites: %d/%d = %.1f%%",
      n_owo_vip, n_cluster, 100 * vip_proportion
    )
  )

vip_cluster_table <- cluster_assignments %>%
  distinct(feature, cluster_number) %>%
  mutate(is_owo_vip = feature %in% owo_vip_metabolites) %>%
  with(table(cluster_number, is_owo_vip))

vip_cluster_fisher <- fisher.test(vip_cluster_table)
vip_cluster_test_label <- sprintf(
  "Fisher's exact test (cluster x key-metabolite status): P = %.2e",
  vip_cluster_fisher$p.value
)

coefficient_trajectory_plot <- ggplot(
  trend_clustered,
  aes(x = time, y = beta, group = feature, colour = cluster)
) +
  geom_hline(
    yintercept = 0, linetype = 2, linewidth = 0.45, colour = "grey70"
  ) +
  geom_line(alpha = 0.16, linewidth = 0.4) +
  geom_point(alpha = 0.20, size = 0.7) +
  geom_line(
    data = cluster_centroids,
    aes(
      x = time, y = mean_beta, group = cluster, colour = cluster
    ), inherit.aes = FALSE, linewidth = 1.1
  ) +
  geom_point(
    data = cluster_centroids,
    aes(x = time, y = mean_beta, colour = cluster),
    inherit.aes = FALSE, size = 2
  ) +
  geom_text(
    data = vip_cluster_annotations,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.5,
    size = 4.3, colour = "black"
  ) +
  facet_wrap(
    ~cluster, nrow = 1,
    labeller = labeller(cluster = cluster_labels)
  ) +
  scale_colour_manual(values = cluster_colours) +
  scale_x_continuous(
    breaks = seq_along(expected_periods),
    labels = expected_periods,
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  labs(
    x = "Timepoint",
    y = "Beta (OWO vs normal weight)",
    title = vip_cluster_test_label,
    tag = 'd'
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    legend.position = "none",
    panel.background = element_blank(), 
    plot.background = element_blank(), 
    legend.background = element_blank(), 
    axis.text = element_text(size = 11, color = 'black'),
    axis.title = element_text(size = 12, color = 'black'),
    legend.text = element_text(size = 10, color = 'black'),
    legend.title = element_text(size = 11, color = 'black'),
    strip.text = element_text(size = 12, color = 'black'),
    plot.title.position = "plot",
    plot.title = element_text(size = 10.5, color = "black", hjust = 1),
    plot.tag = element_text(size = 22, color = 'black'),
    plot.tag.position = c(0.01, 0.98)
  )

ggsave(
  file.path(figure_dir, "6.2 owo_metabolite_coefficient_trajectories_v2.pdf"),
  coefficient_trajectory_plot,
  width = 8, height = 3.2
)


# Figure 3e: pBMI and gestational-week contributions to OWO-related metabolites ---------------

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(grid)
})

## 1. Input data 
load("result/PART2/3.3_variance_icc.rdata")
load("result/PART2/4.1_OPLS_DA_longitudinal/4.1_OPLS_DA_40metab_id_thsbc.rdata")
load("result/PART2/6_metabolites_trend/6.2_metabolites_coeffi_trend.RData")

figure_dir <- "figure/PART2"
result_dir <- "result/PART2"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

cluster_levels <- paste0("Cluster ", 1:3)
cluster_colours <- c(
  "Cluster 1" = "#D8745A",
  "Cluster 2" = "#FCC575",
  "Cluster 3" = "#A8C5A0"
)
component_colours <- c(
  "pBMI" = "#f2e1a2",
  "Gestational week" = "#a3c1e0"
)

metabolite_trajectory_clusters <- cluster_assignments %>%
  select(feature, cluster_number, cluster) %>%
  distinct(feature, .keep_all = TRUE)

## 2. Prepare the 40 prioritized OWO-associated metabolites
plot_data <- variance_results %>%
  filter(feature %in% unique(owo_metabolites), status == "ok") %>%
  inner_join(metabolite_trajectory_clusters, by = "feature") %>%
  mutate(
    p_interaction_fdr = p.adjust(p_interaction_lrt, method = "fdr"),
    interaction_sig = if_else(
      !is.na(p_interaction_fdr) & p_interaction_fdr < 0.05,
      "*", ""
    ),
    shapley_pbmi_pct = 100 * shapley_pbmi,
    shapley_wk_pct = 100 * shapley_wk,
    metabolite = sugg_cmpd_name,
    metabolite = metabolite %>%
      str_replace_all(fixed("<U+00A0>"), " ") %>%
      str_replace_all(fixed("<U+03B1>"), "alpha"),
    cluster = factor(cluster, levels = cluster_levels)
  ) %>%
  arrange(cluster, desc(shapley_pbmi), feature)

plot_data <- plot_data %>%
  mutate( order = row_number(), x= order)

shapley_long <- plot_data %>%
  select(
    feature, metabolite, cluster_number, cluster, x,
    shapley_pbmi_pct, shapley_wk_pct
  ) %>%
  pivot_longer(
    cols = c(shapley_pbmi_pct, shapley_wk_pct),
    names_to = "component", values_to = "shapley_pct"
  ) %>%
  mutate(
    component = factor(
      component,
      levels = c("shapley_pbmi_pct", "shapley_wk_pct"),
      labels = c("pBMI", "Gestational week")
    )
  )

cluster_ranges <- plot_data %>%
  group_by(cluster_number, cluster) %>%
  summarise(
    xmin = min(x) - 0.5,
    xmax = max(x) + 0.5,
    xmid = mean(c(xmin, xmax)),
    .groups = "drop"
  )

## 3. Upper lollipop plot 
y_max <- max(c(plot_data$shapley_pbmi_pct, plot_data$shapley_wk_pct), na.rm = TRUE)
x_limits <- c(min(plot_data$x) - 0.7, max(plot_data$x) + 0.7)
interaction_limit <- max(abs(plot_data$interaction_beta), na.rm = TRUE)
if (!is.finite(interaction_limit) || interaction_limit == 0) interaction_limit <- 1

fig3e_lollipop <- ggplot() +
  geom_segment(
    data = plot_data,
    aes( x = x, xend = x,  y = shapley_wk_pct, yend = shapley_pbmi_pct),
    colour = "grey70", linewidth = 0.55
  ) +
  geom_point(
    data = shapley_long,
    aes(x = x, y = shapley_pct, colour = component),
    size = 3
  ) +
  scale_colour_manual(values = component_colours, guide = "none") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.05, 0))
  ) +
  coord_cartesian(ylim = c(min(0, min(plot_data$shapley_wk_pct)), 1.15 * y_max)) +
  labs(
    x = NULL,
    y = expression("Shapley contribution (" * R^2 * ", %)"),
    tag = "e"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    plot.tag = element_text(size = 18),
    plot.tag.position = c(0.005, 0.995),
    plot.margin = margin(t = 5, r = 5, b = -2, l = 5)
  )

## 4. Lower interaction heatmap, cluster strip and metabolite names 

fig3e_annotation <- ggplot() +
  geom_rect(
    data = cluster_ranges,
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 0.5, fill = cluster),
    inherit.aes = FALSE, colour = NA
  ) +
  geom_text(
    data = cluster_ranges,
    aes(x = xmid, y = 0.25, label = cluster),
    inherit.aes = FALSE, size = 4, colour = "black"
  ) +
  scale_fill_manual(values = cluster_colours, guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = plot_data,
    aes(x = x, y = 0.75, fill = interaction_beta),
    width = 1, height = 0.5,
    colour = NA
  ) +
  geom_text(
    data = plot_data %>% filter(interaction_sig == "*"),
    aes(x = x, y = 0.75, label = interaction_sig),
    colour = "white", size = 6.5, vjust = 0.72
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-interaction_limit, interaction_limit), guide = "none"
  ) +
  annotate(
    "text", x = min(plot_data$x) - 0.70, y = 0.75,
    label = "Beta", angle = 90, size = 3.5
  ) +
  scale_x_continuous(
    limits = x_limits,
    breaks = plot_data$x,
    labels = plot_data$metabolite,
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5, colour = "black", size = 9
    ),
    axis.line.x = element_line(colour = "black", linewidth = 0.45),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.4),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(2.2, "pt"),
    legend.position = "none",
    plot.margin = margin(t = -2, r = 5, b = 4, l = 5)
  )

## 5. Extract and inset the two legends into the upper-right open region
legend_source <- ggplot() +
  geom_point(
    data = shapley_long,
    aes(x = 1, y = 1, colour = component), size = 3.5
  ) +
  geom_tile(
    data = plot_data,
    aes(x = 1, y = 1, fill = interaction_beta)
  ) +
  scale_colour_manual(
    name = NULL, values = component_colours,
    guide = guide_legend(
      order = 1, direction = "horizontal", nrow = 1,
      override.aes = list(size = 3)
    )
  ) +
  scale_fill_gradient2(
    name = "Interaction beta",
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-interaction_limit, interaction_limit),
    breaks = c(-interaction_limit, 0, interaction_limit),
    labels = scales::label_number(accuracy = 0.001),
    guide = guide_colourbar(
      order = 2, direction = "horizontal", title.position = "top",
      barwidth = unit(28, "mm"), barheight = unit(3, "mm")
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key.height = unit(3.5, "mm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

legend_grob <- cowplot::get_legend(legend_source)
legend_panel <- cowplot::ggdraw(legend_grob)

fig3e_lollipop_inset <- fig3e_lollipop +
  inset_element(
    legend_panel,
    left = 0.74, bottom = 0.6, right = 0.98, top = 0.94,
    align_to = "panel", clip = FALSE
  )

fig3e_metab_wk_r2_interaction <-
  fig3e_lollipop_inset / plot_spacer() / fig3e_annotation +
  plot_layout(heights = c(16, 0.15, 1.4))

## 6. Export 
ggsave(
  file.path(figure_dir, "fig3e_metab_wk_r2_interaction.pdf"),
  fig3e_metab_wk_r2_interaction,
  width = 10, height = 6
)


# Figure 3f: clinical relevance ---------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

# Input 
thsbc_env <- new.env(parent = emptyenv())
hbc_env <- new.env(parent = emptyenv())

load(
  "result/PART2/5.1_association_metab_clinical_traits_thsbc.RData",
  envir = thsbc_env
)
load(
  "result/PART2/5.2_association_metab_clinical_traits_hbc.RData",
  envir = hbc_env
)

if (!exists("score_association_results", envir = thsbc_env, inherits = FALSE)) {
  stop("THSBC score_association_results was not found in the 5.1 result file.")
}
if (!exists("score_association_results", envir = hbc_env, inherits = FALSE)) {
  stop("HBC score_association_results was not found in the 5.2 result file.")
}

# Harmonize clinical-trait categories 
clinical_group_levels <- c(
  "Blood pressure", "Blood glucose", "Liver", "Renal", "Thyroid"
)

clinical_group_colors <- c(
  "Blood pressure" = "#F8766D",
  "Blood glucose" = "#A3A500",
  "Liver" = "#00BF7D",
  "Renal" = "#00B0F6",
  "Thyroid" = "#BC80BD"
)

add_clinical_group <- function(data) {
  data %>%
    mutate(
      clinical_group = case_when(
        clinical_trait %in% c("SBP", "DBP") ~ "Blood pressure",
        clinical_trait %in% c("FPG", "OGTT_0h", "OGTT_1h", "OGTT_2h") ~
          "Blood glucose",
        clinical_trait %in% c(
          "ALB", "TP", "a_HBDH", "A/G Ratio", "CHE", "GGT", "ALT",
          "AST", "IBIL", "ALP", "GLB", "LAC", "LDH", "ADA", "DBIL",
          "TBIL"
        ) ~ "Liver",
        clinical_trait %in% c("Cr", "Ccr", "BUN", "UA", "CysC") ~ "Renal",
        clinical_trait %in% c("FT3", "FT4", "Tg", "TPOAb", "TSH") ~ "Thyroid",
        TRUE ~ NA_character_
      ),
      clinical_group = factor(clinical_group, levels = clinical_group_levels)
    ) %>%
    filter(!is.na(clinical_group))
}

prepare_pc1_results <- function(results, cohort) {
  required_columns <- c(
    "clinical_trait", "feature_type", "metabolite", "clinical_period",
    "metabolite_period", "n_obs", "beta", "se", "lower_ci", "upper_ci",
    "p_value", "status"
  )
  missing_columns <- setdiff(required_columns, names(results))
  if (length(missing_columns) > 0L) {
    stop(
      cohort, " result is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  results %>%
    filter(
      feature_type == "PC1 score",
      metabolite == "metabolite_score",
      metabolite_period == "T1",
      status == "ok"
    ) %>%
    distinct(clinical_trait, .keep_all = TRUE) %>%
    mutate(
      p_fdr_figure = p.adjust(p_value, method = "fdr"),
      log10_fdr = -log10(pmax(p_fdr_figure, .Machine$double.xmin)),
      cohort = cohort
    ) %>%
    add_clinical_group()
}

thsbc_pc1 <- prepare_pc1_results(
  get("score_association_results", envir = thsbc_env),
  cohort = "THSBC"
)
hbc_pc1 <- prepare_pc1_results(
  get("score_association_results", envir = hbc_env),
  cohort = "HBC"
)

if (anyDuplicated(thsbc_pc1$clinical_trait)) {
  stop("Duplicated THSBC clinical traits remain after data preparation.")
}
if (anyDuplicated(hbc_pc1$clinical_trait)) {
  stop("Duplicated HBC clinical traits remain after data preparation.")
}

# Left panel: THSBC PC1 associations
thsbc_plot_data <- thsbc_pc1 %>%
  mutate(
    label = if_else(p_fdr_figure < 0.05, clinical_trait, NA_character_),
    significant_fdr = p_fdr_figure < 0.05
  )

fig3f_left <- ggplot( thsbc_plot_data, aes(x = beta, y = log10_fdr, fill = clinical_group)) +
  geom_hline(
    yintercept = -log10(0.05),
    color = "darkgrey", linetype = "dashed", linewidth = 0.5
  ) +
  geom_point(
    aes(size = abs(beta)), shape = 21, color = "black", stroke = 0.3
  ) +
  geom_text_repel(
    data = thsbc_plot_data %>% filter(!is.na(label)),
    aes(label = label),
    size = 3.5, box.padding = 0.25, point.padding = 0.15, direction="both",
    min.segment.length = 0.5, max.overlaps = Inf, seed = 111,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = clinical_group_colors,
    breaks = clinical_group_levels,
    drop = FALSE
  ) +
  scale_size_continuous(range = c(2.2, 5.2), guide = "none") +
  labs(
    x = "Effect size (THSBC)",
    y = expression(-log[10](italic(P)[FDR])),
    fill = NULL, tag = 'f'
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, color = NA, size = 4)
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(14, "pt"),
    legend.position = c(0.2, 0.85),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    plot.tag.position = c(0.02, 0.97),
    plot.tag = element_text(size = 22, colour = "black")
  )

# Right panel: cross-cohort consistency 
consistency_data <- thsbc_pc1 %>%
  select(
    clinical_trait, clinical_group,
    beta_thsbc = beta, p_thsbc = p_value, p_fdr_thsbc = p_fdr_figure
  ) %>%
  inner_join(
    hbc_pc1 %>%
      select(
        clinical_trait,
        beta_hbc = beta, p_hbc = p_value, p_fdr_hbc = p_fdr_figure
      ),
    by = "clinical_trait"
  ) %>%
  mutate(
    direction_consistent = beta_thsbc * beta_hbc > 0,
    replicated_fdr =
      p_fdr_thsbc < 0.05 & p_fdr_hbc < 0.05 & direction_consistent,
    hbc_significant_fdr = p_fdr_hbc < 0.05,
    point_stroke = if_else(hbc_significant_fdr, 1.2, 0.35),
    label = if_else(direction_consistent, clinical_trait, NA_character_)
  )

if (nrow(consistency_data) < 3L) {
  stop("Fewer than three shared clinical traits are available for correlation.")
}

correlation_test <- cor.test(
  consistency_data$beta_thsbc,
  consistency_data$beta_hbc,
  method = "spearman"
)
correlation_r <- unname(correlation_test$estimate)
correlation_p <- correlation_test$p.value

correlation_label <- if (correlation_p < 0.001) {
  sprintf("italic(r) == %.3f*','~~italic(P) < 0.001", correlation_r)
} else {
  sprintf(
    "italic(r) == %.3f*','~~italic(P) == %.3f",
    correlation_r, correlation_p
  )
}

fig3f_right <- ggplot( consistency_data, aes(x = beta_thsbc, y = beta_hbc, fill = clinical_group)) +
  geom_hline(
    yintercept = 0, color = "darkgrey", linetype = "dashed", linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0, color = "darkgrey", linetype = "dashed", linewidth = 0.5
  ) +
  geom_point(
    shape = 21, size = 4, color = "black",
    stroke = consistency_data$point_stroke
  ) +
  geom_text_repel(
    data = consistency_data %>% filter(!is.na(label), clinical_trait != "FPG"),
    aes(label = label),
    size = 3.5, box.padding = 0.3, point.padding = 0.2,
    min.segment.length = 0.5, max.overlaps = Inf, seed = 111,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = consistency_data %>% filter(clinical_trait == "FPG"),
    aes(label = label),
    nudge_x = -0.02, nudge_y = 0.015,
    size = 3.5, box.padding = 0.3, point.padding = 0.2,
    min.segment.length = 0.5, max.overlaps = Inf, seed = 111,
    show.legend = FALSE
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    hjust = -0.1, vjust = 2,
    label = correlation_label, parse = TRUE, size = 4.5
  ) +
  scale_fill_manual(
    values = clinical_group_colors,
    breaks = clinical_group_levels,
    drop = FALSE
  ) +
  labs( x = "Effect size (THSBC)", y = "Effect size (HBC)" ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black")
  )

fig3f <- fig3f_left + fig3f_right +
  plot_layout(ncol = 2, widths = c(1, 1))

# Output
dir.create("figure/PART2", recursive = TRUE, showWarnings = FALSE)
ggsave(
  "figure/PART2/fig3f_metabolites_clinical.pdf",
  plot = fig3e, width = 8, height = 4
)


# Figure 3g: determinants of OWO-related serum metabolites ------------------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

# Input and output 
determinant_dir <- "result/PART2/7_determinant_metabolites"
figure_dir <- "figure/PART2"
result_dir <- "result/PART2"

individual_file <- file.path(
  determinant_dir,
  "7-11_individual_metabolite_marginal_r2.csv"
)
matrix_file <- file.path(
  determinant_dir,
  "7-06_matrix_dbrda_period_specific.csv"
)
metabolite_info_file <-
  "data/01_thsbc_cleaned_data/metabolites_1292_991.RData"
class_colour_file <- "data/project_col_set.Rdata"

required_files <- c(
  individual_file, matrix_file, metabolite_info_file, class_colour_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing input file(s): ", paste(missing_files, collapse = ", "))
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

individual_results <- readr::read_csv(individual_file, show_col_types = FALSE)
matrix_results <- readr::read_csv(matrix_file, show_col_types = FALSE)
load(metabolite_info_file)
load(class_colour_file)

if (!exists("info_metab")) {
  stop("Object 'info_metab' was not found in ", metabolite_info_file)
}
if (!exists("metab_class_col")) {
  stop("Object 'metab_class_col' was not found in ", class_colour_file)
}

# Configuration 
period_levels <- c("T1", "T2", "T3")
group_levels <- c(
  "Demographic factors",
  "Behavioral factors",
  "Dietary factors",
  "Socioeconomic factors",
  "Genetics",
  "Gut microbiota"
)

group_labels <- c(
  demographics = "Demographic factors",
  lifestyle = "Behavioral factors",
  diet = "Dietary factors",
  social = "Socioeconomic factors",
  genetics = "Genetics",
  microbiota = "Gut microbiota"
)

group_colors <- c(
  "Demographic factors" = "#b2e1d4",
  "Behavioral factors" = "#f6e2b3",
  "Dietary factors" = "#a9b5e4",
  "Socioeconomic factors" = "#d6c1e0",
  "Genetics" = "#a3c6e5",
  "Gut microbiota" = "#f0b1b5"
)

# Microbiota-versus-genetics dominance classification. A normalized dominance
# score of +/-0.20 corresponds to an approximately 1.5-fold contribution ratio.
# Metabolites with a combined microbiota + genetics contribution below one
# percentage point are kept separate instead of being labelled "comparable".
omics_low_cutoff_percent <- 1
dominance_score_cutoff <- 0.20
dominance_levels <- c(
  "Microbiota-dominant",
  "Comparable",
  "Genetics-dominant",
  "Both low"
)

# Data checks and metabolite annotation
required_individual_columns <- c(
  "period", "metabolite", "group",
  "delta_adjusted_marginal_r2",
  "delta_adjusted_marginal_r2_nonnegative"
)
required_matrix_columns <- c("period", "group", "adjusted_r2", "pvalue")


metabolite_lookup <- info_metab %>%
  transmute(
    metabolite = as.character(Index),
    metabolite_name = coalesce(
      na_if(as.character(sugg_cmpd_name), ""),
      na_if(as.character(Compounds), ""),
      metabolite
    ),
    metab_class = if_else(
      !is.na(class_shrt) & class_shrt %in% names(metab_class_col),
      as.character(class_shrt),
      "Others"
    )
  ) %>%
  distinct(metabolite, .keep_all = TRUE)


# Main-panel data: weighted longitudinal mixed models 
individual_plot_data <- individual_results %>%
  filter(period == "T1-T3", status == "ok") %>%
  mutate(
    group_label = recode(group, !!!group_labels),
    group_label = factor(group_label, levels = group_levels),
    incremental_r2_raw = delta_adjusted_marginal_r2,
    incremental_r2_plot = pmax(
      delta_adjusted_marginal_r2_nonnegative, 0
    ),
    incremental_r2_percent = 100 * incremental_r2_plot
  ) %>%
  left_join(metabolite_lookup, by = "metabolite") %>%
  mutate(
    metabolite_name = coalesce(metabolite_name, metabolite),
    metab_class = coalesce(metab_class, "Others")
  )

omics_dominance <- individual_plot_data %>%
  filter(group %in% c("microbiota", "genetics")) %>%
  select(metabolite, group, incremental_r2_percent) %>%
  pivot_wider(
    names_from = group,
    values_from = incremental_r2_percent,
    values_fill = 0
  ) %>%
  mutate(
    omics_r2_percent = microbiota + genetics,
    dominance_score = if_else(
      omics_r2_percent > 0,
      (microbiota - genetics) / omics_r2_percent,
      0
    ),
    dominance_class = case_when(
      omics_r2_percent < omics_low_cutoff_percent ~ "Both low",
      dominance_score > dominance_score_cutoff ~ "Microbiota-dominant",
      dominance_score < -dominance_score_cutoff ~ "Genetics-dominant",
      TRUE ~ "Comparable"
    ),
    dominance_class = factor(dominance_class, levels = dominance_levels)
  )

metabolite_order <- individual_plot_data %>%
  group_by(metabolite, metabolite_name, metab_class) %>%
  summarise(
    total_r2_percent = sum(incremental_r2_percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(omics_dominance, by = "metabolite") %>%
  arrange(dominance_class, desc(total_r2_percent), metabolite_name, metabolite) %>%
  mutate(
    plot_position = row_number(),
    duplicated_name = duplicated(metabolite_name) |
      duplicated(metabolite_name, fromLast = TRUE),
    metabolite_label = if_else(
      duplicated_name,
      paste0(metabolite_name, " [", metabolite, "]"),
      metabolite_name
    )
  )

dominance_blocks <- metabolite_order %>%
  group_by(dominance_class) %>%
  summarise(
    n_metabolites = n(),
    xmin = min(plot_position) - 0.40,
    xmax = max(plot_position) + 0.40,
    xmid = mean(c(xmin, xmax)),
    .groups = "drop"
  ) %>%
  arrange(dominance_class) %>%
  mutate(
    block_label = paste0(as.character(dominance_class), " (n=", n_metabolites, ")")
  )

dominance_separators <- dominance_blocks %>%
  filter(row_number() < n()) %>%
  transmute(
    xintercept = xmax + 0.10,
    separator_ymax = max(metabolite_order$total_r2_percent) * 0.38
  )

individual_plot_data <- individual_plot_data %>%
  left_join(
    metabolite_order %>%
      select(
        metabolite, metabolite_label, total_r2_percent,
        microbiota, genetics, omics_r2_percent,
        dominance_score, dominance_class, plot_position
      ),
    by = "metabolite"
  ) %>%
  mutate(
    metabolite_label = factor(
      metabolite_label,
      levels = metabolite_order$metabolite_label
    )
  )

class_strip_data <- metabolite_order %>%
  mutate(
    metabolite_label = factor(
      metabolite_label,
      levels = metabolite_order$metabolite_label
    ),
    metab_class = factor(metab_class, levels = names(metab_class_col))
  )

active_classes <- names(metab_class_col)[
  names(metab_class_col) %in% unique(as.character(class_strip_data$metab_class))
]


# Inset data: period-specific matrix dbRDA 
matrix_plot_data <- matrix_results %>%
  filter(period %in% period_levels) %>%
  mutate(
    period = factor(period, levels = period_levels),
    group_label = recode(group, !!!group_labels),
    adjusted_r2_raw = adjusted_r2,
    adjusted_r2_percent = 100 * pmax(adjusted_r2, 0),
    value_label = number(adjusted_r2_percent, accuracy = 0.01)
  )

matrix_group_order <- matrix_plot_data %>%
  group_by(group_label) %>%
  summarise(
    mean_adjusted_r2_percent = mean(adjusted_r2_percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_adjusted_r2_percent), group_label) %>%
  pull(group_label)

matrix_plot_data <- matrix_plot_data %>%
  mutate(group_label = factor(group_label, levels = matrix_group_order))


# Matrix-level dbRDA inset 
matrix_label_offset <- max(matrix_plot_data$adjusted_r2_percent, na.rm = TRUE) * 0.02

fig3g_matrix <- ggplot(
  matrix_plot_data,
  aes(x = adjusted_r2_percent, y = group_label, fill = group_label)
) +
  geom_col(width = 0.68, color = "#333333", linewidth = 0.2) +
  geom_text(
    aes(x = adjusted_r2_percent + matrix_label_offset, label = value_label),
    size = 3.2, hjust = 0
  ) +
  facet_grid(cols = vars(period)) +
  scale_fill_manual(values = group_colors, breaks = group_levels, drop = FALSE) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.16)),
    labels = label_number(accuracy = 0.1)
  ) +
  labs(
    title = "Entire set of 40 metabolites",
    x = expression("Adjusted "*R^2*" (%)"),
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 7) +
  theme(
    plot.background = element_blank(),
    plot.title = element_text(size = 11,  hjust = 0.5, margin = margin(b = 2)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 11),
    axis.title.x = element_text(size = 9),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.spacing.x = unit(0, "pt"),
    legend.position = "none"
  )


# Main stacked bar and metabolite-class strip 
fig3g_bar <- ggplot(
  individual_plot_data,
  aes(x = metabolite_label, y = incremental_r2_percent, fill = group_label)
) +
  geom_col(width = 0.80, color = NA) +
  geom_segment(
    data = dominance_separators,
    aes(
      x = xintercept, xend = xintercept,
      y = 0, yend = separator_ymax
    ),
    inherit.aes = FALSE,
    color = "#555555",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  scale_fill_manual(values = group_colors, breaks = group_levels, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    labels = label_number(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = expression(Delta*" adjusted marginal "*R^2*" (%)"),
    fill = NULL, tag = 'g'
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = c(0.08,0.95),
    legend.justification = "left",
    legend.background = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(5, "mm"),
    axis.ticks.length.x = unit(0, "mm"),
    plot.tag.position = c(0.02, 0.98),
    plot.tag = element_text(size = 18),
    plot.margin = margin(t = 5, r = 5, b = -2, l = 5)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# The decreasing bar heights create an open upper-right region. Place the
# matrix dbRDA panel there without changing the common metabolite x scale.
fig3g_bar_inset <- fig3g_bar +
  inset_element(
    fig3g_matrix,
    left = 0.42, bottom = 0.33, right = 0.99, top = 0.995,
    align_to = "panel",
    clip = FALSE
  )

fig3g_class_strip <- ggplot(
  class_strip_data,
  aes(x = metabolite_label, y = 0.57, fill = metab_class)
) +
  geom_tile(height = 0.14) +
  geom_segment(
    data = dominance_blocks,
    aes(x = xmin, xend = xmax, y = 0.7, yend = 0.7),
    inherit.aes = FALSE,
    color = "#555555",
    linewidth = 0.35
  ) +
  geom_text(
    data = dominance_blocks,
    aes(x = xmid, y = 0.8, label = block_label),
    inherit.aes = FALSE,
    size = 2.8,
    color = "#333333"
  ) +
  scale_fill_manual(
    values = metab_class_col,
    breaks = active_classes,
    drop = FALSE
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.495, 0.90), clip = "off") +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5, size = 9, color = "black"
    ),
    axis.line.x = element_line(color = "black", linewidth = 0.45),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(2.2, "pt"),
    legend.position = "none",
    plot.margin = margin(t = -2, r = 5, b = 4, l = 5)
  )

fig3g <- fig3g_bar_inset / fig3g_class_strip +
  plot_layout(heights = c(16, 1.4))

# Output 
ggsave(
  file.path(figure_dir, "fig3gv2_metabolites_determinant.pdf"),
  plot = fig3g, width = 10, height = 6
)
