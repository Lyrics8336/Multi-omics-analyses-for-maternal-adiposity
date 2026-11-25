####################################################
# Visualization and statistical testing for Figure 4
#####################################################

# Fig 4a: circular heatmap of OWO-associated gut microbiota ------------------------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(MetBrewer)
})

result_file <- "result/PART3/1_diff_microbe.rdata"
abundance_file <- "data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData"
output_dir <- "figure/PART3/graphlan_onlyowo"
fdr_cutoff <- 0.05
period_order <- c("T1", "T2", "T3", "Pooled")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
load(result_file)

stopifnot(exists("res_pbmi_microbe_categorical"))

# Remove rank prefixes from every taxonomic level (d__, p__, c__, ...), not
# only from genus labels. Dots are replaced because GraPhlAn uses them as the
# hierarchy delimiter.
clean_taxon <- function(x) {
  x %>%
    str_remove("^[A-Za-z]__") %>%
    str_replace_all("\\.", "_") %>%
    str_trim()
}

write_graphlan <- function(x, filename) {
  write.table(
    x, file.path(output_dir, filename), sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE, na = ""
  )
}

diff_16s <- res_pbmi_microbe_categorical %>%
  filter(
    omics_type == "16S_g",
    exposure %in% c("OWO vs normal weight"),
    period %in% period_order,
    status == "ok",
    !is.na(OTU), !is.na(name)
  ) %>%
  mutate(
    group = recode( exposure, "OWO vs normal weight" = "OWO" ),
    genus = clean_taxon(name)
  ) %>%
  filter(
    genus != "", genus != "uncultured",
    !str_detect(genus, regex("uncultured|unclassified", ignore_case = TRUE)),
    !str_detect(genus, "^[gG]?[0-9]+$")
  )

# Select the union of genera significant in at least one time-specific or
# pooled model for OWO.
sig <- diff_16s %>% filter(!is.na(p_fdr), p_fdr < fdr_cutoff)
selected_features <- unique(sig$feature)
if (!length(selected_features)) stop("No significant 16S genera at FDR < ", fdr_cutoff)

plot_data <- diff_16s %>% filter(feature %in% selected_features)

## One row per selected feature and a clean five-level bacterial lineage
taxa <- plot_data %>%
  distinct(feature, OTU, genus) %>%
  separate(
    OTU, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus_otu"),
    sep = ";", fill = "right", extra = "drop", remove = FALSE
  ) %>%
  mutate(across(Kingdom:Genus_otu, clean_taxon), Genus = genus) %>%
  filter(!is.na(Phylum), !is.na(Genus), Genus != "") %>%
  distinct(feature, .keep_all = TRUE)

write.table(
    taxa %>% select(Phylum, Class, Order, Family, Genus), 
    file.path(output_dir, "tree_backbone.txt"), sep = ".", quote = FALSE,
    col.names = FALSE, row.names = FALSE, na = ""
  )

# beta rings
beta_limit <- range(plot_data$beta)
palette_21 <- rev( met.brewer("Cassatt1", n = 21) )
breaks1 <- quantile(plot_data$beta[which(plot_data$beta < 0)], 
                    probs = seq(0, 1, length.out = 11), na.rm = TRUE)
breaks2 <- quantile(plot_data$beta[which(plot_data$beta > 0)], 
                    probs = seq(0, 1, length.out = 11), na.rm = TRUE)
breaks2[1] = 0
beta_breaks <- c(breaks1, breaks2)

cell_color <- function(beta, p_fdr) {
  out <- rep("#FFFFFF", length(beta))
  keep <- !is.na(beta) & !is.na(p_fdr) & p_fdr < fdr_cutoff
  bin <- findInterval(beta[keep], beta_breaks, all.inside = TRUE)
  out[keep] <- palette_21[pmin(bin, 21)]
  out
}

make_ring <- function(group_name, period_name, ring_number) {
  cells <- plot_data %>%
    filter(group == group_name, period == period_name) %>%
    transmute(Genus = genus, color = cell_color(beta, p_fdr)) %>%
    right_join(taxa %>% select(Genus), by = "Genus") %>%
    mutate(color = replace_na(color, "#FFFFFF"))

  bind_rows(
    cells %>% transmute(Genus, type = "ring_color", ring_number, value = color),
    cells %>% transmute(Genus, type = "ring_edge_color", ring_number,
                        value = "#000000")
  )
}

# Inner heatmap: OWO T1/T2/T3/Pooled (rings 1-4)
walk2(period_order, 1:4, ~ write_graphlan(make_ring("OWO", .x, .y),
                                          paste0("ring", .y, "_owo_", .x, ".txt")))


# Ring 5: Prevalence across all samples
load(abundance_file)
stopifnot(exists("bact_fung_tab"))
abund_tab <- bact_fung_tab %>% select(matches("^bact_g\\d+$"))
height <- colMeans(abund_tab > 0)[taxa$feature]

prevalence_ring <- tibble(feature = names(height), prevalence = as.numeric(height)*1.5) %>%
  left_join(taxa %>% select(feature, Genus), by = "feature") %>%
  transmute(Genus, type = "ring_height", ring = 5, value = prevalence)
write_graphlan(prevalence_ring, "ring5_prevalence.txt")

load("result/PART3/thsbc_hbc_rep_genus_species_owo_sig_sets.rdata")   
owo_rep_name$rep_genus
rep_genus <- owo_rep_name$rep_genus
length(rep_genus)

point_color <- data.frame(
  V1 = str_remove(rep_genus, 'g__'),
  V2 = 'clade_marker_color') %>% 
  left_join(taxa %>% select(Genus, Phylum), by = c('V1' = 'Genus')) %>% 
  mutate(V3 = case_when(
    Phylum == 'Actinobacteriota' ~ '#6A4099',
    Phylum == 'Bacteroidota' ~ '#F4B336',
    Phylum == 'Firmicutes' ~ '#37AB5A',
    Phylum == 'Proteobacteria' ~ '#88a4c9',
    TRUE ~ '#D3D3D3'
  )) %>% select(V1,V2,V3)

seq <- c(
  'Mediterraneibacter', "Dorea", "Blautia", "Eisenbergiella", "Hungatella", 'Fusicatenibacter', 
  'Anaeromassilibacillus', "Pseudoruminococcus", 'Vescimonas', 'Intestinimonas', "Streptococcus", "Granulicatella", 
  "Weissella", "Megasphaera", "Dialister", "Megamonas", "Phascolarctobacterium", 
  "Odoribacter", "Alistipes", "Prevotella", "Olsenella", "Raoultibacter", "Fusobacterium"
)
point_LETTER <- data.frame(
  V1 = seq,
  V2 = 'annotation',
  V3 = LETTERS[1:23]) %>% 
  mutate(V3 = str_c(V3, ":", V1))
point_LETTER_bg <- data.frame(
  V1 = point_LETTER$V1,
  V2 = "annotation_background_color",
  V3 = "#D3D3D3"
)

point_color %>% rbind(point_LETTER) %>% rbind(point_LETTER_bg) %>%
  write.table(file.path(output_dir, "track0.txt"),  sep = "\t", quote = F,col.names = F,row.names = F, na="")


#####
# cat global_anno.txt track0.txt > track0_all
# graphlan_annotate.py --annot track0_all tree_backbone.txt graphlan.xml
# graphlan_annotate.py --annot ring1_owo_T1.txt graphlan.xml graphlan1.xml
# graphlan_annotate.py --annot ring2_owo_T2.txt graphlan1.xml graphlan2.xml
# graphlan_annotate.py --annot ring3_owo_T3.txt graphlan2.xml graphlan3.xml
# graphlan_annotate.py --annot ring4_owo_Pooled.txt graphlan3.xml graphlan4.xml
# graphlan_annotate.py --annot ring5_prevalence.txt graphlan4.xml graphlan5.xml
# graphlan.py graphlan5.xml fig4a_micro_circle_onlyowo.pdf --dpi 300 --size 4 --pad 0.0


# Fig 4b: bacteria MetaG replication in THSBC and HBC ------------------------------------
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(patchwork)

library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(MetBrewer)

load('data/01_thsbc_cleaned_data/gut_microbe_16s_its_3688_815.RData')
info_bact_fung_thsbc = info_bact_fung
load('data/01_thsbc_cleaned_data/gut_microbe_MetaG_1250_473.RData')
info_mgs_thsbc = info_mgs
load('data/02_hbc_cleaned_data/gut_microbe_MetaG_393_158.RData')
info_mgs_hbc = info_mgs

load("result/PART3/1_diff_microbe.rdata")
load("result/PART3/2_diff_microb_HBC.rdata")

load("result/PART3/1_diff_16s_its_sig_sets_thsbc.rdata")
load("result/PART3/genus_match_mgs_16s.RData")


## Bacterial MetaG validation in THSBC and HBC 
period_order <- c("T1", "T2", "T3", "Pooled")
target_exposure <- "OWO vs normal weight"
cohort_order <- c("THSBC", "HBC")

# THSBC stores MetaG IDs as `mgs_g*` / `mgs_s*`, whereas the HBC metadata
# stores the same reference IDs as `g*` / `s*` and adds the `mgs_` prefix
# when model results are assembled.  Build an explicit cross-cohort taxonomy
# audit so that a shared feature ID is accepted only when both its displayed
# name and complete lineage are identical in the two cohort dictionaries.
standardise_hbc_mgs_id <- function(id) {
  if_else(str_starts(id, "mgs_"), id, paste0("mgs_", id))
}

build_cross_cohort_taxonomy_audit <- function(taxonomic_level) {
  thsbc_taxonomy <- info_mgs_thsbc %>%
    filter(level == taxonomic_level) %>%
    transmute(
      feature = ID,
      name_THSBC = name,
      OTU_THSBC = OTU
    )

  hbc_taxonomy <- info_mgs_hbc %>%
    filter(level == taxonomic_level) %>%
    transmute(
      feature = standardise_hbc_mgs_id(ID),
      name_HBC = name,
      OTU_HBC = OTU
    )

  if (anyDuplicated(thsbc_taxonomy$feature) ||
      anyDuplicated(hbc_taxonomy$feature)) {
    stop("Duplicated ", taxonomic_level,
         " feature IDs were found in a cohort taxonomy dictionary.")
  }

  full_join(thsbc_taxonomy, hbc_taxonomy, by = "feature") %>%
    mutate(
      present_THSBC = !is.na(name_THSBC) & !is.na(OTU_THSBC),
      present_HBC = !is.na(name_HBC) & !is.na(OTU_HBC),
      present_both = present_THSBC & present_HBC,
      name_identical = present_both & name_THSBC == name_HBC,
      OTU_identical = present_both & OTU_THSBC == OTU_HBC,
      taxonomy_match = name_identical & OTU_identical
    )
}

genus_cross_cohort_taxonomy <-
  build_cross_cohort_taxonomy_audit("Genus")
species_cross_cohort_taxonomy <-
  build_cross_cohort_taxonomy_audit("Species")

# Require every taxon shown in the heatmaps to have an estimable beta and
# P value for every displayed period in both metagenomic cohorts. Statistical
# support may still come from either cohort alone.
summarise_cross_cohort_coverage <- function(
    evidence,
    feature_ids,
    feature_column
) {
  observed <- evidence %>%
    transmute(
      feature_key = .data[[feature_column]],
      cohort,
      period,
      usable = !is.na(beta) & !is.na(p_value)
    ) %>%
    filter(usable) %>%
    distinct(feature_key, cohort, period) %>%
    mutate(usable = TRUE)

  expand_grid(
    feature_key = feature_ids,
    cohort = cohort_order,
    period = period_order
  ) %>%
    left_join(observed, by = c("feature_key", "cohort", "period")) %>%
    mutate(usable = replace_na(usable, FALSE)) %>%
    group_by(feature_key, cohort) %>%
    summarise(
      n_periods = sum(usable),
      available = any(usable),
      complete = all(usable),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = cohort,
      values_from = c(n_periods, available, complete),
      names_sep = "_"
    ) %>%
    mutate(
      available_both = available_THSBC & available_HBC,
      complete_both = complete_THSBC & complete_HBC
    )
}

# genus_match_mgs_16s is the sole authority for linking the two taxonomies.
# This allows reliable SILVA-NCBI cross-name matches to enter the validation
# analysis while excluding placeholder or unresolved 16S genera.
required_mapping_columns <- c(
  "ID_16S", "ID_MetaG", "OTU_16S", "OTU_MetaG",
  "name_16S_original", "name_MetaG_original", "common_name",
  "matched", "match_rule", "match_evidence", "name_sim"
)

owo_genus_mapping <- genus_match_mgs_16s %>%
  filter(
    matched,
    ID_16S %in% sig_16s_sets$sigid_16s_owo
  ) %>%
  select(all_of(required_mapping_columns)) %>%
  arrange(common_name)

genus_taxonomy_audit <- owo_genus_mapping %>%
  transmute(
    ID_16S,
    feature = ID_MetaG,
    expected_name = name_MetaG_original,
    expected_OTU = OTU_MetaG,
    common_name
  ) %>%
  left_join(genus_cross_cohort_taxonomy, by = "feature") %>%
  mutate(
    expected_taxonomy_match =
      taxonomy_match &
      expected_name == name_THSBC &
      expected_name == name_HBC &
      expected_OTU == OTU_THSBC &
      expected_OTU == OTU_HBC
  )

if (anyNA(genus_taxonomy_audit$expected_taxonomy_match) ||
    any(!genus_taxonomy_audit$expected_taxonomy_match)) {
  failed_genus <- genus_taxonomy_audit %>%
    filter(is.na(expected_taxonomy_match) | !expected_taxonomy_match) %>%
    pull(feature)
  stop(
    "Cross-cohort genus taxonomy mismatch for: ",
    paste(failed_genus, collapse = ", ")
  )
}

message(
  "Cross-cohort taxonomy audit: all ", nrow(genus_taxonomy_audit),
  " mapped OWO-related genera have identical IDs, names and lineages."
)

# Define the discovery direction from the T1 THSBC 16S association estimate
discovery_direction <- res_pbmi_microbe_categorical %>%
  filter(
    omics_type == "16S_g",
    exposure == target_exposure,
    period == "T1",
    status == "ok",
    feature %in% owo_genus_mapping$ID_16S
  ) %>%
  transmute(
    ID_16S = feature,
    beta_16S_T1 = beta,
    direction = case_when(
      beta > 0 ~ "Enriched in OWO",
      beta < 0 ~ "Depleted in OWO",
      TRUE ~ NA_character_
    )
  )

owo_genus_mapping <- owo_genus_mapping %>%
  left_join(discovery_direction, by = "ID_16S")

collect_genus_evidence <- function(results, cohort_name) {
  results %>%
    filter(
      omics_type == "MetaG_g",
      exposure == target_exposure,
      period %in% period_order,
      status == "ok",
      !is.na(beta),
      !is.na(p_value),
      feature %in% owo_genus_mapping$ID_MetaG
    ) %>%
    select(feature, period, beta, p_value) %>%
    inner_join(
      owo_genus_mapping,
      by = c("feature" = "ID_MetaG")
    ) %>%
    mutate(
      cohort = cohort_name,
      direction_concordant = beta * beta_16S_T1 > 0,
      concordant_nominal = p_value < 0.05 & direction_concordant,
      discordant_nominal = p_value < 0.05 & !direction_concordant
    )
}

genus_evidence <- bind_rows(
  collect_genus_evidence(res_pbmi_microbe_categorical, "THSBC"),
  collect_genus_evidence(res_hbc_microbe_pbmi_categorical, "HBC")
)

genus_coverage <- summarise_cross_cohort_coverage(
  evidence = genus_evidence,
  feature_ids = owo_genus_mapping$ID_16S,
  feature_column = "ID_16S"
) %>%
  rename(ID_16S = feature_key)

# A genus is displayed only when both cohorts have complete data for all four
# periods. Metagenomic support still requires a directionally concordant
# nominal association in either THSBC or HBC, not significance in both.
genus_validation <- owo_genus_mapping %>%
  select(
    ID_16S, ID_MetaG, name_16S_original, name_MetaG_original,
    common_name, OTU_16S, OTU_MetaG, match_rule, match_evidence,
    name_sim, beta_16S_T1, direction
  ) %>%
  left_join(
    genus_evidence %>%
      group_by(ID_16S) %>%
      summarise(
        has_concordant_nominal = any(concordant_nominal, na.rm = TRUE),
        has_discordant_nominal = any(discordant_nominal, na.rm = TRUE),
        supporting_cohorts = paste(
          sort(unique(cohort[concordant_nominal])),
          collapse = "+"
        ),
        .groups = "drop"
      ),
    by = "ID_16S"
  ) %>%
  left_join(genus_coverage, by = "ID_16S") %>%
  mutate(
    across(
      c(has_concordant_nominal, has_discordant_nominal),
      ~ replace_na(.x, FALSE)
    ),
    supported =
      complete_both &
      has_concordant_nominal &
      !has_discordant_nominal,
    supporting_cohorts = na_if(supporting_cohorts, "")
  )

validated_genus_map <- genus_validation %>%
  filter(supported) %>%
  arrange(
    factor(direction, levels = c("Depleted in OWO", "Enriched in OWO")),
    common_name
  )

if (nrow(validated_genus_map) == 0L ||
    any(!validated_genus_map$complete_both)) {
  stop("Every displayed genus must have complete data in THSBC and HBC.")
}

message(
  "Mapped OWO-related genera: ", nrow(owo_genus_mapping),
  "; metagenomically supported genera: ", nrow(validated_genus_map),
  " (", sum(validated_genus_map$direction == "Enriched in OWO"),
  " enriched, ",
  sum(validated_genus_map$direction == "Depleted in OWO"), " depleted)."
)

build_genus_wide <- function(results) {
  results %>%
    filter(
      omics_type == "MetaG_g",
      exposure == target_exposure,
      period %in% period_order,
      status == "ok",
      feature %in% validated_genus_map$ID_MetaG
    ) %>%
    select(ID_MetaG = feature, period, beta, p_value) %>%
    inner_join(
      validated_genus_map %>%
        select(ID_MetaG, taxon = common_name),
      by = "ID_MetaG"
    ) %>%
    select(taxon, period, beta, p_value) %>%
    pivot_wider(
      id_cols = taxon,
      names_from = period,
      values_from = c(beta, p_value),
      names_sep = "_"
    )
}

pdata_thsbc_mgs_g <- build_genus_wide(res_pbmi_microbe_categorical)
pdata_hbc_mgs_g <- build_genus_wide(res_hbc_microbe_pbmi_categorical)

row_info <- validated_genus_map %>%
  transmute(
    taxon = common_name,
    OTU = OTU_MetaG,
    direction,
    phylum = str_extract(OTU_MetaG, "(?<=p__)[^|;]+"),
    phylum = recode(
      phylum,
      "Actinobacteria" = "Actinobacteriota",
      "Bacteroidetes" = "Bacteroidota"
    ),
    phylum = if_else(
      phylum %in% c(
        "Actinobacteriota", "Bacteroidota", "Firmicutes",
        "Proteobacteria"
      ),
      phylum,
      "Others"
    )
  ) %>%
  mutate(
    direction = factor(
      direction,
      levels = c("Depleted in OWO", "Enriched in OWO")
    )
  ) %>%
  arrange(direction, taxon)

stopifnot(nrow(row_info) > 0L, !anyDuplicated(row_info$taxon))

make_study_matrix <- function(data, value_prefix, row_names) {
  wanted <- paste0(value_prefix, "_", period_order)
  missing_columns <- setdiff(wanted, names(data))
  if (length(missing_columns) > 0L) data[missing_columns] <- NA_real_

  data %>%
    select(taxon, all_of(wanted)) %>%
    right_join(tibble(taxon = row_names), by = "taxon") %>%
    slice(match(row_names, taxon)) %>%
    column_to_rownames("taxon") %>%
    as.matrix()
}

htbeta <- cbind(
  make_study_matrix(pdata_thsbc_mgs_g, "beta", row_info$taxon),
  make_study_matrix(pdata_hbc_mgs_g, "beta", row_info$taxon)
)
htpv <- cbind(
  make_study_matrix(pdata_thsbc_mgs_g, "p_value", row_info$taxon),
  make_study_matrix(pdata_hbc_mgs_g, "p_value", row_info$taxon)
)
colnames(htbeta) <- colnames(htpv) <- rep(period_order, 2)

if (anyNA(htbeta) || anyNA(htpv)) {
  stop("Missing genus-level beta or P values remain in the plotted matrix.")
}

beta_limit <- max(abs(htbeta), na.rm = TRUE)
if (!is.finite(beta_limit) || beta_limit == 0) beta_limit <- 1
col_fun <- colorRamp2(c(-beta_limit, 0, beta_limit),
                      c("#5A5A83", "white", "#B1615C"))

row_anno <- rowAnnotation(
  Direction = row_info$direction,
  Phylum = row_info$phylum,
  col = list(
    Direction = c("Enriched in OWO" = "#F4C2C2",
                  "Depleted in OWO" = "#A3C4E0"),
    Phylum = c("Actinobacteriota" = "#6B6CA3",
               "Bacteroidota" = "#D9A34D",
               "Firmicutes" = "#80A46E",
               "Proteobacteria" = "#87BCBD",
               "Others" = "#D3D3D3")
  ),
  show_annotation_name = FALSE,
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(
    Direction = list(title = NULL),
    Phylum = list(title = "Phylum")
  )
)

dir.create("figure/PART3", recursive = TRUE, showWarnings = FALSE)
pdf("figure/PART3/fig4b_thsbc_hbc_bact_MetaG_heatmap.pdf",
    width = 7, height = max(5, 0.28 * nrow(row_info) + 2.2))
ht <- Heatmap(
  htbeta,
  name = "Beta",
  col = col_fun,
  na_col = "#F2F2F2",
  left_annotation = row_anno,
  column_split = factor(rep(c("THSBC", "HBC"), each = length(period_order)),
                        levels = c("THSBC", "HBC")),
  column_labels = rep(period_order, 2),
  column_names_side = "bottom",
  column_names_rot = 30,
  column_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 15),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_labels = str_remove(rownames(htbeta), "^g__"),
  row_names_gp = gpar(fontsize = 11, fontface = "italic"),
  row_names_max_width = unit(7, "cm"),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1.5),
  border_gp = gpar(col = "black", lwd = 2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(htpv[i, j])) {
      if (htpv[i, j] < 0.01) {
        grid.text("*", x, y, vjust = 0.8, gp = gpar(fontsize = 20))
      } else if (htpv[i, j] < 0.05) {
        grid.text("+", x, y, gp = gpar(fontsize = 14))
      }
    }
  },
  heatmap_legend_param = list(
    title = "Beta",
    at = c(-beta_limit, 0, beta_limit),
    labels = round(c(-beta_limit, 0, beta_limit), 2)
  )
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legends = FALSE)
grid.text("+  P < 0.05     *  P < 0.01",
          x = unit(0.99, "npc"), y = unit(0.01, "npc"),
          just = c("right", "bottom"), gp = gpar(fontsize = 9))
dev.off()


## Species-level MetaG heatmap and vertical combination 
# Restrict species to descendants of the metagenomically supported genera
species_taxonomy_all <- info_mgs_thsbc %>%
  filter(level == "Species") %>%
  transmute(
    ID_species = ID,
    species_name = name,
    OTU,
    genus_metag = paste0("g__", str_extract(OTU, "(?<=g__)[^|;]+"))
  ) %>%
  inner_join(
    validated_genus_map %>%
      select(
        ID_16S, ID_MetaG, name_16S_original, name_MetaG_original,
        parent_common_name = common_name, direction
    ),
    by = c("genus_metag" = "name_MetaG_original")
  )

species_taxonomy_audit <- species_taxonomy_all %>%
  left_join(
    species_cross_cohort_taxonomy,
    by = c("ID_species" = "feature")
  ) %>%
  mutate(
    expected_taxonomy_match =
      taxonomy_match &
      species_name == name_THSBC &
      species_name == name_HBC &
      OTU == OTU_THSBC &
      OTU == OTU_HBC
  )

# Only exact cross-cohort species mappings are eligible for validation.  The
# later coverage rule additionally requires estimates at every displayed
# period in both cohorts.
species_taxonomy <- species_taxonomy_audit %>%
  filter(!is.na(expected_taxonomy_match), expected_taxonomy_match) %>%
  select(all_of(names(species_taxonomy_all)))

message(
  "Cross-cohort taxonomy audit: ", nrow(species_taxonomy), " of ",
  nrow(species_taxonomy_all),
  " descendant species have identical IDs, names and lineages."
)

collect_species_evidence <- function(results, cohort_name) {
  results %>%
    filter(
      omics_type == "MetaG_s",
      exposure == target_exposure,
      period %in% period_order,
      status == "ok",
      !is.na(beta),
      !is.na(p_value),
      feature %in% species_taxonomy$ID_species
    ) %>%
    select(feature, period, beta, p_value) %>%
    inner_join(
      species_taxonomy,
      by = c("feature" = "ID_species")
    ) %>%
    mutate(
      cohort = cohort_name,
      direction_concordant = case_when(
        direction == "Enriched in OWO" ~ beta > 0,
        direction == "Depleted in OWO" ~ beta < 0,
        TRUE ~ FALSE
      ),
      concordant_nominal = p_value < 0.05 & direction_concordant,
      discordant_nominal = p_value < 0.05 & !direction_concordant
    )
}

species_evidence <- bind_rows(
  collect_species_evidence(res_pbmi_microbe_categorical, "THSBC"),
  collect_species_evidence(res_hbc_microbe_pbmi_categorical, "HBC")
)

species_coverage <- summarise_cross_cohort_coverage(
  evidence = species_evidence,
  feature_ids = species_taxonomy$ID_species,
  feature_column = "feature"
) %>%
  rename(ID_species = feature_key)

species_validation <- species_taxonomy %>%
  left_join(
    species_evidence %>%
      group_by(feature) %>%
      summarise(
        has_concordant_nominal = any(concordant_nominal, na.rm = TRUE),
        has_discordant_nominal = any(discordant_nominal, na.rm = TRUE),
        supporting_cohorts = paste(
          sort(unique(cohort[concordant_nominal])),
          collapse = "+"
        ),
        .groups = "drop"
      ),
    by = c("ID_species" = "feature")
  ) %>%
  left_join(species_coverage, by = "ID_species") %>%
  mutate(
    across(
      c(has_concordant_nominal, has_discordant_nominal),
      ~ replace_na(.x, FALSE)
    ),
    supported =
      complete_both &
      has_concordant_nominal &
      !has_discordant_nominal,
    supporting_cohorts = na_if(supporting_cohorts, "")
  )

validated_species_map <- species_validation %>%
  filter(supported) %>%
  arrange(
    factor(direction, levels = c("Depleted in OWO", "Enriched in OWO")),
    species_name
  )

if (nrow(validated_species_map) == 0L ||
    any(!validated_species_map$complete_both)) {
  stop("Every displayed species must have complete data in THSBC and HBC.")
}

if (any(!validated_species_map$ID_species %in%
        species_taxonomy$ID_species)) {
  stop("A displayed species failed the cross-cohort taxonomy audit.")
}

message(
  "Metagenomically supported species: ", nrow(validated_species_map),
  " (", sum(validated_species_map$direction == "Enriched in OWO"),
  " enriched, ",
  sum(validated_species_map$direction == "Depleted in OWO"), " depleted)."
)

build_species_wide <- function(results) {
  results %>%
    filter(
      omics_type == "MetaG_s",
      exposure == target_exposure,
      period %in% period_order,
      status == "ok",
      feature %in% validated_species_map$ID_species
    ) %>%
    select(ID_species = feature, period, beta, p_value) %>%
    inner_join(
      validated_species_map %>%
        select(ID_species, taxon = species_name),
      by = "ID_species"
    ) %>%
    select(taxon, period, beta, p_value) %>%
    pivot_wider(
      id_cols = taxon,
      names_from = period,
      values_from = c(beta, p_value),
      names_sep = "_"
    )
}

pdata_thsbc_mgs_s <- build_species_wide(res_pbmi_microbe_categorical)
pdata_hbc_mgs_s <- build_species_wide(res_hbc_microbe_pbmi_categorical)

species_row_info <- validated_species_map %>%
  transmute(
    taxon = species_name,
    OTU,
    genus = genus_metag,
    direction,
    phylum = str_extract(OTU, "(?<=p__)[^|;]+"),
    phylum = recode(
      phylum,
      "Actinobacteria" = "Actinobacteriota",
      "Bacteroidetes" = "Bacteroidota"
    ),
    phylum = if_else(
      phylum %in% c(
        "Actinobacteriota", "Bacteroidota", "Firmicutes",
        "Proteobacteria"
      ),
      phylum,
      "Others"
    )
  ) %>%
  mutate(
    direction = factor(
      direction,
      levels = c("Depleted in OWO", "Enriched in OWO")
    )
  ) %>%
  arrange(direction, taxon)

stopifnot(
  nrow(species_row_info) > 0L,
  !anyDuplicated(species_row_info$taxon)
)

htbeta_species <- cbind(
  make_study_matrix(pdata_thsbc_mgs_s, "beta", species_row_info$taxon),
  make_study_matrix(pdata_hbc_mgs_s, "beta", species_row_info$taxon)
)
htpv_species <- cbind(
  make_study_matrix(pdata_thsbc_mgs_s, "p_value", species_row_info$taxon),
  make_study_matrix(pdata_hbc_mgs_s, "p_value", species_row_info$taxon)
)
colnames(htbeta_species) <- colnames(htpv_species) <- rep(period_order, 2)

if (anyNA(htbeta_species) || anyNA(htpv_species)) {
  stop("Missing species-level beta or P values remain in the plotted matrix.")
}

combined_beta_limit <- max(abs(c(htbeta, htbeta_species)), na.rm = TRUE)
if (!is.finite(combined_beta_limit) || combined_beta_limit == 0) {
  combined_beta_limit <- 1
}
combined_col_fun <- colorRamp2(
  c(-0.52, 0, 0.72),
  c("#5A5A83", "white", "#B1615C")
)

species_row_anno <- rowAnnotation(
  Direction = species_row_info$direction,
  Phylum = species_row_info$phylum,
  col = list(
    Direction = c(
      "Enriched in OWO" = "#F4C2C2",
      "Depleted in OWO" = "#A3C4E0"
    ),
    Phylum = c(
      "Actinobacteriota" = "#6B6CA3",
      "Bacteroidota" = "#D9A34D",
      "Firmicutes" = "#80A46E",
      "Proteobacteria" = "#87BCBD",
      "Others" = "#D3D3D3"
    )
  ),
  show_annotation_name = FALSE,
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(
    Direction = list(title = NULL),
    Phylum = list(title = "Phylum")
  )
)

ht_species <- Heatmap(
  htbeta_species,
  name = "Beta_species",
  show_heatmap_legend = FALSE,
  col = combined_col_fun,
  na_col = "#F2F2F2",
  left_annotation = species_row_anno,
  column_split = factor(
    rep(c("THSBC", "HBC"), each = length(period_order)),
    levels = c("THSBC", "HBC")
  ),
  column_title = NULL,
  column_labels = rep(period_order, 2),
  column_names_side = "bottom",
  column_names_rot = 30,
  column_names_gp = gpar(fontsize = 12),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_labels = rownames(htbeta_species) %>%
    str_remove("^s__") %>%
    str_replace("_", " "),
  row_names_gp = gpar(fontsize = 12, fontface = "italic"),
  row_names_max_width = unit(7, "cm"),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1.5),
  border_gp = gpar(col = "black", lwd = 2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(htpv_species[i, j])) {
      if (htpv_species[i, j] < 0.01) {
        grid.text("*", x, y, vjust = 0.8, gp = gpar(fontsize = 20))
      } else if (htpv_species[i, j] < 0.05) {
        grid.text("+", x, y, gp = gpar(fontsize = 14))
      }
    }
  },
  heatmap_legend_param = list(
    title = "Beta",
    at = c(-0.52, 0, 0.72),
    labels = round(c(-0.52, 0, 0.72), 2)
  )
)

# Recreate the genus panel with the shared scale and hide its period labels;
# period labels are displayed once below the species panel.
ht_genus_combined <- Heatmap(
  htbeta,
  name = "Beta",
  col = combined_col_fun,
  na_col = "#F2F2F2",
  left_annotation = row_anno,
  column_split = factor(
    rep(c("THSBC", "HBC"), each = length(period_order)),
    levels = c("THSBC", "HBC")
  ),
  column_labels = rep(period_order, 2),
  show_column_names = FALSE,
  column_title_gp = gpar(fontsize = 15),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_labels = str_remove(rownames(htbeta), "^g__"),
  row_names_gp = gpar(fontsize = 12, fontface = "italic"),
  row_names_max_width = unit(7, "cm"),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1.5),
  border_gp = gpar(col = "black", lwd = 2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(htpv[i, j])) {
      if (htpv[i, j] < 0.01) {
        grid.text("*", x, y, vjust = 0.8, gp = gpar(fontsize = 16))
      } else if (htpv[i, j] < 0.05) {
        grid.text("+", x, y, gp = gpar(fontsize = 12))
      }
    }
  },
  heatmap_legend_param = list(
    title = "Beta",
    at = c(-0.52, 0, 0.72),
    labels = round(c(-0.52, 0, 0.72), 2)
  )
)

pdf(
  "figure/PART3/fig4b_thsbc_hbc_bact_MetaG_genus_species_heatmap.pdf",
  width = 6.8,
  height = 14
)
combined_ht <- ht_genus_combined %v% ht_species
draw(
  combined_ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)
grid.text(
  "+  P < 0.05\n*  P < 0.01",
  x = unit(0.95, "npc"), y = unit(0.1, "npc"),
  just = c("right", "bottom"), gp = gpar(fontsize = 12)
)
dev.off()


## save sets 
validated_genus_map <- validated_genus_map %>%
  mutate(
    direction = factor(
      direction,
      levels = c("Depleted in OWO", "Enriched in OWO")
    )
  ) %>%
  arrange(direction, common_name)

validated_species_map <- validated_species_map %>%
  mutate(
    direction = factor(
      direction,
      levels = c("Depleted in OWO", "Enriched in OWO")
    )
  ) %>%
  arrange(direction, species_name)

# Store platform-specific IDs and names together with harmonized display names
owo_rep_name <- list(
  rep_genus = validated_genus_map$name_16S_original,
  rep_genus_16s = validated_genus_map$name_16S_original,
  rep_genus_mgs = validated_genus_map$name_MetaG_original,
  rep_genus_common = validated_genus_map$common_name,
  rep_genus_neg = validated_genus_map %>%
    filter(direction == "Depleted in OWO") %>%
    pull(name_16S_original),
  rep_genus_pos = validated_genus_map %>%
    filter(direction == "Enriched in OWO") %>%
    pull(name_16S_original),
  rep_genus_mgs_neg = validated_genus_map %>%
    filter(direction == "Depleted in OWO") %>%
    pull(name_MetaG_original),
  rep_genus_mgs_pos = validated_genus_map %>%
    filter(direction == "Enriched in OWO") %>%
    pull(name_MetaG_original),
  rep_species = validated_species_map$species_name,
  rep_species_neg = validated_species_map %>%
    filter(direction == "Depleted in OWO") %>%
    pull(species_name),
  rep_species_pos = validated_species_map %>%
    filter(direction == "Enriched in OWO") %>%
    pull(species_name)
)

owo_rep_id <- list(
  rep_genus = validated_genus_map$ID_16S,
  rep_genus_mgs = validated_genus_map$ID_MetaG,
  rep_genus_pos = validated_genus_map %>%
    filter(direction == "Enriched in OWO") %>%
    pull(ID_16S),
  rep_genus_neg = validated_genus_map %>%
    filter(direction == "Depleted in OWO") %>%
    pull(ID_16S),
  rep_genus_mgs_pos = validated_genus_map %>%
    filter(direction == "Enriched in OWO") %>%
    pull(ID_MetaG),
  rep_genus_mgs_neg = validated_genus_map %>%
    filter(direction == "Depleted in OWO") %>%
    pull(ID_MetaG),
  rep_species = validated_species_map$ID_species,
  rep_species_pos = validated_species_map %>%
    filter(direction == "Enriched in OWO") %>%
    pull(ID_species),
  rep_species_neg = validated_species_map %>%
    filter(direction == "Depleted in OWO") %>%
    pull(ID_species)
)

dir.create("result/PART3", recursive = TRUE, showWarnings = FALSE)
save(
  owo_rep_name,
  owo_rep_id,
  validated_genus_map,
  validated_species_map,
  genus_validation,
  species_validation,
  genus_coverage,
  species_coverage,
  genus_taxonomy_audit,
  species_taxonomy_audit,
  file = "result/PART3/thsbc_hbc_rep_genus_species_owo_sig_sets.rdata"
)


# Fig 4c: key 16S genera and clinical traits in THSBC ------------------------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggh4x)
  library(ggnewscale)
})

# Input 
load("result/PART3/3.1_association_micro_clinical_traits_thsbc.RData")
load("result/PART3/thsbc_hbc_rep_genus_species_owo_sig_sets.rdata")
load("result/PART3/genus_match_mgs_16s.RData")

# Display order and annotations 
clinical_traits <- c(
  "DBP", "SBP",
  "FPG", "OGTT_0h", "OGTT_1h", "OGTT_2h",
  "GGT", "ALP", "ALT", "AST", "CHE", "LAC", "LDH",
  "UA",
  "FT3"
)

clinical_labels <- c(
  "DBP", "SBP",
  "FPG", "OGTT 0h", "OGTT 1h", "OGTT 2h",
  "GGT", "ALP", "ALT", "AST", "CHE", "LAC", "LDH",
  "UA",
  "FT3"
)

clinical_group_levels <- c(
  "Blood pressure", "Blood glucose", "Liver", "Renal", "Thyroid"
)

clinical_group_colors <- c(
  "Blood pressure" = "#F8766D",
  "Blood glucose" = "#C4B800",
  "Liver" = "#5AB897",
  "Renal" = "#58B6D9",
  "Thyroid" = "#C995C7"
)

direction_colors <- c(
  "Depleted" = "#A3C4E0",
  "Enriched" = "#F4C2C2"
)

# Use the saved metagenomically supported genus set and verify agreement with the clinical-association results
selected_genus_ids <- unique(owo_rep_id$rep_genus)

stopifnot(
  length(selected_genus_ids) == 23L,
  setequal(selected_genus_ids, unique(owo_16s)),
  all(selected_genus_ids %in% unique(results_16s$microbe))
)

genus_key <- microbiota_annotation %>%
  filter(
    microbiome_type == "16S_g",
    microbe %in% selected_genus_ids
  ) %>%
  transmute(
    microbe,
    name_16S_analysis = microbe_name
  ) %>%
  left_join(
    genus_match_mgs_16s %>%
      transmute(
        microbe = ID_16S,
        name_16S_original,
        common_name,
        matched,
        match_rule
      ) %>%
      distinct(microbe, .keep_all = TRUE),
    by = "microbe"
  ) %>%
  left_join(
    validated_genus_map %>%
      select(microbe = ID_16S, direction) %>%
      distinct(microbe, .keep_all = TRUE),
    by = "microbe"
  ) %>%
  mutate(
    direction = recode(
      direction,
      "Depleted in OWO" = "Depleted",
      "Enriched in OWO" = "Enriched"
    )
  )

missing_genus_ids <- setdiff(selected_genus_ids, genus_key$microbe)
if (length(missing_genus_ids) > 0L) {
  stop(
    "The following current 16S genera are absent from the plotting key: ",
    paste(missing_genus_ids, collapse = ", ")
  )
}
if (any(is.na(genus_key$common_name) | str_trim(genus_key$common_name) == "")) {
  stop(
    "A common 16S genus name is unavailable for: ",
    paste(
      genus_key$microbe[
        is.na(genus_key$common_name) | str_trim(genus_key$common_name) == ""
      ],
      collapse = ", "
    )
  )
}
if (any(is.na(genus_key$direction))) {
  stop(
    "OWO direction is unavailable for: ",
    paste(genus_key$microbe[is.na(genus_key$direction)], collapse = ", ")
  )
}
if (anyDuplicated(genus_key$common_name)) {
  stop("The current 16S genus set contains duplicated common names.")
}

# Order both sections alphabetically by the displayed common genus name.
genus_key <- genus_key %>%
  mutate(direction = factor(direction, levels = c("Depleted", "Enriched"))) %>%
  arrange(direction, common_name)

depleted_genera <- genus_key %>%
  filter(direction == "Depleted") %>%
  pull(microbe)
enriched_genera <- genus_key %>%
  filter(direction == "Enriched") %>%
  pull(microbe)
selected_genera <- c(depleted_genera, enriched_genera)

# A discrete y axis runs from bottom to top; reverse the desired display order.
genus_levels_bottom_to_top <- rev(genus_key$common_name)

background_data <- tibble(
  clinical_group = factor(clinical_group_levels, levels = clinical_group_levels),
  xmin = c(0.5, 2.5, 6.5, 13.5, 14.5),
  xmax = c(2.5, 6.5, 13.5, 14.5, 15.5)
)

# Prepare bubble-heatmap data 
plot_data <- results_16s %>%
  filter(
    microbiome_type == "16S_g",
    microbiome_period == "T1",
    clinical_trait %in% clinical_traits,
    microbe %in% selected_genera,
    status == "ok"
  ) %>%
  distinct(microbe, clinical_trait, .keep_all = TRUE) %>%
  left_join(
    genus_key %>%
      select(
        microbe, common_name, direction,
        name_16S_analysis, name_16S_original, matched, match_rule
      ),
    by = "microbe"
  ) %>%
  mutate(
    direction = factor(direction, levels = c("Depleted", "Enriched")),
    common_name = factor(
      common_name,
      levels = genus_levels_bottom_to_top
    ),
    clinical_trait = factor(clinical_trait, levels = clinical_traits),
    p_fdr_class = case_when(
      p_fdr < 0.01 ~ "< 0.01",
      p_fdr < 0.05 ~ "< 0.05",
      TRUE ~ "Non-Sig"
    ),
    p_fdr_class = factor(
      p_fdr_class,
      levels = c("Non-Sig", "< 0.05", "< 0.01")
    ),
    beta_plot = if_else(p_value < 0.05, beta, NA_real_)
  ) %>%
  arrange(direction, common_name, clinical_trait)

expected_rows <- length(selected_genera) * length(clinical_traits)
if (nrow(plot_data) != expected_rows) {
  missing_pairs <- tidyr::crossing(
    microbe = selected_genera,
    clinical_trait = clinical_traits
  ) %>%
    anti_join(
      plot_data %>%
        transmute(
          microbe = microbe,
          clinical_trait = as.character(clinical_trait)
        ),
      by = c("microbe", "clinical_trait")
    )

  stop(
    "Expected ", expected_rows, " genus-trait results but found ",
    nrow(plot_data), ". Missing pairs: ",
    paste(
      paste(missing_pairs$microbe, missing_pairs$clinical_trait, sep = "/"),
      collapse = ", "
    )
  )
}

beta_limit <- ceiling(max(abs(plot_data$beta), na.rm = TRUE) * 100) / 100

# Plot 
fig4c <- ggplot(plot_data, aes(x = clinical_trait, y = common_name)) +
  geom_rect( data = background_data,
    aes( xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
      fill = clinical_group
    ), inherit.aes = FALSE, alpha = 0.1, color = NA
  ) +
  scale_fill_manual(
    values = clinical_group_colors,
    breaks = clinical_group_levels,
    name = NULL,
    guide = guide_legend(
      order = 1, nrow = 2, byrow = TRUE,
      override.aes = list(alpha = 0.18, color = NA)
    )
  ) +
  ggnewscale::new_scale_fill() +
  geom_point(
    aes(fill = beta_plot, size = p_fdr_class), shape = 21, color = "black"
  ) +
  facet_grid2(
    direction ~ ., scales = "free_y", space = "free_y",
    strip = strip_themed(
      background_y = elem_list_rect(fill = direction_colors)
    )
  ) +
  scale_x_discrete(
    limits = clinical_traits,
    labels = setNames(clinical_labels, clinical_traits),
    drop = FALSE
  ) +
  # Drop unused genus levels within each free-y facet so depleted and enriched
  # genera appear only in their own panels.
  scale_y_discrete(drop = TRUE) +
  scale_size_manual(
    values = c("Non-Sig" = 1.8, "< 0.01" = 5.2, "< 0.05" = 3.8),
    name = expression(italic(P)[FDR]),
    drop = FALSE,
    guide = guide_legend(
      order = 2, direction = "horizontal", title.position = "left",
      override.aes = list(shape = 21, fill = "white", color = "black")
    )
  ) +
  scale_fill_gradient2(
    low = "#276FB0", mid = "white", high = "#B82530",
    midpoint = 0, limits = c(-beta_limit, beta_limit),
    breaks = c(-0.1, 0, 0.1),
    labels = c("-0.1", "0.0", "0.1"),
    oob = scales::squish, na.value = "white",
    name = "Beta",
    guide = guide_colorbar(
      order = 3, direction = "horizontal", title.position = "left",
      barwidth = unit(3.2, "cm"), barheight = unit(0.38, "cm")
    )
  ) +
  labs(x = NULL, y = NULL, tag = 'c') +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(1.5, "mm"),
    axis.text.y = element_text( size = 12, color = "black", hjust = 1, face = "italic" ),
    axis.text.x = element_text(
      size = 12, color = "black", angle = 40, hjust = 1, vjust = 1
    ),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.y = element_text(size = 15, angle = -90, color = "black"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.justification = "center",
    legend.margin = margin(0, 0, 1, 0),
    legend.spacing.x = unit(2.5, "mm"),
    legend.key.height = unit(4.5, "mm"),
    legend.key.width = unit(5, "mm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.background = element_rect(fill = "white", color = NA),
    plot.tag.position = c(0.05, 0.98),
    plot.tag = element_text(size = 20) 
  )

# Output 
dir.create("figure/PART3", recursive = TRUE, showWarnings = FALSE)
ggsave(
  "figure/PART3/fig4c_thsbc_16s_clinical.pdf",
  plot = fig4c, width = 9.5, height = 6
)
