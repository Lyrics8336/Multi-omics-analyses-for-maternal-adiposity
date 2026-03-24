

# fig 4a ------------------------------------------------------------------
# graphlan
# cat global_anno.txt track0.txt > track0_all
# graphlan_annotate.py --annot track0_all tree1_backbone.txt graphlan.xml
# graphlan_annotate.py --annot ring1.txt graphlan.xml graphlan1.xml
# graphlan_annotate.py --annot ring2.txt graphlan1.xml graphlan2.xml
# graphlan_annotate.py --annot ring3.txt graphlan2.xml graphlan3.xml
# graphlan_annotate.py --annot ring4.txt graphlan3.xml graphlan4.xml
# graphlan_annotate.py --annot label5.txt graphlan4.xml graphlan5.xml
# 
# graphlan_annotate.py --annot ring6.txt graphlan5.xml graphlan6.xml
# graphlan_annotate.py --annot ring7.txt graphlan6.xml graphlan7.xml
# graphlan_annotate.py --annot ring8.txt graphlan7.xml graphlan8.xml
# graphlan_annotate.py --annot ring9.txt graphlan8.xml graphlan9.xml
# graphlan_annotate.py --annot bar10.txt graphlan9.xml graphlan10.xml
# graphlan.py graphlan10.xml graphlan_tree_v5.pdf --dpi 300 --size 4 --pad 0.0




# fig 4b ------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(forcats)
library(readxl);library(writexl)
library(haven)

library(ggplot2)
library(ggbreak)
library(cowplot)

## genus
htotu <- read_xlsx('pub_code/data/Figure 4.xlsx', sheet = 'figure4b_genus_otu')
htpv <- read_xlsx('pub_code/data/Figure 4.xlsx', sheet = 'figure4b_genus_pv')
htbeta <- read_xlsx('pub_code/data/Figure 4.xlsx', sheet = 'figure4b_genus_beta')

col_fun = colorRamp2(c(-0.5, 0, 0.7), 
                     c("#5A5A83", "white", "#B1615C"))
met_gp = rowAnnotation(dir = htotu$dir,
                       phylum = htotu$OTU,
                       show_annotation_name = F,
                       col = list(dir = c("Enriched in OWO" = "#f4c2c2",
                                          "Depleted in OWO" = '#a3c4e0'),
                                  phylum = c("Actinobacteriota" = "#6B6CA3",
                                             "Bacteroidota" = "#D9A34D",
                                             "Firmicutes" = "#80A46E",
                                             "Proteobacteria" = '#87bcbd',
                                             "Others" = '#D3D3D3'),
                                  width = unit(2, "mm")),
                       #gp = gpar(color = 'white'),
                       annotation_legend_param = list(
                         dir = list(
                           title = "",
                           title_gp = gpar(col = "black", fontsize = 12),
                           labels_gp = gpar(col = "black", fontsize = 10),
                           title_position = "topleft"
                         ),
                         phylum = list(
                           title = "Phylum",
                           title_gp = gpar(col = "black", fontsize = 12),
                           labels_gp = gpar(col = "black", fontsize = 10),
                           title_position = "topleft"
                         )
                       ))

main1 = Heatmap(as.matrix(htbeta), 
                col = col_fun,
                row_dend_width = unit(2.5, "cm"),
                row_dend_gp = gpar(lwd = 2),
                show_row_dend = F,
                show_row_names = T,
                row_names_gp = gpar(fontsize = 11, col = ifelse(rownames(htotu) %in% rep_genus, 'black','grey')),
                row_names_max_width = unit(8, "cm"),
                column_split = factor(rep(c("THSBC", "HBC"), each = 4), levels = c("THSBC", "HBC")), 
                column_labels = c("T1", "T2", "T3", "Pooled", "T1", "T2", "T3", "Pooled"),
                #column_labels = c("T1", "T2", "T3", "T1", "T2", "T3", "Pooled", "HBC_Pooled"),
                # column_labels = c("THSBC_T1", "HBC_T1", "THSBC_T2", "HBC_T2", "THSBC_T3", "HBC_T3", "THSBC_Pooled", "HBC_Pooled"),
                column_names_side = "bottom",
                column_names_gp = gpar(fontsize = 10),
                column_names_rot = 30,
                column_title = NULL,
                border = TRUE,
                rect_gp = gpar(col = "white", lwd = 1.5),
                border_gp = gpar(col = "black", lwd = 2),
                cluster_columns = FALSE,  cluster_rows = FALSE,
                cell_fun = function(j, i, x, y, w, h, fill){
                  if(htpv[i, j] < 0.01) {
                    grid.text("*", x, y, vjust = 0.8,
                              gp = gpar(fontsize = 27))
                  } else if(htpv[i, j] < 0.05){
                    grid.text("+", x, y, #vjust = 0.8,
                              gp = gpar(fontsize = 14))
                  }},
                show_heatmap_legend = T,
                heatmap_legend_param = list(title = 'Beta'),
                left_annotation = met_gp)
draw(main1)

## species
htotu <- read_xlsx('pub_code/data/Figure 4.xlsx', sheet = 'figure4b_spec_otu')
htpv <- read_xlsx('pub_code/data/Figure 4.xlsx', sheet = 'figure4b_spec_pv')
htbeta <- read_xlsx('pub_code/data/Figure 4.xlsx', sheet = 'figure4b_spec_beta')


col_fun = colorRamp2(c(-0.5, 0, 0.7), 
                     c("#5A5A83", "white", "#B1615C"))
met_gp = rowAnnotation(dir = htotu$dir,
                       phylum = htotu$OTU,
                       show_annotation_name = F,
                       col = list(dir = c("Enriched in OWO" = "#f4c2c2",
                                          "Depleted in OWO" = '#a3c4e0'),
                                  phylum = c("Actinobacteriota" = "#6B6CA3",
                                             "Bacteroidota" = "#D9A34D",
                                             "Firmicutes" = "#80A46E",
                                             "Proteobacteria" = '#87bcbd',
                                             "Others" = '#D3D3D3'),
                                  width = unit(2, "mm")),
                       #gp = gpar(color = 'white'),
                       annotation_legend_param = list(
                         dir = list(
                           title = "",
                           title_gp = gpar(col = "black", fontsize = 12),
                           labels_gp = gpar(col = "black", fontsize = 10),
                           title_position = "topleft"
                         ),
                         phylum = list(
                           title = "Phylum",
                           title_gp = gpar(col = "black", fontsize = 12),
                           labels_gp = gpar(col = "black", fontsize = 10),
                           title_position = "topleft"
                         )
                       ))

main2 = Heatmap(as.matrix(htbeta), 
                col = col_fun,
                row_dend_width = unit(2.5, "cm"),
                row_dend_gp = gpar(lwd = 2),
                show_row_dend = F,
                show_row_names = T,
                row_names_gp = gpar(fontsize = 11, col = ifelse(rownames(htotu) %in% rep_spec, 'black','grey')),
                row_names_max_width = unit(8, "cm"),
                column_split = factor(rep(c("THSBC", "HBC"), each = 4), levels = c("THSBC", "HBC")), 
                column_labels = c("T1", "T2", "T3", "Pooled", "T1", "T2", "T3", "Pooled"),
                #column_labels = c("T1", "T2", "T3", "T1", "T2", "T3", "Pooled", "HBC_Pooled"),
                # column_labels = c("THSBC_T1", "HBC_T1", "THSBC_T2", "HBC_T2", "THSBC_T3", "HBC_T3", "THSBC_Pooled", "HBC_Pooled"),
                column_names_side = "bottom",
                column_names_gp = gpar(fontsize = 10),
                column_names_rot = 30,
                column_title = NULL,
                border = TRUE,
                rect_gp = gpar(col = "white", lwd = 1.5),
                border_gp = gpar(col = "black", lwd = 2),
                cluster_columns = FALSE,  cluster_rows = FALSE,
                cell_fun = function(j, i, x, y, w, h, fill){
                  if(htpv[i, j] < 0.01) {
                    grid.text("*", x, y, vjust = 0.8,
                              gp = gpar(fontsize = 27))
                  } else if(htpv[i, j] < 0.05){
                    grid.text("+", x, y, #vjust = 0.8,
                              gp = gpar(fontsize = 14))
                  }},
                show_heatmap_legend = T,
                heatmap_legend_param = list(title = 'Beta'),
                left_annotation = met_gp)
draw(main2)


ht_list <- main1 %v% main2
draw(ht_list, ht_gap = unit(0.5, "cm"),
     column_title = "THSBC    HBC",  
     column_title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
     heatmap_legend_side = "bottom",merge_legend = TRUE,
     annotation_legend_side = "bottom") 


# fig 4c ------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(forcats)
library(readxl);library(writexl)
library(haven)
library(furrr)
library(survey)
library(ggrepel) # For labels
library(pheatmap)
library(ggforce)
library(ggh4x)
library(MetBrewer)
library(ggtext)

tmp <- read_xlsx('pub_code/data/Figure 4.xlsx', sheet = 'figure4c')

ggplot(tmp, aes(x=outcome, y=name))+
  facet_grid2(dir~., scales = 'free', space = 'free',
              strip = strip_themed(background_y = elem_list_rect(fill = col_dir))
  ) +
  annotate("rect",xmin=-Inf,xmax=2.5,ymin=-Inf,ymax=Inf,alpha=0.1,fill="#F8766D")+
  annotate("rect",xmin=2.5,xmax=6.5,ymin=-Inf,ymax=Inf,alpha=0.1,fill="#A3A500")+
  annotate("rect",xmin=6.5,xmax=13.5,ymin=-Inf,ymax=Inf,alpha=0.1,fill="#00BF7D")+
  annotate("rect",xmin=13.5,xmax=15.5,ymin=-Inf,ymax=Inf,alpha=0.1,fill="#00B0F6")+
  annotate("rect",xmin=15.5,xmax=16.5,ymin=-Inf,ymax=Inf,alpha=0.1,fill="#E76BF3")+
  geom_point(aes(fill = beta_0, size = pv_class), shape=21)  +
  labs(fill = "Beta") +
  scale_fill_gradient2(low="#276fb0", mid='white', high = '#b82530',
                       limits = c(-0.13, 0.19),
                       na.value = 'transparent') +
  theme_minimal() +  
  theme(axis.text.y = element_text(size = 12, color = 'black', hjust = 1),
        panel.grid.major = element_line(colour = "grey90",size=0.2),        
        panel.grid.minor = element_blank(),     
        axis.text.x = element_text(color = 'black', size = 12, angle = 40, 
                                   vjust = 1, hjust = 1),
        axis.title = element_blank(),
        strip.background = element_rect(fill = 'grey90', color = 'transparent'),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(size = 12, angle = -90),
        legend.position = 'top',
        plot.background = element_rect(fill = 'white', color = 'transparent'))










