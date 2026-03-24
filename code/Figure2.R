


# fig 2a ------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(data.table)
library(readxl);library(writexl)
library(haven)

library(vegan)
library(compositions)

library(ggplot2)
library(ggpubr)
library(ggExtra)
library(patchwork)
library(ComplexHeatmap)

pdata_heatmap <- read_xlsx('pub_code/data/Figure 2.xlsx', sheet = 'figure2a_heat')
pdata_bar <- read_xlsx('pub_code/data/Figure 2.xlsx', sheet = 'figure2a_bar')

var_gp <- columnAnnotation(fa = c('Anthropometric','Anthropometric','Anthropometric',
                                  'Demography','Diet',
                                  'Diet', 'Demography', 'Anthropometric',
                                  'Diet','Demography',
                                  'Diet','Diet',
                                  'Lifestyle','Lifestyle',
                                  'Diet','Lifestyle'),
                           col = list(fa = c("Anthropometric" = "#82b181", 
                                             "Diet" = "#87bcbd", 
                                             "Lifestyle" = "#bbaacd", 
                                             "Demography" = "#dab736")),
                           gp = gpar(col = "white"), 
                           show_legend = FALSE)

heat1 <- Heatmap(as.matrix(pdata_heatmap), 
                 col = structure(c("white", "#CD5C5C"), names = c("non-sig", "sig")),
                 show_row_dend = F,
                 show_row_names = T,
                 row_names_gp = gpar(fontsize = 14),
                 row_names_max_width = unit(25, "cm"),
                 row_labels = rev(c("Metabolites", "Bact", "Fung", "MetaG_g", "MetaG_s")),
                 column_names_side = "bottom", row_names_side = 'left',
                 column_names_gp = gpar(fontsize = 14, hjust = -0.1),
                 column_names_rot = 45,
                 column_title = NULL,
                 rect_gp = gpar(col = "white", lwd = 3),
                 cluster_columns = FALSE, cluster_rows = FALSE, 
                 show_heatmap_legend = F,
                 heatmap_legend_param = list(title = 'Statistics'),
                 bottom_annotation = var_gp)

figure2a_bar <- ggplot()+  #柱状堆积图  
  geom_col(data = pdata_bar,aes(variable, R2, fill=type)) + 
  labs(x=NULL,y=NULL,fill=NULL) +  
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.2,0.4,0.6,0.8,1.0)) +  
  scale_fill_manual(values = rev(c("#7EABCA","#FFB2B9", "#95C6C6","#ebdcb2", "#d4cfdf")))+ 
  theme_classic()+  
  theme(axis.text.y = element_text(color = "black",size=12),        
        axis.text.x = element_blank(),
        legend.position = 'top',
        legend.text = element_text(color = "black",size=12), 
        plot.margin = margin(t = 2, r = 12, b = 6, l = 55)) 

pdf(file = 'figure2a.pdf', width = 8, height = 7.5)
library(gridExtra)
ht_grob = grid.grabExpr(draw(heat1, padding = unit(c(0.01, 0.01, 0.01, 0.01), 'cm')))
grid.arrange(figure2a_bar, ht_grob, ncol=1, heights=c(0.60, 0.40)) 
dev.off()


# fig 2b ------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(data.table)
library(forcats)
library(readxl);library(writexl)

library(vegan) # For dist
library(pairwiseAdonis)  # For pairwise Adonis
library(patchwork) # For combine plots

library(ggExtra) # For marginal plot
library(ggpubr)

df_pca <- read_xlsx('pub_code/data/Figure 2.xlsx', sheet = 'figure2b')

# 绘图
col_group = c('Underweight'="#c0b9e6", 'Normal'="#b0d2e8",'OWO'="#f7d4bc")
col_period = shades::saturation(c('T1'='#be8a7f','T2'='#a2b98f','T3'='#e2c083'), 0.6)
# 主图
ggplot(df_pca, aes(x = PC1, y = PC2,color = BMI_prep_cut3)) +
  geom_point(aes(shape = period), size = 1.5, alpha = 0.8)+
  stat_ellipse(aes(fill = BMI_prep_cut3), geom = "polygon",
               level = 0.95, linetype = 2, linewidth= 0.8, alpha = 0.2) +
  scale_x_continuous(limits = c(-25, 25))+
  scale_y_continuous(limits = c(-15, 25))+
  scale_color_manual(values = col_group, guide = "none") +
  scale_fill_manual(values = col_group, guide = "none") +
  labs(x=sprintf('PC1 (%.1f%%)', PCA_df2$contr[1]),
       y=sprintf('PC2 (%.1f%%)', PCA_df2$contr[2])) +
  theme_test(base_size = 15) + 
  theme(
    panel.background = element_blank(), 
    plot.background = element_blank(), 
    panel.grid = element_blank(),
    legend.position = c(0.92, 0.1),
    legend.title = element_blank(),
    legend.background = element_blank(), 
    legend.key.size = unit(14, "pt"), 
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12, color = 'black'),
    plot.margin = margin(3,3,3,3)
  ) 

ggplot(data=df_pca,
       aes(x= BMI_prep_cut3, y=PC1, fill=BMI_prep_cut3, color=BMI_prep_cut3)) +
  facet_wrap(. ~ period, nrow = 1, ncol = 3) +
  geom_boxplot(size= 0.8, width = 0.3, alpha = 0.2, outlier.shape = NA, 
               position = position_dodge(0.4)) +
  stat_compare_means(aes(group = BMI_prep_cut3), method = "wilcox.test", 
                     label.x.npc = 0.5,label.y.npc = 0.9,
                     ref.group = 'Normal',
                     label = "p.signif", show.legend = F, 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                     size = 5) +
  labs(y = "PCA1", x = "") +
  scale_color_manual(values = col_group, guide = "none") +
  scale_fill_manual(values = col_group, guide = "none") + 
  scale_y_continuous(limits = c(-20, 20))+
  theme_bw() +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8),
        panel.background = element_blank(),
        axis.title = element_text(color = "black",size = 14),
        axis.text = element_text(color = "black",size = 12),
        axis.text.x = element_text(color = "black",angle = 45,size = 12, hjust = 1),
        strip.text = element_text(color = "black",size = 14),
        legend.title = element_blank(),
        legend.key.size = unit(12, "pt"), 
        legend.text = element_text(size = 10),
        legend.position = 'none')



# fig 2c ------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(data.table)
library(forcats)
library(readxl);library(writexl)
library(ggplot2)
library(ggridges)
library(ggpubr)

pdata <- read_xlsx('pub_code/data/Figure 2.xlsx', sheet = 'figure2c')

p1 = ggplot(pdata, aes(x = richness, y = BMI_prep_cut3)) +  
  geom_density_ridges(aes(fill = BMI_prep_cut3),
                      color = 'white', alpha = 0.9,
                      rel_min_height = 0.02,
                      quantile_lines = TRUE, quantiles = 2) +  
  scale_fill_manual(values = c("#c0b9e6", "#b0d2e8","#f7d4bc")) +
  facet_grid(micro_type ~ ., switch = "y") +
  labs(y = "", x = "Richness (T1)")+
  scale_x_continuous(expand = c(0.02,0), limits = c(0, 240))+
  theme_classic2() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y =  element_blank(),
        axis.text.y =  element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 14),
        strip.placement = "outside",  # 确保标签在绘图区域外
        strip.text.y.left = element_text(
          angle = 90, vjust = -0.8, hjust = 0.5,
          size = 13, color = "black"),
        strip.background = element_rect(fill = "transparent", color = NA),
        legend.position = 'none',
        plot.margin = margin(1, 1, 1, 1, unit = "pt"))

p2 = ggplot(pdata, aes(x = shannon, y = BMI_prep_cut3)) +  
  geom_density_ridges(aes(fill = BMI_prep_cut3),
                      color = 'white', alpha = 0.9,
                      rel_min_height = 0.02,
                      quantile_lines = TRUE, quantiles = 2) +  
  scale_fill_manual(values = c("#c0b9e6", "#b0d2e8","#f7d4bc")) +
  facet_grid(micro_type ~ ., switch = "y") +
  labs(y = "", x = "Shannon (T1)")+
  scale_x_continuous(expand = c(0.02,0), limits = c(-0.5, 4.5))+
  theme_classic2() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks.y =  element_blank(),
        axis.text.y =  element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 14),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_blank(),
        legend.key.size = unit(14, "pt"), 
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.margin = margin(1, 1, 1, 1, unit = "pt"))

ggarrange(p1, p2, ncol = 2)


# fig 2d ------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(data.table)
library(forcats)
library(readxl);library(writexl)
library(ggplot2)
library(ggridges)
library(ggpubr)

bact_fb <- read_xlsx('pub_code/data/Figure 2.xlsx', sheet = 'figure2d_bact')
mgs_fb <- read_xlsx('pub_code/data/Figure 2.xlsx', sheet = 'figure2d_mgs')

p3 <- ggplot(bact_fb, aes(x = fb_ratio, y = BMI_prep_cut3)) +  
  geom_density_ridges(aes(fill = BMI_prep_cut3),
                      color = 'white', alpha = 0.9,
                      rel_min_height = 0.08,
                      quantile_lines = TRUE, quantiles = 2) +  
  scale_fill_manual(values = c("#c0b9e6", "#b0d2e8","#f7d4bc")) +
  facet_grid(period ~ ., switch = "y") +
  labs(y = "", x = "log(F/B ratio) (16S rRNA)")+
  annotate("text", x = 8, y = 4.2, label = "Ptrend=0.001") +
  annotate("text", x = 8, y = 1, label = "*", size = 5) +
  annotate("text", x = 8, y = 3, label = "*", size = 5) +
  scale_x_continuous(expand = c(0.02,0), limits = c(0, 10))+
  theme_classic2() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks.y =  element_blank(),
        axis.text.y =  element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 14),
        strip.placement = "outside",  # 确保标签在绘图区域外
        strip.text.y.left = element_text(
          angle = 90, vjust = -0.8, hjust = 0.5,
          size = 13, color = "black"),
        strip.background = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(14, "pt"), 
        legend.text = element_text(size = 12),
        legend.position = 'none')

p4 <- ggplot(mgs_fb, aes(x = fb_ratio, y = BMI_prep_cut3)) +  
  geom_density_ridges(aes(fill = BMI_prep_cut3),
                      color = 'white', alpha = 0.9,
                      rel_min_height = 0.08,
                      quantile_lines = TRUE, quantiles = 2) +  
  scale_fill_manual(values = c("#c0b9e6", "#b0d2e8","#f7d4bc")) +
  facet_grid(period ~ ., switch = "y") +
  labs(y = "", x = "log(F/B ratio) (Metagenome)") +
  annotate("text", x = 12, y = 4.2, label = "Ptrend=0.001") +
  annotate("text", x = 11, y = 1, label = "*", size = 5) +
  annotate("text", x = 11, y = 3, label = "*", size = 5) +
  scale_x_continuous(expand = c(0.02,0), limits = c(-0.5, 14))+
  theme_classic2() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks.y =  element_blank(),
        axis.text.y =  element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 14),
        strip.placement = "outside",  # 确保标签在绘图区域外
        strip.text.y.left = element_blank(),
        strip.background = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(14, "pt"), 
        legend.text = element_text(size = 12),
        legend.position = 'none')

ggarrange(p3, p4, ncol = 2)

# fig 2e ------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(data.table)
library(readxl);library(writexl)
library(haven)

library(vegan)
library(compositions)

library(ggplot2)
library(ggpubr)


ggplot(data=pdata[which(!is.na(pdata$time_group2)),],
       aes(x= time_group2, y=distance, fill=time_group2, color=time_group2)) +
  facet_wrap(. ~ class, nrow = 1, ncol = 5) +
  geom_violin(width = 0.5, alpha = 0.2, position = position_dodge(0.4)) +
  geom_boxplot(size= 0.8, width = 0.3, alpha = 0.1, outlier.shape = NA, 
               position = position_dodge(0.4)) +
  stat_compare_means(aes(group = time_group2), method = "wilcox.test", label.x.npc = 0.5,
                     label = "p.signif", show.legend = F, 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                     size = 6, label.y.npc = 0.9) +
  labs(y = "Dissimilarity based on distance", x = "") +
  scale_color_manual(values = c('#E6956F', '#788FCE')) +  
  scale_fill_manual(values = c('#E6956F', '#788FCE')) +   
  theme_bw() +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8),
        panel.background = element_blank(),
        axis.title = element_text(color = "black",size = 14),
        axis.text = element_text(color = "black",size = 12),
        axis.text.x = element_text(color = "black",size = 12, angle = 45, hjust=1),
        strip.text = element_text(color = "black",size = 14),
        legend.title = element_blank(),
        legend.key.size = unit(12, "pt"), 
        legend.text = element_text(size = 10),
        legend.position = 'none')



