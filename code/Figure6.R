# fig 6a ------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(forcats)
library(readxl);library(writexl)
library(data.table)
library(haven)
library(ggplot2)
library(patchwork)
library(scales)

pdata_summary <- read_xlsx('pub_code/data/Figure 6.xlsx', sheet = 'figure6a_summary')
pdata <- read_xlsx('pub_code/data/Figure 6.xlsx', sheet = 'figure6a_raw')

pred_stand <- read_xlsx('pub_code/data/Figure 6.xlsx', sheet = 'figure6a_stand')
pred_gene <- read_xlsx('pub_code/data/Figure 6.xlsx', sheet = 'figure6a_gene')
pred_micro <- read_xlsx('pub_code/data/Figure 6.xlsx', sheet = 'figure6a_micro')
pred_metab <- read_xlsx('pub_code/data/Figure 6.xlsx', sheet = 'figure6a_metab')
pred_comb <- read_xlsx('pub_code/data/Figure 6.xlsx', sheet = 'figure6a_comb')


p1 <- ggplot() +
  geom_bar(data = pdata_summary, aes(x = type, y = mean_R2, fill = type), 
           stat = "identity", color = 'black', width = 0.6) +
  #scale_fill_manual(values = c( "#0072b5", "#e18727", "#20854e", "#7876b1" ))+
  scale_fill_manual(values = c( "#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1"))+
  #ggsci::scale_fill_nejm(alpha = 0.8) +
  geom_errorbar(data = pdata_summary, aes(x = type, ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, color = "black") +
  geom_jitter(data = pdata, aes(x = type, y = R2), 
              position = position_jitter(width = 0.1), 
              size = 2, alpha = 0.4, fill = "gray") +
  labs(x = "", y = "Out-of-sample R²")+ 
  scale_y_continuous(expand = c(0, 0), limits=c(0, 0.5),breaks=seq(0, 0.5,0.1))+
  theme_bw() +
  theme(plot.background = element_blank(),  # 移除图表背景
        panel.background = element_blank(),  # 移除绘图区域背景
        legend.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(0.95), colour = 'black'),
        axis.text.x = element_text(angle = 40, hjust = 1), 
        legend.position = 'none')


cols = c("#bc3c29", "#0072b5", "#e18727", "#20854e", "#7876b1")
pfunc_point <- function(data, col, xlab, ylab, title){
  
  cor = cor.test(data$BMI_prep_raw, data$predict_y_raw)$estimate %>% round(., 3)
  pv = cor.test(data$BMI_prep_raw, data$predict_y_raw)$p.value*5
  
  plot <- ggplot(data, aes(x = predict_y_raw, y = BMI_prep_raw)) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    geom_point(color = col, alpha = 0.2, size = 2, shape = 21, fill = "white", stroke = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = col, fill = alpha(col, 0.2)) +
    labs(x = xlab, y = ylab, title = title)+ 
    scale_x_continuous(limits = c(16, 28), breaks=seq(16, 28, 4))+
    scale_y_continuous(limits = c(15, 35), breaks=seq(15, 35, 5))+
    #scale_x_continuous(breaks = seq(20, 22.5, 1))+
    annotate("text", x = 16, y = 34,
             label = paste0("Pearson's r = ", cor, '\n','P = ', scientific(pv, digits = 2)),
             # label = bquote(atop("Pearson's "*italic("r")*" = "~.(cor), 
             #                     italic("P")*" = "*.(scientific(pv, digits = 2)), "")),
             hjust = 0, vjust = 1, size = 3.5) +  # 使用bquote和scientific函数
    theme_classic() +
    theme(plot.background = element_blank(),  # 移除图表背景
          panel.background = element_blank(),  # 移除绘图区域背景
          legend.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = rel(1.2)),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = rel(0.95), colour = 'black'),
          legend.position = 'none')
  
  return(plot)
}

p2 = pfunc_point(pred_stand, "#0072b5", '', expression("Measured pBMI (kg/m2)"), 'Basic')
p3 = pfunc_point(pred_gene, "#bc3c29", '', '', 'Genomics')
p4 = pfunc_point(pred_micro, "#e18727", expression("Inferred pBMI (kg/m2)"), '', 'Microbiome')
p5 = pfunc_point(pred_metab, "#20854e", '', '', 'Metabolomics')
p6 = pfunc_point(pred_comb, "#7876b1", '', '', 'Multi-omics')

all <- (p1+ p2 + p3 + p4 + p5 + p6) + plot_layout(ncol = 6) &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))
all



# fig 6b ------------------------------------------------------------------

rm(list = ls())

library(ggradar)

pdata1 <- read_xlsx('pub_code/data/Figure 6.xlsx', sheet = 'figure6b_raw')
pdata2 <- read_xlsx('pub_code/data/Figure 6.xlsx', sheet = 'figure6b_omic')

ggradar(pdata1, #font.radar = "Calibri",
        background.circle.colour = "transparent",   # 背景颜色         
        values.radar = c("0", "10%", "20%"),       # 不显示刻度标签 
        grid.min = 0, grid.mid = 10, grid.max = 20,
        group.line.width = 0.8,               # 线宽          
        group.point.size = 2,               # 数据点大小          
        #group.colours = c("grey", "#8bcee5", "#25ab60"),
        group.colours = c("grey", "#f5a3a3", "#a4d7f3"),       
        fill = TRUE, fill.alpha = 0.2,
        
        plot.title = "Phenotypic vs. Genetic Mismatch of pBMI",
        
        grid.label.size = 4,
        gridline.label.offset = 0.1,
        
        axis.label.size = 4,                # 轴标签字体大小          
        axis.label.offset = 1.1,          
        axis.line.colour = "grey",            # 轴线颜色          
        grid.line.width = 0.4,  
        
        legend.text.size = 10,
        legend.position = c(0.95, 0.2),
        
        gridline.min.linetype = "solid",      # 最小网格线          
        gridline.max.linetype = "solid",      # 最大网格线          
        gridline.max.colour = "black"
) +  
  theme(plot.background = element_blank(),  
        legend.background = element_blank(),  
        legend.box.background = element_blank(),  
        panel.background = element_blank(),   
        plot.margin = margin(0, 0, 0, 0))


ggradar(pdata2, #font.radar = "Calibri",
        background.circle.colour = "transparent",   # 背景颜色         
        values.radar = c("0", "15%", "30%"),       # 不显示刻度标签 
        grid.min = 0, grid.mid = 15, grid.max = 30,
        group.line.width = 0.8,               # 线宽          
        group.point.size = 2,               # 数据点大小          
        group.colours = c("grey", "#fcbb71", "#8faadc"),     
        fill = TRUE, fill.alpha = 0.2,
        
        plot.title = "Phenotypic vs. Genetic Mismatch of Predicted pBMI",
        
        grid.label.size = 4,
        gridline.label.offset = 0.1,
        
        axis.label.size = 4,                # 轴标签字体大小          
        axis.label.offset = 1.1,          
        axis.line.colour = "grey",            # 轴线颜色          
        grid.line.width = 0.4,  
        
        legend.text.size = 10,
        legend.position = c(0.98, 0.2),
        
        gridline.min.linetype = "solid",      # 最小网格线          
        gridline.max.linetype = "solid",      # 最大网格线          
        gridline.max.colour = "black"
) +  
  theme(plot.background = element_blank(),  
        legend.background = element_blank(),  
        panel.background = element_blank(),   
        plot.margin = margin(0, 0, 0, 0))



# fig 6d ------------------------------------------------------------------

rm(list = ls())

pdata <- read_xlsx('pub_code/data/Figure 6.xlsx', sheet = 'figure6d')

pright <- pdata %>% filter(popclass %in% c('N_N', 'N_O')) %>% 
  mutate(popclass = factor(popclass, c('N_N', 'N_O')))
rightplot <- ggplot(pright, aes(Frequency, outcome,group = popclass, fill = popclass)) + 
  geom_bar(stat = "identity", position = "dodge", color = "white", width =0.8)+
  scale_x_continuous(expand = c(0,0), limits = c(0,75),breaks = seq(0,75,15))+ 
  scale_y_discrete(expand = c(0.05,0.05))+
  #scale_fill_manual(values = c('#aebfce', '#ceb5b9')) + 
  scale_fill_manual(values = c( '#a1c6e7','#b9d9eb'),
                    name = "",  
                    labels = c("Observed non-OWO & Predicted non-OWO", "Observed non-OWO & Predicted OWO")) + 
  labs(subtitle = "Non-OWO", x = "", y = "") + 
  theme_classic() + 
  theme(plot.background = element_blank(),  
        panel.background = element_blank(),  
        #text = element_text(family = "Calibri"),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        axis.text.y= element_text(size = 10,hjust = 0.5, colour = 'black'),
        axis.text.x = element_text(size = 10, colour = 'black'),
        axis.title = element_text(size = rel(1.2)),
        plot.subtitle = element_text(hjust = 0.5, vjust = -1, size = 12),
        plot.margin = margin(r = 15, b=5),
        legend.title = element_text(size = 10),
        legend.key.size = unit(13, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.position = c(0.7, 0.3)
  )

pleft <- pdata %>% filter(popclass %in% c('O_N', 'O_O')) %>% 
  mutate(popclass = factor(popclass, c('O_O', 'O_N')))
leftplot <- ggplot(pleft, aes(-Frequency, outcome,group = popclass, fill = popclass)) + 
  geom_bar(stat = "identity", position = "dodge", color = "white", width =0.8)+
  scale_y_discrete(expand = c(0.05,0.05), position = "right") +
  scale_x_continuous(expand = c(0,0), limits = c(-75, 0),breaks = seq(-75,0,15),
                     labels = c('75', '60', '45', '30', '15', '0'))+ 
  scale_fill_manual(values = c('#f4a261', '#f8c3b1'),
                    name = "",  
                    # labels = c("Matched", "Mismatched")
                    labels = c("Observed OWO & Predicted OWO", "Observed OWO & Predicted non-OWO")
  ) + 
  labs(subtitle = "OWO", x = "", y = "") + 
  theme_classic() + 
  theme(plot.background = element_blank(),  # 移除图表背景
        panel.background = element_blank(),  # 移除绘图区域背景
        legend.background = element_blank(),
        #text = element_text(family = "Calibri"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, colour = 'black'),
        axis.title = element_text(size = rel(1.2)),
        plot.subtitle = element_text(hjust = 0.5, vjust = -1, size = 12),
        legend.position = c(0.4, 0.3),
        legend.title = element_text(size = 10),
        legend.key.size = unit(13, "pt"),
        plot.margin = margin(l = 15, b=5),
        legend.margin = margin(0, 0, 0, 0)
  )

library(cowplot)
combined_plot <- ggdraw() +
  draw_plot(leftplot, x = 0, width = 0.46) +  
  draw_plot(rightplot, x = 0.42, width = 0.58)  +
  draw_label("Frequency of health outcome (%)", 
             x = 0.5, y = 0.03, hjust = 0.5, vjust = 0.5, size = 13)  
combined_plot



