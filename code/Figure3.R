# fig 3a ------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(data.table)
library(forcats)
library(readxl);library(writexl)
library(haven)
library(furrr)
library(survey)

library(circlize) # For circle plot
library(ComplexHeatmap) # For heatmap
library(RColorBrewer) # For color
library(ggrepel) # For labels

library(pheatmap)


pdata <- read_xlsx('pub_code/data/Figure 3.xlsx', sheet = 'figure3a')


## heatmap
beta_mat <- pdata %>% select(starts_with("beta")) %>% as.matrix()
beta_mat <- beta_mat[, c(4,3,2,1)]
rownames(beta_mat) <- pdata$Compounds
beta_col = colorRamp2(c(min(beta_mat), 0, max(beta_mat)), c("#003366", "white", "#990033"))


circos.clear()

circos.par(
  "track.height" = 0.2,
  start.degree = 90,
  clock.wise = T,
  gap.after = 90,
  track.margin = c(0.005, 0.005),
  cell.padding = c(0, 0.8, 0, 0.8),
  points.overflow.warning = T,
  canvas.xlim = c(-0.8, 0.8), 
  canvas.ylim = c(-0.8, 0.8)
)

# circos.initialize(pdata$sectors, x = pdata$x)
circos.heatmap.initialize(beta_mat, split = pdata$y, cluster = F)


### track1: heatmap
circos.heatmap(beta_mat,
               col = beta_col,
               bg.border = "white",
               cell.border = "white",
               cluster=F,
               rownames.side = "outside", rownames.cex = 0.55,
               track.height=0.2)

seq = which(pdata$apv_T1 < 0.05)
circos.points(seq-0.5, 1-0.5, pch = "*", col = "black", cex = 0.6)

seq = which(pdata$apv_T2 < 0.05)
circos.points(seq-0.5, 2-0.5, pch = "*", col = "black", cex = 0.6)

seq = which(pdata$apv_T3 < 0.05)
circos.points(seq-0.5, 3-0.5, pch = "*", col = "black", cex = 0.6)

seq = which(pdata$apv_TA < 0.05)
circos.points(seq-0.5, 4-0.5, pch = "*", col = "black", cex = 0.6)

### track2: bar-point
plot_circos_track <- function(temp_value, color, ylim_mult = c(0.8, 1.1), track_height = 0.12) {
  circos.track(
    ylim = range(temp_value, na.rm = TRUE) * ylim_mult,
    bg.border = "black",
    # bg.col = NA,
    track.height = track_height,
    panel.fun = function(x, y) {
      y = temp_value[CELL_META$subset]
      circos.yaxis(
        side = "left",
        at = c(0 * min(temp_value),
               1,
               round(max(temp_value, na.rm = TRUE), 2)),
        sector.index = get.all.sector.index()[1],
        labels.cex = 0.4, #tick = F,
        labels.niceFacing = F
      )
      
      circos.segments(0, 1, 120, 1, col = 'black', lty = "dotted")
      
      circos.lines(
        x = seq_along(y) - 0.5,
        y = temp_value[CELL_META$subset],
        type = "h",
        col = color,
        lwd = 1.3
      )
      
      circos.points(
        x = seq_along(y) - 0.5,
        y =  temp_value[CELL_META$subset],
        pch = 16,
        cex = 0.6,
        col = color
      )
      
      # circos.text(x = seq_along(y) - 0.5, 
      #             y = 2 + mm_y(4), 
      #             labels = pdata$Compounds,
      #             facing = "clockwise",
      #             niceFacing = TRUE,
      #             cex = 0.5
      # )
      
      
    }
  )
  
}

plot_circos_track(pdata$vip, '#B98FAA')

### track3: group info
class_col = c("AAM" = "#94CEC5", "FA"= "#F6F6BA",
              "OAD" = "#C1BFDB", "BSD"= "#EC8479",
              "GP" = "#86B1CF", "NM" = "#F4B56C",
              "BA" = "#EEDD82", "HC" = "#FCCDE5",
              "HHC" = "#BC80BD", "Others" = "#D9D9D9")

class_col = metab_class_col[!names(metab_class_col) %in% c('CV', 'CHM')]

# library(RColorBrewer)
# col_group_type <- structure(brewer.pal(length(unique(pdata$class1)), "Set3"), names = unique(pdata$class1))
circos.heatmap(pdata$class1, col = class_col, track.height = 0.01)

text(0, 0, "120 serum\nmetabolites", cex = 1)

circos.clear()


dev.off()


## bar plot
num1 = pdata %>% filter(apv_T1 < 0.05) %>% group_by(class1) %>% summarise(Up = sum(beta_T1 >0), Down = sum(beta_T1 <0)) %>% mutate(time = 'T1')
num2 = pdata %>% filter(apv_T2 < 0.05) %>% group_by(class1) %>% summarise(Up = sum(beta_T2 >0), Down = sum(beta_T2 <0)) %>% mutate(time = 'T2')
num3 = pdata %>% filter(apv_T3 < 0.05) %>% group_by(class1) %>% summarise(Up = sum(beta_T3 >0), Down = sum(beta_T3 <0)) %>% mutate(time = 'T3')
num4 = pdata %>% filter(apv_TA < 0.05) %>% group_by(class1) %>% summarise(Up = sum(beta_TA >0), Down = sum(beta_TA <0)) %>% mutate(time = 'Pooled')

pdata_bar <- Reduce(rbind, list(num1, num2, num3, num4)) %>% 
  mutate(time = factor(time, levels = c('Pooled', 'T3', 'T2', 'T1'))) %>% 
  mutate(class1 = factor(class1, levels = c('Others', 'BA', 'AA', 'HHC', 'HetC', 'BSD',
                                            'GP', 'NTM', 'FA', 'OAD', 'AAM'))) 

Up = ggplot(pdata_bar, aes(x=time, y=Up, fill=class1))+
  geom_bar(stat = 'identity',width = 0.6) +
  scale_fill_manual(values  = class_col) +
  scale_y_continuous(expand = c(0, 0.02), limits = c(0, 80))+
  scale_x_discrete(expand = c(0.08, 0.08))+
  labs(x = "", y = "") +
  theme_bw() +guides(fill = "none")+
  theme(plot.background = element_blank(),  # 移除图表背景
        panel.background = element_blank(),  # 移除绘图区域背景
        legend.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 13, colour = 'black'),
        axis.text.x = element_text(size = 13, colour = 'black', face = "bold")
  )

Down = ggplot(pdata_bar, aes(x=time, y=Down, fill=class1))+
  geom_bar(stat = 'identity',width = 0.8) +
  scale_fill_manual(values  = class_col) +
  scale_y_continuous(expand = c(0, 0.02), limits = c(0, 20))+
  scale_x_discrete(expand = c(0.08, 0.08))+
  labs(x = "", y = "") +
  theme_bw() +guides(fill = "none")+
  theme(plot.background = element_blank(),  # 移除图表背景
        panel.background = element_blank(),  # 移除绘图区域背景
        legend.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.text.x = element_blank())


library(cowplot)
ggdraw() +
  draw_plot(Up, y = 0, height = 0.5) +  
  draw_plot(Down, y = 0.48, height = 0.5)  +
  draw_label("Number of significant metabolites", 
             x = 0.1, y = 0.5, angle = 90, vjust = -0.8, size = 16) +
  draw_label("Up", x = 0.55, y = 0.46, size = 12) +
  draw_label("Down", x = 0.55, y = 0.95, size = 12)

# fig 3b ------------------------------------------------------------------

rm(list = ls())
library(ggsankeyfier)
library(patchwork)
library(RColorBrewer)


library(ggsankey)
library(ggplot2)
library(cols4all)
library(cowplot)

dot <- read_xlsx('pub_code/data/Figure 3.xlsx', sheet = 'figure3b_dot')
df <- read_xlsx('pub_code/data/Figure 3.xlsx', sheet = 'figure3b_df')

p1 <- ggplot() + 
  geom_point(data= dot,      
             aes(x = p_trans, y = num_y,        
                 size = hits, color = hit_ratio)) +# 气泡大小及颜色设置 
  theme_bw() + 
  labs(x ="-log(P value)", y="", size = 'Count',color = 'Hit ratio')+
  scale_x_continuous(expand = c(0.1,0.1),limits = c(1.2,7.5), breaks = seq(2,7,1))+
  scale_y_continuous(expand = c(0,0),limits = c(0,29))+
  scale_colour_distiller(palette ="Reds", direction =1) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 13),     
        axis.text = element_text(size = 11, color = "black"),  
        axis.text.y = element_blank(),     
        axis.ticks.y = element_blank(),     
        legend.title = element_text(size = 13),     
        legend.text = element_text(size = 11))


#c4a_gui()
mycol<- c4a('rainbow_wh_rd', 25)
#绘图：
p2 <- ggplot(df, aes(x = x, next_x= next_x, node= node, next_node= next_node,
                     fill= node, label= node)) +
  geom_sankey(flow.alpha = 0.5,
              flow.fill = 'grey',
              flow.color = 'grey80', #条带描边色
              node.fill = mycol, #节点填充色
              smooth= 8,
              width= 0.08) +
  geom_sankey_text(size = 3.5, color= "black", hjust = 1)+
  theme_void() +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(), 
        legend.position = 'none',
        plot.margin = unit(c(0,5,0,0), units="cm"))


ggdraw() + draw_plot(p2) + draw_plot(p1, scale = 0.5, x = 0.33, y=-0.48, width=0.8, height = 1.85)





# fig 3c ------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(forcats)
library(readxl);library(writexl)
library(haven)
library(furrr)

pdata <- read_xlsx('pub_code/data/Figure 3.xlsx', sheet = 'figure3c')

cor.test(pdata$beta_sl, pdata$beta_hbc)
ggplot(pdata, aes(x=beta_sl, y=beta_hbc, fill = col_group))+
  geom_smooth(method = "lm", color = "black", fill = "gray", alpha = 0.3)+
  geom_point(size=3.5, pch=21, stroke = ifelse(pdata$sig_group == 'null', 0.3, 1.2)) +
  scale_fill_manual(values = cols, name = NULL)+ 
  annotate("text", x = 0.05, y = 0.75, size = 4.5, 
           label = "paste(italic(r), ' = 0.590, ', italic(P), ' < 0.001 (all shared metabolites)')", parse = TRUE)+
  #scale_x_continuous(limits=c(-0.4, 0.4))+
  #scale_y_continuous(limits=c(-0.65, 0.65))+
  theme_bw() +
  theme(plot.background = element_blank(),  # 移除图表背景
        panel.background = element_blank(),  # 移除绘图区域背景
        legend.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.25, 0.75),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour = 'black'))+
  guides(fill = guide_legend(ncol = 2))+
  geom_hline(aes(yintercept=0),color="darkgrey",linetype="dashed",size=0.5)+
  geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+ 
  xlab("Effect size (THSBC discovery)")+ylab("Effect size (HBC Validation)") 


# fig 3d ------------------------------------------------------------------

rm(list = ls())

df <- read_xlsx('pub_code/data/Figure 3.xlsx', sheet = 'figure3d')

df_unique = df %>% select(ft_name, diff, Cluster) %>% unique()
count = table(df_unique$diff, df_unique$Cluster)
ratio = count[2,]/ count[1,]
df_text = data.frame(Cluster = colnames(count), 
                     count0 = 47, count1 = count[2,]) %>% 
  mutate(ratio = count1/count0*100)

cluster_label = paste0(paste0('Cluster ',names(table(df$Cluster))), 
                       ' (n=',table(df$Cluster)/3,')')
names(cluster_label) = names(table(df$Cluster))

class_col = c('#6c91a7','#dfc858','#a1b383','#b491a2')
names(class_col) = names(cluster_label)
ggplot(df, aes(x = time, y = value))+
  facet_wrap(.~Cluster, labeller = labeller(Cluster=cluster_label))+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5, color='grey80')+
  geom_line(aes(group = ID, color = Cluster), alpha = 0.04)+
  geom_line(aes(group = Cluster, color = Cluster), data=cluster_traj, linewidth=0.8)+
  geom_text(aes(label = sprintf('Key DEMs: %d/%d (%.2f%%)', count1, count0, ratio)), 
            x = 3.1, y = 0.75, data = df_text, hjust = 1, size = 3)+
  scale_color_manual(values = class_col)+
  scale_x_continuous(expand = c(0,0.2), breaks = 1:3, labels = c('T1','T2','T3'))+
  labs(y = 'Standerdized value')+
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10.2, color = "black", 
                                    vjust = 1.9, hjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 10.2, color = "black"), 
        axis.text.y = element_text(size = 9.2, color = "grey10", 
                                   vjust = 0.5, hjust = 0.5, angle = 0),
        axis.ticks.length = unit(1.5, 'mm'),
        strip.text = element_text(size = 11, color = 'black'),
        strip.background = element_rect(fill = 'grey90'),
        legend.position = 'none')




# fig 3e ------------------------------------------------------------------

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

pdata <- read_xlsx('pub_code/data/Figure 3.xlsx', sheet = 'figure3e_thsbc')
pdata_both <- read_xlsx('pub_code/data/Figure 3.xlsx', sheet = 'figure3e_both')


p1 = ggplot(pdata, aes(x = beta, y = logq_discovery))+
  geom_point(pch=21, color="black", stroke=0.3,aes(size=abs(beta), fill=group))+
  geom_hline(aes(yintercept=1.30103),color="darkgrey",linetype="dashed",size=0.5)+
  geom_text_repel(aes(label=labels), size=3.5, direction="both", force)+
  #scale_x_continuous(limits=c(-0.01, 0.04))+
  labs(x = "Effect size", y = "-log(P value)") +
  theme_bw() +
  guides(size = FALSE,
         fill = guide_legend(override.aes = list(shape = 22, color = NA, size = 4))) +
  theme(plot.background = element_blank(),  
        panel.background = element_blank(),  
        panel.grid = element_blank(),
        legend.background = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour = 'black'),
        panel.spacing = unit(0,'mm'),
        legend.key.size = unit(14, "pt"),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.8)
  ) 
p1 

cor.test(pdata_both$thsbc_beta, pdata_both$hbc_beta)

p2 = ggplot(pdata_both, aes(x=thsbc_beta, y=hbc_beta))+
  geom_point(size=4,pch=21,color="black", aes(fill=group),
             stroke = ifelse(pdata_both$sig_group == 'null', 0.3, 1.2))+
  geom_hline(aes(yintercept=0),color="darkgrey",linetype="dashed",size=0.5)+
  geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+ 
  geom_text_repel(aes(label=labels),size=3.5, direction="both",force)+
  annotate("text", x = 0.5, y = 2.5, size = 3.5,
           label = "paste(italic(r), ' = 0.634, ', italic(P), ' < 0.001')", parse = TRUE)+
  labs(x = "Effect size (THSBC)", y = "Effect size (HBC)") +
  theme_bw() +
  theme(plot.background = element_blank(),  
        panel.background = element_blank(),  
        panel.grid = element_blank(),
        legend.background = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour = 'black'),
        panel.spacing = unit(0,'mm'),
        legend.position = 'none'
  ) 
p2 

library(patchwork)
(p1+ p2) + plot_layout(ncol = 2)



