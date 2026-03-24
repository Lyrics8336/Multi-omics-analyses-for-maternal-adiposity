# fig 5a ------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(data.table)
library(forcats)
library(readxl);library(writexl)
library(haven)
library(ggplot2)

pdata_pbmi <- read_xlsx('pub_code/data/Figure 5.xlsx', sheet = 'figure5a_logic')
pdata_pbmi_linear <- read_xlsx('pub_code/data/Figure 5.xlsx', sheet = 'figure5a_linear')

p1 <- ggplot(pdata_pbmi, aes(color = subclass, y = beta, ymin = lci, ymax = uci, x = outcome))+
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 0, xmax = 4.5, alpha = 0.3, fill = "#fef4e8")+ #添加背景色块  
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 4.5, xmax = 9.5, alpha = 0.3, fill = "#bae7fa")+ #添加背景色块  
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 9.5, xmax = 10.5, alpha = 0.5, fill = "#dddcc7")+
  geom_pointrange(size=2, fatten = 1, position = position_dodge(width = 0.3))+  
  geom_text(
    data = subset(pdata_pbmi, pv_sig != "ns"), 
    aes(y = uci-0.15, x = outcome, label = pv_sig, color = subclass),
    vjust = -0.2, hjust = 0.5, size = 5, fontface = 'bold',
    position = position_dodge(width = 0.3)      # 与 pointrange 对齐
  ) +
  geom_hline(yintercept = 1, color = "black", size = 0.8, linetype = 3) +
  scale_color_manual(values = c("#f7d4bc", "#c0b9e6"), name = 'pBMI group')+ #设置颜色  
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), limits = c(-0.05, 4.05))+ #设置X轴刻度 
  labs(y = "Odds ratio", x = '') +
  theme_bw()+ 
  theme(
    plot.background = element_blank(),  # 移除图表背景
    panel.background = element_blank(),  # 移除绘图区域背景
    legend.background = element_blank(),  # 移除图例背景
    panel.grid = element_blank(), 
    axis.title.y = element_text(size = 13, vjust = 1.5),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(color = "black", size = 12), 
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = c(0.8, 0.8)
  ) 
p1

p3 <- ggplot(pdata_pbmi_linear, aes(color = subclass, y = beta, ymin = lci, ymax = uci, x = outcome))+
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 0, xmax = 3.5, alpha = 0.3, fill = "#bae7fa")+ #添加背景色块  
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 3.5, xmax = 4.5, alpha = 0.3, fill = "#dddcc7")+ #添加背景色块  
  #annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1.5, alpha = 0.5, fill = "#fef4e8")+
  geom_pointrange(size=2, fatten = 1, position = position_dodge(width = 0.3))+  
  geom_text(
    data = subset(pdata_pbmi_linear, pv_sig != "ns"), 
    aes(y = uci-0.03, x = outcome, label = pv_sig, color = subclass),
    vjust = -0.2, hjust = 0.5, size = 5, fontface = 'bold',
    position = position_dodge(width = 0.3)      # 与 pointrange 对齐
  ) +
  geom_hline(yintercept = 0, color = "black", size = 0.8, linetype = 3) +
  scale_color_manual(values = c("#f7d4bc", "#c0b9e6"))+ #设置颜色  
  scale_x_discrete(expand = c(0,0))+
  #scale_y_continuous(expand = c(0,0), limits = c(-0.02, 0.25))+ #设置X轴刻度 
  labs(y = "Coefficient", x = '') +
  theme_bw()+ 
  theme(
    plot.background = element_blank(),  # 移除图表背景
    panel.background = element_blank(),  # 移除绘图区域背景
    legend.background = element_blank(),  # 移除图例背景
    panel.grid = element_blank(), 
    axis.title.y = element_text(size = 13, vjust = 1.5),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(color = "black", size = 12), 
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = 'none'
  ) 
p3
p1 + p3 + plot_layout(widths = c(10, 4))


# fig 5b ------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(data.table)
library(forcats)
library(readxl);library(writexl)
library(haven)
library(ggtext)

pdata_gdmhdp <- read_xlsx('pub_code/data/Figure 5.xlsx', sheet = 'figure5b')

#文本旋转角度
myAngle <- c(seq(85, 85 - 360, length.out = 31)[1:13],
             seq(90, 90 - 360, length.out = 31)[14:31])


# myAngle[26:42] <- seq(88, 88 + 360 - 360 / 42, length.out = 42)[26:42]
mysize = c(seq(4,3,length.out = 11), rep(3, 2),
           seq(4,3,length.out = 11), rep(3, 7))

ggplot(pdata_gdmhdp, aes(x=seq))+
  annotate('text', x=22.5, y=c(12,22,32,42,52),
           label=c('10','20','30','40','50'), color="black")+
  geom_hline(yintercept=seq(0,50,by=10),colour="grey70",linewidth=0.3)+
  geom_col(aes(y = prop, fill= outcome), width = 0.9, alpha=0.8)+
  geom_text(aes(y = prop, label = label1),  
            angle = myAngle, vjust = 0.2, hjust = 0, size = 4 ) +
  geom_text(aes(y = prop, label = label2),color = 'white',
            angle = myAngle, vjust = 0.5, hjust = 1.2, size = mysize) +
  scale_fill_manual(values=c("#788FCE","#E6956F"))+
  scale_y_continuous(limits=c(-10,85), breaks=seq(0,80,20))+
  coord_polar(clip = 'off')+ 
  theme_void()+ 
  theme(panel.grid=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.background=element_rect(fill=NA,color=NA),
        plot.background=element_rect(fill=NA,color=NA),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position = c(0.65, 0.35))








# fig 5c ------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(data.table)
library(forcats)
library(readxl);library(writexl)
library(haven)
library(ggtext)


library(tidyverse)
library(viridis)
library(patchwork)
library(circlize)
library(networkD3)
library(htmlwidgets)
library(webshot2)

tmp <- read_xlsx('pub_code/data/Figure 5.xlsx', sheet = 'figure5c')

col <- c(
  rep("#FFB2B9",10),
  rep("#95C6C6",1),
  rep("#ebdcb2",3),
  rep("#d4cfdf",6),
  c("#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF")
)

ColourScal ='d3.scaleOrdinal() .range(["#FFB2B9", "#FFB2B9", "#FFB2B9", "#FFB2B9", "#FFB2B9", "#FFB2B9", 
"#FFB2B9", "#FFB2B9", "#FFB2B9", "#FFB2B9", "#95C6C6", "#ebdcb2", 
"#ebdcb2", "#ebdcb2", "#d4cfdf", "#d4cfdf", "#d4cfdf", "#d4cfdf", 
"#d4cfdf", "#d4cfdf", "#FDE725FF", "#B4DE2CFF", "#6DCD59FF", 
"#35B779FF", "#1F9E89FF", "#26828EFF"])'

## draw and save
p <- sankeyNetwork(Links = tmp, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "prop", NodeID = "name",
                   sinksRight = FALSE, colourScale = ColourScal, 
                   nodeWidth = 40, fontSize = 13, nodePadding = 5)



# fig 5d ------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(data.table)
library(forcats)
library(readxl);library(writexl)

library(circlize)


dir_forw_data <- read_xlsx('pub_code/data/Figure 5.xlsx', sheet = 'figure5d_forw')
dir_back_data <- read_xlsx('pub_code/data/Figure 5.xlsx', sheet = 'figure5d_back')


# 绘制"back"方向的条形图
plot_back <- ggplot(dir_back_data, aes(x = n, y = outcome, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) + 
  scale_x_continuous(limits = c(0, 16), expand = c(0.01, 0.01),breaks = seq(0,16,2))+
  theme_bw() +
  labs(x = "", y = "Backward") +
  theme(plot.background = element_blank(),  # 移除图表背景
        panel.background = element_blank(),  # 移除绘图区域背景
        legend.background = element_blank(),  # 移除图例背景
        strip.background = element_blank(), 
        panel.grid = element_blank(),
        legend.position = c(0.9, 0.3), 
        legend.title = element_blank(), 
        axis.text = element_text(color = 'black', hjust = 1, size = 12),
        axis.title.y = element_text(hjust = 0.5, size = 14)) 

# 绘制"forw"方向的条形图
plot_forw <- ggplot(dir_forw_data, aes(x = n, y = outcome, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +
  scale_x_continuous(limits = c(0, 16), expand = c(0.01, 0.01),breaks = seq(0,16,2))+
  theme_bw() +
  labs(x = "", y = "Forward") +
  theme(plot.background = element_blank(),  # 移除图表背景
        panel.background = element_blank(),  # 移除绘图区域背景
        legend.background = element_blank(),  # 移除图例背景
        strip.background = element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none", 
        axis.text = element_text(color = 'black', hjust = 1, size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5, size = 14)) 

# 使用patchwork包将两个图表纵向拼接
library(patchwork)
plot_forw / plot_back





# fig 5e ------------------------------------------------------------------

rm(list = ls())


Prepare_Micros <- function(input, dir){
  # input = forw_pbmi_owo_m36; dir = 'forw'
  
  data = input %>% arrange(type, name, Compounds)
  
  value <- 'prop.mediated'
  if(dir == 'forw'){
    name_x <- 'name'
    name_y <- 'Compounds'
    
    sector_name <- c(data[[name_x]], data[[name_y]])
    sec_source <- c(data$type, rep('Met', length(data[[name_y]]))) 
  } else if(dir == 'back'){
    name_x <- 'Compounds'
    name_y <- 'name'
    
    sector_name <- c(data[[name_x]], data[[name_y]])
    sec_source <- c(rep('Met', length(data[[name_x]])), data$type) 
  }
  
  
  
  dat <- data  %>%  mutate(trans = ifelse(pv_adj_FDR < 0.1, 0.1, 0.8)) %>% 
    select('from'=name_x, 'to'=name_y, 'value'=value, trans)
  dat$from <- as.character(dat$from)
  dat$to <- as.character(dat$to)
  
  
  gaps <- data.frame(
    name = sector_name,
    source = sec_source
  )
  gaps <- unique(gaps) %>% mutate(
    color=case_match(source, 'Bact'~'#FFB2B9', 'MetaG_g'~'#ebdcb2', 'MetaG_s'~'#d4cfdf', 'Met'~'#7EABCA')
  )
  
  n_micro_x <- dim(gaps  %>%  filter(source=='MetaG_s'))[1]
  n_micro_y <- dim(gaps  %>%  filter(source=='Bact'))[1]
  
  bgap <- 10
  gaps$gap <- 1
  gaps$gap[n_micro_x] <- bgap       # 根据个数更改
  gaps$gap[n_micro_x+n_micro_y] <- bgap 
  
  return(list('dat'=dat, 'gaps'=gaps))
}

# 自定义函数，用于将长文本分割成多行
split_text <- function(text, max_length = 32) {
  if (nchar(text) <= max_length) {
    return(text)
  } else {
    # 将文本分割成多个部分，每个部分不超过max_length长度
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
    # 添加最后一行
    lines <- c(lines, current_line)
    return(paste(lines, collapse = "\n"))
  }
}

Draw_chrod <- function(dat, gaps, dir, save_file, start_deg=85, hw=20,big_gap = 15){
  
  gaps <- unique(gaps)
  gap <- setNames(gaps$gap, gaps$name)
  grid.col <- setNames(gaps$color, gaps$name)
  group <- gaps  %>%  select(name, source)  %>%  arrange(source)
  group <- setNames(group$source, group$name)
  
  if(dir == 'forw'){
    arrow_color <- dat  %>%  left_join(gaps, by=c('from'='name'))  %>%  select(color)
    arrow_color <- arrow_color[[1]]
  } else if(dir == 'back'){
    arrow_color <- dat  %>%  left_join(gaps, by=c('to'='name'))  %>%  select(color)
    arrow_color <- arrow_color[[1]]
  }
  
  pdf(save_file, width = hw, height = hw)
  circos.clear()
  circos.par(start.degree = start_deg, 
             # gap.degree = 2, 
             # gap.after=gap
             canvas.xlim = c(-1.5, 1.5), # 扩大 X 轴画布范围
             canvas.ylim = c(-1.5, 1.5),  # 扩大 Y 轴画布范围
             points.overflow.warning = FALSE
             
  )
  
  chordDiagram(dat,
               big.gap = big_gap, small.gap = 2,
               grid.col = grid.col,
               col = arrow_color,
               #group = group,
               #transparency = 0.25,
               transparency = dat$trans,
               directional = 1,
               direction.type = c("arrows", "diffHeight"), 
               diffHeight  = -0.02,
               annotationTrack = "grid", 
               annotationTrackHeight = c(0.08, 0.1),
               link.arr.type = "big.arrow", 
               link.arr.length = 0.03,
               #  link.sort = TRUE, 
               #  link.decreasing = TRUE,
               link.largest.ontop = TRUE,
               preAllocateTracks = list(
                 track.height = 0.1,
                 track.margin = c(0.01, 0)
               ))
  
  circos.trackPlotRegion(
    track.index = 2, 
    bg.border = NA, 
    panel.fun = function(x, y) {
      # xlim = get.cell.meta.data("xlim")
      # sector.index = get.cell.meta.data("sector.index")
      # sector.index = gsub("[a-z]+_", "", get.cell.meta.data("sector.index"))
      sector.index = sapply(CELL_META$sector.index, split_text)
      
      # Add names to the sector. 
      # circos.text(
      #   x = mean(xlim), 
      #   y = 0.5, 
      #   col = "white",
      #   labels = sector.index, 
      #   facing = "bending", 
      #   cex = 1,
      #   niceFacing = TRUE
      # )
      circos.text(
        x = CELL_META$xcenter,  # 文字切向位置，按角度给
        y = 1.2,                # 文字径向位置
        labels = sector.index, # 需要显示的字
        facing = "clockwise",   # 字体排列方式
        niceFacing = TRUE, 
        cex = 2,        # 字体大小
        col = "black",  # 字体颜色
        adj = c(0, 0.5))   # 字体位置
    }
  )
  
  dev.off()
  print(save_file)
}



forw_pbmi_GDM <- read_xlsx('pub_code/data/Figure 5.xlsx', sheet = 'figure5e_forw')
back_pbmi_GDM <- read_xlsx('pub_code/data/Figure 5.xlsx', sheet = 'figure5e_back')

## pbmi-forw-GDM
DATA <- Prepare_Micros(forw_pbmi_GDM, 'forw')

Draw_chrod(DATA$dat, DATA$gaps, start_deg=270, dir = 'forw',
           save_file = ' GDM_forw_micro_metab.pdf')
## pbmi-back-GDM
back_pbmi_GDM$name[which(back_pbmi_GDM$name == 'g__Blautia'&back_pbmi_GDM$type=='MetaG_g')] <- 'g__Blautia_mg'
DATA <- Prepare_Micros(back_pbmi_GDM, 'back')

Draw_chrod(DATA$dat, DATA$gaps, start_deg=270, dir = 'back',
           save_file = 'GDM_back_micro_metab.pdf')









