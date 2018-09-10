### common libs and dependencies ### 
# especially for plotting #

# dependencies
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library("gplots")
library(gridExtra)
library(matrixStats)
library("WGCNA")
library(grid)
# library("PerformanceAnalytics")
# library("lattice")
# library(VennDiagram)
# library(grDevices)


col_paletteB = colorRampPalette(brewer.pal(9,"Blues"))
col_paletteR = colorRampPalette(brewer.pal(9,"Reds"))
RdBu = brewer.pal(9, "RdBu") 
getPalette = colorRampPalette(rev(RdBu))
RdBu1024 = colorRampPalette(rev(RdBu))(1024)
YlGnBu = brewer.pal(9, "YlGnBu") 
getPalette2 = colorRampPalette(YlGnBu)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
set1 = brewer.pal(9,"Set1")
set2 = brewer.pal(9,"Set2")

col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals) ))

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())

theme0 = function(...) theme( legend.position = "none",
                               panel.background = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               #panel.margin = unit(0,"null"),
                               axis.ticks = element_blank(),
                               axis.text.x = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               axis.ticks.length = unit(0,"null"),
                               axis.ticks.margin = unit(0,"null"),
                               panel.border=element_rect(color=NA),...)

colors = c(
  "#C1A72F",
  "#FAD2D9",
  "#ED2891",
  "#F6B667",
  "#104A7F",
  "#9EDDF9",
  "#3953A4",
  "#007EB5",
  "#B2509E",
  "#97D1A9",
  "#ED1C24",
  "#F8AFB3",
  "#EA7075",
  "#754C29",
  "#D49DC7",
  "#CACCDB",
  "#D3C3E0",
  "#A084BD",
  "#542C88",
  "#D97D25",
  "#6E7BA2",
  "#E8C51D",
  "#7E1918",
  "#DAF1FC",
  "#00A99D",
  "#BBD642",
  "#00AEEF",
  "#BE1E2D",
  "#F9ED32",
  "#CEAC8F",
  "#FBE3C7",
  "#F89420",
  "#009444")    
color.names = c("ACC",
                "BLCA",
                "BRCA",
                "CESC",
                "CHOL",
                "COAD",
                "DLBC",
                "ESCA",
                "GBM",
                "HNSC",
                "KICH",
                "KIRC",
                "KIRP",
                "LAML",
                "LGG",
                "LIHC",
                "LUAD",
                "LUSC",
                "MESO",
                "OV",
                "PAAD",
                "PCPG",
                "PRAD",
                "READ",
                "SARC",
                "SKCM",
                "STAD",
                "TGCT",
                "THCA",
                "THYM",
                "UCEC",
                "UCS",
                "UVM")
names(colors) = color.names