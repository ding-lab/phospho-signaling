# Yige Wu @ WashU 2018 Jan
# draw a grid showing the distribution of correlated kinase-substrate pairs are distributed across oncogenic pathways

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(pheatmap)



# set variables -----------------------------------------------------------
reg_nonNA <- 20
fdr_thres <- c(0.05, 0.1); names(fdr_thres) <- c("kinase", "phosphatase")
SELF = "trans"
num_top <- 10
# enzyme_type <- "kinase"
enzyme_type <- "phosphatase"

for (enzyme_type in c("kinase")) {
  # inputs ------------------------------------------------------------------
  regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                     "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
  ## filtering
  regression <- regression[regression$enzyme_type == enzyme_type & regression$SELF == SELF,]
  regression <- markSigKS(regression = regression, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
  
  # get table of how many pairs are annotated to different pathways ---------
  tab_pairs <- data.frame(table(regression[, c("GENE", "Cancer", "regulated")]))
  tab_pairs %>%
    head()
  
  tab_pairs.wratio <- merge(tab_pairs[tab_pairs$regulated == "TRUE",], 
                            tab_pairs[tab_pairs$regulated == "FALSE", c("GENE", "Cancer", "Freq")],
                            by = c("GENE", "Cancer"), suffixes = c(".regulatedT", ".regulatedF"), all.y = T)
  tab_pairs.wratio %>%
    head()
  tab_pairs.wratio$ratio.regulatedT <- tab_pairs.wratio$Freq.regulatedT/(tab_pairs.wratio$Freq.regulatedT +  tab_pairs.wratio$Freq.regulatedF)
  tab_pairs.wratio %>%
    head()
  ## merge into SDs of kinase phosphoprotein
  tab_GENE_sd <- unique(regression[, c("GENE", "sd_pho_kin", "sd_pho_sub", "Cancer", "Size")])  
  tab_GENE_sd$id <- paste0(tab_GENE_sd$GENE, ":", tab_GENE_sd$Cancer)
  tab_GENE_sd <- tab_GENE_sd[order(tab_GENE_sd$sd_pho_kin, decreasing = T),]
  tab_GENE_sd <- tab_GENE_sd[order(tab_GENE_sd$Size, decreasing = T),]
  tab_GENE_sd %>% head(n = 50)
  tab_GENE_sd <- tab_GENE_sd[!duplicated(tab_GENE_sd$id),]
  tab_GENE_sd %>% tail()
  
  # Plot range of the number of substrate phoshposites per kinase -----------------------------------
  tab2p <- tab_pairs.wratio
  tab2p <- tab2p[tab2p$Freq.regulatedF > 0,]
  tab2p <- merge(tab2p, tab_GENE_sd, by = c("GENE", "Cancer"), all = T)
  tab2p %>% head()
  
  lm_eqn = function(df){
    m = lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(P)~"="~p, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3),
                          p = format(summary(m)$coefficients[2,4], digits = 3)))
    as.character(as.expression(eq));                 
  }
  
  ## filtering
  tab2p <- tab2p[tab2p$ratio.regulatedT > 0,]
  tab2p$y <- tab2p$ratio.regulatedT
  tab2p$x <- tab2p$sd_pho_kin
  tab2p$Cancer <- order_cancer(tab2p$Cancer)
  eq <- ddply(tab2p,.(Cancer),lm_eqn)
  p <- ggplot(data = tab2p)
  p <- p + geom_point(mapping = aes(x = x, y = y, fill = Cancer, group = Cancer))
  p <- p + scale_fill_manual(values = color_cancers2)
  p <- p + geom_smooth(mapping = aes(x = x, y = y), method = "lm", se=FALSE, color="black", formula = y ~ x)
  p = p + geom_text(data=eq, aes(x = 1, y = 0.4, label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
  p <- p + facet_grid(Cancer~.)
  p <- p + theme_bw()
  p <- p + ylab("SD(kinase phosphoprotein)") + xlab("% correlated substrate per kinase")
  p
  fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_sd_pho_kin_vs_ratio.pdf',sep = "")
  ggsave(filename = fn, width = 6, height = 8)
  
  tab2p$y <- tab2p$ratio.regulatedT
  tab2p$x <- tab2p$Size
  eq <- ddply(tab2p,.(Cancer),lm_eqn)
  
  p <- ggplot(data = tab2p)
  p <- p + geom_point(mapping = aes(x = x, y = y, fill = Cancer, group = Cancer))
  p <- p + scale_fill_manual(values = color_cancers2)
  p <- p + geom_smooth(mapping = aes(x = x, y = y), method = "lm", se=FALSE, color="black", formula = y ~ x)
  p = p + geom_text(data=eq, aes(x = 60, y = 0.3,label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
  p <- p + facet_grid(Cancer~.)
  p <- p + ylab("number of samples with available data") + xlab("% correlated substrate per kinase")
  p <- p + theme_bw()
  p
  fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_Size_vs_ratio.pdf',sep = "")
  ggsave(filename = fn, width = 6, height = 8)
  
  tab2p$y <- tab2p$ratio.regulatedT
  tab2p$x <- tab2p$sd_pho_sub
  eq <- ddply(tab2p,.(Cancer), lm_eqn)
  p <- ggplot(data = tab2p)
  p <- p + geom_point(mapping = aes(x = x, y = y, fill = Cancer, group = Cancer))
  p <- p + scale_fill_manual(values = color_cancers2)
  p <- p + geom_smooth(mapping = aes(x = x, y = y), method = "lm", se=FALSE, color="black", formula = y ~ x)
  p = p + geom_text(data=eq, aes(x = 1, y = 0.3,label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
  p <- p + facet_grid(Cancer~.)
  p <- p + ylab("SD(substrate phosphorylation)") + xlab("% correlated substrate per kinase")
  p <- p + theme_bw()
  p
  fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_sd_pho_sub_vs_ratio.pdf',sep = "")
  ggsave(filename = fn, width = 6, height = 8)

}





