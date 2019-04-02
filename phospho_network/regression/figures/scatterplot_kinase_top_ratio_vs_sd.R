# Yige Wu @ WashU 2018 Jan
# draw a grid showing the distribution of correlated kinase-substrate pairs are distributed across oncogenic pathways

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)

source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(pheatmap)
library(plyr)
library(grDevices)

lm_eqn = function(df){
  m = lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(P)~"="~p, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(summary(m)$coefficients[2,4], digits = 3)))
  as.character(as.expression(eq));                 
}

lm_eqn1 = function(df){
  m = lm(y ~ x, df);
  eq1 <- substitute(italic(y) == a + b %.% italic(x)*",", 
                    list(a = format(coef(m)[1], digits = 2), 
                         b = format(coef(m)[2], digits = 2), 
                         r2 = format(summary(m)$r.squared, digits = 3),
                         p = format(summary(m)$coefficients[2,4], digits = 3)))
  string2return <- as.character(as.expression(eq1))
  return(string2return)
}

lm_eqn2 = function(df){
  m = lm(y ~ x, df);
  eq2 <- substitute(~~italic(r)^2~"="~r2*","~~italic(P)~"="~p,
                    list(a = format(coef(m)[1], digits = 2),
                         b = format(coef(m)[2], digits = 2),
                         r2 = format(summary(m)$r.squared, digits = 3),
                         p = format(summary(m)$coefficients[2,4], digits = 3)))
  string2return <- as.character(as.expression(eq2))
  return(string2return)
}


# set variables -----------------------------------------------------------
reg_nonNA <- 20
# fdr_thres <- c(0.05, 0.1); names(fdr_thres) <- c("kinase", "phosphatase")
fdr_thres <- c(0.05, 0.05); names(fdr_thres) <- c("kinase", "phosphatase")
SELF = "trans"
num_top <- 10
# enzyme_type <- "kinase"
enzyme_type <- "phosphatase"

# input regression table --------------------------------------------------
regression_sup <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                       "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                        data.table = F)
regression_sup <- regression_sup %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression_sup <- adjust_regression_by_nonNA(regression = regression_sup, reg_nonNA = 20, reg_sig = reg_sig)
regression_sup <- annotate_ks_source(regression = regression_sup)

# for (SELF_tmp in c("cis", "trans")) {
for (SELF_tmp in c("trans")) {
    
  # first plot site-level comparison with known ones ------------------------
  for (enzyme_type_tmp in c("kinase")) {
    regression <- regression_sup %>%
      filter(enzyme_type == enzyme_type_tmp) %>%
      filter(SELF == SELF_tmp)
    
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
    write.table(x = tab_pairs.wratio, file = paste0(makeOutDir(resultD = resultD), enzyme_type_tmp, "_", SELF_tmp, "_correlated_ratio.txt"), quote = F, sep = "\t", row.names = F)
    ## merge into SDs of kinase phosphoprotein
    tab_GENE_sd <- unique(regression[, c("GENE", ifelse(SELF_tmp == "cis", "sd_pro_kin", "sd_pho_kin"), "sd_pho_sub", "Cancer", "Size")])  
    tab_GENE_sd$id <- paste0(tab_GENE_sd$GENE, ":", tab_GENE_sd$Cancer)
    tab_GENE_sd <- tab_GENE_sd[order(tab_GENE_sd[, ifelse(SELF_tmp == "cis", "sd_pro_kin", "sd_pho_kin")], decreasing = T),]
    tab_GENE_sd <- tab_GENE_sd[order(tab_GENE_sd$Size, decreasing = T),]
    tab_GENE_sd %>% head(n = 50)
    tab_GENE_sd <- tab_GENE_sd[!duplicated(tab_GENE_sd$id),]
    tab_GENE_sd %>% tail()
    
    tab_GENE_num_sub <- data.frame(table(unique(regression[, c("Cancer" , "GENE", "SUB_GENE")])[, c("Cancer", "GENE")]))
    
    # Plot range of the number of substrate phoshposites per kinase -----------------------------------
    tab2p <- tab_pairs.wratio
    tab2p <- tab2p[tab2p$Freq.regulatedT >= 5,]
    tab2p <- merge(tab2p, tab_GENE_sd, by = c("GENE", "Cancer"), all.x = T)
    tab2p <- merge(tab2p, tab_GENE_num_sub, by =  c("GENE", "Cancer"), all.x = T)
    tab2p %>% head()

    
    ## filtering
    tab2p <- tab2p[tab2p$ratio.regulatedT > 0,]
    tab2p$Cancer <- order_cancer(tab2p$Cancer)
    tab2p$y <- tab2p$ratio.regulatedT
    if (SELF_tmp == "trans") {
      x_var <- "sd_pho_kin"
    }
    if (SELF_tmp == "cis") {
      x_var <- "sd_pro_kin"
    }
    tab2p$x <- tab2p[, x_var]
    tab2p$x_is.top <- F
    tab2p$x_is.top <- get_top_QT_by_id(value_vector = tab2p$x, id_vector = tab2p$Cancer, qt = 0.75)
    tab2p$y_is.top <- get_top_QT_by_id(value_vector = tab2p$y, id_vector = tab2p$Cancer, qt = 0.75)
    
    eq1 <- ddply(tab2p,.(Cancer),lm_eqn1)
    eq1$x <- sapply(eq1$Cancer, FUN = function(cancer, tab2p) quantile(x = tab2p$x[tab2p$Cancer == cancer], probs = 0.9, na.rm = T), tab2p = tab2p)
    eq2 <- ddply(tab2p,.(Cancer),lm_eqn2)
    eq2$x <- eq1$x
    p <- ggplot(data = tab2p)
    p <- p + geom_smooth(formula = y ~ x, mapping = aes(x = x, y = y, alpha = 0.5), method = "lm", se=FALSE, color="blue", alpha = 0.5)
    p <- p + geom_point(mapping = aes(x = x, y = y))
    p <- p + scale_fill_manual(values = color_cancers2)
    p = p + geom_text(data=eq1, aes(x = 0.8*x, y = 0.95*max(tab2p$y, na.rm = T), label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
    p = p + geom_text(data=eq2, aes(x = 0.8*x, y = 0.9*max(tab2p$y, na.rm = T), label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
    p = p + geom_text_repel(data = tab2p[tab2p$y_is.top & tab2p$x_is.top,], 
                            mapping = aes(y = y, x = x, label = GENE), color = "black", segment.color = "black",
                            force = 3, segment.size = 1, segment.alpha = 0.5, size = 3, alpha = 0.8)
    p <- p + ylim(c(0,1.2*max(tab2p$y, na.rm = T)))
    p <- p + facet_grid(.~Cancer, scales = "free")
    p <- p + xlab(paste0("SD(", ifelse(SELF_tmp == "trans", "kinase phosphoprotein", "kinase protein"), ")")) + ylab("% correlated substrate per kinase")
    p = p + theme_bw()# + theme_nogrid()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), axis.ticks = element_blank(),
                   axis.title.y = element_text(size = 10), axis.title.x = element_text(size = 10),
                   strip.text.y = element_text(size = 10, angle = 0),
                   strip.background = element_rect(fill = "white", color = "white"),
                   panel.spacing.y = unit(0, "lines"),
                   panel.spacing.x = unit(0, "lines"))
    p
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, "_", x_var, '_vs_ratio.pdf',sep = "")
    ggsave(filename = fn, width = 15, height = 5)
    write.table(x = tab2p, file = paste0(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, "_for_subtype_plotting.txt"), quote = F, sep = "\t", row.names = F)
    
    if (SELF_tmp == "cis") {
      tab2p <- tab2p[tab2p$Cancer %in% c("BRCA", "UCEC", "CCRCC"),]
      tab2p$x_is.top <- F
      tab2p$x_is.top <- get_top_QT_by_id(value_vector = tab2p$x, id_vector = tab2p$Cancer, qt = 0.75)
      tab2p$y_is.top <- get_top_QT_by_id(value_vector = tab2p$y, id_vector = tab2p$Cancer, qt = 0.75)
      write.table(x = tab2p, file = paste0(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, "_for_subtype_plotting.txt"), quote = F, sep = "\t", row.names = F)
      
      
      eq1 <- ddply(tab2p,.(Cancer),lm_eqn1)
      eq1$x <- sapply(eq1$Cancer, FUN = function(cancer, tab2p) quantile(x = tab2p$x[tab2p$Cancer == cancer], probs = 0.9, na.rm = T), tab2p = tab2p)
      eq2 <- ddply(tab2p,.(Cancer),lm_eqn2)
      eq2$x <- eq1$x
      p <- ggplot(data = tab2p)
      p <- p + geom_smooth(formula = y ~ x, mapping = aes(x = x, y = y, alpha = 0.5), method = "lm", se=FALSE, color="blue", alpha = 0.5)
      p <- p + geom_point(mapping = aes(x = x, y = y))
      p <- p + scale_fill_manual(values = color_cancers2)
      p = p + geom_text(data=eq1, aes(x = x, y = 0.4, label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
      p = p + geom_text(data=eq2, aes(x = x, y = 0.3, label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
      
      p = p + geom_text_repel(data = tab2p[tab2p$y_is.top & tab2p$x_is.top,], 
                              mapping = aes(y = y, x = x, label = GENE), color = "black", segment.color = "black",
                              force = 3, segment.size = 1, segment.alpha = 0.5, size = 3, alpha = 0.8)
      p <- p + ylim(c(0,1.2))
      p <- p + facet_grid(.~Cancer, scales = "free")
      p <- p + xlab(paste0("SD(", ifelse(SELF_tmp == "trans", "kinase phosphoprotein", "kinase protein"), ")")) + ylab("% correlated substrate per kinase")
      p = p + theme_bw()# + theme_nogrid()
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), axis.ticks = element_blank(),
                     axis.title.y = element_text(size = 10), axis.title.x = element_text(size = 10),
                     strip.text.y = element_text(size = 10, angle = 0),
                     strip.background = element_rect(fill = "white", color = "white"),
                     panel.spacing.y = unit(0, "lines"),
                     panel.spacing.x = unit(0, "lines"))
      p
      fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, "_", x_var, '_vs_ratio_BRCA_UCEC_CCRCC.pdf',sep = "")
      ggsave(filename = fn, width = 8, height = 3, useDingbats = F)
    }
    
    tab2p$y <- tab2p$ratio.regulatedT
    tab2p$x <- tab2p$Size
    tab2p$x_is.top <- F
    tab2p$x_is.top <- get_top_QT_by_id(value_vector = tab2p$x, id_vector = tab2p$Cancer, qt = 0.75)
    tab2p$y_is.top <- get_top_QT_by_id(value_vector = tab2p$y, id_vector = tab2p$Cancer, qt = 0.75)
    
    eq1 <- ddply(tab2p,.(Cancer),lm_eqn1)
    eq1$x <- sapply(eq1$Cancer, FUN = function(cancer, tab2p) quantile(x = tab2p$x[tab2p$Cancer == cancer], probs = 0.9, na.rm = T), tab2p = tab2p)
    eq2 <- ddply(tab2p,.(Cancer),lm_eqn2)
    eq2$x <- eq1$x
    p <- ggplot(data = tab2p)
    p <- p + geom_smooth(formula = y ~ x, mapping = aes(x = x, y = y, alpha = 0.5), method = "lm", se=FALSE, color="blue", alpha = 0.5)
    p <- p + geom_point(mapping = aes(x = x, y = y))
    p <- p + scale_fill_manual(values = color_cancers2)
    p = p + geom_text(data=eq1, aes(x = 0.8*x, y = 0.4*max(tab2p$y, na.rm = T), label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
    p = p + geom_text(data=eq2, aes(x = 0.8*x, y = 0.3*max(tab2p$y, na.rm = T), label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
    p = p + geom_text_repel(data = tab2p[tab2p$y_is.top & tab2p$x_is.top,], 
                            mapping = aes(y = y, x = x, label = GENE), color = "black", segment.color = "black",
                            force = 3, segment.size = 1, segment.alpha = 0.5, size = 3, alpha = 0.8)
    p <- p + ylim(c(0,1.2*max(tab2p$y, na.rm = T)))
    p <- p + facet_grid(.~Cancer, scales = "free")
    p <- p + xlab("#samples with kinase phosphoprotein detected") + ylab("% correlated substrate per kinase")
    p = p + theme_bw()# + theme_nogrid()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), axis.ticks = element_blank(),
                   axis.title.y = element_text(size = 10), axis.title.x = element_text(size = 10),
                   strip.text.y = element_text(size = 10, angle = 0),
                   strip.background = element_rect(fill = "white", color = "white"),
                   panel.spacing.y = unit(0, "lines"),
                   panel.spacing.x = unit(0, "lines"))
    p
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_Size_vs_ratio.pdf',sep = "")
    ggsave(filename = fn, width = 15, height = 5)
    
    tab2p$y <- tab2p$ratio.regulatedT
    tab2p$x <- tab2p$Freq
    tab2p$x_is.top <- F
    tab2p$x_is.top <- get_top_QT_by_id(value_vector = tab2p$x, id_vector = tab2p$Cancer, qt = 0.75)
    tab2p$x_is.bottom <- get_bottom_QT_by_id(value_vector = tab2p$x, id_vector = tab2p$Cancer, qt = 0.25)
    tab2p$y_is.top <- get_top_QT_by_id(value_vector = tab2p$y, id_vector = tab2p$Cancer, qt = 0.75)
    
    eq1 <- ddply(tab2p,.(Cancer),lm_eqn1)
    eq1$x <- sapply(eq1$Cancer, FUN = function(cancer, tab2p) quantile(x = tab2p$x[tab2p$Cancer == cancer], probs = 0.9, na.rm = T), tab2p = tab2p)
    eq2 <- ddply(tab2p,.(Cancer),lm_eqn2)
    eq2$x <- eq1$x
    p <- ggplot(data = tab2p)
    p <- p + geom_smooth(formula = y ~ x, mapping = aes(x = x, y = y, alpha = 0.5), method = "lm", se=FALSE, color="blue", alpha = 0.5)
    p <- p + geom_point(mapping = aes(x = x, y = y))
    p <- p + scale_fill_manual(values = color_cancers2)
    p = p + geom_text(data=eq1, aes(x = 0.8*x, y = 0.4*max(tab2p$y, na.rm = T), label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
    p = p + geom_text(data=eq2, aes(x = 0.8*x, y = 0.3*max(tab2p$y, na.rm = T), label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
    p = p + geom_text_repel(data = tab2p[tab2p$y_is.top & tab2p$x_is.bottom,], 
                            mapping = aes(y = y, x = x, label = GENE), color = "black", segment.color = "black",
                            force = 3, segment.size = 1, segment.alpha = 0.5, size = 3, alpha = 0.8)
    p <- p + ylim(c(0,1.2*max(tab2p$y, na.rm = T)))
    p <- p + facet_grid(.~Cancer, scales = "free")
    p <- p + xlab("#samples with kinase phosphoprotein detected") + ylab("% correlated substrate per kinase")
    p = p + theme_bw()# + theme_nogrid()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), axis.ticks = element_blank(),
                   axis.title.y = element_text(size = 10), axis.title.x = element_text(size = 10),
                   strip.text.y = element_text(size = 10, angle = 0),
                   strip.background = element_rect(fill = "white", color = "white"),
                   panel.spacing.y = unit(0, "lines"),
                   panel.spacing.x = unit(0, "lines"))
    p
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_num_sub_vs_ratio.pdf',sep = "")
    ggsave(filename = fn, width = 15, height = 5)

    
  }
}

# tab2p$y <- tab2p$ratio.regulatedT
# tab2p$x <- tab2p$sd_pho_sub
# eq <- ddply(tab2p,.(Cancer), lm_eqn)
# p <- ggplot(data = tab2p)
# p <- p + geom_point(mapping = aes(x = x, y = y, fill = Cancer, group = Cancer))
# p <- p + scale_fill_manual(values = color_cancers2)
# p <- p + geom_smooth(mapping = aes(x = x, y = y), method = "lm", se=FALSE, color="black", formula = y ~ x)
# p = p + geom_text(data=eq, aes(x = 1, y = 0.3,label=V1), parse = TRUE, inherit.aes=FALSE, size = 3)
# p <- p + facet_grid(Cancer~.)
# p <- p + ylab("SD(substrate phosphorylation)") + xlab("% correlated substrate per kinase")
# p <- p + theme_bw()
# p
# fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_sd_pho_sub_vs_ratio.pdf',sep = "")
# ggsave(filename = fn, width = 6, height = 8)



