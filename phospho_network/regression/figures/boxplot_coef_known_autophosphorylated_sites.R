# Yige Wu @ WashU 2019 Feb
# venn diagram showing whether regulated enzyme-substrate-phosphosite pairs are also the sites previously reported

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# set variables -----------------------------------------------------------
reg_nonNA <- 20
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
SELF = "trans"

# inputs regression result ------------------------------------------------------------------

# regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
#                                    "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
# ## filtering
# regression <- regression[regression$enzyme_type_tmp == enzyme_type_tmp & regression$SELF == SELF,]
# regression <- markSigKS(regression = regression, sig_thres = fdr_thres[enzyme_type_tmp], enzyme_type_tmp = enzyme_type_tmp)

regression_sup <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                       "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                        data.table = F)
regression_sup <- regression_sup %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression_sup <- adjust_regression_by_nonNA(regression = regression_sup, reg_nonNA = 20, reg_sig = reg_sig)
regression_sup <- annotate_ks_source(regression = regression_sup)

regression_sup %>%
  filter(SELF == "cis") %>%
  filter(enzyme_type == "kinase") %>%
  filter(fdr_sig == T) %>%
  filter(coef_pro_kin < 0) %>%
  filter(is.direct == T) %>%
  select(GENE) %>%
  unique() %>%
  nrow()

for (SELF_tmp in c("cis")) {
  # first plot site-level comparison with known ones ------------------------
  for (enzyme_type_tmp in c("kinase")) {
    regression <- regression_sup %>%
      filter(enzyme_type == enzyme_type_tmp) %>%
      filter(SELF == SELF_tmp) %>%
      mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD))
    
    regression <- regression %>%
      filter(pair_pro %in% regression$pair_pro[regression$regulated & regression$is.direct])
    
    tab2p <- regression
    tab2p$is.top <- F
    tab2p$is.top[tab2p$regulated] <- get_top_QT_by_id(value_vector = tab2p$coef_pro_kin[tab2p$regulated], id_vector = tab2p$pair_pro_cancer[tab2p$regulated], qt = 0.9)
    tab2p$Cancer <- order_cancer(tab2p$Cancer)
    pos <- position_jitter(width = 0.2, seed = 1)
    
    p <- ggplot()
    # p = p + geom_boxplot(data = tab2p, mapping = aes(x = GENE, y = coef_pro_kin, fill = Cancer), color = NA, alpha = 0.4, notch = F)
    p = p + geom_point(data = tab2p, mapping = aes(y = GENE, x = coef_pro_kin, shape = is.direct, color = regulated, group = is.direct), 
                       fill = "black", stroke = 0, alpha = 0.6, size = 2) #position = pos, 
    p <- p + scale_color_manual(values = c("TRUE" = set1[1], "FALSE" = "grey50"))
    p <- p + scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 17))
    p = p + geom_text_repel(data = tab2p[tab2p$is.direct & tab2p$is.top,],
                            mapping = aes(y = GENE, x = coef_pro_kin, label = phosphosite), color = "black", segment.color = "black",
                            force = 2, segment.size = 2, segment.alpha = 0.8, size = 3, alpha = 0.8)
    p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey50", alpha = 0.5)
    p <- p + xlab(label = paste0("beta correlation coefficient"))
    p <- p + ylab(label = paste0("autophosphorylated kinases"))
    p <- p + facet_grid(.~Cancer, scales = "free", space = "fixed")
    p = p + theme_bw()# + theme_nogrid()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), axis.ticks = element_blank(),
                   axis.title.y = element_text(size = 10), axis.title.x = element_text(size = 10),
                   strip.text.y = element_text(size = 10, angle = 0),
                   strip.background = element_rect(fill = "white", color = "white"),
                   panel.spacing.y = unit(0, "lines"),
                   panel.spacing.x = unit(0, "lines"))
    p
    fn <- paste0(makeOutDir(resultD = resultD), "autophosphorylated_kinase_coefs.pdf")
    ggsave(filename = fn, width = 12, height = 8)
    
    tab2p <- regression
    tab2p$is.top <- F
    tab2p$is.top[tab2p$regulated] <- get_top_QT_by_id(value_vector = tab2p$coef_pro_kin[tab2p$regulated], id_vector = tab2p$pair_pro_cancer[tab2p$regulated], qt = 0.9)
    tab2p$Cancer <- order_cancer(tab2p$Cancer)
    tab2p <- tab2p[tab2p$GENE %in% driver_genes$Gene,]
    pos <- position_jitter(width = 0.2, seed = 1)
    pos_up <- position_nudge(x = 0, y = 0.2)
    pos_down <- position_nudge(x = 0, y = -0.2)
    
    p <- ggplot()
    p = p + geom_point(data = tab2p[tab2p$is.direct,], mapping = aes(y = GENE, x = coef_pro_kin, shape = is.direct, color = regulated, group = is.direct), 
                       position = pos_up, fill = "black", stroke = 0, alpha = 0.6, size = 2)
    p = p + geom_point(data = tab2p[!tab2p$is.direct,], mapping = aes(y = GENE, x = coef_pro_kin, shape = is.direct, color = regulated, group = is.direct), 
                       position = pos_down, fill = "black", stroke = 0, alpha = 0.6, size = 2) #  
    p <- p + scale_color_manual(values = c("TRUE" = set1[1], "FALSE" = "grey50"))
    p <- p + scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 17))
    p = p + geom_text_repel(data = tab2p[tab2p$is.direct & tab2p$is.top,],
                            mapping = aes(y = GENE, x = coef_pro_kin, label = phosphosite), color = "black", segment.color = "black",
                            min.segment.length = unit(0, 'lines'), 
                            force = 4, segment.size = 0.5, segment.alpha = 0.8, size = 2.5, alpha = 0.8, position = pos_up)
    p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey50", alpha = 0.5)
    p <- p + xlab(label = paste0("beta correlation coefficient"))
    p <- p + ylab(label = paste0("autophosphorylated kinases"))
    p <- p + facet_grid(.~Cancer, scales = "free", space = "fixed")
    p = p + theme_bw()# + theme_nogrid()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), axis.ticks = element_blank(),
                   axis.title.y = element_text(size = 10), axis.title.x = element_text(size = 10),
                   strip.text.y = element_text(size = 10, angle = 0),
                   strip.background = element_rect(fill = "white", color = "white"),
                   panel.spacing.y = unit(0, "lines"),
                   panel.spacing.x = unit(0, "lines"))
    p
    fn <- paste0(makeOutDir(resultD = resultD), "autophosphorylated_kinase_coefs_driver_gene_highlighted.pdf")
    ggsave(filename = fn, width = 8, height = 3)
    stop()
  }
}

stop("")

## no difference in correlation coeffcient overall
wilcox.test(regression$coef_pro_kin[regression$regulated & regression$is.direct], regression$coef_pro_kin[regression$regulated & !regression$is.direct])

## no difference in correlation coeffcient in each cancer type
for (cancer_tmp in unique(regression$Cancer)) {
  wilcox_tmp <- wilcox.test(regression$coef_pro_kin[regression$regulated & regression$is.direct & regression$Cancer == cancer_tmp], regression$coef_pro_kin[regression$regulated & !regression$is.direct & regression$Cancer == cancer_tmp])
  print(wilcox_tmp)
}

## no difference in correlation coeffcient in each cancer type in each gene
for (cancer_tmp in unique(regression$Cancer)) {
  for (gene_tmp in unique(regression$GENE)) {
    regression_tmp <- regression %>%
      filter(regulated == T) %>%
      filter(Cancer == cancer_tmp) %>%
      filter(GENE == gene_tmp)
    if (nrow(regression_tmp[regression_tmp$is.direct,]) < 4) {
      next()
    }
    wilcox_tmp <- wilcox.test(regression_tmp$coef_pro_kin[regression_tmp$is.direct], regression_tmp$coef_pro_kin[!regression_tmp$is.direct])
    if (wilcox_tmp$p.value < 0.05) {
      print(wilcox_tmp)
    }
  }
}

## 75% quantile
regression$above_75 <- F
regression$above_75[regression$regulated] <- get_top_QT_by_id(value_vector = regression$coef_pro_kin[regression$regulated], id_vector = regression$pair_pro_cancer[regression$regulated], qt = 0.75)
nrow(regression[regression$above_75 & regression$is.direct,])
nrow(regression[regression$regulated & regression$is.direct,])
nrow(regression[regression$above_75 & regression$is.direct,])/nrow(regression[regression$regulated & regression$is.direct,])

## 90% quantile
regression$above_90 <- F
regression$above_90[regression$regulated] <- get_top_QT_by_id(value_vector = regression$coef_pro_kin[regression$regulated], id_vector = regression$pair_pro_cancer[regression$regulated], qt = 0.9)
nrow(regression[regression$above_90 & regression$is.direct,])
nrow(regression[regression$regulated & regression$is.direct,])
nrow(regression[regression$above_90 & regression$is.direct,])/nrow(regression[regression$regulated & regression$is.direct,])

## 50% quantile
regression$above_50 <- F
regression$above_50[regression$regulated] <- get_top_QT_by_id(value_vector = regression$coef_pro_kin[regression$regulated], id_vector = regression$pair_pro_cancer[regression$regulated], qt = 0.5)
nrow(regression[regression$above_50 & regression$is.direct,])
nrow(regression[regression$regulated & regression$is.direct,])
nrow(regression[regression$above_50 & regression$is.direct,])/nrow(regression[regression$regulated & regression$is.direct,])

## 95% quantile
regression$above_95 <- F
regression$above_95[regression$regulated] <- get_top_QT_by_id(value_vector = regression$coef_pro_kin[regression$regulated], id_vector = regression$pair_pro_cancer[regression$regulated], qt = 0.95)
nrow(regression[regression$above_95 & regression$is.direct,])
nrow(regression[regression$regulated & regression$is.direct,])
nrow(regression[regression$above_95 & regression$is.direct,])/nrow(regression[regression$regulated & regression$is.direct,])
test <- regression[regression$above_95 & regression$is.direct,]
test <- regression[regression$above_90 & regression$is.direct,]

for (SELF_tmp in c("cis")) {
  # first plot site-level comparison with known ones ------------------------
  for (enzyme_type_tmp in c("kinase")) {
    regression <- regression_sup %>%
      filter(enzyme_type == enzyme_type_tmp) %>%
      filter(SELF == SELF_tmp) %>%
      mutate(phosphosite = paste0(SUB_GENE, "\n", SUB_MOD_RSD))
    
    regression <- regression %>%
      filter(pair_pro_cancer %in% regression$pair_pro_cancer[regression$regulated & regression$is.direct])
    
    regression %>%
      select(GENE) %>%
      unique() %>%
      nrow()
    
    regression %>%
      select(Cancer) %>%
      unique() %>%
      nrow()
    
    regression %>%
      filter(is.direct == T & regulated == T) %>%
      unique() %>%
      nrow()
    
    regression %>%
      filter(is.direct == T & regulated == T) %>%
      select(pair) %>%
      unique() %>%
      nrow()
    
    regression %>%
      filter(is.direct == T & regulated == T) %>%
      select(pair) %>%
      unique()
  }
}
