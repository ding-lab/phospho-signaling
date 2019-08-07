# Yige Wu @ WashU 2018 Jan
# draw volcano plots for regression result (3 cancer types together)

# Yige Wu @ WashU 2018 Apr
## draw volcano plots 

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
library(ggrepel)
## input enzyme substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv"))
op_enzyme <- fread(input = paste0(ppnD, "diffexp/tables/integrate_enzyme_KSEAordiffexp_sub_FC_plus_regression/omnipath_enzyme.txt"), data.table = F)
color_map <- c("firebrick1", "black")
names(color_map) <- c("TRUE", "FALSE")
color_direction <- c(set1[1], set1[2], "black")
names(color_direction) <- c("up", "down", "NA")
fdr_pp <- 0.2


# kinases -----------------------------------------------------------------
for (cancer in c("BRCA", "OV", "CO")) {
  ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE,":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  
  for (enzyme_type in c("kinase")) {
    subdir1 <- paste0(makeOutDir(), enzyme_type, "/")
    dir.create(subdir1)
    
    for (self in c("trans")) {
      df0 <- ksea_diffexp[ksea_diffexp$enzyme_type == enzyme_type & !is.na(ksea_diffexp$FDR_pho_kin),]
      df1 <- df0[-log10(df0$FDR_pho_kin) < 30,]
      df1$coef_pho_kin_filtered = remove_outliers(df1$coef_pho_kin, out_thres = 1.5)
      df <- unique(df1[,c("coef_pho_kin_filtered", "FDR_pho_kin", "pair")])
      df$pass_FDR <- ifelse(df$coef_pho_kin_filtered > 0 & df$FDR_pho_kin < 0.1, TRUE, FALSE)
      
      for (text_cutoff in c(1)) {
        ## get KSEA scores integrated with differential phosphorylation results
        p = ggplot(df)
        p = p + geom_point(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin), 
                               color = pass_FDR), stroke = 0, alpha = 0.05)
        p = p + theme_bw() + theme_nogrid()
        p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
        p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
        p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
        p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
        p = p + scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey")) 
        p = p + theme(legend.position="none")
        p
        subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
        dir.create(subdir2)
        fn = paste(subdir2, "volcano_trans_", cancer , "_", enzyme_type, '_substrate_textcutoff', text_cutoff, '.pdf',sep ="")
        ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
      }
    }
  }
}
for (cancer in c("BRCA", "OV", "CO")) {
  ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE,":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  
  for (enzyme_type in c("kinase")) {
    subdir1 <- paste0(makeOutDir(), enzyme_type, "/")
    dir.create(subdir1)
    
    for (self in c("cis")) {
      df0 <- ksea_diffexp[ksea_diffexp$enzyme_type == enzyme_type & !is.na(ksea_diffexp$FDR_pro_kin),]
      df1 <- df0[-log10(df0$FDR_pro_kin) < 30,]
      df1$coef_pro_kin_filtered = remove_outliers(df1$coef_pro_kin, out_thres = 1.5)
      df <- unique(df1[,c("coef_pro_kin_filtered", "FDR_pro_kin", "pair")])
      df$pass_FDR <- ifelse(df$coef_pro_kin_filtered > 0 & df$FDR_pro_kin < 0.1, TRUE, FALSE)
      
      for (text_cutoff in c(1)) {
        ## get KSEA scores integrated with differential phosphorylation results
        p = ggplot(df)
        p = p + geom_point(aes(x=coef_pro_kin_filtered, y=-log10(FDR_pro_kin), 
                               color = pass_FDR), stroke = 0, alpha = 0.1)
        p = p + theme_bw() + theme_nogrid()
        p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
        p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
        p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
        p = p + labs(x = paste0("Coefficient for ", enzyme_type, " protein level"), y="-log10(FDR)")
        p = p + scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey")) 
        p = p + theme(legend.position="none")
        p
        subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
        dir.create(subdir2)
        fn = paste(subdir2, "volcano_cis_", cancer , "_", enzyme_type, '_substrate_textcutoff', text_cutoff, '.pdf',sep ="")
        ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
      }
    }
  }
}


# phosphotase -----------------------------------------------------------------

for (cancer in c("BRCA", "OV", "CO")) {
# for (cancer in c("BRCA")) {
  ksea_diffexp <- fread(input = paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE,":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  
  for (enzyme_type in c("phosphotase")) {
    subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
    dir.create(subdir1)
    
    for (self in c("trans")) {
      df0 <- ksea_diffexp[ksea_diffexp$SELF == self & ksea_diffexp$aa_reported & ksea_diffexp$enzyme_type == enzyme_type & !is.na(ksea_diffexp$FDR_pho_kin),]
      # df1 <- df0[-log10(df0$FDR_pho_kin) < 30,]
      df1 <- df0
      df1$coef_pho_kin_filtered = remove_outliers(df1$coef_pho_kin, out_thres = 1.5)
      df <- unique(df1[,c("coef_pho_kin_filtered", "FDR_pho_kin.aa_reported", "P_pho_kin", "pair")])
      df <- df[!is.na(df$coef_pho_kin_filtered) & !is.na(df$FDR_pho_kin.aa_reported),]
      df$pass_FDR <- ifelse(df$coef_pho_kin_filtered < 0 & df$FDR_pho_kin.aa_reported < fdr_pp, TRUE, FALSE)
      df$pass_pvalue <- ifelse(df$coef_pho_kin_filtered < 0 & df$P_pho_kin < fdr_pp, TRUE, FALSE)
      
      for (text_cutoff in c(-log10(fdr_pp))) {
        ## get KSEA scores integrated with differential phosphorylation results
        p = ggplot(df)
        p = p + geom_point(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin.aa_reported), 
                               color = pass_FDR), stroke = 0, alpha = 0.5)
        p = p + theme_bw() + theme_nogrid()
        p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
        p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
        p = p + geom_text_repel(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin.aa_reported), 
                                    label= ifelse(-log10(FDR_pho_kin.aa_reported) > text_cutoff & coef_pho_kin_filtered > 0, as.character(pair), NA)), 
                                size=2, color = "grey", alpha = 0.5)
        p = p + geom_text_repel(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin.aa_reported), 
                                    label= ifelse(-log10(FDR_pho_kin.aa_reported) > text_cutoff & coef_pho_kin_filtered < 0, as.character(pair), NA)),
                                size=3, color = "black", alpha = 1, force = 2)
        p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
        p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
        p = p + scale_color_manual(values = c("TRUE" = "#377EB8", "FALSE" = "grey", "NA" = "grey")) 
        p = p + theme(legend.position="none")
        p
        subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
        dir.create(subdir2)
        fn = paste(subdir2, "volcano_trans_", cancer , "_", enzyme_type, '_substrate_textcutoff', text_cutoff, '.pdf',sep ="")
        ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
      }
    }
  }
}

df_cans <- NULL
for (cancer in c("BRCA", "OV")) {
  # for (cancer in c("BRCA")) {
  ksea_diffexp <- fread(input = paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE,":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  
  for (enzyme_type in c("phosphotase")) {
    subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
    dir.create(subdir1)
    
    for (self in c("trans")) {
      df0 <- ksea_diffexp[ksea_diffexp$SELF == self & ksea_diffexp$aa_reported & ksea_diffexp$enzyme_type == enzyme_type & !is.na(ksea_diffexp$FDR_pho_kin),]
      # df1 <- df0[-log10(df0$FDR_pho_kin) < 30,]
      df1 <- df0
      df1$coef_pho_kin_filtered = remove_outliers(df1$coef_pho_kin, out_thres = 1.5)
      df <- unique(df1[,c("coef_pho_kin_filtered", "FDR_pho_kin.aa_reported", "P_pho_kin", "pair")])
      df <- df[!is.na(df$coef_pho_kin_filtered) & !is.na(df$FDR_pho_kin.aa_reported),]
      df$pass_FDR <- ifelse(df$coef_pho_kin_filtered < 0 & df$FDR_pho_kin.aa_reported < fdr_pp, TRUE, FALSE)
      df$pass_pvalue <- ifelse(df$coef_pho_kin_filtered < 0 & df$P_pho_kin < fdr_pp, TRUE, FALSE)
      df$Cancer <- cancer
    }
  }
  df_cans <- rbind(df_cans, df)
}

for (text_cutoff in c(-log10(fdr_pp))) {
  ## get KSEA scores integrated with differential phosphorylation results
  p = ggplot(df_cans)
  p = p + geom_point(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin.aa_reported), 
                         color = pass_FDR), stroke = 0, alpha = 0.5)
  p = p + theme_bw() + theme_nogrid()
  p = p + facet_grid(Cancer~., scales = "free", space = "fixed")
  p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
  p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
  p = p + geom_text_repel(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin.aa_reported), 
                              label= ifelse(-log10(FDR_pho_kin.aa_reported) > text_cutoff & coef_pho_kin_filtered > 0, as.character(pair), NA)), 
                          size=2, color = "grey", alpha = 0.5)
  p = p + geom_text_repel(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin.aa_reported), 
                              label= ifelse(-log10(FDR_pho_kin.aa_reported) > text_cutoff & coef_pho_kin_filtered < 0, as.character(pair), NA)),
                          size=2, color = "black", alpha = 1, force = 3)
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(panel.spacing.y=unit(0, "lines"), strip.text.y = element_text(size = 12))
  p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
  p = p + scale_color_manual(values = c("TRUE" = "#377EB8", "FALSE" = "grey", "NA" = "grey")) 
  p = p + theme(legend.position="none")
  p
  subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
  dir.create(subdir2)
  fn = paste(subdir2, "volcano_trans_", "cancers" , "_", enzyme_type, '_substrate_textcutoff', text_cutoff, '.pdf',sep ="")
  ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
}
