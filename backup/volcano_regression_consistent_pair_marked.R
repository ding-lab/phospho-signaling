# Yige Wu @ WashU 2018 Apr
## draw volcano plots 

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
library(ggrepel)
## input enzyme substrate table
ptms_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv")
op_enzyme <- fread(input = paste0(resultD, "diffexp/tables/integrate_enzyme_KSEAordiffexp_sub_FC_plus_regression/omnipath_enzyme.txt"), data.table = F)
color_map <- c("firebrick1", "black")
names(color_map) <- c("TRUE", "FALSE")
color_direction <- c(set1[1], set1[2], "black")
names(color_direction) <- c("up", "down", "NA")


# kinases -----------------------------------------------------------------
for (cancer in c("BRCA", "OV", "CO")) {
  ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$SELF <- ifelse(as.vector(ksea_diffexp$GENE == ksea_diffexp$SUB_GENE), TRUE, FALSE)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE,":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  
  for (enzyme_type in c("kinase")) {
    subdir1 <- paste0(makeOutDir(), enzyme_type, "/")
    dir.create(subdir1)
    df0 <- ksea_diffexp[ksea_diffexp$enzyme_type == enzyme_type & !is.na(ksea_diffexp$FDR_pho_kin) & !ksea_diffexp$SELF,]
    
    df1 <- df0[!is.na(df0$FDR_pho_kin) & -log10(df0$FDR_pho_kin) < 30,]
    df1$coef_pho_kin_filtered = remove_outliers(df1$coef_pho_kin, out_thres = 1.5)
    df1$substrate_direction_filtered <- "NA"
    df1$substrate_direction_filtered[which(df1$consistent)] <- df1$substrate_direction[which(df1$consistent)]
    
    df <- unique(df1[,c("coef_pho_kin_filtered", "FDR_pho_kin", "substrate_direction_filtered", "pair")])
    
    ## write a summary about the overlap btw regulated pair and consistently high/low pairs
    df2 <- df1[df1$FDR_pho_kin < 0.1 & df1$coef_pho_kin > 0,]
    df3 <- df1[df1$FDR_pho_kin < 0.1 & df1$coef_pho_kin > 0.1,]
    df2.5 <- df2[!is.na(df2$consistent),]
    df4 <- df2[df2$substrate_direction_filtered == "up",]
    df5<- df2[df2$substrate_direction_filtered == "down",]
    df6 <- df4[df4$substrate_FC > 2,]
    df7 <- df5[df5$substrate_FC < 0.5,]
    
    tn = paste(subdir1, cancer , "_", enzyme_type, '_substrate.txt',sep ="")
    sink(file = tn)
    cat(paste0("starting from ", length(unique(df0$pair)), " known trans ", enzyme_type, " phosphosite relationships\n"))
    cat(paste0("tested the co-phosphorylation of ", length(unique(df1$pair)), " ", enzyme_type, "-phosphosite pairs\n"))
    cat(paste0(length(unique(df2$pair)), " of above have FDR<0.1, beta>0; ", length(unique(df3$pair)), " of above have FDR<0.1, beta>0.1; ", "\n"))
    cat(paste0(length(unique(df2.5$pair)), " of above(FDR<0.1, beta>0) are tested for consistent fold change in ", enzyme_type, "-phosphosite pairs", "\n"))
    cat(paste0(length(unique(df4$pair)), " of above(FDR<0.1, beta>0) are ", enzyme_type, "-high-substrate-high pairs", "\n"))
    cat(paste0(length(unique(df6$pair)), " of above ", 
               enzyme_type, "-high-substrate-high pairs are significantly differentially phosphorylated in their phosphosite(substrate_log2FC>1)", "\n"))
    print(df6$pair)
    cat(paste0("\n", length(unique(df5$pair)), " of above(FDR<0.1, beta>0) are ", enzyme_type, "-low-substrate-low pairs", "\n"))
    cat(paste0(length(unique(df7$pair)), " of above ", 
               enzyme_type, "-low-substrate-low pairs are significantly differentially phosphorylated in their phosphosite(substrate_log2FC<-1)", "\n"))
    print(df7$pair)
    
    sink()
    closeAllConnections()
    for (text_cutoff in c(4)) {
      ## get KSEA scores integrated with differential phosphorylation results
      p = ggplot(df)
      p = p + geom_point(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin), 
                             color = substrate_direction_filtered, alpha = ifelse(substrate_direction_filtered == "NA", 0.05, 0.5)),
                         stroke = 0)
      p = p + theme_bw() + theme_nogrid()
      p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
      p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
      p = p + geom_text_repel(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin), 
                                  label= ifelse(substrate_direction_filtered != "NA" & -log10(FDR_pho_kin)>text_cutoff, as.character(pair), NA), 
                                  color = substrate_direction_filtered), size=3, alpha=0.5)
      p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
      p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
      p = p + scale_color_manual(values = color_direction) 
      p = p + theme(legend.position="none")
      p
      subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
      dir.create(subdir2)
      fn = paste(subdir2, "volcano_", cancer , "_", enzyme_type, '_substrate_textcutoff', text_cutoff, '.pdf',sep ="")
      ggsave(file=fn, height=6, width=6, useDingbats=FALSE)
      
    }
  }
}


# phosphotase -------------------------------------------------------------


for (cancer in c("BRCA", "OV", "CO")) {
  ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$SELF <- ifelse(as.vector(ksea_diffexp$GENE == ksea_diffexp$SUB_GENE), TRUE, FALSE)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE,":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  
  for (enzyme_type in c("phosphotase")) {
    subdir1 <- paste0(makeOutDir(), enzyme_type, "/")
    dir.create(subdir1)
    df01 <- ksea_diffexp[ksea_diffexp$enzyme_type == enzyme_type & !is.na(ksea_diffexp$KINASE) & !ksea_diffexp$SELF,]
    df0 <- ksea_diffexp[ksea_diffexp$enzyme_type == enzyme_type & !ksea_diffexp$SELF,]
    
    df1 <- df0[!is.na(df0$FDR_pho_kin) & -log10(df0$FDR_pho_kin) < 30,]
    df1$coef_pho_kin_filtered = remove_outliers(df1$coef_pho_kin, out_thres = 1.5)
    df1$substrate_direction_filtered <- "NA"
    df1$substrate_direction_filtered[which(df1$consistent)] <- df1$substrate_direction[which(df1$consistent)]
  
    ## write a summary about the overlap btw regulated pair and consistently high/low pairs
    df2 <- df1[df1$FDR_pho_kin < 0.1 & df1$coef_pho_kin < 0,]
    df3 <- df1[df1$FDR_pho_kin < 0.1 & df1$coef_pho_kin < -0.1,]
    df2.5 <- df2[!is.na(df2$consistent),]
    df4 <- df2[df2$substrate_direction_filtered == "up",]
    df5<- df2[df2$substrate_direction_filtered == "down",]
    df6 <- df4[df4$substrate_FC > 2,]
    df7 <- df5[df5$substrate_FC < 0.5,]
    
    tn = paste(subdir1, cancer , "_", enzyme_type, '_substrate.txt',sep ="")
    sink(file = tn)
    cat(paste0("starting from ", length(unique(df01$pair)), " known trans ", enzyme_type, " phosphosite relationships\n"))
    cat(paste0("expanding to test all phosphosites per substrate: tested ", length(unique(df0$pair)), " trans ", enzyme_type, " phosphosite relationships\n"))
    cat(paste0("tested the co-phosphorylation of ", length(unique(df1$pair)), " ", enzyme_type, "-phosphosite pairs in tumors\n"))
    cat(paste0(length(unique(df2$pair)), " of above have FDR<0.1, beta<0; ", length(unique(df3$pair)), " of above have FDR<0.1, beta<-0.1; ", "\n"))
    print(df2$pair)
    cat("\n")
    cat(paste0(length(unique(df2.5$pair)), " of above(FDR<0.1, beta<0) are tested for consistent fold change in ", enzyme_type, "-phosphosite pairs", "\n"))
    cat(paste0(length(unique(df4$pair)), " of above(FDR<0.1, beta<0) are ", enzyme_type, "-low-substrate-high pairs", "\n"))
    cat(paste0(length(unique(df6$pair)), " of above ", 
               enzyme_type, "-low-substrate-high pairs are significantly differentially phosphorylated in their phosphosite(substrate_log2FC>1)", "\n"))
    print(df6$pair)
    cat(paste0(length(unique(df5$pair)), " of above(FDR<0.1, beta<0) are ", enzyme_type, "-high-substrate-low pairs", "\n"))
    cat(paste0(length(unique(df7$pair)), " of above ", 
               enzyme_type, "-high-substrate-low pairs are significantly differentially phosphorylated in their phosphosite(substrate_log2FC<-1)", "\n"))
    print(df7$pair)
    
    df2 <- df1[df1$FDR_pho_kin.tn < 0.1 & df1$coef_pho_kin.tn < 0,]
    df3 <- df1[df1$FDR_pho_kin.tn < 0.1 & df1$coef_pho_kin.tn < -0.1,]
    df2.5 <- df2[!is.na(df2$consistent),]
    df4 <- df2[df2$substrate_direction_filtered == "up",]
    df5<- df2[df2$substrate_direction_filtered == "down",]
    df6 <- df4[df4$substrate_FC > 2,]
    df7 <- df5[df5$substrate_FC < 0.5,]
    cat(paste0("\nwhen tested the co-phosphorylation in tumors+normals: \n", length(unique(df2$pair)), " of above have FDR<0.1, beta<0; ", length(unique(df3$pair)), " of above have FDR<0.1, beta<-0.1; ", "\n"))
    print(df2$pair)
    cat(paste0(length(unique(df2.5$pair)), " of above(FDR<0.1, beta<0) are tested for consistent fold change in ", enzyme_type, "-phosphosite pairs", "\n"))
    cat(paste0(length(unique(df4$pair)), " of above(FDR<0.1, beta<0) are ", enzyme_type, "-low-substrate-high pairs", "\n"))
    cat(paste0(length(unique(df6$pair)), " of above ", 
               enzyme_type, "-low-substrate-high pairs are significantly differentially phosphorylated in their phosphosite(substrate_log2FC>1)", "\n"))
    print(df6$pair)
    cat("\n")
    cat(paste0(length(unique(df5$pair)), " of above(FDR<0.1, beta<0) are ", enzyme_type, "-high-substrate-low pairs", "\n"))
    cat(paste0(length(unique(df7$pair)), " of above ", 
               enzyme_type, "-high-substrate-low pairs are significantly differentially phosphorylated in their phosphosite(substrate_log2FC<-1)", "\n"))
    print(df7$pair)
    sink()
    closeAllConnections()
    for (text_cutoff in c(1)) {
      df <- unique(df1[,c("coef_pho_kin", "FDR_pho_kin", "coef_pho_kin.tn", "FDR_pho_kin.tn", "substrate_direction_filtered", "pair")])
      p = ggplot(df[!is.na(remove_outliers(df$coef_pho_kin, out_thres = 1.5)) & -log10(df$FDR_pho_kin.tn) < 10,])
      p = p + geom_point(aes(x=coef_pho_kin, y=-log10(FDR_pho_kin), 
                             color = substrate_direction_filtered, alpha = ifelse(substrate_direction_filtered == "NA", 0.05, 0.5)),
                         stroke = 0)
      p = p + theme_bw() + theme_nogrid()
      p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
      p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
      p = p + geom_text_repel(aes(x=coef_pho_kin, y=-log10(FDR_pho_kin), 
                                  label= ifelse(-log10(FDR_pho_kin)>text_cutoff & coef_pho_kin < 0, as.character(pair), NA), 
                                  color = substrate_direction_filtered,
                                  alpha = 1), size=4)
      p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
      p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
      p = p + scale_color_manual(values = color_direction) 
      p = p + theme(legend.position="none")
      p
      subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
      dir.create(subdir2)
      fn = paste(subdir2, "volcano_", cancer , "_", enzyme_type, '_substrate_textcutoff', text_cutoff, '_expanded.pdf',sep ="")
      ggsave(file=fn, height=6, width=6, useDingbats=FALSE)
      
      p = ggplot(df[!is.na(remove_outliers(df$coef_pho_kin.tn, out_thres = 1.5)) & -log10(df$FDR_pho_kin.tn) < 10,])
      p = p + geom_point(aes(x=coef_pho_kin.tn, y=-log10(FDR_pho_kin.tn), 
                             color = substrate_direction_filtered, alpha = 0.9),
                         stroke = 0)
      p = p + theme_bw() + theme_nogrid()
      p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
      p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
      p = p + geom_text_repel(aes(x=coef_pho_kin.tn, y=-log10(FDR_pho_kin.tn), 
                                  label= ifelse((-log10(FDR_pho_kin.tn) > text_cutoff & coef_pho_kin.tn < 0) | (substrate_direction_filtered != "NA"), as.character(pair), NA), 
                                  color = substrate_direction_filtered,
                                  alpha = 1), size=4)
      p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
      p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
      p = p + scale_color_manual(values = color_direction) 
      p = p + theme(legend.position="none")
      p
      subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
      dir.create(subdir2)
      fn = paste(subdir2, "volcano_", cancer , "_", enzyme_type, '_substrate_textcutoff', text_cutoff, '_expanded_tumors&normals.pdf',sep ="")
      ggsave(file=fn, height=6, width=6, useDingbats=FALSE)
    }
  }
}


