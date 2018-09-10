# Yige Wu @ WashU 2018 Apr
## check regulated pairs with mutational impact from kinases

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

library(ggrepel)
devtools::install_github("slowkow/ggrepel")
library(ggplot2)

# inputs ------------------------------------------------------------------
## input druggable gene list
# drug_genes <- fread()
color_direction <- c(set1[1], set1[2], "black")
names(color_direction) <- c("up", "down", "NA")
clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)
enzyme_type <- 'kinase'
ksea_diffexp <- fread(input = paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                                     '3can', "_KSEA_enzyme_diffexp_substrate_aa_reported.txt"), data.table = F)
ksea_diffexp$pair <- paste0(ksea_diffexp$GENE, ":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
ksea_diffexp <- markSigSiteCan(ksea_diffexp, sig_thres = 0.1, enzyme_type = enzyme_type)
cancers_sort <- c("BRCA", "OV", "CO")

# trans-regulated pairs overlap with driver mutations --------------
# for (cancer in c("BRCA", "OV", "CO")) {
for (cancer in c("OV")) { 
# for (cancer in c("CO")) {
  ## input super table
  for (enzyme_type in c("kinase")) {
    tmp <- ksea_diffexp
    tmp <- tmp[tmp$fdr_sig & tmp$coef_sig & tmp$SELF == "trans" & tmp$Cancer == cancer,]
    tmp <- tmp[!is.na(tmp[, paste0('uniq_', cancer)]) & tmp[, paste0('uniq_', cancer)] ,]
    tmp <- tmp[ (tmp$enzyme.driver_type != "notdriver" | tmp$substrate.driver_type != "notdriver"),]
    # tmp <- ksea_diffexp[ksea_diffexp$FDR_pho_kin < 0.2 & !is.na(ksea_diffexp$FDR_pho_kin) & (ksea_diffexp$enzyme.driver_type != "notdriver" | ksea_diffexp$substrate.driver_type != "notdriver"),]
    # tmp <- tmp[(!is.na(tmp$p_value) & tmp$p_value < 0.2) | (!is.na(tmp$p_value.pro) & tmp$p_value.pro < 0.2),]
    
    if (nrow(tmp) > 0) {
      ## create output directory
      outDir1 <- paste0(makeOutDir(resultD = resultD), cancer, "/")
      dir.create(outDir1)
      #for (i in 1:nrow(tmp)) {
      for (i in c(1)) {
        # enzyme <- as.character(tmp[i, "GENE"])
        # substrate <- as.character(tmp[i, "SUB_GENE"])
        enzyme <- 'IKBKB'
        substrate <- 'TP53'
        # for (gene in c("MAP2K4")) {
        # for (rsd in c("S268")) {
        for (rsd in unique(tmp$SUB_MOD_RSD[tmp$SUB_GENE == substrate])) {
          sup_3can <- NULL
          for (cancer2 in c("BRCA", "OV", "CO")) {
            ## input protein level
            pro_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PRO_formatted_normalized_noControl.txt",sep=""), data.table = F)
            pro_sub <- pro_data[pro_data$Gene == substrate,]
            pro_kin <- pro_data[pro_data$Gene == enzyme,]
            rm(pro_data)
            ## input phospho level
            pho_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
            pho_data <- pho_data[pho_data$Gene == substrate,]
            
            ## input phospho level
            phog_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
            phog_data <- phog_data[phog_data$Gene == enzyme,]
            
            if (nrow(pro_sub) > 0 & nrow(pho_data) > 0 & nrow(phog_data)) {
              pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
              pho_data <- pho_data[pho_head$SUB_MOD_RSD == rsd,]
              
              if ( nrow(pho_data) > 0) {
                pro.m <- melt(pro_sub)
                colnames(pro.m)[ncol(pro.m)] <- "pro_sub"
                
                pro.m2 <- melt(pro_kin)
                colnames(pro.m2)[ncol(pro.m2)] <- "pro_kin"
                
                pho.m <- melt(pho_data)
                colnames(pho.m)[ncol(pho.m)] <- "pho_sub"
                
                phog.m <- melt(phog_data)
                colnames(phog.m)[ncol(phog.m)] <- "pho_kin"
                
                sup_tab <- merge(pro.m[, c("variable", "pro_sub")], pho.m[, c("variable", "pho_sub")], all = T)
                sup_tab <- merge(sup_tab, phog.m[, c("variable", "pho_kin")], all = T)
                sup_tab <- merge(sup_tab, pro.m2[, c("variable", "pro_kin")], all = T)
                sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
                
                ## input maf
                maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/CPTAC2_", cancer2, "_prospective.v1.3.somatic.variants.031918.maf"), data.table = F)
                maf <- maf[maf$Hugo_Symbol == enzyme & maf$Variant_Classification != "Silent",]
                if (nrow(maf) > 0) {
                  maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
                  maf$aa_change <- paste0(maf$Hugo_Symbol, ":", maf$HGVSp_Short)
                  maf$is.upstream <- ifelse(maf$Hugo_Symbol == enzyme, TRUE, FALSE)
                  sup_tab <- merge(sup_tab, maf[, c("partID", "Variant_Classification", "aa_change", "is.upstream")], all.x = T) 
                } else {
                  sup_tab$aa_change <- NA
                  sup_tab$is.upstream <- NA
                  sup_tab$Variant_Classification <- NA
                }
                sup_tab$cancer <- cancer2
                sup_3can <- rbind(sup_3can, sup_tab)
              }
            }
          }
          ## make deeper directory for outpus
          outDir2 <- paste0(outDir1, enzyme, '_', substrate, '_', rsd, "/")
          dir.create(outDir2)
          
          ##  add label
          sup_3can$cancer <- factor(sup_3can$cancer, levels = cancers_sort)
          if (length(which(is.na(sup_3can$is.upstream))) > 0) {
            sup_3can$is.upstream[is.na(sup_3can$is.upstream)] <- "NA"
          }  
          
          ## scatterplot
          p = ggplot(sup_3can, aes(x=pho_kin, y=pho_sub))
          p = p + geom_point(aes(color = Variant_Classification, 
                                 alpha = ifelse(is.na(Variant_Classification), 0.8, 1)), stroke = 0)
          p <- p + facet_grid(.~cancer, scales = "free_y", space = "fixed")
          p = p + geom_text_repel(aes(pho_kin, pho_sub, label= ifelse(Variant_Classification != "Silent", as.character(aa_change), NA), color = Variant_Classification, segment.color = Variant_Classification),
                                  force = 1, 
                                  segment.size = 0.5, segment.alpha = 0.2, 
                                  size=1.5,alpha=0.8)
          p = p + geom_smooth(mapping = aes(x = pho_kin, y = pho_sub), method = "glm", se=FALSE, color="grey", formula = y ~ x, linetype = 2, size = 0.5)
          p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
          p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
          p = p + labs(x = paste0(enzyme, " phosphorylation abundance(log2 ratio)"), 
                       y = paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio"))
          p = p + theme_nogrid()
          p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                        axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                        axis.text.y = element_text(colour="black", size=6))#element_text(colour="black", size=14))
          p = p + theme(title = element_text(size = 8))
          p
          fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_sub+pho_kin~pho_sub.pdf")
          ggsave(file=fn, height=3, width=6, useDingbats=FALSE)
          
          ## boxplot for substrate phosphorylation
          p = ggplot(sup_3can, aes(x=cancer, y=pho_sub, color = is.upstream, label= as.character(aa_change)))
          p = p + geom_point(aes(shape = ifelse(!is.na(is.upstream) & (is.upstream == TRUE), "b", "a")), 
                             position = position_jitterdodge(), stroke = 0, alpha = 0.5)
          p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
          p = p + geom_text_repel(aes(segment.color = Variant_Classification),
                                  force = 1,
                                  segment.size = 0.5, segment.alpha = 0.2,
                                  size=1.5,alpha=0.8)
          p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio"))
          p = p + theme_nogrid()
          p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                        axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                        axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
          p = p + theme(title = element_text(size = 8))
          p
          fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pho_sub.pdf")
          ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
          
          ## boxplot for substrate protein
          p = ggplot(sup_3can, aes(x=cancer, y=pro_sub, color = is.upstream, label= as.character(aa_change)))
          p = p + geom_point(aes(shape = ifelse(!is.na(is.upstream) & (is.upstream == TRUE), "b", "a")), 
                             position = position_jitterdodge(), stroke = 0, alpha = 0.5)
          p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
          p = p + geom_text_repel(aes(segment.color = Variant_Classification),
                                  force = 1,
                                  segment.size = 0.5, segment.alpha = 0.2,
                                  size=1.5,alpha=0.8)
          p = p + labs(y = paste0(substrate, " protein abundance(log2 ratio"))
          p = p + theme_nogrid()
          p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                        axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                        axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
          p = p + theme(title = element_text(size = 8))
          p
          fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_sub.pdf")
          ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
          
          ## boxplot for kinase protein
          p = ggplot(sup_3can, aes(x=cancer, y=pro_kin, color = is.upstream, label= as.character(aa_change)))
          p = p + geom_point(aes(shape = ifelse(!is.na(is.upstream) & (is.upstream == TRUE), "b", "a")), 
                             position = position_jitterdodge(), stroke = 0, alpha = 0.5)
          p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
          p = p + geom_text_repel(aes(segment.color = Variant_Classification),
                                  force = 1,
                                  segment.size = 0.5, segment.alpha = 0.2,
                                  size=1.5,alpha=0.8)
          p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio"))
          p = p + theme_nogrid()
          p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                        axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                        axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
          p = p + theme(title = element_text(size = 8))
          p
          fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_kin.pdf")
          ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
        }
      }
    }
  }
}

for (cancer in c("CO")) {
  ## input super table
  for (enzyme_type in c("kinase")) {
    tmp <- ksea_diffexp
    tmp <- tmp[tmp$fdr_sig & tmp$coef_sig & tmp$SELF == "trans" & tmp$Cancer == cancer,]
    tmp <- tmp[!is.na(tmp[, paste0('uniq_', cancer)]) & tmp[, paste0('uniq_', cancer)] ,]
    tmp <- tmp[ (tmp$enzyme.driver_type != "notdriver" | tmp$substrate.driver_type != "notdriver"),]
    # tmp <- ksea_diffexp[ksea_diffexp$FDR_pho_kin < 0.2 & !is.na(ksea_diffexp$FDR_pho_kin) & (ksea_diffexp$enzyme.driver_type != "notdriver" | ksea_diffexp$substrate.driver_type != "notdriver"),]
    # tmp <- tmp[(!is.na(tmp$p_value) & tmp$p_value < 0.2) | (!is.na(tmp$p_value.pro) & tmp$p_value.pro < 0.2),]
    
    if (nrow(tmp) > 0) {
      ## create output directory
      outDir1 <- paste0(makeOutDir(resultD = resultD), cancer, "/")
      dir.create(outDir1)
      for (i in 1:nrow(tmp)) {
      # for (i in c(1)) {
        # enzyme <- as.character(tmp[i, "GENE"])
        # substrate <- as.character(tmp[i, "SUB_GENE"])
        enzyme <- 'IKBKB'
        substrate <- 'TP53'
        # for (gene in c("MAP2K4")) {
        # for (rsd in c("S268")) {
        for (rsd in unique(tmp$SUB_MOD_RSD[tmp$SUB_GENE == substrate])) {
          sup_3can <- NULL
          for (cancer2 in c("BRCA", "OV", "CO")) {
            ## input protein level
            pro_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PRO_formatted_normalized_noControl.txt",sep=""), data.table = F)
            pro_sub <- pro_data[pro_data$Gene == substrate,]
            pro_kin <- pro_data[pro_data$Gene == enzyme,]
            rm(pro_data)
            ## input phospho level
            pho_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
            pho_data <- pho_data[pho_data$Gene == substrate,]
            
            ## input phospho level
            phog_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
            phog_data <- phog_data[phog_data$Gene == enzyme,]
            
            if (nrow(pro_sub) > 0 & nrow(pho_data) > 0 & nrow(phog_data)) {
              pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
              pho_data <- pho_data[pho_head$SUB_MOD_RSD == rsd,]
              
              if ( nrow(pho_data) > 0) {
                pro.m <- melt(pro_sub)
                colnames(pro.m)[ncol(pro.m)] <- "pro_sub"
                
                pro.m2 <- melt(pro_kin)
                colnames(pro.m2)[ncol(pro.m2)] <- "pro_kin"
                
                pho.m <- melt(pho_data)
                colnames(pho.m)[ncol(pho.m)] <- "pho_sub"
                
                phog.m <- melt(phog_data)
                colnames(phog.m)[ncol(phog.m)] <- "pho_kin"
                
                sup_tab <- merge(pro.m[, c("variable", "pro_sub")], pho.m[, c("variable", "pho_sub")], all = T)
                sup_tab <- merge(sup_tab, phog.m[, c("variable", "pho_kin")], all = T)
                sup_tab <- merge(sup_tab, pro.m2[, c("variable", "pro_kin")], all = T)
                sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
                
                ## input maf
                maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/CPTAC2_", cancer2, "_prospective.v1.3.somatic.variants.031918.maf"), data.table = F)
                maf <- maf[maf$Hugo_Symbol == enzyme & maf$Variant_Classification != "Silent",]
                if (nrow(maf) > 0) {
                  maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
                  maf$aa_change <- paste0(maf$Hugo_Symbol, ":", maf$HGVSp_Short)
                  maf$is.upstream <- ifelse(maf$Hugo_Symbol == enzyme, TRUE, FALSE)
                  sup_tab <- merge(sup_tab, maf[, c("partID", "Variant_Classification", "aa_change", "is.upstream")], all.x = T) 
                } else {
                  sup_tab$aa_change <- NA
                  sup_tab$is.upstream <- NA
                  sup_tab$Variant_Classification <- NA
                }
                sup_tab$cancer <- cancer2
                sup_3can <- rbind(sup_3can, sup_tab)
              }
            }
          }
          ## make deeper directory for outpus
          outDir2 <- paste0(outDir1, enzyme, '_', substrate, '_', rsd, "/")
          dir.create(outDir2)
          
          ##  add label
          sup_3can$cancer <- factor(sup_3can$cancer, levels = cancers_sort)
          if (length(which(is.na(sup_3can$is.upstream))) > 0) {
            sup_3can$is.upstream[is.na(sup_3can$is.upstream)] <- "NA"
          }  
          
          ## scatterplot
          p = ggplot(sup_3can, aes(x=pho_kin, y=pho_sub))
          p = p + geom_point(aes(color = Variant_Classification, 
                                 alpha = ifelse(is.na(Variant_Classification), 0.8, 1)), stroke = 0)
          p <- p + facet_grid(.~cancer, scales = "free_y", space = "fixed")
          p = p + geom_text_repel(aes(pho_kin, pho_sub, label= ifelse(Variant_Classification != "Silent", as.character(aa_change), NA), color = Variant_Classification, segment.color = Variant_Classification),
                                  force = 1, 
                                  segment.size = 0.5, segment.alpha = 0.2, 
                                  size=1.5,alpha=0.8)
          p = p + geom_smooth(mapping = aes(x = pho_kin, y = pho_sub), method = "glm", se=FALSE, color="grey", formula = y ~ x, linetype = 2, size = 0.5)
          p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
          p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
          p = p + labs(x = paste0(enzyme, " phosphorylation abundance(log2 ratio)"), 
                       y = paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio"))
          p = p + theme_nogrid()
          p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                        axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                        axis.text.y = element_text(colour="black", size=6))#element_text(colour="black", size=14))
          p = p + theme(title = element_text(size = 8))
          p
          fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_sub+pho_kin~pho_sub.pdf")
          ggsave(file=fn, height=3, width=6, useDingbats=FALSE)
          
          ## boxplot for substrate phosphorylation
          p = ggplot(sup_3can, aes(x=cancer, y=pho_sub, color = is.upstream, label= as.character(aa_change)))
          p = p + geom_point(aes(shape = ifelse(!is.na(is.upstream) & (is.upstream == TRUE), "b", "a")), 
                             position = position_jitterdodge(), stroke = 0, alpha = 0.5)
          p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
          p = p + geom_text_repel(aes(segment.color = Variant_Classification),
                                  force = 1,
                                  segment.size = 0.5, segment.alpha = 0.2,
                                  size=1.5,alpha=0.8)
          p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio"))
          p = p + theme_nogrid()
          p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                        axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                        axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
          p = p + theme(title = element_text(size = 8))
          p
          fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pho_sub.pdf")
          ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
          
          ## boxplot for substrate protein
          p = ggplot(sup_3can, aes(x=cancer, y=pro_sub, color = is.upstream, label= as.character(aa_change)))
          p = p + geom_point(aes(shape = ifelse(!is.na(is.upstream) & (is.upstream == TRUE), "b", "a")), 
                             position = position_jitterdodge(), stroke = 0, alpha = 0.5)
          p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
          p = p + geom_text_repel(aes(segment.color = Variant_Classification),
                                  force = 1,
                                  segment.size = 0.5, segment.alpha = 0.2,
                                  size=1.5,alpha=0.8)
          p = p + labs(y = paste0(substrate, " protein abundance(log2 ratio"))
          p = p + theme_nogrid()
          p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                        axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                        axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
          p = p + theme(title = element_text(size = 8))
          p
          fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_sub.pdf")
          ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
          
          ## boxplot for kinase protein
          p = ggplot(sup_3can, aes(x=cancer, y=pro_kin, color = is.upstream, label= as.character(aa_change)))
          p = p + geom_point(aes(shape = ifelse(!is.na(is.upstream) & (is.upstream == TRUE), "b", "a")), 
                             position = position_jitterdodge(), stroke = 0, alpha = 0.5)
          p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
          p = p + geom_text_repel(aes(segment.color = Variant_Classification),
                                  force = 1,
                                  segment.size = 0.5, segment.alpha = 0.2,
                                  size=1.5,alpha=0.8)
          p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio"))
          p = p + theme_nogrid()
          p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                        axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                        axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
          p = p + theme(title = element_text(size = 8))
          p
          fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_kin.pdf")
          ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
        }
      }
    }
  }
}


# plot shared -------------------------------------------------------------
for (enzyme_type in c("kinase")) {
  tmp <- ksea_diffexp
  tmp <- tmp[tmp$fdr_sig & tmp$coef_sig & tmp$SELF == "trans",]
  tmp <- tmp[!is.na(tmp$shared3can) & tmp$shared3can,]
  # tmp <- ksea_diffexp[ksea_diffexp$FDR_pho_kin < 0.2 & !is.na(ksea_diffexp$FDR_pho_kin) & (ksea_diffexp$enzyme.driver_type != "notdriver" | ksea_diffexp$substrate.driver_type != "notdriver"),]
  # tmp <- tmp[(!is.na(tmp$p_value) & tmp$p_value < 0.2) | (!is.na(tmp$p_value.pro) & tmp$p_value.pro < 0.2),]
  
  if (nrow(tmp) > 0) {
    ## create output directory
    outDir1 <- paste0(makeOutDir(resultD = resultD), "shared3can/")
    dir.create(outDir1)
    for (i in 1:nrow(tmp)) {
      # for (i in c(1)) {
      enzyme <- as.character(tmp[i, "GENE"])
      substrate <- as.character(tmp[i, "SUB_GENE"])

      # for (gene in c("MAP2K4")) {
      # for (rsd in c("S268")) {
      for (rsd in unique(tmp$SUB_MOD_RSD[tmp$SUB_GENE == substrate])) {
        sup_3can <- NULL
        for (cancer2 in c("BRCA", "OV", "CO")) {
          ## input protein level
          pro_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PRO_formatted_normalized_noControl.txt",sep=""), data.table = F)
          pro_sub <- pro_data[pro_data$Gene == substrate,]
          pro_kin <- pro_data[pro_data$Gene == enzyme,]
          rm(pro_data)
          ## input phospho level
          pho_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
          pho_data <- pho_data[pho_data$Gene == substrate,]
          
          ## input phospho level
          phog_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
          phog_data <- phog_data[phog_data$Gene == enzyme,]
          
          if (nrow(pro_sub) > 0 & nrow(pho_data) > 0 & nrow(phog_data)) {
            pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
            pho_data <- pho_data[pho_head$SUB_MOD_RSD == rsd,]
            
            if ( nrow(pho_data) > 0) {
              pro.m <- melt(pro_sub)
              colnames(pro.m)[ncol(pro.m)] <- "pro_sub"
              
              pro.m2 <- melt(pro_kin)
              colnames(pro.m2)[ncol(pro.m2)] <- "pro_kin"
              
              pho.m <- melt(pho_data)
              colnames(pho.m)[ncol(pho.m)] <- "pho_sub"
              
              phog.m <- melt(phog_data)
              colnames(phog.m)[ncol(phog.m)] <- "pho_kin"
              
              sup_tab <- merge(pro.m[, c("variable", "pro_sub")], pho.m[, c("variable", "pho_sub")], all = T)
              sup_tab <- merge(sup_tab, phog.m[, c("variable", "pho_kin")], all = T)
              sup_tab <- merge(sup_tab, pro.m2[, c("variable", "pro_kin")], all = T)
              sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
              
              ## input maf
              maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/CPTAC2_", cancer2, "_prospective.v1.3.somatic.variants.031918.maf"), data.table = F)
              maf <- maf[maf$Hugo_Symbol == enzyme & maf$Variant_Classification != "Silent",]
              if (nrow(maf) > 0) {
                maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
                maf$aa_change <- paste0(maf$Hugo_Symbol, ":", maf$HGVSp_Short)
                maf$is.upstream <- ifelse(maf$Hugo_Symbol == enzyme, TRUE, FALSE)
                sup_tab <- merge(sup_tab, maf[, c("partID", "Variant_Classification", "aa_change", "is.upstream")], all.x = T) 
              } else {
                sup_tab$aa_change <- NA
                sup_tab$is.upstream <- NA
                sup_tab$Variant_Classification <- NA
              }
              sup_tab$cancer <- cancer2
              sup_3can <- rbind(sup_3can, sup_tab)
            }
          }
        }
        ## make deeper directory for outpus
        outDir2 <- paste0(outDir1, enzyme, '_', substrate, '_', rsd, "/")
        dir.create(outDir2)
        
        ##  add label
        sup_3can$cancer <- factor(sup_3can$cancer, levels = cancers_sort)
        if (length(which(is.na(sup_3can$is.upstream))) > 0) {
          sup_3can$is.upstream[is.na(sup_3can$is.upstream)] <- "NA"
        }  
        
        ## scatterplot
        p = ggplot(sup_3can, aes(x=pho_kin, y=pho_sub))
        if (length(which(!is.na(sup_3can$Variant_Classification))) > 0) {
          p = p + geom_point(aes(color = Variant_Classification, 
                                 alpha = ifelse(is.na(Variant_Classification), 0.8, 1)), stroke = 0)
          p = p + geom_text_repel(aes(pho_kin, pho_sub, label= ifelse(Variant_Classification != "Silent", as.character(aa_change), NA), color = Variant_Classification, segment.color = Variant_Classification),
                                  force = 1, 
                                  segment.size = 0.5, segment.alpha = 0.2, 
                                  size=1.5,alpha=0.8)
        } else {
          p = p + geom_point(stroke = 0, alpha = 0.8, color = 'grey')
        }
        p = p + geom_smooth(mapping = aes(x = pho_kin, y = pho_sub), method = "glm", se=FALSE, color="grey", formula = y ~ x, linetype = 2, size = 0.5)
        p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
        p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
        p <- p + facet_grid(.~cancer, scales = "free_y", space = "fixed")
        p = p + labs(x = paste0(enzyme, " phosphorylation abundance(log2 ratio)"), 
                     y = paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio"))
        p = p + theme_nogrid()
        p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                      axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                      axis.text.y = element_text(colour="black", size=6))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p
        fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_sub+pho_kin~pho_sub.pdf")
        ggsave(file=fn, height=3, width=6, useDingbats=FALSE)
        
        ## boxplot for substrate phosphorylation
        p = ggplot(sup_3can, aes(x=cancer, y=pho_sub, color = is.upstream, label= as.character(aa_change)))
        p = p + geom_point(aes(shape = ifelse(!is.na(is.upstream) & (is.upstream == TRUE), "b", "a")), 
                           position = position_jitterdodge(), stroke = 0, alpha = 0.5)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        if (length(which(!is.na(sup_3can$Variant_Classification))) > 0) {
         p = p + geom_text_repel(aes(segment.color = Variant_Classification),
                                  force = 1,
                                  segment.size = 0.5, segment.alpha = 0.2,
                                  size=1.5,alpha=0.8)
        }
        p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio"))
        p = p + theme_nogrid()
        p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                      axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p
        fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pho_sub.pdf")
        ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
        
        ## boxplot for substrate protein
        p = ggplot(sup_3can, aes(x=cancer, y=pro_sub, color = is.upstream, label= as.character(aa_change)))
        p = p + geom_point(aes(shape = ifelse(!is.na(is.upstream) & (is.upstream == TRUE), "b", "a")), 
                           position = position_jitterdodge(), stroke = 0, alpha = 0.5)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        p = p + geom_text_repel(aes(segment.color = Variant_Classification),
                                force = 1,
                                segment.size = 0.5, segment.alpha = 0.2,
                                size=1.5,alpha=0.8)
        p = p + labs(y = paste0(substrate, " protein abundance(log2 ratio"))
        p = p + theme_nogrid()
        p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                      axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p
        fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_sub.pdf")
        ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
        
        ## boxplot for kinase protein
        p = ggplot(sup_3can, aes(x=cancer, y=pro_kin, color = is.upstream, label= as.character(aa_change)))
        p = p + geom_point(aes(shape = ifelse(!is.na(is.upstream) & (is.upstream == TRUE), "b", "a")), 
                           position = position_jitterdodge(), stroke = 0, alpha = 0.5)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        p = p + geom_text_repel(aes(segment.color = Variant_Classification),
                                force = 1,
                                segment.size = 0.5, segment.alpha = 0.2,
                                size=1.5,alpha=0.8)
        p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio"))
        p = p + theme_nogrid()
        p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                      axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p
        fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_kin.pdf")
        ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
      }
    }
  }
}


