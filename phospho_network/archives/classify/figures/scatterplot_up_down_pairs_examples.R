# Yige Wu @ WashU 2018 Jul
# plot example of up and down regulated kinase substrate pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# variables ---------------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
reg_sig <- 0.05
diff_sig <- 0.2
diff_log2fc <- 1
color_tn <- c("#FB9A99", "#B2DF8A")
names(color_tn) <- c("tumor", "normal")

# inputs -------------------------------------------------------------------
sup_cans_tab <- NULL
for (cancer in cancers_sort) {
  sup_tab <- fread(paste0(ppnD, "kinase_activity/tables/fisher_es_pairs/", 
                          cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  sup_tab$Cancer <- cancer
  sup_cans_tab <- rbind(sup_cans_tab, sup_tab)
}
clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)

# summarize kinase-substrate pairs ----------------------------------------
enzyme_type <- "kinase"
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == enzyme_type,]
sup_cans_tab_pk <- data.frame(sup_cans_tab_pk)
sup_cans_tab_pk <- markSigSiteCan(sup_cans_tab_pk, sig_thres = reg_sig, enzyme_type = enzyme_type)

## annotate kinase substrate regulation
sup_cans_tab_pk$regulated <- (sup_cans_tab_pk$coef_sig & sup_cans_tab_pk$fdr_sig)

## annotate up and down regulation
sup_cans_tab_pk$diffexp_type <- "other"
diffexp_sig <- ((!is.na(sup_cans_tab_pk$KSEA_pvalue) & sup_cans_tab_pk$KSEA_pvalue < diff_sig) | (!is.na(sup_cans_tab_pk$diffexp_log2FC & abs(sup_cans_tab_pk$diffexp_log2FC) > diff_log2fc)))
sup_cans_tab_pk$diffexp_type[diffexp_sig & sup_cans_tab_pk$enzyme_direction == "up" & sup_cans_tab_pk$substrate_direction == "up"] <- "up"
sup_cans_tab_pk$diffexp_type[diffexp_sig & sup_cans_tab_pk$enzyme_direction == "down" & sup_cans_tab_pk$substrate_direction == "down"] <- "down"


# up regulated kinase-substrate pairs -------------------------------------
tmp1 <- sup_cans_tab_pk
tmp1 <- tmp1[tmp1$diffexp_type == "up",]
tmp1 <- tmp1[tmp1$known_enzyme_phosphosite & !is.na(tmp1$known_enzyme_phosphosite) & !is.na(tmp1$Source) & tmp1$Source != "NetKIN",]
tmp1 <- tmp1[order(tmp1$kin_up_sub_up, decreasing = T),]
tmp1 <- tmp1[!is.na(tmp1$sig_BRCA) & !is.na(tmp1$sig_OV) & !is.na(tmp1$sig_CO),]
tmp1 <- tmp1[tmp1$kin_up_sub_up > 0,]

for (enzyme_type in c("kinase")) {
  tmp <- tmp1

  if (nrow(tmp) > 0) {
    ## create output directory
    outDir1 <- paste0(makeOutDir(resultD = resultD), "up/")
    dir.create(outDir1)
    for (enzyme in unique(tmp$GENE)) {
      for (substrate in unique(tmp$SUB_GENE[tmp$GENE == enzyme])[1]) {
      # for (substrate in c("MCM4")) {
        for (rsd in unique(tmp$SUB_MOD_RSD[tmp$SUB_GENE == substrate])) {
          sup_3can <- NULL
          for (cancer2 in c("BRCA", "OV", "CO")) {
            ## input susbstrate phospho level
            pho_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
            pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
            pho_t <- pho_data[pho_data$Gene == substrate & pho_head$SUB_MOD_RSD == rsd,]
            
            pho_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PHO_formatted_normalized_Control.txt",sep=""), data.table = F)
            pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
            pho_n <- pho_data[pho_data$Gene == substrate & pho_head$SUB_MOD_RSD == rsd,]
            
            rm(pho_data)
            
            ## input kinase phospho level
            phog_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
            phog_t <- phog_data[phog_data$Gene == enzyme,]
            
            phog_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_collapsed_PHO_formatted_normalized_replicate_averaged_Normal.txt",sep=""), data.table = F)
            phog_n <- phog_data[phog_data$Gene == enzyme,]
            
            rm(phog_data)
            
            if (nrow(pho_t) > 0 & nrow(pho_n) > 0 & nrow(phog_t) & nrow(phog_n)) {
              pho_t.m <- melt(pho_t)
              colnames(pho_t.m)[ncol(pho_t.m)] <- "pho_sub"
              pho_t.m$tumor_normal <- "tumor"
              
              pho_n.m <- melt(pho_n)
              colnames(pho_n.m)[ncol(pho_n.m)] <- "pho_sub"
              pho_n.m$tumor_normal <- "normal"
              
              phog_t.m <- melt(phog_t)
              colnames(phog_t.m)[ncol(phog_t.m)] <- "phog_kin"
              phog_t.m$tumor_normal <- "tumor"
              
              phog_n.m <- melt(phog_n)
              colnames(phog_n.m)[ncol(phog_n.m)] <- "phog_kin"
              phog_n.m$tumor_normal <- "normal"
              
              tab_t <- merge(pho_t.m[, c("variable", "pho_sub", "tumor_normal")], phog_t.m[, c("variable", "phog_kin", "tumor_normal")], all = T)
              tab_n <- merge(pho_n.m[, c("variable", "pho_sub", "tumor_normal")], phog_n.m[, c("variable", "phog_kin", "tumor_normal")], all = T)
              sup_tab <- rbind(tab_t, tab_n)
              sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
              
              sup_tab$cancer <- cancer2
              sup_3can <- rbind(sup_3can, sup_tab)
            }
          }
          ## make deeper directory for outpus
          outDir2 <- paste0(outDir1, enzyme, '_', substrate, '_', rsd, "/")
          dir.create(outDir2)
          
          ##  add label
          sup_3can$cancer <- factor(sup_3can$cancer, levels = cancers_sort)
          
          ## scatterplot
          p = ggplot(sup_3can, aes(x=phog_kin, y=pho_sub))
          p = p + geom_point(aes(color = tumor_normal), stroke = 0, alpha = 0.8, shape = 16)
          p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
          p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
          p <- p + facet_grid(.~cancer, scales = "free_y", space = "fixed")
          p = p + scale_color_manual(values = color_tn) 
          p = p + labs(x = paste0(enzyme, " phosphorylation abundance(log2 ratio)"), 
                       y = paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio"))
          p = p + theme_nogrid()
          p = p + theme(axis.title = element_text(size=6),
                        axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                        axis.text.y = element_text(colour="black", size=6))
          p = p + theme(title = element_text(size = 8),
                        legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                        legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                        legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                        legend.title = element_text(size = 5), legend.justification = c(0, 1), legend.position = c(0, 1))
          p
          fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_sub+pho_kin~pho_sub.pdf")
          ggsave(file=fn, height=2.5, width=6, useDingbats=FALSE)
        }
      }
    }
  }
}

# down regulated kinase-substrate pairs -------------------------------------
tmp2 <- sup_cans_tab_pk
tmp2 <- tmp2[tmp2$diffexp_type == "down",]
tmp2 <- tmp2[tmp2$known_enzyme_phosphosite & !is.na(tmp2$known_enzyme_phosphosite) & !is.na(tmp2$Source) & tmp2$Source != "NetKIN",]
## filter by data available in 3 cancers
pair_can_count <- data.frame(table(unique(tmp2[, c("pair", "Cancer")])[, c("pair")]))
pair_can_count <- pair_can_count[order(pair_can_count$Freq, decreasing = T),]

tmp2 <- tmp2[tmp2$pair %in% pair_can_count$Var1[pair_can_count$Freq == length(cancers_sort)],]

for (enzyme_type in c("kinase")) {
  tmp <- tmp2
  if (nrow(tmp) > 0) {
    ## create output directory
    outDir1 <- paste0(makeOutDir(resultD = resultD), "down/")
    dir.create(outDir1)
    for (enzyme in unique(tmp$GENE)) {
      for (substrate in unique(tmp$SUB_GENE[tmp$GENE == enzyme])) {
        for (rsd in unique(tmp$SUB_MOD_RSD[tmp$SUB_GENE == substrate])) {
          sup_3can <- NULL
          for (cancer2 in c("BRCA", "OV", "CO")) {
            ## input susbstrate phospho level
            pho_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
            pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
            pho_t <- pho_data[pho_data$Gene == substrate & pho_head$SUB_MOD_RSD == rsd,]
            
            pho_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PHO_formatted_normalized_Control.txt",sep=""), data.table = F)
            pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
            pho_n <- pho_data[pho_data$Gene == substrate & pho_head$SUB_MOD_RSD == rsd,]
            
            rm(pho_data)
            
            ## input kinase phospho level
            phog_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
            phog_t <- phog_data[phog_data$Gene == enzyme,]
            
            phog_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_collapsed_PHO_formatted_normalized_replicate_averaged_Normal.txt",sep=""), data.table = F)
            phog_n <- phog_data[phog_data$Gene == enzyme,]
            
            rm(phog_data)
            
            if (nrow(pho_t) > 0 & nrow(pho_n) > 0 & nrow(phog_t) & nrow(phog_n)) {
              pho_t.m <- melt(pho_t)
              colnames(pho_t.m)[ncol(pho_t.m)] <- "pho_sub"
              pho_t.m$tumor_normal <- "tumor"
              
              pho_n.m <- melt(pho_n)
              colnames(pho_n.m)[ncol(pho_n.m)] <- "pho_sub"
              pho_n.m$tumor_normal <- "normal"
              
              phog_t.m <- melt(phog_t)
              colnames(phog_t.m)[ncol(phog_t.m)] <- "phog_kin"
              phog_t.m$tumor_normal <- "tumor"
              
              phog_n.m <- melt(phog_n)
              colnames(phog_n.m)[ncol(phog_n.m)] <- "phog_kin"
              phog_n.m$tumor_normal <- "normal"
              
              tab_t <- merge(pho_t.m[, c("variable", "pho_sub", "tumor_normal")], phog_t.m[, c("variable", "phog_kin", "tumor_normal")], all = T)
              tab_n <- merge(pho_n.m[, c("variable", "pho_sub", "tumor_normal")], phog_n.m[, c("variable", "phog_kin", "tumor_normal")], all = T)
              sup_tab <- rbind(tab_t, tab_n)
              sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
              
              sup_tab$cancer <- cancer2
              sup_3can <- rbind(sup_3can, sup_tab)
            }
          }
          
          if (length(unique(sup_3can$cancer)) == length(cancers_sort)) {
            ## make deeper directory for outpus
            outDir2 <- paste0(outDir1, enzyme, '_', substrate, '_', rsd, "/")
            dir.create(outDir2)
            
            ##  add label
            sup_3can$cancer <- factor(sup_3can$cancer, levels = cancers_sort)
            
            ## scatterplot
            p = ggplot(sup_3can, aes(x=phog_kin, y=pho_sub))
            p = p + geom_point(aes(color = tumor_normal), stroke = 0, alpha = 0.8, shape = 16)
            p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
            p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
            p <- p + facet_grid(.~cancer, scales = "free_y", space = "fixed")
            p = p + scale_color_manual(values = color_tn) 
            p = p + labs(x = paste0(enzyme, " phosphorylation abundance(log2 ratio)"), 
                         y = paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio"))
            p = p + theme_nogrid()
            p = p + theme(axis.title = element_text(size=6),
                          axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                          axis.text.y = element_text(colour="black", size=6))
            p = p + theme(title = element_text(size = 8),
                          legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                          legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                          legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                          legend.title = element_text(size = 5), legend.justification = c(0, 1), legend.position = c(0, 1))
            p
            fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_sub+pho_kin~pho_sub.pdf")
            ggsave(file=fn, height=2.5, width=6, useDingbats=FALSE)
          }
        }
      }
    }
  }
}

# regulated kinase-substrate pairs -------------------------------------

for (enzyme_type in c("kinase")) {
  tmp <- sup_cans_tab_pk
  tmp <- tmp[tmp$SELF == 'trans',]
  # tmp <- tmp[tmp$known_enzyme_phosphosite & !is.na(tmp$known_enzyme_phosphosite) & !is.na(tmp$Source) & tmp$Source != "NetKIN",]
  tmp <- tmp[tmp$shared3can & !is.na(tmp$shared3can),]
  
  if (nrow(tmp) > 0) {
    ## create output directory
    outDir1 <- paste0(makeOutDir(resultD = resultD), "regulated/")
    dir.create(outDir1)
    for (enzyme in unique(tmp$GENE)) {
      for (substrate in unique(tmp$SUB_GENE[tmp$GENE == enzyme])) {
        for (rsd in unique(tmp$SUB_MOD_RSD[tmp$SUB_GENE == substrate])) {
          sup_3can <- NULL
          for (cancer2 in c("BRCA", "OV", "CO")) {
            ## input susbstrate phospho level
            pho_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
            pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
            pho_t <- pho_data[pho_data$Gene == substrate & pho_head$SUB_MOD_RSD == rsd,]
            
            pho_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_PHO_formatted_normalized_Control.txt",sep=""), data.table = F)
            pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
            pho_n <- pho_data[pho_data$Gene == substrate & pho_head$SUB_MOD_RSD == rsd,]
            
            rm(pho_data)
            
            ## input kinase phospho level
            phog_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
            phog_t <- phog_data[phog_data$Gene == enzyme,]
            
            phog_data <- fread(input = paste(cptac_sharedD, cancer2,"/",prefix[cancer2], "_collapsed_PHO_formatted_normalized_replicate_averaged_Normal.txt",sep=""), data.table = F)
            phog_n <- phog_data[phog_data$Gene == enzyme,]
            
            rm(phog_data)
            
            if (nrow(pho_t) > 0 & nrow(pho_n) > 0 & nrow(phog_t) & nrow(phog_n)) {
              pho_t.m <- melt(pho_t)
              colnames(pho_t.m)[ncol(pho_t.m)] <- "pho_sub"
              pho_t.m$tumor_normal <- "tumor"
              
              pho_n.m <- melt(pho_n)
              colnames(pho_n.m)[ncol(pho_n.m)] <- "pho_sub"
              pho_n.m$tumor_normal <- "normal"
              
              phog_t.m <- melt(phog_t)
              colnames(phog_t.m)[ncol(phog_t.m)] <- "phog_kin"
              phog_t.m$tumor_normal <- "tumor"
              
              phog_n.m <- melt(phog_n)
              colnames(phog_n.m)[ncol(phog_n.m)] <- "phog_kin"
              phog_n.m$tumor_normal <- "normal"
              
              tab_t <- merge(pho_t.m[, c("variable", "pho_sub", "tumor_normal")], phog_t.m[, c("variable", "phog_kin", "tumor_normal")], all = T)
              tab_n <- merge(pho_n.m[, c("variable", "pho_sub", "tumor_normal")], phog_n.m[, c("variable", "phog_kin", "tumor_normal")], all = T)
              sup_tab <- rbind(tab_t, tab_n)
              sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
              
              sup_tab$cancer <- cancer2
              sup_3can <- rbind(sup_3can, sup_tab)
            }
          }
          
          if (length(unique(sup_3can$cancer)) == length(cancers_sort)) {
            ## make deeper directory for outpus
            outDir2 <- paste0(outDir1, enzyme, '_', substrate, '_', rsd, "/")
            dir.create(outDir2)
            
            ##  add label
            sup_3can$cancer <- factor(sup_3can$cancer, levels = cancers_sort)
            
            ## scatterplot
            p = ggplot(sup_3can, aes(x=phog_kin, y=pho_sub))
            p = p + geom_point(aes(color = tumor_normal), stroke = 0, alpha = 0.8, shape = 16)
            p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
            p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
            p <- p + facet_grid(.~cancer, scales = "free_y", space = "fixed")
            p = p + scale_color_manual(values = color_tn) 
            p = p + labs(x = paste0(enzyme, " phosphorylation abundance(log2 ratio)"), 
                         y = paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio"))
            p = p + theme_nogrid()
            p = p + theme(axis.title = element_text(size=6),
                          axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                          axis.text.y = element_text(colour="black", size=6))
            p = p + theme(title = element_text(size = 8),
                          legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                          legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                          legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                          legend.title = element_text(size = 5), legend.justification = c(0, 1), legend.position = c(0, 1))
            p
            fn = paste0(outDir2, enzyme, "_", substrate, "_", rsd, "_pro_sub+pho_kin~pho_sub.pdf")
            ggsave(file=fn, height=2.5, width=6, useDingbats=FALSE)
          }
        }
      }
    }
  }
}
