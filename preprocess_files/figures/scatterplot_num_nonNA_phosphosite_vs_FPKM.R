# Yige Wu @ WashU 2018 Apr
# venn diagram showing the phosphosites with at least 5 non-NA values

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
# devtools::install_github("jolars/eulerr")
library(eulerr)
library(ggrepel)

# inputs ------------------------------------------------------------------
## input enzyme-substrate pairs examined
en_sub_cans <- NULL
for (enzyme_type in c("phosphatase", "kinase")) {
  sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/", 
                                          enzyme_type, "_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
  sup_cans_tab_en$GENE <- as.vector(sup_cans_tab_en$KINASE)
  sup_cans_tab_en$SUB_GENE <- as.vector(sup_cans_tab_en$SUBSTRATE)
  en_sub_can <- unique(sup_cans_tab_en[, c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer")])
  en_sub_can$enzyme_type <- enzyme_type
  en_sub_cans <- rbind(en_sub_cans, en_sub_can)
}
## input RPKM

# set variables -----------------------------------------------------------
# num_nonNA2test <- seq(from = 5, to = 30, by = 5)
num_nonNA2test <- c(5,25)
num_sites2test <- c(1)
subdir_phosphosite <- paste0(makeOutDir(resultD = resultD), "phosphosite/")
dir.create(subdir_phosphosite)
subdir_phosphopro <- paste0(makeOutDir(resultD = resultD), "phosphopro/")
dir.create(subdir_phosphopro)
subdir_sites_per_phosphopro <- paste0(makeOutDir(resultD = resultD), "sites_per_phosphopro/")
dir.create(subdir_sites_per_phosphopro)
pos <- position_jitter(width = 0.5, seed = 1)
## threshold to time SD below median
pro_sd_thres <- 2
fpkm_sd_thres <- 1

# plot the distribution of median fpkm abundance per gene with outlier threshold------------------------------
fpkm_median_bottom <- NULL
fpkm_medians_cans <- NULL
for (cancer in cancers_sort) {
  fpkm_can <- readRDS(file = paste0("./cptac2p/analysis_results/expression_matrices/tables/get_mRNA_quantile/cptac2p_mRNA_score.RDS"))
  fpkm_can <- fpkm_can[[cancer]]
  fpkm_can$log2FPKM <- log2(as.vector(fpkm_can$score) + 0.0001)
  fpkm_mat <- dcast(fpkm_can, formula = Gene ~ Participant.ID, value.var = "log2FPKM")
  fpkm_medians <- data.frame(gene = as.vector(fpkm_mat$Gene),
                             fpkm_median = rowMedians(as.matrix(fpkm_mat[, -1]), na.rm = T),
                             cancer = cancer)
  
  fpkm_medians_median <- median(x = fpkm_medians$fpkm_median, na.rm = T)
  fpkm_medians_sd <- sd(x = fpkm_medians$fpkm_median, na.rm = T)

  fpkm_median_bottom[cancer] <- (fpkm_medians_median - fpkm_sd_thres*fpkm_medians_sd)
  
  fpkm_medians$bottom <- (fpkm_medians$fpkm_median < fpkm_median_bottom[cancer])
  fpkm_medians_cans <- rbind(fpkm_medians_cans, fpkm_medians)
  
  p <- ggplot()
  p <- p + geom_density(data = fpkm_medians, mapping = aes(x = fpkm_median))
  p <- p + geom_vline(xintercept = fpkm_median_bottom[cancer], linetype = 2)
  p <- p + geom_text(mapping = aes(x = fpkm_median_bottom[cancer], y = 0.1, label = signif(fpkm_median_bottom[cancer], digits = 2)), nudge_x = -1)
  p <- p + theme_nogrid()
  p <- p + xlim(c(-5, 5))
  p
  ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer,  "_fpkm_abundance_low_outlier_cutoff_", fpkm_sd_thres, "sd.pdf"),
         width = 5, height = 4)
}
## check whether have any driver gene filtered out
print(fpkm_medians_cans$gene[fpkm_medians_cans$bottom & fpkm_medians_cans$gene %in% loadGeneList(gene_type = "driver", cancer = "BRCA", is.soft.limit = "")])
print(fpkm_medians_cans$gene[fpkm_medians_cans$bottom & fpkm_medians_cans$gene %in% loadGeneList(gene_type = "driver", cancer = "OV", is.soft.limit = "")])
print(fpkm_medians_cans$gene[fpkm_medians_cans$bottom & fpkm_medians_cans$gene %in% loadGeneList(gene_type = "driver", cancer = "CO", is.soft.limit = "")])
print(fpkm_medians_cans$gene[fpkm_medians_cans$bottom & fpkm_medians_cans$gene %in% loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")])



# plot the distribution of median protein abundance per gene with outlier threshold------------------------------
protein_median_bottom <- NULL
protein_medians_cans <- NULL
for (cancer in cancers_sort) {
  protein <- loadProteinNormalizedTumor(cancer = cancer)
  samples <- colnames(protein)
  samples <- samples[!(samples == "Gene")]
  protein_medians <- rowMedians(as.matrix(protein[, samples]), na.rm = T)
  protein_medians <- data.frame(gene = as.vector(protein$Gene),
                                protein_median = protein_medians,
                                cancer = cancer)

  protein_medians_median <- median(x = protein_medians$protein_median, na.rm = T)
  protein_medians_sd <- sd(x = protein_medians$protein_median, na.rm = T)

  protein_median_bottom[cancer] <- (protein_medians_median - pro_sd_thres*protein_medians_sd)
  
  protein_medians$bottom <- (protein_medians$protein_median < protein_median_bottom[cancer])
  protein_medians_cans <- rbind(protein_medians_cans, protein_medians)
  
  p <- ggplot()
  p <- p + geom_density(data = protein_medians, mapping = aes(x = protein_median))
  p <- p + geom_vline(xintercept = protein_median_bottom[cancer], linetype = 2)
  p <- p + geom_text(mapping = aes(x = protein_median_bottom[cancer], y = 0.1, label = signif(protein_median_bottom[cancer], digits = 2)), nudge_x = -1)
  p <- p + theme_nogrid()
  p <- p + xlim(c(-5, 5))
  p
  ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer,  "_protein_abundance_low_outlier_cutoff_", pro_sd_thres, "sd.pdf"),
         width = 5, height = 4)
}
## check whether have any driver gene filtered out
print(protein_medians_cans$gene[protein_medians_cans$bottom & protein_medians_cans$gene %in% loadGeneList(gene_type = "driver", cancer = "BRCA", is.soft.limit = "")])
print(protein_medians_cans$gene[protein_medians_cans$bottom & protein_medians_cans$gene %in% loadGeneList(gene_type = "driver", cancer = "OV", is.soft.limit = "")])
print(protein_medians_cans$gene[protein_medians_cans$bottom & protein_medians_cans$gene %in% loadGeneList(gene_type = "driver", cancer = "CO", is.soft.limit = "")])
print(protein_medians_cans$gene[protein_medians_cans$bottom & protein_medians_cans$gene %in% loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")])

# plot all phosphosites ---------------------------------------------------
## input phosphosite data from 3 cancers
gene_cat <- "all"
subdir_phosphosite_cat <- paste0(subdir_phosphosite, gene_cat, "/")
dir.create(subdir_phosphosite_cat)
subdir_phosphopro_cat <- paste0(subdir_phosphopro, gene_cat, "/")
dir.create(subdir_phosphopro_cat)
subdir_sites_per_phosphopro_cat <- paste0(subdir_sites_per_phosphopro, gene_cat, "/")
dir.create(subdir_sites_per_phosphopro_cat)

for (num_nonNA in num_nonNA2test) {
  for (num_sites in num_sites2test) {
    pho_sites_cans <- list()
    pho_sites_all <- NULL
    pro_all <- NULL
    pro_cans <- list()
    for (cancer in cancers_sort) {
      phog_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
      
      
      pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
      samples <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
      pho_data <- pho_data[rowSums(!is.na(pho_data[, samples])) >= num_nonNA,]
      pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
      pho_sites <- unique(pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")])
      pho_sites$site <- paste0(pho_sites$SUBSTRATE, "_", pho_sites$SUB_MOD_RSD)
      pho_sites_cans[[cancer]] <- unique(as.vector(pho_sites$site))
      pho_sites_all <- unique(c(pho_sites_all, as.vector(pho_sites$site)))
      
      pro_tab <- data.frame(table(pho_sites$SUBSTRATE))
      pro_tab <- pro_tab[pro_tab$Freq >= num_sites,]
      pro_cans[[cancer]] <- as.vector(pro_tab$Var1)
      pro_all <- unique(c(pro_all, as.vector(pro_tab$Var1)))
    }
    
    ## scatterplot for sites phosphoprotein
    dat_pro <- data.frame(gene = pro_all)
    dat_phosite <- data.frame(site = pho_sites_all)
    dat_phosite$gene <- str_split_fixed(string = dat_phosite$site, pattern = "_", 2)[,1]
    for (cancer in cancers_sort) {
      sites_can <- data.frame(site = pho_sites_cans[[cancer]])
      sites_can$gene <- str_split_fixed(sites_can$site, pattern = "_", 2)[,1]
      sites_can_gene <- data.frame(table(sites_can$gene))
      dat_pro <- merge(dat_pro, sites_can_gene, by.x = c("gene"), by.y = c("Var1"), all.x = T)
      colnames(dat_pro)[ncol(dat_pro)] <- paste0("sites_per_pro_", cancer)
      tmp <- as.vector(dat_pro[, paste0("sites_per_pro_", cancer)])
      tmp[is.na(tmp)] <- 0
      dat_pro[, paste0("sites_per_pro_", cancer)] <- tmp
      dat_phosite[, paste0("detected_", cancer)] <- FALSE
      dat_phosite[dat_phosite$site %in% pho_sites_cans[[cancer]], paste0("detected_", cancer)]  <- TRUE
    }
    dat_phosite <- merge(dat_phosite, dat_pro, all.x = T, by = c("gene"))
    dat_phosite <- unique(dat_phosite)

    dat_phosite <- merge(dat_phosite, dat_pro, all.x = T)
    
    for (cancer in cancers_sort)
    ## per cancer, number of phosphosites detected per protein vs median protein abundance
    
    ## per cancer number of phosphosites detected per protein vs median log2FPKM
    
    ## per cancer median log2FPKM vs median protein abundance, number of phosphosite detected
    
    
    
    ## BRCA vs OV number of phosphosites detected per protein
    tab2p <- dat_phosite
    detected_site_cat <- vector(mode = "character", length = nrow(dat_phosite))
    detected_site_cat[tab2p$detected_BRCA == TRUE] <- "detected_BRCA"
    detected_site_cat[tab2p$detected_OV == TRUE] <- "detected_OV"
    detected_site_cat[tab2p$detected_BRCA == TRUE & tab2p$detected_OV == TRUE] <- "detected_BRCA_OV"
    tab2p$detected_site_cat <- detected_site_cat
    tab2p <- tab2p[!is.na(tab2p$detected_site_cat) & tab2p$detected_site_cat != "",]
    fn = paste(subdir_sites_per_phosphopro_cat, 'sites_per_pro_BRCA_OV_', num_sites, "phosphosites_", num_nonNA, 'nonNA_scatterplot','.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 5, width = 12, useDingbats = FALSE)
    p <- ggplot()
    p <- p + geom_point(data = tab2p, mapping = aes(x = sites_per_pro_BRCA, y = sites_per_pro_OV, color = detected_site_cat), alpha = 0.3, position = pos)
    p <- p + coord_fixed()
    p <- p + geom_abline(slope = 1, linetype = 2, color = "grey")
    grid.draw(p)
    dev.off()
    ## BRCA vs CO number of phosphosites detected per protein
    

    
    
  }
}

# substrate phosphosites in pairs -----------------------------------------
gene_cat <- "substrate"
subdir_phosphosite_cat <- paste0(subdir_phosphosite, gene_cat, "/")
dir.create(subdir_phosphosite_cat)
subdir_phosphopro_cat <- paste0(subdir_phosphopro, gene_cat, "/")
dir.create(subdir_phosphopro_cat)
subdir_sites_per_phosphopro_cat <- paste0(subdir_sites_per_phosphopro, gene_cat, "/")
dir.create(subdir_sites_per_phosphopro_cat)

for (num_nonNA in c(5)) {
  for (num_sites in num_sites2test) {
    pho_sites_cans <- list()
    pho_sites_all <- NULL
    pro_all <- NULL
    pro_cans <- list()
    for (cancer in cancers_sort) {
      pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
      samples <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
      pho_data <- pho_data[(rowSums(!is.na(pho_data[, samples])) >= num_nonNA) & (pho_data$Gene %in% en_sub_cans$SUB_GENE[en_sub_cans$Cancer == cancer]),]
      pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
      pho_sites <- unique(pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")])
      pho_sites$site <- paste0(pho_sites$SUBSTRATE, "_", pho_sites$SUB_MOD_RSD)
      pho_sites_cans[[cancer]] <- as.vector(pho_sites$site)
      pho_sites_all <- unique(c(pho_sites_all, as.vector(pho_sites$site)))
      
      pro_tab <- data.frame(table(pho_sites$SUBSTRATE))
      pro_tab <- pro_tab[pro_tab$Freq >= num_sites,]
      pro_cans[[cancer]] <- as.vector(pro_tab$Var1)
      pro_all <- c(pro_all, as.vector(pro_tab$Var1))
    }
    
    ## venn diagram for phosphosites
    dat <- data.frame(site = pho_sites_all)
    dat$gene <- str_split_fixed(string = dat$site, pattern = "_", 2)[,1]
    for (cancer in cancers_sort) {
      dat[, cancer] <- FALSE
      dat[dat$site %in% pho_sites_cans[[cancer]], cancer]  <- TRUE
      dat[, paste0("bottom_protein_", cancer)] <- (dat$gene %in% protein_medians_cans$gene[protein_medians_cans$cancer == cancer & protein_medians_cans$bottom])
      dat[, paste0("bottom_FPKM_", cancer)] <- (dat$gene %in% fpkm_medians_cans$gene[fpkm_medians_cans$cancer == cancer & fpkm_medians_cans$bottom])
      
    }
    
    for (cancer in c("BRCA")) {
      ## check whether phosphosites not detected in OV and CO is because low abundance
      fn = paste(subdir_phosphosite_cat, 'substrate_in_pairs_phosphosites_', num_nonNA, 'nonNA_', cancer,'_detected_overlap_bottom',pro_sd_thres, 'sd_protein_in_others_venn.pdf',sep ="")
      grid.newpage()
      pdf(fn, height = 12, width = 12, useDingbats = FALSE)
      fit <- euler(combinations = dat[, c(cancers_sort, paste0("bottom_protein_", cancers_sort[cancers_sort != cancer]))], input = "disjoint", shape = 'circle')
      p <-plot(fit, quantities = list(fontsize = 20), fills = color_cancers2, labels = F, legend = list(fontsize = 20))
      grid.draw(p)
      dev.off()
      
      fn = paste(subdir_phosphosite_cat, 'substrate_in_pairs_phosphosites_', num_nonNA, 'nonNA_', cancer,'_detected_overlap_bottom',fpkm_sd_thres, 'sd_log2FPKM_in_others_venn.pdf',sep ="")
      grid.newpage()
      pdf(fn, height = 12, width = 12, useDingbats = FALSE)
      fit <- euler(combinations = dat[, c(cancers_sort, paste0("bottom_logfpkm_", cancers_sort[cancers_sort != cancer]))], input = "disjoint", shape = 'circle')
      p <-plot(fit, quantities = list(fontsize = 20), fills = color_cancers2, labels = F, legend = list(fontsize = 20))
      grid.draw(p)
      dev.off()
    }

    ## venn diagram for phosphoprotein
    dat <- data.frame(site = pro_all)
    for (cancer in cancers_sort) {
      dat[, cancer] <- FALSE
      dat[dat$site %in% pro_cans[[cancer]], cancer]  <- TRUE
    }
    fn = paste(subdir_phosphopro_cat, 'substrate_in_pairs_phosphoproteins_', num_sites, "phosphosites_", num_nonNA, 'nonNA_venn','.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 12, useDingbats = FALSE)
    fit <- euler(combinations = dat, input = "disjoint", shape = 'circle')
    
    p <-plot(fit, quantities = list(fontsize = 50), fills = color_cancers2, labels = F, legend = list(fontsize = 50))
    grid.draw(p)
    dev.off()
    
    ## scatterplot for sites phosphoprotein
    dat_pro <- data.frame(gene = pro_all)
    dat_phosite <- data.frame(site = pho_sites_all)
    dat_phosite$gene <- str_split_fixed(string = dat_phosite$site, pattern = "_", 2)[,1]
    for (cancer in cancers_sort) {
      sites_can <- data.frame(site = pho_sites_cans[[cancer]])
      sites_can$gene <- str_split_fixed(sites_can$site, pattern = "_", 2)[,1]
      sites_can_gene <- data.frame(table(sites_can$gene))
      dat_pro <- merge(dat_pro, sites_can_gene, by.x = c("gene"), by.y = c("Var1"), all.x = T)
      colnames(dat_pro)[ncol(dat_pro)] <- paste0("sites_per_pro_", cancer)
      tmp <- as.vector(dat_pro[, paste0("sites_per_pro_", cancer)])
      tmp[is.na(tmp)] <- 0
      dat_pro[, paste0("sites_per_pro_", cancer)] <- tmp
      dat_phosite[, paste0("detected_", cancer)] <- FALSE
      dat_phosite[dat_phosite$site %in% pho_sites_cans[[cancer]], paste0("detected_", cancer)]  <- TRUE
    }
    dat_phosite <- merge(dat_phosite, dat_pro, all.x = T, by = c("gene"))
    dat_phosite <- unique(dat_phosite)
    
    ## OV number of phosphosites detected vs FPKM median
    # tab2p <- dat_phosite[dat_phosite$sites_per_pro_BRCA > 0 | dat_phosite$sites_per_pro_OV >0,]
    # tab2p <- unique(tab2p[, c("sites_per_pro_OV", "gene")])
    # fpkm_list <- readRDS(file = paste0("./cptac2p/analysis_results/expression_matrices/tables/get_mRNA_quantile/cptac2p_mRNA_score.RDS"))
    # fpkm_can <- fpkm_list[["OV"]]
    # rm(fpkm_list)
    # fpkm_can <- fpkm_can[fpkm_can$Gene %in% tab2p$gene,]
    # fpkm_can$log2fpkm <- log2(fpkm_can$score + 0.001)
    # tab2p$median_log2fpkm <- sapply(X = tab2p$gene, FUN = function(gene, fpkm_tab) median(fpkm_tab$log2fpkm[fpkm_tab$Gene == gene]),
    #                                 fpkm_tab = fpkm_can)
    # rm(fpkm_can)
    # 
    # 
    # fn = paste(subdir_sites_per_phosphopro_cat, 'sites_per_pro_vs_median_log2FPKM_OV_', num_sites, "phosphosites_", num_nonNA, 'nonNA_scatterplot','.pdf',sep ="")
    # grid.newpage()
    # pdf(fn, height = 5, width = 12, useDingbats = FALSE)
    # 
    # grid.draw(p)
    # dev.off()
    # 
    # p <- ggplot()
    # p <- p + geom_point(data = tab2p, mapping = aes(x = median_log2fpkm, y = sites_per_pro_OV, color = (sites_per_pro_OV == 0)), alpha = 0.3, position = pos)
    # p
    # 
    # p <- ggplot()
    # p <- p + geom_density(data = tab2p, mapping = aes(x = median_log2fpkm, group = (sites_per_pro_OV == 0), color = (sites_per_pro_OV == 0)), alpha = 0.3, position = pos)
    # p
    # 
    # protein <- loadProteinNormalizedTumor(cancer = "OV")
    # samples <- colnames(protein); samples <- samples[samples != "Gene"]
    # protein <- protein[protein$Gene %in% tab2p$gene,]
    # protein_exp <- data.frame(gene = protein$Gene,
    #                           median_pro = rowMedians(as.matrix(protein[,samples]), na.rm = T))
    # tab2p <- merge(tab2p, protein_exp, all.x = T)
    # 
    # p <- ggplot()
    # p <- p + geom_point(data = tab2p, mapping = aes(x = median_pro, y = sites_per_pro_OV, color = (sites_per_pro_OV == 0)), alpha = 0.3, position = pos)
    # p
    # 
    # p <- ggplot()
    # p <- p + geom_density(data = tab2p, mapping = aes(x = median_pro, group = (sites_per_pro_OV == 0), color = (sites_per_pro_OV == 0)), alpha = 0.3, position = pos)
    # p
  }
}

# enzyme phosphosites in pairs --------------------------------------------
for (enzyme_type in unique(en_sub_cans$enzyme_type)) {
  for (num_nonNA in num_nonNA2test) {
    for (num_sites in num_sites2test) {
      pho_sites_cans <- list()
      pho_sites_all <- NULL
      pro_all <- NULL
      pro_cans <- list()
      for (cancer in cancers_sort) {
        pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
        samples <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
        pho_data <- pho_data[(rowSums(!is.na(pho_data[, samples])) >= num_nonNA) & (pho_data$Gene %in% en_sub_cans$GENE[en_sub_cans$Cancer == cancer & en_sub_cans$enzyme_type == enzyme_type]),]
        pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
        pho_sites <- unique(pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")])
        pho_sites$site <- paste0(pho_sites$SUBSTRATE, "_", pho_sites$SUB_MOD_RSD)
        pho_sites_cans[[cancer]] <- as.vector(pho_sites$site)
        pho_sites_all <- unique(c(pho_sites_all, as.vector(pho_sites$site)))
        
        pro_tab <- data.frame(table(pho_sites$SUBSTRATE))
        pro_tab <- pro_tab[pro_tab$Freq >= num_sites,]
        pro_cans[[cancer]] <- as.vector(pro_tab$Var1)
        pro_all <- unique(c(pro_all, as.vector(pro_tab$Var1)))
      }
      
      ## venn diagram for phosphosites
      dat <- data.frame(site = pho_sites_all)
      for (cancer in cancers_sort) {
        dat[, cancer] <- FALSE
        dat[dat$site %in% pho_sites_cans[[cancer]], cancer]  <- TRUE
      }
      fn = paste(makeOutDir(resultD = resultD),  enzyme_type, '_in_pairs_phosphosites_', num_nonNA, 'nonNA_venn','.pdf',sep ="")
      grid.newpage()
      pdf(fn, height = 12, width = 12, useDingbats = FALSE)
      fit <- euler(combinations = dat, input = "disjoint", shape = 'circle')
      
      p <-plot(fit, quantities = list(fontsize = 50), fills = color_cancers2, labels = F, legend = list(fontsize = 50))
      grid.draw(p)
      dev.off()
      
      ## venn diagram for phosphoprotein
      dat <- data.frame(site = pro_all)
      for (cancer in cancers_sort) {
        dat[, cancer] <- FALSE
        dat[dat$site %in% pro_cans[[cancer]], cancer]  <- TRUE
        print(table(dat[, cancer]))
      }
      fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_in_pairs_phosphoproteins_', num_sites, "phosphosites_", num_nonNA, 'nonNA_venn','.pdf',sep ="")
      grid.newpage()
      pdf(fn, height = 12, width = 12, useDingbats = FALSE)
      fit <- euler(combinations = dat, input = "disjoint", shape = 'circle')
      
      p <-plot(fit, quantities = list(fontsize = 50), fills = color_cancers2, labels = F, legend = list(fontsize = 50))
      grid.draw(p)
      dev.off()
    }
  }
}




