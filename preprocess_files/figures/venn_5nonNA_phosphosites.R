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
fpkm_list <- readRDS(file = paste0("./cptac2p/analysis_results/expression_matrices/tables/get_mRNA_quantile/cptac2p_mRNA_score.RDS"))

# set variables -----------------------------------------------------------
num_nonNA2test <- seq(from = 5, to = 80, by = 5)
num_sites2test <- c(1)
subdir_phosphosite <- paste0(makeOutDir(resultD = resultD), "phosphosite/")
dir.create(subdir_phosphosite)
subdir_phosphopro <- paste0(makeOutDir(resultD = resultD), "phosphopro/")
dir.create(subdir_phosphopro)
subdir_sites_per_phosphopro <- paste0(makeOutDir(resultD = resultD), "sites_per_phosphopro/")
dir.create(subdir_sites_per_phosphopro)
pos <- position_jitter(width = 0.5, seed = 1)

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
    
    ## venn diagram for phosphosites
    dat <- data.frame(site = pho_sites_all)
    for (cancer in cancers_sort) {
      dat[, cancer] <- FALSE
      dat[dat$site %in% pho_sites_cans[[cancer]], cancer]  <- TRUE
      # print(table(dat[, cancer]))
    }

    fn = paste(subdir_phosphosite_cat, 'all_phosphosites_', num_nonNA, 'nonNA_venn','.pdf',sep ="")
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
    }

    fn = paste(subdir_phosphopro_cat, 'all_phosphoproteins_', num_sites, "phosphosites_", num_nonNA, 'nonNA_venn','.pdf',sep ="")
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
      
      dat_phosite[, paste0("detected_", cancer)] <- FALSE
      dat_phosite[dat_phosite$site %in% pho_sites_cans[[cancer]], paste0("detected_", cancer)]  <- TRUE
    }

    dat_phosite <- merge(dat_phosite, dat_pro, all.x = T)
    ## BRCA vs OV
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
  }
}

# plot phosphosite in driver genes ----------------------------------------
driver_genes <- fread("./TCGA_data/reference_files/Consensus.Genelist.full.txt", data.table = F, header = T)
for (driver_can in c("BRCA", "OV", "COADREAD", "PANCAN")) {
  drivers <- as.vector(driver_genes$Gene[driver_genes$`Cancer type` == driver_can])
  for (num_nonNA in num_nonNA2test[1]) {
    for (num_sites in num_sites2test) {
      pho_sites_cans <- list()
      pho_sites_all <- NULL
      pro_all <- NULL
      pro_cans <- list()
      for (cancer in cancers_sort) {
        phog_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
        
        
        pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
        samples <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
        pho_data <- pho_data[rowSums(!is.na(pho_data[, samples])) >= num_nonNA & pho_data$Gene %in% drivers,]
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
      
      ## venn diagram for phosphosites
      dat <- data.frame(site = pho_sites_all)
      for (cancer in cancers_sort) {
        dat[, cancer] <- FALSE
        dat[dat$site %in% pho_sites_cans[[cancer]], cancer]  <- TRUE
        print(table(dat[, cancer]))
      }
      fn = paste(makeOutDir(resultD = resultD), driver_can, '_driver_phosphosites_', num_nonNA, 'nonNA_venn','.pdf',sep ="")
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
      }
      fn = paste(makeOutDir(resultD = resultD), driver_can, '_driver_phosphoproteins_', num_sites, "phosphosites_", num_nonNA, 'nonNA_venn','.pdf',sep ="")
      grid.newpage()
      pdf(fn, height = 12, width = 12, useDingbats = FALSE)
      fit <- euler(combinations = dat, input = "disjoint", shape = 'circle')
      # p <-plot(fit, quantities = list(fontsize = 20), fills = color_cancers2, labels = list(fontsize = 20))
      p <-plot(fit, quantities = list(fontsize = 50), fills = color_cancers2, labels = F, legend = list(fontsize = 50))
      grid.draw(p)
      dev.off()
      
      print(paste0(dat$site[dat$BRCA & dat$OV & dat$CO], collapse = ","))
      print(paste0(dat$site[dat$BRCA & !dat$OV & !dat$CO], collapse = ","))
      print(paste0(dat$site[!dat$BRCA & dat$OV & !dat$CO], collapse = ","))
      print(paste0(dat$site[!dat$BRCA & !dat$OV & dat$CO], collapse = ","))
      
      print(paste0(dat$site[dat$BRCA & dat$OV & !dat$CO], collapse = ","))
      print(paste0(dat$site[dat$BRCA & !dat$OV & dat$CO], collapse = ","))
      print(paste0(dat$site[!dat$BRCA & dat$OV & dat$CO], collapse = ","))
      
    }
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
    for (cancer in cancers_sort) {
      dat[, cancer] <- FALSE
      dat[dat$site %in% pho_sites_cans[[cancer]], cancer]  <- TRUE
      print(table(dat[, cancer]))
    }
    fn = paste(subdir_phosphosite_cat, 'substrate_in_pairs_phosphosites_', num_nonNA, 'nonNA_venn','.pdf',sep ="")
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
    
    ## OV number of phosphosites detected vs FPKM
    tab2p <- dat_phosite
    tab2p <- unique(tab2p[, c("sites_per_pro_OV", "gene")])
    fpkm_can <- fpkm_list[["OV"]]
    fpkm_can <- fpkm_can[fpkm_can$Gene %in% tab2p$gene,]
    
    
    fn = paste(subdir_sites_per_phosphopro_cat, 'sites_per_pro_BRCA_OV_', num_sites, "phosphosites_", num_nonNA, 'nonNA_scatterplot','.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 5, width = 12, useDingbats = FALSE)
    p <- ggplot()
    p <- p + geom_point(data = tab2p, mapping = aes(x = sites_per_pro_BRCA, y = sites_per_pro_OV, color = detected_site_cat), alpha = 0.3, position = pos)
    p <- p + coord_fixed()
    p <- p + geom_abline(slope = 1, linetype = 2, color = "grey")
    p <- p + geom_text_repel(data = unique(tab2p[!is.na(tab2p$text), c("sites_per_pro_BRCA", "sites_per_pro_OV", "text")]), 
                             mapping = aes(x = sites_per_pro_BRCA, y = sites_per_pro_OV, label = text))
    grid.draw(p)
    dev.off()
    
    ## BRCA vs OV
    # tab2p <- dat_phosite
    # detected_site_cat <- vector(mode = "character", length = nrow(dat_phosite))
    # detected_site_cat[tab2p$detected_BRCA == TRUE] <- "detected_BRCA"
    # detected_site_cat[tab2p$detected_OV == TRUE] <- "detected_OV"
    # detected_site_cat[tab2p$detected_BRCA == TRUE & tab2p$detected_OV == TRUE] <- "detected_BRCA_OV"
    # tab2p$detected_site_cat <- detected_site_cat
    # tab2p <- tab2p[!is.na(tab2p$detected_site_cat) & tab2p$detected_site_cat != "",]
    # tab2p <- unique(tab2p[, c("sites_per_pro_BRCA", "sites_per_pro_OV", "detected_site_cat", "gene")])
    # tab2p$text <- NA
    # text_row <- (tab2p$sites_per_pro_BRCA >= quantile(x = tab2p$sites_per_pro_BRCA, probs = 0.99, na.rm = T) | tab2p$sites_per_pro_OV >= quantile(x = tab2p$sites_per_pro_OV, probs = 0.99, na.rm = T))
    # tab2p$text[text_row] <- as.vector(tab2p$gene[text_row])
    # fn = paste(subdir_sites_per_phosphopro_cat, 'sites_per_pro_BRCA_OV_', num_sites, "phosphosites_", num_nonNA, 'nonNA_scatterplot','.pdf',sep ="")
    # grid.newpage()
    # pdf(fn, height = 5, width = 12, useDingbats = FALSE)
    # p <- ggplot()
    # p <- p + geom_point(data = tab2p, mapping = aes(x = sites_per_pro_BRCA, y = sites_per_pro_OV, color = detected_site_cat), alpha = 0.3, position = pos)
    # p <- p + coord_fixed()
    # p <- p + geom_abline(slope = 1, linetype = 2, color = "grey")
    # p <- p + geom_text_repel(data = unique(tab2p[!is.na(tab2p$text), c("sites_per_pro_BRCA", "sites_per_pro_OV", "text")]), 
    #                          mapping = aes(x = sites_per_pro_BRCA, y = sites_per_pro_OV, label = text))
    # grid.draw(p)
    # dev.off()
    
    # tab2p <- dat_phosite
    # detected_site_cat <- vector(mode = "character", length = nrow(dat_phosite))
    # detected_site_cat[tab2p$detected_BRCA == TRUE] <- "detected_BRCA"
    # detected_site_cat[tab2p$detected_OV == TRUE] <- "detected_OV"
    # detected_site_cat[tab2p$detected_BRCA == TRUE & tab2p$detected_OV == TRUE] <- "detected_BRCA_OV"
    # tab2p$detected_site_cat <- detected_site_cat
    # tab2p <- tab2p[!is.na(tab2p$detected_site_cat) & tab2p$detected_site_cat != "",]
    # tab2p <- tab2p[tab2p$sites_per_pro_BRCA <= quantile(x = tab2p$sites_per_pro_BRCA, probs = 0.75, na.rm = T) & !is.na(tab2p$sites_per_pro_BRCA),]
    # fn = paste(subdir_sites_per_phosphopro_cat, 'sites_per_pro_BRCA_OV_', num_sites, "phosphosites_", num_nonNA, 'nonNA_scatterplot_faceted','.pdf',sep ="")
    # grid.newpage()
    # pdf(fn, height = 5, width = 12, useDingbats = FALSE)
    # p <- ggplot()
    # p <- p + geom_point(data = tab2p, mapping = aes(x = sites_per_pro_BRCA, y = sites_per_pro_OV, color = detected_site_cat), alpha = 0.3, position = pos)
    # p <- p + facet_grid(.~detected_site_cat, scales = "fixed", space = "fixed")
    # p <- p + coord_fixed()
    # p <- p + geom_abline(slope = 1, linetype = 2, color = "grey")
    # p <- p + geom_text_repel(data = unique(tab2p[tab2p$sites_per_pro_BRCA >= quantile(x = tab2p$sites_per_pro_BRCA, probs = 0.99, na.rm = T) & !is.na(tab2p$sites_per_pro_BRCA), c("sites_per_pro_BRCA", "sites_per_pro_OV", "gene")]), 
    #                          mapping = aes(x = sites_per_pro_BRCA, y = sites_per_pro_OV, label = gene))
    # grid.draw(p)
    # dev.off()
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




