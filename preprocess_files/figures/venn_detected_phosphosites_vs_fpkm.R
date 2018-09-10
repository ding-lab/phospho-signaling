# Yige Wu @ WashU 2018 Apr
# venn diagram showing the phosphosites with at least 5 non-NA values

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
# devtools::install_github("jolars/eulerr")
library(eulerr)

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

# set variables -----------------------------------------------------------
num_nonNA2test <- seq(from = 5, to = 30, by = 5)
num_sites2test <- c(1)

# plot all phosphosites ---------------------------------------------------
## input phosphosite data from 3 cancers
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
      print(table(dat[, cancer]))
    }
    fn = paste(makeOutDir(resultD = resultD), 'all_phosphosites_', num_nonNA, 'nonNA_venn','.pdf',sep ="")
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
    fn = paste(makeOutDir(resultD = resultD), 'all_phosphoproteins_', num_sites, "phosphosites_", num_nonNA, 'nonNA_venn','.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 12, useDingbats = FALSE)
    fit <- euler(combinations = dat, input = "disjoint", shape = 'circle')
    
    p <-plot(fit, quantities = list(fontsize = 50), fills = color_cancers2, labels = F, legend = list(fontsize = 50))
    grid.draw(p)
    dev.off()
  }
}
