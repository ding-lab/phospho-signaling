# Yige Wu @ WashU 2018 Apr
# filter the regression results by number of non-NA values and adjust FDR accordingly

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source('cptac2p_analysis/phospho_network/phospho_network_shared.R')



# inputs regression table------------------------------------------------------------------
if (!exists("table2changenonNA")) {
  stop("input regression table first!")
}
# table2changenonNA <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/gather_regression_super_table/regression_cptac2p_cptac3_tumor.txt"), data.table = F)


# annotate with source for kinase-substrate relations ----------------------------------------------------
table2annotate <- table2changenonNA
table2annotate$GENE <- table2annotate$KINASE
table2annotate$SUB_GENE <- table2annotate$SUBSTRATE
source("./cptac2p_analysis/phospho_network/regression/tables/annotate_regression_with_source.R")
source("./cptac2p_analysis/phospho_network/regression/tables/annotate_regression_with_enzyme_type.R")
table2changenonNA <- table2annotate

# set variables -----------------------------------------------------------
if (!exists("reg_nonNA")) {
  stop("input reg_nonNAfirst!")
}

# business ------------------------------------------------------------------
tab_tmp <- table2changenonNA[table2changenonNA$Size >= reg_nonNA,]

## adjust p-values to FDR
name = c("pro_kin","pro_sub","pho_kin")
for (cancer_tmp in unique(tab_tmp$Cancer)) {
  for (enzyme_type_tmp in unique(tab_tmp$enzyme_type)) {
    for(self_tmp in c(TRUE,FALSE)) {
      for(coln_tmp in name) {#adjust pvalues for each variable
        row <- (tab_tmp$self==self_tmp) & (tab_tmp$Cancer==cancer_tmp) & (tab_tmp$enzyme_type==enzyme_type_tmp)
        tab_tmp[row, paste("FDR_",coln_tmp,sep = "")] <-p.adjust(tab_tmp[row,paste("P_",coln_tmp,sep = "")],method = "fdr")
      }
    }
  }
}

## clean up
tab_tmp <- tab_tmp[!(tab_tmp$GENE %in% c("APC", "AXIN1", "CCNE1", "CCNA2")),]
table2changenonNA <- tab_tmp

# for (reg_nonNA in seq(from = 15, to = 25, by = 5)) {
#   tab_tmp <- table2changenonNA[table2changenonNA$Size >= reg_nonNA,]
#   
#   ## adjust p-values to FDR
#   name = c("pro_kin","pro_sub","pho_kin")
#   for (cancer in unique(tab_tmp$Cancer)) {
#     for (enzyme_type in unique(tab_tmp$enzyme_type)) {
#       for(self in c(TRUE,FALSE)) {
#         for(coln in name) {#adjust pvalues for each variable
#           row <- (tab_tmp$self==self) & (tab_tmp$Cancer==cancer) & (tab_tmp$enzyme_type==enzyme_type)
#           tab_tmp[row, paste("FDR_",coln,sep = "")] <-p.adjust(tab_tmp[row,paste("P_",coln,sep = "")],method = "fdr")
#         }
#       }
#     }
#   }
#   ## clean up
#   tab_tmp <- tab_tmp[!(tab_tmp$GENE %in% c("APC", "AXIN1", "CCNE1", "CCNA2")),]
#   
#   write.table(x = tab_tmp, file = paste0(makeOutDir(resultD = resultD),
#                                          "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"),
#               row.names = F, quote = F, sep = "\t")
# }
#   
# ## get common detected phosphosites
# num_nonNA_common <- reg_nonNA
# pho_sites_cans <- list()
# pho_sites_all <- NULL
# for (cancer in cancers_sort) {
#   pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
#   samples <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
#   pho_data <- pho_data[rowSums(!is.na(pho_data[, samples])) >= num_nonNA_common,]
#   pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
#   pho_sites <- unique(pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")])
#   pho_sites$site <- paste0(pho_sites$SUBSTRATE, "_", pho_sites$SUB_MOD_RSD)
#   pho_sites_cans[[cancer]] <- unique(as.vector(pho_sites$site))
#   pho_sites_all <- unique(c(pho_sites_all, as.vector(pho_sites$site)))
# }
# dat <- data.frame(site = pho_sites_all)
# for (cancer in cancers_sort) {
#   dat[, cancer] <- FALSE
#   dat[dat$site %in% pho_sites_cans[[cancer]], cancer]  <- TRUE
# }
# pho_site_common <- dat[dat$BRCA & dat$OV & dat$CO,]
# 
# ## filter by phosphosite
# tab_tmp$site <- paste0(tab_tmp$SUB_GENE, "_", tab_tmp$SUB_MOD_RSD)
# tab_tmp <- tab_tmp[tab_tmp$site %in% pho_site_common$site,]
# 
# write.table(x = tab_tmp, file = paste0(makeOutDir(resultD = resultD), 
#                                        enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
#                                        "_reg_nonNA", reg_nonNA, "_commonsite_nonNA", num_nonNA_common,".txt"), row.names = F, quote = F, sep = "\t")


# split table for smaller supplementary tables ----------------------------
# table2changenonNA <- NULL
# for (enzyme_type in c("kinase", "phosphatase")) {
#   tab_tmp <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_manualAdded_protein_level/",
#                                           enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor.txt"), data.table = F)
#   
#   for (cancer in cancers2process) {
#     for (self in c("trans", "cis")) {
#       tab2w <- tab_tmp[tab_tmp$Cancer == cancer & tab_tmp$SELF == self,]
#       tab2w$self <- NULL
#       tab2w$pair <- NULL
#       
#       write.table(x = tab2w, file = paste0(makeOutDir(resultD = resultD), enzyme_type, "_substrate_regression_", cancer, "_", self , "_tumor.txt"), quote = F, row.names = F, sep = "\t")
#       write.table(x = tab2w, file = paste0(makeOutDir(resultD = resultD), enzyme_type, "_substrate_regression_", cancer, "_", self , "_tumor.csv"), quote = F, row.names = F, sep = ",")
#       
#       if (self == "trans") {
#         tab2w <- tab2w[tab2w$FDR_pho_kin < 0.1,]
#         write.table(x = tab2w, file = paste0(makeOutDir(resultD = resultD), enzyme_type, "_substrate_regression_", cancer, "_", self , "_FDR0.1", "_tumor.csv"), quote = F, row.names = F, sep = ",")
#         
#         tab2w <- tab2w[tab2w$FDR_pho_kin < 0.05,]
#         write.table(x = tab2w, file = paste0(makeOutDir(resultD = resultD), enzyme_type, "_substrate_regression_", cancer, "_", self , "_FDR0.05", "_tumor.csv"), quote = F, row.names = F, sep = ",")
#         
#       }
#     }
#   }
#   table2changenonNA <- rbind(table2changenonNA, tab_tmp)
# }

