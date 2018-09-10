# Yige Wu @ WashU 2018 Jan
# integrate differential expression data with KSEA results

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

## input enzyme substrate table
ptms_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_networkin_enzyme_substrate_site_level_union.csv")
op <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_genename_annotated.csv")
op$enzyme_type <- NA
op$enzyme_type[op$modification == "phosphorylation"] <- "kinase"
op$enzyme_type[op$modification == "dephosphorylation"] <- "phosphotase"
op_enzyme <- unique(op[op$enzyme_type == "kinase" | op$enzyme_type == "phosphotase", c("GENE", "enzyme_type")])
write.table(x = op_enzyme,
            file = paste0(makeOutDir(), "omnipath_enzyme.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")


# inputs ------------------------------------------------------------------
for (m.cutoff in 2:2) {
  for (NetworKIN.cutoff in 5:5) {
    ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$networkin_score >= NetworKIN.cutoff,]
    for (sig in c(0.2)) {
      ## input differential expression
      Pho.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "PHO_diffexp_3can_paired_FDR", sig, ".txt"),
                               data.table = F)
      Pho.diffexp3can <- cbind(Pho.diffexp3can, formatPhosphosite(phosphosite_vector = Pho.diffexp3can$Phosphosite, Pho.diffexp3can$Gene))
      
      for (cancer in c("BRCA", "OV", "CO")) {
        ## input KSEA score
        ksea <- fread(input = paste0("/Users/yigewu/Box\ Sync/", "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/OmniPath_NetworKIN_", 
                                     "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "/", cancer, "/KSEA Kinase Scores.csv"),
                      data.table = F)
        ksea <- ksea[ksea$p.value < 0.05 & ksea$m >= m.cutoff,]
        pho <- Pho.diffexp3can[Pho.diffexp3can$Cancer == cancer,]
        ptms_table <- merge(ptms_site_pairs_sup, ksea[, c('z.score', "p.value", "Kinase.Gene")], by.x = c("GENE"), by.y = c("Kinase.Gene"))
        colnames(ptms_table)[(-1):0+ncol(ptms_table)] <- c("kinase_zscore", "kinase_pvalue")
        ptms_table <- merge(ptms_table, pho[,c("SUBSTRATE", "SUB_MOD_RSD", "Fold Change")], by.x = c("SUB_GENE","SUB_MOD_RSD"), by.y = c("SUBSTRATE", "SUB_MOD_RSD"))
        colnames(ptms_table)[ncol(ptms_table)] <- c("substrate_FC")
        ptms_table <- merge(ptms_table, op_enzyme, by = c("GENE"), all.x = T)
        ptms_table$enzyme_type[is.na(ptms_table$enzyme_type)] <- "kinase"
        ptms_table$enzyme_direction <- ifelse(ptms_table$kinase_zscore > 0, "up", "down")
        ptms_table$substrate_direction <- ifelse(ptms_table$substrate_FC > 1, "up", "down")
        ptms_table$same_direction <- (as.vector(ptms_table$enzyme_direction) == as.vector(ptms_table$substrate_direction))
        ptms_table$consistent <- ifelse((ptms_table$enzyme_type == "kinase" & ptms_table$same_direction) | (ptms_table$enzyme_type == "phosphotase" & !ptms_table$same_direction), TRUE, FALSE)
        
        ## annoate drivers
        tsg <- loadGeneList(id = "tsg", cancer = cancer, is.soft.limit = "soft")
        oncogene <- loadGeneList(id = "oncogene", cancer = cancer, is.soft.limit = "soft")
        ptms_table$enzyme.driver_type <- "NA"
        ptms_table$substrate.driver_type <- "NA"
        
        ptms_table$enzyme.driver_type[ptms_table$GENE %in% tsg] <- "tsg"
        ptms_table$enzyme.driver_type[ptms_table$GENE %in% oncogene] <- "oncogene"
        ptms_table$substrate.driver_type[ptms_table$SUB_GENE %in% tsg] <- "tsg"
        ptms_table$substrate.driver_type[ptms_table$SUB_GENE %in% oncogene] <- "oncogene"
        
        
        write.table(x = ptms_table,
                    file = paste0(makeOutDir(), cancer, "_KSEA_enzyme_diffexp_substrate.txt"),
                    row.names = F, col.names = T, quote = F, sep = "\t")
        
        ## write out consistent directed enzyme-substrate pairs
        ptms_table_consistent <- ptms_table[ptms_table$consistent,]
        write.table(x = ptms_table_consistent,
                    file = paste0(makeOutDir(), cancer, "_KSEA_enzyme_diffexp_substrate_consistent.txt"),
                    row.names = F, col.names = T, quote = F, sep = "\t")
        
        ## write out consistent enzyme-substrate pairs and also in concord with its driver role
        ptms_table_consistent_driver <- ptms_table_consistent[((ptms_table_consistent$substrate_direction == "up" & ptms_table_consistent$substrate.driver_type == "oncogene") | (ptms_table_consistent$enzyme_direction == "up" & ptms_table_consistent$enzyme.driver_type == "oncogene")) | ((ptms_table_consistent$substrate_direction == "down" & ptms_table_consistent$substrate.driver_type == "tsg") | (ptms_table_consistent$enzyme_direction == "down" & ptms_table_consistent$enzyme.driver_type == "tsg")),]
        write.table(x = ptms_table_consistent_driver,
                    file = paste0(makeOutDir(), cancer, "_KSEA_enzyme_diffexp_substrate_consistent_driver.txt"),
                    row.names = F, col.names = T, quote = F, sep = "\t")
      }
    }
  }
}



