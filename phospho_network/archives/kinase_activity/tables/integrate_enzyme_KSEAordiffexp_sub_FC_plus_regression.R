# Yige Wu @ WashU 2018 Apr
# integrate differential expression data with KSEA results

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

## input enzyme substrate table
ptms_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv")
op <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_genename_annotated.csv")
op$enzyme_type <- NA
op$enzyme_type[op$modification == "phosphorylation"] <- "kinase"
op$enzyme_type[op$modification == "dephosphorylation"] <- "phosphotase"
op_kinase <- unique(op[op$enzyme_type == "kinase" & !is.na(op$enzyme_type), c("GENE", "enzyme_type")])
op_phosphotase <- unique(op[op$enzyme_type == "phosphotase" & !is.na(op$enzyme_type), c("GENE", "enzyme_type")])
op_enzyme <- rbind(op_phosphotase, op_kinase[!(op_kinase$GENE %in% op_phosphotase$GENE),])
write.table(x = op_enzyme,
            file = paste0(makeOutDir(), "omnipath_enzyme.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")



# inputs ------------------------------------------------------------------
for (m.cutoff in 2:2) {
  for (NetworKIN.cutoff in 5:5) {
    ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$networkin_score >= NetworKIN.cutoff,]
    for (sig in c(0.2)) {
      for (cancer in c("BRCA", "OV", "CO", "UCEC")) {
        ## input KSEA scores
        df <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA_genomic_matched/OmniPath_NetworKIN_", 
                                   "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "/", cancer, "/KSEA Kinase Scores.csv"),
                    data.table = F)
        df <- df[df$p.value < 0.05 & df$m >= m.cutoff,]
        df$enzyme_direction <- ifelse(df$z.score > 0, "up", "down")
        df$method <- "KSEA"
        df$method_enzyme_direction <- paste0("KSEA-", df$enzyme_direction)
        
        ## input differential expression
        pho <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", cancer, "_PHO_diffexp_paired_FDR", sig, "_wilcoxon.txt"),
                                 data.table = F)
        pho <- cbind(pho, formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene))
        
        ## input differential expression for global phosphorylation
        phog <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", cancer, "_collapsed_PHO_diffexp_paired_FDR", sig, "_wilcoxon.txt"),
                     data.table = F)
        pho_kinase <- phog[phog$Gene %in% ptms_site_pairs_sup$GENE,]
        pho_kinase$log2FC <- log2(pho_kinase$`Fold Change`)
        
        ## merge KSEA and differential phosphorylation tgt
        df <- rbind(df[, c("Kinase.Gene", "enzyme_direction", "method_enzyme_direction", "log2FC", "method")], 
                    data.frame(Kinase.Gene = pho_kinase$Gene, 
                               enzyme_direction = pho_kinase$diffexp_type, 
                               method_enzyme_direction = paste0("diffexp-", pho_kinase$diffexp_type),
                               log2FC = pho_kinase$log2FC,
                               method = "diffexp"))
        colnames(df)[which(colnames(df) == "log2FC")] <- "enzyme_log2FC"
        
        ## add KSEA score into the enzyme-substrate table
        ptms_table <- merge(ptms_site_pairs_sup, df, by.x = c("GENE"), by.y = c("Kinase.Gene"), all.x = T)
        ptms_table$known_enzyme_phosphosite <- TRUE
        
        ## add in the regression result for just tumors
        reg_kinase <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/preprocess_files/tables/regression_omnipath&newworkin_protein_level/", "kinase", "_substrate_regression_cptac2p_", cancer, "_", "tumor", ".txt") , data.table = F)
        reg_phosphotase <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/preprocess_files/tables/regression_omnipath&newworkin_protein_level/", "phosphotase", "_substrate_regression_cptac2p_", cancer, "_", "tumor", ".txt") , data.table = F)
        reg_tab <- rbind(reg_kinase, reg_phosphotase)
        reg_cancer <- reg_tab[reg_tab$SELF=="trans" & reg_tab$Cancer == cancer, c("KINASE", "SUBSTRATE", "SUB_MOD_RSD", "FDR_pho_kin", "FDR_pro_sub", "coef_pho_kin", "coef_pro_sub", "Size", "P_pho_kin")]
        reg_cancer <- unique(reg_cancer)
        ptms_table <- merge(ptms_table, reg_cancer, by.x = c("GENE", "SUB_GENE", "SUB_MOD_RSD"), by.y = c("KINASE", "SUBSTRATE", "SUB_MOD_RSD"), all = T)

        ## add in the regression result for tumors + normals
        reg_kinase <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/preprocess_files/tables/regression_omnipath&newworkin_protein_level/", "kinase", "_substrate_regression_cptac2p_", cancer, "_", "normal+tumor", ".txt") , data.table = F)
        reg_phosphotase <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/preprocess_files/tables/regression_omnipath&newworkin_protein_level/", "phosphotase", "_substrate_regression_cptac2p_", cancer, "_", "normal+tumor", ".txt") , data.table = F)
        reg_tab <- rbind(reg_kinase, reg_phosphotase)
        reg_cancer <- reg_tab[reg_tab$SELF=="trans" & reg_tab$Cancer == cancer, c("KINASE", "SUBSTRATE", "SUB_MOD_RSD", "FDR_pho_kin", "FDR_pro_sub", "coef_pho_kin", "coef_pro_sub", "Size", "P_pho_kin")]
        reg_cancer <- unique(reg_cancer)
        ptms_table <- merge(ptms_table, reg_cancer, by.x = c("GENE", "SUB_GENE", "SUB_MOD_RSD"), by.y = c("KINASE", "SUBSTRATE", "SUB_MOD_RSD"), all = T, suffixes = c("",".tn"))
        
        ## add differential phosphorylation result
        ptms_table <- merge(ptms_table, pho[,c("SUBSTRATE", "SUB_MOD_RSD", "Fold Change")], by.x = c("SUB_GENE","SUB_MOD_RSD"), by.y = c("SUBSTRATE", "SUB_MOD_RSD"), all.x = T)
        colnames(ptms_table)[ncol(ptms_table)] <- c("substrate_diff_FC")
        
        ## add in the fold change for all substrates
        fc <- fread(input = paste0("/Users/yigewu/Box\ Sync/", "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/generate_table4KSEA/", cancer, "_phosphposite_FC_for_KSEA.csv"),
                    data.table = F)
        fc <- unique(fc[,c("Gene", "Residue.Both", "FC")])
        ptms_table <- merge(ptms_table, fc, by.x = c("SUB_GENE", "SUB_MOD_RSD"), by.y = c("Gene", "Residue.Both"), all.x = T)
        colnames(ptms_table)[ncol(ptms_table)] <- c("substrate_FC")
        
        ##
        ptms_table <- unique(ptms_table)
        
        ## annotate enzyme as kinase/phosphotase, etc
        ptms_table <- merge(ptms_table, op_enzyme, by = c("GENE"), all.x = T)
        ptms_table$enzyme_type[is.na(ptms_table$enzyme_type)] <- "kinase"
        
        ## annotate substrate differential phosphorylation
        tmp <- rep(c("NA"), nrow(ptms_table))
        tmp[is.na(ptms_table$substrate_FC)] <- "sub_FC_unavailable" 
        tmp[!is.na(ptms_table$substrate_FC)] <- paste0("FDR>=", sig)
        tmp[!is.na(ptms_table$substrate_diff_FC)] <- paste0("FDR<", sig)
        ptms_table$substrate_diffexp <- tmp
        
        ## annotate substrate FC directions
        tmp <- rep(c("NA"), nrow(ptms_table))
        tmp[!is.na(ptms_table$substrate_FC) & ptms_table$substrate_FC > 1] <- "up"
        tmp[!is.na(ptms_table$substrate_FC) & ptms_table$substrate_FC < 1] <- "down"
        ptms_table$substrate_direction <- tmp
        
        ptms_table$substrate_FC_diffexp <- paste0(ptms_table$substrate_direction, " (", ptms_table$substrate_diffexp, ")")
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
        ptms_table_consistent <- ptms_table[!is.na(ptms_table$consistent) & ptms_table$consistent,]
        write.table(x = ptms_table_consistent,
                    file = paste0(makeOutDir(), cancer, "_KSEA_enzyme_diffexp_substrate_consistent.txt"),
                    row.names = F, col.names = T, quote = F, sep = "\t")
        
        ## write out consistent enzyme-substrate pairs involvng
        ptms_table_driver <- ptms_table_consistent[((ptms_table_consistent$substrate.driver_type == "oncogene") | (ptms_table_consistent$enzyme.driver_type == "oncogene")) | ((ptms_table_consistent$substrate.driver_type == "tsg") | (ptms_table_consistent$enzyme.driver_type == "tsg")),]
        write.table(x = ptms_table_driver,
                    file = paste0(makeOutDir(), cancer, "_KSEA_enzyme_diffexp_substrate_consistent_w.driver.txt"),
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
