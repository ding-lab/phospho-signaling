# Yige Wu @ WashU 2018 Apr
# integrate differential expression data with KSEA results

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# inputs ------------------------------------------------------------------
## input enzyme substrate table
ptms_site_pairs_sup <- read.csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv"))
op <- read.csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv"))
k_s_psp <- load_ks_table(protein = "kinase")
psp_cis <- k_s_psp[as.vector(k_s_psp$GENE) == as.vector(k_s_psp$SUB_GENE),]
psp_cis$pair <- paste0(psp_cis$GENE, ":", psp_cis$SUB_GENE, ':', psp_cis$SUB_MOD_RSD)

## annotate kinase and phosphotase
op_enzyme <- fread(input = "./cptac2p/analysis_results/phospho_network/diffexp/tables/integrate_enzyme_KSEAordiffexp_sub_FC_plus_regression/omnipath_enzyme.txt", data.table = F)

## annotate the amino acid specificity of kinase/phosphotase
ptms_site_pairs_sup$residue_type <- substr(ptms_site_pairs_sup$SUB_MOD_RSD, start = 1, stop = 1)
aa_tab <- unique(ptms_site_pairs_sup[ptms_site_pairs_sup$residue_type %in% c("S", "T", "Y"), c("residue_type", "GENE")])
aa_tab$aa_reported <- TRUE
# write.table(x = op_enzyme,
#             file = paste0(makeOutDir(), "omnipath_enzyme.txt"),
#             row.names = F, col.names = T, quote = F, sep = "\t")

## create output directory
dir.create(path = ppnD, "kinase_activity")
dir.create(path = ppnD, "kinase_activity/tables")

## input mutational impact result
# mutimpact <- fread(input = paste0("./PhosphoDrug/StatTests/12142017/correctedPvalue_merged"), data.table = F)
mutimpact <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/genoalt/tables/merge_mutation_impact/mutation_impact.txt"), data.table = F)
mutimpact$id <- paste(mutimpact$Cancer, mutimpact$Mutated_Gene, mutimpact$Substrate_Gene, mutimpact$Phosphosite, sep = ":")
mutimpact <- mutimpact[order(-mutimpact$p_value),]
mutimpact <- mutimpact[!duplicated(mutimpact$id),]
tmp <- formatPhosphosite(phosphosite_vector = mutimpact$Phosphosite, gene_vector = mutimpact$Substrate_Gene)
names(tmp) <- c("SUB_MOD_RSD", "transcript", "SUBSTRATE")
tmp$SUB_MOD_RSD <- toupper(tmp$SUB_MOD_RSD)
mutimpact <- cbind(mutimpact, tmp)
mutimpact <- mutimpact[!is.na(mutimpact$SUB_MOD_RSD),]

# meat ------------------------------------------------------------------
for (m.cutoff in 2:2) {
  for (NetworKIN.cutoff in 5:5) {
    for (sig in c(0.2)) {
      ptms_3can_aa_reported <- NULL
      for (cancer in c("BRCA", "OV", "CO")) {
      # for (cancer in c("OV", "CO")) {
        ptms_table <- ptms_site_pairs_sup[ptms_site_pairs_sup$networkin_score >= NetworKIN.cutoff,]
        ptms_table$known_enzyme_phosphosite <- TRUE
        
        ## add in the regression result for just tumors
        # reg_kinase <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/preprocess_files/tables/regression_omnipath&newworkin_protein_level/", "kinase", "_substrate_regression_cptac2p_", cancer, "_", "tumor", ".txt") , data.table = F)
        reg_kinase <- fread(input = paste0(preprocess_filesD, "tables/regression_omnipath&newworkin_protein_level/", "kinase", "_substrate_regression_cptac2p_", cancer, "_", "tumor", ".txt"),
                            data.table = F)
        
        # reg_phosphotase <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/preprocess_files/tables/regression_omnipath&newworkin_protein_level/", "phosphotase", "_substrate_regression_cptac2p_", cancer, "_", "tumor", ".txt") , data.table = F)
        reg_phosphotase <- fread(input = paste0(preprocess_filesD, "tables/regression_omnipath&newworkin_protein_level/", "phosphotase", "_substrate_regression_cptac2p_", cancer, "_", "tumor", ".txt"),
                                 data.table = F)
        
        reg_tab <- rbind(reg_kinase, reg_phosphotase)
        reg_cancer <- reg_tab[reg_tab$Cancer == cancer, c("KINASE", "SUBSTRATE", "SUB_MOD_RSD", 
                                                                                  "FDR_pho_kin", "FDR_pro_sub", "FDR_pro_kin", 
                                                                                  "coef_pho_kin", "coef_pro_sub", "coef_pro_kin", 
                                                                                  "Size", "P_pho_kin", "P_pro_kin", "SELF")]
        reg_cancer <- unique(reg_cancer)
        ptms_table <- merge(ptms_table, reg_cancer, by.x = c("GENE", "SUB_GENE", "SUB_MOD_RSD"), by.y = c("KINASE", "SUBSTRATE", "SUB_MOD_RSD"), all = T)
        
        ## add in the regression result for tumors + normals
        # reg_kinase <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/preprocess_files/tables/regression_omnipath&newworkin_protein_level/", "kinase", "_substrate_regression_cptac2p_", cancer, "_", "normal+tumor", ".txt") , data.table = F)
        reg_kinase <- fread(input = paste0(preprocess_filesD, "tables/regression_omnipath&newworkin_protein_level/", "kinase", "_substrate_regression_cptac2p_", cancer, "_", "normal+tumor", ".txt") , data.table = F)
        
        # reg_phosphotase <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/preprocess_files/tables/regression_omnipath&newworkin_protein_level/", "phosphotase", "_substrate_regression_cptac2p_", cancer, "_", "normal+tumor", ".txt") , data.table = F)
        reg_phosphotase <- fread(input = paste0(preprocess_filesD, "tables/regression_omnipath&newworkin_protein_level/", "phosphotase", "_substrate_regression_cptac2p_", cancer, "_", "normal+tumor", ".txt") , data.table = F)
        
        reg_tab <- rbind(reg_kinase, reg_phosphotase)
        reg_cancer <- reg_tab[reg_tab$Cancer == cancer, c("KINASE", "SUBSTRATE", "SUB_MOD_RSD", 
                                                                                  "FDR_pho_kin", "FDR_pro_sub", "FDR_pro_kin", 
                                                                                  "coef_pho_kin", "coef_pro_sub", "coef_pro_kin", 
                                                                                  "Size", "P_pho_kin", "P_pro_kin")]
        reg_cancer <- unique(reg_cancer)
        ptms_table <- merge(ptms_table, reg_cancer, by.x = c("GENE", "SUB_GENE", "SUB_MOD_RSD"), by.y = c("KINASE", "SUBSTRATE", "SUB_MOD_RSD"), all = T, suffixes = c("",".tn"))
        ptms_table$known_enzyme_phosphosite[is.na(ptms_table$known_enzyme_phosphosite)] <- FALSE
        ptms_table$SELF[is.na(ptms_table$SELF)] <- "trans"
        
        ## input KSEA scores
        df <- read.csv(paste0(resultD, "phospho_network/classify/tables/runKSEA_genomic_matched/OmniPath_NetworKIN_", 
                                   "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "/", cancer, "/KSEA Kinase Scores.csv"))
        df <- df[df$m >= m.cutoff,]
        df$KSEA_enzyme_direction <- ifelse(df$z.score > 0, "up", "down")
        df$KSEA_log2FC <- df$log2FC
        df$KSEA_pvalue <- df$p.value

        ## input differential expression
        # pho <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", cancer, "_PHO_diffexp_paired_FDR", sig, "_wilcoxon.txt"),
        #                          data.table = F)
        pho <- read.table(paste0(resultD, "phospho_network/diffexp/tables/differential_expression_paired/", cancer, "_PHO_diffexp_paired_FDR", sig, "_wilcoxon.txt"),
                     header = T, sep = "\t")
        pho <- cbind(pho, formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene))
        
        ## input differential expression for global phosphorylation
        # phog <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", cancer, "_collapsed_PHO_diffexp_paired_FDR", sig, "_wilcoxon.txt"),
        #              data.table = F)
        phog <- read.table(paste0(resultD, "phospho_network/diffexp/tables/differential_expression_paired/", cancer, "_collapsed_PHO_diffexp_paired_FDR", sig, "_wilcoxon.txt"),
                           header = T, sep = "\t")
        pho_kinase <- phog[phog$Gene %in% ptms_site_pairs_sup$GENE,]
        pho_kinase$log2FC <- log2(pho_kinase[, colnames(pho_kinase)[grepl(pattern = "Fold", x = colnames(pho_kinase))]])
        
        ## merge KSEA and differential phosphorylation tgt
        df <- merge(df[, c("Kinase.Gene", "KSEA_enzyme_direction", "KSEA_log2FC", "KSEA_pvalue")], 
                    data.frame(Kinase.Gene = pho_kinase$Gene, 
                               diffexp_enzyme_direction = pho_kinase$diffexp_type,
                               diffexp_log2FC = pho_kinase$log2FC), all = T)
        ## annotate enzyme phospho fold change direction
        df$enzyme_direction <- NA
        df$enzyme_direction[!is.na(df$KSEA_enzyme_direction) & df$KSEA_pvalue < sig] <- df$KSEA_enzyme_direction[!is.na(df$KSEA_enzyme_direction) & df$KSEA_pvalue < sig]
        df$enzyme_direction[is.na(df$KSEA_enzyme_direction)] <- as.vector(df$diffexp_enzyme_direction)[is.na(df$KSEA_enzyme_direction)]
        
        ## add enzyme fold changes into the enzyme-substrate table
        ptms_table <- merge(ptms_table, df, by.x = c("GENE"), by.y = c("Kinase.Gene"), all.x = T)
        
        ## add in the fold change for all substrates
        # fc <- fread(input = paste0("/Users/yigewu/Box\ Sync/", "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/generate_table4KSEA/", cancer, "_phosphposite_FC_for_KSEA.csv"),
        #             data.table = F)
        fc <- read.table(paste0(resultD, "phospho_network/classify/tables/generate_table4KSEA/", cancer, "_phosphposite_FC_for_KSEA.csv"),
                         header = T, sep = ",")
        fc <- unique(fc[,c("Gene", "Residue.Both", "FC")])
        fc$log2FC <- log2(fc$FC)
        ptms_table <- merge(ptms_table, fc[,c("Gene", "Residue.Both", "log2FC")], by.x = c("SUB_GENE", "SUB_MOD_RSD"), by.y = c("Gene", "Residue.Both"), all.x = T)
        colnames(ptms_table)[ncol(ptms_table)] <- c("substrate_log2FC")
        
        ## annotate substrate phoshp fold change direction
        ptms_table <- unique(ptms_table)
        ptms_table$substrate_direction <- NA
        ptms_table$substrate_direction[ptms_table$substrate_log2FC < (-1)] <- "down"
        ptms_table$substrate_direction[ptms_table$substrate_log2FC > 1] <- "up"
        
        ## annotate enzyme as kinase/phosphotase, etc
        ptms_table <- merge(ptms_table, op_enzyme, by = c("GENE"), all.x = T)
        ptms_table$enzyme_type[is.na(ptms_table$enzyme_type)] <- "kinase"
        
        ## annotate the consistency of phospho FC betwen enzyme and substrate
        ptms_table$same_direction <- (as.vector(ptms_table$enzyme_direction) == as.vector(ptms_table$substrate_direction))
        ptms_table$consistent <- ifelse((ptms_table$enzyme_type == "kinase" & ptms_table$same_direction) | (ptms_table$enzyme_type == "phosphotase" & !ptms_table$same_direction), TRUE, FALSE)
        
        ## annotate mutation impact
        mutimpact_can <- mutimpact[mutimpact$Cancer == cancer,]
        mutimpact_can_pro <- mutimpact_can[mutimpact_can$SUB_MOD_RSD == "PROTEIN",]
        mutimpact_can_pho <- mutimpact_can[mutimpact_can$SUB_MOD_RSD != "PROTEIN",]
        
        ptms_table <- merge(ptms_table, mutimpact_can_pho[,c("Mutated_Gene", "SUBSTRATE", "SUB_MOD_RSD", 
                                                       "Mutated_Sample_Size", "Non_Mutated_Sample_Size", 
                                                       "Fold_Change", "p_value")], 
                      by.x = c("GENE", "SUB_GENE", "SUB_MOD_RSD"),
                      by.y = c("Mutated_Gene", "SUBSTRATE", "SUB_MOD_RSD"), all.x = T)
        ptms_table <- merge(ptms_table, mutimpact_can_pro[,c("Mutated_Gene", "SUBSTRATE", 
                                                             "Mutated_Sample_Size", "Non_Mutated_Sample_Size", 
                                                             "Fold_Change", "p_value")], 
                            by.x = c("GENE", "SUB_GENE"),
                            by.y = c("Mutated_Gene", "SUBSTRATE"), all.x = T, suffixes = c("", ".pro"))
        
        
        ## annoate drivers
        tsg <- loadGeneList(gene_type = "tsg", cancer = cancer, is.soft.limit = "soft")
        oncogene <- loadGeneList(gene_type = "oncogene", cancer = cancer, is.soft.limit = "soft")
        ptms_table$enzyme.driver_type <- "notdriver"
        ptms_table$substrate.driver_type <- "notdriver"
        
        ptms_table$enzyme.driver_type[ptms_table$GENE %in% tsg] <- "tsg"
        ptms_table$enzyme.driver_type[ptms_table$GENE %in% oncogene] <- "oncogene"
        ptms_table$substrate.driver_type[ptms_table$SUB_GENE %in% tsg] <- "tsg"
        ptms_table$substrate.driver_type[ptms_table$SUB_GENE %in% oncogene] <- "oncogene"
        
        ## annotate whether the ptm by this enzyme to each type amino acid has been reported
        ptms_table$residue_type <- substr(ptms_table$SUB_MOD_RSD, start = 1, stop = 1)
        ptms_table <- merge(ptms_table, aa_tab, all.x = T)
        ptms_table$aa_reported[is.na(ptms_table$aa_reported)] <- FALSE
        
        ## re-adjust the FDR within the pairs with reported amino acid ptm
        ptms_table$FDR_pho_kin.aa_reported <- NA
        ptms_table$FDR_pro_kin.aa_reported <- NA
        tmp <- !is.na(ptms_table$P_pho_kin) & ptms_table$aa_reported & ptms_table$enzyme_type == "kinase"
        ptms_table$FDR_pho_kin.aa_reported[tmp] <- p.adjust(ptms_table$P_pho_kin[tmp], method = "fdr")
        tmp <- !is.na(ptms_table$P_pho_kin) & ptms_table$aa_reported & ptms_table$enzyme_type == "phosphotase"
        ptms_table$FDR_pho_kin.aa_reported[tmp] <- p.adjust(ptms_table$P_pho_kin[tmp], method = "fdr")
        tmp <- !is.na(ptms_table$P_pro_kin) & ptms_table$aa_reported & ptms_table$enzyme_type == "kinase"
        ptms_table$FDR_pro_kin.aa_reported[tmp] <- p.adjust(ptms_table$P_pro_kin[tmp], method = "fdr")
        
        ptms_table$pair <- paste0(ptms_table$GENE, ":", ptms_table$SUB_GENE, ':', ptms_table$SUB_MOD_RSD)
        ptms_table$known_enzyme_phosphosite[ptms_table$pair %in% psp_cis$pair] <- TRUE
        
        write.table(x = ptms_table,
                    file = paste0(ppnD, "/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/",
                                  cancer, "_KSEA_enzyme_diffexp_substrate.txt"),
                    row.names = F, col.names = T, quote = F, sep = "\t")
        
        ptms_table_aa_reported <- ptms_table[ptms_table$aa_reported,]
        ptms_table_aa_reported$Cancer <- cancer
        ptms_3can_aa_reported <- rbind(ptms_3can_aa_reported, ptms_table_aa_reported)
        write.table(x = ptms_table_aa_reported,
                    file = paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/",
                                  cancer, "_KSEA_enzyme_diffexp_substrate_aa_reported.txt"),
                    row.names = F, col.names = T, quote = F, sep = "\t")

        ## write out consistent directed enzyme-substrate pairs
        ptms_table_consistent <- ptms_table_aa_reported[!is.na(ptms_table$consistent) & ptms_table$consistent,]
        write.table(x = ptms_table_consistent,
                    file = paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/",
                                  cancer, "_KSEA_enzyme_diffexp_substrate_consistent.txt"),
                    row.names = F, col.names = T, quote = F, sep = "\t")

        ## write out consistent enzyme-substrate pairs involvng
        ptms_table_driver <- ptms_table_consistent[((ptms_table_consistent$substrate.driver_type == "oncogene") | (ptms_table_consistent$enzyme.driver_type == "oncogene")) | ((ptms_table_consistent$substrate.driver_type == "tsg") | (ptms_table_consistent$enzyme.driver_type == "tsg")),]
        write.table(x = ptms_table_driver,
                    file = paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/",
                                  cancer, "_KSEA_enzyme_diffexp_substrate_consistent_w.driver.txt"),
                    row.names = F, col.names = T, quote = F, sep = "\t")

        ## write out consistent enzyme-substrate pairs and also in concord with its driver role
        ptms_table_consistent_driver <- ptms_table_consistent[((ptms_table_consistent$substrate_direction == "up" & ptms_table_consistent$substrate.driver_type == "oncogene") | (ptms_table_consistent$enzyme_direction == "up" & ptms_table_consistent$enzyme.driver_type == "oncogene")) | ((ptms_table_consistent$substrate_direction == "down" & ptms_table_consistent$substrate.driver_type == "tsg") | (ptms_table_consistent$enzyme_direction == "down" & ptms_table_consistent$enzyme.driver_type == "tsg")),]
        write.table(x = ptms_table_consistent_driver,
                    file = paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/",
                                  cancer, "_KSEA_enzyme_diffexp_substrate_consistent_driver.txt"),
                    row.names = F, col.names = T, quote = F, sep = "\t")
        
        print(cancer)
      }
    }
  }
}
write.table(x = ptms_3can_aa_reported,
            file = paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                          "3can", "_KSEA_enzyme_diffexp_substrate_aa_reported.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")

