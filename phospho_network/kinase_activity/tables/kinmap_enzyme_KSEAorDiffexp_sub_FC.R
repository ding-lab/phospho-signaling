# Yige Wu @ WashU 2018 Mar
# barplot displaying KSEA score for significantly enriched kinase/phosphotase

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# test --------------------------------------------------------------------
## input omnipath to annotate kinase/phosphotase
op <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_genename_annotated.csv")
op$enzyme_type <- NA
op$enzyme_type[op$modification == "phosphorylation"] <- "kinase"
op$enzyme_type[op$modification == "dephosphorylation"] <- "phosphotase"
op_kinase <- unique(op[op$enzyme_type == "kinase" & !is.na(op$enzyme_type), c("GENE", "enzyme_type")])
op_phosphotase <- unique(op[op$enzyme_type == "phosphotase" & !is.na(op$enzyme_type), c("GENE", "enzyme_type")])
op_enzyme <- rbind(op_phosphotase, op_kinase[!(op_kinase$GENE %in% op_phosphotase$GENE),])

## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv")

for (m.cutoff in 2:2) {
  for (NetworKIN.cutoff in 5:5) {
    ptms_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv")
    ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$networkin_score >= NetworKIN.cutoff,]
    for (cancer in c("BRCA", "OV", "CO")) {
    # for (cancer in c("BRCA")) {
      for (reg_sig in c(0.05)) {
        ## get KSEA scores integrated with differential phosphorylation results
        ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/diffexp/tables/integrate_enzyme_KSEAordiffexp_sub_FC_plus_regression/", 
                                             cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
        ksea_diffexp$reg_issig <- (ksea_diffexp$FDR_pho_kin < reg_sig)
        ksea_diffexp_direction <- data.frame(table(ksea_diffexp[, c("GENE", "substrate_FC_diffexp")]))
        ksea_diffexp_reg <- data.frame(table(ksea_diffexp[, c("GENE", "substrate_FC_diffexp", "reg_issig")]))
        
        ## get a list of enzymes with KSEA or differential expression data
        df <- unique(ksea_diffexp[!is.na(ksea_diffexp$method), c("GENE", "method", colnames(ksea_diffexp)[grepl(pattern = "enzyme", x = colnames(ksea_diffexp))])])
        
        df <- merge(df, ksea_diffexp_direction, by = c("GENE"), all.x = T)
        colnames(df)[1] <- "Kinase.Gene" 
        
        df$Freq2p <- df$Freq; df$Freq2p[grepl(x = df$substrate_FC_diffexp, pattern = "down")] <- (-df$Freq2p[grepl(x = df$substrate_FC_diffexp, pattern = "down")])
        df$Freq2p <- as.numeric(df$Freq2p)
        df$substrate_direction <- ifelse(grepl(pattern = "up", x = df$substrate_FC_diffexp), "up", "down")
        df <- unique(df)
        
        for (enzyme_type in c("kinase")) {
          df_e <- df[df$Kinase.Gene %in% op_enzyme$GENE[op_enzyme$enzyme_type == enzyme_type],]
          df_e$Kinase.Gene <- factor(df_e$Kinase.Gene, levels = as.vector(df_e$Kinase.Gene)[order(as.vector(abs(df_e$enzyme_log2FC)))])
          df_e$enzyme_direction <- factor(df_e$enzyme_direction, levels = c("up", "down"))
          
          df_e_kin <- unique(df_e[, c("Kinase.Gene", "method", colnames(df_e)[grepl(pattern = "enzyme", x = colnames(df_e))])])
          df_e_kin$method <- factor(df_e_kin$method, levels = c("diffexp", "KSEA"))
          
          for (enzyme_direction in c("up")) {
            df_e_kin_d <- df_e_kin[df_e_kin$enzyme_direction == enzyme_direction,]
            df_e_kin_d$method <- factor(df_e_kin_d$method, levels = c("diffexp", "KSEA"))
            df_e_kin_d <- df_e_kin_d[!duplicated(df_e_kin_d$Kinase.Gene),]
            tn = paste(makeOutDir(), cancer , "_", enzyme_type, '_', enzyme_direction, '.txt',sep ="")
            sink(file = tn)
            for (i in 1:nrow(df_e_kin_d)) {
              enzyme <- as.character(df_e_kin_d[i, "Kinase.Gene"])
              method <- as.character(df_e_kin_d[i, "method"])
              if (method != "diffexp") {
                cat(paste0("@ 3 : ", 20*(df_e_kin_d[i, "enzyme_log2FC"]), " : red : red : 3\n"))
              } else {
                cat(paste0("@ 3 : ", 20*(df_e_kin_d[i, "enzyme_log2FC"]), " : red : yellow : 3\n"))
              }
              cat(paste0(enzyme, "\n\n"))
            }
            sink()
            closeAllConnections()
          }
          for (enzyme_direction in c("down")) {
            df_e_kin_d <- df_e_kin[df_e_kin$enzyme_direction == enzyme_direction,]
            df_e_kin_d$method <- factor(df_e_kin_d$method, levels = c("diffexp", "KSEA"))
            df_e_kin_d <- df_e_kin_d[!duplicated(df_e_kin_d$Kinase.Gene),]
            tn = paste(makeOutDir(), cancer , "_", enzyme_type, '_', enzyme_direction, '.txt',sep ="")
            sink(file = tn)
            for (i in 1:nrow(df_e_kin_d)) {
              enzyme <- as.character(df_e_kin_d[i, "Kinase.Gene"])
              method <- as.character(df_e_kin_d[i, "method"])
              if (method != "diffexp") {
                cat(paste0("@ 3 : ", 20*(df_e_kin_d[i, "enzyme_log2FC"]), " : blue : blue : 3\n"))
              } else {
                cat(paste0("@ 3 : ", 20*(df_e_kin_d[i, "enzyme_log2FC"]), " : blue : yellow : 3\n"))
              }
              cat(paste0(enzyme, "\n\n"))
            }
            sink()
            closeAllConnections()
          }
        }
        
        for (enzyme_type in c("phosphotase")) {
          df_e <- df[df$Kinase.Gene %in% op_enzyme$GENE[op_enzyme$enzyme_type == enzyme_type],]
          df_e$enzyme_log2FC[df_e$method == "KSEA"] <- -(df_e$enzyme_log2FC[df_e$method == "KSEA"])
          df_e$Kinase.Gene <- factor(df_e$Kinase.Gene, levels = as.vector(df_e$Kinase.Gene)[order(as.vector(abs(df_e$enzyme_log2FC)))])
          df_e$enzyme_direction <- factor(df_e$enzyme_direction, levels = c("up", "down"))
          
          df_e_kin <- unique(df_e[, c("Kinase.Gene", "method", colnames(df_e)[grepl(pattern = "enzyme", x = colnames(df_e))])])
          df_e_kin$method <- factor(df_e_kin$method, levels = c("diffexp", "KSEA"))
        }
      }
    }
  }
}
