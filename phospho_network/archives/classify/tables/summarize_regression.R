# Yige Wu @ WashU 2018 Jan
# summarize the prevalence of up-regulated and down-regulated ks pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# set variables -----------------------------------------------------------
datatypes <- c("PRO", "PHO", "collapsed_PHO")
fdr_thres <- c(0.1, 0.05)
names(fdr_thres) <- c("phosphatase", "kinase")

# inputs ------------------------------------------------------------------
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))
ptms_site_pairs_sup <- read_csv(paste0("./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv"))

# get actual genes with summary numbers-----------------------------------------------------
for (enzyme_type in c("phosphatase")) {
  regression <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/phosphatase_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
  regression$pair <- paste0(regression$KINASE, ":", regression$SUBSTRATE, ":", regression$SUB_MOD_RSD)
  length(unique(regression$pair))
  fdr_thres_en <- fdr_thres[enzyme_type] 
  
  es_table <- ptms_site_pairs_sup[ptms_site_pairs_sup$enzyme_type == enzyme_type,]
  es_table$pair_pro <- paste0(es_table$GENE, ":", es_table$SUB_GENE)
  for (sig in c(1)) {
    sum_table <- NULL

    for (cancer in c("BRCA", "OV", "CO")) {
      sink(paste0(makeOutDir(resultD = resultD), cancer,"_", enzyme_type, '_substrate_pairs_summary_log_fdr', fdr_thres_en, ".txt"), append=FALSE, split=FALSE)
      cat(paste0("-", cancer, " summary:\n"))
      
      cat(paste0("original known relationships: ", length(unique(es_table$pair[es_table$SUB_MOD_RSD != "NULL"])), "(site-lvel), ", length(unique(es_table$pair_pro)), "(protien-level), ", ":\n"))
      ## summarize protein/phosphoprotein number input to differential expression analysis (all, kinase, substrate)
      phog <- fread(paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_", "collapsed_PHO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                    data.table = F)
      tmp <- data.frame(datatype = "collapsed_PHO",
                        cancer = cancer,
                        diffexp_type = "all",
                        gene_cat = c(enzyme_type, paste0("not-", enzyme_type)),
                        num_gene = as.matrix(table(phog$Gene %in% es_table$GENE))[c("TRUE", "FALSE"),1])
      tmp$num_site <- tmp$num_gene
      sum_table <- rbind(sum_table, tmp)
      
      pho <- fread(paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                   data.table = F)
      tmp <- data.frame(datatype = "PHO",
                        cancer = cancer,
                        diffexp_type = "all",
                        gene_cat = c("substrate", "not-substrate"),
                        num_gene = as.matrix(table(unique(pho$Gene) %in% es_table$SUB_GENE))[c("TRUE", "FALSE"),1],
                        num_site = as.matrix(table(pho$Gene %in% es_table$SUB_GENE))[c("TRUE", "FALSE"),1])
      sum_table <- rbind(sum_table, tmp)
      
      pro <- fread(paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                   data.table = F)
      tmp <- data.frame(datatype = "PRO",
                        cancer = cancer,
                        diffexp_type = "all",
                        gene_cat = c(enzyme_type, paste0("not-", enzyme_type)),
                        num_gene = as.matrix(table(pro$Gene %in% es_table$GENE))[c("TRUE", "FALSE"),1])
      tmp$num_site <- tmp$num_gene
      sum_table <- rbind(sum_table, tmp)
      
      sink()
      closeAllConnections()
    }
  }
}


# measured in all 3 cancers -----------------------------------------------
enzyme_type <- "kinase"
sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin_protein_level/", enzyme_type, "_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
sup_cans_tab_en$pair <- paste0(sup_cans_tab_en$KINASE, ":", sup_cans_tab_en$SUBSTRATE, ":", sup_cans_tab_en$SUB_MOD_RSD)
colnames(sup_cans_tab_en)[1:2] <- c("GENE", "SUB_GENE")
pair_cans_count  <- data.frame(table(sup_cans_tab_en$pair))
pair_cans_count3 <- pair_cans_count[pair_cans_count$Freq == length(cancers_sort),]
nrow(pair_cans_count3)


fdr_pk <- 0.05
sup_cans_tab_en <- markSigSiteCan(sup_cans_tab_en, sig_thres = fdr_pk, enzyme_type = enzyme_type)
## annotate kinase substrate regulation
sup_cans_tab_en$regulated <- (sup_cans_tab_en$coef_sig & sup_cans_tab_en$fdr_sig)
length(unique(sup_cans_tab_en$pair[sup_cans_tab_en$uniq_BRCA]))
