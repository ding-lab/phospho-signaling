# Yige Wu @ WashU 2018 Jan
# annotate the regression results in tumor samples and normal samples
# annotate the regression results with differential expression results
# annotate the regression results with genomic alteration correlation results

source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()


# overlap DE phosphosites and trans-regulated phosphosites ----------------
# Pho.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression/", "Pho.diffexp3can.txt"),
#                          data.table = F)
Pho.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "PHO_diffexp_3can_paired_FDR0.05.txt"),
                         data.table = F)
Pho.diffexp3can <- cbind(Pho.diffexp3can, formatPhosphosite(phosphosite_vector = Pho.diffexp3can$Phosphosite, Pho.diffexp3can$Gene))
Pho.diffexp3can <- Pho.diffexp3can[, !(colnames(Pho.diffexp3can) %in% c("Gene", "Phosphosite"))]
colnames(Pho.diffexp3can)[colnames(Pho.diffexp3can) == "Fold Change"] <- "FC_sub"

Phog.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "collapsed_PHO_diffexp_3can_paired_FDR0.05.txt"),
                          data.table = F)
Phog.diffexp3can$KINASE <- Phog.diffexp3can$Gene
Phog.diffexp3can <- Phog.diffexp3can[, !(colnames(Phog.diffexp3can) %in% c("Gene", "Phosphosite"))]
colnames(Phog.diffexp3can)[colnames(Phog.diffexp3can) == "Fold Change"] <- "FC_kin"

Pro.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "PRO_diffexp_3can_paired_FDR0.05.txt"),
                         data.table = F)
Pro.diffexp3can$KINASE <- Pro.diffexp3can$Gene
Pro.diffexp3can <- Pro.diffexp3can[, !(colnames(Pro.diffexp3can) %in% c("Gene", "Phosphosite"))]
colnames(Pro.diffexp3can)[colnames(Pro.diffexp3can) == "Fold Change"] <- "FC_kin"

## input mutational impact result
mutimpact <- fread(input = "PhosphoDrug/StatTests/12142017/correctedPvalue_merged", data.table = F)
mutimpact$id <- paste(mutimpact$Cancer, mutimpact$Mutated_Gene, mutimpact$Substrate_Gene, mutimpact$Phosphosite, sep = ":")
mutimpact <- mutimpact[order(-mutimpact$adjusted_pvalue),]
mutimpact <- mutimpact[!duplicated(mutimpact$id),]
tmp <- formatPhosphosite(phosphosite_vector = mutimpact$Phosphosite, gene_vector = mutimpact$Substrate_Gene)
names(tmp) <- c("SUB_MOD_RSD", "transcript", "SUBSTRATE")
tmp$SUB_MOD_RSD <- toupper(tmp$SUB_MOD_RSD)
mutimpact <- cbind(mutimpact, tmp)

#for( protein in c("kinase", "phosphotase")) {
for( protein in c("kinase")) {
  tn = paste0(resultD, "regression/tables/regression_cptac2p/", protein,"_substrate_regression_cptac2p_3can_", "tumor", ".txt")
  table_tumor <- fread(tn, data.table = F)
  table_tumor <- markSigSiteCan(regression = table_tumor, sig_thres = sig, protein_type = protein)
  
  tn = paste0(resultD, "regression/tables/regression_cptac2p/", protein,"_substrate_regression_cptac2p_3can_", "normal", ".txt")
  table_normal <- fread(tn, data.table = F)
  table_normal <- markSigSiteCan(regression = table_normal, sig_thres = sig, protein_type = protein)
  
  columns <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","FDR_pro_kin","FDR_pro_sub","FDR_pho_kin","coef_pro_kin",
               "coef_pro_sub","coef_pho_kin","Cancer","transcript","Size",
               "SELF","pair","fdr_sig","coef_sig" ,
               "sig_BRCA", "sig_OV","sig_CO","shared3can", "uniq_BRCA","uniq_OV","uniq_CO" )
  table_tn <- merge(table_tumor[,columns], table_normal[,columns], by = c("KINASE", "SUBSTRATE","SUB_MOD_RSD","Cancer","SELF","pair","transcript"),
                    all.x = T,  suffixes = c(".t", ".n"))
  
  ## add differential expression
  table_tn <- merge(table_tn, Pho.diffexp3can[,c("SUBSTRATE", "SUB_MOD_RSD", "Cancer", "transcript", "FC_sub")], by = c("SUBSTRATE", "SUB_MOD_RSD", "Cancer", "transcript"), all.x = T)
  table_tn_cis <- table_tn[table_tn$SELF == "cis",]
  table_tn_trans <- table_tn[table_tn$SELF == "trans",]
  
  table_tn_cis <- merge(table_tn_cis, Pro.diffexp3can[,c("KINASE", "Cancer", "FC_kin")], by = c("KINASE", "Cancer"), all.x = T)
  table_tn_trans <- merge(table_tn_trans, Phog.diffexp3can[,c("KINASE", "Cancer", "FC_kin")], by = c("KINASE", "Cancer"), all.x = T)
  table_tn <- rbind(table_tn_cis, table_tn_trans)
  
  table_tn$ks_diffexp <- ((table_tn$FC_kin>1) == (table_tn$FC_sub>1))
  table_tn$ks_diffexp_type <- "no_data+insig"
  table_tn$ks_diffexp_type[!is.na(table_tn$ks_diffexp) & (table_tn$ks_diffexp) == FALSE] <- "inconsistent"
  table_tn$ks_diffexp_type[table_tn$FC_kin>1 & table_tn$FC_sub>1] <- "up" 
  table_tn$ks_diffexp_type[table_tn$FC_kin<1 & table_tn$FC_sub<1] <- "down" 

  ## add genomic alteration info
  table_tn <- merge(table_tn, mutimpact[,c("Cancer", "Mutated_Gene", "SUBSTRATE", "SUB_MOD_RSD", "transcript",
                                       "Mutated_Sample_Size", "Non_Mutated_Sample_Size", "Fold_Change",
                                       "adjusted_pvalue")], 
                by.x = c("KINASE", "SUBSTRATE", "SUB_MOD_RSD","transcript", "Cancer"), 
                by.y = c("Mutated_Gene", "SUBSTRATE", "SUB_MOD_RSD", "transcript","Cancer"), all.x = T)

  ## reformat the table, fill in NA
  table_tn2p <- table_tn
  table_tn2p$Cancer <- factor(table_tn2p$Cancer, levels = c("BRCA", "OV", "CO"))
  table_tn2p$regulated.t <- as.character(table_tn2p$fdr_sig.t & table_tn2p$coef_sig.t)
  table_tn2p$regulated.n <- as.character(table_tn2p$fdr_sig.n & table_tn2p$coef_sig.n)
  table_tn2p$regulated.t[is.na(table_tn2p$regulated.t)] <- "no_data"
  table_tn2p$regulated.n[is.na(table_tn2p$regulated.n)] <- "no_data"
  table_tn2p$regulated.n <- factor(table_tn2p$regulated.n, levels = c("TRUE", "FALSE", "no_data"))
  table_tn2p$regulated.t <- factor(table_tn2p$regulated.t, levels = c("TRUE", "FALSE", "no_data"))
  table_tn2p$ks_diffexp_type <- factor(table_tn2p$ks_diffexp_type, levels = c("up", "down", "inconsistent" ,"no_data+insig"))
  
  ## tumor only classification
  table_tn2p$class.t <- "uncorrelated"
  table_tn2p$class.t[table_tn2p$regulated.t == "TRUE"] <- "regulated"
  table_tn2p$class.t[table_tn2p$regulated.t == "FALSE" & table_tn2p$ks_diffexp_type == "up"] <- "consistently high"
  table_tn2p$class.t[table_tn2p$regulated.t == "FALSE" & table_tn2p$ks_diffexp_type == "down"] <- "consistently low"
  
  ## tumor+normal classification
  table_tn2p$class.tn <- "unclassified"
  table_tn2p$class.tn[table_tn2p$regulated.t == "FALSE" & table_tn2p$regulated.n == "FALSE"] <- "uncorrelated"
  table_tn2p$class.tn[table_tn2p$ks_diffexp_type == "up"] <- "consistently high"
  table_tn2p$class.tn[table_tn2p$ks_diffexp_type == "down"] <- "consistently low"
  table_tn2p$class.tn[table_tn2p$regulated.t == "TRUE" & table_tn2p$regulated.n == "TRUE" ] <- "regulated"
  table_tn2p$class.tn[table_tn2p$regulated.t == "FALSE" & table_tn2p$regulated.n == "TRUE" ] <- "dis-regulated"
  table_tn2p$class.tn[table_tn2p$regulated.t == "TRUE" & table_tn2p$regulated.n == "FALSE" ] <- "neo-regulated"
  
  tn = paste0(resultDnow, protein,"_substrate_regression_cptac2p_3can_plus_normal_plus_diffexp_plus_mutimpact.txt")
  write.table(table_tn2p, file=tn, quote=F, sep = '\t', row.names = FALSE)
}

