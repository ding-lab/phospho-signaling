# Yige Wu @ WashU 2018 Jan
# integrate differential expression data with regression results and mutational impact results

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
## input mutational impact result
mutimpact <- fread(input = "PhosphoDrug/StatTests/12142017/correctedPvalue_merged", data.table = F)
mutimpact$id <- paste(mutimpact$Cancer, mutimpact$Mutated_Gene, mutimpact$Substrate_Gene, mutimpact$Phosphosite, sep = ":")
mutimpact <- mutimpact[order(-mutimpact$adjusted_pvalue),]
mutimpact <- mutimpact[!duplicated(mutimpact$id),]
tmp <- formatPhosphosite(phosphosite_vector = mutimpact$Phosphosite, gene_vector = mutimpact$Substrate_Gene)
names(tmp) <- c("SUB_MOD_RSD", "transcript", "SUBSTRATE")
tmp$SUB_MOD_RSD <- toupper(tmp$SUB_MOD_RSD)
mutimpact <- cbind(mutimpact, tmp)
mutimpact <- mutimpact[!is.na(mutimpact$SUB_MOD_RSD),]

## input kinase-substrate table and phosphotase-substrate table
ks_table <- load_ks_table("kinase"); ks_table$upstream <- "kinase"
ps_table <- load_ks_table("phosphotase"); ps_table$upstream <- "phosphotase"
us_table <- unique(rbind(ks_table[,c("GENE", "SUB_GENE", "upstream")], ps_table[,c("GENE", "SUB_GENE", "upstream")]))
colnames(us_table) <- c("KINASE", "SUBSTRATE", "upstream")
us_table$SELF <- ifelse(us_table$KINASE == as.vector(us_table$SUBSTRATE), "cis", "trans")
us_table_trans <- us_table[us_table$SELF == "trans",]
us_table_cis <- us_table[us_table$SELF == "cis",]
### some consistently high/low ks pairs might be lost in regression analysis because the later requires at least 5 samples with both kinase and substrate data
### but we probably wouldn't want those pairs if they didn't even fullfill this requirement
### only use pairs when both the upstream and the substrate is both available in one cancer type

# merge based on regression table -------------------------------------------------------------------
## input regression results
tn = paste0(resultD,"regression/tables/kinase_substrate_regression_cptac2p_3can.txt")
table_kinase <- fread(tn, data.table = F)
table_kinase$upstream <- "kinase"
tn = paste0(resultD,"regression/tables/phosphotase_substrate_regression_cptac2p_3can.txt")
table_phosphotase <- fread(tn, data.table = F)
table_phosphotase$upstream <- "phosphotase"
regression <- rbind(table_kinase, table_phosphotase)

for (sig in c(0.1)) {
  ## input differential expression results
  Pho.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "PHO_diffexp_3can_paired_FDR", sig, ".txt"),
                           data.table = F)
  Pho.diffexp3can <- cbind(Pho.diffexp3can, formatPhosphosite(phosphosite_vector = Pho.diffexp3can$Phosphosite, Pho.diffexp3can$Gene))
  Pho.diffexp3can <- Pho.diffexp3can[, !(colnames(Pho.diffexp3can) %in% c("Gene", "Phosphosite"))]
  colnames(Pho.diffexp3can)[colnames(Pho.diffexp3can) == "Fold Change"] <- "FC_sub"
  
  Phog.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "collapsed_PHO_diffexp_3can_paired_FDR", sig, ".txt"),
                            data.table = F)
  Phog.diffexp3can$KINASE <- Phog.diffexp3can$Gene
  Phog.diffexp3can <- Phog.diffexp3can[, !(colnames(Phog.diffexp3can) %in% c("Gene", "Phosphosite"))]
  colnames(Phog.diffexp3can)[colnames(Phog.diffexp3can) == "Fold Change"] <- "FC_kin"
  
  Pro.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "PRO_diffexp_3can_paired_FDR", sig, ".txt"),
                           data.table = F)
  Pro.diffexp3can$KINASE <- Pro.diffexp3can$Gene
  Pro.diffexp3can <- Pro.diffexp3can[, !(colnames(Pro.diffexp3can) %in% c("Gene", "Phosphosite"))]
  colnames(Pro.diffexp3can)[colnames(Pro.diffexp3can) == "Fold Change"] <- "FC_kin"
  
  ## merge based on the initial kinase/phosphotase-substrate table
  us_diffexp <- merge(us_table, 
                      Pho.diffexp3can[, c("FC_sub", "SUBSTRATE", "SUB_MOD_RSD", "Cancer", "transcript")], 
                      by = c("SUBSTRATE"))
  us_trans <- us_diffexp[us_diffexp$SELF == "trans",]
  us_cis <- us_diffexp[us_diffexp$SELF == "cis",]
  us_cis <- merge(us_cis, Pro.diffexp3can[,c("KINASE", "Cancer", "FC_kin")], by = c("KINASE", "Cancer"), all.x = T)
  us_trans <- merge(us_trans, Phog.diffexp3can[,c("KINASE", "Cancer", "FC_kin")], by = c("KINASE", "Cancer"), all.x = T)
  us_diffexp <- rbind(us_cis, us_trans)
  us_diffexp$upstream_direction <- sapply(1:nrow(us_diffexp), function(n) ifelse(is.na(us_diffexp$FC_kin[n]), NA, ifelse(us_diffexp$FC_kin[n] > 1, "up", "down")))
  us_diffexp$downstream_direction <- sapply(1:nrow(us_diffexp), function(n) ifelse(is.na(us_diffexp$FC_sub[n]), NA, ifelse(us_diffexp$FC_sub[n] > 1, "up", "down")))
  resultDnow <- makeOutDir()
  write.table(us_diffexp, row.names = F, quote=F, sep = '\t', file=paste0(resultDnow, "kinase_phosphotase_substrate_diffexp.txt"))
  
  
  ## merge regression with substrate phosphosite differential expression
  regression_diffexp <- merge(regression, 
                              Pho.diffexp3can[, c("FC_sub", "SUBSTRATE", "SUB_MOD_RSD", "Cancer", "transcript")], 
                              by = c("SUBSTRATE", "SUB_MOD_RSD", "Cancer", "transcript"), all.x = T)
  
  ## merge regression with kinase differential expression
  regression_cis <- regression_diffexp[!is.na(regression_diffexp$SELF) & regression_diffexp$SELF == "cis",]
  regression_cis <- merge(regression_cis, Pro.diffexp3can[,c("KINASE", "Cancer", "FC_kin")], by = c("KINASE", "Cancer"), all.x = T)
  regression_trans <- regression_diffexp[!is.na(regression_diffexp$SELF) & regression_diffexp$SELF == "trans",]
  regression_trans <- merge(regression_trans, Phog.diffexp3can[,c("KINASE", "Cancer", "FC_kin")], by = c("KINASE", "Cancer"), all.x = T)
  ## annotate FC_kin and FC_sub with whether data is available
  regression_diffexp <- NULL
  for (cancer in c("BRCA", "OV", "CO")) {
    pro <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                 data.table = F)
    phog <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "collapsed_PHO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                  data.table = F)

    ## annotate kinases in cis pairs with data availability
    regression_cis_can <- regression_cis[regression_cis$Cancer == cancer,]
    regression_cis_can_na <- regression_cis_can[is.na(regression_cis_can$FC_kin),]
    regression_cis_can_nonna <- regression_cis_can[!is.na(regression_cis_can$FC_kin),]
    regression_cis_can_na$FC_kin <- ifelse(regression_cis_can_na$KINASE %in% pro$Gene, 1, NA)
    
    ## annotate kinases in cis pairs with data availability
    regression_trans_can <- regression_trans[regression_trans$Cancer == cancer,]
    regression_trans_can_na <- regression_trans_can[is.na(regression_trans_can$FC_kin),]
    regression_trans_can_nonna <- regression_trans_can[!is.na(regression_trans_can$FC_kin),]
    regression_trans_can_na$FC_kin <- ifelse(regression_trans_can_na$KINASE %in% phog$Gene, 1, NA)
    regression_diffexp <- rbindlist(list(regression_diffexp,
                                         regression_cis_can_nonna, regression_trans_can_nonna, 
                                         regression_cis_can_na, regression_trans_can_na))
  }
  ## annotate ks_diffexp_type columns with up/down or whether data is available for kinase/substrate or insignificant result
  regression_diffexp$consistent_pair <- NA
  regression_diffexp$consistent_pair[((regression_diffexp$FC_kin>1 & regression_diffexp$FC_sub>1) | (regression_diffexp$FC_kin<1 & regression_diffexp$FC_sub<1)) & regression_diffexp$upstream == "kinase"] <- TRUE
  regression_diffexp$consistent_pair[((regression_diffexp$FC_kin>1 & regression_diffexp$FC_sub<1) | (regression_diffexp$FC_kin<1 & regression_diffexp$FC_sub>1)) & regression_diffexp$upstream == "phosphotase"] <- TRUE
  
  ## merge with mutational impact
  regression_diffexp <- merge(regression_diffexp, mutimpact[,c("Cancer", "Mutated_Gene", "SUBSTRATE", "SUB_MOD_RSD", "transcript",
                                                               "Mutated_Sample_Size", "Non_Mutated_Sample_Size", "Fold_Change",
                                                               "adjusted_pvalue")], 
                              by.x = c("KINASE", "SUBSTRATE", "SUB_MOD_RSD","transcript", "Cancer"), 
                              by.y = c("Mutated_Gene", "SUBSTRATE", "SUB_MOD_RSD", "transcript","Cancer"), all.x = T, sort = F)
  resultDnow <- makeOutDir()
  write.table(regression_diffexp, row.names = F, quote=F, sep = '\t', file=paste0(resultDnow, "union_diffexp_regression_genoalt.txt"))
}


