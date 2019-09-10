# Yige Wu @ WashU 2018 Jan
# summarize the prevalence of different kinds of ks pairs


# inputs ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

protein <- "kinase"
table_tn <- fread(input = paste0(ppnD, "regression/tables/annotate_tumor_normal/", protein, "_substrate_regression_cptac2p_3can_plus_normal_plus_diffexp_plus_mutimpact.txt"),
                  data.table = F)

# summary only by tumor ---------------------------------------------------
sum_table <- data.frame(table(table_tn2p[, c("class.t", "SELF","Cancer")]))
sum_table <- sum_table[sum_table$Freq > 0,]
sum_table <- sum_table[order(sum_table$class.t, sum_table$SELF, sum_table$Cancer),]
sum_table

# summary by tumor and normal ---------------------------------------------------
sum_table <- data.frame(table(table_tn2p[, c("class.tn", "SELF","Cancer")]))
sum_table <- sum_table[sum_table$Freq > 0,]
sum_table <- sum_table[order(sum_table$class.tn, sum_table$SELF, sum_table$Cancer),]
sum_table

sum_table <- data.frame(table(table_tn2p[, c("class.t", "class.tn", "SELF","Cancer")]))
sum_table <- sum_table[sum_table$Freq > 0,]
sum_table

sum_table <- data.frame(table(table_tn2p[, c("regulated.n", "regulated.t", "ks_diffexp_type", "SELF", "Cancer")]))
sum_table <- sum_table[sum_table$Freq > 0,]
sum_table <- sum_table[order(sum_table$SELF, sum_table$regulated.n, sum_table$regulated.t, sum_table$ks_diffexp_type, sum_table$Cancer),]
sum_table

#for (class.t in unique(table_tn2p$class.t)) {
# for (class.t in c("uncorrelated")) {
#  for (class.tn in unique(table_tn2p$class.tn)) {
#   for (class.tn in c("dis-regulated")) {
#     if (!(class.t %in% c("uncorrelated", "unclassified") & class.tn %in% c("uncorrelated", "unclassified"))) {
#       kspairList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = 1, prefix = prefix, 
#                                 pairlist = unique(table_tn2p$pair[table_tn2p$class.t == class.t & table_tn2p$class.tn == class.tn & table_tn2p$SELF == self]), 
#                                 id = paste0(class.t, "-", class.tn))
#     }
#   }
# }

# table_tRnD <- table_tn[!is.na(table_tn$fdr_sig.n) & !is.na(table_tn$fdr_sig.t) & table_tn$fdr_sig.n & !table_tn$fdr_sig.t,]
# table_tRnD_diffexp <- table_tRnD[!is.na(table_tRnD$diffexp_type),]
# write.table(table_tRnD, file=tn, quote=F, sep = '\t', row.names = FALSE)

