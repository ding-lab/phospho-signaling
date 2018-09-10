# Yige Wu @ WashU 2018 July
## output table for cytoscape visualization

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
sig_thres <- 0.1


# inputs ------------------------------------------------------------------
mutimpact <- fread(input = paste0(ppnD, "genoalt/tables/merge_mutation_impact/mutation_impact.txt"), data.table = F)
mutimpact_sig <- mutimpact[mutimpact$p_value < sig_thres,]

# load cancer hallmark pathways -------------------------------------------
load("~/Box Sync/pan3can_shared_data/Gene_family/2015-08-01_Gene_Set.RData")
keywords <- c("cell cycle", "PI3K","Akt","MAPK","adherens","Mismatch","EMT", "apoptosis","immunological","stromal","transmembrane","receptors","integrin",
              "TGFβ","LKB1","AMPK","TSC","mTOR","Ras","Notch","Wnt", "catenin",  "p53","RTK","erbb", "Spliceosome")
all_keggs <- names(KEGG)
keggs <- NULL
for (key in keywords) {
  keggs <- c(keggs, all_keggs[grepl(key, all_keggs, ignore.case = T)])
}
keggs <- unique(keggs)

## annotate cancer hallmark pathways

## annotate kinases to pathways
sub_tab <- data.frame(Substrate_Gene = unique(mutimpact_sig$Substrate_Gene))
for (kegg in keggs) {
  sub_tab[, kegg] <- (sub_tab$Substrate_Gene %in% KEGG[[kegg]])
}
sub_tab$hallmark_count <- rowSums(sub_tab[, keggs])
hallmark_assign <- vector(mode = "character", length = nrow(sub_tab))
for (i in 1:nrow(sub_tab)) {
  if (sub_tab$hallmark_count[i] == 0) {
    hallmark_assign[i] <- "other"
  }
  if (sub_tab$hallmark_count[i] > 0) {
    hallmark_assign[i] <- keggs[min(which(sub_tab[i, keggs] == T))]
  }
}
sub_tab$hallmark_assign <- hallmark_assign

mutimpact_sig <- merge(mutimpact_sig, sub_tab[, c("Substrate_Gene", "hallmark_assign")], by = c('Substrate_Gene'), all.x = T)

## only keep substrates at hallmark pathways
mutimpact_sig_hallmark <- mutimpact_sig[mutimpact_sig$hallmark_assign != "other",]
mutimpact_sig_hallmark$hallmark_assign <- NULL
mutimpact_sig_hallmark$log10_pvalue <- -log10(mutimpact_sig_hallmark$p_value)
## write network file
write.table(x = mutimpact_sig_hallmark, file = paste0(makeOutDir(resultD = resultD), "mutimpact_sig_hallmark.txt"), 
            quote = F, row.names = F, col.names = T, sep = "\t")
for (cancer in unique(mutimpact_sig_hallmark$Cancer)) {
  mutimpact_sig_hallmark_can <- mutimpact_sig_hallmark[mutimpact_sig_hallmark$Cancer == cancer,]
  write.table(x = mutimpact_sig_hallmark_can, file = paste0(makeOutDir(resultD = resultD), cancer,  "_mutimpact_sig_hallmark.txt"), 
              quote = F, row.names = F, col.names = T, sep = "\t")
}
  



# write node file ---------------------------------------------------------
enzyme_tab <- data.frame(Substrate_Gene = unique(mutimpact_sig$Mutated_Gene))
for (kegg in keggs) {
  enzyme_tab[, kegg] <- (enzyme_tab$Substrate_Gene %in% KEGG[[kegg]])
}
enzyme_tab$hallmark_count <- rowSums(enzyme_tab[, keggs])
hallmark_assign <- vector(mode = "character", length = nrow(enzyme_tab))
for (i in 1:nrow(enzyme_tab)) {
  if (enzyme_tab$hallmark_count[i] == 0) {
    hallmark_assign[i] <- "other"
  }
  if (enzyme_tab$hallmark_count[i] > 0) {
    hallmark_assign[i] <- keggs[min(which(enzyme_tab[i, keggs] == T))]
  }
}
enzyme_tab$hallmark_assign <- hallmark_assign


tmp <- sub_tab[,c("Substrate_Gene", "hallmark_assign")]
colnames(tmp) <- c("shared_name", "hallmark_assign")
node_tab <- tmp
tmp <- enzyme_tab[,c("Substrate_Gene", "hallmark_assign")]
colnames(tmp) <- c("shared_name", "hallmark_assign")
node_tab <- rbind(node_tab, tmp)
node_tab$is_enzyme <- (node_tab$shared_name %in% mutimpact_sig_hallmark$Mutated_Gene)
node_tab$pathway <- str_split_fixed(string = node_tab$hallmark_assign, pattern = "\t", 2)[,2]
node_tab$pathway[node_tab$pathway == ""] <- "other"
node_tab$hallmark_assign <- NULL
write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), "mutimpact_sig_hallmark_node_tab.txt"), 
            quote = F, row.names = F, col.names = T, sep = ",")


table(mutimpact_sig_hallmark[mutimpact_sig_hallmark$SELF == "trans" & mutimpact_sig_hallmark$Mutated_Gene == "AKT1", c("substrate_type", "hallmark_assign")])
mutimpact_sig_hallmark[mutimpact_sig_hallmark$Substrate_Gene %in% mutimpact_sig_hallmark$Mutated_Gene & mutimpact_sig_hallmark$SELF == "trans",]


mutimpact_sig[mutimpact_sig$Substrate_Gene %in% mutimpact_sig$Mutated_Gene & mutimpact_sig$SELF == "trans",]
