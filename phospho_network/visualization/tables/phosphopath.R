# Yige Wu @ WashU 2018 Jan
# generate phosphopath input

# inputs ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()
load("~/Box Sync/pan3can_shared_data/Gene_family/2015-08-01_Gene_Set.RData")
keywords <- c("cell cycle", "Hippo", "Notch", "PI3K", "TGF", "Ras", "p53", "Wnt")
dir.create("./cptac2p/cptac_shared/analysis_results/phospho_network/cytoscape/")

cancers <- c("BRCA", "OV", "CO")

# mark cancer specific or shared kinase-substrate pairs -------------------
## formatting
cons_tabs <- list()
pairs_uniq <- NULL
for (cancer in c("BRCA", "OV", "CO")) {
  # for (cancer in c("BRCA")) {
  ksea_diffexp <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  
  cons_tab <- ksea_diffexp[!is.na(ksea_diffexp$consistent) ,]
  cons_tab$pair <- paste0(cons_tab$GENE, ":", cons_tab$SUB_GENE, ":", cons_tab$SUB_MOD_RSD)
  cons_tab <- cons_tab[!duplicated(cons_tab$pair),]
  cons_tab$id <- paste0(cons_tab$pair, ":", cons_tab$enzyme_direction, ":", cons_tab$substrate_direction)
  rownames(cons_tab) <- cons_tab$id
  cons_tabs[[cancer]] <- cons_tab
  pairs_uniq <- unique(rbind(pairs_uniq, cons_tab[cons_tab$consistent, c("GENE", "SUB_GENE", "SUB_MOD_RSD", "pair", "enzyme_type", "enzyme_direction" ,"substrate_direction", "id", "KSEA_log2FC", "substrate_log2FC")]))
}

## mark consistently high/low pairs in each cancer types
table4upset <- data.frame(pairs_uniq)
table4upset <- table4upset[order(table4upset$GENE, table4upset$SUB_GENE, table4upset$SUB_MOD_RSD),]
for (cancer in c("BRCA", "OV", "CO")) {
  cons_tab <- cons_tabs[[cancer]]
  sig_pairs <- cons_tab$id[cons_tab$consistent]
  table4upset[, paste0("sig_", cancer)] <- (table4upset$id %in% sig_pairs)
}

## get the uniprot ids for each protein
genes2uniprot <- NULL
for (cancer in cancers) {
  ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  tmp <- unique(rbind(data.frame(gene_name = ksea_diffexp$GENE, uniprot_id = ksea_diffexp$KIN_ACC_ID),
                                data.frame(gene_name = ksea_diffexp$SUB_GENE, uniprot_id = ksea_diffexp$SUB_ACC_ID)))
  tmp <- tmp[!is.na(tmp$uniprot_id) & tmp$uniprot_id != "NULL",]
  genes2uniprot <- rbind(genes2uniprot, tmp)
  genes2uniprot <- unique(genes2uniprot)
}

# get kegg pathways -------------------------------------------------------
all_keggs <- names(KEGG)
keggs <- NULL
for (key in keywords) {
  keggs <- c(keggs, all_keggs[grepl(key, all_keggs, ignore.case = T)])
}
keggs <- unique(keggs)

# for (kegg_tmp in keggs) {
for (kegg_tmp in c("hsa03020\tRNA polymerase")) {
    ## get consistently high/low pairs
    cons_tab <- table4upset[table4upset$sig_BRCA | table4upset$sig_OV | table4upset$sig_CO,]
    cons_tab <- cons_tab[cons_tab$GENE %in% KEGG[[kegg_tmp]] | cons_tab$SUB_GENE %in% KEGG[[kegg_tmp]],]
    cons_tab <- merge(cons_tab, genes2uniprot, by.x = c("GENE"), by.y = c("gene_name"), all.x = T)
    colnames(cons_tab)[ncol(cons_tab)] <- "enzyme_uniprot_id"
    cons_tab <- merge(cons_tab, genes2uniprot, by.x = c("SUB_GENE"), by.y = c("gene_name"), all.x = T)
    colnames(cons_tab)[ncol(cons_tab)] <- "substrate_uniprot_id"
    
    enzyme_tab <- table4upset[table4upset$SUB_GENE %in% cons_tab$GENE,]
    enzyme_tab <- merge(enzyme_tab, genes2uniprot, by.x = c("SUB_GENE"), by.y = c("gene_name"), all.x = T)
    colnames(enzyme_tab)[ncol(enzyme_tab)] <- "substrate_uniprot_id"
    
    ## write network file
    phos_file = data.frame(uniprot = cons_tab$substrate_uniprot_id,
                           node =  paste(cons_tab$substrate_uniprot_id, cons_tab$SUB_MOD_RSD, sep = '-'),
                           site = cons_tab$SUB_MOD_RSD)
    # phos_file <- rbind(phos_file, data.frame(uniprot = enzyme_tab$substrate_uniprot_id,
    #                                          node =  paste(enzyme_tab$substrate_uniprot_id, enzyme_tab$SUB_MOD_RSD, sep = '-'),
    #                                          site = enzyme_tab$SUB_MOD_RSD))
    phos_file <- rbind(phos_file, data.frame(uniprot = cons_tab$enzyme_uniprot_id,
                                             node =  paste(cons_tab$enzyme_uniprot_id, sep = '-'),
                                             site = ""))
    phos_file <- unique(phos_file)
    tn = paste0(makeOutDir(), kegg_tmp,"_consistent_pairs.phos")
    write_delim(phos_file, tn, delim = '\t', col_names = F)
    
    attr_file <- data.frame(id =paste(cons_tab$substrate_uniprot_id, cons_tab$SUB_MOD_RSD, sep = '-'),
                            fc = cons_tab$substrate_log2FC,
                            substrate_direction = cons_tab$substrate_direction,
                            enzyme_type = "substrate")
    attr_file <- rbind(attr_file,data.frame(id =cons_tab$enzyme_uniprot_id,
                                            fc = cons_tab$KSEA_log2FC,
                                            substrate_direction = cons_tab$substrate_direction,
                                            enzyme_type = cons_tab$enzyme_type))
    attr_file <- unique(attr_file)
    tn = paste0(makeOutDir(), kegg_tmp,"_consistent_pairs_nodeAttr.txt")
    write_delim(attr_file, tn, delim = '\t', col_names = T)
}

protein <- "kinase"
tn = paste0(resultD, "classify/summarize_class/", protein,"_substrate_regression_cptac2p_3can_plus_normal_plus_diffexp_plus_mutimpact.txt")
table_tn2p <- fread(input = paste0(resultD, "regression/tables/annotate_tumor_normal/", protein, "_substrate_regression_cptac2p_3can_plus_normal_and_diffexp.txt"))

library(org.Hs.eg.db)
library(clusterProfiler)
library("pathview")
library(stringr)

# transfer IDs ------------------------------------------------------------
df <- table_tn2p[table_tn2p$regulated.n == "TRUE" & table_tn2p$regulated.t == "FALSE" & table_tn2p$SELF == "trans",]
de_genes <- unique(c(cons_tab$KINASE, cons_tab$SUBSTRATE))
de_ids <- bitr(de_genes, fromType="SYMBOL", toType=c("ENTREZID", "UNIPROT"), OrgDb="org.Hs.eg.db")
de_uniprot_uniq <- de_ids[,c("SYMBOL", "UNIPROT")]
de_uniprot_uniq <- de_uniprot_uniq[!duplicated(de_uniprot_uniq$SYMBOL),]


cancers <- c("BRCA", "OV", "CO")

# mark cancer specific or shared kinase-substrate pairs -------------------
## formatting
cons_tabs <- list()
pairs_uniq <- NULL
for (cancer in c("BRCA", "OV", "CO")) {
  # for (cancer in c("BRCA")) {
  ksea_diffexp <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  
  cons_tab <- ksea_diffexp[!is.na(ksea_diffexp$consistent) ,]
  cons_tab$pair <- paste0(cons_tab$GENE, ":", cons_tab$SUB_GENE, ":", cons_tab$SUB_MOD_RSD)
  cons_tab <- cons_tab[!duplicated(cons_tab$pair),]
  cons_tab$id <- paste0(cons_tab$pair, ":", cons_tab$enzyme_direction, ":", cons_tab$substrate_direction)
  rownames(cons_tab) <- cons_tab$id
  cons_tabs[[cancer]] <- cons_tab
  pairs_uniq <- unique(rbind(pairs_uniq, cons_tab[cons_tab$consistent, c("GENE", "SUB_GENE", "SUB_MOD_RSD", "pair", "enzyme_type", "enzyme_direction" ,"substrate_direction", "id")]))
}

## mark consistently high/low pairs in each cancer types
table4upset <- data.frame(pairs_uniq)
table4upset <- table4upset[order(table4upset$GENE, table4upset$SUB_GENE, table4upset$SUB_MOD_RSD),]
for (cancer in c("BRCA", "OV", "CO")) {
  cons_tab <- cons_tabs[[cancer]]
  sig_pairs <- cons_tab$id[cons_tab$consistent]
  table4upset[, paste0("sig_", cancer)] <- (table4upset$id %in% sig_pairs)
}

## mark unique or shared across cancer types
for (cancer in cancers) {
  rest_cancers <- cancers[cancers != cancer]
  table4upset[, paste0("uniq_", cancer)] <- (table4upset[, paste0("sig_", cancer)] & rowSums(table4upset[, paste0("sig_", rest_cancers)]) == 0)
}
table4upset[, "shared"] <- (rowSums(table4upset[, paste0("sig_", cancers)]) == length(cancers))
table4upset[, paste0("sig_", cancers)] <- table4upset[, paste0("sig_", cancers)]*1
table4upset <- data.frame(table4upset)

# write network file ------------------------------------------------------
df <- merge(df, de_uniprot_uniq, by.x = c('KINASE'), by.y = c('SYMBOL'), all.x = T)
colnames(df)[ncol(df)] <- "KINASE_uniprot"
df <- merge(df, de_uniprot_uniq, by.x = c('SUBSTRATE'), by.y = c('SYMBOL'), all.x = T)
colnames(df)[ncol(df)] <- "SUBSTRATE_uniprot"

phos_file = data.frame(uniprot = cons_tab$SUBSTRATE_uniprot,
                       node =  paste(cons_tab$SUBSTRATE_uniprot, cons_tab$SUB_MOD_RSD, sep = '-'),
                       site = cons_tab$SUB_MOD_RS)
phos_file <- rbind(phos_file, data.frame(uniprot = cons_tab$KINASE_uniprot,
                                         node =  paste(cons_tab$KINASE_uniprot, "protein", sep = '-'),
                                         site = "protein"))
phos_file <- unique(phos_file)
tn = paste0(resultDnow, protein,"_substrate_regression_cptac2p_3can_nRtD.phos")
write_delim(phos_file, tn, delim = '\t', col_names = F)


# write node attribute file -----------------------------------------------

# write edge attribute file -----------------------------------------------
