# Yige Wu @ WashU 2018 Jan
# generate phosphopath input

# inputs ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()

protein <- "kinase"
tn = paste0(resultD, "classify/summarize_class/", protein,"_substrate_regression_cptac2p_3can_plus_normal_plus_diffexp_plus_mutimpact.txt")
table_tn2p <- fread(input = paste0(resultD, "regression/tables/annotate_tumor_normal/", protein, "_substrate_regression_cptac2p_3can_plus_normal_and_diffexp.txt"))

library(org.Hs.eg.db)
library(clusterProfiler)
library("pathview")
library(stringr)

# transfer IDs ------------------------------------------------------------
df <- table_tn2p[table_tn2p$regulated.n == "TRUE" & table_tn2p$regulated.t == "FALSE" & table_tn2p$SELF == "trans",]
de_genes <- unique(c(df$KINASE, df$SUBSTRATE))
de_ids <- bitr(de_genes, fromType="SYMBOL", toType=c("ENTREZID", "UNIPROT"), OrgDb="org.Hs.eg.db")
de_uniprot_uniq <- de_ids[,c("SYMBOL", "UNIPROT")]
de_uniprot_uniq <- de_uniprot_uniq[!duplicated(de_uniprot_uniq$SYMBOL),]


# write network file ------------------------------------------------------
df <- merge(df, de_uniprot_uniq, by.x = c('KINASE'), by.y = c('SYMBOL'), all.x = T)
colnames(df)[ncol(df)] <- "KINASE_uniprot"
df <- merge(df, de_uniprot_uniq, by.x = c('SUBSTRATE'), by.y = c('SYMBOL'), all.x = T)
colnames(df)[ncol(df)] <- "SUBSTRATE_uniprot"

phos_file = data.frame(uniprot = df$SUBSTRATE_uniprot,
                       node =  paste(df$SUBSTRATE_uniprot, df$SUB_MOD_RSD, sep = '-'),
                       site = df$SUB_MOD_RS)
phos_file <- rbind(phos_file, data.frame(uniprot = df$KINASE_uniprot,
                                         node =  paste(df$KINASE_uniprot, "protein", sep = '-'),
                                         site = "protein"))
phos_file <- unique(phos_file)
tn = paste0(resultDnow, protein,"_substrate_regression_cptac2p_3can_nRtD.phos")
write_delim(phos_file, tn, delim = '\t', col_names = F)


# write node attribute file -----------------------------------------------
attr_file <- data.frame(id =paste(df$SUBSTRATE_uniprot, df$SUB_MOD_RSD, sep = '-'),
                        fc = df$FC_sub,
                        label = df$SUB_MOD_RSD)
attr_file <- rbind(attr_file,data.frame(id =df$KINASE_uniprot,
                                        fc = df$FC_kin,
                                        label = df$KINASE),
                   data.frame(id =paste(df$KINASE_uniprot, 'protein', sep = '-'),
                                                                      fc = NA,
                                                                      label = NA))
attr_file <- unique(attr_file)
attr_file$log2fc <- log2(attr_file$fc)
tn = paste0(resultDnow, protein,"_substrate_regression_cptac2p_3can_nRtD_nodeAttr.txt")
write_delim(attr_file, tn, delim = '\t', col_names = T)

attr_file <- attr_file[!is.na(attr_file$fc),]
attr_file <- unique(attr_file)
attr_file$log2fc <- log2(attr_file$fc)
tn = paste0(resultDnow, protein,"_substrate_regression_cptac2p_3can_nRtD_nodeAttr.txt")
write_delim(attr_file, tn, delim = '\t', col_names = T)

# write edge attribute file -----------------------------------------------
df$ks <- paste0(df$KINASE, df$SUBSTRATE)
tn = paste0(resultDnow, protein,"_substrate_regression_cptac2p_3can_nRtD_edgeAttr.txt")
write_delim(df[which(df$ks=="MAPK1MAP2K1")[1],], tn, delim = '\t', col_names = T)
