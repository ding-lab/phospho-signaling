# Yige Wu @ WashU 2018 Jul
# identify enzyme_substrate pairs where different cancers use different phosphosites on the same substrate by the same kinase/phosphatase

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# variables ---------------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
reg_sig <- 0.1
color_cat_man <- c(colors['BRCA'], colors['OV'], colors['COAD'], "#bdbdbd"); names(color_cat_man) <- c("BRCA", "OV", "CO", "other")

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

# inputs -------------------------------------------------------------------
enzyme_type <- "phosphatase"
sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/phosphatase_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
sup_cans_tab_en$pair <- paste0(sup_cans_tab_en$KINASE, ":", sup_cans_tab_en$SUBSTRATE, ":", sup_cans_tab_en$SUB_MOD_RSD)
sup_cans_tab_en$pair_pro <- paste0(sup_cans_tab_en$KINASE, ":", sup_cans_tab_en$SUBSTRATE)
colnames(sup_cans_tab_en)[1:2] <- c("GENE", "SUB_GENE")
sup_cans_tab_en <- markSigSiteCan(sup_cans_tab_en, sig_thres = reg_sig, enzyme_type = enzyme_type)

sup_cans_tab <- sup_cans_tab_en
## annotate kinase substrate regulation
sup_cans_tab$regulated <- (sup_cans_tab$coef_sig & sup_cans_tab$fdr_sig)

## annotate kinases to pathways
sub_tab <- data.frame(sub = unique(sup_cans_tab[,c("SUB_GENE")]))
for (kegg in keggs) {
  sub_tab[, kegg] <- (sub_tab$sub %in% KEGG[[kegg]])
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
hallmark_assign[hallmark_assign != "other"] <- str_split_fixed(string = hallmark_assign[hallmark_assign != "other"], pattern = "\t", 2)[,2]
sub_tab$hallmark_assign_sub <- hallmark_assign
sup_cans_tab <- merge(sup_cans_tab, sub_tab[,c('sub',"hallmark_assign_sub")], by.x = c('SUB_GENE'), by.y = c('sub'), all.x = T)

## annotate kinases to pathways
enzyme_tab <- data.frame(enzyme = unique(sup_cans_tab[,c("GENE")]))
for (kegg in keggs) {
  enzyme_tab[, kegg] <- (enzyme_tab$enzyme %in% KEGG[[kegg]])
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
hallmark_assign[hallmark_assign != "other"] <- str_split_fixed(string = hallmark_assign[hallmark_assign != "other"], pattern = "\t", 2)[,2]
enzyme_tab$hallmark_assign_enz <- hallmark_assign
# sup_cans_tab <- merge(sup_cans_tab, enzyme_tab[,c('enzyme',"hallmark_assign_enz")], by.x = c('GENE'), by.y = c('enzyme'), all.x = T)

# summarize the kinase-substrate pairs ------------------------------------
sup_cans_tab <- data.frame(sup_cans_tab)
sum1 <- data.frame(table(unique(sup_cans_tab[sup_cans_tab$regulated,c("pair", "pair_pro", "Cancer", "regulated")])[,c("pair_pro", "Cancer")]))
sum1 <- sum1[sum1$Freq > 0,]
sum2 <- data.frame(table(sum1$pair_pro))
sum_pair <- data.frame(table(unique(sup_cans_tab[sup_cans_tab$regulated,c("pair",  "Cancer" )])[,c("pair")]))
sum2 <- sum2[sum2$Freq > 1,]
sup_cans_tab2 <- sup_cans_tab[sup_cans_tab$pair_pro %in% sum2$Var1,]


## get a list of enzyme substrate pairs with regression output in all 3 cancers
sup_cans_tab_nona <- sup_cans_tab[!is.na(sup_cans_tab$fdr_sig),]
pairs_wvalue_tab <- data.frame(table(unique(sup_cans_tab_nona[,c("pair", "Cancer")])[,c("pair")]))

# find enzyme substrate pairs with different sites regulated in each cancer ----------------------------
sup_cans_tab2_trans <- sup_cans_tab2[sup_cans_tab2$SELF == "trans",]
sup_cans_tab2_trans_reg <- sup_cans_tab2_trans[sup_cans_tab2_trans$regulated & !is.na(sup_cans_tab2_trans$regulated),]
sup_cans_tab2_trans_reg <- sup_cans_tab2_trans_reg[sup_cans_tab2_trans_reg$pair %in% pairs_wvalue_tab$Var1[pairs_wvalue_tab$Freq == 3],]
sup_cans_tab2_trans_reg <- sup_cans_tab2_trans_reg[order(sup_cans_tab2_trans_reg$pair),]
## remove the shared phosphosites
sup_cans_tab2_trans_reg <- sup_cans_tab2_trans_reg[!(sup_cans_tab2_trans_reg$pair %in% sum_pair$Var1[sum_pair$Freq == 3]),]
sum21 <- data.frame(table(unique(sup_cans_tab2_trans_reg[sup_cans_tab2_trans_reg$regulated,c("pair", "pair_pro", "Cancer", "regulated")])[,c("pair_pro", "Cancer")]))
sum21 <- sum21[sum21$Freq > 0,]
sum22 <- data.frame(table(sum21$pair_pro))
sum22 <- sum22[sum22$Freq > 1,]
sup_cans_tab2_trans_reg2 <- sup_cans_tab2_trans_reg[sup_cans_tab2_trans_reg$pair_pro %in% sum22$Var1,]
sup_cans_tab2_trans_reg_hallmark <- sup_cans_tab2_trans_reg2[sup_cans_tab2_trans_reg2$hallmark_assign_sub != "other" & sup_cans_tab2_trans_reg2$hallmark_assign_enz != "other",]
sup_cans_tab2_trans_reg_hallmark$hallmark_assign_sub <- factor(sup_cans_tab2_trans_reg_hallmark$hallmark_assign_sub, levels = c(keggs, "other"))
sup_cans_tab2_trans_reg_hallmark$x_print <- paste0(sup_cans_tab2_trans_reg_hallmark$Cancer, "-", sup_cans_tab2_trans_reg_hallmark$SUB_MOD_RSD)

# df <- unique(sup_cans_tab2_trans_reg_hallmark[, c("pair_pro", "Cancer", "SUB_MOD_RSD", "pair", "GENE")])
df <- unique(sup_cans_tab2_trans_reg2[, c("pair_pro", "Cancer", "SUB_MOD_RSD", "pair", "GENE")])
df <- merge(df, sum21, all.x = T)
df$bar_len <- 1/df$Freq
df$id <- paste0(df$pair_pro, ":", df$Cancer)
df <- df[order(df$id),]
y_print <- vector(mode = "numeric", length = nrow(df))
for (i in 1:nrow(df)) {
  if (duplicated(df$id)[i]) {
    y_print[i] <- y_print[i-1] + df$bar_len[i]
  } else {
    y_print[i] <- df$bar_len[i]/2
  }
}
df$y_print <- y_print
df$Cancer <- factor(df$Cancer, levels = cancers_sort)
p <- ggplot()
p <- p + geom_bar(data = df, mapping = aes(x = pair_pro, y = bar_len, fill = Cancer), stat = "identity", position = "stack", color = "black")
p <- p + geom_text(data = df, mapping = aes(x = pair_pro, y = y_print, label = SUB_MOD_RSD), size = 2.5)
p <- p + scale_fill_manual(values = color_cat_man)
p <- p + facet_grid(GENE~Cancer, drop=T,shrink = T, space = "free",scales = "free")#, space = "free", scales = "free")
p <- p + coord_flip()
p <- p + theme_nogrid()
p <- p + theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
               axis.title.y = element_blank(), axis.title.x = element_blank(),
               # strip.text.x = element_text(size = 5, angle = 90),
               panel.spacing.y = unit(0, "lines"),
               panel.spacing.x = unit(0, "lines"))
p
fn = paste0(makeOutDir(resultD = resultD),'different_site_same_es_pairs_by_different_cancers.pdf')
ggsave(file=fn, height=6.5, width=6)

sup_cans_tab2_trans4cytoscape <- sup_cans_tab2_trans[sup_cans_tab2_trans$pair %in% df$pair,]
sup_cans_tab2_trans4cytoscape <- sup_cans_tab2_trans4cytoscape[sup_cans_tab2_trans4cytoscape$regulated & !is.na(sup_cans_tab2_trans4cytoscape$regulated),]
sup_cans_tab2_trans4cytoscape$scale_FDR <- -log10(sup_cans_tab2_trans4cytoscape$FDR_pho_kin)
sup_cans_tab2_trans4cytoscape$site <- paste0(sup_cans_tab2_trans4cytoscape$SUB_GENE, ":", sup_cans_tab2_trans4cytoscape$SUB_MOD_RSD)

node_tab <- data.frame(shared_name = sup_cans_tab2_trans4cytoscape$site,
                       target_type = "phosphosite", name_print = sup_cans_tab2_trans4cytoscape$SUB_MOD_RSD)
tmp <- data.frame(shared_name = sup_cans_tab2_trans4cytoscape$GENE,
                  target_type = enzyme_type, name_print = sup_cans_tab2_trans4cytoscape$GENE)
node_tab <- rbind(node_tab, tmp)
tmp <- data.frame(shared_name = unique(sup_cans_tab2_trans4cytoscape$SUB_GENE[!(sup_cans_tab2_trans4cytoscape$SUB_GENE %in% sup_cans_tab2_trans4cytoscape$GENE)]),
                  target_type = "substrate")
tmp$name_print <- tmp$shared_name
node_tab <- rbind(node_tab, tmp)
node_tab <- unique(node_tab)
tmp <- sup_cans_tab2_trans4cytoscape[!duplicated(sup_cans_tab2_trans4cytoscape$site),]
tmp$GENE <- tmp$site
tmp$site <- tmp$SUB_GENE
tmp$Cancer <- "placeholder"
tmp$coef_pho_kin <- 0
sup_cans_tab2_trans4cytoscape <- rbind(sup_cans_tab2_trans4cytoscape, tmp)

write.table(x = sup_cans_tab2_trans4cytoscape, file = paste0(makeOutDir(resultD = resultD), "network_tab4cytoscape.txt"), quote = F, row.names = F, sep = "\t")
write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), "node_tab4cytoscape.txt"), quote = F, row.names = F, sep = "\t")

