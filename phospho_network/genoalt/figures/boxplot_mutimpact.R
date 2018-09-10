# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
library(ggrepel)

# inputs ------------------------------------------------------------------
mutimpact <- fread(input = paste0(ppnD, "genoalt/tables/merge_mutation_impact/mutation_impact.txt"), data.table = F)
driver_genes <- read.delim(file = "./TCGA_data/reference_files/Consensus.Genelist.full.txt")
driver_brca <- loadGeneList(gene_type = "driver", cancer = "BRCA", is.soft.limit = "soft.limit")
driver_co <- loadGeneList(gene_type = "driver", cancer = "CO", is.soft.limit = "soft.limit")
loadGeneList(gene_type = "tsg", cancer = "BRCA", is.soft.limit = "")
oncogenes <- c(loadGeneList(gene_type = "oncogene", cancer = "BRCA", is.soft.limit = "soft.limit"), 
               loadGeneList(gene_type = "oncogene", cancer = "CO", is.soft.limit = "soft.limit"))

# set variables -----------------------------------------------------------
sig_thres <- 0.05
substrate_types_sort <- c("protein", "phosphosite")
color_cat_man <- c(colors['BRCA'], colors['OV'], colors['COAD'], "#bdbdbd"); names(color_cat_man) <- c("BRCA", "OV", "CO", "other")


# reshape inputs ----------------------------------------------------------
mutimpact$log10_pvalue <- -log10(mutimpact$p_value)
mutimpact_sig <- mutimpact[mutimpact$p_value < sig_thres,]
mutimpact_sig_pho <- mutimpact_sig[mutimpact_sig$Phosphosite != "Protein",]
mutimpact_sig_pho$site <- paste0(mutimpact_sig_pho$Substrate_Gene, ":", mutimpact_sig_pho$Phosphosite)
mutimpact_sig_pho_capped <- mutimpact_sig_pho[abs(mutimpact_sig_pho$Fold_Change) < fc2p_pho,]
mutimpact_sig_pro <- mutimpact_sig[mutimpact_sig$substrate_type == "protein",]
fc2p_pho <- quantile(x = mutimpact_sig_pho$Fold_Change, probs = 0.95)
mutimpact_sig_pho <- data.frame(mutimpact_sig_pho)
# plot impacted phosphosites only show oncogene up-regulate substrate/tsg down-regulate----------------------------------------------
tab2p <- mutimpact_sig_pho
tab2p <- tab2p[!is.na(tab2p$enzyme_driver_type) & ((tab2p$enzyme_driver_type == "oncogene" & tab2p$substrate_direction == "up") |(tab2p$enzyme_driver_type == "tsg" & tab2p$substrate_direction == "down")),]
tab2p <- tab2p[((tab2p$Mutated_Gene == "AKT1") & (tab2p$Substrate_Gene %in% driver_genes$Gene)) | (tab2p$Mutated_Gene != "AKT1"),]

lim = max(abs(max(tab2p$Fold_Change, na.rm = T)), abs(min(tab2p$Fold_Change, na.rm = T)))
cap <- min(lim, fc2p_pho)
tab2p <- tab2p[abs(tab2p$Fold_Change) > 1,]
tab2p$Fold_Change_capped <- tab2p$Fold_Change
tab2p$Fold_Change_capped[tab2p$Fold_Change > cap] <- cap
tab2p$Fold_Change_capped[tab2p$Fold_Change < -cap] <- (-cap)
tab2p$x <- paste0(tab2p$Cancer, "-", tab2p$Mutated_Gene)
tab2p$site_fc <- paste0(tab2p$site, "\n(FC=", signif(tab2p$Fold_Change, digits = 1), ")")
tab2p$text <- as.vector(tab2p$site_fc)
tab2p$text[abs(tab2p$Fold_Change) < cap] <- as.vector(tab2p$site[abs(tab2p$Fold_Change) < cap])

p <- ggplot()
p <- p + geom_point(data=tab2p, aes(x = x, y = Fold_Change_capped, fill = Cancer, color = Cancer, size = log10_pvalue), 
                    shape = 16, stroke = 0, alpha = 0.6)
p <- p + scale_fill_gradientn(name= "log2 Fold Change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
p <- p + geom_text_repel(data = tab2p, aes(x = x, y = Fold_Change_capped, label = text, color = Cancer), size = 2, force = 2)
## label the site with the max or min fold changes
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + scale_x_discrete(breaks=tab2p$x,
                          labels=tab2p$Mutated_Gene)
p <- p + scale_color_manual(values = color_cat_man)
p <- p + scale_fill_manual(values = color_cat_man)
p <- p + theme_bw()
p <- p + theme_nogrid()
p <- p + xlab('kinase')+ylab("log2FC(substrate phosphorylation")
p <- p + theme(axis.title=element_text(size=10))
p <- p + theme( axis.text.y = element_text(colour="black", size=10),
                panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 5))
p <- p + theme(axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
# p <- p + ylim(c(-cap, cap))
p
fn = paste(makeOutDir(resultD = resultD), 'enzyme_mutations_affect_substrates_right_direction_p_value', sig_thres,'.pdf',sep ="")
ggsave(filename = fn, width = 6, height = 3)





# plot impacted phosphosites only show highest/lowest fold change----------------------------------------------
tab2p <- mutimpact_sig_pho
## mark the row with the min and max fold change for each mutated gene
minmax <- vector(mode = "logical", length = nrow(tab2p))
for (mut_gene in unique(tab2p$Mutated_Gene)) {
  tab_mut <- tab2p[tab2p$Mutated_Gene == mut_gene,]
  min_fc <- min(tab_mut$Fold_Change)
  max_fc <- max(tab_mut$Fold_Change)
  minmax[tab2p$Fold_Change == min_fc] <- TRUE
  minmax[tab2p$Fold_Change == max_fc] <- TRUE
}
tab2p$minmax <- minmax
tab2p$x <- paste0(tab2p$Cancer, "-", tab2p$Mutated_Gene)

lim = max(abs(max(tab2p$Fold_Change)), abs(min(tab2p$Fold_Change)))
cap <- min(lim, fc2p_pho)
tab2p$Fold_Change_capped <- tab2p$Fold_Change
tab2p$Fold_Change_capped[tab2p$Fold_Change > cap] <- cap
tab2p$Fold_Change_capped[tab2p$Fold_Change < -cap] <- (-cap)
tab2p <- tab2p[tab2p$minmax,]
tab2p$site_fc <- paste0(tab2p$site, "\n(FC=", signif(tab2p$Fold_Change, digits = 1), ")")
tab2p$text <- as.vector(tab2p$site_fc)
tab2p$text[abs(tab2p$Fold_Change) < cap] <- as.vector(tab2p$site[abs(tab2p$Fold_Change) < cap])

p <- ggplot()
p <- p + geom_point(data=tab2p, aes(x = x, y = Fold_Change_capped, fill = Fold_Change_capped, color = Cancer, size = log10_pvalue), 
                    shape = 21, stroke = 1, alpha = 0.6)
p <- p + scale_fill_gradientn(name= "log2 Fold Change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
p <- p + geom_text_repel(data = tab2p, aes(x = x, y = Fold_Change_capped, label = text, color = Cancer), size = 2, force = 3)
## label the site with the max or min fold changes
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + scale_x_discrete(breaks=tab2p$x,
                          labels=tab2p$Mutated_Gene)
p <- p + scale_color_manual(values = color_cat_man)
p <- p + theme_bw()
p <- p + theme_nogrid()
p <- p + xlab('kinase')+ylab("log2FC(substrate phosphorylation")
p <- p + theme(axis.title=element_text(size=10))
p <- p + theme( axis.text.y = element_text(colour="black", size=10),
                panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 5))
p <- p + theme(axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
p <- p + ylim(c(-cap, cap))
p
fn = paste(makeOutDir(resultD = resultD), 'enzyme_mutations_affect_substrates_capped', signif(fc2p_pho, digits = 1), '_p_value', sig_thres,'.pdf',sep ="")
ggsave(filename = fn, width = 7, height = 3)

# plot impacted phosphosites only show substrates in driver genes----------------------------------------------
tab2p <- mutimpact_sig_pho
# tab2p <- tab2p[tab2p$Substrate_Gene %in% driver_genes$Gene,]
tab2p <- tab2p[tab2p$Substrate_Gene %in% c(driver_brca, driver_co) | (tab2p$Substrate_Gene %in% driver_genes$Gene & abs(tab2p$Fold_Change) > 5),]
# tab2p <- tab2p[tab2p$Substrate_Gene %in% c(driver_brca, driver_co),]
lim = max(abs(max(tab2p$Fold_Change)), abs(min(tab2p$Fold_Change)))
cap <- min(lim, fc2p_pho)
tab2p$Fold_Change_capped <- tab2p$Fold_Change
tab2p$Fold_Change_capped[tab2p$Fold_Change > cap] <- cap
tab2p$Fold_Change_capped[tab2p$Fold_Change < -cap] <- (-cap)
tab2p$x <- paste0(tab2p$Cancer, "-", tab2p$Mutated_Gene)
tab2p$site_fc <- paste0(tab2p$site, "\n(FC=", signif(tab2p$Fold_Change, digits = 1), ")")
tab2p$text <- as.vector(tab2p$site_fc)
tab2p$text[abs(tab2p$Fold_Change) < cap] <- as.vector(tab2p$site[abs(tab2p$Fold_Change) < cap])

p <- ggplot()
p <- p + geom_point(data=tab2p, aes(x = x, y = Fold_Change_capped, fill = Fold_Change_capped, color = Cancer, size = log10_pvalue), 
                    shape = 21, stroke = 1, alpha = 0.6)
p <- p + scale_fill_gradientn(name= "log2 Fold Change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
p <- p + geom_text_repel(data = tab2p, aes(x = x, y = Fold_Change_capped, label = text, color = Cancer), size = 2, force = 2)
## label the site with the max or min fold changes
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + scale_x_discrete(breaks=tab2p$x,
                          labels=tab2p$Mutated_Gene)
p <- p + scale_color_manual(values = color_cat_man)
p <- p + theme_bw()
p <- p + theme_nogrid()
p <- p + xlab('kinase')+ylab("log2FC(substrate phosphorylation")
p <- p + theme(axis.title=element_text(size=10))
p <- p + theme( axis.text.y = element_text(colour="black", size=10),
                panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 5))
p <- p + theme(axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
# p <- p + ylim(c(-cap, cap))
p
fn = paste(makeOutDir(resultD = resultD), 'enzyme_mutations_affect_substrates_in_driver_genes_p_value', sig_thres,'.pdf',sep ="")
ggsave(filename = fn, width = 6, height = 3)



# plot impacted proteins only show substrates in driver genes----------------------------------------------
tab2p <- mutimpact_sig_pro
# tab2p <- tab2p[tab2p$Substrate_Gene %in% driver_genes$Gene,]
tab2p <- tab2p[tab2p$Substrate_Gene %in% c(driver_brca, driver_co) | (tab2p$Substrate_Gene %in% driver_genes$Gene & abs(tab2p$Fold_Change) > 5),]
# tab2p <- tab2p[tab2p$Substrate_Gene %in% c(driver_brca, driver_co),]
lim = max(abs(max(tab2p$Fold_Change)), abs(min(tab2p$Fold_Change)))
cap <- min(lim, fc2p_pho)
tab2p$Fold_Change_capped <- tab2p$Fold_Change
tab2p$Fold_Change_capped[tab2p$Fold_Change > cap] <- cap
tab2p$Fold_Change_capped[tab2p$Fold_Change < -cap] <- (-cap)
tab2p$x <- paste0(tab2p$Cancer, "-", tab2p$Mutated_Gene)
tab2p$site_fc <- paste0(tab2p$Substrate_Gene, "\n(FC=", signif(tab2p$Fold_Change, digits = 1), ")")
tab2p$text <- as.vector(tab2p$site_fc)
tab2p$text[abs(tab2p$Fold_Change) < cap] <- as.vector(tab2p$site[abs(tab2p$Fold_Change) < cap])

p <- ggplot()
p <- p + geom_point(data=tab2p, aes(x = x, y = Fold_Change_capped, fill = Fold_Change_capped, color = Cancer, size = log10_pvalue), 
                    shape = 21, stroke = 1, alpha = 0.6)
p <- p + scale_fill_gradientn(name= "log2 Fold Change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
p <- p + geom_text_repel(data = tab2p, aes(x = x, y = Fold_Change_capped, label = text, color = Cancer), size = 2, force = 2)
## label the site with the max or min fold changes
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + scale_x_discrete(breaks=tab2p$x,
                          labels=tab2p$Mutated_Gene)
p <- p + scale_color_manual(values = color_cat_man)
p <- p + theme_bw()
p <- p + theme_nogrid()
p <- p + xlab('kinase')+ylab("log2FC(substrate phosphorylation")
p <- p + theme(axis.title=element_text(size=10))
p <- p + theme( axis.text.y = element_text(colour="black", size=10),
                panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 5))
p <- p + theme(axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
# p <- p + ylim(c(-cap, cap))
p
fn = paste(makeOutDir(resultD = resultD), 'enzyme_mutations_affect_substrates_proteins_in_driver_genes_p_value', sig_thres,'.pdf',sep ="")
ggsave(filename = fn, width = 6, height = 3)
