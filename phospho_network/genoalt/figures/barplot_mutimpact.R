# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
sig_thres <- 0.05
substrate_types_sort <- c("protein", "phosphosite")
substrate_types_sort2 <- c("protein", "phosphosite-only", "phosphosite-protein")
color_cat_man <- c(colors['BRCA'], colors['OV'], colors['COAD'], "#bdbdbd"); names(color_cat_man) <- c("BRCA", "OV", "CO", "other")
cancers_sort <- c("BRCA", "OV", "CO")

# inputs ------------------------------------------------------------------
mutimpact <- fread(input = paste0(ppnD, "genoalt/tables/merge_mutation_impact/mutation_impact.txt"), data.table = F)
# mutimpact <- fread(input = paste0("./PhosphoDrug/PhosphoDrug_shared_data/analysis_results/mutationimpact/07232018/correctedPvalueDriver"), data.table = F)
mutimpact$sig <- ifelse(mutimpact$p_value < sig_thres, paste0("P_value<", sig_thres), paste0("P_value>=", sig_thres))
mutimpact_sig <- mutimpact[mutimpact$p_value < sig_thres,]
mutimpact_sig_pro <- mutimpact_sig[mutimpact_sig$substrate_type == "protein",]
mutimpact_sig_pro_oncogene <- mutimpact_sig_pro[!is.na(mutimpact_sig_pro$enzyme_driver_type) & mutimpact_sig_pro$enzyme_driver_type == "oncogene" & mutimpact_sig_pro$substrate_direction == "up",]
mutimpact_sig_pro_tsg <- mutimpact_sig_pro[!is.na(mutimpact_sig_pro$enzyme_driver_type) & mutimpact_sig_pro$enzyme_driver_type == "tsg" & mutimpact_sig_pro$substrate_direction == "down",]
mutimpact_sig_pho <- mutimpact_sig[mutimpact_sig$substrate_type == "phosphosite",]
mutimpact_sig_pho_oncogene <- mutimpact_sig_pho[!is.na(mutimpact_sig_pho$enzyme_driver_type) & mutimpact_sig_pho$enzyme_driver_type == "oncogene" & mutimpact_sig_pho$substrate_direction == "up",]
mutimpact_sig_pho_tsg <- mutimpact_sig_pho[!is.na(mutimpact_sig_pho$enzyme_driver_type) & mutimpact_sig_pho$enzyme_driver_type == "tsg" & mutimpact_sig_pho$substrate_direction == "down",]


# show enzymes mutations affecting phosphosites  ------------------------
# mutimpact_sig_pho_tab <- data.frame(table(mutimpact_sig_pho[, c("Mutated_Gene", "Cancer", "substrate_direction")]))
mutimpact_sig_pho_tab <- data.frame(table(mutimpact[mutimpact$Phosphosite != "Protein", c("Mutated_Gene", "Cancer", "sig")]))
tab2p <- mutimpact_sig_pho_tab
tab2p <- tab2p[tab2p$Mutated_Gene %in% tab2p$Mutated_Gene[tab2p$sig == paste0("P_value<", sig_thres) & tab2p$Freq > 0],]
tab2p <- tab2p[order(tab2p$sig, tab2p$Freq, decreasing = T),]
tab2p <- tab2p[order(tab2p$Cancer, decreasing = F),]
color_cat <- as.vector(tab2p$Cancer)
color_cat[which(tab2p$sig == paste0("P_value>=", sig_thres))] <- "other"
tab2p$color_cat <- color_cat
tab2p$color_cat <- factor(tab2p$color_cat, levels = c("other", cancers_sort))
tab2p$Mutated_Gene <- factor(tab2p$Mutated_Gene, levels = unique(tab2p$Mutated_Gene[tab2p$sig == paste0("P_value<", sig_thres)]))
p <- ggplot()
# p <- p + geom_bar(data=tab2p, aes(x = Mutated_Gene, y = Freq, fill = substrate_direction, color = Cancer), stat="identity", position='stack')
p <- p + geom_bar(data=tab2p, aes(x = Mutated_Gene, y = Freq, fill = color_cat), stat="identity", position='stack', color = "#000000")
p <- p + geom_text(data=tab2p[tab2p$Freq > 0 & tab2p$sig == paste0("P_value<", sig_thres),], aes(x = Mutated_Gene, y = Freq, color = color_cat, label = Freq), nudge_y = 0.2)
# p <- p + scale_fill_manual(values = c("down" = "#1F78B4", "up" = "#E31A1C"))
p <- p + scale_fill_manual(values = color_cat_man)
p <- p + scale_color_manual(values = color_cat_man)
p <- p + scale_y_log10()
p <- p + theme_nogrid()
p <- p + theme(axis.title.y = element_blank())
p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
p <- p + theme(axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
p
fn = paste(makeOutDir(resultD = resultD), 'enzyme_mutations_affect_phosphosites_self_included_p_value', sig_thres,'.pdf',sep ="")
ggsave(filename = fn, width = 5, height = 3)


# show enzymes mutations affecting phosphosites  ------------------------
# mutimpact_sig_pro_tab <- data.frame(table(mutimpact_sig_pro[, c("Mutated_Gene", "Cancer", "substrate_direction")]))
mutimpact_sig_pro_tab <- data.frame(table(mutimpact[mutimpact$substrate_type == "protein", c("Mutated_Gene", "Cancer", "sig")]))
tab2p <- mutimpact_sig_pro_tab
tab2p <- tab2p[tab2p$Mutated_Gene %in% tab2p$Mutated_Gene[tab2p$sig == paste0("P_value<", sig_thres) & tab2p$Freq > 0],]
tab2p <- tab2p[order(tab2p$sig, tab2p$Freq, decreasing = T),]
tab2p <- tab2p[order(tab2p$Cancer, decreasing = F),]
color_cat <- as.vector(tab2p$Cancer)
color_cat[which(tab2p$sig == paste0("P_value>=", sig_thres))] <- "other"
tab2p$color_cat <- color_cat
tab2p$color_cat <- factor(tab2p$color_cat, levels = c("other", cancers_sort))
tab2p$Mutated_Gene <- factor(tab2p$Mutated_Gene, levels = unique(tab2p$Mutated_Gene[tab2p$sig == paste0("P_value<", sig_thres)]))
p <- ggplot()
# p <- p + geom_bar(data=tab2p, aes(x = Mutated_Gene, y = Freq, fill = substrate_direction, color = Cancer), stat="identity", position='stack')
p <- p + geom_bar(data=tab2p, aes(x = Mutated_Gene, y = Freq, fill = color_cat), stat="identity", position='stack', color = "#000000")
p <- p + geom_text(data=tab2p[tab2p$Freq > 0 & tab2p$sig == paste0("P_value<", sig_thres),], aes(x = Mutated_Gene, y = Freq, color = color_cat, label = Freq), nudge_y = 0.2)
# p <- p + scale_fill_manual(values = c("down" = "#1F78B4", "up" = "#E31A1C"))
p <- p + scale_fill_manual(values = color_cat_man)
p <- p + scale_color_manual(values = color_cat_man)
p <- p + scale_y_log10()
p <- p + theme_nogrid()
p <- p + theme(axis.title.y = element_blank())
p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
p <- p + theme(axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
p
fn = paste(makeOutDir(resultD = resultD), 'enzyme_mutations_affect_proteins_self_included_p_value', sig_thres,'.pdf',sep ="")
ggsave(filename = fn, width = 5, height = 3)

