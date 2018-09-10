# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
sig_thres <- 0.1
substrate_types_sort <- c("protein", "phosphosite")
substrate_types_sort2 <- c("protein", "phosphosite-only", "phosphosite-protein")

# inputs ------------------------------------------------------------------
mutimpact <- fread(input = paste0(ppnD, "genoalt/tables/merge_mutation_impact/mutation_impact.txt"), data.table = F)
mutimpact_sig <- mutimpact[mutimpact$p_value < sig_thres,]
mutimpact_sig.m <- data.frame(table(mutimpact_sig[, c("Mutated_Gene", "Cancer", "SELF", "substrate_type")]))
mutimpact_sig.m$substrate_type <- factor(mutimpact_sig.m$substrate_type, levels = substrate_types_sort)

# show protein/phosphosites affected------------------------
p <- ggplot()
p <- p + geom_tile(data = mutimpact_sig.m, mapping = aes(x = Cancer, y = Mutated_Gene, fill = Freq))
p = p + scale_fill_gradientn(name= "", na.value=NA, colours=col_paletteR(max(mutimpact_sig.m$Freq)), lim = c(1, max(mutimpact_sig.m$Freq)))
p = p + geom_text(data = mutimpact_sig.m[mutimpact_sig.m$Freq > 0,], mapping = aes(x=Cancer, y=Mutated_Gene, label = Freq), size = 5)
p <- p + facet_grid(.~SELF + substrate_type, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p <- p + theme_nogrid()
p <- p + theme(panel.spacing.x = unit(0, "lines"))
p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
p
fn = paste(makeOutDir(resultD = resultD), 'driver_druggable_enzyme_mutations_affect_substrate_numbers.pdf',sep ="")
ggsave(filename = fn, width = 5, height = 5)


# show phosphosites by whether protein is affected ------------------------
mutimpact_sig.m2 <- merge(mutimpact_sig.m, mutimpact_sig.m[mutimpact_sig.m$substrate_type == "protein",],
                          by = c("Mutated_Gene", "Cancer", "SELF"), all.x = T, suffixes = c("", ".pro"))
mutimpact_sig.m2 <- data.frame(mutimpact_sig.m2)

tmp <- as.vector(mutimpact_sig.m2$substrate_type)
tmp[tmp == "phosphosite" & mutimpact_sig.m2$Freq.pro > 0] <- "phosphosite-protein"
tmp[tmp == "phosphosite" & mutimpact_sig.m2$Freq.pro == 0] <- "phosphosite-only"
tmp <- factor(tmp, levels = substrate_types_sort2)
mutimpact_sig.m2$substrate_type <- tmp


p <- ggplot()
p <- p + geom_tile(data = mutimpact_sig.m2, mapping = aes(x = Cancer, y = Mutated_Gene, fill = Freq))
p = p + scale_fill_gradientn(name= "", na.value=NA, colours=col_paletteR(max(mutimpact_sig.m$Freq)), lim = c(1, max(mutimpact_sig.m$Freq)))
p = p + geom_text(data = mutimpact_sig.m2[mutimpact_sig.m2$Freq > 0,], mapping = aes(x=Cancer, y=Mutated_Gene, label = Freq), size = 4)
p <- p + facet_grid(.~SELF + substrate_type, drop=T, space = "free", scales = "fixed")#, space = "free", scales = "free")
p <- p + theme_nogrid()
p <- p + theme(panel.spacing.x = unit(0, "lines"))
p <- p + theme(strip.text.x = element_text(size = 5))
p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
p <- p + theme(axis.text.x = element_text(size = 8, face = "bold", angle = 45))
p
fn = paste(makeOutDir(resultD = resultD), 'driver_druggable_enzyme_mutations_affect_substrate_numbers2.pdf',sep ="")
ggsave(filename = fn, width = 5, height = 3)

