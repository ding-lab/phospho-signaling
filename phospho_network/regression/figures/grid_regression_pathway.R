# Yige Wu @ WashU 2018 Jan
# draw a grid showing the distribution of correlated kinase-substrate pairs are distributed across oncogenic pathways

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")


# set variables -----------------------------------------------------------
reg_nonNA <- 20
fdr_thres <- c(0.05, 0.2); names(fdr_thres) <- c("kinase", "phosphatase")
enzyme_type <- "kinase"
SELF = "trans"

# inputs ------------------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                   "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
regression <- markSigKS(regression = regression, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)

nrow(regression)


# annotate to pathways ----------------------------------------------------
table2pathway <- regression
source("./cptac2p_analysis/phospho_network/regression/tables/annotate_regression_with_pathway.R")
nrow(table2pathway)
## will duplicate a lot from original table



# get table of how many pairs are annotated to different pathways ---------
tab_pathway_per_pair <- melt(table2pathway[, c("pair", "Cancer", "GENE.path", "SUB_GENE.path", "regulated")], id.vars = c("pair", "Cancer", "regulated"))
tab_pathway_per_pair %>%
  head()
colnames(tab_pathway_per_pair) <- c("pair", "Cancer", "regulated", "variable", "pathway")
tab_pathway_per_pair$role_in_pair <- ifelse(tab_pathway_per_pair$variable == "GENE.path", enzyme_type, "substrate_phosphosite")
tab_pathway_pairs <- data.frame(table(tab_pathway_per_pair[, c("Cancer", "regulated", "role_in_pair", "pathway")]))
tab_pathway_pairs %>%
  head()

tab_regulated_pairs <- data.frame(table(tab_pathway_per_pair[, c("Cancer", "regulated")]))
tab_regulated_pairs %>%
  head()

tab_pathway_pairs.wratio <- merge(tab_pathway_pairs[tab_pathway_pairs$regulated == "TRUE", c("Cancer","role_in_pair", "pathway", "Freq")], 
                                  tab_regulated_pairs[tab_pathway_pairs$regulated == "TRUE", c("Cancer", "Freq")],
                                  by = c("Cancer"), suffixes = c(".pathway", ".regulated"), all.y = T)
tab_pathway_pairs.wratio %>%
  head()
tab_pathway_pairs.wratio$ratio_in_regulated_pair <- tab_pathway_pairs.wratio$Freq.pathway/tab_pathway_pairs.wratio$Freq.regulated
tab_pathway_pairs.wratio %>%
  head()


# compile data frame for plotting -----------------------------------------
tab2p <- tab_pathway_pairs.wratio
tab2p$y <- tab2p$Cancer
tab2p$x <- tab2p$role_in_pair
tab2p$fill <- tab2p$ratio_in_regulated_pair

## filtering
tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
tab2p <- tab2p[tab2p$pathway != "other",]
tab2p %>%
  head()
boxplot(tab2p$fill)

## ordering
tab2p$y <- order_cancer(tab2p$y)
tab2p$pathway <- order_pathway(tab2p$pathway)


# plotting ----------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = x, y = y, color = fill, size = fill), alpha = 0.6, shape = 16)
# p <- p + scale_color_manual(values = sig_colors)
p = p + scale_color_gradientn(name= "% pairs with kinase/substrae\nannotated to each pathway", na.value=NA, colours=col_paletteR(100), limit=c(0,0.15))
p <- p + facet_grid(. ~ pathway, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p <- p + xlab("Affected Protein") + ylab("Mutated Gene")
p <- p + guides(colour = guide_legend(title = "P-value"))
p <- p + theme_bw()
p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
p <- p + theme(strip.text.x = element_text(angle = 90))
p
