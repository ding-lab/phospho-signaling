# Yige Wu @ WashU 2018 Jan
# draw volcano plots for regression result (3 cancer types together)

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggrepel)

# set variables -----------------------------------------------------------
reg_nonNA <- 20
fdr_thres <- c(0.05, 0.2); names(fdr_thres) <- c("kinase", "phosphatase")

# inputs ------------------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                   "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)


# kinases-trans-all-cancer -----------------------------------------------------------------
## set variables
### for filtering
enzyme_type <- "kinase"
SELF = "trans"
### for plotting
text_cutoff <- -log10(fdr_pk)
y_cap <- 15
x_cap <- 2

tab2p <- regression
tab2p <- markSigKS(regression = tab2p, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
tab2p <- tab2p[tab2p$FDR_pho_kin > 0,]
tab2p <- tab2p[tab2p$SELF == SELF & tab2p$enzyme_type == enzyme_type,]
## take out extra not significant ones
tab2p_keep <- tab2p[tab2p$regulated,]
tab2p_sample <- tab2p[!(tab2p$regulated),]
tab2p_sample <- tab2p_sample[sample(x = 1:nrow(tab2p_sample), size = 20000, replace = F),]
tab2p <- rbind(tab2p_keep, tab2p_sample)
tab2p$coef_capped <- remove_outliers_cap(x = tab2p$coef_pho_kin, cap = x_cap)
tab2p$log_FDR <- -log10(tab2p$FDR_pho_kin)
tab2p$log_FDR_capped <- remove_outliers_cap(x = tab2p$log_FDR, cap = y_cap)
tab2p$y <- tab2p$log_FDR_capped
tab2p$x <- tab2p$coef_capped

## add a column for color
tab2p$point_color <- "other"
tab2p$point_color[tab2p$regulated & !is.na(tab2p$regulated)] <- tab2p$Cancer[tab2p$regulated& !is.na(tab2p$regulated)]

## order
tab2p <- tab2p[order(tab2p$log_FDR, decreasing = T),]
tab2p$Cancer <- order_cancer(tab2p$Cancer)
tab2p %>%
  head()

p = ggplot(data = tab2p)
p = p + geom_point(mapping = aes(x=x, y= y, color = point_color), stroke = 0, alpha = 0.5, shape = 16, size = 1)
p = p + scale_color_manual(values = color_cancers2)
p = p + facet_wrap(facets = c("Cancer"), scales = "fixed", nrow = 2)
p = p + geom_hline(yintercept = text_cutoff, color = 'black', linetype = 2)
p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
p = p + theme_minimal() + theme_nogrid()
# p = p + geom_text_repel(data = tab2p[1:top_num,], mapping = aes(x=x, y=y, label= as.character(pair),), size= 2, alpha = 0.7, force = 1)
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"), panel.background = element_rect(color = NA),
              strip.text.y = element_text(size = 12), strip.background.x = element_rect(fill = "white", color = "white"))
p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
p = p + theme(legend.position = "none")
p
fn = paste(makeOutDir(resultD = resultD), "volcano_trans_", "cancers" , "_", enzyme_type,'.pdf',sep ="")
ggsave(file=fn, height=4, width=5.5, useDingbats=FALSE)

stop()

# phosphatase all cancer -----------------------------------------------------------------

for (text_cutoff in c(-log10(fdr_pp))) {
  ## get KSEA scores integrated with differential phosphorylation results
  tab2p <- regression[regression$SELF == "trans",]
  tab2p <- tab2p[sample(x = 1:nrow(tab2p), size = nrow(tab2p), replace = F),]
  tab2p$coef_capped <- remove_outliers(x = as.vector(tab2p$coef_pho_kin), out_thres = 2, na.rm = T)
  p = ggplot()
  p = p + geom_point(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y=-log10(FDR_pho_kin),
                         color = Cancer), stroke = 0, alpha = 0.8, shape = 16)
  p = p + geom_point(data = tab2p[!(tab2p$regulated),], mapping = aes(x=coef_capped, y=-log10(FDR_pho_kin)), color = "grey", stroke = 0, alpha = 0.2, shape = 16)
  p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
  p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
  p = p + theme_bw() + theme_nogrid()
  # p = p + geom_text_repel(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y=-log10(FDR_pho_kin),
  #                             label= as.character(pair),
  #                             color = Cancer), size= 1, alpha = 0.7, force = 1)
  p = p + scale_color_manual(values = color_cancers2)
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(panel.spacing.y=unit(0, "lines"), strip.text.y = element_text(size = 12))
  p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
  p
  subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
  dir.create(subdir2)
  fn = paste(subdir2, "volcano_trans_", "cancers" , "_", enzyme_type, '_substrate_textcutoff', text_cutoff, '.pdf',sep ="")
  ggsave(file=fn, height=4, width=5, useDingbats=FALSE)

}
