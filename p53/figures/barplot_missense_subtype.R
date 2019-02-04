# Yige Wu @ WashU 2019 Jan
## test mutation impact on protein/phosphorylation within kinase-substrate pairs or protein complex pairs



# source ------------------------------------------------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/p53/TP53_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")

# set variables -----------------------------------------------------------
gene <- "TP53"
## the position to draw
position_range <- c(170, 300)

# inputs ------------------------------------------------------------------

## input by cancer
for (cancer in c("BRCA")) {
  ## input MAF
  maf <- loadMaf(cancer = cancer, maf_files = maf_files)
  maf <- maf[maf$Hugo_Symbol == gene,]
  maf$position <- str_split_fixed(string = str_split_fixed(string = maf$HGVSp_Short, pattern = 'p.[A-Z]', 2)[,2], pattern = '\\*|[A-Z]', n = 2)[,1]
  maf$position <- str_split_fixed(string = maf$position, pattern = '\\_', n = 2)[,1]
  maf$position <- str_split_fixed(string = maf$position, pattern = '[a-z]', n = 2)[,1]
  maf$position <- as.numeric(as.vector(maf$position))
  maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  
  ## only take missenses
  maf <- maf[maf$Variant_Classification == "Missense_Mutation",]
  maf$cancer <- cancer
  sup_tab <-  maf[, c("position", "cancer", "partID")]
  sup_tab$subtype <- partID2pam50(patientID_vector = sup_tab$partID, pam50_map = loadPAM50Map())

  tab2p %>% head()
  tab2p <- sup_tab[, c("position", "subtype")]
  tab2p <- tab2p[tab2p$position >= position_range[1] & tab2p$position <= position_range[2],]
  p <- ggplot()
  p <- p + geom_histogram(data=tab2p, aes(x = position, fill = subtype, group = subtype), position='stack', color = NA, binwidth = 1)
  p <- p + geom_vline(yintercept = c(175, 248, 273), linetype = 2, color = "grey50", alpha = 0.5)
  p <- p + coord_flip()
  p <- p + scale_x_continuous(breaks = c(175, 248,273))
  # p <- p + scale_fill_manual(values = color_cat_man)
  p <- p + theme_bw()
  p <- p + theme_nogrid()
  p <- p + xlab('position of TP53 missense')+ylab("number of missenses")
  p <- p + theme(axis.title=element_text(size=10))
  p <- p + theme( axis.text.y = element_text(colour="black", size=10), 
                  panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 5))
  p
  fn = paste0(makeOutDir(resultD = resultD), cancer, "_missense_barplot.pdf")
  ggsave(file=fn, height=5, width=5)
  
}
