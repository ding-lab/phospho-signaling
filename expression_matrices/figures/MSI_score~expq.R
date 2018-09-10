# Yige Wu @ WashU March 2018
# plot MSI score against expression quantile for MMR genes

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')


# input -------------------------------------------------------------------
## input MMR genes
mmr_gene <- read_delim("~/Box Sync/MSI_CPTAC/Data/mmr_gene.txt","\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
mmr_gene <- mmr_gene$X1

## input MSI scores
msi_score <- read_delim("~/Box Sync/MSI_CPTAC/Data/MSIsensor_Score_qgao/CPTAC.MSI.score.tsv",
                        "\t", escape_double = FALSE, trim_ws = TRUE)

## input clinical info
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180307.txt"), sep = "\t")
clinical <- data.frame(clinical)

# UCEC -----------------------------------------------------------
## input expression quantile and take UCEC samples and intent gene set
cancer <- "UCEC"
expq <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_expression_quantile/cptac3_protein_quantile.RDS"))
expq <- expq[[cancer]]
expq <- expq[expq$Gene %in% mmr_gene,]

## transform sample IDs to patient IDs
expq$Participant.ID <- sampID2partID(sampleID_vector = expq$Specimen.ID, sample_map = clinical)

## merge expq with MSI score
msi_expq <- merge(msi_score, expq, by.x = c("Sample"), by.y = c("Participant.ID"))
msi_expq$MSI_status <- ifelse(msi_expq$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))

## scatterplots
library(ggrepel)
p <- ggplot(data = msi_expq)
p = p + geom_vline(xintercept = 3.5, color = 'grey')
p <- p + geom_point(mapping = aes(x = Score, y = qt, color = Gene), alpha=0.5, stroke = 0)
p = p + geom_text_repel(mapping = aes(x = Score, y = qt, label= ifelse(Score > 3.5 & qt < 0.25, as.character(Gene), NA), color = Gene),
                        force = 1, segment.color = '#cccccc', segment.size = 0.5, segment.alpha = 0.2, 
                        size=1.5,alpha=0.8)
p <- p + theme_nogrid()
p = p + theme(legend.position = "none")
p
resultDnow <- makeOutDir()
fn = paste0(resultDnow, 'MSI_score~protein_quantile_in_', cancer, ".pdf")
ggsave(file=fn, height=4, width=5, useDingbats=FALSE)

## barplot
msi_expq$QT_rank <- .bincode(msi_expq$qt, c(0,0.25,0.5,0.75,1))
msi_expq$QT_range <- "NA"
msi_expq$QT_range[msi_expq$QT_rank==1] <- "<25%"
msi_expq$QT_range[msi_expq$QT_rank==2] <- "25%-50%"
msi_expq$QT_range[msi_expq$QT_rank==3] <- "50%-75%"
msi_expq$QT_range[msi_expq$QT_rank==4] <- ">75%"
df <- data.frame(table(msi_expq[, c("Gene", "QT_range", "MSI_status")]))
df <- df[df$Freq > 0,]
df$QT_range <- factor(df$QT_range, levels = c("NA", ">75%", "50%-75%","25%-50%","<25%") )

p <- ggplot(data = df)
p <- p + geom_bar(mapping = aes(x = Gene, y = Freq, fill = QT_range, group = MSI_status), 
                         stat="identity", position='dodge')
p <- p + scale_fill_manual(values = col_paletteB(5))
p <- p + theme_bw() + theme_nogrid()
p <- p + xlab("TSGs")+ylab("number of Deletions")
p <- p + theme(axis.title=element_text(size=10))
p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), 
               axis.text.y = element_text(colour="black", size=10))
p <- p + ggtitle("XHMM TSG(+Rahman) Deletions gene expression quantile")
p
resultDnow <- makeOutDir()
fn = paste0(resultDnow, 'MSI_status~protein_quantile_in_', cancer, ".pdf")
ggsave(file=fn, height=4, width=5, useDingbats=FALSE)



fn = paste0(resultDnow, 'MSI_score~protein_quantile_fisher_test_in_', cancer, ".txt")
sink(file = fn, append=FALSE, split=FALSE)
fisher.test(table(msi_expq$Score>3.5, msi_expq$qt<0.25))
sink()


