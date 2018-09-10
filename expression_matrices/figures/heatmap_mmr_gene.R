# Yige Wu @ WashU March 2018
# plot expression heatmap for MMR genes or some other MSI genes

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/expression_matrices/figures/violin_mmr_gene.R')

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

# COAD --------------------------------------------------------------------
## input expression score(normalized log ratios) and take CO samples and intent gene set
cancer <- "CO"
exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_expression_quantile/cptac2p_protein_score.RDS"))
exps <- exps[[cancer]]
exps <- exps[exps$Gene %in% mmr_gene,]

## transform sample IDs to patient IDs
exps$Participant.ID <- sampID2partID(sampleID_vector = exps$Specimen.ID, sample_map = clinical)
## merge exps with MSI score
msi_exps <- merge(msi_score, exps, by.x = c("Sample"), by.y = c("Participant.ID"))
msi_exps$MSI_status <- ifelse(msi_exps$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
resultDnow <- makeOutDir()

myvjust = 0.7
myhjust = 0
lim = max(abs(max(msi_exps$score, na.rm = T)),abs(min(msi_exps$score, na.rm = T)))
cap <- min(2, lim)
msi_exps$score_capped <- as.numeric(msi_exps$score)
msi_exps$score_capped[msi_exps$score_capped > cap] <- cap
msi_exps$score_capped[msi_exps$score_capped < (-cap)] <- (-cap)
genes <- as.vector(msi_exps$Gene)
tmp <- unlist(mann_results[[cancer]])
sig_genes <- as.vector(str_split_fixed(string = names(tmp[tmp<0.1]), pattern = "\\.", 2)[,1])
genes[genes %in% sig_genes] <- paste0(genes[genes %in% sig_genes], "*")
msi_exps$Gene2p <- genes
tmp <- unique(msi_exps$Gene2p)
msi_exps$Gene2p <- factor(msi_exps$Gene2p, levels = c(tmp[!grepl(pattern = "MSH3", x = tmp)], tmp[grepl(pattern = "MSH3", x = tmp)]))
tmp <- msi_exps[msi_exps$Gene == "MSH3",]
msi_exps$Sample2p <- factor(msi_exps$Sample, levels = c(as.vector(tmp$Sample)[order(tmp$score)]))
p = ggplot(msi_exps)
p = p + geom_tile(aes(x=Sample2p, y=Gene2p, fill=score_capped), color=NA)#, linetype="blank") 
p = p + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p = p + scale_fill_gradientn(name= "protein abundance(log2 ratio)", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
p = p + theme_bw() + theme_nogrid()
p = p + theme(axis.title = element_blank(), axis.text.x = element_blank())
p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
p
fn = paste0(resultDnow, cancer, '_MMR_genes_protein_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
ggsave(file=fn, height=4, width=6, useDingbats=FALSE)




# UCEC --------------------------------------------------------------------
## input expression score(normalized log ratios) and take CO samples and intent gene set
cancer <- "UCEC"
exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_expression_quantile/cptac3_protein_score.RDS"))
exps <- exps[[cancer]]
exps <- exps[exps$Gene %in% mmr_gene,]

## transform sample IDs to patient IDs
exps$Participant.ID <- sampID2partID(sampleID_vector = exps$Specimen.ID, sample_map = clinical)
## merge exps with MSI score
msi_exps <- merge(msi_score, exps, by.x = c("Sample"), by.y = c("Participant.ID"))
msi_exps$MSI_status <- ifelse(msi_exps$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
resultDnow <- makeOutDir()

myvjust = 0.7
myhjust = 0
lim = max(abs(max(msi_exps$score, na.rm = T)),abs(min(msi_exps$score, na.rm = T)))
cap <- min(2, lim)
msi_exps$score_capped <- as.numeric(msi_exps$score)
msi_exps$score_capped[msi_exps$score_capped > cap] <- cap
msi_exps$score_capped[msi_exps$score_capped < (-cap)] <- (-cap)
genes <- as.vector(msi_exps$Gene)
tmp <- unlist(mann_results[[cancer]])
sig_genes <- as.vector(str_split_fixed(string = names(tmp[tmp<0.1]), pattern = "\\.", 2)[,1])
genes[genes %in% sig_genes] <- paste0(genes[genes %in% sig_genes], "*")
msi_exps$Gene2p <- genes
tmp <- unique(msi_exps$Gene2p)
msi_exps$Gene2p <- factor(msi_exps$Gene2p, levels = c(tmp[!grepl(pattern = "MLH1", x = tmp)], tmp[grepl(pattern = "MLH1", x = tmp)]))
tmp <- msi_exps[msi_exps$Gene == "MLH1",]
msi_exps$Sample2p <- factor(msi_exps$Sample, levels = c(as.vector(tmp$Sample)[order(tmp$score)]))

p = ggplot(msi_exps)
p = p + geom_tile(aes(x=Sample2p, y=Gene2p, fill=score_capped), color=NA)#, linetype="blank") 
p = p + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p = p + scale_fill_gradientn(name= "protein abundance(log2 ratio)", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
p = p + theme_bw() + theme_nogrid()
p = p + theme(axis.title = element_blank(), axis.text.x = element_blank())
p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
p
fn = paste0(resultDnow, cancer, '_MMR_genes_protein_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
ggsave(file=fn, height=4, width=6, useDingbats=FALSE)


