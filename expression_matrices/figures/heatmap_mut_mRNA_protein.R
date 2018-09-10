# Yige Wu @ WashU March 2018
# plot heatmap showing mutations, gene expression and protein expression for MMR genes

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')


# inputs ------------------------------------------------------------------
## input somatic mutations
mmr_somatic <- fread("MSI_CPTAC/Analysis/Wen-Wei/MutationTable_CPTAC.txt", data.table = F)

## input MSI scores
msi_score <- read_delim("~/Box Sync/MSI_CPTAC/Data/MSIsensor_Score_qgao/CPTAC.MSI.score.tsv",
                        "\t", escape_double = FALSE, trim_ws = TRUE)

## input UCEC clinical info with PCR results and IHC testing
ucec_clinical <- fread(input = "./MSI_CPTAC/Data/clinical_info/UCEC_report_Nov1_2017.csv", data.table = F)
colnames(ucec_clinical)[1] <- "Sample"

## input COAD clinical info with PCR results and IHC testing
co_clinical <- read_xls(path = "./MSI_CPTAC/Data/clinical_info/CPTAC Colon Clinical Data_20160915.xls")
co_clinical <- co_clinical[!is.na(co_clinical$`Microsatellite Instability (Abnormal @ >33% loci tested)`),]

## input MMR genes
mmr_gene <- read_delim("~/Box Sync/MSI_CPTAC/Data/mmr_gene.txt","\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
mmr_gene <- mmr_gene$X1

## input merged CPTAC somatic maf
test <- fread(input = "./MSI_CPTAC/Data/mergerdMMR_CPTAC.somatic.maf", data.table = F, nrows = 10)
somatic <- fread(input = "./MSI_CPTAC/Data/MAF/mergerdCPTAC.somatic.maf", data.table = F, header = F, col.names = colnames(test))
somatic <- somatic[somatic$Hugo_Symbol %in% mmr_gene, c("Hugo_Symbol", "Variant_Classification", "CaseID", "MSIsensorScore", "MSIstatus", "Cancer_Type")]
test <- fread(input = "./MSI_CPTAC/Data/mergerdMMR_CPTAC.germline.maf", data.table = F, nrows = 10)
germline <- fread(input = "./MSI_CPTAC/Data/MAF/mergerdCPTAC.germline.maf", data.table = F, header = F, col.names = colnames(test))
germline <- germline[germline$Hugo_Symbol %in% mmr_gene & germline$ExAC_AF < 0.01, c("Hugo_Symbol", "Variant_Classification", "CaseID", "MSIsensorScore", "MSIstatus", "Cancer_Type")]

# expression quantile, order by protein, UCEC --------------------------------------------------------------------
for (msi_thres in c(0.8, 3.5)) {
  for (cancer in c("UCEC")) {
    samples.msih <- as.vector(unique(msi_score$Sample[msi_score$Cancer  == cancer & msi_score$Score >= msi_thres]))
    df_somatic_m <- somatic[somatic$CaseID %in% msi_score$Sample & somatic$Cancer_Type == cancer,]
    df_somatic_m$mutation_type <- "somatic"
    df_germline_m <- germline[germline$CaseID %in% msi_score$Sample & germline$Cancer_Type == cancer,]
    df_germline_m$mutation_type <- "germline"
    df <- rbind(df_somatic_m, df_germline_m)
    colnames(df) <- c("Gene", "Variant_Classification", "Sample", "Score", "MSIstatus", "Cancer_Type", "mutation_type")
    
    ## input expression quantile
    msi_expq_list <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_expq_list.RDS"))
    df_mRNA <- msi_expq_list[[cancer]][["mRNA"]]
    df <- merge(df, df_mRNA[, c("Sample", "Gene", "qt")], by = c("Sample", "Gene"), all = T)
    df_protein <- msi_expq_list[[cancer]][["protein"]]
    df <- merge(df, df_protein[, c("Sample", "Gene", "qt")], by = c("Sample", "Gene"), all = T)
    df$MSI_status <- ifelse(df$Sample %in% samples.msih, paste0("MSI-H_tumors"), paste0("other_tumors"))
    df <- df[!(is.na(df$Variant_Classification) & is.na(df$qt.x) & is.na(df$qt.y)),]
    df <- df[!is.na(df$Gene) & !is.na(df$Sample),]
    
    ## correct MSI score
    df$Score <- NULL
    df <- merge(df, msi_score[, c("Sample", "Score")], by = c("Sample"), all.x = T)
    
    ## order the genes
    df_protein_decast <- dcast(df_protein , Gene~Sample, value.var = "qt")
    df_protein_decast_msih <- df_protein_decast[,colnames(df_protein_decast) %in% samples.msih]
    df_protein_decast_mss <- df_protein_decast[,!(colnames(df_protein_decast) %in% c("Gene", samples.msih))]
    tmp <- as.vector(df_protein_decast$Gene)[order(rowSums(df_protein_decast_msih, na.rm = T), decreasing = T)]
    tmp <- c(unique(as.vector(df$Gene))[!(unique(df$Gene) %in% tmp)], tmp)
    df$Gene <- factor(df$Gene, levels = tmp)
    
    ## order the samples
    tmp1 <- colnames(df_protein_decast_msih)[order(colSums(df_protein_decast_msih, na.rm = T), decreasing = F)]
    tmp11 <- colnames(df_protein_decast_mss)[order(colSums(df_protein_decast_mss, na.rm = T), decreasing = F)]
    tmp2 <- unique(as.vector(df$Sample))[!(unique(as.vector(df$Sample)) %in% c(tmp1, tmp11))]
    df$Sample <- factor(df$Sample, levels = c(tmp1, tmp11, tmp2))
    myvjust = 0.7
    myhjust = -0.2
    
    df_score <- unique(df[!is.na(df$Score), c("Sample", "Score", "MSI_status")])
    df_score$Sample <- factor(df_score$Sample, levels = c(tmp1, tmp11, tmp2))
    df$Variant_Classification <- factor(df$Variant_Classification, levels = c("Missense_Mutation", "Splice_Site", 
                                                                              "Frame_Shift_Del", "Frame_Shift_Ins", 
                                                                              "In_Frame_Del", "In_Frame_Ins", 
                                                                              "Nonstop_Mutation", "Nonsense_Mutation", NA))
    
    p1 = ggplot()
    p1 = p1 + geom_bar(data = df_score, mapping = aes(x = Sample, y = Score, fill = MSI_status), stat="identity")
    p1 = p1 + facet_grid(.~MSI_status, drop=T, space = "free", scales = "free")
    p1 = p1 + geom_hline(yintercept = 0.8, linetype = 2, color = "grey")
    p1 = p1 + geom_hline(yintercept = 3.5, linetype = 2, color = "grey")
    p1 = p1 + theme_bw() + theme_nogrid()
    p1 = p1 + theme(axis.text.x = element_blank())
    p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
    p1 = p1 + theme(axis.ticks=element_blank(),legend.position="right", 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p1 = p1 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    p1_zoom <- p1 + coord_cartesian(ylim = c(0, 1.5))
    
    tmp <- ucec_clinical[ucec_clinical$Sample %in% df_score$Sample, c("Sample", "MSI")]
    df_pcr <- merge(df_score, tmp, by = c("Sample"), all.x = T)
    colnames(df_pcr)[ncol(df_pcr)] <- "PCR_MSI_testing"
    df_pcr$PCR_MSI_testing[df_pcr$PCR_MSI_testing == ""] <- NA
    df_pcr$Sample <- factor(df_pcr$Sample, levels = c(tmp1, tmp11, tmp2))
    p12 <-  ggplot(df_pcr)
    p12 = p12 + geom_tile(aes(x=Sample, y="PCR_MSI_testing", fill=PCR_MSI_testing), color=NA)#, linetype="blank") 
    p12 = p12 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p12 <- p12 + scale_fill_manual(values = c("MSI indeterminate" = "grey", "MSI low" = "#FDBF6F", "MSI stable" = "#33A02C", "NA" = NA))
    p12 = p12 + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                      axis.title.y = element_blank(), axis.title.x = element_blank(),
                      strip.background = element_blank(), strip.text.x = element_blank())
    p12 = p12 + theme(axis.ticks=element_blank(),legend.position="right", 
                      legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p12 = p12 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    ## add IHC results
    tmp <- ucec_clinical[ucec_clinical$Sample %in% df_score$Sample, c("Sample", "MLH1", "MSH2", "MSH6", "PMS2", "P53")]
    tmp <- melt(tmp, id.vars =  c("Sample"))
    df_ihc <- merge(df_score, tmp, by = c("Sample"), all.x = T)
    colnames(df_ihc)[ncol(df_ihc) + (-1:0)] <- c("Gene", "IHC")
    df_ihc$IHC[df_ihc$IHC == ""] <- NA
    df_ihc <- df_ihc[!is.na(df_ihc$Gene),]
    df_ihc$Sample <- factor(df_ihc$Sample, levels = c(tmp1, tmp11, tmp2))
    p13 <-  ggplot(df_ihc)
    p13 = p13 + geom_tile(aes(x=Sample, y=Gene, fill=IHC), color=NA)
    p13 = p13 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p13 <- p13 + scale_fill_manual(values = c("Loss of nuclear expression" = "#377EB8", 
                                              "Overexpression" = "#E41A1C", 
                                              "Intact nuclear expression" = "#33A02C", 
                                              "Normal" = "#B2DF8A",
                                              "Cannot be determined" = "grey",
                                              "NA" = NA))
    p13 = p13 + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 7),
                      axis.title.y = element_blank(), axis.title.x = element_blank(),
                      strip.background = element_blank(), strip.text.x = element_blank())
    p13 = p13 + theme(axis.ticks=element_blank(),legend.position="right", 
                      legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p13 = p13 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    ## add MLH promoter methylation info
    tmp <- ucec_clinical[ucec_clinical$Sample %in% df_score$Sample, c("Sample", "MLH1_PROMOTER_HYPERMETHYLATION")]
    tmp <- melt(tmp, id.vars =  c("Sample"))
    df_mlh1_met <- merge(df_score, tmp, by = c("Sample"), all.x = T)
    colnames(df_mlh1_met)[ncol(df_mlh1_met) + (-1:0)] <- c("Gene", "mlh1_met")
    df_mlh1_met$mlh1_met[df_mlh1_met$mlh1_met == ""] <- NA
    df_mlh1_met <- df_mlh1_met[!is.na(df_mlh1_met$Gene),]
    df_mlh1_met$Sample <- factor(df_mlh1_met$Sample, levels = c(tmp1, tmp11, tmp2))
    p14 <-  ggplot(df_mlh1_met)
    p14 = p14 + geom_tile(aes(x=Sample, y=Gene, fill=mlh1_met), color=NA)
    p14 = p14 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p14 <- p14 + scale_fill_manual(values = c("Present" = "#E41A1C", 
                                              "Absent" = "#33A02C", 
                                              "Cannot be determined" = "grey",
                                              "NA" = NA))
    p14 = p14 + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 5),
                      axis.title.y = element_blank(), axis.title.x = element_blank(),
                      strip.background = element_blank(), strip.text.x = element_blank())
    p14 = p14 + theme(axis.ticks=element_blank(),legend.position="right", 
                      legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p14 = p14 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))

    ## plot mRNA
    p2 = ggplot(df)
    p2 = p2 + geom_tile(aes(x=Sample, y=Gene, fill=qt.x), color=NA)
    p2 = p2 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")
    p2 = p2 + geom_point(aes(x=Sample, y=Gene, color = mutation_type, shape = Variant_Classification), size = 1.5, alpha = 0.6)
    p2 = p2 + scale_shape_manual(values = 0:7)
    p2 <- p2 + scale_color_manual(values = c("somatic" = "#6A3D9A",
                                             "germline" = "#33A02C",
                                             "NA" = NA))
    p2 = p2 + scale_fill_gradientn(name= "mRNA abundance quantile", na.value=NA, colours=RdBu1024, limit=c(0,1))
    p2 = p2 + theme_bw() + theme_nogrid()
    p2 = p2 + theme(axis.text.x = element_blank())
    p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                    strip.background = element_blank(), strip.text.x = element_blank())
    p2 = p2 + theme(axis.ticks=element_blank(), 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p2 = p2 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    ## plot protein
    p3 = ggplot(df)
    p3 = p3 + geom_tile(aes(x=Sample, y=Gene, fill=qt.y), color=NA)#, linetype="blank") 
    p3 = p3 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p3 = p3 + geom_point(aes(x=Sample, y=Gene, color = mutation_type, shape = Variant_Classification), size = 1.5, alpha = 0.6)
    p3 = p3 + scale_shape_manual(values = 0:7)
    p3 <- p3 + scale_color_manual(values = c("somatic" = "#6A3D9A",
                                             "germline" = "#33A02C",
                                             "NA" = NA))
    p3 = p3 + scale_fill_gradientn(name= "protein abundance quantile", na.value=NA, colours=RdBu1024, limit=c(0,1))
    p3 = p3 + theme_bw() + theme_nogrid()
    p3 = p3 + theme(axis.text.x = element_text(size = 5, angle = 90, hjust = myhjust, vjust = myvjust))
    p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                    strip.background = element_blank(), strip.text.x = element_blank())
    p3 = p3 + theme(axis.ticks=element_blank(), 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p3 = p3 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    rm(msi_expq_list)
    
    resultDnow <- makeOutDir()
    gb1 <- ggplot_build(p1)
    gb2 <- ggplot_build(p2)
    gb3 <- ggplot_build(p3)
    gb1_zoom <- ggplot_build(p1_zoom)
    gb12 <- ggplot_build(p12)
    
    gA <- ggplot_gtable(gb1)
    gB <- ggplot_gtable(gb2)
    gC <- ggplot_gtable(gb3)
    gA_zoom <- ggplot_gtable(gb1_zoom)
    gA2 <- ggplot_gtable(gb12)
    
    
    g <- gtable:::rbind_gtable(gA, gA_zoom, "last")
    g <- gtable:::rbind_gtable(g, gA2, "last")
    
    gb13 <- ggplot_build(p13)
    gA3 <- ggplot_gtable(gb13)
    n13 <- length(gb13$layout$panel_ranges[[1]]$y.labels)
    g <- gtable:::rbind_gtable(g, gA3, "last")
    
    gb14 <- ggplot_build(p14)
    gA4 <- ggplot_gtable(gb14)
    n14 <- length(gb14$layout$panel_ranges[[1]]$y.labels)
    g <- gtable:::rbind_gtable(g, gA4, "last")
    
    g <- gtable:::rbind_gtable(g, gB, "last")
    g <- gtable:::rbind_gtable(g, gC, "last")
    panels <- g$layout$t[grep("panel", g$layout$name)]
    g$layout$name[grep("panel", g$layout$name)]
    
    ## adjust height of each panel
    n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
    n12 <- length(gb12$layout$panel_ranges[[1]]$y.labels)
    
    n1_zoom <- length(gb1_zoom$layout$panel_ranges[[1]]$y.labels)
    n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
    n3 <- length(gb3$layout$panel_ranges[[1]]$y.labels)
    g$heights[panels] <- unit(x = c(n1*3, n1*3, n1_zoom*5, n1_zoom*5, n12*3, n12*3, n13, n13, n14, n14,
                                    n2, n2, n3, n3), units = "null")
    fn = paste0(resultDnow, "ordered_by_protein/", cancer,  '_MMR_somatic&germline&mRNA&protein_abundance_quantile_between_MSI-H_and_rest_cutoff', msi_thres, ".pdf")
    grid.newpage()
    pdf(fn, height = nrow(df_score)/4, width = nrow(df_score)/7)
    grid.draw(g)
    dev.off()
  }
}





# expression quantile, order by protein, CO --------------------------------------------------------------------
for (msi_thres in c(0.8, 3.5)) {
  for (cancer in c("CO")) {
    samples.msih <- as.vector(unique(msi_score$Sample[msi_score$Cancer  == "COAD" & msi_score$Score >= msi_thres]))
    df_somatic_m <- somatic[somatic$CaseID %in% msi_score$Sample & somatic$Cancer_Type == "COAD",]
    df_somatic_m$mutation_type <- "somatic"
    df_germline_m <- germline[germline$CaseID %in% msi_score$Sample & germline$Cancer_Type == "COAD",]
    df_germline_m$mutation_type <- "germline"
    df <- rbind(df_somatic_m, df_germline_m)
    colnames(df) <- c("Gene", "Variant_Classification", "Sample", "Score", "MSIstatus", "Cancer_Type", "mutation_type")
    
    ## input expression quantile
    msi_expq_list <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_expq_list.RDS"))
    df_mRNA <- msi_expq_list[[cancer]][["mRNA"]]
    df <- merge(df, df_mRNA[, c("Sample", "Gene", "qt")], by = c("Sample", "Gene"), all = T)
    df_protein <- msi_expq_list[[cancer]][["protein"]]
    df <- merge(df, df_protein[, c("Sample", "Gene", "qt")], by = c("Sample", "Gene"), all = T)
    df$MSI_status <- ifelse(df$Sample %in% samples.msih, paste0("MSI-H_tumors"), paste0("other_tumors"))
    df <- df[!(is.na(df$Variant_Classification) & is.na(df$qt.x) & is.na(df$qt.y)),]
    df <- df[!is.na(df$Gene) & !is.na(df$Sample),]
    
    ## correct MSI score
    df$Score <- NULL
    df <- merge(df, msi_score[, c("Sample", "Score")], by = c("Sample"), all.x = T)
    
    ## order the genes
    df_protein_decast <- dcast(df_protein , Gene~Sample, value.var = "qt")
    df_protein_decast_msih <- df_protein_decast[,colnames(df_protein_decast) %in% samples.msih]
    df_protein_decast_mss <- df_protein_decast[,!(colnames(df_protein_decast) %in% c("Gene", samples.msih))]
    tmp <- as.vector(df_protein_decast$Gene)[order(rowSums(df_protein_decast_msih, na.rm = T), decreasing = T)]
    tmp <- c(unique(as.vector(df$Gene))[!(unique(df$Gene) %in% tmp)], tmp)
    df$Gene <- factor(df$Gene, levels = tmp)
    
    ## order the samples
    tmp1 <- colnames(df_protein_decast_msih)[order(colSums(df_protein_decast_msih, na.rm = T), decreasing = F)]
    tmp11 <- colnames(df_protein_decast_mss)[order(colSums(df_protein_decast_mss, na.rm = T), decreasing = F)]
    tmp2 <- unique(as.vector(df$Sample))[!(unique(as.vector(df$Sample)) %in% c(tmp1, tmp11))]
    df$Sample <- factor(df$Sample, levels = c(tmp1, tmp11, tmp2))
    myvjust = 0.7
    myhjust = -0.2
    
    df_score <- unique(df[!is.na(df$Score), c("Sample", "Score", "MSI_status")])
    df_score$Sample <- factor(df_score$Sample, levels = c(tmp1, tmp11, tmp2))
    df$Variant_Classification <- factor(df$Variant_Classification, levels = c("Missense_Mutation", "Splice_Site", 
                                                                              "Frame_Shift_Del", "Frame_Shift_Ins", 
                                                                              "In_Frame_Del", "In_Frame_Ins", 
                                                                              "Nonstop_Mutation", "Nonsense_Mutation", NA))
    
    p1 = ggplot()
    p1 = p1 + geom_bar(data = df_score, mapping = aes(x = Sample, y = Score, fill = MSI_status), stat="identity")
    p1 = p1 + facet_grid(.~MSI_status, drop=T, space = "free", scales = "free")
    p1 = p1 + geom_hline(yintercept = 0.8, linetype = 2, color = "grey")
    p1 = p1 + geom_hline(yintercept = 3.5, linetype = 2, color = "grey")
    p1 = p1 + theme_bw() + theme_nogrid()
    p1 = p1 + theme(axis.text.x = element_blank())
    p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
    p1 = p1 + theme(axis.ticks=element_blank(),legend.position="right", 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p1 = p1 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    p1_zoom <- p1 + coord_cartesian(ylim = c(0, 1.5))
    
    tmp <- co_clinical[co_clinical$`Participant ID` %in% df_score$Sample, c("Participant ID", "Microsatellite Instability (Abnormal @ >33% loci tested)")]
    colnames(tmp) <- c("Sample", "PCR_MSI_testing")
    df_pcr <- merge(df_score, tmp, by = c("Sample"), all.x = T)
    df_pcr$PCR_MSI_testing[df_pcr$PCR_MSI_testing == "Unknown" | df_pcr$PCR_MSI_testing == "Not Reported/ Unknown"] <- NA
    df_pcr$Sample <- factor(df_pcr$Sample, levels = c(tmp1, tmp11, tmp2))
    p12 <-  ggplot(df_pcr)
    p12 = p12 + geom_tile(aes(x=Sample, y="PCR_MSI_testing", fill=PCR_MSI_testing), color=NA)#, linetype="blank") 
    p12 = p12 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p12 <- p12 + scale_fill_manual(values = c("Yes" = "#E41A1C", "No" = "#33A02C", "Not Tested" = "grey", "NA" = NA))
    p12 = p12 + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                      axis.title.y = element_blank(), axis.title.x = element_blank(),
                      strip.background = element_blank(), strip.text.x = element_blank())
    p12 = p12 + theme(axis.ticks=element_blank(),legend.position="right", 
                      legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p12 = p12 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    ## add IHC results
    tmp <- co_clinical[co_clinical$`Participant ID` %in% df_score$Sample, c("Participant ID", "MLH1", "MSH2", "MSH6", "PMS2")]
    colnames(tmp)[1] <- "Sample"
    tmp <- data.frame(tmp)
    tmp <- melt(tmp, id.vars =  c("Sample"))
    df_ihc <- merge(df_score, tmp, by = c("Sample"), all.x = T)
    colnames(df_ihc)[ncol(df_ihc)-1] <- c("Gene")
    df_ihc$IHC <- str_split_fixed(string = df_ihc$value, pattern = "-", 2)[,2]
    df_ihc$IHC[df_ihc$IHC == ""] <- NA
    df_ihc <- df_ihc[!is.na(df_ihc$Gene),]
    df_ihc$Sample <- factor(df_ihc$Sample, levels = c(tmp1, tmp11, tmp2))
    p13 <-  ggplot(df_ihc)
    p13 = p13 + geom_tile(aes(x=Sample, y=Gene, fill=IHC), color=NA)
    p13 = p13 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p13 <- p13 + scale_fill_manual(values = c("expressed" = "#E41A1C", 
                                              "not expressed" = "#1F78B4", 
                                              "not tested" = "grey",
                                              "NA" = NA))
    p13 = p13 + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 7),
                      axis.title.y = element_blank(), axis.title.x = element_blank(),
                      strip.background = element_blank(), strip.text.x = element_blank())
    p13 = p13 + theme(axis.ticks=element_blank(),legend.position="right", 
                      legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p13 = p13 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    p2 = ggplot(df)
    p2 = p2 + geom_tile(aes(x=Sample, y=Gene, fill=qt.x), color=NA)
    p2 = p2 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")
    p2 = p2 + geom_point(aes(x=Sample, y=Gene, color = mutation_type, shape = Variant_Classification), size = 1.5, alpha = 0.6)
    p2 = p2 + scale_shape_manual(values = 0:7)
    p2 <- p2 + scale_color_manual(values = c("somatic" = "#6A3D9A",
                                             "germline" = "#33A02C",
                                             "NA" = NA))
    p2 = p2 + scale_fill_gradientn(name= "mRNA abundance quantile", na.value=NA, colours=RdBu1024, limit=c(0,1))
    p2 = p2 + theme_bw() + theme_nogrid()
    p2 = p2 + theme(axis.text.x = element_blank())
    p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                    strip.background = element_blank(), strip.text.x = element_blank())
    p2 = p2 + theme(axis.ticks=element_blank(), 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p2 = p2 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    p3 = ggplot(df)
    p3 = p3 + geom_tile(aes(x=Sample, y=Gene, fill=qt.y), color=NA)#, linetype="blank") 
    p3 = p3 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p3 = p3 + geom_point(aes(x=Sample, y=Gene, color = mutation_type, shape = Variant_Classification), size = 1.5, alpha = 0.6)
    p3 = p3 + scale_shape_manual(values = 0:7)
    p3 <- p3 + scale_color_manual(values = c("somatic" = "#6A3D9A",
                                             "germline" = "#33A02C",
                                             "NA" = NA))
    p3 = p3 + scale_fill_gradientn(name= "protein abundance quantile", na.value=NA, colours=RdBu1024, limit=c(0,1))
    p3 = p3 + theme_bw() + theme_nogrid()
    p3 = p3 + theme(axis.text.x = element_text(size = 5, angle = 90, hjust = myhjust, vjust = myvjust))
    p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                    strip.background = element_blank(), strip.text.x = element_blank())
    p3 = p3 + theme(axis.ticks=element_blank(), 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p3 = p3 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    rm(msi_expq_list)
    
    resultDnow <- makeOutDir()
    gb1 <- ggplot_build(p1)
    gb2 <- ggplot_build(p2)
    gb3 <- ggplot_build(p3)
    gb1_zoom <- ggplot_build(p1_zoom)
    gb12 <- ggplot_build(p12)
    
    gA <- ggplot_gtable(gb1)
    gB <- ggplot_gtable(gb2)
    gC <- ggplot_gtable(gb3)
    gA_zoom <- ggplot_gtable(gb1_zoom)
    gA2 <- ggplot_gtable(gb12)
    
    
    g <- gtable:::rbind_gtable(gA, gA_zoom, "last")
    g <- gtable:::rbind_gtable(g, gA2, "last")
    
    gb13 <- ggplot_build(p13)
    gA3 <- ggplot_gtable(gb13)
    n13 <- length(gb13$layout$panel_ranges[[1]]$y.labels)
    g <- gtable:::rbind_gtable(g, gA3, "last")
    
    g <- gtable:::rbind_gtable(g, gB, "last")
    g <- gtable:::rbind_gtable(g, gC, "last")
    panels <- g$layout$t[grep("panel", g$layout$name)]
    g$layout$name[grep("panel", g$layout$name)]
    
    ## adjust height of each panel
    n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
    n12 <- length(gb12$layout$panel_ranges[[1]]$y.labels)
    
    n1_zoom <- length(gb1_zoom$layout$panel_ranges[[1]]$y.labels)
    n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
    n3 <- length(gb3$layout$panel_ranges[[1]]$y.labels)
    g$heights[panels] <- unit(x = c(n1*3, n1*3, n1_zoom*5, n1_zoom*5, n12*3, n12*3, n13, n13, n2, n2, n3, n3), units = "null")
    fn = paste0(resultDnow, "ordered_by_protein/", cancer,  '_MMR_somatic&mRNA&protein_abundance_quantile_between_MSI-H_and_rest_cutoff', msi_thres, ".pdf")
    grid.newpage()
    pdf(fn, height = nrow(df_score)/6, width = nrow(df_score)/7)
    grid.draw(g)
    dev.off()
  }
}

# expression score, order by protein, UCEC --------------------------------------------------------------------
for (msi_thres in c(0.8, 3.5)) {
  for (cancer in c("UCEC")) {
    samples.msih <- as.vector(unique(msi_score$Sample[msi_score$Cancer  == cancer & msi_score$Score >= msi_thres]))
    df_somatic_m <- somatic[somatic$CaseID %in% msi_score$Sample & somatic$Cancer_Type == cancer,]
    df_somatic_m$mutation_type <- "somatic"
    df_germline_m <- germline[germline$CaseID %in% msi_score$Sample & germline$Cancer_Type == cancer,]
    df_germline_m$mutation_type <- "germline"
    df <- rbind(df_somatic_m, df_germline_m)
    colnames(df) <- c("Gene", "Variant_Classification", "Sample", "Score", "MSIstatus", "Cancer_Type", "mutation_type")
    
    ## input expression quantile
    msi_exps_list <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_exps_list.RDS"))
    df_mRNA <- msi_exps_list[[cancer]][["mRNA"]]
    cap_mRNA <- 5
    df_mRNA$score_capped <- log2(df_mRNA$score + 0.001)
    df_mRNA$score_capped[df_mRNA$score_capped > cap_mRNA] <- cap_mRNA
    df_mRNA$score_capped[df_mRNA$score_capped < (-cap_mRNA)] <- (-cap_mRNA)
    df <- merge(df, df_mRNA[, c("Sample", "Gene", "score_capped")], by = c("Sample", "Gene"), all = T)
    df_protein <- msi_exps_list[[cancer]][["protein"]]
    cap_protein <- 2
    df_protein$score_capped <- df_protein$score
    df_protein$score_capped[df_protein$score_capped > cap_protein] <- cap_protein
    df_protein$score_capped[df_protein$score_capped < (-cap_protein)] <- (-cap_protein)
    df <- merge(df, df_protein[, c("Sample", "Gene", "score_capped")], by = c("Sample", "Gene"), all = T)
    df$MSI_status <- ifelse(df$Sample %in% samples.msih, paste0("MSI-H_tumors"), paste0("other_tumors"))
    df <- df[!(is.na(df$Variant_Classification) & is.na(df$score_capped.x) & is.na(df$score_capped.y)),]
    df <- df[!is.na(df$Gene) & !is.na(df$Sample),]
    
    ## correct MSI score
    df$Score <- NULL
    df <- merge(df, msi_score[, c("Sample", "Score")], by = c("Sample"), all.x = T)
    
    ## order the genes
    df_protein_decast <- dcast(df_protein , Gene~Sample, value.var = "score_capped")
    df_protein_decast_msih <- df_protein_decast[,colnames(df_protein_decast) %in% samples.msih]
    df_protein_decast_mss <- df_protein_decast[,!(colnames(df_protein_decast) %in% c("Gene", samples.msih))]
    tmp <- as.vector(df_protein_decast$Gene)[order(rowSums(df_protein_decast_msih, na.rm = T), decreasing = T)]
    tmp <- c(unique(as.vector(df$Gene))[!(unique(df$Gene) %in% tmp)], tmp)
    df$Gene <- factor(df$Gene, levels = tmp)
    
    ## order the samples
    tmp1 <- colnames(df_protein_decast_msih)[order(colSums(df_protein_decast_msih, na.rm = T), decreasing = F)]
    tmp11 <- colnames(df_protein_decast_mss)[order(colSums(df_protein_decast_mss, na.rm = T), decreasing = F)]
    tmp2 <- unique(as.vector(df$Sample))[!(unique(as.vector(df$Sample)) %in% c(tmp1, tmp11))]
    df$Sample <- factor(df$Sample, levels = c(tmp1, tmp11, tmp2))
    myvjust = 0.7
    myhjust = -0.2
    
    df_score <- unique(df[!is.na(df$Score), c("Sample", "Score", "MSI_status")])
    df_score$Sample <- factor(df_score$Sample, levels = c(tmp1, tmp11, tmp2))
    df$Variant_Classification <- factor(df$Variant_Classification, levels = c("Missense_Mutation", "Splice_Site", 
                                                                              "Frame_Shift_Del", "Frame_Shift_Ins", 
                                                                              "In_Frame_Del", "In_Frame_Ins", 
                                                                              "Nonstop_Mutation", "Nonsense_Mutation", NA))
    
    p1 = ggplot()
    p1 = p1 + geom_bar(data = df_score, mapping = aes(x = Sample, y = Score, fill = MSI_status), stat="identity")
    p1 = p1 + facet_grid(.~MSI_status, drop=T, space = "free", scales = "free")
    p1 = p1 + geom_hline(yintercept = 0.8, linetype = 2, color = "grey")
    p1 = p1 + geom_hline(yintercept = 3.5, linetype = 2, color = "grey")
    p1 = p1 + theme_bw() + theme_nogrid()
    p1 = p1 + theme(axis.text.x = element_blank())
    p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
    p1 = p1 + theme(axis.ticks=element_blank(),legend.position="right", 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p1 = p1 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    p1_zoom <- p1 + coord_cartesian(ylim = c(0, 1.5))
    
    tmp <- ucec_clinical[ucec_clinical$Sample %in% df_score$Sample, c("Sample", "MSI")]
    df_pcr <- merge(df_score, tmp, by = c("Sample"), all.x = T)
    colnames(df_pcr)[ncol(df_pcr)] <- "PCR_MSI_testing"
    df_pcr$PCR_MSI_testing[df_pcr$PCR_MSI_testing == ""] <- NA
    df_pcr$Sample <- factor(df_pcr$Sample, levels = c(tmp1, tmp11, tmp2))
    p12 <-  ggplot(df_pcr)
    p12 = p12 + geom_tile(aes(x=Sample, y="PCR_MSI_testing", fill=PCR_MSI_testing), color=NA)#, linetype="blank") 
    p12 = p12 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p12 <- p12 + scale_fill_manual(values = c("MSI indeterminate" = "grey", "MSI low" = "#FDBF6F", "MSI stable" = "#33A02C", "NA" = NA))
    p12 = p12 + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                      axis.title.y = element_blank(), axis.title.x = element_blank(),
                      strip.background = element_blank(), strip.text.x = element_blank())
    p12 = p12 + theme(axis.ticks=element_blank(),legend.position="right", 
                      legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p12 = p12 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    ## add IHC results
    tmp <- ucec_clinical[ucec_clinical$Sample %in% df_score$Sample, c("Sample", "MLH1", "MSH2", "MSH6", "PMS2", "P53")]
    tmp <- melt(tmp, id.vars =  c("Sample"))
    df_ihc <- merge(df_score, tmp, by = c("Sample"), all.x = T)
    colnames(df_ihc)[ncol(df_ihc) + (-1:0)] <- c("Gene", "IHC")
    df_ihc$IHC[df_ihc$IHC == ""] <- NA
    df_ihc <- df_ihc[!is.na(df_ihc$Gene),]
    df_ihc$Sample <- factor(df_ihc$Sample, levels = c(tmp1, tmp11, tmp2))
    p13 <-  ggplot(df_ihc)
    p13 = p13 + geom_tile(aes(x=Sample, y=Gene, fill=IHC), color=NA)
    p13 = p13 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p13 <- p13 + scale_fill_manual(values = c("Loss of nuclear expression" = "#377EB8", 
                                              "Overexpression" = "#E41A1C", 
                                              "Intact nuclear expression" = "#33A02C", 
                                              "Normal" = "#B2DF8A",
                                              "Cannot be determined" = "grey",
                                              "NA" = NA))
    p13 = p13 + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 7),
                      axis.title.y = element_blank(), axis.title.x = element_blank(),
                      strip.background = element_blank(), strip.text.x = element_blank())
    p13 = p13 + theme(axis.ticks=element_blank(),legend.position="right", 
                      legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p13 = p13 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    ## add MLH promoter methylation info
    tmp <- ucec_clinical[ucec_clinical$Sample %in% df_score$Sample, c("Sample", "MLH1_PROMOTER_HYPERMETHYLATION")]
    tmp <- melt(tmp, id.vars =  c("Sample"))
    df_mlh1_met <- merge(df_score, tmp, by = c("Sample"), all.x = T)
    colnames(df_mlh1_met)[ncol(df_mlh1_met) + (-1:0)] <- c("Gene", "mlh1_met")
    df_mlh1_met$mlh1_met[df_mlh1_met$mlh1_met == ""] <- NA
    df_mlh1_met <- df_mlh1_met[!is.na(df_mlh1_met$Gene),]
    df_mlh1_met$Sample <- factor(df_mlh1_met$Sample, levels = c(tmp1, tmp11, tmp2))
    p14 <-  ggplot(df_mlh1_met)
    p14 = p14 + geom_tile(aes(x=Sample, y=Gene, fill=mlh1_met), color=NA)
    p14 = p14 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p14 <- p14 + scale_fill_manual(values = c("Present" = "#E41A1C", 
                                              "Absent" = "#33A02C", 
                                              "Cannot be determined" = "grey",
                                              "NA" = NA))
    p14 = p14 + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 5),
                      axis.title.y = element_blank(), axis.title.x = element_blank(),
                      strip.background = element_blank(), strip.text.x = element_blank())
    p14 = p14 + theme(axis.ticks=element_blank(),legend.position="right", 
                      legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p14 = p14 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    df$shape <- factor(mod(e1 = 1:nrow(df), e2 = 10))
    p2 = ggplot(df)
    p2 = p2 + geom_tile(aes(x=Sample, y=Gene, fill=score_capped.x), color=NA)
    p2 = p2 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")
    p2 = p2 + geom_point(data = df, aes(x=Sample, y=Gene, color = mutation_type, shape = Variant_Classification), size = 1.5, alpha = 0.6)
    p2 = p2 + scale_shape_manual(values = 0:7)
    p2 <- p2 + scale_color_manual(values = c("somatic" = "#6A3D9A",
                                             "germline" = "#33A02C",
                                             "NA" = NA))
    p2 = p2 + scale_fill_gradientn(name= "mRNA abundance quantile", na.value=NA, colours=RdBu1024, limit=c(-cap_mRNA,cap_mRNA))
    p2 = p2 + theme_bw() + theme_nogrid()
    p2 = p2 + theme(axis.text.x = element_blank())
    p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                    strip.background = element_blank(), strip.text.x = element_blank())
    p2 = p2 + theme(axis.ticks=element_blank(), 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p2 = p2 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    p3 = ggplot(df)
    p3 = p3 + geom_tile(aes(x=Sample, y=Gene, fill=score_capped.y), color=NA)#, linetype="blank") 
    p3 = p3 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p3 = p3 + geom_point(aes(x=Sample, y=Gene, color = mutation_type, shape = Variant_Classification), size = 1.5, alpha = 0.6)
    p3 = p3 + scale_shape_manual(values = 0:7)
    p3 <- p3 + scale_color_manual(values = c("somatic" = "#6A3D9A",
                                             "germline" = "#33A02C",
                                             "NA" = NA))
    p3 = p3 + scale_fill_gradientn(name= "protein abundance quantile", na.value=NA, colours=RdBu1024, limit=c(-cap_protein, cap_protein))
    p3 = p3 + theme_bw() + theme_nogrid()
    p3 = p3 + theme(axis.text.x = element_text(size = 5, angle = 90, hjust = myhjust, vjust = myvjust))
    p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                    strip.background = element_blank(), strip.text.x = element_blank())
    p3 = p3 + theme(axis.ticks=element_blank(), 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p3 = p3 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    rm(msi_expq_list)
    
    resultDnow <- makeOutDir()
    gb1 <- ggplot_build(p1)
    gb2 <- ggplot_build(p2)
    gb3 <- ggplot_build(p3)
    gb1_zoom <- ggplot_build(p1_zoom)
    gb12 <- ggplot_build(p12)
    
    gA <- ggplot_gtable(gb1)
    gB <- ggplot_gtable(gb2)
    gC <- ggplot_gtable(gb3)
    gA_zoom <- ggplot_gtable(gb1_zoom)
    gA2 <- ggplot_gtable(gb12)
    
    
    g <- gtable:::rbind_gtable(gA, gA_zoom, "last")
    g <- gtable:::rbind_gtable(g, gA2, "last")
    
    gb13 <- ggplot_build(p13)
    gA3 <- ggplot_gtable(gb13)
    n13 <- length(gb13$layout$panel_ranges[[1]]$y.labels)
    g <- gtable:::rbind_gtable(g, gA3, "last")
    
    gb14 <- ggplot_build(p14)
    gA4 <- ggplot_gtable(gb14)
    n14 <- length(gb14$layout$panel_ranges[[1]]$y.labels)
    g <- gtable:::rbind_gtable(g, gA4, "last")
    
    g <- gtable:::rbind_gtable(g, gB, "last")
    g <- gtable:::rbind_gtable(g, gC, "last")
    panels <- g$layout$t[grep("panel", g$layout$name)]
    g$layout$name[grep("panel", g$layout$name)]
    
    ## adjust height of each panel
    n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
    n12 <- length(gb12$layout$panel_ranges[[1]]$y.labels)
    
    n1_zoom <- length(gb1_zoom$layout$panel_ranges[[1]]$y.labels)
    n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
    n3 <- length(gb3$layout$panel_ranges[[1]]$y.labels)
    g$heights[panels] <- unit(x = c(n1*3, n1*3, n1_zoom*5, n1_zoom*5, n12*3, n12*3, n13, n13, n14, n14,
                                    n2, n2, n3, n3), units = "null")
    fn = paste0(resultDnow, "ordered_by_protein/", cancer,  '_MMR_somatic&germline&mRNA&protein_abundance_score_between_MSI-H_and_rest_cutoff', msi_thres, ".pdf")
    grid.newpage()
    pdf(fn, height = nrow(df_score)/4, width = nrow(df_score)/7)
    grid.draw(g)
    dev.off()
  }
}





# expression score, order by protein, CO --------------------------------------------------------------------
for (msi_thres in c(0.8, 3.5)) {
  for (cancer in c("CO")) {
    samples.msih <- as.vector(unique(msi_score$Sample[msi_score$Cancer  == "COAD" & msi_score$Score >= msi_thres]))
    df_somatic_m <- somatic[somatic$CaseID %in% msi_score$Sample & somatic$Cancer_Type == "COAD",]
    df_somatic_m$mutation_type <- "somatic"
    df_germline_m <- germline[germline$CaseID %in% msi_score$Sample & germline$Cancer_Type == "COAD",]
    df_germline_m$mutation_type <- "germline"
    df <- rbind(df_somatic_m, df_germline_m)
    colnames(df) <- c("Gene", "Variant_Classification", "Sample", "Score", "MSIstatus", "Cancer_Type", "mutation_type")
    
    ## input expression quantile
    msi_exps_list <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_exps_list.RDS"))
    df_mRNA <- msi_exps_list[[cancer]][["mRNA"]]
    cap_mRNA <- 5
    df_mRNA$score_capped <- log2(df_mRNA$score + 0.001)
    df_mRNA$score_capped[df_mRNA$score_capped > cap_mRNA] <- cap_mRNA
    df_mRNA$score_capped[df_mRNA$score_capped < (-cap_mRNA)] <- (-cap_mRNA)
    df <- merge(df, df_mRNA[, c("Sample", "Gene", "score_capped")], by = c("Sample", "Gene"), all = T)
    df_protein <- msi_exps_list[[cancer]][["protein"]]
    cap_protein <- 2
    df_protein$score_capped <- df_protein$score
    df_protein$score_capped[df_protein$score_capped > cap_protein] <- cap_protein
    df_protein$score_capped[df_protein$score_capped < (-cap_protein)] <- (-cap_protein)
    df <- merge(df, df_protein[, c("Sample", "Gene", "score_capped")], by = c("Sample", "Gene"), all = T)
    df$MSI_status <- ifelse(df$Sample %in% samples.msih, paste0("MSI-H_tumors"), paste0("other_tumors"))
    df <- df[!(is.na(df$Variant_Classification) & is.na(df$score_capped.x) & is.na(df$score_capped.y)),]
    df <- df[!is.na(df$Gene) & !is.na(df$Sample),]
    
    ## correct MSI score
    df$Score <- NULL
    df <- merge(df, msi_score[, c("Sample", "Score")], by = c("Sample"), all.x = T)
    
    ## order the genes
    df_protein_decast <- dcast(df_protein , Gene~Sample, value.var = "score_capped")
    df_protein_decast_msih <- df_protein_decast[,colnames(df_protein_decast) %in% samples.msih]
    df_protein_decast_mss <- df_protein_decast[,!(colnames(df_protein_decast) %in% c("Gene", samples.msih))]
    tmp <- as.vector(df_protein_decast$Gene)[order(rowSums(df_protein_decast_msih, na.rm = T), decreasing = T)]
    tmp <- c(unique(as.vector(df$Gene))[!(unique(df$Gene) %in% tmp)], tmp)
    df$Gene <- factor(df$Gene, levels = tmp)
    
    ## order the samples
    tmp1 <- colnames(df_protein_decast_msih)[order(colSums(df_protein_decast_msih, na.rm = T), decreasing = F)]
    tmp11 <- colnames(df_protein_decast_mss)[order(colSums(df_protein_decast_mss, na.rm = T), decreasing = F)]
    tmp2 <- unique(as.vector(df$Sample))[!(unique(as.vector(df$Sample)) %in% c(tmp1, tmp11))]
    df$Sample <- factor(df$Sample, levels = c(tmp1, tmp11, tmp2))
    myvjust = 0.7
    myhjust = -0.2
    
    df_score <- unique(df[!is.na(df$Score), c("Sample", "Score", "MSI_status")])
    df_score$Sample <- factor(df_score$Sample, levels = c(tmp1, tmp11, tmp2))
    df$Variant_Classification <- factor(df$Variant_Classification, levels = c("Missense_Mutation", "Splice_Site", 
                                                                              "Frame_Shift_Del", "Frame_Shift_Ins", 
                                                                              "In_Frame_Del", "In_Frame_Ins", 
                                                                              "Nonstop_Mutation", "Nonsense_Mutation", NA))
    p1 = ggplot()
    p1 = p1 + geom_bar(data = df_score, mapping = aes(x = Sample, y = Score, fill = MSI_status), stat="identity")
    p1 = p1 + facet_grid(.~MSI_status, drop=T, space = "free", scales = "free")
    p1 = p1 + geom_hline(yintercept = 0.8, linetype = 2, color = "grey")
    p1 = p1 + geom_hline(yintercept = 3.5, linetype = 2, color = "grey")
    p1 = p1 + theme_bw() + theme_nogrid()
    p1 = p1 + theme(axis.text.x = element_blank())
    p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
    p1 = p1 + theme(axis.ticks=element_blank(),legend.position="right", 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p1 = p1 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    p1_zoom <- p1 + coord_cartesian(ylim = c(0, 1.5))
    
    tmp <- co_clinical[co_clinical$`Participant ID` %in% df_score$Sample, c("Participant ID", "Microsatellite Instability (Abnormal @ >33% loci tested)")]
    colnames(tmp) <- c("Sample", "PCR_MSI_testing")
    df_pcr <- merge(df_score, tmp, by = c("Sample"), all.x = T)
    df_pcr$PCR_MSI_testing[df_pcr$PCR_MSI_testing == "Unknown" | df_pcr$PCR_MSI_testing == "Not Reported/ Unknown"] <- NA
    df_pcr$Sample <- factor(df_pcr$Sample, levels = c(tmp1, tmp11, tmp2))
    p12 <-  ggplot(df_pcr)
    p12 = p12 + geom_tile(aes(x=Sample, y="PCR_MSI_testing", fill=PCR_MSI_testing), color=NA)#, linetype="blank") 
    p12 = p12 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p12 <- p12 + scale_fill_manual(values = c("Yes" = "#E41A1C", "No" = "#33A02C", "Not Tested" = "grey", "NA" = NA))
    p12 = p12 + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                      axis.title.y = element_blank(), axis.title.x = element_blank(),
                      strip.background = element_blank(), strip.text.x = element_blank())
    p12 = p12 + theme(axis.ticks=element_blank(),legend.position="right", 
                      legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p12 = p12 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    ## add IHC results
    tmp <- co_clinical[co_clinical$`Participant ID` %in% df_score$Sample, c("Participant ID", "MLH1", "MSH2", "MSH6", "PMS2")]
    colnames(tmp)[1] <- "Sample"
    tmp <- data.frame(tmp)
    tmp <- melt(tmp, id.vars =  c("Sample"))
    df_ihc <- merge(df_score, tmp, by = c("Sample"), all.x = T)
    colnames(df_ihc)[ncol(df_ihc)-1] <- c("Gene")
    df_ihc$IHC <- str_split_fixed(string = df_ihc$value, pattern = "-", 2)[,2]
    df_ihc$IHC[df_ihc$IHC == ""] <- NA
    df_ihc <- df_ihc[!is.na(df_ihc$Gene),]
    df_ihc$Sample <- factor(df_ihc$Sample, levels = c(tmp1, tmp11, tmp2))
    p13 <-  ggplot(df_ihc)
    p13 = p13 + geom_tile(aes(x=Sample, y=Gene, fill=IHC), color=NA)
    p13 = p13 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p13 <- p13 + scale_fill_manual(values = c("expressed" = "#E41A1C", 
                                              "not expressed" = "#1F78B4", 
                                              "not tested" = "grey",
                                              "NA" = NA))
    p13 = p13 + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 7),
                      axis.title.y = element_blank(), axis.title.x = element_blank(),
                      strip.background = element_blank(), strip.text.x = element_blank())
    p13 = p13 + theme(axis.ticks=element_blank(),legend.position="right", 
                      legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p13 = p13 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    p2 = ggplot(df)
    p2 = p2 + geom_tile(aes(x=Sample, y=Gene, fill=score_capped.x), color=NA)
    p2 = p2 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")
    p2 = p2 + geom_point(aes(x=Sample, y=Gene, color = mutation_type, shape = Variant_Classification), size = 1.5, alpha = 0.6)
    p2 = p2 + scale_shape_manual(values = 0:7)
    p2 <- p2 + scale_color_manual(values = c("somatic" = "#6A3D9A",
                                             "germline" = "#33A02C",
                                             "NA" = NA))
    p2 = p2 + scale_fill_gradientn(name= "mRNA abundance quantile", na.value=NA, colours=RdBu1024, limit=c(-cap_mRNA,cap_mRNA))
    p2 = p2 + theme_bw() + theme_nogrid()
    p2 = p2 + theme(axis.text.x = element_blank())
    p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                    strip.background = element_blank(), strip.text.x = element_blank())
    p2 = p2 + theme(axis.ticks=element_blank(), 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p2 = p2 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    p3 = ggplot(df)
    p3 = p3 + geom_tile(aes(x=Sample, y=Gene, fill=score_capped.y), color=NA)#, linetype="blank") 
    p3 = p3 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p3 = p3 + geom_point(aes(x=Sample, y=Gene, color = mutation_type, shape = Variant_Classification), size = 1.5, alpha = 0.6)
    p3 = p3 + scale_shape_manual(values = 0:7)
    p3 <- p3 + scale_color_manual(values = c("somatic" = "#6A3D9A",
                                             "germline" = "#33A02C",
                                             "NA" = NA))
    p3 = p3 + scale_fill_gradientn(name= "protein abundance quantile", na.value=NA, colours=RdBu1024, limit=c(-cap_protein, cap_protein))
    p3 = p3 + theme_bw() + theme_nogrid()
    p3 = p3 + theme(axis.text.x = element_text(size = 5, angle = 90, hjust = myhjust, vjust = myvjust))
    p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                    strip.background = element_blank(), strip.text.x = element_blank())
    p3 = p3 + theme(axis.ticks=element_blank(), 
                    legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
    p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
    p3 = p3 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
    
    rm(msi_expq_list)
    
    resultDnow <- makeOutDir()
    gb1 <- ggplot_build(p1)
    gb2 <- ggplot_build(p2)
    gb3 <- ggplot_build(p3)
    gb1_zoom <- ggplot_build(p1_zoom)
    gb12 <- ggplot_build(p12)
    
    gA <- ggplot_gtable(gb1)
    gB <- ggplot_gtable(gb2)
    gC <- ggplot_gtable(gb3)
    gA_zoom <- ggplot_gtable(gb1_zoom)
    gA2 <- ggplot_gtable(gb12)
    
    
    g <- gtable:::rbind_gtable(gA, gA_zoom, "last")
    g <- gtable:::rbind_gtable(g, gA2, "last")
    
    gb13 <- ggplot_build(p13)
    gA3 <- ggplot_gtable(gb13)
    n13 <- length(gb13$layout$panel_ranges[[1]]$y.labels)
    g <- gtable:::rbind_gtable(g, gA3, "last")
    
    g <- gtable:::rbind_gtable(g, gB, "last")
    g <- gtable:::rbind_gtable(g, gC, "last")
    panels <- g$layout$t[grep("panel", g$layout$name)]
    g$layout$name[grep("panel", g$layout$name)]
    
    ## adjust height of each panel
    n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
    n12 <- length(gb12$layout$panel_ranges[[1]]$y.labels)
    
    n1_zoom <- length(gb1_zoom$layout$panel_ranges[[1]]$y.labels)
    n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
    n3 <- length(gb3$layout$panel_ranges[[1]]$y.labels)
    g$heights[panels] <- unit(x = c(n1*3, n1*3, n1_zoom*5, n1_zoom*5, n12*3, n12*3, n13, n13, n2, n2, n3, n3), units = "null")
    fn = paste0(resultDnow, "ordered_by_protein/", cancer,  '_MMR_somatic&germline&mRNA&protein_abundance_score_between_MSI-H_and_rest_cutoff', msi_thres, ".pdf")
    grid.newpage()
    pdf(fn, height = nrow(df_score)/6, width = nrow(df_score)/7)
    grid.draw(g)
    dev.off()
  }
}


# expression quantile, order by mRNA --------------------------------------
for (cancer in c( "CO", "UCEC")) {
  df_somatic_m <- melt(mmr_somatic, id.vars = c("Hugo_Symbol"))
  colnames(df_somatic_m) <- c("Gene", "Sample", "somatic_tier1_count")
  if (cancer == "CO") {
    df_somatic_m <- df_somatic_m[grepl(x = as.vector(df_somatic_m$Sample), pattern = "CO"),]
    df_somatic_m <- merge(df_somatic_m, msi_score[msi_score$Cancer  == "COAD",], by = c("Sample"), all = T)
    samples.msih <- as.vector(unique(msi_score$Sample[msi_score$Cancer  == "COAD" & msi_score$Score >= 3.5]))
  } else {
    df_somatic_m <- df_somatic_m[!(grepl(x = as.vector(df_somatic_m$Sample), pattern = "CO")),]
    df_somatic_m <- merge(df_somatic_m, msi_score[msi_score$Cancer  == cancer,], by = c("Sample"), all = T)
    samples.msih <- as.vector(unique(msi_score$Sample[msi_score$Cancer  == cancer & msi_score$Score >= 3.5]))
  }
  df <- df_somatic_m
  df <- df[df$Sample %in% msi_score$Sample,]
  
  ## input expression quantile
  msi_expq_list <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_expq_list.RDS"))
  df_mRNA <- msi_expq_list[[cancer]][["mRNA"]]
  df <- merge(df, df_mRNA[, c("Sample", "Gene", "qt")], by = c("Sample", "Gene"), all = T)
  df_protein <- msi_expq_list[[cancer]][["protein"]]
  df <- merge(df, df_protein[, c("Sample", "Gene", "qt")], by = c("Sample", "Gene"), all = T)
  df$MSI_status <- ifelse(df$Sample %in% samples.msih, paste0("MSI-H_tumors"), paste0("other_tumors"))
  df <- df[!(is.na(df$somatic_tier1_count) & is.na(df$qt.x) & is.na(df$qt.y)),]
  df <- df[!is.na(df$Gene) & !is.na(df$Sample),]
  
  ## order the genes
  df_mRNA_decast <- dcast(df_mRNA , Gene~Sample, value.var = "qt")
  df_mRNA_decast_msih <- df_mRNA_decast[,colnames(df_mRNA_decast) %in% samples.msih]
  df_mRNA_decast_mss <- df_mRNA_decast[,!(colnames(df_mRNA_decast) %in% c("Gene", samples.msih))]
  tmp <- as.vector(df_mRNA_decast$Gene)[order(rowSums(df_mRNA_decast_msih, na.rm = T), decreasing = T)]
  tmp <- c(unique(as.vector(df$Gene))[!(unique(df$Gene) %in% tmp)], tmp)
  df$Gene <- factor(df$Gene, levels = tmp)
  
  ## order the samples
  tmp1 <- colnames(df_mRNA_decast_msih)[order(colSums(df_mRNA_decast_msih, na.rm = T), decreasing = F)]
  tmp11 <- colnames(df_mRNA_decast_mss)[order(colSums(df_mRNA_decast_mss, na.rm = T), decreasing = F)]
  tmp2 <- unique(as.vector(df$Sample))[!(unique(as.vector(df$Sample)) %in% c(tmp1, tmp11))]
  df$Sample <- factor(df$Sample, levels = c(tmp1, tmp11, tmp2))
  myvjust = 0.7
  myhjust = -0.2
  
  p2 = ggplot(df)
  p2 = p2 + geom_tile(aes(x=Sample, y=Gene, fill=qt.x), color=NA)#, linetype="blank") 
  p2 = p2 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p2 = p2 + geom_text(aes(x=Sample, y=Gene, label = ifelse(somatic_tier1_count > 0 , somatic_tier1_count, NA)), size = 1.5)
  p2 = p2 + scale_fill_gradientn(name= "mRNA abundance quantile", na.value=NA, colours=RdBu1024, limit=c(0,1))
  p2 = p2 + theme_bw() + theme_nogrid()
  p2 = p2 + theme(axis.text.x = element_blank())
  p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
  p2 = p2 + theme(axis.ticks=element_blank(),legend.position="right", legend.text = element_text(size = 5))
  p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  
  p3 = ggplot(df)
  p3 = p3 + geom_tile(aes(x=Sample, y=Gene, fill=qt.y), color=NA)#, linetype="blank") 
  p3 = p3 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p3 = p3 + geom_text(aes(x=Sample, y=Gene, label = ifelse(somatic_tier1_count > 0 , somatic_tier1_count, NA)), size = 1.5)
  p3 = p3 + scale_fill_gradientn(name= "protein abundance quantile", na.value=NA, colours=RdBu1024, limit=c(0,1))
  p3 = p3 + theme_bw() + theme_nogrid()
  # p3 = p3 + theme(axis.text.x = element_blank())
  p3 = p3 + theme(axis.text.x = element_text(size = 3, angle = 90, hjust = myhjust, vjust = myvjust))
  p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
  p3 = p3 + theme(axis.ticks=element_blank(),legend.position="right", legend.text = element_text(size = 5))
  p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  
  rm(msi_expq_list)
  
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cancer, '_MMR_somatic&mRNA&protein_abundance_quantile_between_MSI-H_and_rest_in_', cancer, ".v2.pdf")
  grid.newpage()
  # pdf(fn, height = 14*(len+10)/(top+10), width = 8*(len+10)/(top+10))
  pdf(fn, height = 6, width = 7)
  plot123 <- rbind(ggplotGrob(p2), ggplotGrob(p3), size = "last")
  title <- textGrob(paste0(cancer, " somatic mutations and gene/protein abundance in MMR genes"),gp=gpar(fontsize=16))
  padding <- unit(0,"lines")
  plottgt <- gtable_add_rows(plot123, 
                             heights = grobHeight(title) + padding,
                             pos = 0)
  plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
  grid.draw(plottgt)
  dev.off()
}


# expression value, order by mRNA --------------------------------------------------------------------

for (cancer in c( "CO", "UCEC")) {
  df_somatic_m <- melt(mmr_somatic, id.vars = c("Hugo_Symbol"))
  colnames(df_somatic_m) <- c("Gene", "Sample", "somatic_tier1_count")
  if (cancer == "CO") {
    df_somatic_m <- df_somatic_m[grepl(x = as.vector(df_somatic_m$Sample), pattern = "CO"),]
    df_somatic_m <- merge(df_somatic_m, msi_score[msi_score$Cancer  == "COAD",], by = c("Sample"), all = T)
    samples.msih <- as.vector(unique(msi_score$Sample[msi_score$Cancer  == "COAD" & msi_score$Score >= 3.5]))
  } else {
    df_somatic_m <- df_somatic_m[!(grepl(x = as.vector(df_somatic_m$Sample), pattern = "CO")),]
    df_somatic_m <- merge(df_somatic_m, msi_score[msi_score$Cancer  == cancer,], by = c("Sample"), all = T)
    samples.msih <- as.vector(unique(msi_score$Sample[msi_score$Cancer  == cancer & msi_score$Score >= 3.5]))
  }
  df <- df_somatic_m
  df <- df[df$Sample %in% msi_score$Sample,]
  
  ## input expression quantile
  msi_exps_list <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_exps_list.RDS"))
  df_mRNA <- msi_exps_list[[cancer]][["mRNA"]]
  cap_mRNA <- 5
  df_mRNA$score_capped <- df_mRNA$log2RPKM
  df_mRNA$score_capped[df_mRNA$score_capped > cap_mRNA] <- cap_mRNA
  df_mRNA$score_capped[df_mRNA$score_capped < (-cap_mRNA)] <- (-cap_mRNA)
  df <- merge(df, df_mRNA[, c("Sample", "Gene", "score_capped")], by = c("Sample", "Gene"), all = T)
  df_protein <- msi_exps_list[[cancer]][["protein"]]
  cap_protein <- 2
  df_protein$score_capped <- df_protein$score
  df_protein$score_capped[df_protein$score_capped > cap_protein] <- cap_protein
  df_protein$score_capped[df_protein$score_capped < (-cap_protein)] <- (-cap_protein)
  df <- merge(df, df_protein[, c("Sample", "Gene", "score_capped")], by = c("Sample", "Gene"), all = T)
  df$MSI_status <- ifelse(df$Sample %in% samples.msih, paste0("MSI-H_tumors"), paste0("other_tumors"))
  df <- df[!(is.na(df$somatic_tier1_count) & is.na(df$score_capped.x) & is.na(df$score_capped.y)),]
  df <- df[!is.na(df$Gene) & !is.na(df$Sample),]
  
  ## order the genes
  df_mRNA_decast <- dcast(df_mRNA , Gene~Sample, value.var = "score_capped")
  df_mRNA_decast_msih <- df_mRNA_decast[,colnames(df_mRNA_decast) %in% samples.msih]
  df_mRNA_decast_mss <- df_mRNA_decast[,!(colnames(df_mRNA_decast) %in% c("Gene", samples.msih))]
  tmp <- as.vector(df_mRNA_decast$Gene)[order(rowSums(df_mRNA_decast_msih, na.rm = T), decreasing = T)]
  tmp <- c(unique(as.vector(df$Gene))[!(unique(df$Gene) %in% tmp)], tmp)
  df$Gene <- factor(df$Gene, levels = tmp)
  
  ## order the samples
  tmp1 <- colnames(df_mRNA_decast_msih)[order(colSums(df_mRNA_decast_msih, na.rm = T), decreasing = F)]
  tmp11 <- colnames(df_mRNA_decast_mss)[order(colSums(df_mRNA_decast_mss, na.rm = T), decreasing = F)]
  tmp2 <- unique(as.vector(df$Sample))[!(unique(as.vector(df$Sample)) %in% c(tmp1, tmp11))]
  df$Sample <- factor(df$Sample, levels = c(tmp1, tmp11, tmp2))
  myvjust = 0.7
  myhjust = -0.2
  
  p2 = ggplot(df)
  p2 = p2 + geom_tile(aes(x=Sample, y=Gene, fill=score_capped.x), color=NA)#, linetype="blank") 
  p2 = p2 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p2 = p2 + geom_text(aes(x=Sample, y=Gene, label = ifelse(somatic_tier1_count > 0 , somatic_tier1_count, NA)), size = 1.5)
  p2 = p2 + scale_fill_gradientn(name= "mRNA abundance (log2RPKM)", na.value=NA, colours=RdBu1024, limit=c(-cap_mRNA,cap_mRNA))
  p2 = p2 + theme_bw() + theme_nogrid()
  p2 = p2 + theme(axis.text.x = element_blank())
  p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
  p2 = p2 + theme(axis.ticks=element_blank(),legend.position="right", legend.text = element_text(size = 5))
  p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  
  p3 = ggplot(df)
  p3 = p3 + geom_tile(aes(x=Sample, y=Gene, fill=score_capped.y), color=NA)#, linetype="blank") 
  p3 = p3 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p3 = p3 + geom_text(aes(x=Sample, y=Gene, label = ifelse(somatic_tier1_count > 0 , somatic_tier1_count, NA)), size = 1.5)
  p3 = p3 + scale_fill_gradientn(name= "protein abundance(log2 ratios normalized)", na.value=NA, colours=RdBu1024, limit=c((-cap_protein),cap_protein))
  p3 = p3 + theme_bw() + theme_nogrid()
  # p3 = p3 + theme(axis.text.x = element_blank())
  p3 = p3 + theme(axis.text.x = element_text(size = 3, angle = 90, hjust = myhjust, vjust = myvjust))
  p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
  p3 = p3 + theme(axis.ticks=element_blank(),legend.position="right", legend.text = element_text(size = 5))
  p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  
  rm(msi_exps_list)
  
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cancer, '_MMR_somatic&mRNA&protein_abundance_value_between_MSI-H_and_rest_in_', cancer, ".v2.pdf")
  grid.newpage()
  # pdf(fn, height = 14*(len+10)/(top+10), width = 8*(len+10)/(top+10))
  pdf(fn, height = 6, width = 7)
  plot123 <- rbind(ggplotGrob(p2), ggplotGrob(p3), size = "last")
  title <- textGrob(paste0(cancer, " somatic mutations and gene/protein abundance in MMR genes"),gp=gpar(fontsize=16))
  padding <- unit(0,"lines")
  plottgt <- gtable_add_rows(plot123, 
                             heights = grobHeight(title) + padding,
                             pos = 0)
  plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
  grid.draw(plottgt)
  dev.off()
}



# CO ----------------------------------------------------------------------
for (cancer in c("CO")) {
  df_protein <- msi_exps_list[[cancer]][["protein"]][, c("Gene", "score", "MSI_status")] ; df_protein$expression_type <- "protein"
  tmp <- mmr_somatic[, colnames(mmr_somatic) %in% unique(msi_exps_list[[cancer]][["mRNA"]]$Sample)]
  df_somatic <- data.frame(Gene = mmr_somatic$Hugo_Symbol)
  df_somatic <- cbind(df_somatic, tmp)
  df_somatic_m <- melt(df_somatic, id.vars = c("Gene"))
  colnames(df_somatic_m) <- c("Gene", "Sample", "somatic_tier1_count")
  df_somatic_m <- merge(df_somatic_m, msi_score, by = c("Sample"), all.x = T)
  df_somatic_m$MSI_status <- ifelse(df_somatic_m$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
  # df_somatic_m <- df_somatic_m[df_somatic_m$somatic_tier1_count > 0,]
  df <- df_somatic_m
  
  p1 = ggplot(df)
  p1 = p1 + geom_tile(aes(x=Sample, y=Gene, fill=somatic_tier1_count), color=NA)#, linetype="blank") 
  p1 = p1 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p1 = p1 + scale_fill_gradientn(name= "somatic_tier1_count", na.value=NA, colours=c("white",RdBu1024[600:1024]))
  p1 = p1 + theme_bw() + theme_nogrid()
  p1 = p1 + theme(axis.text.x = element_blank())
  p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
  p1 = p1 + theme(axis.ticks=element_blank(), legend.position = "right")
  p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))

  cap <- 5
  df_mRNA <- msi_exps_list[[cancer]][["mRNA"]]
  df_mRNA$log2RPKM_capped <- df_mRNA$log2RPKM
  df_mRNA$log2RPKM_capped[df_mRNA$log2RPKM > cap] <- cap
  df_mRNA$log2RPKM_capped[df_mRNA$log2RPKM < (-cap)] <- -cap
  df <- merge(df, df_mRNA[, c("Sample", "Gene", "log2RPKM_capped")], by = c("Sample", "Gene"), all.x = T)
  p2 = ggplot(df)
  p2 = p2 + geom_tile(aes(x=Sample, y=Gene, fill=log2RPKM_capped), color=NA)#, linetype="blank") 
  p2 = p2 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p2 = p2 + scale_fill_gradientn(name= "mRNA abundance (log2RPKM)", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p2 = p2 + theme_bw() + theme_nogrid()
  p2 = p2 + theme(axis.text.x = element_blank())
  p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
  p2 = p2 + theme(axis.ticks=element_blank(),legend.position="right")
  p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))

  cap <- 2
  df_protein <- msi_exps_list[[cancer]][["protein"]]
  df_protein$score_capped <- df_protein$score
  df_protein$score_capped[df_protein$score > cap] <- cap
  df_protein$score_capped[df_protein$score < (-cap)] <- -cap
  df <- merge(df, df_protein[, c("Sample", "Gene", "score_capped")], by = c("Sample", "Gene"), all.x = T)
  p3 = ggplot(df)
  p3 = p3 + geom_tile(aes(x=Sample, y=Gene, fill=score_capped), color=NA)#, linetype="blank") 
  p3 = p3 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p3 = p3 + scale_fill_gradientn(name= "protein abundance (log2 ratio normalized)", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p3 = p3 + theme_bw() + theme_nogrid()
  p3 = p3 + theme(axis.text.x = element_blank())
  p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
  p3 = p3 + theme(axis.ticks=element_blank(),legend.position="right")
  p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))

  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cancer, '_MMR_somatic&mRNA&protein_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
  grid.newpage()
  # pdf(fn, height = 14*(len+10)/(top+10), width = 8*(len+10)/(top+10))
  pdf(fn, height = 6, width = 5)
  plot123 <- rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last")
  title <- textGrob(paste0(cancer, " somatic mutations and gene/protein abundance in MMR genes"),gp=gpar(fontsize=16))
  padding <- unit(0,"lines")
  plottgt <- gtable_add_rows(plot123, 
                             heights = grobHeight(title) + padding,
                             pos = 0)
  plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
  grid.draw(plottgt)
  dev.off()
}

for (cancer in c("CO")) {
  df_protein <- msi_exps_list[[cancer]][["protein"]][, c("Gene", "score", "MSI_status")] ; df_protein$expression_type <- "protein"
  tmp <- mmr_somatic[, colnames(mmr_somatic) %in% unique(msi_exps_list[[cancer]][["mRNA"]]$Sample)]
  df_somatic <- data.frame(Gene = mmr_somatic$Hugo_Symbol)
  df_somatic <- cbind(df_somatic, tmp)
  df_somatic_m <- melt(df_somatic, id.vars = c("Gene"))
  colnames(df_somatic_m) <- c("Gene", "Sample", "somatic_tier1_count")
  df_somatic_m <- merge(df_somatic_m, msi_score, by = c("Sample"), all.x = T)
  df_somatic_m$MSI_status <- ifelse(df_somatic_m$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
  # df_somatic_m <- df_somatic_m[df_somatic_m$somatic_tier1_count > 0,]
  df <- df_somatic_m
  cap <- 5
  df_mRNA <- msi_exps_list[[cancer]][["mRNA"]]
  df_mRNA$log2RPKM_capped <- df_mRNA$log2RPKM
  df_mRNA$log2RPKM_capped[df_mRNA$log2RPKM > cap] <- cap
  df_mRNA$log2RPKM_capped[df_mRNA$log2RPKM < (-cap)] <- -cap
  df <- merge(df, df_mRNA[, c("Sample", "Gene", "log2RPKM_capped")], by = c("Sample", "Gene"), all.x = T)
  p2 = ggplot(df)
  p2 = p2 + geom_tile(aes(x=Sample, y=Gene, fill=log2RPKM_capped), color=NA)#, linetype="blank") 
  p2 = p2 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p2 = p2 + geom_text(aes(x=Sample, y=Gene, label = ifelse(somatic_tier1_count > 0 , somatic_tier1_count, NA)), size = 1.5)
  p2 = p2 + scale_fill_gradientn(name= "mRNA abundance (log2RPKM)", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p2 = p2 + theme_bw() + theme_nogrid()
  p2 = p2 + theme(axis.text.x = element_blank())
  p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
  p2 = p2 + theme(axis.ticks=element_blank(),legend.position="right", legend.text = element_text(size = 5))
  p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  
  cap <- 2
  df_protein <- msi_exps_list[[cancer]][["protein"]]
  df_protein$score_capped <- df_protein$score
  df_protein$score_capped[df_protein$score > cap] <- cap
  df_protein$score_capped[df_protein$score < (-cap)] <- -cap
  df <- merge(df, df_protein[, c("Sample", "Gene", "score_capped")], by = c("Sample", "Gene"), all.x = T)
  p3 = ggplot(df)
  p3 = p3 + geom_tile(aes(x=Sample, y=Gene, fill=score_capped), color=NA)#, linetype="blank") 
  p3 = p3 + facet_grid(.~MSI_status, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p3 = p3 + geom_text(aes(x=Sample, y=Gene, label = ifelse(somatic_tier1_count > 0 , somatic_tier1_count, NA)), size = 1.5)
  p3 = p3 + scale_fill_gradientn(name= "protein abundance (log2 ratio normalized)", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p3 = p3 + theme_bw() + theme_nogrid()
  p3 = p3 + theme(axis.text.x = element_blank())
  p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5))
  p3 = p3 + theme(axis.ticks=element_blank(),legend.position="right", legend.text = element_text(size = 5))
  p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cancer, '_MMR_somatic&mRNA&protein_abundance_between_MSI-H_and_rest_in_', cancer, ".v2.pdf")
  grid.newpage()
  # pdf(fn, height = 14*(len+10)/(top+10), width = 8*(len+10)/(top+10))
  pdf(fn, height = 6, width = 8)
  plot123 <- rbind(ggplotGrob(p2), ggplotGrob(p3), size = "last")
  title <- textGrob(paste0(cancer, " somatic mutations and gene/protein abundance in MMR genes"),gp=gpar(fontsize=16))
  padding <- unit(0,"lines")
  plottgt <- gtable_add_rows(plot123, 
                             heights = grobHeight(title) + padding,
                             pos = 0)
  plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
  grid.draw(plottgt)
  dev.off()
}





