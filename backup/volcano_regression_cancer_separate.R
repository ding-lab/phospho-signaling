# Yige Wu @ WashU 2018 Jan
# draw volcano plots for regression result

# choose kinase/phosphotase, cancer , significance level, model -----------------------------------------------
protein <- "kinase"
# protein <- "phosphotase"
sig <- 0.05
out_thres <- 1.5

# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes; it should be one directory out of the working directory
tn = paste0(resultD,"regression/generate/",protein,"_substrate_regression_cptac2p_3can.txt")
df <- read.delim(file = tn)


# process input -----------------------------------------------------------
df <- markSigCan(df, sig_thres = sig)
df$sig_3can <- "other"
df$sig_3can[rowSums(df[,c('sig_BRCA', 'sig_OV', 'sig_CO')]) == 3]  <- "shared_by_3_cancers"
df$sig_3can[rowSums(df[,c('sig_BRCA', 'sig_OV', 'sig_CO')]) == 1]  <- "cancer_type_unique"
color_map <- c("firebrick1", "dodgerblue", "black")
names(color_map) <- c("shared_by_3_cancers", "cancer_type_unique", "other")

# set parameters for plotting
if ( protein == "kinase") {
  plot_fdr_scale <- 5
}
if ( protein == "phosphotase") {
  plot_fdr_scale <- 2
}


# cis volcano plotting module -------------------------------------------------
table_cis_outlier_removed_m = df[df$SELF=="cis",]
table_cis_outlier_removed_m$coef_pro_kin_filtered = remove_outliers(table_cis_outlier_removed_m$coef_pro_kin, out_thres = out_thres)
table_cis_outlier_removed_m = table_cis_outlier_removed_m[!is.na(table_cis_outlier_removed_m$coef_pro_kin_filtered),]
p = ggplot(table_cis_outlier_removed_m,aes(x=coef_pro_kin, y=-log10(FDR_pro_kin)))
p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.05 ,  stroke = 0 , color = "grey")
p = p + geom_text(aes(label= ifelse(-log10(FDR_pro_kin)>plot_fdr_scale, as.character(pair), NA), color = sig_3can),size=1.5,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + geom_vline(xintercept = 0, color = 'grey')
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for kinase protein expression", y="-log10(FDR)")
p = p + scale_color_manual(values = color_map)
p
fn = paste(resultD,'regression/plot/Cis_',protein,'_substrate_volcano_cptac2p_3can.pdf',sep ="")
ggsave(file=fn, height=6, width=15, useDingbats=FALSE)

# trans volcano plotting module -------------------------------------------------
table_trans_outlier_removed_m = df[df$SELF=="trans",]
table_trans_outlier_removed_m$coef_pho_kin_filtered = remove_outliers(table_trans_outlier_removed_m$coef_pho_kin, out_thres = out_thres)
table_trans_outlier_removed_m = table_trans_outlier_removed_m[!is.na(table_trans_outlier_removed_m$coef_pho_kin_filtered),]
p = ggplot(table_trans_outlier_removed_m,aes(x=coef_pho_kin, y=-log10(FDR_pho_kin)))
p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.05 ,  stroke = 0 )
p = p + geom_text(aes(label= ifelse(-log10(FDR_pho_kin)>plot_fdr_scale, as.character(pair), NA), color = sig_3can),size=1.5,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + geom_vline(xintercept = 0, color = 'grey')
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p = p + labs(x = "Coefficient for kinase phosphorylation level", y="-log10(FDR)")
p = p + scale_color_manual(values = color_map)
p
fn = paste(resultD,'regression/plot/Trans_',protein,'_substrate_volcano_cptac2p_3can.pdf',sep ="")
ggsave(file=fn, height=6, width=15, useDingbats=FALSE)

