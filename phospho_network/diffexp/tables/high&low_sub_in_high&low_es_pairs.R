# Yige Wu @ WashU 2018 Jul
# output the up-regulated substrates of up-regulated kinases in all 3 cancers, whether there's difference in the downstream signaling

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(dplyr)
# source("https://bioconductor.org/biocLite.R")
# biocLite("RDAVIDWebService")
#sudo ln -f -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib

library("RDAVIDWebService")



# variables ---------------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
reg_sig <- 0.05
diff_sig <- 0.2
diff_log2fc <- 1

# inputs -------------------------------------------------------------------
sup_cans_tab <- NULL
for (cancer in cancers_sort) {
  sup_tab <- fread(paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                          cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  sup_tab$Cancer <- cancer
  sup_cans_tab <- rbind(sup_cans_tab, sup_tab)
}

# classify ----------------------------------------------------------------
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == "kinase",]
sup_cans_tab_pk <- markSigSiteCan(sup_cans_tab_pk, sig_thres = reg_sig, enzyme_type = "kinase")
sup_cans_tab_pp <- sup_cans_tab[sup_cans_tab$enzyme_type == "phosphotase",]
sup_cans_tab_pp <- markSigSiteCan(sup_cans_tab_pp, sig_thres = reg_sig, enzyme_type = "phosphotase")
sup_cans_tab <- rbind(sup_cans_tab_pk, sup_cans_tab_pp)
## annotate kinase substrate regulation
sup_cans_tab$regulated <- (sup_cans_tab$coef_sig & sup_cans_tab$fdr_sig)
## annotate up and down regulation
sup_cans_tab$diffexp_type <- "other"
diffexp_sig <- ((!is.na(sup_cans_tab$KSEA_pvalue) & sup_cans_tab$KSEA_pvalue < diff_sig) | (!is.na(sup_cans_tab$diffexp_log2FC & abs(sup_cans_tab$diffexp_log2FC) > diff_log2fc)))
sup_cans_tab$diffexp_type[diffexp_sig & sup_cans_tab$enzyme_direction == "up" & sup_cans_tab$substrate_direction == "up"] <- "up"
sup_cans_tab$diffexp_type[diffexp_sig & sup_cans_tab$enzyme_direction == "down" & sup_cans_tab$substrate_direction == "down"] <- "down"
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == "kinase",]
sup_cans_tab_pp <- sup_cans_tab[sup_cans_tab$enzyme_type == "phosphotase",]

## pairs with any regulation
length(unique(sup_cans_tab_pk$pair[sup_cans_tab_pk$regulated | sup_cans_tab_pk$diffexp_type != "other"]))
length(unique(sup_cans_tab_pk$pair[sup_cans_tab_pk$regulated]))
length(unique(sup_cans_tab_pk$pair[sup_cans_tab_pk$diffexp_type != "other"]))
length(unique(sup_cans_tab_pp$pair[sup_cans_tab_pp$regulated | sup_cans_tab_pp$diffexp_type != "other"]))


# get uniprot IDs ---------------------------------------------------------
## summarize the Uniprot IDs for the substrates
tmp <- sup_cans_tab[!is.na(sup_cans_tab$SUB_ACC_ID) & sup_cans_tab$SUB_ACC_ID != "NULL",]
sub_acc_ids <- unique(tmp[, c("SUB_GENE", "SUB_ACC_ID")])

tmp <- sup_cans_tab[!is.na(sup_cans_tab$KIN_ACC_ID) & sup_cans_tab$KIN_ACC_ID != "NULL",]
kin_acc_ids <- unique(tmp[, c("GENE", "KIN_ACC_ID")])

## input gene names with uniprot IDs
op <- fread(input = "cptac2p/analysis_results/preprocess_files/tables/convert_omnipath2genomic_coord/op_w.phosphosite.availability.txt", data.table = F)
op <- data.frame(op)
op_sub_ids <- unique(op[, c("SUB_GENE", "substrate")])
colnames(op_sub_ids) <- c("SUB_GENE", "SUB_ACC_ID")

op_kin_ids <- unique(op[, c("GENE", "enzyme")])
colnames(op_kin_ids) <- c("GENE", "KIN_ACC_ID")


## refill the super table with substrate uniprot IDs
sub_acc_ids <- unique(rbind(sub_acc_ids, op_sub_ids))
sup_cans_tab$SUB_ACC_ID <- NULL
sup_cans_tab <- merge(sup_cans_tab, sub_acc_ids, all.x = T)

kin_acc_ids <- unique(rbind(kin_acc_ids, op_kin_ids))
sup_cans_tab$KIN_ACC_ID <- NULL
sup_cans_tab <- merge(sup_cans_tab, kin_acc_ids, all.x = T)

write.table(x = kin_acc_ids, file = paste0("./cptac2p/analysis_results/enzyme_acc_ids.txt"), quote = F, row.names = F, col.names = T)

# summarize kinase-substrate pairs ----------------------------------------
sum_tab3 <- data.frame(table(unique(sup_cans_tab_pk[, c("pair", "GENE", "diffexp_type", "Cancer")])[, c("GENE","diffexp_type", "Cancer")]))
tab_pk_cans_trans <- sup_cans_tab_pk[sup_cans_tab_pk$SELF == "trans" & !is.na(sup_cans_tab_pk$SELF),]

# config DAVID Web Service------------------------------------------------------------
david<-DAVIDWebService(email="yigewu@wustl.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
david$is.connected()

## get all kinases uniprots
for (diffexp_type in c("up", "down")) {
  dir1 <- paste0(makeOutDir(resultD = resultD), diffexp_type, "/")
  dir.create(dir1)
  
  kins_diffexp <- as.vector(unique(sum_tab3$GENE[sum_tab3$diffexp_type == diffexp_type & sum_tab3$Freq >= 5]))
  kins_diffexp_uniprot <- kin_acc_ids$KIN_ACC_ID[kin_acc_ids$GENE %in% kins_diffexp]
  tmp <- addList(david, kins_diffexp_uniprot, idType="UNIPROT_ACCESSION", listName= paste0(diffexp_type, "_", enzyme_type), listType="Gene")
  fn <- paste0(dir1, diffexp_type, "_",enzyme_type, "_DAVID.txt")
  getFunctionalAnnotationTableFile(david, fn)
}


for (diffexp_type in c("up", "down")) {
  ## get the list of kinase-substrate pairs shared across all cancers
  df <- sum_tab3[sum_tab3$diffexp_type == diffexp_type & sum_tab3$Freq > 1,]
  df_count2 <- data.frame(table(df$GENE))
  df <- df[df$GENE %in% df_count2$Var1[df_count2$Freq == length(cancers_sort)],]
  tab_pk_cans_diffexp <- tab_pk_cans_trans[tab_pk_cans_trans$diffexp_type == diffexp_type,]
  tab_pk_cans_other <- tab_pk_cans_trans[tab_pk_cans_trans$diffexp_type == "other",]
  
  dir1 <- paste0(makeOutDir(resultD = resultD), diffexp_type, "/")
  dir.create(dir1)
  

  
  for (enzyme in unique(df$GENE)) {
    tab_pk_cans_diffexp_e <- tab_pk_cans_diffexp[tab_pk_cans_diffexp$GENE == enzyme,]
    tab_pk_cans_other_e <- tab_pk_cans_other[tab_pk_cans_other$GENE == enzyme & !is.na(tab_pk_cans_other$substrate_log2FC),]
    
    dir2 <- paste0(dir1, enzyme, "/")
    dir.create(dir2)
    
    for (cancer in cancers_sort){
      dir3 <- paste0(dir2, cancer, "/")
      dir.create(dir3)
      tab_pk_can_diffexp_e <- tab_pk_cans_diffexp_e[tab_pk_cans_diffexp_e$Cancer == cancer,]
      tab_pk_can_other_e <- tab_pk_cans_other_e[tab_pk_cans_other_e$Cancer == cancer,]
      
      ## output regulated phosphosites
      tmp <- tab_pk_can_diffexp_e[is.na(tab_pk_can_diffexp_e$SUB_ACC_ID),]
      regulated_sites <- unique(tab_pk_can_diffexp_e[!is.na(tab_pk_can_diffexp_e$SUB_ACC_ID), c("SUB_ACC_ID", "SUB_MOD_RSD")])
      write.table(x = regulated_sites, file = paste0(dir3,
                                                     diffexp_type, "_", enzyme_type, "_substrate_sites_in", cancer, ".txt"), 
                  quote = F, sep = '\t', col.names = F, row.names = F)

      regulated_subs <- unique(tab_pk_can_diffexp_e[!is.na(tab_pk_can_diffexp_e$SUB_ACC_ID), c("SUB_ACC_ID")])
      write.table(x = regulated_subs, file = paste0(dir3, diffexp_type, "_", enzyme_type, "_substrate_proteins_in", cancer, ".txt"), quote = F, sep = '\t', col.names = F, row.names = F)
      
      ## write out functional annotation for regulated proteins
      FG <- addList(david, regulated_subs, idType="UNIPROT_ACCESSION", listName= paste0(diffexp_type, "_", enzyme, "_", cancer), listType="Gene")
      fn <- paste0(dir3, "FwDAVID.txt")
      getFunctionalAnnotationTableFile(david, fn)

      ## output regulated phosphosites
      unregulated_sites <- unique(tab_pk_can_other_e[!is.na(tab_pk_can_other_e$SUB_ACC_ID), c("SUB_ACC_ID", "SUB_MOD_RSD")])
      write.table(x = unregulated_sites, file = paste0(dir3, "unregulated", "_", enzyme_type, "_substrate_sites_in", cancer, ".txt"), quote = F, sep = '\t', col.names = F, row.names = F)
      
      unregulated_subs <- unique(tab_pk_can_other_e[!is.na(tab_pk_can_other_e$SUB_ACC_ID), c("SUB_ACC_ID")])
      write.table(x = unregulated_subs, file = paste0(dir3, "unregulated", "_", enzyme_type, "_substrate_proteins_in", cancer, ".txt"), quote = F, sep = '\t', col.names = F, row.names = F)
      
      ## write out functional annotation for regulated proteins
      BG <- addList(david, unregulated_subs, idType="UNIPROT_ACCESSION", listName= paste0("unregulated","_", enzyme, "_", cancer), listType="Gene")
      fn <- paste0(dir3, "BgDAVID.txt")
      getFunctionalAnnotationTableFile(david, fn)
      
    }
  }
}
