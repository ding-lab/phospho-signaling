# Yige Wu @ WashU 2018 Jan
# summarize the prevalence of up-regulated and down-regulated ks pairs

# inputs ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()

datatypes <- c("PRO", "PHO", "collapsed_PHO")


# get summary numbers -----------------------------------------------------
for (sig in c(0.1)) {
  upstream <- "kinase"
  ks_table <- load_ks_table(upstream)
  sum_table <- NULL

  Pho.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "PHO_diffexp_3can_paired_FDR", sig, ".txt"),
                           data.table = F)
  Pho.diffexp3can <- cbind(Pho.diffexp3can, formatPhosphosite(phosphosite_vector = Pho.diffexp3can$Phosphosite, Pho.diffexp3can$Gene))
  Pho.diffexp3can <- Pho.diffexp3can[, !(colnames(Pho.diffexp3can) %in% c("Gene", "Phosphosite"))]
  colnames(Pho.diffexp3can)[colnames(Pho.diffexp3can) == "Fold Change"] <- "FC_sub"
  
  Phog.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "collapsed_PHO_diffexp_3can_paired_FDR", sig, ".txt"),
                            data.table = F)
  Phog.diffexp3can$KINASE <- Phog.diffexp3can$Gene
  Phog.diffexp3can <- Phog.diffexp3can[, !(colnames(Phog.diffexp3can) %in% c("Gene", "Phosphosite"))]
  colnames(Phog.diffexp3can)[colnames(Phog.diffexp3can) == "Fold Change"] <- "FC_kin"
  
  Pro.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "PRO_diffexp_3can_paired_FDR", sig, ".txt"),
                           data.table = F)
  Pro.diffexp3can$KINASE <- Pro.diffexp3can$Gene
  Pro.diffexp3can <- Pro.diffexp3can[, !(colnames(Pro.diffexp3can) %in% c("Gene", "Phosphosite"))]
  colnames(Pro.diffexp3can)[colnames(Pro.diffexp3can) == "Fold Change"] <- "FC_kin"
  
  for (cancer in c("BRCA", "OV", "CO")) {
    ## summarize protein/phosphoprotein number input to differential expression analysis (all, kinase, substrate)
    phog <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "collapsed_PHO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                  data.table = F)
    tmp <- data.frame(datatype = "collapsed_PHO",
                      cancer = cancer,
                      diffexp_type = "all",
                      gene_cat = c("kinase", "not-kinase"),
                      num_gene = as.matrix(table(phog$Gene %in% ks_table$GENE))[c("TRUE", "FALSE"),1])
    tmp$num_site <- tmp$num_gene
    sum_table <- rbind(sum_table, tmp)
    
    pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                 data.table = F)
    tmp <- data.frame(datatype = "PHO",
                      cancer = cancer,
                      diffexp_type = "all",
                      gene_cat = c("substrate", "not-substrate"),
                      num_gene = as.matrix(table(unique(pho$Gene) %in% ks_table$SUB_GENE))[c("TRUE", "FALSE"),1],
                      num_site = as.matrix(table(pho$Gene %in% ks_table$SUB_GENE))[c("TRUE", "FALSE"),1])
    sum_table <- rbind(sum_table, tmp)
    
    pro <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                 data.table = F)
    tmp <- data.frame(datatype = "PRO",
                      cancer = cancer,
                      diffexp_type = "all",
                      gene_cat = c("kinase", "not-kinase"),
                      num_gene = as.matrix(table(pro$Gene %in% ks_table$GENE))[c("TRUE", "FALSE"),1])
    tmp$num_site <- tmp$num_gene
    sum_table <- rbind(sum_table, tmp)
    
    ## summarize protein/phosphoprotein number output by differential expression analysis
    phog_can <- Phog.diffexp3can[Phog.diffexp3can$Cancer == cancer,]
    pho_can <- Pho.diffexp3can[Pho.diffexp3can$Cancer == cancer,]
    pro_can <- Pro.diffexp3can[Pro.diffexp3can$Cancer == cancer,]
    
    for (diffexp_type in c("up", "down")) {
      phog_diff <- phog_can[phog_can$diffexp_type == diffexp_type,]
      tmp <- data.frame(datatype = "collapsed_PHO",
                        cancer = cancer,
                        diffexp_type = diffexp_type,
                        gene_cat = c("kinase", "not-kinase"),
                        num_gene = as.matrix(table(phog_diff$KINASE %in% ks_table$GENE))[c("TRUE", "FALSE"),1])
      tmp$num_site <- tmp$num_gene
      sum_table <- rbind(sum_table, tmp)
      
      pho_diff <- pho_can[pho_can$diffexp_type == diffexp_type,]
      tmp <- data.frame(datatype = "PHO",
                        cancer = cancer,
                        diffexp_type = diffexp_type,
                        gene_cat = c("substrate", "not-substrate"),
                        num_gene = as.matrix(table(unique(pho_diff$SUBSTRATE) %in% ks_table$SUB_GENE))[c("TRUE", "FALSE"),1],
                        num_site = as.matrix(table(pho_diff$SUBSTRATE %in% ks_table$SUB_GENE))[c("TRUE", "FALSE"),1])
      sum_table <- rbind(sum_table, tmp)
      
      pro_diff <- pro_can[pro_can$diffexp_type == diffexp_type,]
      tmp <- data.frame(datatype = "PRO",
                        cancer = cancer,
                        diffexp_type = diffexp_type,
                        gene_cat = c("kinase", "not-kinase"),
                        num_gene = as.matrix(table(pro_diff$KINASE %in% ks_table$GENE))[c("TRUE", "FALSE"),1])
      tmp$num_site <- tmp$num_gene
      sum_table <- rbind(sum_table, tmp)
      
      ## summarize kinase-substrate pairs that are consistently high/low
      ks_can <- ks_table[ks_table$GENE %in% phog$Gene & ks_table$SUB_GENE %in% pho$Gene,]
      tmp <- data.frame(datatype = "trans",
                        cancer = cancer,
                        diffexp_type = diffexp_type,
                        gene_cat = c("consistent", "not-consistent"),
                        num_gene = as.matrix(table(ks_can$GENE %in% phog_diff$KINASE & ks_can$SUB_GENE %in% pho_diff$SUBSTRATE))[c("TRUE", "FALSE"),1])
      tmp$num_site <- NA
      sum_table <- rbind(sum_table, tmp)
    }
  }
}


# get actual genes in summary numbers-----------------------------------------------------
for (sig in c(0.1)) {
  upstream <- "kinase"
  ks_table <- load_ks_table(upstream)
  sum_table <- NULL
  
  Pho.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "PHO_diffexp_3can_paired_FDR", sig, ".txt"),
                           data.table = F)
  Pho.diffexp3can <- cbind(Pho.diffexp3can, formatPhosphosite(phosphosite_vector = Pho.diffexp3can$Phosphosite, Pho.diffexp3can$Gene))
  Pho.diffexp3can <- Pho.diffexp3can[, !(colnames(Pho.diffexp3can) %in% c("Gene", "Phosphosite"))]
  colnames(Pho.diffexp3can)[colnames(Pho.diffexp3can) == "Fold Change"] <- "FC_sub"
  
  Phog.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "collapsed_PHO_diffexp_3can_paired_FDR", sig, ".txt"),
                            data.table = F)
  Phog.diffexp3can$KINASE <- Phog.diffexp3can$Gene
  Phog.diffexp3can <- Phog.diffexp3can[, !(colnames(Phog.diffexp3can) %in% c("Gene", "Phosphosite"))]
  colnames(Phog.diffexp3can)[colnames(Phog.diffexp3can) == "Fold Change"] <- "FC_kin"
  
  Pro.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", "PRO_diffexp_3can_paired_FDR", sig, ".txt"),
                           data.table = F)
  Pro.diffexp3can$KINASE <- Pro.diffexp3can$Gene
  Pro.diffexp3can <- Pro.diffexp3can[, !(colnames(Pro.diffexp3can) %in% c("Gene", "Phosphosite"))]
  colnames(Pro.diffexp3can)[colnames(Pro.diffexp3can) == "Fold Change"] <- "FC_kin"
  
  for (cancer in c("BRCA", "OV", "CO")) {
    sink(paste0(makeOutDir(), cancer,"_", upstream, '_substrate_pairs_summary_log.txt'), append=FALSE, split=FALSE)
    cat(paste0("-", cancer, " summary:\n"))
    
    ## summarize protein/phosphoprotein number input to differential expression analysis (all, kinase, substrate)
    phog <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "collapsed_PHO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                  data.table = F)
    tmp <- data.frame(datatype = "collapsed_PHO",
                      cancer = cancer,
                      diffexp_type = "all",
                      gene_cat = c("kinase", "not-kinase"),
                      num_gene = as.matrix(table(phog$Gene %in% ks_table$GENE))[c("TRUE", "FALSE"),1])
    tmp$num_site <- tmp$num_gene
    sum_table <- rbind(sum_table, tmp)
    
    pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                 data.table = F)
    tmp <- data.frame(datatype = "PHO",
                      cancer = cancer,
                      diffexp_type = "all",
                      gene_cat = c("substrate", "not-substrate"),
                      num_gene = as.matrix(table(unique(pho$Gene) %in% ks_table$SUB_GENE))[c("TRUE", "FALSE"),1],
                      num_site = as.matrix(table(pho$Gene %in% ks_table$SUB_GENE))[c("TRUE", "FALSE"),1])
    sum_table <- rbind(sum_table, tmp)
    
    pro <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_replicate_averaged_imputed.txt"),
                 data.table = F)
    tmp <- data.frame(datatype = "PRO",
                      cancer = cancer,
                      diffexp_type = "all",
                      gene_cat = c("kinase", "not-kinase"),
                      num_gene = as.matrix(table(pro$Gene %in% ks_table$GENE))[c("TRUE", "FALSE"),1])
    tmp$num_site <- tmp$num_gene
    sum_table <- rbind(sum_table, tmp)
    
    ## summarize protein/phosphoprotein number output by differential expression analysis
    phog_can <- Phog.diffexp3can[Phog.diffexp3can$Cancer == cancer,]
    pho_can <- Pho.diffexp3can[Pho.diffexp3can$Cancer == cancer,]
    pro_can <- Pro.diffexp3can[Pro.diffexp3can$Cancer == cancer,]
    
    for (diffexp_type in c("up", "down")) {
      phog_diff <- phog_can[phog_can$diffexp_type == diffexp_type,]
      tmp <- data.frame(datatype = "collapsed_PHO",
                        cancer = cancer,
                        diffexp_type = diffexp_type,
                        gene_cat = c("kinase", "not-kinase"),
                        num_gene = as.matrix(table(phog_diff$KINASE %in% ks_table$GENE))[c("TRUE", "FALSE"),1])
      tmp$num_site <- tmp$num_gene
      sum_table <- rbind(sum_table, tmp)
      ## record diffexp genes
      cat(paste0("--", diffexp_type, " regulated collapsed_PHO kinases:\n"))
      cat(phog_diff$KINASE[phog_diff$KINASE %in% ks_table$GENE])
      cat(paste0("\n"))
      
      pho_diff <- pho_can[pho_can$diffexp_type == diffexp_type,]
      tmp <- data.frame(datatype = "PHO",
                        cancer = cancer,
                        diffexp_type = diffexp_type,
                        gene_cat = c("substrate", "not-substrate"),
                        num_gene = as.matrix(table(unique(pho_diff$SUBSTRATE) %in% ks_table$SUB_GENE))[c("TRUE", "FALSE"),1],
                        num_site = as.matrix(table(pho_diff$SUBSTRATE %in% ks_table$SUB_GENE))[c("TRUE", "FALSE"),1])
      sum_table <- rbind(sum_table, tmp)
      for (gene in phog_diff$KINASE[phog_diff$KINASE %in% ks_table$GENE]) {
        subs_all_tmp <- unique(as.vector(ks_table$SUB_GENE)[ks_table$GENE==gene])
        subs_diff_tmp <- subs_all_tmp[subs_all_tmp %in% as.vector(pho_diff$SUBSTRATE)]
        cat(paste0("---", gene, "(", length(subs_diff_tmp), "substrates", diffexp_type, "):"))
        cat(subs_diff_tmp)
        cat(paste0("\n"))
      }
      
      pro_diff <- pro_can[pro_can$diffexp_type == diffexp_type,]
      tmp <- data.frame(datatype = "PRO",
                        cancer = cancer,
                        diffexp_type = diffexp_type,
                        gene_cat = c("kinase", "not-kinase"),
                        num_gene = as.matrix(table(pro_diff$KINASE %in% ks_table$GENE))[c("TRUE", "FALSE"),1])
      tmp$num_site <- tmp$num_gene
      sum_table <- rbind(sum_table, tmp)

      ## summarize kinase-substrate pairs that are consistently high/low
      ks_can <- ks_table[ks_table$GENE %in% phog$Gene & ks_table$SUB_GENE %in% pho$Gene,]
      tmp <- data.frame(datatype = "trans",
                        cancer = cancer,
                        diffexp_type = diffexp_type,
                        gene_cat = c("consistent", "not-consistent"),
                        num_gene = as.matrix(table(ks_can$GENE %in% phog_diff$KINASE & ks_can$SUB_GENE %in% pho_diff$SUBSTRATE))[c("TRUE", "FALSE"),1])
      tmp$num_site <- NA
      sum_table <- rbind(sum_table, tmp)
    }
    sink()
    closeAllConnections()
  }

}


