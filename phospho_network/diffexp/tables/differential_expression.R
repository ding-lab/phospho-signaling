# Yige Wu @ WashU 2018 Feb
# differential expression of protein and phosphoprotein (missing values imputed) using SAM
## TODO: get tumor normal pair info

# inputs ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')

resultDnow <- makeOutDir()
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180127.txt"), sep = "\t")
library(samr)


# function ----------------------------------------------------------------
getSAMfiltered <- function(expression, identifier) {
  ## construct vector differentiating tumor and normal
  expression_t <- get_tumor(expression, clinical.m = clinical)
  expression_n <- get_normal(expression, clinical.m = clinical)
  y = c(rep(2, ncol(expression_t)), rep(1, ncol(expression_n)))
  x = cbind(expression_t, expression_n)
  
  ## keep only the tumors and normals
  samfit<-SAM(x,y,resp.type="Two class unpaired", nperms = 100,
              geneid = expression[, identifier], genenames = expression[, "Gene"],
              fdr.output = 0.05,
              logged2 = T, testStatistic = "wilcoxon")
  
  ## filter by fold change
  siggenes.table <- samfit$siggenes.table
  genes.up.filtered <- as.data.frame(siggenes.table$genes.up)
  genes.up.filtered <- genes.up.filtered[as.numeric(as.vector(genes.up.filtered$`Fold Change`)) > 2,]
  genes.up.filtered$diffexp_type <- "up"
  genes.lo.filtered <- as.data.frame(siggenes.table$genes.lo)
  genes.lo.filtered <- genes.lo.filtered[as.numeric(as.vector(genes.lo.filtered$`Fold Change`)) < 1/2,]
  genes.lo.filtered$diffexp_type <- "down"
  genes.diffexp <- rbind(genes.up.filtered, genes.lo.filtered)
  return(genes.diffexp)
}

# paired tumor vs normal ---------------------------------------------------------------


# unpaired tumor vs normal ---------------------------------------------------
Pro.diffexp3can <- NULL
Pho.diffexp3can <- NULL
for (cancer in c("BRCA", "OV", "CO")) {
  ## input protein and phosphoprotein data
  Pro <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_replicate_averaged_imputed.txt"),
               data.table = F)
  Pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_replicate_averaged_imputed.txt"),
               data.table = F)

  ## combine tables
  Pro.diffexp <- getSAMfiltered(Pro, 'Gene')
  Pro.diffexp$Cancer <- cancer
  Pro.diffexp3can <- rbindlist(list(Pro.diffexp3can, Pro.diffexp))
  Pho.diffexp <- getSAMfiltered(Pho, 'Phosphosite')
  Pho.diffexp$Cancer <- cancer
  Pho.diffexp3can <- rbindlist(list(Pho.diffexp3can, Pho.diffexp))
}
colnames(Pro.diffexp3can)[1:2] <- c("Gene", "Phosphosite"); Pro.diffexp3can$Phosphosite <- "Protein"
colnames(Pho.diffexp3can)[1:2] <- c("Gene", "Phosphosite")
write.table(Pro.diffexp3can, row.names = F, quote=F, sep = '\t', file=paste0(resultDnow, "/", "Pro.diffexp3can.txt"))
write.table(Pho.diffexp3can, row.names = F, quote=F, sep = '\t', file=paste0(resultDnow, "/", "Pho.diffexp3can.txt"))




