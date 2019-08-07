# Yige Wu @ WashU Apr 2019
## generate a gene*sample CNV table from individual gene-annotated segmentation files for CPTAC3 ccRCC discovery cohort

# source ------------------------------------------------------------------
source("~/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R")

# set variables -----------------------------------------------------------

# input CPTAC3 case dat file ----------------------------------------------
cptac3_case_dat <- fread("../CPTAC3.catalog/CPTAC3.cases.dat", data.table = F, sep = " ")
tmp <- str_split_fixed(n = 4, pattern = "\t", string = cptac3_case_dat$`# Case	Disease	Cohort	AnalysisBatch`)
## some cases within Y2.b2 have weird separator, ignore for now
tmp[tmp[,3] == "",]
cptac3_case_dat$Case <- tmp[,1]
cptac3_case_dat$Disease <- tmp[,2]
cptac3_case_dat$Cohort <- tmp[,3]
cptac3_case_dat$AnalysisBatch <- tmp[,4]
cptac3_case_dat <- as.data.frame(cptac3_case_dat)
rows2correct <- (cptac3_case_dat$Cohort == "" & cptac3_case_dat$AnalysisBatch != "")
tmp <- cptac3_case_dat$AnalysisBatch[rows2correct]
tmp_corrected <- str_split_fixed(string = tmp, pattern = "\t", n = 2)
tmp_corrected
cptac3_case_dat[rows2correct, "Cohort"] <- tmp_corrected[,1]
cptac3_case_dat[rows2correct, "AnalysisBatch"] <- tmp_corrected[,2]
cptac3_case_dat$`# Case	Disease	Cohort	AnalysisBatch` <- NULL
colnames(cptac3_case_dat) <- c("# Case", "Disease", "Cohort", "AnalysisBatch")

## write table
write.table(x = cptac3_case_dat, file = paste0(makeOutDir(), "CPTAC3.cases.dat"), sep = "\t", row.names = F, quote = F)
