# Yige Wu @ WashU Mar 2018 
## annotate Omni kinase-substrate pairs from UniprotKB IDs to gene names

# source -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
library(UniProt.ws)


# inputs ------------------------------------------------------------------
ptms <- read_delim("~/Box Sync/pan3can_shared_data/Phospho_databases/OmniPath/ptms.txt",
                   "\t", escape_double = FALSE, trim_ws = TRUE)

## homo sapiens taxids: 9606
up <- UniProt.ws(taxId=9606)
res <- select(up,keys = c("P17252"),columns = c("GENES", "HGNC"),keytype = "UNIPROTKB")
res <- select(up,keys = c("PRKACA"),columns = c("UNIPROTKB", "HGNC"),keytype = "GENES")
