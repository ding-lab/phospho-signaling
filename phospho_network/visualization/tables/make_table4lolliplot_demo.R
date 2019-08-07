# Yige Wu @ WashU 2019 July
## make table for inputting into Protein Paint

# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)

# set variables -----------------------------------------------------------
cancers2process <- c("UCEC")
genes2process <- c("KEAP1", "NFE2L2")

# loop around cancer types ------------------------------------------------
for (gene_tmp in genes2process) {
  for (cancer_type_tmp in cancers2process) {
    maf_tmp <- loadMaf(cancer = cancer_type_tmp)
    
    tab4jude <- maf_tmp %>%
      filter(Hugo_Symbol == gene_tmp) %>%
      mutate(text4proteinpaint = HGVSp_Short) %>%
      mutate(Position = str_split_fixed(string = HGVSp_Short, pattern = '[A-Z]|[a-z]|\\*', n = 4)[, 3]) %>%
      mutate(judeClass = substr(x = Variant_Classification, start = 1, stop = 1)) %>%
      select(text4proteinpaint, Position, judeClass)

    write.table(x = tab4jude, file = paste0(makeOutDir(), "tab4jude_", gene_tmp, "_", cancer_type_tmp, ".txt"), 
                row.names = F, quote = F, sep = ";", col.names = F)
    
  }
  
} 

