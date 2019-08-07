# Yige Wu @ WashU 2018 Nov
# parse complex information from corum, signor and reactome paired gene table

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
## input signor pair table
signor_complex_pair_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_signor_complex/signor_complex_pair_tab.txt", data.table = F)
signor_complex_pair_tab$pair_pro <- paste0(signor_complex_pair_tab$geneA, ":", signor_complex_pair_tab$geneB)
which(duplicated(signor_complex_pair_tab$pair_pro))
signor_complex_pair_uniq <- signor_complex_pair_tab[!duplicated(signor_complex_pair_tab$pair_pro), c("geneA", "geneB", "pair_pro")]
signor_complex_pair_uniq$Source <- "Signor"

## input corum pair table
corum_complex_pair_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_corum_complex/corum_complex_pair_tab.txt", data.table = F)
corum_complex_pair_tab$pair_pro <- paste0(corum_complex_pair_tab$geneA, ":", corum_complex_pair_tab$geneB)
corum_complex_pair_uniq <- corum_complex_pair_tab[!duplicated(corum_complex_pair_tab$pair_pro), c("geneA", "geneB", "pair_pro")]
corum_complex_pair_uniq$Source <- "Corum"

# input corum and reactome complex table from Jessica ---------------------
reactome_corum_tab <- read_csv(file = "./Ding_Lab/Projects_Current/PPI/PPI_shared_data/protein_complex/MERGED_complexesGeneNames_v2.csv")
reactome_tab <- reactome_corum_tab %>%
  filter(Database == "REACTOME")

# reshape the reactome complex table to pair table ------------------------
reactome_pairs <- sapply(X = reactome_tab$geneNameIdentified, FUN = function(x) {
  str_tmp <- x
  genes_tmp <- str_split(string = x, pattern = ";")[[1]]
  pairs_df_tmp <- expand.grid(genes_tmp, genes_tmp)
  pairs_df_tmp <- pairs_df_tmp %>%
    filter(Var1 != Var2)
  pairs_tmp <- paste0(pairs_df_tmp$Var1, ":", pairs_df_tmp$Var2)
  return(pairs_tmp)
})
reactome_pairs %>% head
reactome_pairs %>% tail()

reactome_pairs <- unlist(reactome_pairs)
reactome_pairs <- reactome_pairs[reactome_pairs != ":"]
reactome_pair_tab <- data.frame(str_split_fixed(string = reactome_pairs, pattern = ":", n = 2))
colnames(reactome_pair_tab) <- c("geneA", "geneB")
reactome_pair_tab$pair_pro <- reactome_pairs
reactome_pair_tab$Source <- "Reactome"
reactome_pair_tab %>% nrow()
reactome_pair_tab_uniq <- unique(reactome_pair_tab)
reactome_pair_tab_uniq %>% nrow()

# manual adding SETD2 related proteins ------------------------------------
## source: https://www.uniprot.org/uniprot/Q9BYW2
### FACT complex
SETD2_genes <- c("SSRP1", # https://reactome.org/content/detail/R-HSA-112417
                 "SUPT16H", 
                 "POLR2A", # Plays a role in chromatin structure modulation during elongation by coordinating recruitment of the FACT complex and by interacting with hyperphosphorylated POLR2A (PubMed:23325844)
                 "MSH6", # Acts as a key regulator of DNA mismatch repair in G1 and early S phase by generating H3K36me3, a mark required to recruit MSH6 subunit of the MutS alpha complex: early recruitment of the MutS alpha complex to chromatin to be replicated allows a quick identification of mismatch DNA to initiate the mismatch repair reaction (PubMed:23622243)
                 "RAD51", # Required for DNA double-strand break repair in response to DNA damage: acts by mediating formation of H3K36me3, promoting recruitment of RAD51 and DNA repair via homologous recombination (HR) (PubMed:24843002).
                 "DNMT3A", # H3K36me3 also plays an essential role in the maintenance of a heterochromatic state, by recruiting DNA methyltransferase DNMT3A (PubMed:27317772).
                 "FGFR2", # Required for endoderm development by promoting embryonic stem cell differentiation toward endoderm: acts by mediating formation of H3K36me3 in distal promoter regions of FGFR3, leading to regulate transcription initiation of FGFR3 (By similarity)
                 "STAT1", # also mediates methylation of other proteins, such as tubulins and STAT1 (PubMed:27518565, PubMed:28753426)
                 "TUBA1B", # Trimethylates 'Lys-40' of alpha-tubulins such as TUBA1B (alpha-TubK40me3); alpha-TubK40me3 is required for normal mitosis and cytokinesis and may be a specific tag in cytoskeletal remodeling (PubMed:27518565).
                 "ISG15") # Involved in interferon-alpha-induced antiviral defense by mediating both monomethylation of STAT1 at 'Lys-525' and catalyzing H3K36me3 on promoters of some interferon-stimulated genes (ISGs) to activate gene transcription (PubMed:28753426).
SETD2_pair_tab <- data.frame(geneA = c(rep("SETD2", length(SETD2_genes)), SETD2_genes),
                             geneB = c(SETD2_genes, rep("SETD2", length(SETD2_genes))))

# manual adding TP53 related proteins -------------------------------------
TP53_genes <- c("TP53BP1", "CREBBP", "KAT5", "SIRT1", "BAX", "PUMA", "NOXA", "CDKN1A", "SLC7A11", "GLS2", "TIGAR", "DRAM1", "FOXR2", "LIF", "THBS1", "ICAM1", "SCO2", "SERPINE1")
TP53_pair_tab <- data.frame(geneA = c(rep("TP53", length(TP53_genes)), TP53_genes),
                             geneB = c(TP53_genes, rep("TP53", length(TP53_genes))))


# combine all manually added pairs ----------------------------------------
manual_pair_tab <- rbind(SETD2_pair_tab, TP53_pair_tab) %>%
  unique() %>%
  mutate(pair_pro = paste0(geneA, ":", geneB)) %>%
  mutate(Source = "Manual_Uniprot")

# combine all complex sources  ---------------------------------------------------------------
sup_complex_pair_uniq <- unique(rbind(reactome_pair_tab_uniq, corum_complex_pair_uniq, signor_complex_pair_uniq))
sup_complex_pair_uniq <- unique(rbind(sup_complex_pair_uniq, manual_pair_tab[!(manual_pair_tab$pair_pro %in% sup_complex_pair_uniq$pair_pro),]))
write.table(x = sup_complex_pair_uniq, file = paste0(makeOutDir(resultD = resultD), "reactome_corum_signor_complex_pair_uniq_v1.txt"), quote = F, sep = "\t", row.names = F)
