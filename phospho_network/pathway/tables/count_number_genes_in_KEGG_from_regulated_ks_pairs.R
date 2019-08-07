# Yige Wu @ WashU Mar 2019
## calculate the number of proteins invovled in regulated kinase-substrate pairs

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")

# set variables -----------------------------------------------------------
omnipath_tab <- annotate_ks_source(omnipath_tab)
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")

# input KEGG --------------------------------------------------------------
load(paste0(pan3can_shared_dataD, "Gene_family/2015-08-01_Gene_Set.RData"))

# input regression result table --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
regression %>% nrow()
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[omnipath_tab$is.direct])
regression %>% nrow()
regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)

# KEGG  ---------------------------------------------------------------
all_kegg_names <- names(KEGG)

kegg_regulated_summary <- NULL
for (cancer_tmp in unique(regression$Cancer)) {
  num_pairs_tested <- vector(mode = "numeric", length = length(all_kegg_names))
  num_pairs_regulated <- num_pairs_tested
  count <- 0
  for (kegg_name_tmp in all_kegg_names) {
    regression_tested_tmp <- regression %>%
      filter(Cancer == cancer_tmp) %>%
      filter(SELF == "trans") %>%
      filter(GENE %in% KEGG[[kegg_name_tmp]]) %>%
      filter(SUB_GENE %in% KEGG[[kegg_name_tmp]])
    
    regression_regulated_tmp <- regression %>%
      filter(regulated == T) %>%
      filter(Cancer == cancer_tmp) %>%
      filter(SELF == "trans") %>%
      filter(GENE %in% KEGG[[kegg_name_tmp]]) %>%
      filter(SUB_GENE %in% KEGG[[kegg_name_tmp]])
    
    count <- count + 1
    num_pairs_tested[count] <- regression_tested_tmp %>%
      nrow()
    num_pairs_regulated[count] <- regression_regulated_tmp %>%
      nrow()
  }
  
  summary_tab_tmp <- data.frame(kegg_name = all_kegg_names,
                                num_pairs_tested = num_pairs_tested,
                                num_pairs_regulated = num_pairs_regulated,
                                cancer = cancer_tmp)
  kegg_regulated_summary <- rbind(kegg_regulated_summary, summary_tab_tmp)
}

kegg_regulated_summary$regulated_ratio <- kegg_regulated_summary$num_pairs_regulated/kegg_regulated_summary$num_pairs_tested

kegg_regulated_summary2rank <- kegg_regulated_summary %>%
  filter(num_pairs_tested > 10) %>%
  filter(num_pairs_regulated > 5)


# react both kinase-substrate in pathway  ---------------------------------------------------------------
all_react_names <- names(REACT)
react_regulated_summary <- NULL
for (cancer_tmp in unique(regression$Cancer)) {
  num_pairs_tested <- vector(mode = "numeric", length = length(all_react_names))
  num_pairs_regulated <- num_pairs_tested
  count <- 0
  for (react_name_tmp in all_react_names) {
    regression_tested_tmp <- regression %>%
      filter(Cancer == cancer_tmp) %>%
      filter(SELF == "trans") %>%
      filter(GENE %in% REACT[[react_name_tmp]]) %>%
      filter(SUB_GENE %in% REACT[[react_name_tmp]])
    
    regression_regulated_tmp <- regression %>%
      filter(regulated == T) %>%
      filter(Cancer == cancer_tmp) %>%
      filter(SELF == "trans") %>%
      filter(GENE %in% REACT[[react_name_tmp]]) %>%
      filter(SUB_GENE %in% REACT[[react_name_tmp]])
    
    count <- count + 1
    num_pairs_tested[count] <- regression_tested_tmp %>%
      nrow()
    num_pairs_regulated[count] <- regression_regulated_tmp %>%
      nrow()
  }
  
  summary_tab_tmp <- data.frame(react_name = all_react_names,
                                num_pairs_tested = num_pairs_tested,
                                num_pairs_regulated = num_pairs_regulated,
                                cancer = cancer_tmp)
  react_regulated_summary <- rbind(react_regulated_summary, summary_tab_tmp)
}
react_regulated_summary$regulated_ratio <- react_regulated_summary$num_pairs_regulated/react_regulated_summary$num_pairs_tested
react_regulated_summary2rank <- react_regulated_summary %>%
  filter(num_pairs_regulated >= 5)


react_name_tmp <- "453279\tMitotic G1-G1/S phases"
react_name_tmp <- "380270\tRecruitment of mitotic centrosome proteins and complexes"
react_name_tmp <- "975138\tTRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation"
react_name_tmp <- "128021\tCytokine Signaling in Immune system"
react_name_tmp <- "5218921\tVEGFR2 mediated cell proliferation"
react_name_tmp <- "166208\tmTORC1-mediated signalling"
react_name_tmp <- "380972\tEnergy dependent regulation of mTOR by LKB1-AMPK"
react_name_tmp <- "165159\tmTOR signalling"
react_name_tmp <- "109704\tPI3K Cascade"
react_name_tmp <- "68882\tMitotic Anaphase"
react_name_tmp <- "68877\tMitotic Prometaphase"
react_name_tmp <- "187687\tSignalling to ERKs"
react_name_tmp <- "1236394\tSignaling by ERBB4"
react_name_tmp <- "3700989\tTranscriptional Regulation by TP53"
react_name_tmp <- "109704\tPI3K Cascade"
react_name_tmp <- "198753\tERK/MAPK targets"
react_name_tmp <- "450282\tMAPK targets/ Nuclear events mediated by MAP kinases"
react_name_tmp <- "177929\tSignaling by EGFR"
react_name_tmp <- "1280215\tCytokine Signaling in Immune system"
react_name_tmp <- "3700989\tTranscriptional Regulation by TP53"
react_name_tmp <- "1430728\tMetabolism"
react_name_tmp <- "453274\tMitotic G2-G2/M phases"
react_name_tmp <- "167044\tSignalling to RAS"
react_name_tmp <- "449147\tSignaling by Interleukins"
react_name_tmp <- "2219528\tPI3K/AKT Signaling in Cancer"
react_name_tmp <- "179812\tGRB2 events in EGFR signaling"
react_name_tmp <- "179812\tSHC1 events in EGFR signaling"
react_name_tmp <- "453276\tRegulation of mitotic cell cycle"
react_name_tmp <- "2219530\tConstitutive Signaling by Aberrant PI3K in Cancer"

cancer_tmp <- "CCRCC"
# cancer_tmp <- "UCEC"
# cancer_tmp <- "CO"
# cancer_tmp <- "OV"
# cancer_tmp <- "BRCA"

regression_react_tmp <- regression %>%
  filter(regulated == T) %>%
  filter(Cancer == cancer_tmp) %>%
  filter(SELF == "trans") %>%
  filter(GENE %in% REACT[[react_name_tmp]]) %>%
  filter(SUB_GENE %in% REACT[[react_name_tmp]])

# react only-substrate in pathway  ---------------------------------------------------------------
all_react_names <- names(REACT)
substrate_in_react_summary <- NULL
for (cancer_tmp in unique(regression$Cancer)) {
  num_pairs_tested <- vector(mode = "numeric", length = length(all_react_names))
  num_pairs_regulated <- num_pairs_tested
  count <- 0
  for (react_name_tmp in all_react_names) {
    regression_tested_tmp <- regression %>%
      filter(Cancer == cancer_tmp) %>%
      filter(SELF == "trans") %>%
      filter(!(GENE %in% REACT[[react_name_tmp]])) %>%
      filter(SUB_GENE %in% REACT[[react_name_tmp]])
    
    regression_regulated_tmp <- regression %>%
      filter(regulated == T) %>%
      filter(Cancer == cancer_tmp) %>%
      filter(SELF == "trans") %>%
      filter(!(GENE %in% REACT[[react_name_tmp]])) %>%
      filter(SUB_GENE %in% REACT[[react_name_tmp]])
    
    count <- count + 1
    num_pairs_tested[count] <- regression_tested_tmp %>%
      nrow()
    num_pairs_regulated[count] <- regression_regulated_tmp %>%
      nrow()
  }
  
  summary_tab_tmp <- data.frame(react_name = all_react_names,
                                num_pairs_tested = num_pairs_tested,
                                num_pairs_regulated = num_pairs_regulated,
                                cancer = cancer_tmp)
  substrate_in_react_summary <- rbind(substrate_in_react_summary, summary_tab_tmp)
}
substrate_in_react_summary$regulated_ratio <- substrate_in_react_summary$num_pairs_regulated/substrate_in_react_summary$num_pairs_tested
substrate_in_react_summary2rank <- substrate_in_react_summary %>%
  filter(num_pairs_regulated >= 5)

react_name_tmp <- "1234158\tRegulation of gene expression by Hypoxia-inducible Factor"
cancer_tmp <- "CCRCC"

regression_react_tmp <- regression %>%
  filter(regulated == T) %>%
  filter(Cancer == cancer_tmp) %>%
  filter(SELF == "trans") %>%
  filter((GENE %in% REACT[[react_name_tmp]]) | (SUB_GENE %in% REACT[[react_name_tmp]]))

# react substrate or kinase in pathway  ---------------------------------------------------------------
all_react_names <- names(REACT)
substrate_or_kinase_in_react_summary <- NULL
for (cancer_tmp in unique(regression$Cancer)) {
  num_pairs_tested <- vector(mode = "numeric", length = length(all_react_names))
  num_pairs_regulated <- num_pairs_tested
  count <- 0
  for (react_name_tmp in all_react_names) {
    regression_tested_tmp <- regression %>%
      filter(Cancer == cancer_tmp) %>%
      filter(SELF == "trans") %>%
      filter((GENE %in% REACT[[react_name_tmp]]) | (SUB_GENE %in% REACT[[react_name_tmp]]))
    
    
    regression_regulated_tmp <- regression %>%
      filter(regulated == T) %>%
      filter(Cancer == cancer_tmp) %>%
      filter(SELF == "trans") %>%
      filter((GENE %in% REACT[[react_name_tmp]]) | (SUB_GENE %in% REACT[[react_name_tmp]]))
    
    count <- count + 1
    num_pairs_tested[count] <- regression_tested_tmp %>%
      nrow()
    num_pairs_regulated[count] <- regression_regulated_tmp %>%
      nrow()
  }
  
  summary_tab_tmp <- data.frame(react_name = all_react_names,
                                num_pairs_tested = num_pairs_tested,
                                num_pairs_regulated = num_pairs_regulated,
                                cancer = cancer_tmp)
  substrate_or_kinase_in_react_summary <- rbind(substrate_or_kinase_in_react_summary, summary_tab_tmp)
}
substrate_or_kinase_in_react_summary$regulated_ratio <- substrate_or_kinase_in_react_summary$num_pairs_regulated/substrate_or_kinase_in_react_summary$num_pairs_tested
substrate_or_kinase_in_react_summary$regulated_ratio_string <- paste0(signif(substrate_or_kinase_in_react_summary$regulated_ratio, digits = 2), "(", substrate_or_kinase_in_react_summary$num_pairs_regulated, "/", substrate_or_kinase_in_react_summary$num_pairs_tested, ")")
substrate_or_kinase_in_react_summary2rank <- substrate_or_kinase_in_react_summary %>%
  filter(num_pairs_regulated >= 5)

react_name_tmp <- "1234158\tRegulation of gene expression by Hypoxia-inducible Factor"
react_name_tmp <- "166208\tmTORC1-mediated signalling"
react_name_tmp <- "71406\tPyruvate metabolism and Citric Acid (TCA) cycle"
react_name_tmp <- "442755\tActivation of NMDA receptor upon glutamate binding and postsynaptic events"
react_name_tmp <- "1445148\tTranslocation of GLUT4 to the plasma membrane"
react_name_tmp <- "70171\tGlycolysis"
react_name_tmp <- "380972\tEnergy dependent regulation of mTOR by LKB1-AMPK"
react_name_tmp <- "75105\tFatty Acyl-CoA Biosynthesis"
react_name_tmp <- "380953\tRegulation of Rheb GTPase activity by AMPK"
react_name_tmp <- "1592230\tMitochondrial biogenesis"
react_name_tmp <- "4839726\tChromatin organization"
react_name_tmp <- "3247509\tChromatin modifying enzymes"
react_name_tmp <- "2993913\tClearance of Nuclear Envelope Membranes from Chromatin"
react_name_tmp <- "2559584\tFormation of Senescence-Associated Heterochromatin Foci (SAHF)"
react_name_tmp <- "194138\tSignaling by VEGF"
react_name_tmp <- "194313	VEGF ligand-receptor interactions"
react_name_tmp <- "1059683	Interleukin-6 signaling"
react_name_tmp <- "1640170	Cell Cycle"
react_name_tmp <- "1280215	Cytokine Signaling in Immune system"
react_name_tmp <- "380108	Chemokine receptors bind chemokines"
react_name_tmp <- "912526	Interleukin receptor SHC signaling"
react_name_tmp <- "389513	CTLA4 inhibitory signaling"
react_name_tmp <- "170834	Signaling by TGF-beta Receptor Complex"
react_name_tmp <- "975138\tTRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation"
react_name_tmp <- "109704\tPI3K Cascade"
react_name_tmp <- "199418	Negative regulation of the PI3K/AKT network"
react_name_tmp <- "1250342	PI3K events in ERBB4 signaling"
react_name_tmp <- "2219530	Constitutive Signaling by Aberrant PI3K in Cancer"
react_name_tmp <- "380972	Energy dependent regulation of mTOR by LKB1-AMPK"
react_name_tmp <- "5654699	SHC-mediated cascade:FGFR2"
react_name_tmp <- "5655304	Signaling by FGFR2 mutants"
react_name_tmp <- "5654700	FRS-mediated FGFR2 signaling"
react_name_tmp <- "1250347	SHC1 events in ERBB4 signaling"
react_name_tmp <- "167044	Signalling to RAS"
react_name_tmp <- "4839735	AXIN mutants destabilize the destruction complex, activating WNT signaling"
react_name_tmp <- "73929	Base-Excision Repair, AP Site Formation"
react_name_tmp <- "110302	Formation of transcription-coupled NER (TC-NER) repair complex"
react_name_tmp <- "68952	DNA replication initiation"
react_name_tmp <- "3700989	Transcriptional Regulation by TP53"
react_name_tmp <- "453276	Regulation of mitotic cell cycle"
react_name_tmp <- "201681	TCF dependent signaling in response to WNT"
react_name_tmp <- "195721	Signaling by Wnt"
react_name_tmp <- "180336	SHC1 events in EGFR signaling"
react_name_tmp <- "167044	Signalling to RAS"
react_name_tmp <- "442742	CREB phosphorylation through the activation of Ras"
react_name_tmp <- "113501	Inhibition of replication initiation of damaged DNA by RB1/E2F1"
react_name_tmp <- "1280215	Cytokine Signaling in Immune system"
react_name_tmp <- "449147	Signaling by Interleukins"
react_name_tmp <- "1059683	Interleukin-6 signaling"
react_name_tmp <- "168180	TRAF6 Mediated Induction of proinflammatory cytokines"
react_name_tmp <- "389513	CTLA4 inhibitory signaling"
react_name_tmp <- "389948	PD-1 signaling"
react_name_tmp <- "2173791	TGF-beta receptor signaling in EMT (epithelial to mesenchymal transition)"
react_name_tmp <- "389357	CD28 dependent PI3K/Akt signaling"
react_name_tmp <- "442755	Activation of NMDA receptor upon glutamate binding and postsynaptic events"
react_name_tmp <- "69563	p53-Dependent G1 DNA Damage Response"
react_name_tmp <- "69541	Stabilization of p53"
react_name_tmp <- "1963642	PI3K events in ERBB2 signaling"
react_name_tmp <- "109869	RAF/MAP kinase cascade"
react_name_tmp <- "110049	MEK activation"
react_name_tmp <- "187687	Signalling to ERKs"
react_name_tmp <- "69278	Cell Cycle, Mitotic"
react_name_tmp <- "453274	Mitotic G2-G2/M phases"
react_name_tmp <- "157118	Signaling by NOTCH"
react_name_tmp <- "69273	Cyclin A/B1 associated events during G2/M transition"
react_name_tmp <- "73888	Homologous Recombination Repair"
react_name_tmp <- "1227986	Signaling by ERBB2"
react_name_tmp <- "1251932	PLCG1 events in ERBB2 signaling"
react_name_tmp <- "1059683	Interleukin-6 signaling"
react_name_tmp <- "432142	Platelet sensitization by LDL"
react_name_tmp <- "933542	TRAF6 mediated NF-kB activation"
react_name_tmp <- "389357	CD28 dependent PI3K/Akt signaling"
react_name_tmp <- "199418	Negative regulation of the PI3K/AKT network"
react_name_tmp <- "1250347	SHC1 events in ERBB4 signaling"
react_name_tmp <- "1236394	Signaling by ERBB4"
react_name_tmp <- "187706	Signalling to p38 via RIT and RIN"
react_name_tmp <- "450302	activated TAK1 mediates p38 MAPK activation"
react_name_tmp <- "69541	Stabilization of p53"
react_name_tmp <- "73888	Homologous Recombination Repair"
react_name_tmp <- "174143	APC/C-mediated degradation of cell cycle proteins"
react_name_tmp <- "5218921	VEGFR2 mediated cell proliferation"
react_name_tmp <- "194138	Signaling by VEGF"
react_name_tmp <- "110049	MEK activation"
react_name_tmp <- "110056	ERK1 activation"
cancer_tmp <- "CCRCC"
# cancer_tmp <- "UCEC"
# cancer_tmp <- "CO"
# cancer_tmp <- "OV"
# cancer_tmp <- "BRCA"

regression_react_tmp <- regression %>%
  filter(regulated == T) %>%
  filter(Cancer == cancer_tmp) %>%
  filter(SELF == "trans") %>%
  filter((GENE %in% REACT[[react_name_tmp]]) | (SUB_GENE %in% REACT[[react_name_tmp]]))

regression_gene_tmp <- regression %>%
  filter(regulated == T) %>%
  filter(Cancer == cancer_tmp) %>%
  filter(SELF == "trans") %>%
  filter(SUB_GENE == "CDKN2B")

regression_gene_tmp <- regression %>%
  filter(SUB_GENE %in% c("STAT1", "STAT3", "STAT5A", "STAT5B")) %>%
  # filter(GENE == "PFKFB3") %>%
  # filter(GENE == "AURKB") %>%
  # filter(GENE == "AURKA") %>%
  # filter(GENE == "WEE1") %>%
  # filter(GENE == "CHEK1") %>%
  # filter(GENE == "CHEK2") %>%
  # filter(GENE == "CDK1") %>%
  # filter(GENE == "CDK4") %>%
  # filter(GENE == "CDK6") %>%
  # filter(GENE == "PLK1") %>%
  filter(regulated == T) %>%
  filter(Cancer == cancer_tmp) %>%
  filter(SELF == "trans") 


