# Yige Wu @ WashU 2019 July
# test the trans correlation between individual kinase/phosphatase phosphosite with substrate phosphsoites
## to limit the search and facilitate interpretation, we only use function-known phosophosites on kinase/phosphatase

# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)

# input version number ----------------------------------------------------
version_num <- 1

# set variables -----------------------------------------------------------
least_samples <- 20# least number of samples with complete data for each model

# input regulatory sites --------------------------------------------------
regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/Resources/Phospho_databases/PhosphositePlus/May_28_2019/Regulatory_sites", data.table = F)
regulatory_sites <- regulatory_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1]) %>%
  mutate(phosphosite = paste0(GENE, "_", SUB_MOD_RSD))
disease_sites <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/Resources/Phospho_databases/PhosphositePlus/Mar_04_2019/Disease-associated_sites", data.table = F)
disease_sites <- disease_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1]) %>%
  mutate(phosphosite = paste0(GENE, "_", SUB_MOD_RSD)) %>%
  filter(DISEASE %in% c("kidney cancer", "clear cell kidney cancer"))


# input cancer driver kinases ---------------------------------------------
kinase_driver_table <- readxl::read_xlsx(path = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/resources/Gene_Lists/Fleuren2016/nrc.2015.18-s1.xlsx", sheet = "Suppl Table 1 - Kinase list")
kinases_genomics_altered <- kinase_driver_table %>%
  filter(grepl(x = `Mutational drivers tumor typed`, pattern = "GBM") | grepl(x = `CNA cancer drivers tumor typed`, pattern = "GBM") | grepl(x = `Oncogenic fusion drivers tumor typed`, pattern = "GBM")) %>%
  select(`Gene symbol`)
kinases_genomics_altered <- unlist(kinases_genomics_altered$`Gene symbol`)
kinases_genomics_altered

# input ccRCC important genes ---------------------------------------------
ks_tab <- load_es_pro_with_predicted_table()
kinases_manually_curated <- intersect(c(SMGs[["GBM"]], CNGs[["GBM"]]), ks_tab$GENE)
kinases_manually_curated
genes2test <- unique(c(kinases_genomics_altered, kinases_manually_curated))
genes2test

# run regression -----------------------------------------------------------------
cancer <- "GBM"
pipeline_type <- "PGDAC"
sample_type <- "tumor"
norm_type <- "d4_d6"
cptac_phase <- "CPTAC3"

pro_data <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "d4")
pho_data <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "d6")

pho_data <- data.frame(pho_data)
pro_data <- data.frame(pro_data)

pho_data %>% head()
pho_rsd_split <- data.frame(SUBSTRATE = pho_data$Gene, SUB_MOD_RSD = pho_data$Phosphosite)
colnames(pho_rsd_split) <- c("SUBSTRATE", "SUB_MOD_RSD")

if (any(is.na(pho_rsd_split$SUBSTRATE) | is.na(pho_rsd_split$SUB_MOD_RSD))) {
  stop("pho_rsd_split weird")
}

sub_mod_rsd <- as.vector(pho_rsd_split$SUB_MOD_RSD)
samples <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene"))]
samples <- intersect(samples, colnames(pro_data))

# make table for pairs to test --------------------------------------------
pairs2test_tab <- load_es_pro_with_predicted_table()
pairs2test_tab <- merge(pairs2test_tab, pho_rsd_split, by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"))
pairs2test_tab <- merge(pairs2test_tab, 
                        pho_rsd_split %>%
                          mutate(ENZ_MOD_RSD = SUB_MOD_RSD) %>%
                          select(SUBSTRATE, ENZ_MOD_RSD), 
                        by.x = c("GENE"), by.y = c("SUBSTRATE"))
pairs2test_tab <- pairs2test_tab %>%
  mutate(ENZ_phosphosite = paste0(GENE, "_", ENZ_MOD_RSD)) %>%
  mutate(SUB_phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  # filter(ENZ_phosphosite %in% regulatory_sites$phosphosite | ENZ_phosphosite %in% disease_sites$phosphosite) %>%
  filter(SUB_GENE %in% pro_data$Gene) %>%
  filter(SUB_GENE != GENE)

# for demo purpose just test SMG------------------------------
pairs2test_tab <- pairs2test_tab %>%
  filter(GENE %in% genes2test)

# initiate ----------------------------------------------------------------
# calculate the length of trans table
npairs2test <- nrow(pairs2test_tab)

# looping over kinases for trans pairs -----------------------------------------------------------------
# initiating the table for trans
vec_char <- vector(mode = "character", length = npairs2test)
vec_num <- vector(mode = "numeric", length = npairs2test) + NA

FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
Cancer <- vec_char;
Transcript_sub <- vec_char; Transcript_enz <- vec_char;
model <- vec_char;
Size <- vec_num;
P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
sd_pro_kin <- vec_num;sd_pro_sub <- vec_num;sd_pho_kin <- vec_num; sd_pho_sub <- vec_num;

for (j in 1:npairs2test){
  print(paste0(j, "/", npairs2test))
  
  enzyme <- pairs2test_tab[j, "GENE"]
  enz_mod_rsd <- pairs2test_tab[j, "ENZ_MOD_RSD"]
  substrate <- pairs2test_tab[j, "SUB_GENE"]
  sub_mod_rsd <- pairs2test_tab[j, "SUB_MOD_RSD"]
  
  pho_sub <- pho_data[pho_data$Gene == substrate & pho_data$Phosphosite == sub_mod_rsd, samples]
  pho_enz <- pho_data[pho_data$Gene == enzyme & pho_data$Phosphosite == enz_mod_rsd, samples]
  pro_sub <- pro_data[pro_data$Gene == substrate,samples]
  
  #prepare regression data for model2
  df4regression <- data.frame(t(rbind(pho_sub,pro_sub,pho_enz)))
  colnames(df4regression) <- c("pho_sub","pro_sub","pho_enz")
  data_complete <- df4regression[complete.cases(df4regression),]
  size <- nrow(data_complete)
  if(size > least_samples ){
    fit2 <- lm(pho_sub ~ pro_sub + pho_enz, data = df4regression)
    
    P_pro_sub[j] <- c(coef(summary(fit2))[2,4])
    coef_pro_sub[j] <- fit2$coefficients[2]
    P_pho_kin[j] <- c(coef(summary(fit2))[3,4])
    coef_pho_kin[j] <- fit2$coefficients[3]
    Size[j] <- size
    sd_pro_sub[j] <- sd(data_complete$pro_sub)
    sd_pho_kin[j] <- sd(data_complete$pho_enz)
    sd_pho_sub[j] <- sd(data_complete$pho_sub)
  }
}
table_trans <- data.frame(pairs2test_tab,
                          FDR_pro_kin,FDR_pro_sub,FDR_pho_kin,
                          coef_pro_kin,coef_pro_sub,coef_pho_kin,
                          Size,
                          P_pro_kin,P_pro_sub,P_pho_kin,
                          sd_pro_kin,sd_pro_sub, sd_pho_kin,sd_pho_sub)
table_trans <- table_trans %>%
  filter(!is.na(Size) == T) %>%
  mutate(model = "pho_sub~pro_sub+pho_enz_site") %>%
  mutate(Cancer = cancer) %>%
  mutate(SELF = "trans") %>%
  mutate(pair = paste(ENZ_phosphosite, SUB_GENE, SUB_MOD_RSD, sep = ":"))

table_trans <- annotate_enzyme_type(regression = table_trans, kinases = kinases, phosphatases = phosphatases)
table_trans <- table_trans %>%
  filter(enzyme_type != "")

## write out
tn = paste0(makeOutDir(), "regression_", cptac_phase , "_", cancer, "_", sample_type, "_", pipeline_type, "_", norm_type, ".txt")
write.table(table_trans, file=tn, quote=F, sep = '\t', row.names = FALSE)

table_trans_filtered <- table_trans %>%
  filter(Size >= 20) 
table_trans_filtered$FDR_pho_kin <- FDR_by_id_columns(p_vector = table_trans_filtered$P_pho_kin, id_columns = c("Cancer", "SELF"), df = table_trans_filtered)
table_trans_filtered$FDR_pro_sub <- FDR_by_id_columns(p_vector = table_trans_filtered$P_pro_sub, id_columns = c("Cancer", "SELF"), df = table_trans_filtered)
table_trans_filtered$ENZ_phosphosite.is_regulatory <- (table_trans_filtered$ENZ_phosphosite %in% regulatory_sites$phosphosite)
table_trans_filtered$ENZ_phosphosite.is_disease_related <- (table_trans_filtered$ENZ_phosphosite %in% disease_sites$phosphosite)
table_trans_filtered$SUB_phosphosite.is_regulatory <- (table_trans_filtered$SUB_phosphosite %in% regulatory_sites$phosphosite)
table_trans_filtered$SUB_phosphosite.is_disease_related <- (table_trans_filtered$SUB_phosphosite %in% disease_sites$phosphosite)
# 
## write out
tn = paste0(makeOutDir(), "regression_", cptac_phase , "_", cancer, "_", sample_type, "_", pipeline_type, "_", norm_type, "_nonNA20", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_num, ".txt")
write.table(table_trans_filtered, file=tn, quote=F, sep = '\t', row.names = FALSE)
