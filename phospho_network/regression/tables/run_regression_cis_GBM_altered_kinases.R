# Yige Wu @ WashU Oct 2019
# look at correlations of kinase and downstream substrates phosphorylation status
# pho_sub~pro_sub+pho_kin(callapsed)

# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)

# set variables -----------------------------------------------------------
least_samples <- 5# least number of samples with complete data for each model

# input version number ----------------------------------------------------
version_num <- 1

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
kinase_cis <- genes2test

# loop -----------------------------------------------------------------
cancer <- "GBM"
pipeline_type <- "PGDAC"
sample_type <- "tumor"
norm_type <- "d4_d6"
cptac_phase <- "CPTAC3"

pro_data <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "d4")
pho_data <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "d6")

pho_rsd_split <- data.frame(SUBSTRATE = pho_data$Gene, SUB_MOD_RSD = pho_data$Phosphosite)
colnames(pho_rsd_split) <- c("SUBSTRATE", "SUB_MOD_RSD")

if (any(is.na(pho_rsd_split$SUBSTRATE) | is.na(pho_rsd_split$SUB_MOD_RSD))) {
  stop()
}

transcripts <- NA
sub_mod_rsd <- as.vector(pho_rsd_split$SUB_MOD_RSD)
samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("Gene"))]
samples <- intersect(samples, colnames(pho_data))
samples <- intersect(samples, colnames(pho_gdata))

# initiate ----------------------------------------------------------------
# calculate the length of cis table
ncis <- 0
for (gene in kinase_cis) {
  ncis <- ncis + length(which(pho_rsd_split$SUBSTRATE==gene))
}

# looping over kinases for cis pairs -----------------------------------------------------------------
# initiating the table for cis
vec_char <- vector(mode = "character", length = ncis)
vec_num <- vector(mode = "numeric", length = ncis) + NA
KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
Size <- vec_num;P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
sd_pro_kin <- vec_num;sd_pro_sub <- vec_num;sd_pho_kin <- vec_num; sd_pho_sub <- vec_num;

count <- 0
for (kinase in kinase_cis){
  #for (kinase in "ERBB2"){#test
  # find enzyme_type expression level for the kinase
  pro_kin <- pro_data[pro_data$Gene == kinase, samples]
  
  substrate <- kinase
  if(nrow(pro_kin) != 0){
    s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
    
    # go through all the phosphosites
    for (i in s_pho_table) {
      
      # find phosphorylation level
      pho_sub <- pho_data[i, samples]
      
      #prepare regression data for model1
      data1 <- data.frame(t(rbind(pho_sub,pro_kin)))
      colnames(data1) <- c("pho_sub","pro_kin")
      data_complete <- data1[complete.cases(data1),]
      size <- nrow(data_complete)
      if( size > least_samples ){
        fit1 <- lm(pho_sub ~ pro_kin,data = data1)
        
        count <- count + 1
        KINASE[count] <- kinase
        SUBSTRATE[count] <- substrate
        SUB_MOD_RSD[count] <- sub_mod_rsd[i]
        transcript[count] <- transcripts[i]
        P_pro_kin[count] <- c(coef(summary(fit1))[2,4])
        coef_pro_kin[count] <- fit1$coefficients[2]
        Size[count] <- size
        sd_pro_kin[count] <- sd(data_complete$pro_kin)
        sd_pho_sub[count] <- sd(data_complete$pho_sub)
      }
    }
  }
}

table_cis <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,
                        FDR_pro_kin,FDR_pro_sub,FDR_pho_kin,
                        coef_pro_kin,coef_pro_sub,coef_pho_kin,
                        Cancer,transcript,model,Size,
                        P_pro_kin,P_pro_sub,P_pho_kin,
                        sd_pro_kin,sd_pro_sub,sd_pho_kin,sd_pho_sub)


# integrate table from all the models(need to repetite again for another cancer dataset) --------------------------------------------------------
table_cis$model <- "pho_sub~pro_kin"

# combine table
table <- table_cis[!is.na(table_cis$P_pro_kin),]
# mark cancer and self-regulation
table$Cancer <- cancer
table$self <- as.character(table$KINASE) == as.character(table$SUBSTRATE)
table$SELF <- "trans"; table$SELF[table$self] <- "cis"
table$pair <- paste(table$KINASE,table$SUBSTRATE,table$SUB_MOD_RSD,sep = ":")

## write out
tn = paste0(makeOutDir(), "regression_cis_regulation_", cptac_phase , "_", cancer, "_", sample_type, "_", pipeline_type, "_", norm_type, ".txt")
write.table(table, file=tn, quote=F, sep = '\t', row.names = FALSE)

table_filtered <- table %>%
  filter(Size >= 20) 


table_filtered$FDR_pro_kin <- FDR_by_id_columns(p_vector = table_filtered$P_pro_kin, id_columns = c("Cancer", "SELF"), df = table_filtered)

## write out
tn = paste0(makeOutDir(), "regression_cis_regulation_", cptac_phase , "_", cancer, "_", sample_type, "_", pipeline_type, "_", norm_type, "_nonNA20", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_num, ".txt")
write.table(table_filtered, file=tn, quote=F, sep = '\t', row.names = FALSE)


