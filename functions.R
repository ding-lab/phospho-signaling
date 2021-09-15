# Yige Wu @ WashU 2020 Apr

# make output directory ---------------------------------------------------
dir_analysis_result <- paste0(dir_base, "/analysis_results/")
makeOutDir_katmai = function(path_script) {
  folders <- strsplit(x = path_script, split = "\\/")[[1]]
  folder_num <- which(folders == "phospho-signaling_analysis") + 1
  dir_analysis_resultnow <- paste(strsplit(paste(folders[folder_num:length(folders)], collapse = "/"), split = "\\.")[[1]][1], sep = "/")
  dir_analysis_resultnow <- paste0(dir_analysis_result, dir_analysis_resultnow, "/")
  dir.create(dir_analysis_resultnow)
  dir_analysis_resultnow_son <- dir_analysis_resultnow
  dirs2make <- NULL
  while (!dir.exists(dir_analysis_resultnow_son)) {
    tmp <- strsplit(dir_analysis_resultnow_son, split = "\\/")[[1]]
    dir_analysis_resultnow_parent <-paste(tmp[-length(tmp)], collapse = "/")
    dir.create(dir_analysis_resultnow_parent)
    dir.create(dir_analysis_resultnow_son)
    dir.create(dir_analysis_resultnow)
    if (!dir.exists(dir_analysis_resultnow_son)) {
      dirs2make[length(dirs2make) + 1] <- dir_analysis_resultnow_son
    }
    dir_analysis_resultnow_son <- dir_analysis_resultnow_parent
  }
  
  if (length(dirs2make) > 0){
    for (i in 1:length(dirs2make)) {
      dir.create(dirs2make[i])
    }
  } 
  return(dir_analysis_resultnow)
}

makeOutDir = function() {
  folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
  folder_num <- which(folders == "phospho-signaling_analysis") + 1
  dir_analysis_resultnow <- paste(strsplit(paste(folders[folder_num:length(folders)], collapse = "/"), split = "\\.")[[1]][1], sep = "/")
  dir_analysis_resultnow <- paste0(dir_analysis_result, dir_analysis_resultnow, "/")
  dir.create(dir_analysis_resultnow)
  dir_analysis_resultnow_son <- dir_analysis_resultnow
  dirs2make <- NULL
  while (!dir.exists(dir_analysis_resultnow_son)) {
    tmp <- strsplit(dir_analysis_resultnow_son, split = "\\/")[[1]]
    dir_analysis_resultnow_parent <-paste(tmp[-length(tmp)], collapse = "/")
    dir.create(dir_analysis_resultnow_parent)
    dir.create(dir_analysis_resultnow_son)
    dir.create(dir_analysis_resultnow)
    if (!dir.exists(dir_analysis_resultnow_son)) {
      dirs2make[length(dirs2make) + 1] <- dir_analysis_resultnow_son
    }
    dir_analysis_resultnow_son <- dir_analysis_resultnow_parent
  }
  
  if (length(dirs2make) > 0){
    for (i in 1:length(dirs2make)) {
      dir.create(dirs2make[i])
    }
  } 
  return(dir_analysis_resultnow)
}

get_mutation_class_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat & maf$Variant_Classification != "Silent",]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  }, value.var = "Variant_Classification", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

generate_variant_aachange_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  # maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- gsub(x = maf$Tumor_Sample_Barcode, pattern = "_T", replacement = "")
  maf$HGVSp_Shorter <- gsub(x = maf$HGVSp_Short, pattern = "p.", replacement = "")
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  }, value.var = "HGVSp_Shorter", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  mut_mat <- mut_mat[genes4mat,]
  return(mut_mat)
}

get_mutation_class_sim_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat & maf$Variant_Classification != "Silent",]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class_sim <- plyr::mapvalues(x = unique(x), 
                                         from = c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site", "Missense_Mutation", "In_Frame_Ins", "In_Frame_Del"),
                                         to = c("Truncation", "Truncation", "Truncation", "Truncation", "Missense", "In_Frame_Ins", "In_Frame_Del"))
    variant_class <- paste0(sort(variant_class_sim), collapse = ",")
    
    return(variant_class)
  }, value.var = "Variant_Classification", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

FDR_by_id_columns <- function(p_vector, id_columns, df) {
  ## make a replicate the id columns and make it charactors
  df_id <- matrix(data = as.character(unlist(df[, id_columns])), nrow = nrow(df), ncol = length(id_columns), byrow = F)
  df_id <- data.frame(df_id); colnames(df_id) <- id_columns
  df_id_combo <- data.frame(table(df_id))
  df_id_combo <- df_id_combo[df_id_combo$Freq > 0,]
  ## give a number for each combo
  df_id_combo$combo_id <- 1:nrow(df_id_combo)
  df_id_combo
  df_id <- merge(df_id, df_id_combo[, c(id_columns, "combo_id")], all.x = T)
  if (any(is.na(df_id$combo_id))) {
    stop()
  }
  
  ## for every combo of values in id columns, adjust a set of pvalues
  fdr_vector <- vector(mode = "numeric", length = length(p_vector)) + NA
  for (i in 1:length(unique(df_id$combo_id))) {
    row2adjust <- (df_id$combo_id == i & !is.na(p_vector))
    fdr_vector[row2adjust] <- p.adjust(p = p_vector[row2adjust], method = "fdr")
  }
  return(fdr_vector)
}
