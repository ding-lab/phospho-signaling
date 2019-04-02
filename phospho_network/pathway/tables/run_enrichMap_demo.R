

# set parameters ----------------------------------------------------------
params <- list()
params$pval_thresh <- 0.05
params$fdr_thresh <- 0.05
params$working_dir <- "Downloads"
params$expression_file <- "Supplementary_Table6_TCGA_OV_RNAseq_expression.txt"

# 2.1 Load in required libraries ------------------------------------------


#install required R and bioconductor packages
tryCatch(expr = { library("RCurl")}, 
         error = function(e) {  install.packages("RCurl")}, 
         finally = library("RCurl"))

#use library
tryCatch(expr = { library("limma")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("limma")}, 
         finally = library("limma"))
tryCatch(expr = { library("Biobase")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("Biobase")}, 
         finally = library("Biobase"))
tryCatch(expr = { library("ggplot2")}, 
         error = function(e) { install.packages("ggplot2")}, 
         finally = library("ggplot2"))

#For creating json and communicating with cytoscape
tryCatch(expr = { library("httr")}, 
         error = function(e) { install.packages("httr")}, 
         finally = library("httr"))
tryCatch(expr = { library("RJSONIO")}, 
         error = function(e) { install.packages("RJSONIO")}, 
         finally = library("RJSONIO"))


# 2.2 Configurable Parameters ---------------------------------------------
#path to GSEA jar 
# In order to run GSEA automatically you need to speciry the path to the gsea jar file.
gsea_jar <- params$gsea_jar
gsea_jar <- "Downloads/gsea-3.0.jar"
#Gsea takes a long time to run.  If you have already run GSEA manually or previously there is no need to re-run GSEA.  Make sure the 
# gsea results are in the current directory and the notebook will be able to find them and use them.
run_gsea = params$run_gsea

#navigate to the directory where you put the downloaded protocol files.
working_dir <- params$working_dir

# leave blank if you want the notebook to discover the gsea directory for itself
#gsea_directory = paste(working_dir,"Mesen_vs_Immuno.GseaPreranked.1497635459262",sep="/") 
gsea_directory = params$gsea_directory

# analysis_name <- params$analysis_name
analysis_name <- "test"

# rnk_file <- params$rnk_file
rnk_file <- "Supplementary_Table2_MesenvsImmuno_RNASeq_ranks.rnk"
expression_file <- params$expression_file
# classes_file <- params$classes_file


# 2.3 Download the latest pathway definition file -------------------------
gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"

#list all the files on the server
filenames = getURL(gmt_url)
tc = textConnection(filenames)
contents = readLines(tc)
close(tc)

#get the gmt that has all the pathways and does not include terms inferred from electronic annotations(IEA)
#start with gmt file that has pathways only
rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",
              contents, perl = TRUE)
gmt_file = unlist(regmatches(contents, rx))

dest_gmt_file <- file.path(working_dir,paste("Supplementary_Table3_",gmt_file,sep="") )

download.file(
  paste(gmt_url,gmt_file,sep=""),
  destfile=dest_gmt_file
)


# 2.4 Run GSEA ------------------------------------------------------------
run_gsea <- T
if(run_gsea){
  command <- paste("/Users/yigewu/miniconda3/bin/java  -Xmx1G -cp",gsea_jar,  "xtools.gsea.GseaPreranked -gmx", dest_gmt_file, "-rnk" ,file.path(working_dir,rnk_file), "-collapse false -nperm 1000 -permute gene_set -scoring_scheme weighted -rpt_label ",analysis_name,"  -num 100 -plot_top_x 20 -rnd_seed 12345  -set_max 200 -set_min 15 -zip_report false -out" ,working_dir, "-gui false > gsea_output.txt",sep=" ")
  system(command)
  system("java -h")
  system("java -version")
  system("which java")

}
working_dir
dest_gmt_file
rnk_file


# 2.5 Get the name of the GSEA output directory ---------------------------
gsea_directory <- ""
if(gsea_directory == ""){
  gsea_directories <- list.files(path = working_dir, pattern = "\\.GseaPreranked")
  
  #get the details on the files
  details = file.info(file.path(getwd(),working_dir,gsea_directories))
  #order according to newest to oldest
  details = details[with(details, order(as.POSIXct(mtime),decreasing = TRUE)), ]
  
  #use the newest file:
  gsea_output_dir <- row.names(details)[1]
} else {
  gsea_output_dir <- gsea_directory
}
gsea_output_dir


# 2.6 Launch Cytoscape ----------------------------------------------------


# 2.7 Set up connection from R to cytoscape -------------------------------
# Basic settings
port.number = 1234
base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")

#print(base.url)

version.url = paste(base.url, "version", sep="/")
cytoscape.open = TRUE

tryCatch(expr = { GET(version.url)}, 
         error = function(e) { return (cytoscape.open = FALSE)}, finally =function(r){ return(cytoscape.open = TRUE)})

if(!cytoscape.open){
  #try and launch cytoscape
  print("Cytoscape is not open.  Please launch cytoscape.")
} else{
  cytoscape.version =  GET(version.url)
  cy.version = fromJSON(rawToChar(cytoscape.version$content))
  print(cy.version)
  
}

# 2.8 Create an Enrichment map --------------------------------------------
#use easy cyRest library to communicate with cytoscape.
tryCatch(expr = { library("RCy3")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("RCy3")}, finally = library("RCy3"))

#defined threshold for GSEA enrichments (need to be strings for cyrest call)
pvalue_gsea_threshold <- params$pval_thresh
qvalue_gsea_threshold <- params$fdr_thresh

similarity_threshold <- "0.375"
similarity_metric = "COMBINED"

cur_model_name <- analysis_name

gsea_results_path <- file.path(gsea_output_dir,"edb")
gsea_results_filename <- file.path(gsea_results_path,"results.edb")

#although there is a gmt file in the gsea edb results directory it have been filtered to 
#contain only genes represented in the expression set.  If you use this fltered file you 
#will get different pathway connectivity depending on the dataset being used.  We recommend 
#using original gmt file used for the gsea analysis and not the filtered one in the results directory.
gmt_gsea_file <- file.path(getwd(),dest_gmt_file)
gsea_ranks_file <- file.path(gsea_results_path,list.files(gsea_results_path,pattern=".rnk"))

#######################################
#create EM
current_network_name <- paste(cur_model_name,pvalue_gsea_threshold,qvalue_gsea_threshold,sep="_")

em_command = paste('enrichmentmap build analysisType="gsea" gmtFile=',gmt_gsea_file,
                   'pvalue=',pvalue_gsea_threshold, 'qvalue=',qvalue_gsea_threshold,
                   'similaritycutoff=',similarity_threshold,
                   'coefficients=',similarity_metric,'ranksDataset1=', 
                   gsea_ranks_file,'enrichmentsDataset1=',gsea_results_filename, 
                   'filterByExpressions=false',
                   'expressionDataset1=',file.path(getwd(),working_dir,expression_file),
                   sep=" ")

commandsGET
#enrichment map command will return the suid of newly created network.
response <- commandsGET(em_command)

current_network_suid <- 0
#enrichment map command will return the suid of newly created network unless it Failed.  
#If it failed it will contain the word failed
if(grepl(pattern="Failed", response)){
  paste(response)
} else {
  current_network_suid <- response
}
response <- renameNetwork(title=current_network_name, 
                          network = as.numeric(current_network_suid),base.url)

response


# 2.9 Get a screen shot of the initial network. ---------------------------
output_network_file <- file.path(getwd(),"initial_screenshot_network.png")

url_png <- paste(base.url,"networks",current_network_suid,"views/first.png", sep="/")
response <- GET(url=url_png)
writeBin(response$content, output_network_file)
htmltools::img(src = knitr::image_uri(output_network_file), 
               alt = 'Initial Enrichment Map', 
               style = 'margin:0px auto;display:block')
