#The script contains the main functions used by the NOODAI app. The Monet analysis is done explicitly in this script as it is wrapped in an additional function.

# 
#     Copyright © 2025, Empa, Tiberiu Totu and Rafael Riudavets Puig.
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#     Contact: tiberiu.totu@proton.me




#setClass(
#  "NOODAI",
#  representation(
#    token = "character",
#    working_dir = "character",
#    BioGRID_data_file = "ANY",
#    STRING_data_file = "ANY",
#    IntAct_data_file = "ANY",
#    file_DEA_names = "character",
#    phenotype_names = "character",
#    phenotype_comparison = "character",
#    Use_precompiled_database = "numeric",
#    LookUp_table_file = "ANY",
#    Results_index = "character",
#    edge_file_path = "ANY",
#    monet_path = "ANY",
#    Monet_method_string = "character",
#    tmp_bin_folder = "ANY",
#    CPDB_databases = "ANY",
#    MONET_background_file = "ANY",
#    CPDB_database_file = "ANY",
#    files_edges_path = "ANY",
#    centralities_file = "ANY",
#    TF_Database = "ANY",
#    file_extension = "ANY",
#    Kinome_database = "ANY",
#    BioMart_Dataset = "character",
#    Client_email = "ANY",
#    weight_penalty = "numeric",
#    circos = "list",
#    circos_aux = "list",
#    pathways = "list",
#    cytoscape_modules = "list",
#    status = "character"
#  )
#)








#' Construct a multi-omics integrative network and compute node centralities
#'
#' Integrative_Network_analysis function performs integrative network analysis using multiple omics datasets and user-defined parameters. 
#' It loads interaction databases (BioGRID, STRING, IntAct), prepares phenotype-specific data, 
#' computes centrality metrics, and saves the results to a results directory. At the same time, it updates the status of a NOODAI object throughout the process.
#' A pre-compiled interaction database can be used.
#'
#' @param working_dir Path to the working directory containing input files and scripts.
#' @param NOODAI_object A custom object containing NOODAI-related data and state.
#' @param BioGRID_data_file Optional path to the BioGRID data file. Defaults to internal 4.4.218 version database if NULL.
#' @param STRING_data_file Optional path to the STRING data file. Defaults to internal 11.5 version database if NULL.
#' @param IntAct_data_file Optional path to the IntAct data file. Defaults to internal 245 version database if NULL.
#' @param file_DEA_names A vector of full path file names containing differentially expressed analysis (DEA) results.
#' @param phenotype_names A vector of phenotype identifiers to be used in the analysis.
#' @param phenotype_comparison A vector describing phenotype comparisons to perform.
#' @param splicing_file_name Name of the file containing alternative splicing results. Default is "None".
#' @param Use_precompiled_database Logical. If TRUE, use a precompiled database for analysis.
#' @param LookUp_table_file Optional path to an interaction lookup table. Defaults to internal file if NULL.
#' @param Results_index A unique identifier used to label the output directory and result files.
#' @param BioMart_Dataset The BioMart dataset name for gene annotations (e.g., "hsapiens_gene_ensembl").
#' @param weight_penalty A numeric value controlling the penalty applied to network edge weights.
#'
#' @return This function does not return an R object and it calls an internal script 'Network_analysis_BioGRID_STRING_IntAct.R'. It creates a results directory with network analysis outputs and a zip archive of results.
#'
#' @examples
#' Integrative_Network_analysis(
#'   working_dir = "/path/to/scripts",
#'   NOODAI_object = my_object,
#'   file_DEA_names = c("/path/to/scripts/Results_001/tmp/DEA_file1.xlsx", "/path/to/scripts/Results_001/tmp/DEA_file2.xlsx"),
#'   phenotype_names = c("Cancer", "Control"),
#'   phenotype_comparison = c("CancerVsControl"),
#'   Use_precompiled_database = TRUE,
#'   Results_index = "001",
#'   BioMart_Dataset = "hsapiens_gene_ensembl",
#'   weight_penalty = 0.1
#' )

Integrative_Network_analysis <- function(
  working_dir,
  NOODAI_object,
  BioGRID_data_file = NULL,
  STRING_data_file = NULL,
  IntAct_data_file = NULL,
  file_DEA_names,
  phenotype_names,
  phenotype_comparison,
  splicing_file_name = "None",
  Use_precompiled_database,
  LookUp_table_file = NULL,
  Results_index,
  BioMart_Dataset,
  weight_penalty
) {
  
  ###Network analysis
  
  #Set the working directory as the one where the scripts are found and the results directory as the one newly created
  setwd(working_dir)
  Results_dir <- paste0(getwd(),"/Results_",Results_index)

  NOODAI_object@status = "Network analysis"
  NOODAI_object_path = file.path(Results_dir, "NOODAI_object.RDS")
  if (file.exists(NOODAI_object_path)) {file.remove(NOODAI_object_path)}
  saveRDS(
    NOODAI_object,
    NOODAI_object_path
  )
  
  #Load the corresponding databases if they were not provided
  if(is.null(BioGRID_data_file)){BioGRID_data_file <- paste0(working_dir,"/Databases/BIOGRID-MV-Physical-4.4.218.mitab.txt")}
  if(is.null(STRING_data_file)){STRING_data_file <- paste0(working_dir,"/Databases/9606.protein.links.full.v11.5.txt")}
  if(is.null(IntAct_data_file)){IntAct_data_file <- paste0(working_dir,"/Databases/IntAct_24022022_07Confidence.txt")}
  if(is.null(LookUp_table_file)){LookUp_table_file <- paste0(working_dir,"/Databases/Interaction_Lookup_table.txt")}
  if(is.null(BioMart_Dataset)){BioMart_Dataset <- 'hsapiens_gene_ensembl'}
  BiomaRT_selected_organisms <- paste0(working_dir,"/Databases/BiomaRT_selected_organisms.txt")
  Uniprot_Ncbi_map <- paste0(working_dir,"/Databases/Unip_NCBI_map.txt")
  ChEBI_map <- paste0(working_dir,"/Databases/ChEBI_entries_map.txt")
  
  #Run the network centrality computation
  source('Network_analysis_BioGRID_STRING_IntAct.R')
  
  status <- NULL
  attempt <- 1
  while(is.null(status) && attempt <= 5) {
    attempt <- attempt + 1
    try(
      status <- Network_analysis(
        working_dir,
        Results_dir,
        BioGRID_data_file,
        STRING_data_file,
        IntAct_data_file,
        file_DEA_names,
        phenotype_names,
        phenotype_comparison,
        splicing_file_name,
        Use_precompiled_database,
        LookUp_table_file,
        BioMart_Dataset,
        BiomaRT_selected_organisms,
        Uniprot_Ncbi_map,
        ChEBI_map,
        weight_penalty
      )
    )
    gc()
  }
  
  #Save the results as a zip file in the results folder
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip.zip'))) {
      file.remove(paste0(Results_dir,'/ResultsZip.zip'))
  }

  utils::zip(
    zipfile = paste0(Results_dir,'/ResultsZip.zip'),
    files  = sapply(
      Results_dir,
      FUN = function(x) {
        strsplit(x, "/home/omics/NOODAI_shiny/")[[1]][2]
      }
    )
  )
}


#' Perform MONET module decomposition analysis using edge files
#'
#' This function runs the MONET module detection tool on precomputed edge files representing biological networks.
#' It prepares temporary directories, manages input/output paths, executes MONET via command-line system calls,
#' and saves the results to a specified results folder. It updates the status of a NOODAI object throughout the process.
#'
#' @param working_dir Path to the working directory containing input data and scripts.
#' @param NOODAI_object A custom NOODAI object that holds workflow status and metadata.
#' @param edge_file_path Optional relative path to a folder containing network edge files. 
#'        Defaults to 'edge_files_PPINetworks/Symbol' if NULL.
#' @param monet_path Optional full path to the MONET binary. If NULL, defaults to '/home/omics/.monet/monet'.
#' @param Monet_method_string A string specifying the MONET method and parameters to be used (e.g., '--method=k1').
#' @param tmp_bin_folder Optional path to a temporary folder for intermediate MONET files.
#'        Defaults to a new folder under the current results directory.
#' @param Results_index A string used to uniquely label the output directory and zip archive.
#'
#' @return This function does not return an R object. It writes MONET module results to a subdirectory 
#' under the results folder and creates a zip archive of the full analysis.
#'
#' @details
#' The function copies input edge files into a temporary binary folder, executes MONET on each file,
#' and copies the results into a structured output directory. If MONET fails to detect modules, 
#' the function halts and updates the NOODAI object status accordingly. 
#' Docker containerization is automatically appended to the MONET command.
#'
#' @examples
#' MONET_analysis(
#'   working_dir = "/path/to/scripts",
#'   NOODAI_object = my_object,
#'   edge_file_path = "edge_files_PPINetworks/Symbol",
#'   monet_path = "/home/omics/.monet/monet",
#'   Monet_method_string = "--method=k1",
#'   tmp_bin_folder = "/tmp/monet_tmp",
#'   Results_index = "001"
#' )

MONET_analysis <- function(
  working_dir,
  NOODAI_object,
  edge_file_path = NULL,
  monet_path = NULL,
  Monet_method_string,
  tmp_bin_folder = NULL,
  Results_index
) {
  
  #Set the working directory as the one where the scripts are found and the results directory as the one newly created
  setwd(working_dir)
  Results_dir <- paste0(getwd(), "/Results_", Results_index)

  
  NOODAI_object@status = "MONET analysis"
  NOODAI_object_path = file.path(Results_dir, "NOODAI_object.RDS")
  if (file.exists(NOODAI_object_path)) {file.remove(NOODAI_object_path)}
  saveRDS(
    NOODAI_object,
    NOODAI_object_path
  )

  #Select the preferred tmp folder and use the default Symbol edge files if none are provided
  if(is.null(edge_file_path)) {
    edge_file_path <- paste0("edge_files_PPINetworks/Symbol")
  }
  if(is.null(tmp_bin_folder)) {
    tmp_bin_folder <- paste0(Results_dir,"/tmp_Monet")
  }
  #######################################################################################################################
  if(is.null(monet_path)) {
    monet_path <- paste0("/home/omics/.monet/monet")

  }
  
  if(dir.exists(paste0(Results_dir,"/",edge_file_path))) {
    
    
    dir.create(paste0(Results_dir,"/MONET_analysis"), recursive = T)
    save_dir <- paste0(Results_dir,"/MONET_analysis")
    
    system(paste("mkdir",tmp_bin_folder))
    dir_MONET <- paste0(Results_dir,"/",edge_file_path,"/*")
    dir_MONET_wsl <- dir_MONET
    system(paste("cp -r", dir_MONET_wsl, tmp_bin_folder))
    system(paste("mkdir", paste0(tmp_bin_folder, "/MONET")))
    Monet_method_string = paste0(Monet_method_string, " --container=docker")
    
    ff <- list.files(paste0(Results_dir, "/", edge_file_path))
    for (i in 1:length(ff)){

      command_monet <- paste0(
        monet_path,
        " --input=",
        tmp_bin_folder,
        "/",
        ff[i],
        " ",
        Monet_method_string,
        " --output=",
        paste0(tmp_bin_folder,"/MONET")
      )
      
      system(command_monet)

    }
    
    system(paste0("cp -r ", paste0(tmp_bin_folder,"/MONET")," ", save_dir))
    system(paste0("rm -r ", tmp_bin_folder))
    
  }else{

    stop("Edge file directory does not exist")

  }

  monet_dir = file.path(save_dir, "MONET")
  if (length(list.files(monet_dir)) == 0) {
    NOODAI_object@status = "Stopped - No modules detected."
    if (file.exists(NOODAI_object_path)) {file.remove(NOODAI_object_path)}
    saveRDS(
      NOODAI_object,
      NOODAI_object_path
    )
    stop("No modules detected")
  }
  
  #Save the results as a zip file in the results folder
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip.zip'))) {
    file.remove(paste0(Results_dir,'/ResultsZip.zip'))
  }

  utils::zip(
    zipfile = paste0(Results_dir,'/ResultsZip.zip'),
    files  = sapply(
      Results_dir,
      FUN = function(x) {
        strsplit(x,"/home/omics/NOODAI_shiny/")[[1]][2]
      }
    )
  )
  gc()
}




#' Extract enriched pathways for MONET-decomposed network modules
#'
#' This function identifies and extracts signaling pathways associated with MONET network modules.
#' It uses ConsensusPathDB (CPDB) annotations to map genes from detected modules to known pathways. 
#' The function manages database paths, updates the NOODAI object status, and exports the results 
#' into the working directory structure.
#'
#' @param working_dir Path to the main project or scripts directory.
#' @param NOODAI_object A custom NOODAI object used for tracking workflow status and state.
#' @param CPDB_databases A character vector indicating which CPDB pathway databases to query (e.g., '"KEGG"', '"Reactome"').
#' @param MONET_background_file Optional path to a background gene set file used for enrichment. 
#'        Defaults to 'Background_total.xlsx' within the results folder if NULL.
#' @param phenotype_names A vector of phenotype identifiers.
#' @param phenotype_comparison A vector describing phenotype comparisons performed.
#' @param CPDB_database_file Optional path to the ConsensusPathDB reference file (e.g., 'CPDB_pathways_genes.tab').
#'        If NULL, a default file in the '/Databases' folder is used.
#' @param Results_index A unique identifier used to build the results directory path.
#' @param BioMart_Dataset The BioMart dataset name to use for gene annotations (e.g., '"hsapiens_gene_ensembl"').
#'
#' @return This function does not return an R object. It saves pathway enrichment results in the 
#' MONET analysis directory and archives them in a zip archive under the results folder.
#'
#' @details
#' This function calls an internal script ('MONET_pathways_extraction_CPDB.R') to perform 
#' enrichment using CPDB annotations. Gene identifier mapping is done using UniProt, ChEBI, and BioMart data.
#' The results include enriched pathways per module and comparison, which are written into structured 
#' directories under 'Results_<index>/MONET_analysis/MONET'.
#'
#' @examples
#' MONET_pathways(
#'   working_dir = "/path/to/scripts",
#'   NOODAI_object = my_object,
#'   CPDB_databases = c("KEGG", "Reactome"),
#'   phenotype_names = c("Cancer", "Control"),
#'   phenotype_comparison = c("CancerVsControl"),
#'   Results_index = "001",
#'   BioMart_Dataset = "hsapiens_gene_ensembl"
#' )

#Extract the pathways associated with the identified modules for each comparison of interest.
MONET_pathways <- function(
  working_dir,
  NOODAI_object,
  CPDB_databases,
  MONET_background_file = NULL,
  phenotype_names,
  phenotype_comparison,
  CPDB_database_file = NULL,
  Results_index,
  BioMart_Dataset
) {
  
  #Set the working directory as the one where the scripts are found and the results directory as the one newly created
  wkd <- working_dir
  setwd(working_dir)
  
  if(is.null(CPDB_database_file)) {
    CPDB_database_file <- paste0(working_dir,"/Databases/CPDB_pathways_genes.tab")
  }
  if(is.null(BioMart_Dataset)) {
    BioMart_Dataset <- 'hsapiens_gene_ensembl'
  }
  BiomaRT_selected_organisms <- paste0(working_dir,"/Databases/BiomaRT_selected_organisms.txt")
  Uniprot_Ncbi_map <- paste0(working_dir,"/Databases/Unip_NCBI_map.txt")
  ChEBI_map <- paste0(working_dir,"/Databases/ChEBI_entries_map.txt")
  
  Results_dir <- paste0(getwd(),"/Results_",Results_index)
  working_dir <- paste0(Results_dir,"/MONET_analysis/MONET")

  NOODAI_object@status = "MONET pathway extraction"
  NOODAI_object_path = file.path(Results_dir, "NOODAI_object.RDS")
  if (file.exists(NOODAI_object_path)) {file.remove(NOODAI_object_path)}
  saveRDS(
    NOODAI_object,
    NOODAI_object_path
  )
  
  if(is.null(MONET_background_file)) {
    MONET_background_file <- paste0(Results_dir,"/Background_total.xlsx")
  }

  cat(file = stderr(), "I'm here!\n")
  
  source("MONET_pathways_extraction_CPDB.R")
  
  status <- NULL
  attempt <- 1
  while(is.null(status) && attempt <= 5) {
    attempt <- attempt + 1
    try(
      status <- MONET_Pathways_extraction(
        working_dir,
        CPDB_database_file,
        CPDB_databases,
        MONET_background_file,
        phenotype_names,
        phenotype_comparison,
        BioMart_Dataset,
        BiomaRT_selected_organisms,
        Uniprot_Ncbi_map,
        ChEBI_map
      )
    )
  }

  setwd(wkd)
  Results_dir <- paste0(getwd(),"/Results_",Results_index)

  #Save the results as a zip file in the results folder
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip.zip'))) {
    file.remove(paste0(Results_dir,'/ResultsZip.zip'))
  }
  utils::zip(
    zipfile = paste0(Results_dir,'/ResultsZip.zip'),
    files  = sapply(
      Results_dir, FUN = function(x) {
        strsplit(x,"/home/omics/NOODAI_shiny/")[[1]][2]
      }
    )
  )
  gc()
}




#' Generate network visualizations, summary reports, and Cytoscape networks
#'
#' This function generates key graphical outputs and summary files following MONET module detection and signaling pathway enrichment analysis. It includes the
#' creation of circos plots, pathway diagrams, centrality visualizations, and optional Cytoscape networks. It also 
#' compiles results into organized folders and produces an HTML report using R Markdown.
#'
#' @param working_dir Path to the main working directory containing input data and scripts.
#' @param NOODAI_object A custom NOODAI object used for status tracking and results aggregation.
#' @param phenotype_names A vector of phenotype identifiers used in the analysis.
#' @param phenotype_comparison A vector of phenotype comparisons used for group analysis (e.g., '"TumorVsControl"').
#' @param files_edges_path Optional path to network edge files. Defaults to 'edge_files_PPINetworks/Symbol'.
#' @param centralities_file Optional path to the CINNA centrality result file. Defaults to 'PPINetworks_centralities_values_CINNA_Total.xlsx'.
#' @param TF_Database Optional path to transcription factor annotation database. Defaults to '_TF.txt' in the Databases folder.
#' @param file_extension A suffix used to identify edge and centrality files (e.g., '"Total"'). Defaults to '"Total"' if NULL.
#' @param Results_index A unique identifier for labeling result directories and files.
#' @param Kinome_database Optional path to a kinome annotation file. Defaults to 'kinome.txt' in the Databases folder.
#' @param BioMart_Dataset BioMart dataset name for gene annotation (e.g., '"hsapiens_gene_ensembl"'). Defaults to this value if NULL.
#' @param Cytoscape_style_file Optional path to a Cytoscape style XML file. If NULL, a default is used.
#' @param file_DEA_names A vector of input omics (DEA) full path file names used in the analysis.
#' @param flag_run_cytoscape Logical or numeric flag (0 or 1) indicating whether Cytoscape network should be generated.
#'
#' @return Returns the updated 'NOODAI_object' with visualizations and modules embedded (e.g., Cytoscape module paths).
#'
#' @details
#' This function executes a series of steps including:
#' - Extracting edges and centrality values for each module.
#' - Creating circos plots and pathway enrichment visualizations.
#' - Generating a PDF/HTML summary report via R Markdown.
#' - Optionally exporting styled Cytoscape networks for each phenotype comparison and module.
#' - Organizing all outputs into 'Main_Results' and 'Additional_Results' folders.
#' - Compressing results into a ZIP archive for download or export.
#'
#' Temporary files, logs, and modules are automatically managed and cleaned up between steps. This function relies
#' on multiple external R scripts for network analysis, visualization, and export, including Extract_edge_files_for_clusters.R,
#' Extract_TopCentralitiesInModuels.R, Circos_plots.R, Monet_cluster_pathways_image_joint.R, and optionally Cytoscape_Network_Generation.R.
#'
#' @examples
#' Circos_and_auxiliary(
#'   working_dir = "/path/to/scripts",
#'   NOODAI_object = my_object,
#'   phenotype_names = c("Cancer", "Control"),
#'   phenotype_comparison = c("CancerVsControl"),
#'   file_DEA_names = c("/path/to/scripts/Results_001/tmp/DEA1.csv", "/path/to/scripts/Results_001/tmp/DEA2.csv"),
#'   Results_index = "001",
#'   flag_run_cytoscape = 1
#' )


#Create the main graphical representations and summary documents
Circos_and_auxiliary <- function(
  working_dir,
  NOODAI_object,
  phenotype_names,
  phenotype_comparison,
  files_edges_path = NULL,
  centralities_file = NULL,
  TF_Database = NULL,
  file_extension = NULL,
  Results_index,
  Kinome_database = NULL,
  BioMart_Dataset,
  Cytoscape_style_file = NULL,
  file_DEA_names,
  flag_run_cytoscape = NULL
) {
  
  #Save separately the directory where the scripts are found as the working directory will change multiple times
  wkd <- working_dir
  
  if(is.null(TF_Database)) {TF_Database <- paste0(working_dir,"/Databases/_TF.txt")}
  if(is.null(BioMart_Dataset)) {BioMart_Dataset <- 'hsapiens_gene_ensembl'}
  BiomaRT_selected_organisms <- paste0(working_dir,"/Databases/BiomaRT_selected_organisms.txt")
  
  setwd(wkd)
  Results_dir <- paste0(getwd(),"/Results_",Results_index)

  NOODAI_object@status = "Plot and report generation"
  NOODAI_object_path = file.path(Results_dir, "NOODAI_object.RDS")
  if (file.exists(NOODAI_object_path)) {file.remove(NOODAI_object_path)}
  saveRDS(
    NOODAI_object,
    NOODAI_object_path
  )
  
  working_dir <- paste0(Results_dir,"/MONET_analysis/MONET")
  
  if(is.null(files_edges_path)) {
    files_edges_path <- paste0(Results_dir,"/edge_files_PPINetworks/Symbol")
  }
  
  if(is.null(file_extension)) {file_extension = "Total"}
  
  source("Extract_edge_files_for_clusters.R")
  
  
  
  status <- NULL
  attempt <- 1
  while( is.null(status) && attempt <= 5 ) {
    attempt <- attempt + 1
    try(
      status <- Extract_edges_tables(
        working_dir,
        phenotype_names,
        phenotype_comparison,
        files_edges_path,
        file_extension
      )
    )
  }
  
  if(is.null(centralities_file)) {
    centralities_file <- paste0(Results_dir,"/PPINetworks_centralities_values_CINNA_Total.xlsx")
  }
  
  setwd(wkd)
  
  
  source("Extract_TopCentralitiesInModuels.R")
  
  status <- NULL
  attempt <- 1
  while(is.null(status) && attempt <= 5) {
    attempt <- attempt + 1
    try(
      status <- Extract_CentralitiesValues(
        working_dir,
        centralities_file,
        phenotype_names,
        file_extension
      )
    )
  }
  
  setwd(wkd)
  
  source("Circos_plots.R")
  
  status <- NULL
  attempt <- 1
  while(is.null(status) && attempt <= 5) {
    attempt <- attempt + 1
    try({
      circos_results <- Cricos_plots(
        working_dir,
        NOODAI_object,
        files_edges_path,
        TF_Database,
        file_extension,
        BioMart_Dataset,
        BiomaRT_selected_organisms)

      status = circos_results[[1]]
      NOODAI_object = circos_results[[2]]
      }
    )
    
  }
  
  setwd(wkd)
  
  
  pathways_dir <- paste0(working_dir,"/Signaling_Pathways")
  
  source("Monet_cluster_pathways_image_joint.R")
  
  status <- NULL
  attempt <- 1
  while(is.null(status) && attempt <= 5) {
    attempt <- attempt + 1
    try({
      # here for pathways
      enrichment_results <- MONET_cluster_pathways_joint_image(
        pathways_dir,
        file_extension,
        NOODAI_object
      )
      
      status = enrichment_results[[1]]
      NOODAI_object = enrichment_results[[2]]
    })
  }
  
  setwd(wkd)
  
  if(is.null(Kinome_database)) {Kinome_database <- paste0(wkd,"/Databases/kinome.txt")}
  
  tmp_omics_names <- list.files(paste0(Results_dir,'/tmp_FilesIni'), recursive = T, full.names = F)
  tmp_omics_names <- basename(tmp_omics_names)
  tmp_omics_names <- strsplit(tmp_omics_names,'\\.')
  tmp_omics_names <- sapply(tmp_omics_names,"[[",1)
  
  rmarkdown::render(
    input  = "Networks_Analyisis_results_PDF_generation.Rmd",
    params = list(
      results_dir = Results_dir,
      centralities_file_name = centralities_file,
      omics_file_names = tmp_omics_names,
      TF_file = TF_Database,
      kinome_file = Kinome_database,
      pathways_folder = pathways_dir,
      BiomaRT_organisms_file = BiomaRT_selected_organisms
    ),
    output_file = paste0(Results_dir,"/Results_Interpretation")
  )

  if(is.null(flag_run_cytoscape)) {flag_run_cytoscape<-0}
 if(flag_run_cytoscape == 1){
  try({

    setwd(wkd)

    if(is.null(Cytoscape_style_file)){
      Cytoscape_style_file <- paste0(wkd,"/Databases/Style_NOODAI_Networks.xml")
    }
    Cytoscape_style_file <- normalizePath(path.expand(Cytoscape_style_file))
    
    Results_dir <- paste0(wkd,"/Results_",Results_index)
    save_dir <- paste0(Results_dir,"/MONET_analysis/MONET/Cytoscape_Networks")
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

    file_DEA_names_short <- tools::file_path_sans_ext(basename(file_DEA_names))

    Centrality_file <- paste0(Results_dir,"/PPINetworks_centralities_values_CINNA_Total.xlsx")
    path_to_Clusters <-  paste0(Results_dir,"/MONET_analysis/MONET/Edges_tables")
    full_network_edge_path <- paste0(Results_dir,"/edge_files_PPINetworks/Symbol")

    source("Cytoscape_Network_Generation.R")

    for (Centrality_file_sheet in phenotype_comparison){
      
      Sys.sleep(10)
      save_dir_aux <- paste0(save_dir,"/",Centrality_file_sheet)
      dir.create(save_dir_aux, recursive = TRUE, showWarnings = FALSE)
      setwd(save_dir_aux)
      Sys.sleep(10)
      Edge_file <- paste0(full_network_edge_path,"/","Total_",Centrality_file_sheet,".txt")
      save_name_network <- paste0("No_Module_",Centrality_file_sheet,"_",Results_index)

      status <- local({
        attempt <- 1
        local_status <- NULL
        while (is.null(local_status) && attempt <= 2) {
          attempt <- attempt + 1
          try({
            local_status <- Cytoscape_networks(
            Centrality_file,
            Centrality_file_sheet,
            file_DEA_names_short,
            Edge_file,
            Edge_file_sheet = NULL,
            Cytoscape_style_file,
            save_name_network
          )
        }, silent = FALSE)
        }
        local_status
      })

      all_edge_files <- list.files(path_to_Clusters, full.names = TRUE)
      Edge_file <- all_edge_files[grepl(Centrality_file_sheet, all_edge_files)]
      Edge_file_sheets <- openxlsx::getSheetNames(Edge_file)
      if(length(Edge_file_sheets) > 10){ Edge_file_sheets <- Edge_file_sheets[1:10] }

      for (Edge_file_sheet in Edge_file_sheets){
        save_name_network <- paste0(Edge_file_sheet,"_",Centrality_file_sheet,"_",Results_index)
        status <- local({
          attempt <- 1
          local_status <- NULL
          while (is.null(local_status) && attempt <= 2) {
            attempt <- attempt + 1
            try({
              local_status <- Cytoscape_networks(
              Centrality_file,
              Centrality_file_sheet,
              file_DEA_names_short,
              Edge_file,
              Edge_file_sheet,
              Cytoscape_style_file,
              save_name_network
            )
          }, silent = FALSE)
          }
          local_status
        })
          
      }

      status <- NULL
      attempt <- 1
      while( is.null(status) && attempt <= 2 ) {
        attempt <- attempt + 1
        try({
          status <- Generate_cytoscape_Legend(
            Centrality_file,
            Centrality_file_sheet,
            file_DEA_names_short,
            "Legend_cytoscape.svg"
          )
        }, silent = FALSE)
      }
      Sys.sleep(5)
    }
    
    modules_ls <- lapply(
      file.path(save_dir, phenotype_comparison),
      function(x) {
        list.files(
          path = x,
          pattern = ".svg",
          full.names = TRUE,
          recursive = TRUE
        )
      }
    )
    names(modules_ls) <- phenotype_comparison
    NOODAI_object@cytoscape_modules <- modules_ls

  })
}

  system(paste0("rm ", Results_dir, "/Main_Results/ -r"))
  system(paste0("rm ", Results_dir, "/Additional_Results/ -r"))
  system(paste0("mkdir ", Results_dir, "/Main_Results"))
  system(paste0("mkdir ", Results_dir, "/Main_Results/Modules_Signaling_Pathways"))
  system(paste0("mkdir ", Results_dir, "/Main_Results/Functional_Enrichment_Plots"))
  system(paste0("mkdir ", Results_dir, "/Main_Results/Circular_Diagrams"))
  system(paste0("mkdir ", Results_dir, "/Additional_Results"))

  system(
    paste0(
      "cp ",
      Results_dir,
      "/Results_Interpretation.html ",
      Results_dir,
      "/Main_Results/"
    )
  )

  system(
    paste0(
      "cp ",
      Results_dir,
      "/PPINetworks_centralities_values_CINNA_Total.xlsx ",
      Results_dir,
      "/Main_Results/"
    )
  )
  
  system(
    paste0(
      "cp ",
      Results_dir,
      "/MONET_analysis/MONET/Signaling_Pathways/Images/* ",
      Results_dir,
      "/Main_Results/Functional_Enrichment_Plots/"
    )
  )

  system(
    paste0(
      "cp ",
      Results_dir,
      "/MONET_analysis/MONET/Signaling_Pathways/*_Total_* ",
      Results_dir,
      "/Main_Results/Modules_Signaling_Pathways/"
    )
  )

  system(
    paste0(
      "cp ",
      Results_dir,
      "/MONET_analysis/MONET/Centralities_modules_links/Images/* ",
      Results_dir,
      "/Main_Results/Circular_Diagrams/"
    )
  )

  system(
    paste0(
      "cp ",
      Results_dir,
      "/* ",
      Results_dir,
      "/Additional_Results/ -r"
    )
  )

  system(
    paste0(
      "rm ",
      Results_dir,
      "/Additional_Results/Additional_Results -r"
    )
  )

  system(
    paste0(
      "rm ",
      Results_dir,
      "/Additional_Results/Main_Results -r"
    )
  )

  system(
    paste0(
      "rm ",
      Results_dir,
      "/Additional_Results/ResultsZip.zip"
    )
  )

  #Save the results as a zip file in the results folder
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip.zip'))) {
    file.remove(paste0(Results_dir,'/ResultsZip.zip'))
  }
  #utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir, extras = '-j')
  utils::zip(
    zipfile = paste0(Results_dir,'/ResultsZip.zip'),
    files  = c(
      sapply(
        paste0(Results_dir,"/Main_Results"),
        FUN = function(x) {
          strsplit(x,"/home/omics/NOODAI_shiny/")[[1]][2]
        }
      ),
      sapply(
        paste0(Results_dir,"/Additional_Results"),
        FUN = function(x) {
          strsplit(x,"/home/omics/NOODAI_shiny/")[[1]][2]
        }
      )
    )
  )
  
   system(paste0("cp ",paste0(Results_dir,'/ResultsZip.zip')," /srv/NOODAI_V200/inst/app/www/Downloads/",Results_index,".zip"))


  NOODAI_object@status = "Finished"
  NOODAI_object_path = file.path(Results_dir, "NOODAI_object.RDS")
  if (file.exists(NOODAI_object_path)) {file.remove(NOODAI_object_path)}
  saveRDS(
    NOODAI_object,
    NOODAI_object_path
  )

  gc()
  return(NOODAI_object)
  
}
