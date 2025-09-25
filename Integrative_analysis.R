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


# setClass(
#   "NOODAI",
#   representation(
#     token = "character",
#     working_dir = "character",
#     BioGRID_data_file = "character",
#     STRING_data_file = "character",
#     IntAct_data_file = "character",
#     file_DEA_names = "character",
#     phenotype_names = "character",
#     phenotype_comparison = "character",
#     Use_precompiled_database = "numeric",
#     LookUp_table_file = "character",
#     Results_index = "character",
#     edge_file_path = "character",
#     monet_path = "character",
#     Monet_method_string = "character",
#     tmp_bin_folder = "character",
#     CPDB_databases = "character",
#     MONET_background_file = "character",
#     CPDB_database_file = "character",
#     files_edges_path = "character",
#     centralities_file = "character",
#     TF_Database = "character",
#     file_extension = "character",
#     Kinome_database = "character",
#     BioMart_Dataset = "character",
#     Client_email = "character",
#     weight_penalty = "numeric",
#     circos = "list",
#     pathways = "list"
#   )
# )


#Perform the network construction based on all the omics datasets that were uploaded as well as the computation of the nodes centralities

Integrative_Network_analysis <- function(
  working_dir,
  NOODAI_object,
  BioGRID_data_file = NULL,
  STRING_data_file = NULL,
  IntAct_data_file = NULL,
  file_DEA_names,
  phenotype_names,
  phenotype_comparison,
  splicing_file_name = "nnn",
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
  #zip(zipfile =  paste0(Results_dir,'/ResultsZip.zip'), files = Results_dir,root = Results_dir,include_directories = FALSE)
  #utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir,extras = '-j')
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


#Perform the MONET analysis considering the edges files for each comparison of interest

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
  #######################################################################################################################
  #server_psw <- paste0(working_dir,"/Auxiliary/server_psw.txt")
  
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
    # monet_path <- paste0("/Users/riur/.monet/monet")
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

      #command_monet <- paste0("echo ",readLines(server_psw)," | sudo -S -i; ",monet_path," --input=",tmp_bin_folder,"/",ff[i]," ",Monet_method_string," --output=",paste0(tmp_bin_folder,"/MONET"))
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
  #utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir, extras = '-j')
  utils::zip(
    zipfile = paste0(Results_dir,'/ResultsZip.zip'),
    files  = sapply(
      Results_dir,
      FUN = function(x) {
        strsplit(x,"/home/omics/NOODAI_shiny/")[[1]][2]
        # strsplit(x, "/Users/riur/Documents/noodai-web-interface/")[[1]][2]
      }
    )
  )
  gc()
}

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
  #utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir, extras = '-j')
  utils::zip(
    zipfile = paste0(Results_dir,'/ResultsZip.zip'),
    files  = sapply(
      Results_dir, FUN = function(x) {
        strsplit(x,"/home/omics/NOODAI_shiny/")[[1]][2]
        # strsplit(x, "/Users/riur/Documents/noodai-web-interface/")[[1]][2]
      }
    )
  )
  gc()
}


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

  setwd(wkd)

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
  #system(paste0("mv ",Results_dir,"/Additional_Results/Main_Results ",Results_dir,"/"))
  
  #Save the results as a zip file in the results folder
  library(utils)
  gc()
  if(file.exists(paste0(Results_dir,'/ResultsZip.zip'))) {
    file.remove(paste0(Results_dir,'/ResultsZip.zip'))
  }


  #utils::zip(
  #  zipfile = paste0(Results_dir,'/ResultsZip.zip'),
  #  files  = c(
  #    sapply(
  #      paste0(Results_dir,"/Main_Results"),
  #      FUN = function(x) {
  #        strsplit(x,"/home/omics/NOODAI_shiny/")[[1]][2]
  #        # strsplit(x,"/Users/riur/Documents/noodai-web-interface/")[[1]][2]
  #      }
  #    ),
  #    sapply(
  #      paste0(Results_dir,"/Additional_Results"),
  #      FUN = function(x) {
  #        strsplit(x,"/home/omics/NOODAI_shiny/")[[1]][2]
  #        # strsplit(x,"/Users/riur/Documents/noodai-web-interface/")[[1]][2]
  #      }
  #    )
  #  )
  #)
  
library(zip)

files_to_zip <- file.path(Results_dir, c("Main_Results", "Additional_Results"))
zip_path <- file.path(Results_dir, "ResultsZip.zip")
if (file.exists(zip_path)) file.remove(zip_path)
zip::zip(
  zipfile = zip_path,
  files   = files_to_zip,
  mode    = "cherry-pick"  # ensures only selected dirs/files are included
)
  
  
   system(paste0("cp ",paste0(Results_dir,'/ResultsZip.zip')," /srv/NOODAI_V200/inst/app/www/Downloads/",Results_index,".zip"))
  
  # saveRDS(
  #   NOODAI_object,
  #   file.path(
  #     Results_dir,
  #     "NOODAI_object.RDS"
  #   )
  # )

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
