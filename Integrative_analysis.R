#The script contains the main functions used by the NOODAI app. The Monet analysis is done explicitly in this script as is wrapped in an additional function.

# 
#     Copyright © 2024, Empa, Tiberiu Totu.
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
#     Contact: tiberiu.totu@empa.ch





#Perform the network construction based on all the omics datasets that were uploaded as well as the computation of the nodes centralities

Integrative_Network_analysis <- function(working_dir,BioGRID_data_file=NULL,STRING_data_file=NULL,IntAct_data_file=NULL,file_DEA_names,
                                         phenotype_names,phenotype_comparison,splicing_file_name="nnn",
                                         Use_precompiled_database,LookUp_table_file=NULL,Results_index,BioMart_Dataset){
  ###Network analysis
  
  #Set the working directory as the one where the scripts are found and the results directory as the one newly created
  setwd(working_dir)
  Results_dir <- paste0(getwd(),"/Results_",Results_index)
  
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
  while( is.null(status) && attempt <= 5 ) {
    attempt <- attempt + 1
    try(
      status <- Network_analysis(working_dir,Results_dir,BioGRID_data_file,
                                 STRING_data_file,IntAct_data_file,file_DEA_names,
                                 phenotype_names,phenotype_comparison,splicing_file_name,
                                 Use_precompiled_database,LookUp_table_file,BioMart_Dataset,
								 BiomaRT_selected_organisms,Uniprot_Ncbi_map,ChEBI_map)
    )
    gc()
  }
  
  #Save the results as a zip file in the results folder
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip.zip'))){file.remove(paste0(Results_dir,'/ResultsZip.zip'))}
  #zip(zipfile =  paste0(Results_dir,'/ResultsZip.zip'), files = Results_dir,root = Results_dir,include_directories = FALSE)
  #utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir,extras = '-j')
  utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = sapply(Results_dir,FUN=function(x){strsplit(x,"/home/omics/Linux_test/")[[1]][2]}))
}


#Perform the MONET analysis considering the edges files for each comparison of interest

MONET_analysis <- function(working_dir,edge_file_path=NULL,monet_path=NULL,Monet_method_string,tmp_bin_folder=NULL,Results_index){
  
  #Set the working directory as the one where the scripts are found and the results directory as the one newly created
  setwd(working_dir)
  Results_dir <- paste0(getwd(),"/Results_",Results_index)
  #######################################################################################################################
  #server_psw <- paste0(working_dir,"/Auxiliary/server_psw.txt")
  
  #Select the preferred tmp folder and use the default Symbol edge files if none are provided
  if(is.null(edge_file_path)){edge_file_path <- paste0("edge_files_PPINetworks/Symbol")}
  if(is.null(tmp_bin_folder)){tmp_bin_folder <- paste0(Results_dir,"/tmp_Monet")}
  #######################################################################################################################
  if(is.null(monet_path)){monet_path <- paste0("/home/omics/.monet/monet")}
  
  if(dir.exists(paste0(Results_dir,"/",edge_file_path))){
    
    
    dir.create(paste0(Results_dir,"/MONET_analysis"))
    save_dir <- paste0(Results_dir,"/MONET_analysis")
    
    system(paste("mkdir",tmp_bin_folder))
    dir_MONET <- paste0(Results_dir,"/",edge_file_path,"/*")
    dir_MONET_wsl <- dir_MONET
    system(paste("cp -r",dir_MONET_wsl,tmp_bin_folder))
    system(paste("mkdir",paste0(tmp_bin_folder,"/MONET")))
    Monet_method_string = paste0(Monet_method_string," --container=docker")
    
    ff <- list.files(paste0(Results_dir,"/",edge_file_path))
    for (i in 1:length(ff)){
      #command_monet <- paste0("echo ",readLines(server_psw)," | sudo -S -i; ",monet_path," --input=",tmp_bin_folder,"/",ff[i]," ",Monet_method_string," --output=",paste0(tmp_bin_folder,"/MONET"))
      command_monet <- paste0(monet_path," --input=",tmp_bin_folder,"/",ff[i]," ",Monet_method_string," --output=",paste0(tmp_bin_folder,"/MONET"))
      
      system(command_monet)
    }
    system(paste0("cp -r ", paste0(tmp_bin_folder,"/MONET")," ", save_dir))
    system(paste0("rm -r ", tmp_bin_folder))
    
  }else{
    stop("Edge file directory does not exist")
  }
  
  #Save the results as a zip file in the results folder
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip.zip'))){file.remove(paste0(Results_dir,'/ResultsZip.zip'))}
  #utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir, extras = '-j')
  utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = sapply(Results_dir,FUN=function(x){strsplit(x,"/home/omics/Linux_test/")[[1]][2]}))
  
  gc()
  
}



#Extract the pathways associated with the identified modules for each comparison of interest.
MONET_pathways <- function(working_dir,CPDB_databases,MONET_background_file=NULL,phenotype_names,
                           phenotype_comparison,CPDB_database_file=NULL,Results_index,BioMart_Dataset){
  
  #Set the working directory as the one where the scripts are found and the results directory as the one newly created
  wkd <- working_dir
  setwd(working_dir)
  
  if(is.null(CPDB_database_file)){CPDB_database_file <- paste0(working_dir,"/Databases/CPDB_pathways_genes.tab")}
  if(is.null(BioMart_Dataset)){BioMart_Dataset <- 'hsapiens_gene_ensembl'}
  BiomaRT_selected_organisms <- paste0(working_dir,"/Databases/BiomaRT_selected_organisms.txt")
  Uniprot_Ncbi_map <- paste0(working_dir,"/Databases/Unip_NCBI_map.txt")
  ChEBI_map <- paste0(working_dir,"/Databases/ChEBI_entries_map.txt")
  
  Results_dir <- paste0(getwd(),"/Results_",Results_index)
  working_dir <- paste0(Results_dir,"/MONET_analysis/MONET")
  
  if(is.null(MONET_background_file)){
    MONET_background_file <- paste0(Results_dir,"/Background_total.xlsx")
  }
  
  source("MONET_pathways_extraction_CPDB.R")
  
  
  status <- NULL
  attempt <- 1
  while( is.null(status) && attempt <= 5 ) {
    attempt <- attempt + 1
    try(
      status <- MONET_Pathways_extraction(working_dir,CPDB_database_file,CPDB_databases,
                                          MONET_background_file,phenotype_names,phenotype_comparison,BioMart_Dataset,BiomaRT_selected_organisms,Uniprot_Ncbi_map,ChEBI_map)
    )
  }
  
  setwd(wkd)
  Results_dir <- paste0(getwd(),"/Results_",Results_index)
  #Save the results as a zip file in the results folder
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip.zip'))){file.remove(paste0(Results_dir,'/ResultsZip.zip'))}
  #utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir, extras = '-j')
  utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = sapply(Results_dir,FUN=function(x){strsplit(x,"/home/omics/Linux_test/")[[1]][2]}))
  
  gc()
  
}


#Create the main graphical representations and summary documents
Circos_and_auxiliary <- function(working_dir,phenotype_names,phenotype_comparison,files_edges_path=NULL,
                                 centralities_file=NULL,TF_Database=NULL,file_extension=NULL,Results_index,
                                 Kinome_database=NULL,BioMart_Dataset){
  
  #Save separately the directory where the scripts are found as the working directory will change multiple times
  wkd <- working_dir
  
  if(is.null(TF_Database)){TF_Database <- paste0(working_dir,"/Databases/_TF.txt")}
  if(is.null(BioMart_Dataset)){BioMart_Dataset <- 'hsapiens_gene_ensembl'}
  BiomaRT_selected_organisms <- paste0(working_dir,"/Databases/BiomaRT_selected_organisms.txt")
  
  setwd(wkd)
  Results_dir <- paste0(getwd(),"/Results_",Results_index)
  
  working_dir <- paste0(Results_dir,"/MONET_analysis/MONET")
  
  if(is.null(files_edges_path)){
    files_edges_path <- paste0(Results_dir,"/edge_files_PPINetworks/Symbol")
  }
  
  if(is.null(file_extension)){file_extension = "Total"}
  
  source("Extract_edge_files_for_clusters.R")
  
  
  
  status <- NULL
  attempt <- 1
  while( is.null(status) && attempt <= 5 ) {
    attempt <- attempt + 1
    try(
      status <- Extract_edges_tables(working_dir,phenotype_names,phenotype_comparison,files_edges_path,file_extension)
    )
  }
  
  if(is.null(centralities_file)){
    centralities_file <- paste0(Results_dir,"/PPINetworks_centralities_values_CINNA_Total.xlsx")
  }
  
  setwd(wkd)
  
  
  source("Extract_TopCentralitiesInModuels.R")
  
  status <- NULL
  attempt <- 1
  while( is.null(status) && attempt <= 5 ) {
    attempt <- attempt + 1
    try(
      status <- Extract_CentralitiesValues(working_dir,centralities_file,phenotype_names,file_extension)
    )
  }
  
  setwd(wkd)
  
  source("Circos_plots.R")
  
  status <- NULL
  attempt <- 1
  while( is.null(status) && attempt <= 5 ) {
    attempt <- attempt + 1
    try(
      status <- Cricos_plots(working_dir,files_edges_path,TF_Database,file_extension,BioMart_Dataset,BiomaRT_selected_organisms)
    )
  }
  
  setwd(wkd)
  
  
  pathways_dir <- paste0(working_dir,"/Signaling_Pathways")
  
  source("Monet_cluster_pathways_image_joint.R")
  
  status <- NULL
  attempt <- 1
  while( is.null(status) && attempt <= 5 ) {
    attempt <- attempt + 1
    try(
      status <- MONET_cluster_pathways_joint_image(pathways_dir,file_extension)
    )
  }
  
  setwd(wkd)
  
  if(is.null(Kinome_database)){Kinome_database <- paste0(wkd,"/Databases/kinome.txt")}
  
  tmp_omics_names <- list.files(paste0(Results_dir,'/tmp_FilesIni'))
  tmp_omics_names <- strsplit(tmp_omics_names,'\\.')
  tmp_omics_names <- sapply(tmp_omics_names,"[[",1)
  
  rmarkdown::render(
    input  = "Networks_Analyisis_results_PDF_generation.Rmd"
    , params = list(
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
  system(paste0("rm ",Results_dir,"/Main_Results/ -r"))
  system(paste0("rm ",Results_dir,"/Additional_Results/ -r"))
  system(paste0("mkdir ",Results_dir,"/Main_Results"))
  system(paste0("mkdir ",Results_dir,"/Main_Results/Modules_Signaling_Pathways"))
  system(paste0("mkdir ",Results_dir,"/Main_Results/Pathway_Images"))
  system(paste0("mkdir ",Results_dir,"/Main_Results/Circular_Diagrams"))
  system(paste0("mkdir ",Results_dir,"/Additional_Results"))
  system(paste0("cp ",Results_dir,"/Results_Interpretation.html ",Results_dir,"/Main_Results/"))
  system(paste0("cp ",Results_dir,"/PPINetworks_centralities_values_CINNA_Total.xlsx ",Results_dir,"/Main_Results/"))
  system(paste0("cp ",Results_dir,"/MONET_analysis/MONET/Signaling_Pathways/Images/* ",Results_dir,"/Main_Results/Pathway_Images/"))
  system(paste0("cp ",Results_dir,"/MONET_analysis/MONET/Signaling_Pathways/*_Total_* ",Results_dir,"/Main_Results/Modules_Signaling_Pathways/"))
  system(paste0("cp ",Results_dir,"/MONET_analysis/MONET/Centralities_modules_links/Images/* ",Results_dir,"/Main_Results/Circular_Diagrams/"))
  system(paste0("cp ",Results_dir,"/* ",Results_dir,"/Additional_Results/ -r"))
  system(paste0("rm ",Results_dir,"/Additional_Results/Additional_Results -r"))
  system(paste0("rm ",Results_dir,"/Additional_Results/Main_Results -r"))
  system(paste0("rm ",Results_dir,"/Additional_Results/ResultsZip.zip"))
  #system(paste0("mv ",Results_dir,"/Additional_Results/Main_Results ",Results_dir,"/"))
  
  #Save the results as a zip file in the results folder
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip.zip'))){file.remove(paste0(Results_dir,'/ResultsZip.zip'))}
  #utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir, extras = '-j')
  utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = c(sapply(paste0(Results_dir,"/Main_Results"),FUN=function(x){strsplit(x,"/home/omics/Linux_test/")[[1]][2]}),
             sapply(paste0(Results_dir,"/Additional_Results"),FUN=function(x){strsplit(x,"/home/omics/Linux_test/")[[1]][2]})) )
  
  gc()
  
  
  
  
}
