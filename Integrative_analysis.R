Integrative_Network_analysis <- function(working_dir,BioGRID_data_file=NULL,STRING_data_file=NULL,IntAct_data_file=NULL,file_DEA_names,
                                         phenotype_names,phenotype_comparison,splicing_file_name="nnn",
                                         Use_precompiled_database,LookUp_table_file=NULL,Results_index){
  ###Network analysis
  
  setwd(working_dir)
  Results_dir <- paste0(getwd(),"/Results_",Results_index)
  
  if(is.null(BioGRID_data_file)){BioGRID_data_file <- paste0(working_dir,"/Databases/BIOGRID-MV-Physical-4.4.218.mitab.txt")}
  if(is.null(STRING_data_file)){STRING_data_file <- paste0(working_dir,"/Databases/9606.protein.links.full.v11.5.txt")}
  if(is.null(IntAct_data_file)){IntAct_data_file <- paste0(working_dir,"/Databases/IntAct_24022022_07Confidence.txt")}
  if(is.null(LookUp_table_file)){LookUp_table_file <- paste0(working_dir,"/Databases/Interaction_Lookup_table.txt")}
  
  source('Network_analysis_BioGRID_STRING_IntAct.R')
  
  status <- NULL
  attempt <- 1
  while( is.null(status) && attempt <= 5 ) {
    attempt <- attempt + 1
    try(
      status <- Network_analysis(working_dir,Results_dir,BioGRID_data_file,
                                 STRING_data_file,IntAct_data_file,file_DEA_names,
                                 phenotype_names,phenotype_comparison,splicing_file_name,
                                 Use_precompiled_database,LookUp_table_file)
    )
    gc()
  }
  
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip'))){file.remove(paste0(Results_dir,'/ResultsZip'))}
  #zip(zipfile =  paste0(Results_dir,'/ResultsZip.zip'), files = Results_dir,root = Results_dir,include_directories = FALSE)
  utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir)
}


MONET_analysis <- function(working_dir,edge_file_path=NULL,monet_path=NULL,Monet_method_string,tmp_bin_folder=NULL,Results_index){
  
  
  setwd(working_dir)
  Results_dir <- paste0(getwd(),"/Results_",Results_index)
  #######################################################################################################################
  #server_psw <- paste0(working_dir,"/Auxiliary/server_psw.txt")
  
  if(is.null(edge_file_path)){edge_file_path <- paste0("edge_files_STRINGBioGRIDIntAct/Symbol")}
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
  
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip'))){file.remove(paste0(Results_dir,'/ResultsZip'))}
  utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir)
  
  gc()
  
}




MONET_pathways <- function(working_dir,CPDB_databases,MONET_background_file=NULL,phenotype_names,
                           phenotype_comparison,CPDB_database_file=NULL,Results_index){
  
  
  wkd <- working_dir
  setwd(working_dir)
  
  if(is.null(CPDB_database_file)){CPDB_database_file <- paste0(working_dir,"/Databases/CPDB_pathways_genes.tab")}
  
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
                                          MONET_background_file,phenotype_names,phenotype_comparison)
    )
  }
  
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip'))){file.remove(paste0(Results_dir,'/ResultsZip'))}
  utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir)
  
  gc()
  
  
}



Circos_and_auxiliary <- function(working_dir,phenotype_names,phenotype_comparison,files_edges_path=NULL,
                                 centralities_file=NULL,TF_Database=NULL,file_extension=NULL,Results_index,
                                 Kinome_database=NULL){
  
  wkd <- working_dir
  
  if(is.null(TF_Database)){TF_Database <- paste0(working_dir,"/Databases/_TF.txt")}
  
  setwd(wkd)
  Results_dir <- paste0(getwd(),"/Results_",Results_index)
  
  working_dir <- paste0(Results_dir,"/MONET_analysis/MONET")
  
  if(is.null(files_edges_path)){
    files_edges_path <- paste0(Results_dir,"/edge_files_STRINGBioGRIDIntAct/Symbol")
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
    centralities_file <- paste0(Results_dir,"/STRINGBioGRIDIntAct_centralities_values_CINNA_Total.xlsx")
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
      status <- Cricos_plots(working_dir,files_edges_path,TF_Database,file_extension)
    )
  }
  
  setwd(wkd)
  
  
  pathways_dir <- paste0(working_dir,"/CPDB_Pathways")
  
  source("Monet_cluster_pathways_image_joint.R")
  
  status <- NULL
  attempt <- 1
  while( is.null(status) && attempt <= 5 ) {
    attempt <- attempt + 1
    try(
      status <- MONET_cluster_pathways_joint_image(pathways_dir,file_extension)
    )
  }
  
  library(utils)
  if(file.exists(paste0(Results_dir,'/ResultsZip'))){file.remove(paste0(Results_dir,'/ResultsZip'))}
  utils::zip(zipfile = paste0(Results_dir,'/ResultsZip.zip'), files  = Results_dir)
  
  gc()
  
  setwd(wkd)
  
  if(is.null(Kinome_database)){Kinome_database <- paste0(wkd,"/Databases/kinome.txt")}
  
  tmp_omics_names <- list.files(paste0(Results_dir,'/tmp_FilesIni'))
  tmp_omics_names <- strsplit(tmp_omics_names,'\\.')
  tmp_omics_names <- sapply(tmp_omics_names,"[[",1)
  
  rmarkdown::render(
    input  = '~/Networks_Analyisis_results_PDF_generation.Rmd'
    , params = list(
      results_dir = Results_dir,
      centralities_file_name = centralities_file,
      omics_file_names = tmp_omics_names,
      TF_file = TF_Database,
      kinome_file = Kinome_database,
      pathways_folder = pathways_dir
    ),
    output_file = paste0(Results_dir,"/Network_Analysis_results")
  )
  
}
