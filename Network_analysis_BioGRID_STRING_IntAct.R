#The function constructs a protein-protein interaction network from the uploaded omics data for each comparison of interest. In addition it computes the main centrality scores for each node, saves the edge files.
# 
#     Copyright © 2023, Empa, Tiberiu Totu.
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







Network_analysis <- function(working_dir,Results_dir,BioGRID_data_file=NULL,STRING_data_file=NULL,IntAct_data_file=NULL,file_DEA_names,phenotype_names,phenotype_comparison,splicing_file_name,Use_precompiled_database,LookUp_table_file=NULL){
  
  splicing_file_index <- grep(splicing_file_name,file_DEA_names)
  
  if(length(splicing_file_index)==0){splicing_file_index <- -1}
  if(length(splicing_file_index)>1){stop("Please provide the name of the DTU file. It should be unique.")}
  
  library(readxl)
  #library(msigdbr)
  #library(AnnotationDbi)
  library(biomaRt)
  #library(org.Hs.eg.db)
  library(centiserve)
  library(CINNA)
  library(igraph)
  library(tidygraph)
  library(stringr)
  library(parallel)
  
  n_cores <- detectCores()
  n_cores_used <- floor(n_cores/2)
  #n_cores_used <- 3
  cluster_1 <- makeCluster(n_cores_used,type = "PSOCK")
  clusterEvalQ(cluster_1, library(centiserve))
  clusterEvalQ(cluster_1,library(CINNA))
  clusterEvalQ(cluster_1,library(igraph))
  clusterEvalQ(cluster_1,library(tidygraph))
  
  
  
  options(timeout = max(1000, getOption("timeout")))
  options(java.parameters = "-Xmx6000m")
  set.seed(21)
  
  setwd(working_dir)
  dir.create(Results_dir)
  edge_file_folder <- paste0(Results_dir,"/edge_files_STRINGBioGRIDIntAct")
  dir.create(edge_file_folder)
  dir.create(paste0(edge_file_folder,'/Uniprot'))
  dir.create(paste0(edge_file_folder,'/Symbol'))
  dir.create(paste0(edge_file_folder,'/PCA_validation'))
  
  #Functions######################################################################
  
  get_symbol_idMART <- function(names,mart){
    MART_sym1 <- getBM(filters = "entrezgene_id",
                       attributes = c("entrezgene_id", "external_gene_name"),
                       values = as.character(names[,1]), mart = mart,useCache = FALSE)
    edge1 <- names[,1]
    edge1 <- MART_sym1[match(edge1, MART_sym1$entrezgene_id),2]
    edge2 <- names[,2]
    MART_sym2 <- getBM(filters = "entrezgene_id",
                       attributes = c("entrezgene_id", "external_gene_name"),
                       values = as.character(names[,2]), mart = mart,useCache = FALSE)
    edge2 <- MART_sym2[match(edge2, MART_sym2$entrezgene_id),2]
    edge <- cbind(edge1,edge2)
    ind <- which(is.na(edge[,1]))
    if(length(ind)>0) {edge <- edge[-ind,]}
    ind <- which(is.na(edge[,2]))
    if(length(ind)>0) {edge <- edge[-ind,]}
    ind <- which(edge[,1]=="")
    if(length(ind)>0) {edge <- edge[-ind,]}
    ind <- which(edge[,2]=="")
    if(length(ind)>0) {edge <- edge[-ind,]}
    return(edge)
  }
  
  get_entrez_idMART <- function(names,mart){
    MART_sym1 <- getBM(filters = "ensembl_peptide_id",
                       attributes = c("ensembl_peptide_id", "entrezgene_id"),
                       values = as.character(names[,1]), mart = mart,useCache = FALSE)
    protein1 <- names[,1]
    protein1 <- MART_sym1[match(protein1, MART_sym1$ensembl_peptide_id),2]
    protein2 <- names[,2]
    MART_sym2 <- getBM(filters = "ensembl_peptide_id",
                       attributes = c("ensembl_peptide_id", "entrezgene_id"),
                       values = as.character(names[,2]), mart = mart,useCache = FALSE)
    protein2 <- MART_sym2[match(protein2, MART_sym2$ensembl_peptide_id),2]
    protein <- cbind(protein1,protein2)
    ind <- which(is.na(protein[,1]))
    if(length(ind)>0) {protein <- protein[-ind,]}
    ind <- which(is.na(protein[,2]))
    if(length(ind)>0) {protein <- protein[-ind,]}
    ind <- which(protein[,1]=="")
    if(length(ind)>0) {protein <- protein[-ind,]}
    ind <- which(protein[,2]=="")
    if(length(ind)>0) {protein <- protein[-ind,]}
    return(protein)
  }
  
  get_entrez_idMART_from_uniprot <- function(names,mart){
    MART_sym1 <- getBM(filters = "uniprot_gn_id",
                       attributes = c("uniprot_gn_id", "entrezgene_id"),
                       values = as.character(names[,1]), mart = mart,useCache = FALSE)
    protein1 <- names[,1]
    protein1 <- MART_sym1[match(protein1, MART_sym1$uniprot_gn_id),2]
    protein2 <- names[,2]
    MART_sym2 <- getBM(filters = "uniprot_gn_id",
                       attributes = c("uniprot_gn_id", "entrezgene_id"),
                       values = as.character(names[,2]), mart = mart,useCache = FALSE)
    protein2 <- MART_sym2[match(protein2, MART_sym2$uniprot_gn_id),2]
    protein <- cbind(protein1,protein2)
    ind <- which(is.na(protein[,1]))
    if(length(ind)>0) {protein <- protein[-ind,]}
    ind <- which(is.na(protein[,2]))
    if(length(ind)>0) {protein <- protein[-ind,]}
    ind <- which(protein[,1]=="")
    if(length(ind)>0) {protein <- protein[-ind,]}
    ind <- which(protein[,2]=="")
    if(length(ind)>0) {protein <- protein[-ind,]}
    return(protein)
  }
  
  get_uniprot_idMART <- function(names,mart){
    MART_sym1 <- getBM(filters = "entrezgene_id",
                       attributes = c("entrezgene_id", "uniprot_gn_id"),
                       values = as.character(names[,1]), mart = mart,useCache = FALSE)
    edge1 <- names[,1]
    edge1 <- MART_sym1[match(edge1, MART_sym1$entrezgene_id),2]
    edge2 <- names[,2]
    MART_sym2 <- getBM(filters = "entrezgene_id",
                       attributes = c("entrezgene_id", "uniprot_gn_id"),
                       values = as.character(names[,2]), mart = mart,useCache = FALSE)
    edge2 <- MART_sym2[match(edge2, MART_sym2$entrezgene_id),2]
    edge <- cbind(edge1,edge2)
    ind <- which(is.na(edge[,1]))
    if(length(ind)>0) {edge <- edge[-ind,]}
    ind <- which(is.na(edge[,2]))
    if(length(ind)>0) {edge <- edge[-ind,]}
    ind <- which(edge[,1]=="")
    if(length(ind)>0) {edge <- edge[-ind,]}
    ind <- which(edge[,2]=="")
    if(length(ind)>0) {edge <- edge[-ind,]}
    return(edge)
  }
  
  get_protein_idMART <- function(names,mart){
    if(length(names)>0){
      aux2 <- getBM(filters = "uniprot_gn_id",
                    attributes = c("entrezgene_id", "uniprot_gn_id"),
                    values = as.character(names), mart = mart,useCache = FALSE)
      return(aux2$entrezgene_id)
    }else {return(names)}
  }
  
  get_protein_idMART1 <- function(names,mart){
    if(length(names)>0){
      aux2 <- getBM(filters = "external_gene_name",
                    attributes = c("entrezgene_id", "external_gene_name"),
                    values = as.character(names), mart = mart,useCache = FALSE)
      return(aux2$entrezgene_id)
    }else {return(names)}
  }
  
  get_dataframefromlist <- function(data,phenotype_comparison){
    increment <- max(lengths(data))
    aux <- lapply(data, `length<-`, increment)
    aux <- do.call(data.frame,aux)
    colnames(aux) <- phenotype_comparison
    return(aux)
  }
  
  extract_edges_NONDIRECTIONAL <- function(data_ini,Lookup_data){
    ind1 <- which(Lookup_data$protein1 %in% data_ini)
    ind2 <- which(Lookup_data$protein2 %in% data_ini)
    ind <- intersect(ind1,ind2)
    aux <- Lookup_data[ind,c(1,2)]
    ###THIS PART DESTROYES AND DIRECTIONALITY######################################################## Remove it to make a directional edge table
    for (j in 1:nrow(aux))
    {
      aux[j, ] = str_sort(aux[j,],numeric = TRUE)
    }
    aux = aux[!duplicated(aux),]
    return(aux)
  }
  
  centrality_extraxtion <- function (aux,save_name,i,mart,phenotype_comparison,Transcriptomics_DTU_total,flag,multi_omics_data,enhancement_flag,file_DEA_names_short){
    
    aux_bck <- aux
    if ((length(which(is.na(aux)))==length(aux)) == FALSE){
      
      isd <- which(is.na(aux[,1]))
      if (length(isd)>0){aux <- aux[-isd,]}
      isd <- which(is.na(aux[,2]))
      if (length(isd)>0){aux <- aux[-isd,]}
      aux <- graph_from_edgelist(aux,directed=FALSE)
      aux_decomposed <- decompose.graph(aux)
      if(max(lengths(aux_decomposed))>10){
        
        length_aux <- lengths(aux_decomposed)
        ind <- which(length_aux == max(length_aux))
        #If the flag is on, then we analyse combined data, hence we do not want networks build from more than 75% DTU data
        if(flag==1){
          while(length(which(as_edgelist(aux_decomposed[[ind]])[,1] %in% Transcriptomics_DTU_total[,i]))/nrow(as_edgelist(aux_decomposed[[ind]]))>0.75){
            aux_decomposed <- aux_decomposed[-ind]
            length_aux <- lengths(aux_decomposed)
            ind <- which(length_aux == max(length_aux))
          }
        }
        
        aux_decomposed <- aux_decomposed[[ind]]
        #isolated = which(degree(aux)==2)
        #aux <- delete.vertices(aux,isolated)
        
        #Write edge tables
        aux_decomposed_edges <- as_edgelist(aux_decomposed)
        aux_uniprot <- get_uniprot_idMART(aux_decomposed_edges,mart)
        aux_uniprot <- cbind(aux_uniprot,rep(1,length(aux_uniprot[,1])))
        save_name1 <- paste(edge_file_folder,'/Uniprot/',save_name,'_',phenotype_comparison[i],'.txt',sep='')
        file.remove(save_name1)
        write.table(aux_uniprot,file=save_name1,row.names=FALSE,col.names=FALSE,sep = '\t',quote = FALSE)
        aux_symbol <- get_symbol_idMART(aux_decomposed_edges,mart)
        aux_symbol <- cbind(aux_symbol,rep(1,length(aux_symbol[,1])))
        save_name1 <- paste(edge_file_folder,'/Symbol/',save_name,'_',phenotype_comparison[i],'.txt',sep='')
        file.remove(save_name1)
        write.table(aux_symbol,file=save_name1,row.names=FALSE,col.names=FALSE,sep = '\t', quote = FALSE)
        #######
        
        cent_m <- proper_centralities(aux_decomposed)
        #Remove very lengthy calculus centralities - do them afterwards with all the cores if necessary
        cent_m <- cent_m[-c(5,6,31,32,33,38,43)]
        ##########################################################################
        nx_per_worker <- floor(length(cent_m)/n_cores_used)
        length_cent_m <- length(cent_m)
        clusterExport(cluster_1, "aux_decomposed",envir=environment())
        clusterExport(cluster_1, "nx_per_worker",envir=environment())
        clusterExport(cluster_1, "length_cent_m",envir=environment())
        clusterExport(cluster_1, "cent_m",envir=environment())
        fsx <- function(k) {
          if ((k+nx_per_worker-1)<=length_cent_m){
            amx <- calculate_centralities(aux_decomposed, include = cent_m[k:(k+nx_per_worker-1)])
          }else{amx <- calculate_centralities(aux_decomposed, include = cent_m[k:length_cent_m])}
          if(length(which(lengths(amx)==0)) >0 ){amx <- amx[-which(lengths(amx)==0)]}
          return(amx)
        }
        steps <- seq(1,length(cent_m),nx_per_worker)
        centralities_values_CINNA <- parSapply(cluster_1, steps, fsx)
        #centralities_values_CINNA <- c(centralities_values_CINNA)
        aax <- list()
        for (ik in 1:length(centralities_values_CINNA)){
          aax <- c(aax,centralities_values_CINNA[[ik]])
        }
        centralities_values_CINNA <- aax
        
        
        if(length(which(lengths(centralities_values_CINNA)==0)) >0 ){
          centralities_values_CINNA1 <- centralities_values_CINNA[-which(lengths(centralities_values_CINNA)==0)]
        }else {centralities_values_CINNA1 <- centralities_values_CINNA}
        
        #Plot the metrics
        a1 <- pca_centralities(centralities_values_CINNA1)
        
        save_name1 <- paste(edge_file_folder,'/PCA_validation/',save_name,'_',phenotype_comparison[i],'.png',sep='')
        save_name_pdf <- paste(edge_file_folder,'/PCA_validation/',save_name,'_',phenotype_comparison[i],'.pdf',sep='')
        
        png(file = save_name1,width = 55, height = 15, units = "cm", res = 800)
        print(a1, newpage = FALSE)
        dev.off()
        
        
        pdf(file = save_name_pdf, width = 50/2.5, height = 16/2.5)
        print(a1, newpage = FALSE)
        dev.off()
        
        
        centralities_values_CINNA1 <- do.call(data.frame,centralities_values_CINNA1)
        #colnames(centralities_values_CINNA1) <- cent_m
        
        specific_metric <- tbl_graph(edges = as.data.frame(aux_decomposed_edges), directed = FALSE)
        specific_metric <- specific_metric %>% activate(nodes) %>% mutate(importance_current = centrality_betweenness_current()) %>% mutate(importance_katz = centrality_katz()) #%>% mutate(importance_random = centrality_random_walk())
        specific_metric <- specific_metric %>% activate(nodes) %>% data.frame()
        specific_metric <- specific_metric[order(match(rownames(centralities_values_CINNA1), specific_metric$name)),]
        
        centralities_values_CINNA1$Current.Flow.Betweenness.Centrality <- specific_metric$importance_current
        centralities_values_CINNA1$Katz.Betweenness <- specific_metric$importance_katz
        centralities_values_CINNA1$Random.Walk.Betweenness <- specific_metric$importance_random
        
        
        axx <- getBM(filters = "entrezgene_id",
                     attributes = c("entrezgene_id", "external_gene_name"),
                     values = as.character(rownames(centralities_values_CINNA1)), mart = mart,useCache = FALSE)
        names_sym <- axx[match(rownames(centralities_values_CINNA1), axx$entrezgene_id),2]
        centralities_values_CINNA1$Symbol <- names_sym
        colnames(centralities_values_CINNA1) <- str_replace_all(colnames(centralities_values_CINNA1),pattern = "\\.", replacement = "_")
		
		if (enhancement_flag){
		vv1 <- data.frame(matrix(data = 0,nrow=nrow(centralities_values_CINNA1),ncol=length(multi_omics_data)))
		rownames(vv1) <- as.character(rownames(centralities_values_CINNA1))
		colnames(vv1) <- file_DEA_names_short
		for (mxk in 1:length(multi_omics_data)){
		arc <- unique(rbind(multi_omics_data[[mxk]][[i]][,1],multi_omics_data[[mxk]][[i]][,2]))
		vv1[which(rownames(vv1) %in% arc),mxk] <- 1
		}
		centralities_values_CINNA1 <- cbind(vv1,centralities_values_CINNA1)
		save_name1 <- paste(Results_dir,'/Background_total.xlsx',sep='')
		xlsx::write.xlsx(unique(as.character(rbind(aux_bck[,1],aux_bck[,2]))),file=save_name1,sheetName = phenotype_comparison[i],col.names=FALSE,row.names=FALSE,append = TRUE)
		
		}
        
        centralities_values_CINNA1 <- centralities_values_CINNA1[order(centralities_values_CINNA1$Current_Flow_Betweenness_Centrality,decreasing=TRUE),]
        ncaux <- ncol(centralities_values_CINNA1)
        centralities_values_CINNA1 <- centralities_values_CINNA1[,c(ncaux,ncaux-2,ncaux-1,1:(ncaux-3))]
        centralities_values_CINNA1 <- cbind(rownames(centralities_values_CINNA1),centralities_values_CINNA1)
        colnames(centralities_values_CINNA1)[1] <- "Entrez ID"
        save_name1 <- paste(Results_dir,'/STRINGBioGRIDIntAct_centralities_values_CINNA_',save_name,'.xlsx',sep='')
        xlsx::write.xlsx(centralities_values_CINNA1,file=save_name1,sheetName = phenotype_comparison[i],col.names=TRUE,row.names=FALSE,append = TRUE)
      }
    }
  }
  
  add_splicing_inforamtion <- function(i,symmetry_value,data_DTU,file_DEA_edges,data_aux_total,splicing_file_index){
    
    if(i<symmetry_value+1){factor<-i}else{factor<- i-symmetry_value}
    
    aux_total_aux1 <- c()
    for (k in 1:length(file_DEA_names)){
      if(k != splicing_file_index){
        aux_total_aux1 <- as.matrix(rbind(aux_total_aux1,as.matrix(file_DEA_edges[[k]][[factor]])))
      }
    }
    
    msx <- file_DEA_edges[[splicing_file_index]][[i]]
    isx1 <- which(msx[,1] %in% aux_total_aux1)
    isx2 <- which(msx[,2] %in% aux_total_aux1)
    isx_v1 <- union(isx1,isx2)
    
    aux_total_aux2 <- c()
    for (k in 1:length(file_DEA_names)){
      if(k != splicing_file_index){
        aux_total_aux2 <- as.matrix(rbind(aux_total_aux2,as.matrix(file_DEA_edges[[k]][[factor+symmetry_value]])))
      }
    }
    
    isx1 <- which(msx[,1] %in% aux_total_aux2)
    isx2 <- which(msx[,2] %in% aux_total_aux2)
    isx_v2 <- union(isx1,isx2)
    isx <- union(isx_v1,isx_v2)
    if(length(isx)>0) {msx <- msx[-isx,]}
    data_aux_total <- as.matrix(rbind(data_aux_total,msx))
    return(data_aux_total)
  }
  
  #END Functions part#############################################################
  
  #Read the databases files#######################################################
  # 
  # cent_m <- c("Closeness.centrality..Latora.","Dangalchev.Closeness.Centrality","Decay.Centrality","Residual.Closeness.Centrality",
  #             "Semi.Local.Centrality","Barycenter.Centrality","Diffusion.Degree","eigenvector.centralities","Markov.Centrality",
  #             "Wiener.Index.Centrality","Average.Distance","Radiality.Centrality","Geodesic.K.Path.Centrality","Current.Flow.Closeness.Centrality",
  #             "Degree.Centrality","Leverage.Centrality","Centroid.value","Laplacian.Centrality","Topological.Coefficient","ClusterRank",
  #             "subgraph.centrality.scores","Eccentricity.Centrality","Information.Centrality","Flow.Betweenness.Centrality","Shortest.Paths.Betweenness.Centrality","Entropy.Centrality")
  
  # cent_m <- c("Closeness centrality (Latora)","eigenvector centralities","Markov Centrality","Shortest-Paths Betweenness Centrality",
  #             "Degree Centrality","Entropy Centrality","Laplacian Centrality","Topological Coefficient","Information Centrality",
  #             "subgraph centrality scores")
  
  mart <- useMart('ENSEMBL_MART_ENSEMBL')
  mart <- useDataset('hsapiens_gene_ensembl', mart)
  
  if(Use_precompiled_database==0){
  
  BioGRID_data <- read.table(BioGRID_data_file,sep="\t",header = TRUE,fill=TRUE)
  
  ind <- which(BioGRID_data$Taxid.Interactor.A %in% "taxid:9606")
  BioGRID_data <- BioGRID_data[ind,]
  
  ind <- which(BioGRID_data$Taxid.Interactor.B %in% "taxid:9606")
  BioGRID_data <- BioGRID_data[ind,]
  
  ind <- which(BioGRID_data$Confidence.Values %in% "-")
  BioGRID_data <- BioGRID_data[-ind,]
  
  aux <- BioGRID_data$Confidence.Values
  aux <- strsplit(aux,split = "[:]")
  aux <- sapply( aux, "[", 2 )
  BioGRID_data$Confidence.Values <- as.numeric(aux)
  
  
  BioGRID_data <- BioGRID_data[,c(1,2)]
  colnames(BioGRID_data) <- c("protein1","protein2")
  aux <- BioGRID_data$protein1
  aux <- strsplit(aux,split = "[:]")
  aux <- sapply( aux, "[", 2 )
  BioGRID_data$protein1 <- aux
  
  aux <- BioGRID_data$protein2
  aux <- strsplit(aux,split = "[:]")
  aux <- sapply( aux, "[", 2 )
  BioGRID_data$protein2 <- aux
  ################################################################################
  
  STRING_data <- read.table(STRING_data_file,sep=" ",header = TRUE)
  
  STRING_data <- STRING_data[STRING_data$combined_score>700,]
  
  aux <- STRING_data$protein1
  aux <- strsplit(aux,split = "[.]")
  aux <- sapply( aux, "[", 2 )
  STRING_data$protein1 <- aux
  
  aux <- STRING_data$protein2
  aux <- strsplit(aux,split = "[.]")
  aux <- sapply( aux, "[", 2 )
  STRING_data$protein2 <- aux
  
  STRING_data <- STRING_data[,c(1,2)]
  
  STRING_data <- as.data.frame(get_entrez_idMART(STRING_data,mart))
  ################################################################################
  IntAct_data <- read.table(IntAct_data_file,sep=" ",header = FALSE)
  
  IntAct_data <- as.data.frame(get_entrez_idMART_from_uniprot(IntAct_data,mart))
  #Create the joint database file between STRING, BioGRID and IntAct
  Lookup_data <- rbind(BioGRID_data,STRING_data,IntAct_data)
  }
  if(Use_precompiled_database==1){
    Lookup_data <- read.table(LookUp_table_file,header=TRUE,colClasses = "character")
  }
  ################################################################################
  
  #Read the protein files names###################################################
  

  
  #Define the analyzed phenotypes
  phenotype_comparison_aux = t(as.data.frame(strsplit(phenotype_comparison, split="vs")))
  rownames(phenotype_comparison_aux) <- phenotype_comparison
  
  symmetry_value <- length(phenotype_comparison)/2 #Half of the number of comparison groups when these are symmetric
  
  #Read the measured data#########################################################
  
  file_DEA <- list(list())
  for (i in 1:length(file_DEA_names)){
    file_DEA[[i]] <- list()
  }
  
  for (j in 1:length(file_DEA_names)){
    for (i in 1:length(phenotype_comparison)){
      
      aux <- as.data.frame(read_excel(file_DEA_names[j],sheet = phenotype_comparison[i]))
      #Get the corresponding entrez ids from the uniprot ones#######################
      if(length(aux)>0){
        file_DEA[[j]][[i]] <- as.character(unique(get_protein_idMART(aux$Protein,mart)))
      }else{
        file_DEA[[j]][[i]] <- NA
      }
      
    }
  }
  ################################################################################
  
  #Construct data frames from the lists###########################################
  
  file_DEA_total <- list()
  for (j in 1:length(file_DEA_names)){
    file_DEA_total[[j]] <- get_dataframefromlist(file_DEA[[j]],phenotype_comparison)
  }
  
  #Edge table creation############################################################
  
  file_DEA_edges <- list(list())
  for (i in 1:length(file_DEA_names)){
    file_DEA_edges[[i]] <- list()
  }
  
  for (j in 1:length(file_DEA_names)){
    for (i in 1:length(phenotype_comparison)){
      file_DEA_edges[[j]][[i]] <- extract_edges_NONDIRECTIONAL(file_DEA_total[[j]][,i],Lookup_data)
    }
  }
  ################################################################################
  
  #Centrality extraction#########################################################
  
  
  file_DEA_names_short <- file_DEA_names
  file_DEA_names_short <- str_split(file_DEA_names,"\\.")
  file_DEA_names_short <- sapply(file_DEA_names_short,"[[",1)
  if(grep('/',file_DEA_names_short[1])==1){
  file_DEA_names_short <- str_split(file_DEA_names_short,"/")
  file_DEA_names_short <- sapply(file_DEA_names_short,"[[",min(lengths(file_DEA_names_short)))
  }
  
  
  for (j in 1:length(file_DEA_names)){
    file.remove(paste(Results_dir,'/STRINGBioGRIDIntAct_centralities_values_CINNA_',file_DEA_names_short[j],'.xlsx',sep=''))
  }
  
    file.remove(paste(Results_dir,'/Background_total.xlsx',sep=''))
    file.remove(paste(Results_dir,'/STRINGBioGRIDIntAct_centralities_values_CINNA_Total.xlsx',sep=''))
  
  for (i in 1:length(phenotype_comparison)){
    
    for (j in 1:length(file_DEA_names)){
      aux <- as.matrix(file_DEA_edges[[j]][[i]])
      colnames(aux) <- c()
      if(splicing_file_index>0){
        centrality_extraxtion(aux,file_DEA_names_short[j],i,mart,phenotype_comparison,file_DEA_total[[splicing_file_index]],0,c(),0,file_DEA_names_short)
      }else{
        centrality_extraxtion(aux,file_DEA_names_short[j],i,mart,phenotype_comparison,c(),0,c(),0,file_DEA_names_short)
        
      }
    }
    
    aux_total <- c()
    for (j in 1:length(file_DEA_names)){
      if(j != splicing_file_index){
        aux_total <- as.matrix(rbind(aux_total,as.matrix(file_DEA_edges[[j]][[i]])))
      }
    }
    
    if(splicing_file_index >0){
      aux_total <- add_splicing_inforamtion(i,symmetry_value,file_DEA_edges[[splicing_file_index]],
                                            file_DEA_edges,aux_total,splicing_file_index)
    }
    colnames(aux_total) <- c()
    ###THIS PART DESTROYES AND DIRECTIONALITY######################################################## Remove it to make a directional edge table - remove it from the initial function as well
    for (jik in 1:nrow(aux_total))
    {
      aux_total[jik, ] = str_sort(aux_total[jik,],numeric = TRUE)
    }
    aux_total = aux_total[!duplicated(aux_total),]
    if(splicing_file_index>0){
      centrality_extraxtion(aux_total,"Total",i,mart,phenotype_comparison,file_DEA_total[[splicing_file_index]],1,file_DEA_edges,1,file_DEA_names_short)
    }else{
      centrality_extraxtion(aux_total,"Total",i,mart,phenotype_comparison,c(),1,file_DEA_edges,1,file_DEA_names_short)
    }
  }
  
  stopCluster(cluster_1)
  
  return(1)
  
  
}
