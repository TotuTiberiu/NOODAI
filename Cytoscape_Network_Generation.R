#The function constructs automatically the protein-protein interaction network. It requires a graphical interface!

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


#The centrality file associated with the respective edge file
Centrality_file <- "C:/work/Files/Conferences/PH_24.08.2024/Live_Demo/Demo/Additional_Results/PPINetworks_centralities_values_CINNA_Total.xlsx"
#The Excel sheet from which to read the centrality values
Centrality_file_sheet <- "M1vsM2a"
#Provide the names of the omics datasets
names_of_omicsDatasets <- c("Omics1","Omics2","Omics3","Omics4")
#The Edge file that contain the network as node-node list of values. It must have the same identifier as the Centrality file!
Edge_file <- "C:/work/Files/Conferences/PH_24.08.2024/Live_Demo/Demo/Additional_Results/edge_files_PPINetworks/Symbol/Total_M1vsM2a.txt"
#Edge_file <- "C:/work/Files/Conferences/PH_24.08.2024/Live_Demo/Demo/Additional_Results/MONET_analysis/MONET/Edges_tables/2024-06-22-225802__M1__result-modules__Total_M1vsM2a_edges.xlsx"
#If the Edge file is an Excel table, provide the appropriate sheet
Edge_file_sheet <- "Cluster_1"
#Working directory where to save the files
working_dir <- "C:/work/Files/Conferences/PH_24.08.2024/Live_Demo/Cytoscape_Networks"
#The path to the Cytoscape executable
Cytoscape_path <- "C:/Programs/Cytoscape_v3.10.1/Cytoscape.exe"
#The Cytoscape style file that contains the color mappings and nodes dimensions mapping
Cytoscape_style_file <- "C:/work/Files/Conferences/PH_24.08.2024/Live_Demo/Style_NOODAI_Networks.xml"
#The name of the created Cytoscape Network
save_name_netowrk <- "Cytoscape_Network"

Cytoscape_networks <- function(Centrality_file,Centrality_file_sheet,names_of_omicsDatasets,Edge_file,
                               Edge_file_sheet=NULL,working_dir,Cytoscape_path,Cytoscape_style_file,save_name_netowrk="Cytoscape_Network"){
  
  #Open Cytoscape 
  tryCatch(
    {
      system(Cytoscape_path,wait = FALSE)
    },
    error = function(cond) {
      stop ("ERROR: Cannot open Cytoscape")
    },
    warning = function(cond) {
      stop ("ERROR: Cannot open Cytoscape")
    }
  )
  
  Sys.sleep(60)
  
  #Library
  if(!"RCy3" %in% installed.packages()){
    install.packages("BiocManager")
    BiocManager::install("RCy3")
  }
  library(RCy3)
  
  if(!"igraph" %in% installed.packages()){
    install.packages("BiocManager")
    BiocManager::install("igraph")
  }
  library(igraph)
  
  #Identify if it is a txt file (full-size netowrk) or a module that is stored in Excel. Read the edges files
  sep_ident <- strsplit(Edge_file,split = "\\.")
  sep_ident <- sep_ident[[1]][lengths(sep_ident)]
  if(sep_ident=="txt"){
    data <- read.table(Edge_file,header = FALSE,sep = "\t",quote = "",fill = TRUE)
  }
  if(sep_ident=="xlsx"){
    data <- openxlsx::read.xlsx(Edge_file,sheet = Edge_file_sheet,colNames = FALSE,rowNames = FALSE)
  }
  
  #Read the centrality data and the omics presence info
  auxiliary_data <- openxlsx::read.xlsx(Centrality_file,sheet = Centrality_file_sheet,colNames = TRUE,rowNames = FALSE)
  
  #Select only the centrality column and omics presence
  selected_col <- c(which(colnames(auxiliary_data) %in% "Symbol"),
                    which(colnames(auxiliary_data) %in% "Current_Flow_Betweenness_Centrality"),
                    which(colnames(auxiliary_data) %in% names_of_omicsDatasets))
  auxiliary_data <- auxiliary_data[,selected_col]
  auxiliary_data$Omics_Total <- apply(auxiliary_data[,which(colnames(auxiliary_data) %in% names_of_omicsDatasets)],MARGIN=1,FUN=function(x){paste(x,collapse = "")})
  
  #Manage if there are more or less than 4 omics datasets
  nr_omics <- length(which(colnames(auxiliary_data) %in% names_of_omicsDatasets))
  
  if(nr_omics<4){
    nr <- 4-length(which(colnames(auxiliary_data) %in% names_of_omicsDatasets))
    auxiliary_data$Omics_Total <- apply(cbind(matrix(0,nrow=nrow(auxiliary_data),ncol=nr),auxiliary_data$Omics_Total),MARGIN=1,FUN=function(x){paste(x,collapse = "")})
  }
  if(nr_omics>4){
    auxiliary_data$Omics_Total <- substring(auxiliary_data$Omics_Total, nchar(auxiliary_data$Omics_Total) - 3, nchar(auxiliary_data$Omics_Total))
  }
  
  #Create the graph (undirected)
  grh <- graph_from_edgelist(as.matrix(data[,c(1,2)]),directed=FALSE)
  
  if(sep_ident=="txt"){
    title_graph <- Centrality_file_sheet
  }else{
    title_graph <- Edge_file_sheet
  }
  
  network <- createNetworkFromIgraph(grh,title = title_graph)
  #Import the preset Cytoscape style
  style <- importVisualStyles(filename = Cytoscape_style_file)
  
  setVisualStyle(style)
  
  Sys.sleep(5)
  
  #Tailor the Centrality and Omics presence info to match the preset style. This means to change the column names to Omics 1-4 and to set the rownames to Symbol values
  ind <- match(unique(auxiliary_data$Symbol),auxiliary_data$Symbol)
  rownames(auxiliary_data) <- auxiliary_data[ind,1]
  auxiliary_data <- auxiliary_data[,-1]
  new_col_names <- c("Omics1","Omics2","Omics3","Omics4")
  if(nr_omics>4){nr_omics<-4}
  colnames(auxiliary_data)[2:(ncol(auxiliary_data)-1)] <- new_col_names[c(1:nr_omics)]
  loadTableData (auxiliary_data,table = "node",table.key.column = "shared name")
  
  #Set a layout
  layoutNetwork("force-directed")
  
  #Save the file
  saveSession(filename = paste0(working_dir,"/",save_name_netowrk))
  Sys.sleep(20)
  
  #Close Cytoscape
  tryCatch(
    {
      system("TASKKILL /F /IM Cytosc* /T",wait = FALSE)
    },
    error = function(cond) {
      stop ("ERROR: Cannot close Cytoscape")
    },
    warning = function(cond) {
      stop ("ERROR: Cannot close Cytoscape")
    }
  )
  
  
}