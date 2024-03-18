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




Kinase_Kinase_networks <- function(working_dir,Cytoscape_dir,Cytoscape_style_file,network_directed){
  
  tryCatch(
    {
      system(Cytoscape_dir,wait = FALSE)
    },
    error = function(cond) {
      stop ("ERROR: Cannot open Cytoscape")
    },
    warning = function(cond) {
      stop ("ERROR: Cannot open Cytoscape")
    }
  )
  
  Sys.sleep(60)
  
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
  
  
  file_name <- paste0(working_dir,"/Kinase_Network_data_phenotype_comaprison_unique_value_no_residues_merged_elements_Prepared_for_visualization.xlsx")
  file_name_auxiliary <- paste0(working_dir,"/Kinase_Network_data_aux_info.xlsx")
  xlsx_sheets <- openxlsx::getSheetNames(file_name)
  
  for (i in xlsx_sheets){
    
  data <- openxlsx::read.xlsx(file_name,sheet = i,colNames = TRUE,rowNames = FALSE)
  auxiliary_data <- openxlsx::read.xlsx(file_name_auxiliary,sheet = i,colNames = TRUE,rowNames = FALSE,)
  
  grh <- graph_from_edgelist(as.matrix(data[,c(ncol(data)-1,ncol(data))]),directed=network_directed)
  
  network <- createNetworkFromIgraph(grh,title = i)
  
  style <- importVisualStyles(filename = Cytoscape_style_file)
  
  setVisualStyle(style)
  
  Sys.sleep(5)
  
  ind <- match(unique(auxiliary_data$Kinase_ID),auxiliary_data$Kinase_ID)
  df <- data.frame (as.integer(auxiliary_data[ind,1]),
                    stringsAsFactors = FALSE)
  rownames(df) <- auxiliary_data[ind,2]
  colnames(df) <- 'Kinase_presence_0_is_not_measured_1_is_measured_and_upreg_2_is_measured'
  loadTableData (df,table = "node",table.key.column = "shared name")
  
  ind <- match(unique(data$Substrate_gene),data$Substrate_gene)
  df <- data.frame (as.integer(data$Substrate_TF[ind]),
                    stringsAsFactors = FALSE)
  rownames(df) <- data$Substrate_gene[ind]
  colnames(df) <- colnames(data)[which(colnames(data) %in% "Substrate_TF")]
  loadTableData (df,table = "node",table.key.column = "shared name")
  
  df <- data.frame (as.integer(data$Database),source=data$Kinase_gene,target=data$Substrate_gene,
                    stringsAsFactors = FALSE)
  
  rownames(df) <- paste0(df$source," (interacts with) ",df$target)
  colnames(df)[1] <- colnames(data)[which(colnames(data) %in% "Database")]
  df <- subset(df,select = 1)
  loadTableData (df,table = "edge",table.key.column = "shared name")
  
  layoutNetwork("force-directed")
  
  }

  saveSession(filename = paste0(working_dir,"/Kinase-Kinase_Networks"))
  Sys.sleep(20)
  
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