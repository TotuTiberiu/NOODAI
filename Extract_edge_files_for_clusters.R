#The function extracts the edges that are found in each identified subnetwork after the MONET decomposition. These are used for the circular diagrams.

# 
#     Copyright © 2025, Empa, Tiberiu Totu.
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


Extract_edges_tables <- function(working_dir,phenotype_names,phenotype_comparison,files_edges_path,file_extension){

library(readxl)
library(stringr)

setwd(working_dir)

files_monet <- list.files(path = getwd(), pattern = paste0("result-modules__",file_extension))
print(files_monet)
dir.create(paste(getwd(),"/Edges_tables",sep=""))
results_dir <- paste(getwd(),"/Edges_tables",sep="")

phenotype_comparison_aux = t(as.data.frame(strsplit(phenotype_comparison, split="vs")))
rownames(phenotype_comparison_aux) <- phenotype_comparison


files_edges <- list.files(path = files_edges_path,pattern = paste0(file_extension,'_'),full.names = TRUE)
print(files_edges)
edges <- list()
for (i in 1:length(phenotype_comparison)){
  ind <- grep(phenotype_comparison[i], files_edges)
  edges[[i]] <- read.table(files_edges[ind])
}
print(length(edges))

for (j in 1:length(phenotype_comparison)){
  
  ind_f <- grep(phenotype_comparison[j], files_monet)
  print(ind_f)
print(getwd())
  data <- read.table(files_monet[ind_f],sep="\t",header = FALSE, fill = TRUE)
  name_index <- str_split(files_monet[ind_f],patter=paste0(file_extension,"_"))
  name_index <- str_split(name_index[[1]][2],patter="[.]")
  name_index <- name_index[[1]][1]
  
  data <- data[,c(3:ncol(data))]
  
  ind <- apply(data,1,function(x) length(which(x=="")))
  keep <- which(ind<=(ncol(data)-10))
  data <- data[keep,]
  
  ss <- str_split(files_monet[ind_f],patter="[.]")
  save_name_edges = paste(results_dir,"/",ss[[1]][1],"_edges.xlsx",sep="")
  if (file.exists(save_name_edges)){file.remove(save_name_edges)}
  
  for (k in 1:nrow(data)){
    
    aux <- data[k,]
    ind <- which(aux=="")
    if (length(ind)>0){ aux <- aux[-ind] }
    
    ind1 <- which(edges[[j]]$V1 %in% aux)
    ind2 <- which(edges[[j]]$V2 %in% aux)
    ind <- intersect(ind1,ind2)
    
    aux_edges <- edges[[j]][ind,]
    
    xlsx::write.xlsx(aux_edges,save_name_edges,sheetName=paste("Cluster_",as.character(k),sep=""),append=TRUE,row.names = FALSE,col.names = FALSE)
    
  }
  
}
return(1)
}
