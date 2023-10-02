Extract_CentralitiesValues <- function(working_dir,centralities_file,phenotype_names,file_extension){

library(readr)
library(stringr)
library(readxl)
  
setwd(working_dir)

sheets <- excel_sheets(centralities_file)

module_files <- list.files(path = getwd(),pattern = paste0("result-modules__",file_extension))

for (i in 1:length(sheets)){
  
  centralities <- read_excel(centralities_file,sheet = sheets[i])
  centralities$Current_Flow_Betweenness_Centrality <- as.numeric(as.character(centralities$Current_Flow_Betweenness_Centrality))
  
  centralities <- centralities[match(sort(centralities$Current_Flow_Betweenness_Centrality,decreasing = TRUE),centralities$Current_Flow_Betweenness_Centrality),match(c("Current_Flow_Betweenness_Centrality","Symbol"),colnames(centralities))]
  
  m_ind <- grep(sheets[i], module_files)
  
  modules <- read.table(module_files[m_ind],header = FALSE,fill = TRUE,sep="\t")
  
  # is <- apply(modules,function(x) {length(x)-length(which(x %in% ""))},MARGIN=1)
  # is <- which(is>=12)
  # modules <- modules[is,]
  
  centralities_moduels <- data.frame(matrix(nrow=nrow(centralities),ncol=4))
  
  for (j in 1:nrow(centralities)){
    
    if(is.na(centralities[j,2]) == 0){
      aux <-  as.character(which(apply(modules,FUN = function(x) {which(centralities[j,2] %in% x)}, MARGIN=1)==1))
      if(length(aux)>0){
    centralities_moduels[j,1] <- as.character(centralities[j,2])
    centralities_moduels[j,2] <- aux
    centralities_moduels[j,3] <- centralities[j,1]
    centralities_moduels[j,4] <- 1
      }
    }
  }
  
  
  aux <- as.character(module_files[m_ind])
  aux <- sapply(strsplit(aux,paste0(file_extension,"_")),"[",2)
  aux <- sapply(strsplit(aux,"\\."),"[",1)
  dir.create(paste0(getwd(),"/Centralities_modules_links"))
  save_name = paste(getwd(),"/Centralities_modules_links/",aux,"_centralities_modules_links.txt",sep="")
  file.remove(save_name)
  colnames(centralities_moduels) <- c("Gene","SumModule","CentralityScore","Weight")
  write.table(centralities_moduels, file = save_name,row.names = FALSE,col.names = TRUE,quote = FALSE,sep="\t")
  
}
return(1)
}