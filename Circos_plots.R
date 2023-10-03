#The function creates the circular diagrams that characterize a full protein-protein interaction network.

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


Cricos_plots <- function(working_dir,edges_dir,TF_Database,file_extension){

library(readr)
library(RColorBrewer)
library(circlize)
library(stringr)
library(ComplexHeatmap)
library(gridBase)
library(gtools)
library(biomaRt)



mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

setwd(paste(working_dir,"/Centralities_modules_links",sep=""))

centralities_modules_files <- list.files(path = getwd(),pattern = "centralities_modules_links.txt")

edge_modules_files <- list.files(path = edges_dir,pattern = paste0(file_extension,"_"))

tf_list <- read.table(TF_Database,sep="\t",fill=TRUE,header=TRUE)

for (i in 1:length(centralities_modules_files)){
  
  
  data <- read.table(centralities_modules_files[i],header = TRUE,fill = TRUE,sep="\t")
  
  
  data_10pc <- data[1:floor(nrow(data)/3),c(1,2,4)]
  
  data <- data_10pc
  
  aux <- str_split(string = centralities_modules_files[i],pattern = "_")
  aux <- aux[[1]][1]
  
  
  
  ind <- grep(aux,edge_modules_files)
  edge_data <- read.table(paste0(edges_dir,"/",edge_modules_files[ind]),header = FALSE,fill = TRUE,sep="\t")
  
  #Remove the centrality score
   # data <- data[1:floor(nrow(data)/10),c(1,2,4)]
  
   ind <- which(data$SumModule>15)
   if (length(ind)>0){data <- data[-ind,]}

   # ams <- colSums(table(data[,c(1,2)]))
   # ind <- which(ams>7)
   # if (length(ind)>0){data <- data[-ind,]}
  
  #data <- data[,c(1,2,4)]
  
  ind <- which(data$Gene %in% tf_list$Symbol)
  
  if(length(ind)>7){ind<-ind[1:7]}

  data_tf <- c()

  for (j in 1:length(ind)){
    ind1 <- which(edge_data$V1 %in% data$Gene[ind[j]])
    ind2 <- which(edge_data$V2 %in% data$Gene[ind[j]])
    aux <- unique(edge_data[union(ind1,ind2),c(1,2)])

    ind1 <- which(data$Gene %in% aux$V1)
    ind2 <- which(data$Gene %in% aux$V2)
    aux2 <- unique(data$Gene[union(ind1,ind2)])
    aux2 <- aux2[-which(aux2 %in% tf_list$Symbol)]
    
    if(length(aux2)>0){
    data_tf <- rbind(data_tf,data.frame(aux2,data$Gene[ind[j]]))
    }

  }
  
  if(length(unique(data_tf[,1]))>30){
    ind <- which(data_10pc[,1] %in% data_tf[,1])
    if(length(ind)>25) {ind <- ind[1:25]}
    data_tf <- data_tf[which(data_tf[,1] %in% data_10pc[ind,1]),]
  }

  if(is.null(data_tf)){next()}
  
  colnames(data_tf) <- c("Protein","TF")
  data_tf <- cbind(data_tf,data.frame(weight=rep(1,nrow(data_tf)), SubModule = sapply(data_tf$Protein, FUN = function(x) {data[which(data[,1] %in% x),2]} )))
  
  data_tf <- data_tf[order(data_tf$Protein),]
  
  color_palette_members <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
  
  color_palette_Group <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(5, "Set1"))
  
  ##########Uncomment to create no grouping plot
  #grid.col <- c(color_palette_Group(length(unique(data_tf[,2]))), color_palette_members(length(unique(data_tf[,1]))))
  #names(grid.col) <- c(unique(data_tf[,2]), as.character(unique(data_tf[,1])))
  ##########################################
  
  #########Comment when commenting the above
  grid.col <- color_palette_members(length(unique(data_tf[,4])))
  names(grid.col) <- as.character(unique(data_tf[,4]))
  grid.col_aux <- grid.col
  grid.col <- sapply(data_tf[,4], FUN = function(x){grid.col[which(names(grid.col) %in% x)]})
  names(grid.col) <- data_tf[,1]

  aux <- as.character(data[which(data$Gene %in% data_tf$TF),2])
  names(aux) <- data[which(data$Gene %in% data_tf$TF),1]

  grid.col <- c(grid.col_aux[match(aux,unique(names(grid.col_aux)))], grid.col[match(unique(names(grid.col)), names(grid.col))])
  #grid.col <- c(rep("darkorange",length(aux)), grid.col[match(unique(names(grid.col)), names(grid.col))])
  names(grid.col) <- c(names(aux),unique(data_tf[,1]))
  ###############Until here#####################
  data_circos <- data_tf[,c(1,2,3)]
  
  circos_plot1 <- function(){
    # plot.new()
    # circle_size = unit(1.3, "snpc")
    # pushViewport(viewport(x = 0.1, y = 0.7, width = circle_size, height = circle_size,
    #                       just = c("left", "center")))
    # par(omi = gridOMI(), new = TRUE)
    circos.clear()
    circos.par(start.degree = 180)
    circos.initialize(sectors = "a", xlim = c(0,length(unique(data_circos[,1]))) )
    par(cex = 1.5, mar = c(0, 0, 0, 0))
    chordDiagram(data_circos,transparency = 0.2,
                 order = c(data_circos[,1],
                           sort((unique(data_circos[,2])),decreasing = TRUE)),
                 annotationTrack = c("name", "grid"), scale = TRUE,
                 grid.col = grid.col,big.gap = 10,
                 annotationTrackHeight = c(0.03, 0.03))
    #Comment to remove the legend##
    aaux1 <- unique(data_tf[,4])
    aaux1 <- sort(aaux1)
    aaux2 <- grid.col_aux[unique(match(unique(data_tf[,4]),unique(names(grid.col_aux))))]
    aaux2 <- aaux2[mixedsort(names(aaux2))]
    lgd_links = Legend(labels = aaux1, legend_gp = gpar(fill=aaux2),
                       title_position = "topleft", title = "Modules",labels_gp = gpar(fontsize=17),title_gp = gpar(fontsize = 19, fontface = "bold"))
   # upViewport()
    # draw(lgd_links, x = circle_size, just = "left")
    # draw(lgd_links, x = circle_size, just = "left")
    draw(lgd_links, x = unit(1, "npc") - unit(3, "mm"), y = unit(13, "mm"),
         just = c("right", "bottom"))
  }
  
  dir.create(paste0(getwd(),"/Images"))
#  save_name <- paste(getwd(),"/Images/",sapply(strsplit(centralities_modules_files[i],"_"),"[",1),"_chordDiagram.png",sep="")
#  if (file.exists(save_name)) {
#    #Delete file if it exists
#    file.remove(save_name)
#  }
#  png(file = save_name,width = 11.5, height = 9, units = "in", res = 600)
#  circos_plot1()
#  dev.off()

  
  save_name <- paste(getwd(),"/Images/",sapply(strsplit(centralities_modules_files[i],"_"),"[",1),"_chordDiagram.pdf",sep="")
  if (file.exists(save_name)) {
    #Delete file if it exists
    file.remove(save_name)
  }
  pdf(file = save_name,
      width = measurements::conv_unit(x = 320, from = "mm", to = "inch"),
      height = measurements::conv_unit(x = 250, from = "mm", to = "inch"))
  circos_plot1()
  dev.off()

  
  
}
return(1)
}