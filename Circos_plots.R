#The function creates the circular diagrams that characterize a full protein-protein interaction network.

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
#     BioGRID_data_file = "ANY",
#     STRING_data_file = "ANY",
#     IntAct_data_file = "ANY",
#     file_DEA_names = "character",
#     phenotype_names = "character",
#     phenotype_comparison = "character",
#     Use_precompiled_database = "numeric",
#     LookUp_table_file = "ANY",
#     Results_index = "character",
#     edge_file_path = "ANY",
#     monet_path = "ANY",
#     Monet_method_string = "character",
#     tmp_bin_folder = "ANY",
#     CPDB_databases = "ANY",
#     MONET_background_file = "ANY",
#     CPDB_database_file = "ANY",
#     files_edges_path = "ANY",
#     centralities_file = "ANY",
#     TF_Database = "ANY",
#     file_extension = "ANY",
#     Kinome_database = "ANY",
#     BioMart_Dataset = "character",
#     Client_email = "ANY",
#     weight_penalty = "numeric",
#     circos = "list",
#     circos_aux = "list",
#     pathways = "list"
#   )
# )

Cricos_plots <- function(
  working_dir,
  NOODAI_object,
  edges_dir,
  TF_Database,
  file_extension,
  BioMart_Dataset,
  BiomaRT_selected_organisms
) {

  library(readr)
  library(RColorBrewer)
  library(circlize)
  library(stringr)
  library(ComplexHeatmap)
  library(gridBase)
  library(gtools)
  library(biomaRt)



  tryCatch(
    {
      mart <- useMart('ENSEMBL_MART_ENSEMBL')
      mart <- useDataset(BioMart_Dataset, mart)
    },
    error = function(e){
      tryCatch(
        {
          mart <- useMart('ENSEMBL_MART_ENSEMBL',host='https://asia.ensembl.org')
          mart <- useDataset(BioMart_Dataset, mart)
        },
        error = function(e){
          mart <- useMart('ENSEMBL_MART_ENSEMBL',host='https://useast.ensembl.org')
          mart <- useDataset(BioMart_Dataset, mart)
        }
      )
    }
  )

  setwd(paste(working_dir, "/Centralities_modules_links", sep = ""))

  centralities_modules_files <- list.files(
    path = getwd(),
    pattern = "centralities_modules_links.txt"
  )

  edge_modules_files <- list.files(
    path = edges_dir,
    pattern = paste0(file_extension, "_")
  )

  BiomaRT_organisms <- read.table(
    BiomaRT_selected_organisms,
    sep = "\t",
    quote = "",
    header = TRUE
  )
  BiomaRT_organisms <- BiomaRT_organisms[which(BiomaRT_organisms$BiomaRT %in% BioMart_Dataset), ]
  BiomaRT_organisms$Name <- stringr::str_replace_all(BiomaRT_organisms$Name, " ", "_")
  
  tf_list <- read.table(
    TF_Database,
    sep = "\t",
    fill = TRUE,
    header = TRUE
  )

  if(length(unique(tf_list$Species)) > 1) {

    tf_list <- tf_list[which(tolower(tf_list$Species) %in% tolower(BiomaRT_organisms$Name)),]

  }

  for (i in 1:length(centralities_modules_files)) {
    
    
    data <- read.table(
      centralities_modules_files[i],
      header = TRUE,
      fill = TRUE,
      sep = "\t"
    )
    
    data_10pc <- data[1:floor(nrow(data) / 3), c(1, 2, 4)]
    
    data <- data_10pc

    # Remove duplicates to avoid issues in table format
    data <- data[!duplicated(data), ]
    
    aux <- str_split(string = centralities_modules_files[i], pattern = "_")
    aux <- aux[[1]][1]
    
    ind <- grep(aux,edge_modules_files)
    edge_data <- read.table(
      paste0(
        edges_dir,
        "/",
        edge_modules_files[ind]
      ),
      header = FALSE,
      fill = TRUE,
      sep = "\t"
    )
    
    ind <- which(data$SumModule > 15)
    if (length(ind) > 0) {data <- data[-ind, ]}

    ind <- which(tolower(data$Gene) %in% tolower(tf_list$Symbol))
    
    if(length(ind) > 7) {ind <- ind[1:7]}

    data_tf <- c()

    for (j in 1:length(ind)) {
      ind1 <- which(edge_data$V1 %in% data$Gene[ind[j]])
      ind2 <- which(edge_data$V2 %in% data$Gene[ind[j]])
      aux <- unique(edge_data[union(ind1,ind2), c(1, 2)])

      ind1 <- which(data$Gene %in% aux$V1)
      ind2 <- which(data$Gene %in% aux$V2)
      aux2 <- unique(data$Gene[union(ind1, ind2)])
      aux2 <- aux2[-which(tolower(aux2) %in% tolower(tf_list$Symbol))]
      
      if(length(aux2) > 0) {
        data_tf <- rbind(
          data_tf,
          data.frame(
            aux2,
            data$Gene[ind[j]]
          )
        )
      }

    }
    
    if(length(unique(data_tf[, 1])) > 30) {
      
      ind <- which(data_10pc[, 1] %in% data_tf[, 1])
      if(length(ind) > 25) {ind <- ind[1:25]}
      data_tf <- data_tf[which(data_tf[,1] %in% data_10pc[ind, 1]), ]

    }

    if(is.null(data_tf)) {next()}
    
    colnames(data_tf) <- c("Protein", "TF")
    data_tf <- cbind(
      data_tf,
      data.frame(
        weight = rep(1,nrow(data_tf)),
        SubModule = sapply(
          data_tf$Protein,
          FUN = function(x) {
            data[which(data[,1] %in% x),2]
          }
        )
      )
    )
    
    data_tf <- data_tf[order(data_tf$Protein), ]
    
    color_palette_members <- grDevices::colorRampPalette(
      c(
        "#7196BE",
        "#BD3737",
        "#CD9B1D",
        "#E8768F",
        "#855C5C",
        "#FFD2DC",
        "#663355",
        "#55CCAA",
        "#CFCF3F",
        "#99DAAF",
        "#CF73B8",
        "#E8890C",
        "#645ED4"
      )
    )
    
    color_palette_Group <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(5, "Set1"))
    
    ##########################################
    
    grid.col <- color_palette_members(length(unique(data_tf[, 4])))
    names(grid.col) <- mixedsort(as.character(unique(data_tf[, 4])), )
    
    vec_cols <- c(
      "#7196BE",
      "#BD3737",
      "#CD9B1D",
      "#E8768F",
      "#855C5C"
    )

    for (imk in 1:5) {
      if(length(which(names(grid.col)[imk] %in% as.character(c(1:5)))) > 0) {
        
        grid.col[imk]<-vec_cols[as.numeric(names(grid.col)[imk])]
        
      }

    }

    grid.col_aux <- grid.col
    grid.col <- sapply(
      data_tf[, 4],
      FUN = function(x) {
        grid.col[which(names(grid.col) %in% x)]
      }
    )

    names(grid.col) <- data_tf[, 1]

    aux <- as.character(
      data[which(data$Gene %in% data_tf$TF), 2]
    )
    names(aux) <- data[which(data$Gene %in% data_tf$TF), 1]

    grid.col <- c(
      grid.col_aux[match(aux,unique(names(grid.col_aux)))], grid.col[match(unique(names(grid.col)), names(grid.col))]
    )

    names(grid.col) <- c(names(aux),unique(data_tf[, 1]))

    ##
    circos_slot_name = strsplit(centralities_modules_files[i], "_")[[1]][1]

    circos_aux_slot_data = list(name = data)
    names(circos_aux_slot_data) = circos_slot_name
    NOODAI_object@circos_aux = c(NOODAI_object@circos_aux, circos_aux_slot_data)

    circos_slot_data = list(name = data_tf)
    names(circos_slot_data) = circos_slot_name
    NOODAI_object@circos = c(NOODAI_object@circos, circos_slot_data)
    
    data_circos <- data_tf[, c(1, 2, 3)]

    circos_plot1 <- function() {
      
      circos.clear()
      circos.par(start.degree = 180)
      circos.initialize(sectors = "a", xlim = c(0,length(unique(data_circos[,1]))))

      par(cex = 1.5, mar = c(0, 0, 0, 0))
      chordDiagram(
        data_circos,
        transparency = 0.2,
        order = c(
          data_circos[, 1],
          sort((unique(data_circos[,2])), decreasing = TRUE)
        ),
        annotationTrack = c("name", "grid"),
        scale = TRUE,
        grid.col = grid.col,
        big.gap = 10,
        annotationTrackHeight = c(0.03, 0.03)
      )

      aaux1 <- unique(data_tf[, 4])
      aaux1 <- sort(aaux1)
      aaux2 <- grid.col_aux[unique(match(unique(data_tf[, 4]), unique(names(grid.col_aux))))]
      aaux2 <- aaux2[mixedsort(names(aaux2))]
      lgd_links = Legend(
        labels = aaux1,
        legend_gp = gpar(fill = aaux2),
        title_position = "topleft",
        title = "Modules",
        labels_gp = gpar(fontsize = 17),
        title_gp = gpar(fontsize = 19, fontface = "bold")
      )

      draw(
        lgd_links,
        x = unit(1, "npc") - unit(3, "mm"),
        y = unit(13, "mm"),
        just = c("right", "bottom")
      )

    }
    
    dir.create(paste0(getwd(), "/Images"), recursive = T)
    
    save_name <- paste(
      getwd(),
      "/Images/",
      sapply(
        strsplit(centralities_modules_files[i], "_"),
        "[",
        1
      ),
      "_chordDiagram.pdf",
      sep = ""
    )

    if (file.exists(save_name)) {
      file.remove(save_name)
    }
    pdf(
      file = save_name,
      onefile = FALSE,
      width = measurements::conv_unit(x = 320, from = "mm", to = "inch"),
      height = measurements::conv_unit(x = 250, from = "mm", to = "inch")
    )
    circos_plot1()
    dev.off()

  }

  # return(1)
  return(list(1, NOODAI_object))

}