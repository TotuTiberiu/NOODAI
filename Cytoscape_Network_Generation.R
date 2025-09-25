#The function constructs automatically the protein-protein interaction network. It requires a graphical interface!

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

Cytoscape_networks <- function(
  Centrality_file,
  Centrality_file_sheet,
  names_of_omicsDatasets,
  Edge_file,
  Edge_file_sheet = NULL,
  Cytoscape_style_file,
  save_name_network
) {
  
  library(RCy3)
  library(igraph)
  library(stringr)
  
  options(rcy3.cytoscape.host = "localhost")
  options(rcy3.cytoscape.port = 1234) 
  
  Sys.sleep(10)
  
  #Identify if it is a txt file (full-size netowrk) or a module that is stored in Excel. Read the edges files
  sep_ident <- strsplit(Edge_file, split = "\\.")
  sep_ident <- sep_ident[[1]][lengths(sep_ident)]
  
  if (sep_ident == "txt") {
    data <- read.table(
      Edge_file,
      header = FALSE,
      sep = "\t",
      quote = "",
      fill = TRUE
    )
  }
  
  if (sep_ident == "xlsx") {
    data <- openxlsx::read.xlsx(
      Edge_file,
      sheet = Edge_file_sheet,
      colNames = FALSE,
      rowNames = FALSE
    )
  }
  
  #Read the centrality data and the omics presence info
  auxiliary_data <- openxlsx::read.xlsx(
    Centrality_file,
    sheet = Centrality_file_sheet,
    colNames = TRUE,
    rowNames = FALSE
  )
  
  #Select only the centrality column and omics presence
  selected_col <- c(
    which(
      colnames(auxiliary_data) %in% "Symbol"
    ),
    which(
      colnames(auxiliary_data) %in% "Current_Flow_Betweenness_Centrality"
    ),
    which(
      colnames(auxiliary_data) %in% names_of_omicsDatasets
    )
  )

  auxiliary_data <- auxiliary_data[, selected_col]
  auxiliary_data$Omics_Total <- apply(
    auxiliary_data[, which(colnames(auxiliary_data) %in% names_of_omicsDatasets)],
    MARGIN = 1,
    FUN = function(x) {
      paste(x, collapse = "")
    }
  )
  
  #Manage if there are more or less than 4 omics datasets
  nr_omics <- length(
    which(
      colnames(auxiliary_data) %in% names_of_omicsDatasets
    )
  )
  
  if (nr_omics < 4) {
    nr <- 4 - length(
      which(
        colnames(auxiliary_data) %in% names_of_omicsDatasets
      )
    )
    auxiliary_data$Omics_Total <- apply(
      cbind(
        matrix(
          0,
          nrow = nrow(auxiliary_data),
          ncol = nr
        ),
        auxiliary_data$Omics_Total
      ),
      MARGIN = 1,
      FUN = function(x) {
        paste(x, collapse = "")
      }
    )
  }

  if (nr_omics > 4) {
    auxiliary_data$Omics_Total <- substring(
      auxiliary_data$Omics_Total,
      nchar(auxiliary_data$Omics_Total) - 3,
      nchar(auxiliary_data$Omics_Total)
    )
  }
  
  #Create the graph (undirected)
  grph_data <- as.matrix(data[, c(1, 2)])
  grph_data <- t(apply(grph_data, 1, sort))
  grph_data <- unique(grph_data)
  grph_data <- apply(grph_data, 2, as.character)
  grh <- graph_from_edgelist(grph_data, directed = FALSE)
  
  network <- createNetworkFromIgraph(
    grh,
    title = save_name_network,
    collection = save_name_network
  )
  Sys.sleep(5)
  #Import the preset Cytoscape style
  imported_styles <- importVisualStyles(filename = Cytoscape_style_file)
  Sys.sleep(5)
  style_name <- imported_styles[1]  # use the first style imported
  setVisualStyle(style_name)
  Sys.sleep(5)
  
  #Tailor the Centrality and Omics presence info to match the preset style.
  #Remove duplicates
  auxiliary_data <- auxiliary_data[!duplicated(auxiliary_data[, 1]), ]
  ind <- match(unique(auxiliary_data$Symbol), auxiliary_data$Symbol)
  rownames(auxiliary_data) <- auxiliary_data[ind, 1]
  auxiliary_data <- auxiliary_data[, -1]
  new_col_names <- names_of_omicsDatasets
  
  if (nr_omics > 4) {
    nr_omics <- 4
  }
  colnames(auxiliary_data)[2:(ncol(auxiliary_data) - 1)] <- new_col_names[c(1:nr_omics)]
  loadTableData(
    auxiliary_data,
    table = "node",
    table.key.column = "shared name"
  )
  Sys.sleep(5)
  
  #Set a layout
  if(sep_ident == "txt") {
    layoutNetwork("force-directed")
  }

  if(sep_ident == "xlsx"){
    layoutNetwork("degree-circle")
  }

  Sys.sleep(5)
  
  #Resize the nodes
  Sys.sleep(5)
  setNodeSizeMapping(
    table.column = "Current_Flow_Betweenness_Centrality",
    table.column.values = c(
      min(
        auxiliary_data[, 1],
        na.rm = TRUE
      ),
      max(
        auxiliary_data[, 1],
        na.rm = TRUE
      )
    ),
    sizes = c(30, 150),
    mapping.type = "c",
    default.size = 20,
    style.name = style_name
  )
  
  Sys.sleep(5)
  save_name_network_s <- sub("^Cluster", "Module", save_name_network)
  #Save the file
  saveSession(filename = paste0(save_name_network_s))
  Sys.sleep(5)
  exportImage(filename = save_name_network_s, type = "SVG", width = 1000, height = 800)
  Sys.sleep(5)
  
  try({
    deleteNetwork(network = save_name_network)
  }, silent = TRUE)
  
  # Delete all networks in the collection
  try({
    networks_in_collection <- getNetworkList()
    for (net in networks_in_collection) {
      net_suid <- getNetworkSuid(net)
      net_collection <- getCollectionNetworkSuid(net_suid)
      if (!is.null(net_collection) && net_collection == save_name_network) {
        deleteNetwork(net)
      }
    }
  }, silent = TRUE)
  #gc()
  
  return(1)
  
}


Generate_cytoscape_Legend <- function(
  Centrality_file,
  Centrality_file_sheet,
  names_of_omicsDatasets,
  save_name
) {
  
  library(grid)
  library(gridSVG)
  auxiliary_data <- openxlsx::read.xlsx(
    Centrality_file,
    sheet = Centrality_file_sheet,
    colNames = TRUE,
    rowNames = FALSE
  )
  
  centrality_range <- auxiliary_data$Current_Flow_Betweenness_Centrality
  centrality_range <- as.character(
    c(
      round(
        min(centrality_range, na.rm = TRUE)
      ),
      round(
        max(centrality_range,na.rm = TRUE)
      )
    )
  )
  size_small <- 5
  size_large <- 17
  used_colors <- c("#8DABCB", "#D9BF1B", "#E8890C", "#AE3333")
  
  if (length(names_of_omicsDatasets) < 5) {
    category_colors <- used_colors[1:length(names_of_omicsDatasets)]
    category_labels <- names_of_omicsDatasets
  } else {
    category_colors <- used_colors
    category_labels <- names_of_omicsDatasets[1:4]
  }
  
  svg(save_name, width = 4.5, height = 6)
  
  grid.newpage()
  
  # Title
  grid.text(
    "Node size ~ Centrality",
    x = unit(0.5, "npc"),
    y = unit(0.98, "npc"),
    gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  # Gradient bar parameters
  n_steps <- 100
  y_pos <- 0.93
  height <- 0.05
  width <- 0.8
  
  # Draw gradient bar (black to white from left to right)
  for (i in 1:n_steps) {
    x_center <- 0.1 + width * (i - 0.5) / n_steps
    shade <- gray(1 - (i - 1) / (n_steps - 1))  # black to white
    
    grid.rect(
      x = unit(x_center, "npc"),
      y = unit(y_pos, "npc"),
      width = unit(width / n_steps, "npc"),
      height = unit(height, "npc"),
      just = "center",
      gp = gpar(fill = shade, col = NA)
    )
  }
  
  grid.text(
    centrality_range[1],
    x = unit(0.1, "npc"),
    y = unit(y_pos - 0.04, "npc"),
    just = "center",
    gp = gpar(fontsize = 10)
  )
  
  grid.text(
    centrality_range[2],
    x = unit(0.9, "npc"),
    y = unit(y_pos - 0.04, "npc"),
    just = "center",
    gp = gpar(fontsize = 10)
  )
  
  y_circle <- y_pos - 0.12
  
  grid.circle(
    x = unit(0.1, "npc"),
    y = unit(y_circle, "npc"),
    r = unit(size_small, "points"),
    gp = gpar(fill = "gray20", col = "black")
  )

  grid.circle(
    x = unit(0.5, "npc"),
    y = unit(y_circle, "npc"),
    r = unit((size_small + size_large)/2, "points"),
    gp = gpar(fill = "gray20", col = "black")
  )

  grid.circle(
    x = unit(0.9, "npc"),
    y = unit(y_circle, "npc"),
    r = unit(size_large, "points"),
    gp = gpar(fill = "gray20", col = "black")
  )
  
  circle_radius <- unit(23, "points")
  circle_radius_npc <- convertUnit(circle_radius, "npc", valueOnly = TRUE)
  spacing <- circle_radius_npc * 2.5 
  
  y_start <- y_circle - 0.15
  
  for (i in seq_along(category_colors)) {
    
    y_cat_pos <- y_start - (i - 1) * spacing
    
    grid.circle(
      x = unit(0.5, "npc"),
      y = unit(y_cat_pos, "npc"),
      r = circle_radius,
      gp = gpar(fill = category_colors[i], col = "black")
    )

    grid.text(
      category_labels[i],
      x = unit(0.5, "npc"),
      y = unit(y_cat_pos - circle_radius_npc - 0.02, "npc"),
      just = "center",
      gp = gpar(fontsize = 10)
    )

  }
  
  dev.off()
  
  return(1)
  
}