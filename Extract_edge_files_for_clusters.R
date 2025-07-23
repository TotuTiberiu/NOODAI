#The function extracts the edges that are found in each identified subnetwork after the MONET decomposition. These are used for the circular diagrams.

# 
#     Copyright © 2025 , Empa, Tiberiu Totu and Rafael Riudavets Puig.
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



#' Extract and Export Module-Specific Edge Tables from MONET Output
#'
#' This function processes MONET module files and edge files for multiple phenotypes,
#' extracting module-specific edges and saving them as Excel files. It matches nodes
#' from MONET clusters with corresponding edges to generate edge tables per module.
#'
#' @param working_dir Character. Directory where MONET module result files are stored and outputs will be saved.
#' @param phenotype_names Character vector. Names of phenotypes.
#' @param phenotype_comparison Character vector. Names of phenotype comparisons (e.g., '"NightVsDay"').
#' @param files_edges_path Character. Directory path where raw edge list files (used for filtering) are located.
#' @param file_extension Character. File extension or identifier used in the filenames (e.g., `"Total"`).
#'
#' @details
#' - For each phenotype comparison, the function loads MONET modules and filters the full edge list,
#'   saving edges that belong to each module into separate Excel sheets.
#' - Empty or sparse modules (based on missing values) are filtered out.
#' - Creates an 'Edges_tables' subdirectory to store output '.xlsx' files.
#'
#' Requires that edge files are named with a pattern like '<file_extension>_...' and MONET files with 'result-modules__<file_extension>...'.
#'
#' @return Integer unity (for pipeline integration).
#'
#' @examples
#' Extract_edges_tables(
#'   working_dir = "Results/Modules/",
#'   phenotype_names = c("Night", "Day"),
#'   phenotype_comparison = c("NightVsDay"),
#'   files_edges_path = "Data/Edges/",
#'   file_extension = "Total"
#' )
#'
#' @import readxl stringr xlsx

Extract_edges_tables <- function(
  working_dir,
  phenotype_names,
  phenotype_comparison,
  files_edges_path,
  file_extension
) {

  library(readxl)
  library(stringr)

  setwd(working_dir)

  files_monet <- list.files(
    path = getwd(),
    pattern = paste0("result-modules__",file_extension)
  )

  print(files_monet)
  dir.create(
    paste(getwd(),"/Edges_tables", sep = ""),
    recursive = T
  )

  results_dir <- paste(getwd(),"/Edges_tables",sep = "")

  phenotype_comparison_aux = t(
    as.data.frame(
      strsplit(
        phenotype_comparison,
        split = "vs"
      )
    )
  )
  rownames(phenotype_comparison_aux) <- phenotype_comparison


  files_edges <- list.files(
    path = files_edges_path,
    pattern = paste0(file_extension, '_'),
    full.names = TRUE
  )

  print(files_edges)
  edges <- list()

  for (i in 1:length(phenotype_comparison)) {
    ind <- grep(phenotype_comparison[i], files_edges)
    edges[[i]] <- read.table(files_edges[ind])

  }

  print(length(edges))

  for (j in 1:length(phenotype_comparison)) {
    
    ind_f <- grep(phenotype_comparison[j], files_monet)
    print(ind_f)
    print(getwd())
    data <- read.table(
      files_monet[ind_f],
      sep = "\t",
      header = FALSE,
      fill = TRUE
    )

    name_index <- str_split(
      files_monet[ind_f],
      pattern = paste0(file_extension, "_")
    )

    name_index <- str_split(
      name_index[[1]][2],
      pattern = "[.]"
    )
    name_index <- name_index[[1]][1]
    
    data <- data[, c(3:ncol(data))]
    
    ind <- apply(data, 1, function(x) length(which(x == "")))
    keep <- which(ind <= (ncol(data) - 10))
    data <- data[keep, ]
    
    ss <- str_split(
      files_monet[ind_f],
      pattern = "[.]"
    )

    save_name_edges = paste(
      results_dir,
      "/",
      ss[[1]][1],
      "_edges.xlsx",
      sep = ""
    )

    if (file.exists(save_name_edges)) {file.remove(save_name_edges)}
    
    for (k in 1:nrow(data)) {
      
      aux <- data[k, ]
      ind <- which(aux == "")
      if (length(ind) > 0) {aux <- aux[-ind]}
      
      ind1 <- which(edges[[j]]$V1 %in% aux)
      ind2 <- which(edges[[j]]$V2 %in% aux)
      ind <- intersect(ind1,ind2)
      
      aux_edges <- edges[[j]][ind, ]
      
      xlsx::write.xlsx(
        aux_edges,
        save_name_edges,
        sheetName = paste(
          "Cluster_",
          as.character(k),
          sep = ""
        ),
        append = TRUE,
        row.names = FALSE,
        col.names = FALSE
      )
      
    }
    
  }
  
  return(1)

}
