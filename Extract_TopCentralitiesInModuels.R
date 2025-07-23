#The function extracts the centrality score for the elements from each module after the MONET decomposition. These are used for the circular diagrams.

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


#' Extract Centrality Metrics and Link to Network Modules
#'
#' This function processes node centrality values (current-flow betweenness centrality) and associates them 
#' with MONET decomposed modules.
#'
#' @param working_dir Character. The path to the working directory where MONET module results are located.
#' @param centralities_file Character. Path to the Excel file containing centrality metrics per phenotype. 
#'        Each sheet should correspond to a phenotype and contain at least the columns 'Symbol' and 
#'        'Current_Flow_Betweenness_Centrality'.
#' @param phenotype_names Character vector. Names of the phenotypes used to match the appropriate Excel sheets 
#'        and MONET module files.
#' @param file_extension Character. Suffix pattern used in MONET module filenames (e.g., '"AcusVsBrancus"').
#'
#' @details
#' For each phenotype sheet in the Excel file, the function:
#' - Reads and ranks genes by current-flow betweenness centrality from an Excel file (one sheet per phenotype).
#' - Finds the module to which each gene belongs based on MONET module files.
#' - Outputs a file listing gene-module assignments with centrality scores and a weight
#'
#' Output files are saved to the 'Centralities_modules_links' folder inside the working directory.
#' Genes not present in any module are omitted from the output.
#'
#' @return Integer of unity.
#'
#' @examples

#' Extract_CentralitiesValues(
#'   working_dir = "/path/to/MONET/results",
#'   centralities_file = "centrality_scores.xlsx",
#'   phenotype_names = c("Acus", "Brancus"),
#'   file_extension = "AcusVsBrancus"
#' )
#'
#' @import readr stringr readxl

Extract_CentralitiesValues <- function(
  working_dir,
  centralities_file,
  phenotype_names,
  file_extension
) {

  library(readr)
  library(stringr)
  library(readxl)
    
  setwd(working_dir)

  sheets <- excel_sheets(centralities_file)

  module_files <- list.files(
    path = getwd(),
    pattern = paste0("result-modules__",file_extension)
  )

  for (i in 1:length(sheets)) {
    
    centralities <- read_excel(centralities_file, sheet = sheets[i])
    centralities$Current_Flow_Betweenness_Centrality <- as.numeric(as.character(centralities$Current_Flow_Betweenness_Centrality))
    centralities$Current_Flow_Betweenness_Centrality <- centralities$Current_Flow_Betweenness_Centrality + c(1:length(centralities$Current_Flow_Betweenness_Centrality)) / length(centralities$Current_Flow_Betweenness_Centrality)
    
    centralities <- centralities[match(sort(centralities$Current_Flow_Betweenness_Centrality,decreasing = TRUE), centralities$Current_Flow_Betweenness_Centrality), match(c("Current_Flow_Betweenness_Centrality", "Symbol"),colnames(centralities))]
    
    m_ind <- grep(sheets[i], module_files)
    
    modules <- read.table(
      module_files[m_ind],
      header = FALSE,
      fill = TRUE,
      sep = "\t"
    )
    
    centralities_moduels <- data.frame(
      matrix(
        nrow = nrow(centralities),
        ncol = 4
      )
    )
    
    for (j in 1:nrow(centralities)) {
      
      if(is.na(centralities[j,2]) == 0) {
        
        aux <-  as.character(
          which(
            apply(
              modules,
              FUN = function(x) {
                which(centralities[j,2] %in% x)
              },
              MARGIN = 1
            ) == 1
          )
        )
        
        if(length(aux) > 0) {
      
          centralities_moduels[j, 1] <- as.character(centralities[j, 2])
          centralities_moduels[j, 2] <- aux
          centralities_moduels[j, 3] <- centralities[j, 1]
          centralities_moduels[j, 4] <- 1

        }
      }
    }
    
    aux <- as.character(module_files[m_ind])
    aux <- sapply(
      strsplit(aux, paste0(file_extension, "_")),
      "[",
      2
    )

    aux <- sapply(
      strsplit(aux, "\\."), 
      "[",
      1
    )

    dir.create(
      paste0(
        getwd(),
        "/Centralities_modules_links"
      ),
      recursive = T
    )

    save_name = paste(
      getwd(),
      "/Centralities_modules_links/",
      aux,
      "_centralities_modules_links.txt",
      sep = ""
    )

    file.remove(save_name)
    colnames(centralities_moduels) <- c("Gene", "SumModule", "CentralityScore", "Weight")
    write.table(
      centralities_moduels,
      file = save_name,
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )
    
  }
  
  return(1)

}