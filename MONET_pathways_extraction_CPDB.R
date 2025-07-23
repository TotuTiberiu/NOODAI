#The function extract the pathways associated with each module identified after the MONET decomposition.
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


#' Perform pathway enrichment analysis for MONET modules using CPDB based databases
#'
#' This function processes MONET module results by mapping module genes to Entrez IDs and 
#' performing pathway enrichment analysis using data from CPDB. 
#' The analysis is parallelized for efficiency and outputs an Excel 
#' file per module with enriched signaling pathways.
#'
#' @description
#' The function reads MONET module output files (named with `result-modules*`), maps gene symbols 
#' to Entrez IDs using BioMart or UniProt/ChEBI pre-mappings, and runs over-representation analysis (ORA) 
#' using 'clusterProfiler'. The results are written to Excel spreadsheets, one per module, in a 
#' 'Signaling_Pathways' subdirectory.
#'
#' The function is optimized for multi-core environments.
#'
#' @param working_dir Character. Path to the working directory containing MONET module outputs and where results will be saved.
#' @param CPDB_database_file Character. Path to the CPDB pathway-to-gene mapping file ('CPDB_pathways_genes.tab').
#' @param CPDB_databases Character vector. Names of the CPDB databases to include (e.g., '"KEGG"', '"Reactome"').
#' @param MONET_background_file Character. Path to Excel file containing background Entrez genes per cluster (one sheet per cluster).
#' @param phenotype_names Character vector. Names of phenotypes used in the comparison (e.g., `c("Control", "Disease")`).
#' @param phenotype_comparison Character. A string like `"DiseaseVsControl"` indicating which groups were compared.
#' @param BioMart_Dataset Character. The BioMart dataset name (e.g., '"hsapiens_gene_ensembl"').
#' @param BiomaRT_selected_organisms Character. Path to a text file listing BioMart-supported organisms and dataset names.
#' @param Uniprot_Ncbi_map Character. Path to UniProt to NCBI gene mapping file.
#' @param ChEBI_map Character. Path to ChEBI to gene mapping file.
#'
#' @return Integer value '1' indicating successful execution. Enrichment results are written to Excel files in the 'Signaling_Pathways/' directory.
#'
#' @details
#' - Input files named 'result-modules*' are assumed to be present in the working directory.
#' - If BioMart mapping fails, fallback mapping is done using UniProt and ChEBI datasets.
#' - Pathway enrichment is performed with 'clusterProfiler::enricher()' using TERM2GENE format.
#' - Only clusters with sufficient genes are processed (minimum 10 genes by default).
#' - Results are saved in '.xlsx' format, one sheet per module.
#'
#' @examples
#' MONET_Pathways_extraction(
#'   working_dir = "path/to/results_dir",
#'   CPDB_database_file = "path/to/CPDB_pathways_genes.tab",
#'   CPDB_databases = c("KEGG", "Reactome"),
#'   MONET_background_file = "path/to/Background_total.xlsx",
#'   phenotype_names = c("Control", "Treated"),
#'   phenotype_comparison = "TreatedVsControl",
#'   BioMart_Dataset = "hsapiens_gene_ensembl",
#'   BiomaRT_selected_organisms = "path/to/BiomaRT_selected_organisms.txt",
#'   Uniprot_Ncbi_map = "path/to/Unip_NCBI_map.txt",
#'   ChEBI_map = "path/to/ChEBI_entries_map.txt"
#' )
#'
#' @import readxl biomaRt stringr parallel org.Hs.eg.db AnnotationDbi msigdbr rbioapi clusterProfiler

MONET_Pathways_extraction <- function(
  working_dir,
  CPDB_database_file,
  CPDB_databases,
  MONET_background_file,
  phenotype_names,
  phenotype_comparison,
  BioMart_Dataset,
  BiomaRT_selected_organisms,
  Uniprot_Ncbi_map,
  ChEBI_map
) {
  
  library(readxl)
  library(biomaRt)
  library(stringr)
  library(parallel)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(msigdbr)
  library(rbioapi)
  library(clusterProfiler)

  BiomaRT_organisms <- read.table(
    BiomaRT_selected_organisms,
    sep = "\t",
    quote = "",
    header = TRUE
  )

  BiomaRT_organisms <- BiomaRT_organisms[which(BiomaRT_organisms$BiomaRT %in% BioMart_Dataset), ]
  flag_biomart <- 0
  
  if(length(BiomaRT_organisms) == 0) {flag_biomart <- 1}
  
  uni_ncbi_mapping <- read.table(
    Uniprot_Ncbi_map,
    sep = "\t",
    header = TRUE,
    quote = ""
  )

  uni_ncbi_mapping <- uni_ncbi_mapping[which(uni_ncbi_mapping$Organism %in% BiomaRT_organisms$Name), c(1, 2, 3)]
  colnames(uni_ncbi_mapping) <- c("Uniprot_ChEBI", "Name", "NCBI_ChEBI")
  
  ChEBI_mapping <- read.table(
    ChEBI_map,
    sep = "\t",
    header = TRUE,
    comment.char = "@",
    quote = ""
  )

  colnames(ChEBI_mapping) <- c("Uniprot_ChEBI", "Name", "NCBI_ChEBI")
  Id_mapp_data <- rbind(uni_ncbi_mapping, ChEBI_mapping)

  setwd(working_dir)
  n_cores <- detectCores()
  n_cores_used <- floor(n_cores / 2)
  cluster_1 <- makeCluster(n_cores_used, type = "PSOCK")
  clusterEvalQ(cluster_1, library(readxl))
  clusterEvalQ(cluster_1, library(biomaRt))
  clusterEvalQ(cluster_1, library(org.Hs.eg.db))
  clusterEvalQ(cluster_1, library(AnnotationDbi))
  clusterEvalQ(cluster_1, library(msigdbr))
  clusterEvalQ(cluster_1, library(rbioapi))
  clusterEvalQ(cluster_1, library(clusterProfiler))
  clusterEvalQ(cluster_1, library(stringr))
  
  options(timeout = max(1000, getOption("timeout")))
  options(java.parameters = "-Xmx6000m")
  set.seed(21)
  
  BiomaRT_organisms <- read.table(
    BiomaRT_selected_organisms,
    sep = "\t",
    quote = "",
    header = TRUE
  )

  BiomaRT_organisms <- BiomaRT_organisms[which(BiomaRT_organisms$BiomaRT %in% BioMart_Dataset), ]
  
  cpdb_db <- read.table(
    CPDB_database_file,
    fill = TRUE,
    sep = "\t",
    header = TRUE,
    quote = ""
  )

  cpdb_db <- cpdb_db[which(tolower(cpdb_db$source) %in% tolower(CPDB_databases)), ]
  cpdb_db <- cpdb_db[which(tolower(cpdb_db$Organism) %in% tolower(BiomaRT_organisms$Name)), ]
  cpdb_db <- cpdb_db[, c(1:4)]
  cpdb_db_aux <- data.frame(matrix(nrow = 1, ncol = 4))
  colnames(cpdb_db_aux) <- colnames(cpdb_db)
  
  for(i in 1:nrow(cpdb_db)) {
    aux <- unlist(str_split(cpdb_db[i,4],','))
    cpdb_db_aux <- rbind(
      cpdb_db_aux,
      data.frame(
        pathway = rep(cpdb_db[i,1],length(aux)),
        external_id = rep(cpdb_db[i,2],length(aux)),
        source = rep(cpdb_db[i,3],length(aux)),
        entrez_gene_ids = aux
      )
    )
  }
  
  cpdb_db <- cpdb_db_aux[2:nrow(cpdb_db_aux), ]
  
  files_monet <- list.files(path = getwd(), pattern = "result-modules")
  
  length_files <- length(files_monet)
  print(head(length_files))
  nx_per_worker <- floor(length_files / n_cores_used)
  if(nx_per_worker == 0) {nx_per_worker <- 1}
  steps <- seq(1, length_files, nx_per_worker)
  att <- 0
  while(att <= 5) {
    att <- att + 1
    try({
      tryCatch(
      {
        mart <- useMart('ENSEMBL_MART_ENSEMBL')
        mart <- useDataset(BioMart_Dataset, mart)
      },
      error = function(e){
        tryCatch(
          {
            mart <- useMart('ENSEMBL_MART_ENSEMBL', host = 'https://asia.ensembl.org')
            mart <- useDataset(BioMart_Dataset, mart)
          },
          error = function(e){
            mart <- useMart('ENSEMBL_MART_ENSEMBL', host = 'https://useast.ensembl.org')
            mart <- useDataset(BioMart_Dataset, mart)
          })
      })
    })
  }

  fsx <- function(files_index) {
    
    if (files_index <= (length_files - nx_per_worker + 1)) {
      
      files_index <- seq(files_index, (files_index + nx_per_worker - 1), by = 1)
    
    } else {
      
      files_index <- seq(files_index, length_files, by = 1)
      
    }
    file_universe_entrez <- MONET_background_file
    
    dir.create(
      paste(
        getwd(),
        "/Signaling_Pathways",
        sep = ""
      ),
      recursive = T
    )
    results_dir <- paste(
      getwd(),
      "/Signaling_Pathways",
      sep = ""
    )
    
    phenotype_comparison_aux = t(
      as.data.frame(
        strsplit(phenotype_comparison, split = "vs")
      )
    )
    rownames(phenotype_comparison_aux) <- phenotype_comparison
    
    for (i in files_index) {
      
      data <- read.table(
        files_monet[i],
        sep = "\t",
        header = FALSE,
        fill = TRUE,
        col.names = paste0("V", seq_len(max(count.fields(files_monet[i]))))
      )
      
      data_resave <- data

      name_index <- str_extract(files_monet[i], "(?<=_)[^_]+(?=\\.txt)")
      
        
      universe_entrez <- read_excel(
        file_universe_entrez,
        sheet = name_index,
        col_names = FALSE
      )

      colnames(universe_entrez) <- c("EntrezID")
      universe_entrez$EntrezID <- as.character(universe_entrez$EntrezID)
      ind <- which(is.na(universe_entrez$EntrezID))
      
      if(length(ind) > 0) {
        universe_entrez <- universe_entrez[-ind, ]
      }
      
      data <- data[, c(3:ncol(data))]
      
      ind <- apply(data, 1, function(x) length(which(x == "")))
      keep <- which(ind <= (ncol(data) - 10))
      if(length(keep) == 0) {next}
      data <- data[keep, ]
      
      data_resave <- data_resave[order(ind), ]
      data_resave[, 1] <- seq_len(nrow(data_resave))
      file.remove(files_monet[i])
      write.table(
        data_resave,
        files_monet[i],
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = "\t"
      )
        
      ss <- str_split(files_monet[i],patter="[.]")
      save_name_CPDB = paste(
        results_dir,
        "/",
        ss[[1]][1],
        "_Signaling.xlsx",
        sep = ""
      )
      
      if (file.exists(save_name_CPDB)) {
        file.remove(save_name_CPDB)
      }
        
      for (j in 1:nrow(data)) {
        
        aux <- data[j, ]
        ind <- which(aux == "")
        if (length(ind) > 0) {
          aux <- aux[-ind]
        }

        if(flag_biomart) {
        aux_converted <- getBM(
          filters = "external_gene_name",
          attributes = c("external_gene_name", "entrezgene_id","uniprot_gn_id"),
          values = as.character(aux),
          mart = mart,
          useCache = FALSE
        )

        } else {
          
          aux_converted <- data.frame(
            entrezgene_id = Id_mapp_data$NCBI_ChEBI[which(Id_mapp_data$Name %in% aux)]
          )
          
        }
        
        KeggReactome_ora_results <- c()
        KeggReactome_ora_results <- enricher(
          gene = aux_converted$entrezgene_id,
          pvalueCutoff = 0.1,
          qvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          universe = universe_entrez$EntrezID,
          TERM2GENE = dplyr::select(
            cpdb_db,
            pathway,
            entrez_gene_ids
          ),
          minGSSize = round(1 * length(aux_converted$entrezgene_id) / 100),
          maxGSSize = 300
        )
        
        if(is.null(KeggReactome_ora_results) == FALSE) {

          ind <- which(KeggReactome_ora_results@result$p.adjust < 1)
          KeggReactome_ora_results <- data.frame(
            Pathway = KeggReactome_ora_results@result$ID[ind],
            Source = cpdb_db$source[match(KeggReactome_ora_results@result$ID[ind], cpdb_db$pathway)],
            External_id = cpdb_db$external_id[match(KeggReactome_ora_results@result$ID[ind], cpdb_db$pathway)],
            padj_values = KeggReactome_ora_results@result$p.adjust[ind],
            pvalue = KeggReactome_ora_results@result$pvalue[ind],
            members_input_overlap = KeggReactome_ora_results@result$geneID[ind],
            geneRatio = KeggReactome_ora_results@result$GeneRatio[ind],
            bcgRatio = KeggReactome_ora_results@result$BgRatio[ind]
          )
          
          KeggReactome_ora_results$members_input_overlap <- str_replace_all(
            KeggReactome_ora_results$members_input_overlap,
            pattern = "/",
            replacement = ";"
          )

        }
        
        if(length(KeggReactome_ora_results) > 0) {

          xlsx::write.xlsx(
            KeggReactome_ora_results,
            save_name_CPDB,
            sheetName = paste("Cluster_", as.character(j), sep = ""),
            append = TRUE
          )
        } 
      }
    }
  }
  
  dir.create('Signaling_Pathways', recursive = T)
  
  clusterExport(cluster_1, "mart", envir = environment())
  clusterExport(cluster_1, "nx_per_worker", envir = environment())
  clusterExport(cluster_1, "length_files", envir = environment())
  clusterExport(cluster_1, "files_monet", envir = environment())
  clusterExport(cluster_1, "cpdb_db", envir = environment())
  clusterExport(cluster_1, "MONET_background_file", envir = environment())
  clusterExport(cluster_1, "phenotype_comparison", envir = environment())
  
  parSapply(cluster_1, steps, fsx)
  
  stopCluster(cluster_1)
  return(1)
  
}

