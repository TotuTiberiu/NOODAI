working_dir <- "~/NOODAI_shiny"
#NOODAI_object <- 
phenotype_names <- c("M1","M2a","M2c")
phenotype_comparison <- c("M1vsM2a","M2avsM1")
files_edges_path = NULL
centralities_file = NULL
TF_Database = NULL
file_extension = NULL
Results_index <- "2025-09-08_1_4_FP8QH8"
Kinome_database = NULL
BioMart_Dataset <- "hsapiens_gene_ensembl"
Cytoscape_style_file = NULL
file_DEA_names <- c("~/NOODAI_shiny/Results_2025-07-21_1_4_J1G794/tmp/DGE.xlsx",
                    "~/NOODAI_shiny/Results_2025-07-21_1_4_J1G794/tmp/DTU.xlsx",
                    "~/NOODAI_shiny/Results_2025-07-21_1_4_J1G794/tmp/PhosphoProteomics.xlsx",
                    "~/NOODAI_shiny/Results_2025-07-21_1_4_J1G794/tmp/Proteomics.xlsx")
flag_run_cytoscape = 1


setClass(
  "NOODAI",
  representation(
    token = "character",
    working_dir = "character",
    BioGRID_data_file = "ANY",
    STRING_data_file = "ANY",
    IntAct_data_file = "ANY",
    file_DEA_names = "character",
    phenotype_names = "character",
    phenotype_comparison = "character",
    Use_precompiled_database = "numeric",
    LookUp_table_file = "ANY",
    Results_index = "character",
    edge_file_path = "ANY",
    monet_path = "ANY",
    Monet_method_string = "character",
    tmp_bin_folder = "ANY",
    CPDB_databases = "ANY",
    MONET_background_file = "ANY",
    CPDB_database_file = "ANY",
    files_edges_path = "ANY",
    centralities_file = "ANY",
    TF_Database = "ANY",
    file_extension = "ANY",
    Kinome_database = "ANY",
    BioMart_Dataset = "character",
    Client_email = "ANY",
    weight_penalty = "numeric",
    circos = "list",
    circos_aux = "list",
    pathways = "list",
    cytoscape_modules = "list",
    status = "character"
  )
)

NOODAI_object = new("NOODAI")