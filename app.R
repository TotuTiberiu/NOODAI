#Main Shiny Interface of the NOODAI software, deployed on the server.


###Modify this
setwd("YOUR_SCRIPTS_DIRECTORY")
###############################



#' NOODAI main interface tailored for local run
#'
#' @description This is the main Rshiny application that one should run to launch a local instance of the Shiny application
#' You should change the working directory path in the first line of the script. In addition, all the github scritps and databases should be available in the same folder.
#' MONET must be pre-configures on machine and set to be run by the Rshiny workers with sudo permisions. The folders must be configured as well to have read/write permissions for the Shiny application.
#' For Cytoscape networks, a headless version of cytoscape should be running at the time the Shiny application is executed. Proper configuration will be dependent on the operating system and available permissions.
#' Full list of required packages for all the functionalities can be found in the Libraries_install.txt file.
#'
#' @noRd 
#'
#' @import shiny reshape2 ggplot2 ggpubr EnvStats moments bslib shinyWidgets archive uuid shinyjs biomaRt purrr dplyr 

js <- "
$(function () {
  $('[data-toggle=tooltip]').tooltip()
})
"

library(shiny)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(moments)
library(bslib)
library(shinyWidgets)
library(archive)
library(uuid)
library(shinyjs)
library(biomaRt)

library(purrr)

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

uniprot_ids_maps <- read.table(
  paste0(
    getwd(),
    "/Databases/Uniprot_FullIds.txt"
  ),
  header = TRUE,
  fill = TRUE
)

biomart_available_datasets <- read.table(
  paste0(
    getwd(),
    "/Databases/BioMart_datasets.txt"
  ),
  header = TRUE,
  fill = TRUE
)

CPDB_databases_list <- read.table(
  paste0(
    getwd(),
    "/Databases/CPDB_databases_list.txt"
  ),
  header = TRUE,
  fill = TRUE
)

Monet_exp <- "--method=[A-Z][0-9] --avgk=?[0-9]+ --linksdir=[a-z]+"

biomart_check <- function() {
  
  tryCatch(
    {
      mart <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL')
      mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart)
      return(1)
    },
    error = function(e){
      tryCatch(
        {
          mart <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL',host='https://asia.ensembl.org')
          mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart)
          return(1)
        },
        error = function(e){
          tryCatch(
            {
              mart <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL',host='https://useast.ensembl.org')
              mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart)
              return(1)
            }, error = function(e){
              return(0)
              stop()
            }
          )
        })
    })
  
}

working_dir <- getwd()

jobs <- list()

save_run_info <- function(pars_ls, out_dir) {

  pars_json = jsonlite::toJSON(pars_ls, auto_unbox = T) %>%
    jsonlite::prettify()

  write(pars_json, file.path(out_dir, "run_info.json"))

}

run_token <- function(
  token,
  working_dir,
  NOODAI_object,
  BioGRID_data_file,
  STRING_data_file,
  IntAct_data_file,
  file_DEA_names,
  phenotype_names,
  phenotype_comparison,
  splicing_file_name,
  Use_precompiled_database,
  LookUp_table_file,
  Results_index,
  edge_file_path,
  monet_path,
  Monet_method_string,
  tmp_bin_folder,
  CPDB_databases,
  MONET_background_file,
  CPDB_database_file,
  files_edges_path,
  centralities_file,
  TF_Database,
  file_extension,
  Kinome_database,
  BioMart_Dataset,
  Client_email,
  weight_penalty,
  Cytoscape_style_file,
  flag_run_cytoscape
) {
  
  source(
    paste0(
      working_dir,
      "/Integrative_analysis.R"
    ),
    local = TRUE,
    echo = TRUE,
    max.deparse.length = 100000
  )
  
  Integrative_Network_analysis(
    working_dir,
    NOODAI_object,
    BioGRID_data_file,
    STRING_data_file,
    IntAct_data_file,
    file_DEA_names,
    phenotype_names,
    phenotype_comparison,
    splicing_file_name,
    Use_precompiled_database,
    LookUp_table_file,
    Results_index,
    BioMart_Dataset,
    weight_penalty
  )
  
  MONET_analysis(
    working_dir,
    NOODAI_object,
    edge_file_path,
    monet_path,
    Monet_method_string,
    tmp_bin_folder,
    Results_index
  )
  
  MONET_pathways(
    working_dir,
    NOODAI_object,
    CPDB_databases,
    MONET_background_file,
    phenotype_names,
    phenotype_comparison,
    CPDB_database_file,
    Results_index,
    BioMart_Dataset
  )
  
  NOODAI_object = Circos_and_auxiliary(
    working_dir = working_dir,
    NOODAI_object = NOODAI_object,
    phenotype_names = phenotype_names,
    phenotype_comparison = phenotype_comparison,
    files_edges_path = files_edges_path,
    centralities_file = centralities_file,
    TF_Database = TF_Database,
    file_extension = file_extension,
    Results_index = Results_index,
    Kinome_database = Kinome_database,
    BioMart_Dataset = BioMart_Dataset,
    Cytoscape_style_file = Cytoscape_style_file,
    file_DEA_names = file_DEA_names,
	flag_run_cytoscape = flag_run_cytoscape
  )

  if(length(Client_email) != 0) {
   # system(
    #  paste0(
    #    "echo ",
    #    "Your analysis associated with the results folder:",
    #    Results_index,
    #    " is done. You can provide this folder ID to the platform or download your results from the following address: omics-oracle.com/www/Downloads/",
    #    Results_index,
    #    ".zip. If the analysis did not work as expected, please check the log files. If the problem persists, please contact the authors. | mutt -s ",
    #    '"NOODAI analysis results are available" ',
    #    Client_email
      )
    )
  }

  return(NOODAI_object)

}

mod_NOODAI_pipeline_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    tags$head(
      tags$script(HTML(js)),
      tags$style(HTML(
        '.info-box {
          min-height: 45px;
          min-width: 45px;
        } 
        .info-box-icon {
          height: 45px;
          line-height: 45px;
          width: 45px;
          line-width: 45px;
        } 
        .info-box-content {
          padding-top: 0px;
          padding-bottom: 0px;
        }'))
    ),
    shiny::fluidRow(
      shinyWidgets::radioGroupButtons(
        inputId = ns("Type_of_analysis"),
        choices = c(
          "Run the full pipeline",
          "Custom algorithms",
          "Results download"
        ),
        status = "primary",
        justified = TRUE
      )
    ),
    conditionalPanel(
      condition = "input.Type_of_analysis == 'Run the full pipeline'",
      ns = ns,
      shinydashboard::box(
        width = 12,
        title = "Input data",
        status = "success",
        shiny::fluidRow(
          shiny::column(
            width = 4,
            fileInput(
              width = "100%",
              inputId = ns("file_DEA_names3"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Omics data archive"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Upload an archive containing at least one omics file. Each omics file needs to be an Excel table with the sheets containing a 'UniProt_ChEBI' column with Uniprot or ChEBI Ids in accordance with the documentation", 
                  "right"
                )
            )
          ),
shiny::column(
  width = 8,
  shiny::div(
    class = "form-group shiny-input-container",
    shiny::tags$label("\u00a0", class = "control-label"),
    shiny::div(
      style = "display: flex; justify-content: flex-end; gap: 5px; flex-wrap: wrap;",
      
      actionButton(
        width = "auto",
        inputId = ns("Demo_values"),
        label = shiny::span(
          shiny::tagList(
            tags$span("Use demo data"),
            tags$span(shiny::icon("info-circle"))
          )
        ) %>% spsComps::bsTooltip(
          "Upload the demo data",
          "right"
        )
      ),
      
      downloadButton(
        width = "auto",
        outputId = ns("Demo_download"),
        icon = NULL,
        label = shiny::span(
          shiny::tagList(
            tags$span("Download demo data"),
            tags$span(shiny::icon("cloud-arrow-down"))
          )
        ) %>% spsComps::bsTooltip(
          "Download the input demo data",
          "left"
        )
      )
    )
  )
)
        ),
        shiny::fluidRow(
          shiny::column(
            width = 6,
            textInput(
              width = "100%",
              inputId = ns("phenotype_comparison"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Condition comparison (contrasts)"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the analyzed conditions (phenotypes) contrasts. They should be the names of the xlsx sheets. Format example: CoolvsHot3,mkeyvsCool,Hot3vsCool",
                  "right"
                )
            )
          )
        )
      ),
      shinydashboard::box(
        status = "success",
        width = 12,
        title = "Input sanity checks",
        collapsible = T,
        shiny::fluidRow(
          shiny::column(
            width = 12,
            htmlOutput(ns("testcheck1"))
          ),
          shiny::column(
            width = 12,
            htmlOutput(ns("testcheck2"))
          ),
          shiny::column(
            width = 12,
            htmlOutput(ns("testcheck3"))
          ),
          shiny::column(
            width = 12,
            htmlOutput(ns("testcheck4"))
          ),
          shiny::column(
            width = 12,
            htmlOutput(ns("testcheck5"))
          ),
          shiny::column(
            width = 12,
            htmlOutput(ns("testcheck6"))
          )
        ),
          shiny::column(
            width = 12,
            htmlOutput(ns("testcheck7"))
          )
      ),
      shinydashboard::box(
        width = 12,
        title = "Network analysis parameters",
        status = "success",
        shiny::fluidRow(
          shiny::column(
            width = 5,
            sliderInput(
              width = "100%",
              inputId = ns("weight_penalty1"), 
              min = 0,
              max = 1,
              value = 0.1,
              step = 0.01,
              ticks = FALSE,
              label = shiny::span(
                shiny::tagList(
                  tags$span("Weight penalty"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Select the NetWalk framework restart probablity for weighted networks. This value will be used only if the input tables have a node weight column named 'Weight'",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 5,
            textInput(
              width = "100%",
              inputId = ns("Monet_method_stringFULL"),
              value = "--method=M1 --avgk=10 --linksdir=undirected",
              label = shiny::span(
                shiny::tagList(
                  tags$span("MONET method"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the MONET analysis method options",
                  "right"
                )
            )
          ),
		  shiny::column(
            width = 2,
            shiny::div(
			 style = "display: flex; flex-direction: column; align-items: flex-end; gap: 8px;font-size: 14px;",
				shiny::div(
				style = "display: flex; justify-content: flex-end; align-items: center; gap: 5px;",
							tags$b("Draw networks"),
							shiny::icon("info-circle")
							) %>% spsComps::bsTooltip(
													"Should Cytoscape networks be generated after MONET? A considerable increase in computation time is expected.",
													"left"
													 ),
				shinyWidgets::switchInput(
						inputId = ns("flag_run_cytoscape1"),
						value = FALSE,
						onLabel = "Yes",
						offLabel = "No",
						size = "small",
						disabled = FALSE,
						inline = FALSE
				)
			)
          )
        )
      ),
      shinydashboard::box(
        width = 12,
        title = "Knowledge databases settings",
        status = "success",
        shiny::fluidRow(
          shiny::column(
            width = 4,
            selectInput(
              width = "100%",
              inputId = ns("Use_precompiled_database1"),
              choices = c("No"="0","Yes"="1"),
              selected = '1',
              label = shiny::span(
                shiny::tagList(
                  tags$span("Use Pre-compiled Interaction file"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Choose if you would like to use a pre-compiled interaction file, formated in accordance with the instructions. A default interaction file is available. For other organisms than the pre-loaded ones, upload a pre-compiled interaction file!",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 4,
            selectInput(
              width = "100%",
              inputId = ns("CPDB_databasesFULL"),
              choices = c(
                "Reactome",
                "Wikipathways",
                "BioCarta",
                "PID",
                "NetPath",
                "HumanCyc",
                "INOH",
                "SMPDB"
              ),
              selected = "Reactome",
              multiple = TRUE,
              label = shiny::span(
                shiny::tagList(
                  tags$span("Pathways databases"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the names of the databases against which the enrichment will be performed. Example: Reactome, Wikipathways",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 4,
            selectInput(
              width = "100%",
              inputId = ns("BioMart_Dataset1"),
              selected = "hsapiens_gene_ensembl",
              multiple = F,
              choices = biomart_available_datasets,
              label = shiny::span(
                shiny::tagList(
                  tags$span("BioMart dataset"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the BioMart dataset used for the mapping of the protein Ids. For other organisms than the pre-loaded ones, provide the pre-formated interaction table! The available biomartID is only for humans. Alternatives: mmusculus_gene_ensembl, rnorvegicus_gene_ensembl, btaurus_gene_ensembl, dmelanogaster_gene_ensembl, drerio_gene_ensembl. See the documentation for the full list",
                  "right"
                )
            )
          )
        ),
        hr(),
        shiny::fluidRow(
          shiny::column(
            width = 2,
            actionButton(
              inputId = ns("AdditionalDatabases"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Add custom databases"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Upload interaction databases (Default ones are available)",
                  "right"
                )
            )
          )
        ),
        br(),
        shiny::fluidRow(
          shiny::column(
            width = 12,
            conditionalPanel(
              condition = "input.AdditionalDatabases%2 == 1 & input.Use_precompiled_database1 == '1'",
              ns = ns,
              fluidRow(
                column(
                  12,
                  span(
                    fileInput(
                      width = "100%",
                      inputId = ns("LookUp_table_file1"),
                      label = shiny::span(
                        shiny::tagList(
                          tags$span("Interaction table file"),
                          tags$span(shiny::icon("info-circle"))
                        )) %>% spsComps::bsTooltip(
                          "Upload the pre-formatted interaction database file. Leave empty for the default version, built by combining STRING v11.5, BioGRID v4.4.218, and IntAct v245",
                          "right"
                        )
                    ),
                    style = "opacity:0.5; input-group{opacity:0.5;}"
                  )
                ),
              )
            ),
            conditionalPanel(
              condition = " input.AdditionalDatabases%2 == 1  & input.Use_precompiled_database1 == '0'",
              ns = ns,
              fluidRow(
                column(
                  4,
                  span(
                    fileInput(
                      width = "100%",
                      inputId = ns("BioGRID_data_file1"),
                      label = shiny::span(
                        shiny::tagList(
                          tags$span("BioGrid database file"),
                          tags$span(shiny::icon("info-circle"))
                        )) %>% spsComps::bsTooltip(
                          "Upload the BioGrid database file. Leave empty for default version 4.4.218",
                          "right"
                        )
                    ),
                    style = "opacity:0.5; input-group{opacity:0.5;}"
                  )
                ),
                column(
                  4,
                  span(
                    fileInput(
                      width = "100%",
                      inputId = ns("STRING_data_file1"),
                      label = shiny::span(
                        shiny::tagList(
                          tags$span("STRING database file"),
                          tags$span(shiny::icon("info-circle"))
                        )) %>% spsComps::bsTooltip(
                          "Upload the STRING database file. Leave empty for default version 11.5",
                          "right"
                        )
                    ),
                    style = "opacity:0.5; input-group{opacity:0.5;}"
                  )
                ),
                column(
                  4,
                  span(
                    fileInput(
                      width = "100%",
                      inputId = ns("IntAct_data_file1"),
                      label = shiny::span(
                        shiny::tagList(
                          tags$span("IntAct database file"),
                          tags$span(shiny::icon("info-circle"))
                        )) %>% spsComps::bsTooltip(
                          "Upload the IntAct database file. Leave empty for default release 245",
                          "right"
                        )
                    ),
                    style = "opacity:0.5; input-group{opacity:0.5;}"
                  )
                )
              )
            )
          )
        )
      ),
      shinydashboard::box(
        width = 12,
        status = "success",
        shiny::fluidRow(
          shiny::column(
            width = 6,
            textInput(
              width = "100%",
              inputId = ns("Client_email"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Email address (Optional)"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "If you would like to be notified when the analysis is completed, you can provide an email address.",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 6,
            shiny::div(
              class = "form-group shiny-input-container",
              shiny::tags$label("\u00a0", class="control-label"),
              shiny::div(
                actionButton(
                  width = "100%",
                  inputId = ns("submit_Full_analysis"),
                  label = shiny::span(
                    shiny::tagList(
                      tags$span("Submit"),
                      tags$span(shiny::icon("info-circle"))
                    )) %>% spsComps::bsTooltip(
                      "Run the full multi-omics integrative analysis.",
                      "right"
                    )
                )
              )
            )
          )
        )
      )
    ),
    conditionalPanel(
      condition = "input.Type_of_analysis == 'Custom algorithms'",
      ns = ns,
      shinydashboard::box(
        width = 12,
        status = "success",
        shiny::fluidRow(
          shiny::column(
            width = 12,
            textInput(
              inputId = ns("Results_index"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Results directory index"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "What is the directory index where your NOODAI analysis results are found? The index is provided after the input data is submitted for the analysis",
                  "right"
                )
            ),
            shiny::strong("Note: this field must be provided before running any of the sections below!")
          )
        )
      ),
      shinydashboard::box(
        width = 12,
        status = "success",
        #Start First Section
        shiny::fluidRow(
          shiny::column(
            width = 4,
            textInput(
              width = "100%",
              inputId = ns("phenotype_comparison1"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Condition comparison (contrasts)"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the analyzed conditions (phenotypes) contrasts. They should be the names of the xlsx sheets. Format example: CoolvsHot3,mkeyvsCool,Hot3vsCool",
                  "right"
                )
            )
          ),
          shiny::column(
            4,
            sliderInput(
              width = "100%",
              inputId = ns("weight_penalty"),
              min = 0,
              max = 1,
              value = 0.1,
              step = 0.01,
              ticks = FALSE,
              label = shiny::span(
                shiny::tagList(
                  tags$span("Weight penalty"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Select the NetWalk framework restart probablity for weighted networks. This value will be used only if the input tables have a node weight column named 'Weight'",
                  "right"
                )
            )
          ),
		  shiny::column(
            width = 4,
            shiny::div(
			 style = "display: flex; flex-direction: column; align-items: flex-end; gap: 8px;font-size: 14px;",
				shiny::div(
				style = "display: flex; justify-content: flex-end; align-items: center; gap: 5px;",
							tags$b("Draw networks"),
							shiny::icon("info-circle")
							) %>% spsComps::bsTooltip(
													"Should Cytoscape networks be generated after MONET? A considerable increase in computation time is expected.",
													"left"
													 ),
				shinyWidgets::switchInput(
						inputId = ns("flag_run_cytoscape"),
						value = FALSE,
						onLabel = "Yes",
						offLabel = "No",
						size = "small",
						disabled = FALSE,
						inline = FALSE
				)
			)
          )
        ),
        shiny::fluidRow(
          shiny::column(
            width = 4,
            selectInput(
              width = "100%",
              inputId = ns("Use_precompiled_database"),
              choices = c("No"="0","Yes"="1"),
              selected ='1',
              label = shiny::span(
                shiny::tagList(
                  tags$span("Use Pre-compiled Interaction file"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Choose if you would like to use a pre-compiled interaction file, formated in accordance with the instructions. A default interaction file is available. For other organisms than the pre-loaded ones, upload a pre-compiled interaction file!",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 4,
            textInput(
              width = "100%",
              inputId = ns("splicing_file_name"),
              value = "DTU",
              label = shiny::span(
                shiny::tagList(
                  tags$span("DTU file"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the name of the splicing file in the archive. If you use splicing data, the contrast must be symmetric! Leave default for no DTU data",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 4,
            selectInput(
              width = "100%",
              inputId = ns("BioMart_Dataset"),
              selected = "hsapiens_gene_ensembl",
              multiple = F,
              choices = biomart_available_datasets,
              label = shiny::span(
                shiny::tagList(
                  tags$span("BioMart dataset"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the BioMart dataset used for the mapping of the protein Ids. For other organisms than the pre-loaded ones, provide the pre-formated interaction table! The available biomartID is only for humans. Alternatives: mmusculus_gene_ensembl, rnorvegicus_gene_ensembl, btaurus_gene_ensembl, dmelanogaster_gene_ensembl, drerio_gene_ensembl. See the documentation for the full list",
                  "right"
                )
            )
          )
        ),
        conditionalPanel(
          condition = "input.Use_precompiled_database == '0'",
          ns = ns,
          shiny::fluidRow(
            shiny::column(
              width = 3,
              fileInput(
                width = "100%",
                inputId = ns("BioGRID_data_file"),
                label = shiny::span(
                shiny::tagList(
                  tags$span("BioGrid database file"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Upload the BioGrid database file. Leave empty for default version 4.4.218",
                  "right"
                )
              )
            ),
            shiny::column(
              width = 3,
              fileInput(
                width = "100%",
                inputId = ns("STRING_data_file"),
                label = shiny::span(
                shiny::tagList(
                  tags$span("STRING database file"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Upload the STRING database file. Leave empty for default version 11.5",
                  "right"
                )
              ),
            ),
            shiny::column(
              width = 3,
              fileInput(
                width = "100%",
                inputId = ns("IntAct_data_file"),
                label = shiny::span(
                shiny::tagList(
                  tags$span("IntAct database file"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Upload the IntAct database file. Leave empty for default release 245",
                  "right"
                )
              )
            ),
            shiny::column(
              width = 3,
              fileInput(
                width = "100%",
                inputId = ns("file_DEA_names1"),
                label = shiny::span(
                shiny::tagList(
                  tags$span("Omics files archive"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Upload an archive containing at least one omics file. Each omics file needs to be an Excel table with the sheets containing a 'UniProt_ChEBI' column with Uniprot or ChEBI Ids in accordance with the documentation",
                  "right"
                )
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.Use_precompiled_database == '1'",
          ns = ns,
          shiny::fluidRow(
            shiny::column(
              width = 6,
              fileInput(
                width = "100%",
                inputId = ns("LookUp_table_file"),
                label = shiny::span(
                shiny::tagList(
                  tags$span("Interaction table file"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Upload the pre-formatted interaction database file. Leave empty for the default version, built by combining STRING v11.5, BioGRID v4.4.218, and IntAct v245",
                  "right"
                )
              )
            ),
            shiny::column(
              width = 6,
              fileInput(
                width = "100%",
                inputId = ns("file_DEA_names2"),
                label = shiny::span(
                shiny::tagList(
                  tags$span("Omics files archive"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Upload an archive containing at least one omics file. Each omics file needs to be an Excel table with the sheets containing a 'UniProt_ChEBI' column with Uniprot or ChEBI Ids in accordance with the documentation",
                  "right"
                )
              )
            )
          )
        ),
        shiny::fluidRow(
          shiny::column(
            offset = 1,
            10,
            actionButton(
              width = "100%",
              inputId = ns("submit_Network_analysis"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Submit"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Run the main network analysis algorithm extracting the top central nodes.",
                  "right"
                )
            )
          )
        ),
        ###End First Section  
      ),
      shinydashboard::box(
        width = 12,
        status = "success",
        ###Start Second Section
        shiny::fluidRow(
          shiny::column(
            width = 3,
            textInput(
              width = "100%",
              inputId = ns("edge_file_path"),
              value = "edge_files_PPINetworks/Symbol",
              label = shiny::span(
                shiny::tagList(
                  tags$span("Edge file path"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the path where the network edges are saved.",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 3,
            textInput(
              width = "100%",
              inputId = ns("Monet_method_string"),
              value = "--method=M1 --avgk=10 --linksdir=undirected",
              label = shiny::span(
                shiny::tagList(
                  tags$span("MONET method"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the MONET analysis method options",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 3,
            textInput(
              width = "100%",
              inputId = ns("tmp_bin_folder"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Temporary folder"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide a temporary bin server. Leave empty for default",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 3,
            textInput(
              width = "100%",
              inputId = ns("monet_path"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("MONET path"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the monet executable path. Developer mode",
                  "right"
                )
            )
          )
        ),
        shiny::fluidRow(
          shiny::column(
            width = 10,
            offset = 1,
            actionButton(
              width = "100%",
              inputId = ns("submit_MONET_analysis"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Submit"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Run the MONET decomposition analysis.",
                  "right"
                )
            )
          )
        ),
        ###End Second Section
      ),
      shinydashboard::box(
        width = 12,
        status = "success",
        ###Start third section
        shiny::fluidRow(
          shiny::column(
            4,
            selectInput(
              width = "100%",
              inputId = ns("CPDB_databases"),
              choices = c(
                "Reactome",
                "Wikipathways",
                "BioCarta",
                "PID",
                "NetPath",
                "HumanCyc",
                "INOH",
                "SMPDB"
              ),
              selected = "Reactome",
              multiple = TRUE,
              label = shiny::span(
                shiny::tagList(
                  tags$span("CPDB databases"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the names of the databases against which the enrichment is to be performed. Example: Reactome, Wikipathways",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 4,
            fileInput(
              width = "100%",
              inputId = ns("MONET_background_file"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("MONET background file"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the MONET background file. Leave empty for default",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 4,
            fileInput(
              width = "100%",
              inputId = ns("CPDB_database_file"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("CPDB database file"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the CPDB database file. Leave empty for default",
                  "right"
                )
            ),
          )
        ),
        shiny::fluidRow(
          shiny::column(
            width = 10,
            offset = 1,
            actionButton(
              width = "100%",
              inputId = ns("submit_MONET_pathways"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Submit"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Extract the pathways associated with each network submodule.",
                  "right"
                )
            )
          )
        ),
        ###End Third Section
      ),
      shinydashboard::box(
        width = 12,
        status = "success",
        ###Start Fourth Section
        shiny::fluidRow(
          shiny::column(
            width = 3,
            textInput(
              width = "100%",
              inputId = ns("files_edges_path"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Edge files directory"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the edges' directory. Leave empty for default",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 3,
            fileInput(
              width = "100%",
              inputId = ns("centralities_file"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Centralities file"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the network centralities file. Leave empty for default",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 3,
            fileInput(
              width = "100%",
              inputId = ns("Kinome_database"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Kinome Dataset"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the Kinome dataset file. Leave empty for default",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 3,
            textInput(
              width = "100%",
              inputId = ns("file_extension"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("File ending"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Provide the file ending specific for a type of analysis. Default to 'Total'. Leave empty for default",
                  "right"
                )
            )
          )
        ),
        shiny::fluidRow(
          shiny::column(
            width = 10,
            offset = 1,
            actionButton(
              width = "100%",
              inputId = ns("submit_Circos_and_auxiliary"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Submit"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "Create the summary plots and report.",
                  "right"
                )
            )
          )
        )
        ###End Fourth Section
      )
    ),
    conditionalPanel(
      condition = "input.Type_of_analysis == 'Results download'",
      ns = ns,
      shinydashboard::box(
        status = "success",
        width = 12,
        shiny::fluidRow(
          shiny::column(
            width = 4,
            textAreaInput(
              width = "100%",
              inputId = ns("Results_index_save"),
              label = shiny::span(
                shiny::tagList(
                  tags$span("Results directory index"),
                  tags$span(shiny::icon("info-circle"))
                )) %>% spsComps::bsTooltip(
                  "The directory where the NOODAI analysis results are found. The index is provided after the input data is submitted for the analysis",
                  "right"
                )
            )
          ),
          shiny::column(
            width = 4,
            div(
              width = "100%",
              style = "margin-top: 35px;",
              downloadButton(
                outputId = ns("downloadData"),
                label = "Save"
              ),
              span(
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                title = "Download the results",
                icon("info-circle")
              )
            )
          ),
          shiny::column(
            width = 2,
            offset = 2,
            div(
              width = "100%",
              style = "margin-top: 35px;",
              downloadButton(
                outputId = ns("Demo_Reults_download"),
                icon = NULL,
                label = shiny::span(
                  shiny::tagList(
                    tags$span("Demo Results"),
                    tags$span(shiny::icon("cloud-arrow-down"))
                  )) %>% spsComps::bsTooltip(
                    "Download the pre-compiled results for the Demo input data.",
                    "right"
                  )
              )
            )
          )
        )
      )
    ),
    conditionalPanel(
      condition = "output.showRes",
      ns = ns,
      shiny::fluidRow(
        shiny::column(
          12,
          align = "center",
          textOutput(ns("text_output_string"))
        ),
        shiny::column(
          12,
          align = "center",
          span(
            htmlOutput(ns("verb")),
            style = "font-weight:bold;"
          )
        ),
        shiny::column(
          12,
          align = "center",
          textOutput(ns("text_output_string2"))
        ),
        shiny::column(
          12,
          align = "center",
          textOutput(ns("text_output_string3"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showResCentralities",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          htmlOutput(ns("CentralitiesOut1"))
        ),
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("CentralitiesOut2"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showResMONET",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          htmlOutput(ns("MONETOut1"))
        ),
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("MONETOut2"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showResPathways",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          htmlOutput(ns("PathwaysOut1"))
        ),
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("PathwaysOut2"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showResImages",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          htmlOutput(ns("ImagesOut1"))
        ),
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("ImagesOut2"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showErrorFormat",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("showErrorFormat1"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showErrorDataset",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("showErrorDataset1"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showErrorSample",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("showErrorSample1"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showErrorContrasts",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("showErrorContrasts1"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showErrorBioMart",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("showErrorBioMart1"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showErrorMonet",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("showErrorMonet1"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showErrorCPDB",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("showErrorCPDB1"))
        )
      )
    ),
    conditionalPanel(
      condition = "output.showErrorBiomartServer",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("showErrorBiomartServer1"))
        )
      )
    ),
	conditionalPanel(
      condition = "output.showWarningInput",
      ns = ns,
      fluidRow(
        column(
          12,
          align = "center",
          verbatimTextOutput(ns("showWarningInput1"))
        )
      )
    ),
    br(),
    fluidRow(
      column(
        offset = 10,
        2,
        span(
          textOutput(ns("text_output_version")),
          style = "color:gray; float:right;"
        )
      )
    )
  )
}
    
mod_NOODAI_pipeline_server <- function(input, output, session) {

      ns <- session$ns

      output$text_output_version <- renderText({"version 2.0.0"})
  
      options(shiny.maxRequestSize=5*2^30)
      
      showVal <- reactiveVal(TRUE)
      showRes1 <- reactiveVal(FALSE)
      showResCentralities1 <- reactiveVal(FALSE)
      showResMONET1 <- reactiveVal(FALSE)
      showResPathways1 <- reactiveVal(FALSE)
      showResImages1 <- reactiveVal(FALSE)
      showMemFull <- reactiveVal(FALSE)
      showErrorFormatMessage <- reactiveVal(FALSE)
      showErrorDatasetMessage <- reactiveVal(FALSE)
      showErrorSampleMessage <- reactiveVal(FALSE)
      showErrorContrastsMessage <- reactiveVal(FALSE)
      showErrorBioMartMessage <- reactiveVal(FALSE)
      showErrorMonetMessage <- reactiveVal(FALSE)
      showErrorCPDBMessage <- reactiveVal(FALSE)
      showErrorBiomartServerMessage <- reactiveVal(FALSE)
	  showWarningLowInput <- reactiveVal(FALSE)
      file_DEA_names5 <- reactiveVal(value=NULL)
	  special_error_messageDataset <- reactiveVal(value = c())
      special_error_messageContrast <- reactiveVal(value = c())
      special_error_messageSample <- reactiveVal(value = c())
  
      #####Start Network analysis#######
  
      Results_index <- reactive({
        x <- input$Results_index
        if(x == "") {

          today <- as.character(Sys.Date())
          x <- paste(today, 1, sep = "_")
          dirs <- dir(working_dir)

          if(length((grep(x, dirs)[1])) > 0) {

            dirs_today <- dirs[grepl(today, dirs)]
            n <- length(dirs_today)
            x <- paste(x, n + 1, paste0(sample(c(0:9, LETTERS), 6, T), collapse=''), sep = "_")

          }
          
        }

        return(x)
      })

      Use_precompiled_database <- reactive({
        x <- input$Use_precompiled_database
        return(as.numeric(x))
      })
  
  
      BioMart_Dataset1 <- reactive({
        x <- input$BioMart_Dataset1
        return(x)
      })
  
      Client_email <- reactive({
        x <- input$Client_email
        x <- gsub(" ", "", x, fixed = TRUE)
        return(x)
      })
      
      BioMart_Dataset <- reactive({
        x <- input$BioMart_Dataset
        return(x)
      })
  
      Use_precompiled_database1 <- reactive({
        x <- input$Use_precompiled_database1
        return(as.numeric(x))
      })
  
      Type_of_analysis <- reactive({
        x <- input$Type_of_analysis
        return(x)
      })
  
  
      AdditionalDatabases <- reactive({
        x <- input$AdditionalDatabases
        return(as.numeric(x))
      })
  
      weight_penalty <- reactive({
        x <- input$weight_penalty
        return(as.numeric(x))
      })
	  
  
      weight_penalty1 <- reactive({
        x <- input$weight_penalty1
        return(as.numeric(x))
      })
	  
	  flag_run_cytoscape1 <- reactive({
        return(as.numeric(isTRUE(input$flag_run_cytoscape1)))
      })
	  
	  flag_run_cytoscape <- reactive({
       return(as.numeric(isTRUE(input$flag_run_cytoscape)))
      })
      
      BioGRID_data_file <- reactive({
        if(length(input$BioGRID_data_file) > 0) {
          x <- input$BioGRID_data_file$datapath
          return(x)
        }
      })
      
      BioGRID_data_file1 <- reactive({
        if(length(input$BioGRID_data_file1) > 0) {
          x <- input$BioGRID_data_file1$datapath
          return(x)
        }
      })
      
      STRING_data_file <- reactive({
        if(length(input$STRING_data_file) > 0) {
          x <- input$STRING_data_file$datapath
          return(x)
        }
      })
      
      STRING_data_file1 <- reactive({
        if(length(input$STRING_data_file1) > 0) {
          x <- input$STRING_data_file1$datapath
          return(x)
        }
      })
      
      IntAct_data_file <- reactive({
        if(length(input$IntAct_data_file) > 0) {
          x <- input$IntAct_data_file$datapath
          return(x)
        }
      })
      
      IntAct_data_file1 <- reactive({
        if(length(input$IntAct_data_file1) > 0) {
          x <- input$IntAct_data_file1$datapath
          return(x)
        }
      })
  
      file_DEA_names1 <- reactive({
        if(length(input$file_DEA_names1) > 0) {

          if(length(working_dir) != 0) {

            tmp_file_loc <- paste0(
              working_dir,
              "/Results_",
              Results_index(),
              "/tmp_FilesIni"
            )

            print(tmp_file_loc)

            unlink("tmp_file_loc", recursive = TRUE)
            dir.create(tmp_file_loc, recursive = T)

            x <- archive::archive_extract(input$file_DEA_names1$datapath, dir = tmp_file_loc)
            return(file.path(tmp_file_loc, x))
          }
        }
      })
  
      file_DEA_names2 <- reactive({
        if(length(input$file_DEA_names2) > 0){

          if(length(working_dir) != 0) {
            
            tmp_file_loc <- paste0(
              working_dir,
              "/Results_",
              Results_index(),
              "/tmp_FilesIni"
            )

            unlink("tmp_file_loc", recursive = TRUE)
            dir.create(tmp_file_loc, recursive = T)

            x <- archive::archive_extract(input$file_DEA_names2$datapath, dir = tmp_file_loc)
            return(file.path(tmp_file_loc, x))

          }
        }
      })
      
      
      file_DEA_names3 <- reactive({
        if(input$Demo_values != 1) {

          if(length(input$file_DEA_names3) > 0) {
            
            if(length(working_dir) != 0) {
              
              tmp_file_loc <- paste0(
                working_dir,
                "/Results_",
                Results_index(),
                "/tmp_FilesIni"
              )

              unlink("tmp_file_loc", recursive=TRUE)
              dir.create(tmp_file_loc, recursive = T)

              x <- archive::archive_extract(input$file_DEA_names3$datapath, dir = tmp_file_loc)
              files <- list.files(pattern = ".xlsx", path = tmp_file_loc, recursive = T, full.names = T)
              # return(file.path(tmp_file_loc, x))
              return(files)
            }
          }
        }
      })
      
      LookUp_table_file <- reactive({
        if(length(input$LookUp_table_file) > 0){
          x <- input$LookUp_table_file$datapath
          return(x)
        }
      })
      
      LookUp_table_file1 <- reactive({
        if(length(input$LookUp_table_file1) > 0) {
          x <- input$LookUp_table_file1$datapath
          return(x)
        }
      })
      
      phenotype_comparison <- reactive({
        if(length(input$phenotype_comparison) > 0) {
          x <- input$phenotype_comparison
          x <- gsub(" ", "", x, fixed = TRUE)
          x <- strsplit(x,',')
          x <- unlist(x)
          return(x)
        }
      })
      
      phenotype_comparison1 <- reactive({
        if(length(input$phenotype_comparison1) > 0) {
          x <- input$phenotype_comparison1
          x <- gsub(" ", "", x, fixed = TRUE)
          x <- strsplit(x,',')
          x <- unlist(x)
          return(x)
        }
      })
      
      splicing_file_name <- reactive({
        if(length(input$splicing_file_name) > 0) {
          x <- input$splicing_file_name
          x <- gsub(" ", "", x, fixed = TRUE)
          return(x)
        }
      })

      ## SANITY CHECKS
      output$testcheck1 = renderText({
        out_text = "Condition contrasts are provided."
        if (length(phenotype_comparison()) == 0) {
          as.character(paste(shiny::icon("xmark"), out_text, sep = " "))
        } else {
          as.character(paste(shiny::icon("check"), out_text, sep = " "))
        }
      })

      output$testcheck3 = renderText({
        out_text = "There are at least 2 condition names."
        if (length(unique(unlist(strsplit(phenotype_comparison(),"vs")))) < 2) {
          as.character(paste(shiny::icon("xmark"), out_text, sep = " "))
        } else {
          as.character(paste(shiny::icon("check"), out_text, sep = " "))
        }
      })

      output$testcheck4 = renderText({
        out_text = "At least 1 omics file is provided."
        if (length(file_DEA_names3()) == 0) {
          as.character(paste(shiny::icon("xmark"), out_text, sep = " "))
        } else {
          as.character(paste(shiny::icon("check"), out_text, sep = " "))
        }
      })

      output$testcheck5 = renderText({
        out_text = "Identifiers are found in the 'UniProt_ChEBI' column."
        checks = c()
        if (length(file_DEA_names3()) > 0) {
          
          for (file in file_DEA_names3()) {
            
            file_content = file %>%
              readxl::excel_sheets() %>%
              set_names() %>%
              map(readxl::read_excel, path = file)

            check_i = any(unlist(lapply(file_content, colnames)) == "UniProt_ChEBI")

            checks = c(checks, check_i)

          }
        
        } else {
          checks = F
        }

        if (!all(checks)) {
          as.character(paste(shiny::icon("xmark"), out_text, sep = " "))
        } else {
          as.character(paste(shiny::icon("check"), out_text, sep = " "))
        }
        
      })

      output$testcheck6 = renderText({
        out_text = "Sheet names coincide with the provided phenotype comparisons."
        checks = c()
        if (length(file_DEA_names3()) > 0) {

          if (length(phenotype_comparison()) > 0) {

            for (file in file_DEA_names3()) {

              sheet_names = file %>%
                readxl::excel_sheets()

              check_i = all(sheet_names %in% phenotype_comparison())
              checks = c(checks, check_i)

            }

          } else {
            checks = FALSE
          }

        } else {
          checks = FALSE
        }

        if (!all(checks)) {
          as.character(paste(shiny::icon("xmark"), out_text, sep = " "))
        } else {
          as.character(paste(shiny::icon("check"), out_text, sep = " "))
        }
      })

      output$testcheck7 = renderText({
        out = "Conditions are: "
        if (length(phenotype_comparison()) > 0) {
          phenotype_names = unique(unlist(strsplit(phenotype_comparison(), "vs")))
          out = paste0("Conditions are: ", paste(phenotype_names, collapse = ", "))
        }
        return(out)
      })

      ##
      
      observeEvent(input$Type_of_analysis, {
        
        if(Type_of_analysis() == 'Run the full pipeline')
        {
          shinyjs::show("Demo_values")
          shinyjs::show("Demo_download")
        }else{
          shinyjs::hide("Demo_values")
          shinyjs::hide("Demo_download")
        }
    
      })
  
      observeEvent(input$Demo_values, {
        
        updateTextInput(session, "phenotype_names", value="M1,M2a,M2c")
        updateTextInput(session, "phenotype_comparison", value="M1vsM2a,M2avsM1,M1vsM2c,M2cvsM1,M2avsM2c,M2cvsM2a")
        updateTextInput(session, "splicing_file_name", value="DTU")
        shinyjs::hide("file_DEA_names3")

        output$testcheck4 = renderText({
          as.character(paste(shiny::icon("check"), "At least 1 omics file is provided.", sep = " "))
        })
        output$testcheck5 = renderText({
          as.character(paste(shiny::icon("check"), "Identifiers are found in the 'UniProt_ChEBI' column.", sep = " "))
        })
        output$testcheck6 = renderText({
          as.character(paste(shiny::icon("check"), "Sheet names coincide with the provided phenotype comparisons", sep = " "))
        })

        if(length(working_dir) != 0) {
          dir.create(paste0(working_dir,"/Results_",Results_index()))
          tmp_file_loc <- paste0(working_dir,"/Results_",Results_index(),"/tmp_FilesIni")
          unlink("tmp_file_loc", recursive = TRUE)        
          dir.create(tmp_file_loc, recursive = T)
          x <- archive::archive_extract(paste0(working_dir,"/Data_Omics_Demo.7z"), dir = tmp_file_loc)
          file_DEA_names5(paste0(tmp_file_loc,"/",x))
        }
        
      })
  
      observeEvent(input$submit_Network_analysis, {

        cat(file = stderr(), "Running network analysis\n")
        
        shinyjs::addClass(id = "submit_Network_analysis", class = "loading dots")
        shinyjs::disable("submit_Full_analysis")
        shinyjs::disable("submit_Network_analysis")
        shinyjs::disable("submit_MONET_pathways")
        shinyjs::disable("submit_Circos_and_auxiliary")
        shinyjs::disable("submit_MONET_analysis")
        
        
        if(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE)) < 300000) {
          showMemFull(TRUE)
        }else{
          
          if(length(working_dir)!=0){
            
            if(biomart_check() == 0) {

              showVal(FALSE)
              showErrorBiomartServerMessage(TRUE)

            } else {
              
              dir.create(paste0(working_dir,"/Results_",Results_index()), recursive = T)

              token <- UUIDgenerate()
              message(paste0("running task for token: ", token))
              
              if(is.null(jobs[[token]])) {

                if(Type_of_analysis() == 'Run the full pipeline') {

                  file_DEA_names <- file_DEA_names3()
                  Use_precompiled_databaseVal <- Use_precompiled_database1()
                  LookUp_table_fileVal <- LookUp_table_file1()
                  BioGRID_data_fileVal <- BioGRID_data_file1()
                  STRING_data_fileVal <- STRING_data_file1()
                  IntAct_data_fileVal <- IntAct_data_file1()
                  phenotype_names <- unique(unlist(strsplit(phenotype_comparison(), "vs")))
                  phenotype_comparison <- phenotype_comparison()
                  BioMart_Dataset <- BioMart_Dataset1()

                } else {
                  
                  Use_precompiled_databaseVal <- Use_precompiled_database()
                  LookUp_table_fileVal <- LookUp_table_file()
                  BioGRID_data_fileVal <- BioGRID_data_file()
                  STRING_data_fileVal <- STRING_data_file()
                  IntAct_data_fileVal <- IntAct_data_file()
                  phenotype_names <- unique(unlist(strsplit(phenotype_comparison1(), "vs")))
                  phenotype_comparison <- phenotype_comparison1()
                  BioMart_Dataset <- BioMart_Dataset()
                  if(Use_precompiled_database() == '0'){file_DEA_names <- file_DEA_names1()}
                  if(Use_precompiled_database() == '1'){file_DEA_names <- file_DEA_names2()}

                }
                
                
                analyses_flag <- 1
                
                
                for (i in 1:length(phenotype_comparison)){
                     if(length(unlist(sapply(phenotype_names,FUN=function(x){grep(x,phenotype_comparison[i])})))!=2){
                        analyses_flag<-0
                        showVal(FALSE)
                        showErrorContrastsMessage(TRUE)
                        showErrorSampleMessage(TRUE)
                        special_error_messageContrast("Cause: There is a mismatch between the conditions names and contrasts!")
                        special_error_messageSample("Cause: There is a mismatch between the conditions names and contrasts!")
                        break()
                    }
                }
                
                if(length(file_DEA_names)==0){
                    analyses_flag<-0
                    showVal(FALSE)
                    showErrorDatasetMessage(TRUE)
                    special_error_messageDataset("Cause: Either the archive is empty or wrongly formatted!")
                }
                if(length(phenotype_comparison)==0){
                    analyses_flag<-0
                    showVal(FALSE)
                    showErrorContrastsMessage(TRUE)
                    special_error_messageContrast("Cause: No contrasts were correctly provided!")
                }
                if(length(grep("vs", phenotype_comparison))!=length(phenotype_comparison)){
                    analyses_flag<-0
                    showVal(FALSE)
                    showErrorContrastsMessage(TRUE)
                    special_error_messageContrast("Cause: 'vs' was not used to separate the conditions in the contrasts!")
                }
                if(length(phenotype_names)<2){
                    analyses_flag<-0
                    showVal(FALSE)
                    showErrorSampleMessage(TRUE)
                    special_error_messageSample("Cause: There are fewer than two condition names!")
                }
                
				
				
				                    if(analyses_flag){
									contrast_uniprot_counts <- setNames(vector("list", length(phenotype_comparison)), phenotype_comparison)
                                      for (j in 1:length(file_DEA_names)){
                                        file_Sheets1 <- readxl::excel_sheets(file_DEA_names[j])
                                        if(length(which(phenotype_comparison %in% file_Sheets1))!=length(phenotype_comparison)){
                                          analyses_flag<-0
                                          showVal(FALSE)
                                          showErrorContrastsMessage(TRUE)
                                          showErrorDatasetMessage(TRUE)
                                          special_error_messageDataset("Cause: The provided contrasts are not found in the Excel tables!")
                                          special_error_messageContrast("Cause: The provided contrasts are not found in the Excel tables!")
                                          break()
                                        }
                                        for (i in 1:length(phenotype_comparison)){
                                          aux <- as.data.frame(readxl::read_excel(file_DEA_names[j],sheet = phenotype_comparison[i]))
                                          #Get the corresponding entrez ids from the uniprot ones#############
                                          if(length(aux)>0){
                                            if(length(which(colnames(aux) %in% "UniProt_ChEBI"))==0){
                                              showVal(FALSE)
                                              showErrorFormatMessage(TRUE)
                                              analyses_flag<-0
                                              break()
                                            }
                                            aux1 <- aux$UniProt_ChEBI[which(aux$UniProt_ChEBI %in% uniprot_ids_maps$UniProt_ChEBI)]
                                            if(length(aux1)<(length(aux$UniProt_ChEBI)/3)){
                                              showVal(FALSE)
                                              showErrorFormatMessage(TRUE)
                                              analyses_flag<-0
                                              break()
                                            }
											 contrast_uniprot_counts[[phenotype_comparison[i]]] <- unique(c(contrast_uniprot_counts[[phenotype_comparison[i]]], aux1))
                                          }
                                        }
                                      }
									  low_input_flag <- any(sapply(contrast_uniprot_counts, function(x) length(x) < 150))
									  if (low_input_flag) {
											showWarningLowInput(TRUE)
											}
                                    }
				
				
				
                
                
                if(length(which(BioMart_Dataset %in% biomart_available_datasets$Dataset)) == 0) {
                  
                  showVal(FALSE)
                  showErrorBioMartMessage(TRUE)
                  analyses_flag <- 0

                }
                
                
                if(analyses_flag) {

                  source(
                    paste0(working_dir,"/Integrative_analysis.R"),
                    local = TRUE,
                    echo = TRUE,
                    max.deparse.length = 100000
                  )

                  jobs[[token]] <- callr::r_bg(
                    Integrative_Network_analysis,
                    supervise = FALSE,
                    stdout = paste0(
                      working_dir,
                      "/Results_",
                      Results_index(),
                      "/out_Centralities.txt"
                    ),
                    package = TRUE,
                    stderr = paste0(
                      working_dir,
                      "/Results_",
                      Results_index(),
                      "/log_Centralities.txt"
                    ),
                    args = list(
                      working_dir = working_dir,
                      BioGRID_data_file = BioGRID_data_fileVal,
                      STRING_data_file = STRING_data_fileVal,
                      IntAct_data_file = IntAct_data_fileVal,
                      file_DEA_names = file_DEA_names,
                      phenotype_names = phenotype_names,
                      phenotype_comparison = phenotype_comparison,
                      splicing_file_name = splicing_file_name(),
                      Use_precompiled_database = Use_precompiled_databaseVal,
                      LookUp_table_file = LookUp_table_fileVal,
                      Results_index = Results_index(),
                      BioMart_Dataset = BioMart_Dataset,
                      weight_penalty = weight_penalty()
                    )
                  )

                  pars_ls = list(
                    NOODAI_module = "Network Analysis",
                    BioGRID_data_file = BioGRID_data_fileVal,
                      STRING_data_file = STRING_data_fileVal,
                      IntAct_data_file = IntAct_data_fileVal,
                      file_DEA_names = file_DEA_names,
                      phenotype_names = phenotype_names,
                      phenotype_comparison = phenotype_comparison,
                      splicing_file_name = splicing_file_name(),
                      Use_precompiled_database = Use_precompiled_databaseVal,
                      LookUp_table_file = LookUp_table_fileVal,
                      Results_index = Results_index(),
                      BioMart_Dataset = BioMart_Dataset,
                      weight_penalty = weight_penalty()
                  )

                  save_run_info(
                    pars_ls,
                    out_dir = file.path(
                      working_dir,
                      paste0("Results_", Results_index())
                    )
                  )
                  
                  showVal(FALSE)
                  showResCentralities1(TRUE)
                }
              }
              
            }
          }
        }

        cat(file = stderr(), "Finished running network analysis\n")
        
      })
  
      # #####End Network analysis#######
      
      # #####Start MONET analysis#######
      
      edge_file_path <- reactive({
        if(length(input$edge_file_path)>0){
          x <- input$edge_file_path
          if(x==""){
            return(NULL)
          }else{
            return(x)
          }
        }
      })
      
      
      Monet_method_string <- reactive({
        if(length(input$Monet_method_string)>0){
          x <- input$Monet_method_string
          return(x)
        }
      })
      
      Monet_method_stringFULL <- reactive({
        if(length(input$Monet_method_stringFULL)>0){
          x <- input$Monet_method_stringFULL
          return(x)
        }
      })
      
      tmp_bin_folder <- reactive({
        if(length(input$tmp_bin_folder)>0){
          x <- input$tmp_bin_folder
          if(x==""){
            return(NULL)
          }else{
            return(x)
          }
        }
      })
      
      
      monet_path <- reactive({
        if(length(input$monet_path)>0){
          x <- input$monet_path
          if(x==""){
            return(NULL)
          }else{
            return(x)
          }
        }
      })
      
      
      observeEvent(input$submit_MONET_analysis, {
        
        shinyjs::addClass(id = "submit_MONET_analysis", class = "loading dots")
        shinyjs::disable("submit_Full_analysis")
        shinyjs::disable("submit_Network_analysis")
        shinyjs::disable("submit_MONET_pathways")
        shinyjs::disable("submit_Circos_and_auxiliary")
        shinyjs::disable("submit_MONET_analysis")
        
        if(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE)) < 300000){
          showMemFull(TRUE)
        }else{
          
          if(length(working_dir)!=0){
            dir.create(paste0(working_dir,"/Results_",Results_index()))
            #con <- file(paste0(working_dir,"/Results_",Results_index(),"/log.txt"))
            #sink(paste0(working_dir,"/Results_",Results_index(),"/log.txt"), append=TRUE,type = c("output", "message"))
            
            
            token <- UUIDgenerate()
            message(paste0("running task for token: ", token))
            
            analyses_flag <- 1
            
            if(length(grep(Monet_exp,Monet_method_string()))==0){
              showVal(FALSE)
              showErrorMonetMessage(TRUE)
              analyses_flag<-0
            }
            
            
            if(analyses_flag){
              if(is.null(jobs[[token]])){
                source(paste0(working_dir,"/Integrative_analysis.R"), local = TRUE,echo=TRUE, max.deparse.length=100000)
                # jobs[[token]] <<- callr::r_bg(
                jobs[[token]] <- callr::r_bg(
                  MONET_analysis,
                  supervise = FALSE,
                  stdout = paste0(
                    working_dir,
                    "/Results_",
                    Results_index(),
                    "/out_MONET.txt"
                  ),
                  package = TRUE,
                  stderr = paste0(
                    working_dir, 
                    "/Results_",
                    Results_index(),
                    "/log_MONET.txt"
                  ),
                  args = list(
                    working_dir = working_dir,
                    edge_file_path = edge_file_path(),
                    monet_path = monet_path(),
                    Monet_method_string = Monet_method_string(),
                    tmp_bin_folder = tmp_bin_folder(),
                    Results_index = Results_index()
                  )
                )
                
                pars_ls = list(
                    NOODAI_module = "Monet",
                    edge_file_path = edge_file_path(),
                    Monet_method_string = Monet_method_string(),
                    Results_index = Results_index()
                )

                save_run_info(
                  pars_ls,
                  out_dir = file.path(
                    working_dir,
                    paste0("Results_", Results_index())
                  )
                )
                
                showVal(FALSE)
                showResMONET1(TRUE)
              }
            }
            
            
          }
        }
        
      })
      
      
      # #####End MONET analysis#######
      
      
      # #####Start MONET pathways analysis#######
      
      CPDB_databases <- reactive({
        if(length(input$CPDB_databases)>0){
          x <- input$CPDB_databases
          return(x)
        }
      })
      
      CPDB_databasesFULL <- reactive({
        if(length(input$CPDB_databasesFULL)>0){
          x <- input$CPDB_databasesFULL
          return(x)
        }
      })
      
      MONET_background_file <- reactive({
        if(length(input$MONET_background_file)>0){
          x <- input$MONET_background_file$datapath
          return(x)
        }
      })
      
      CPDB_database_file <- reactive({
        if(length(input$CPDB_database_file)>0){
          x <- input$CPDB_database_file$datapath
          return(x)
        }
      })
      
      
      
      observeEvent(input$submit_MONET_pathways, {
        
        shinyjs::addClass(id = "submit_MONET_pathways", class = "loading dots")
        shinyjs::disable("submit_Full_analysis")
        shinyjs::disable("submit_Network_analysis")
        shinyjs::disable("submit_MONET_pathways")
        shinyjs::disable("submit_Circos_and_auxiliary")
        shinyjs::disable("submit_MONET_analysis")

        cat(file = stderr(), "Running MONET pathways module\n")
        
        if(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern = TRUE)) < 300000){
          showMemFull(TRUE)
        } else {
          
          if(length(working_dir) !=0 ) {
            
            if(biomart_check() == 0) {
              showVal(FALSE)
              showErrorBiomartServerMessage(TRUE)
            } else {
              
              dir.create(paste0(working_dir,"/Results_",Results_index()))
              
              token <- UUIDgenerate()
              message(paste0("running task for token: ", token))
              
              if(is.null(jobs[[token]])) {
                
                if(Type_of_analysis() == 'Run the full pipeline') {
                  phenotype_names <- unique(unlist(strsplit(phenotype_comparison(), "vs")))
                  phenotype_comparison <- phenotype_comparison()
                } else {
                  phenotype_names <- unique(unlist(strsplit(phenotype_comparison1(), "vs")))
                  phenotype_comparison <- phenotype_comparison1()
                }
                
                analyses_flag <- 1
                if(length(which(BioMart_Dataset() %in% biomart_available_datasets$Dataset)) == 0) {
                  showVal(FALSE)
                  showErrorBioMartMessage(TRUE)
                  analyses_flag <- 0
                }
                
                if(length(which(CPDB_databases() %in% CPDB_databases_list$Database)) == 0) {
                  showVal(FALSE)
                  showErrorCPDBMessage(TRUE)
                  analyses_flag <- 0
                }
                
                if(analyses_flag) {
                  
                  source(
                    paste0(working_dir, "/Integrative_analysis.R"),
                    local = TRUE,
                    echo = TRUE,
                    max.deparse.length = 100000
                  )

                  jobs[[token]] <- callr::r_bg(
                    MONET_pathways,
                    supervise = FALSE,
                    stdout = paste0(working_dir, "/Results_", Results_index(), "/out_Pathways.txt"),
                    package = TRUE,
                    stderr = paste0(working_dir, "/Results_", Results_index(), "/log_Pathways.txt"),
                    args = list(
                      working_dir = working_dir,
                      CPDB_databases = CPDB_databases(),
                      MONET_background_file = MONET_background_file(),
                      phenotype_names = phenotype_names,
                      phenotype_comparison = phenotype_comparison,
                      CPDB_database_file = CPDB_database_file(),
                      Results_index = Results_index(),
                      BioMart_Dataset = BioMart_Dataset()
                    )
                  )

                  pars_ls = list(
                    NOODAI_module = "Monet pathway analysis",
                    CPDB_databases = CPDB_databases(),
                    MONET_background_file = MONET_background_file(),
                    phenotype_names = phenotype_names,
                    phenotype_comparison = phenotype_comparison,
                    CPDB_database_file = CPDB_database_file(),
                    Results_index = Results_index(),
                    BioMart_Dataset = BioMart_Dataset()
                  )

                  save_run_info(
                    pars_ls,
                    out_dir = file.path(
                      working_dir,
                      paste0("Results_", Results_index())
                    )
                  )
                  
                  showVal(FALSE)
                  showResPathways1(TRUE)
                  
                }
              }
              
            }
          }
        }

        cat(file = stderr(), "Finished running MONET pathways module\n")
        
      })
      
      
      # #####End MONET pathways analysis#######
      
      
      # #####Start MONET images generation#######
      
      files_edges_path <- reactive({
        if(length(input$files_edges_path)>0){
          x <- input$files_edges_path
          if(x==""){
            return(NULL)
          }else{
            return(x)
          }
        }
      })
      
      TF_Database <- reactive({
        if(length(input$TF_Database)>0){
          x <- input$TF_Database$datapath
          return(x)
        }
      })
      
      centralities_file <- reactive({
        if(length(input$centralities_file)>0){
          x <- input$centralities_file$datapath
          return(x)
        }
      })
      
      file_extension <- reactive({
        if(length(input$file_extension)>0){
          x <- input$file_extension
          if(x==""){
            return(NULL)
          }else{
            return(x)
          }
        }
      })
      
      Kinome_database <- reactive({
        if(length(input$Kinome_database)>0){
          x <- input$Kinome_database$datapath
          return(x)
        }
      })
      

      observeEvent(input$submit_Circos_and_auxiliary, {
        
        shinyjs::addClass(id = "submit_Circos_and_auxiliary", class = "loading dots")
        shinyjs::disable("submit_Full_analysis")
        shinyjs::disable("submit_Network_analysis")
        shinyjs::disable("submit_MONET_pathways")
        shinyjs::disable("submit_Circos_and_auxiliary")
        shinyjs::disable("submit_MONET_analysis")

        cat(file = stderr(), "Running plotting module\n")

        NOODAI_object = new("NOODAI")
        
        if(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE)) < 300000) {
          showMemFull(TRUE)
        } else {
          
          if(length(working_dir) != 0){
            
            if(biomart_check() == 0) {
              showVal(FALSE)
              showErrorBiomartServerMessage(TRUE)
            } else {
              
              dir.create(paste0(working_dir,"/Results_",Results_index()))
              
              token <- UUIDgenerate()
              message(paste0("running task for token: ", token))
              
              if(is.null(jobs[[token]])) {
                
                if(Type_of_analysis() == 'Run the full pipeline') {
				file_DEA_names <- file_DEA_names3()
                  phenotype_names <- unique(unlist(strsplit(phenotype_comparison(), "vs")))
                  phenotype_comparison <- phenotype_comparison()
                } else {
                  phenotype_names <- unique(unlist(strsplit(phenotype_comparison1(), "vs")))
                  phenotype_comparison <- phenotype_comparison1()
				  if(Use_precompiled_database() == '0'){file_DEA_names <- file_DEA_names1()}
                  if(Use_precompiled_database() == '1'){file_DEA_names <- file_DEA_names2()}
                }
                
                analyses_flag <- 1
                if(length(which(BioMart_Dataset() %in% biomart_available_datasets$Dataset)) == 0) {
                  showVal(FALSE)
                  showErrorBioMartMessage(TRUE)
                  analyses_flag <- 0
                }
                
                if(analyses_flag) {
                  source(
                    paste0(working_dir, "/Integrative_analysis.R"),
                    local = TRUE,
                    echo = TRUE,
                    max.deparse.length = 100000
                  )

                  NOODAI_object@working_dir = working_dir
                  NOODAI_object@phenotype_names = phenotype_names
                  NOODAI_object@phenotype_comparison = phenotype_comparison
                  NOODAI_object@edge_file_path = files_edges_path()
                  NOODAI_object@centralities_file = centralities_file()
                  NOODAI_object@TF_Database = TF_Database()
                  NOODAI_object@file_extension = file_extension()
                  NOODAI_object@Results_index = Results_index()
                  NOODAI_object@Kinome_database = Kinome_database()

                  message(paste0("working dir: ", working_dir))
                  message(paste0("phenotype names: ", phenotype_names))
                  message(paste0("phenotype comparison: ", phenotype_comparison))
                  message(paste0("file edges path: ", files_edges_path()))
                  message(paste0("centralities file: ", centralities_file()))
                  message(paste0("TF database: ", TF_Database()))
                  message(paste0("file extension: ", file_extension()))
                  message(paste0("Results index: ", Results_index()))
                  message(paste0("Kinome database: ", Kinome_database()))

                  jobs[[token]] <- callr::r_bg(
                    Circos_and_auxiliary,
                    supervise = FALSE,
                    stdout = paste0(working_dir, "/Results_", Results_index(), "/out_Images.txt"),
                    package = TRUE,
                    stderr = paste0(working_dir, "/Results_", Results_index(), "/log_Images.txt"),
                    args = list(
                      working_dir = working_dir,
                      NOODAI_object = NOODAI_object,
                      phenotype_names = phenotype_names,
                      phenotype_comparison = phenotype_comparison,
                      files_edges_path = files_edges_path(),
                      centralities_file = centralities_file(),
                      TF_Database = TF_Database(),
                      file_extension = file_extension(),
                      Results_index = Results_index(),
                      Kinome_database = Kinome_database(),
                      BioMart_Dataset = BioMart_Dataset(),
                      Cytoscape_style_file = NULL,
                      file_DEA_names = file_DEA_names,
					  flag_run_cytoscape = flag_run_cytoscape()
                    )
                  )
                  
                  pars_ls = list(
                    NOODAI_module = "Image generation",
                    phenotype_names = phenotype_names,
                    phenotype_comparison = phenotype_comparison,
                    files_edges_path = files_edges_path(),
                    centralities_file = centralities_file(),
                    TF_Database = TF_Database(),
                    file_extension = file_extension(),
                    Results_index = Results_index(),
                    Kinome_database = Kinome_database(),
                    BioMart_Dataset = BioMart_Dataset(),
                    Cytoscape_style_file = NULL,
                    file_DEA_names = file_DEA_names,
					flag_run_cytoscape = flag_run_cytoscape()
                  )

                  save_run_info(
                    pars_ls,
                    out_dir = file.path(
                      working_dir,
                      paste0("Results_", Results_index())
                    )
                  )

                  showVal(FALSE)
                  showResImages1(TRUE)
                  
                }
              }
              
            }
          }
        }

        cat(file = stderr(), "Finished running plotting module\n")
        
      })
      
      #####End MONET images generation#######  
      
      #####Run the full Network analysis####### 
      
      
      observeEvent(input$submit_Full_analysis, {

        cat(file = stderr(), "Running full analysis\n")
        
        shinyjs::addClass(id = "submit_Full_analysis", class = "loading dots")
        shinyjs::disable("submit_Full_analysis")
        shinyjs::disable("submit_Network_analysis")
        shinyjs::disable("submit_MONET_pathways")
        shinyjs::disable("submit_Circos_and_auxiliary")
        shinyjs::disable("submit_MONET_analysis")
        
        if(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE)) < 300000) {
          showMemFull(TRUE)
        } else {
          
          if(length(working_dir) != 0) {
            
            if(biomart_check() == 0) {
              showVal(FALSE)
              showErrorBiomartServerMessage(TRUE)
            } else {
              
              dir.create(paste0(working_dir,"/Results_",Results_index()), recursive = T)

              shinyalert::shinyalert(
                "Thank you for using NOODAI!",
                paste(
                  "This is your results ID:",
                  Results_index(),
                  "Make sure to copy it somewhere safe to be able to access your results later. ",
                  "In addition, please read the additional footnote information.",
                  sep = "\n"),
                type = "success"
              )

              token <- uuid::UUIDgenerate()
              message(paste0("running task for token: ", token))
              # the if statement is to avoid rerunning a job again
              if(is.null(jobs[[token]])){
                # call the job in the background session
                if(Type_of_analysis() == 'Run the full pipeline') {

                  if(input$Demo_values >= 1) {
                    file_DEA_names <- file_DEA_names5()
                  } else {
                    file_DEA_names <- file_DEA_names3()
                  }

                  Use_precompiled_databaseVal <- Use_precompiled_database1()
                  LookUp_table_fileVal <- LookUp_table_file1()
                  BioGRID_data_fileVal <- BioGRID_data_file1()
                  STRING_data_fileVal <- STRING_data_file1()
                  IntAct_data_fileVal <- IntAct_data_file1()
                  phenotype_names <- unique(unlist(strsplit(phenotype_comparison(), "vs")))
                  phenotype_comparison <- phenotype_comparison()
                  BioMart_Dataset <- BioMart_Dataset1()

                } else {

                  Use_precompiled_databaseVal <- Use_precompiled_database()
                  LookUp_table_fileVal <- LookUp_table_file()
                  BioGRID_data_fileVal <- BioGRID_data_file()
                  STRING_data_fileVal <- STRING_data_file()
                  IntAct_data_fileVal <- IntAct_data_file()
                  phenotype_names <- unique(unlist(strsplit(phenotype_comparison1(), "vs")))
                  phenotype_comparison <- phenotype_comparison1()
                  BioMart_Dataset <- BioMart_Dataset()
                  
                  if(Use_precompiled_database() == '0') {
                    file_DEA_names <- file_DEA_names1()
                  }

                  if(Use_precompiled_database() == '1') {
                    file_DEA_names <- file_DEA_names2()
                  }

                }
                
                analyses_flag <- 1
                
				
				
				                                    for (i in 1:length(phenotype_comparison)){
                                      if(length(unlist(sapply(phenotype_names,FUN=function(x){grep(x,phenotype_comparison[i])})))!=2){
                                        analyses_flag<-0
                                        showVal(FALSE)
                                        showErrorContrastsMessage(TRUE)
                                        showErrorSampleMessage(TRUE)
                                        special_error_messageContrast("Cause: There is a mismatch between the conditions names and contrasts!")
                                        special_error_messageSample("Cause: There is a mismatch between the conditions names and contrasts!")
                                        break()
                                      }
                                    }
                                    
                                    if(length(file_DEA_names)==0){
                                      analyses_flag<-0
                                      showVal(FALSE)
                                      showErrorDatasetMessage(TRUE)
                                      special_error_messageDataset("Cause: Either the archive is empty or wrongly formatted!")
                                    }
                                    if(length(phenotype_comparison)==0){
                                      analyses_flag<-0
                                      showVal(FALSE)
                                      showErrorContrastsMessage(TRUE)
                                      special_error_messageContrast("Cause: No contrasts were correctly provided!")
                                    }
                                    if(length(grep("vs", phenotype_comparison))!=length(phenotype_comparison)){
                                      analyses_flag<-0
                                      showVal(FALSE)
                                      showErrorContrastsMessage(TRUE)
                                      special_error_messageContrast("Cause: 'vs' was not used to separate the conditions in the contrasts!")
                                    }
                                    
                                    if(length(phenotype_names)<2){
                                      analyses_flag<-0
                                      showVal(FALSE)
                                      showErrorSampleMessage(TRUE)
                                      special_error_messageSample("Cause: There are fewer than two condition names!")
                                    }
                                    
                                    if(analyses_flag){
									contrast_uniprot_counts <- setNames(vector("list", length(phenotype_comparison)), phenotype_comparison)
                                      for (j in 1:length(file_DEA_names)){
                                        file_Sheets1 <- readxl::excel_sheets(file_DEA_names[j])
                                        if(length(which(phenotype_comparison %in% file_Sheets1))!=length(phenotype_comparison)){
                                          analyses_flag<-0
                                          showVal(FALSE)
                                          showErrorContrastsMessage(TRUE)
                                          showErrorDatasetMessage(TRUE)
                                          special_error_messageContrast("Cause: The provided contrasts are not found in the Excel tables!")
                                          special_error_messageDataset("Cause: The provided contrasts are not found in the Excel tables!")
                                          break()
                                        }
                                        for (i in 1:length(phenotype_comparison)){
                                          aux <- as.data.frame(readxl::read_excel(file_DEA_names[j],sheet = phenotype_comparison[i]))
                                          #Get the corresponding entrez ids from the uniprot ones#####################
                                          if(length(aux)>0){
                                            if(length(which(colnames(aux) %in% "UniProt_ChEBI"))==0){
                                              showVal(FALSE)
                                              showErrorFormatMessage(TRUE)
                                              analyses_flag<-0
                                              break()
                                            }
                                            aux1 <- aux$UniProt_ChEBI[which(aux$UniProt_ChEBI %in% uniprot_ids_maps$UniProt_ChEBI)]
                                            if(length(aux1)<(length(aux$UniProt_ChEBI)/3)){
                                              showVal(FALSE)
                                              showErrorFormatMessage(TRUE)
                                              analyses_flag<-0
                                              break()
                                            }
											contrast_uniprot_counts[[phenotype_comparison[i]]] <- unique(c(contrast_uniprot_counts[[phenotype_comparison[i]]], aux1))
                                          }
                                        }
                                      }
									  low_input_flag <- any(sapply(contrast_uniprot_counts, function(x) length(x) < 150))
									  if (low_input_flag) {
											showWarningLowInput(TRUE)
										}
                                    }
                                    
                                    if(length(which(BioMart_Dataset %in% biomart_available_datasets$Dataset))==0){
                                      showVal(FALSE)
                                      showErrorBioMartMessage(TRUE)
                                      analyses_flag<-0
                                    }
                                    
                                    if(length(which(CPDB_databasesFULL() %in% CPDB_databases_list$Database))==0){
                                      showVal(FALSE)
                                      showErrorCPDBMessage(TRUE)
                                      analyses_flag<-0
                                    }
                                    
                                    if(length(grep(Monet_exp,Monet_method_stringFULL()))==0){
                                      showVal(FALSE)
                                      showErrorMonetMessage(TRUE)
                                      analyses_flag<-0
                                    }
				

                
                if(analyses_flag) {

                  message(paste0("Token is: ", token))
                  message(paste0("Working_dir is: ", working_dir))
                  message(paste0("BioGRID_data_file is: ", BioGRID_data_fileVal))
                  message(paste0("STRING_data_file is: ", STRING_data_fileVal))
                  message(paste0("IntAct_data_file is: ", IntAct_data_fileVal))
                  message(paste0("file_DEA_names is: ", file_DEA_names))
                  message(paste0("phenotype_names is: ", phenotype_names))
                  message(paste0("phenotype_comparison is: ", phenotype_comparison))
                  message(paste0("Use_precompiled_database is: ", Use_precompiled_databaseVal))
                  message(paste0("LookUp_table_file is: ", LookUp_table_fileVal))
                  message(paste0("Results_index is: ", Results_index()))
                  message(paste0("edge_file_path is: ", edge_file_path()))
                  message(paste0("monet_path is: ", monet_path()))
                  message(paste0("Monet_method_string is: ", Monet_method_stringFULL()))
                  message(paste0("tmp_bin_folder is: ", tmp_bin_folder()))
                  message(paste0("CPDB_databases is: ", CPDB_databasesFULL()))
                  message(paste0("MONET_background_file is: ", MONET_background_file()))
                  message(paste0("CPDB_database_file is: ", CPDB_database_file()))
                  message(paste0("files_edges_path is: ", files_edges_path()))
                  message(paste0("centralities_file is: ", centralities_file()))
                  message(paste0("TF_Database is: ", TF_Database()))
                  message(paste0("file_extension is: ", file_extension()))
                  message(paste0("Kinome_database is: ", Kinome_database()))
                  message(paste0("BioMart_Dataset is: ", BioMart_Dataset))
                  message(paste0("Client_email is: ", Client_email()))
                  message(paste0("weight_penalty is: ", weight_penalty1()))

                  NOODAI_object = new("NOODAI")

                  NOODAI_object@token = token
                  NOODAI_object@working_dir = working_dir
                  NOODAI_object@BioGRID_data_file = BioGRID_data_fileVal
                  NOODAI_object@STRING_data_file = STRING_data_fileVal
                  NOODAI_object@IntAct_data_file = IntAct_data_fileVal
                  NOODAI_object@file_DEA_names = file_DEA_names
                  NOODAI_object@phenotype_names = phenotype_names
                  NOODAI_object@phenotype_comparison = phenotype_comparison
                  NOODAI_object@Use_precompiled_database = Use_precompiled_databaseVal
                  NOODAI_object@LookUp_table_file = LookUp_table_fileVal
                  NOODAI_object@Results_index = Results_index()
                  NOODAI_object@edge_file_path = edge_file_path()
                  NOODAI_object@MONET_background_file = MONET_background_file()
                  NOODAI_object@CPDB_database_file = CPDB_database_file()
                  NOODAI_object@files_edges_path = edge_file_path()
                  NOODAI_object@centralities_file = centralities_file()
                  NOODAI_object@TF_Database = TF_Database()
                  NOODAI_object@file_extension = file_extension()
                  NOODAI_object@Kinome_database = Kinome_database()
                  NOODAI_object@BioMart_Dataset = BioMart_Dataset
                  NOODAI_object@Client_email = Client_email()
                  NOODAI_object@weight_penalty = weight_penalty1()
                  NOODAI_object@status = "To start"

                  saveRDS(
                    NOODAI_object,
                    file.path(
                      working_dir,
                      paste0("Results_", Results_index()),
                      "NOODAI_object.RDS"
                    )
                  )
                  
                  jobs[[token]] <- callr::r_bg(
                    run_token,
                    supervise = FALSE,
                    stdout = paste0(
                      working_dir,
                      "/Results_",
                      Results_index(),
                      "/out.txt"
                    ),
                    package = TRUE,
                    stderr = paste0(
                      working_dir,
                      "/Results_",
                      Results_index(),
                      "/log.txt"
                    ),
                    args = list(
                      token = token,
                      working_dir = working_dir,
                      NOODAI_object = NOODAI_object,
                      BioGRID_data_file = BioGRID_data_fileVal,
                      STRING_data_file = STRING_data_fileVal,
                      IntAct_data_file = IntAct_data_fileVal,
                      file_DEA_names = file_DEA_names,
                      phenotype_names = phenotype_names,
                      phenotype_comparison = phenotype_comparison,
                      splicing_file_name = splicing_file_name(),
                      Use_precompiled_database = Use_precompiled_databaseVal,
                      LookUp_table_file = LookUp_table_fileVal,
                      Results_index = Results_index(),
                      edge_file_path = edge_file_path(),
                      monet_path = monet_path(),
                      Monet_method_string = Monet_method_stringFULL(),
                      tmp_bin_folder = tmp_bin_folder(),
                      CPDB_databases = CPDB_databasesFULL(),
                      MONET_background_file = MONET_background_file(),
                      CPDB_database_file = CPDB_database_file(),
                      files_edges_path = files_edges_path(),
                      centralities_file = centralities_file(),
                      TF_Database = TF_Database(),
                      file_extension = file_extension(),
                      Kinome_database = Kinome_database(),
                      BioMart_Dataset = BioMart_Dataset,
                      Client_email = Client_email(),
                      weight_penalty = weight_penalty1(),
                      Cytoscape_style_file = NULL,
					  flag_run_cytoscape = flag_run_cytoscape1()
                    )
                  )

                  pars_ls = list(
                    NOODAI_module = "Full pipeline",
                    BioGRID_data_file = BioGRID_data_fileVal,
                    STRING_data_file = STRING_data_fileVal,
                    IntAct_data_file = IntAct_data_fileVal,
                    file_DEA_names = file_DEA_names,
                    phenotype_names = phenotype_names,
                    phenotype_comparison = phenotype_comparison,
                    splicing_file_name = splicing_file_name(),
                    Use_precompiled_database = Use_precompiled_databaseVal,
                    LookUp_table_file = LookUp_table_fileVal,
                    Results_index = Results_index(),
                    edge_file_path = edge_file_path(),
                    monet_path = monet_path(),
                    Monet_method_string = Monet_method_stringFULL(),
                    tmp_bin_folder = tmp_bin_folder(),
                    CPDB_databases = CPDB_databasesFULL(),
                    MONET_background_file = MONET_background_file(),
                    CPDB_database_file = CPDB_database_file(),
                    files_edges_path = files_edges_path(),
                    centralities_file = centralities_file(),
                    TF_Database = TF_Database(),
                    file_extension = file_extension(),
                    Kinome_database = Kinome_database(),
                    BioMart_Dataset = BioMart_Dataset,
                    Client_email = Client_email(),
                    weight_penalty = weight_penalty1(),
                    Cytoscape_style_file = NULL,
					flag_run_cytoscape = flag_run_cytoscape1()
                  )

                  save_run_info(
                    pars_ls,
                    out_dir = file.path(
                      working_dir,
                      paste0("Results_", Results_index())
                    )
                  )
                  
                  showVal(FALSE)
                  showRes1(TRUE)

                }
              }
            }
          }
        }
        message("Finished with full pipeline! :D")
        
      })
      
      #####End the full Network analysis####### 
      
      #####Download the results#######
      
      
      Results_index_save <- reactive({
        if(length(input$Results_index_save)>0){
          x <- input$Results_index_save
          x <- gsub(" ", "", x, fixed = TRUE)
          if(x==""){
            return(NULL)
          }else{
            return(as.character(x))
          }
        }
      })
      
      observe({
        
        
        useReactiveDownload <- function(input, output, session, df){
          
          downloadHandler(
            filename = function() {
              paste0("Results_",Sys.Date(),".zip")
            },
            content = function(file) {
              ns <- paste0(working_dir,"/Results_",df,"/ResultsZip.zip")
              print(ns)
              file.copy(ns, file)
            },
            contentType = "application/zip"
          )
        }
        
        output$downloadData <- useReactiveDownload(input, output, session, df = Results_index_save())
        
      })
      
      
      observe({
        
        useReactiveDownload <- function(input,output,session,df){
          
          downloadHandler(
            filename = function() {
              paste0("Demo_Input_",Sys.Date(),".zip")
            },
            content = function(file) {
              ns <- df
              print(ns)
              file.copy(ns, file)
            },
            contentType = "application/zip"
          )
        }
        
        output$Demo_download <- useReactiveDownload(input,output,session,df = paste0(working_dir,"/Demo_data.zip"))
        
      })
      
      
      observe({
        
        useReactiveDownload <- function(input,output,session,df){
          
          downloadHandler(
            filename = function() {
              paste0("Demo_Results_",Sys.Date(),".zip")
            },
            content = function(file) {
              ns <- df
              print(ns)
              file.copy(ns, file)
            },
            contentType = "application/zip"
          )
        }
        
        output$Demo_Reults_download <- useReactiveDownload(input,output,session,df = paste0(working_dir,"/Results_Demo_Precompiled.zip"))
        
      })
      
      #####End Download the results#######
      
      
      #####Output#######
      
      output$MemFull <- reactive({
        return(showMemFull() == "TRUE")
      })
      outputOptions(output, "MemFull", suspendWhenHidden = FALSE)
      
      output$text_output_MemFull <- renderText({ paste0("The server is currently at maximum capacity! Please try again in at least 30 minutes!") })
      
      output$panelStatus <- reactive({
        return(showVal() == "TRUE")
      })
      outputOptions(output, "panelStatus", suspendWhenHidden = FALSE)
      
      output$showRes <- reactive({
        return(showRes1() == "TRUE")
      })
      outputOptions(output, "showRes", suspendWhenHidden = FALSE)
      
      output$verb <- renderUI({
        HTML(
          paste0(
            "Your results folder has the following ID: ",
            Results_index(),
            "<br>Additionally, you can download the results from: omics-oracle.com/www/Downloads/",
            Results_index(),
            ".zip"
          )
        )
      })
      output$text_output_string <- renderText({ "Thank you for submitting the analysis request. If you have provided an email address you will be notified when the analysis is completed. Otherwise, we recommend to come back in 20-120 minutes to download the results."})
      output$text_output_string2 <- renderText({ "Please save this information if you did not provide an email address. Your analysis will be automatically removed from the server after 5 days."})
      
      
      output$showResCentralities <- reactive({
        return(showResCentralities1() == "TRUE")
      })
      outputOptions(output, "showResCentralities", suspendWhenHidden = FALSE)
      
      output$CentralitiesOut1 <- renderUI({
        HTML(
          paste0(
            "Your results folder has the following ID: ",
            Results_index(),
            "<br>Additionally, you can download the results from: omics-oracle.com/www/Downloads/",
            Results_index(),
            ".zip"
          )
        )
      })
      output$CentralitiesOut2 <- renderText({ "Thank you for submitting the analysis request. Please come back in at least one hour in order to download the results."})
      
      output$showResMONET <- reactive({
        return(showResMONET1() == "TRUE")
      })
      outputOptions(output, "showResMONET", suspendWhenHidden = FALSE)
      
      output$MONETOut1 <- renderUI({
        HTML(
          paste0(
            "Your results folder has the following ID: ",
            Results_index(),
            "<br>Additionally, you can download the results from: omics-oracle.com/www/Downloads/",
            Results_index(),
            ".zip"
          )
        )
      })
      output$MONETOut2 <- renderText({ "Thank you for submitting the analysis request. Please come back in at least one hour in order to download the results."})
      
      output$showResPathways <- reactive({
        return(showResPathways1() == "TRUE")
      })
      outputOptions(output, "showResPathways", suspendWhenHidden = FALSE)
      
      output$PathwaysOut1 <- renderUI({
        HTML(
          paste0(
            "Your results folder has the following ID: ",
            Results_index(),
            "<br>Additionally, you can download the results from: omics-oracle.com/www/Downloads/",
            Results_index(),
            ".zip"
          )
        )
      })
      output$PathwaysOut2 <- renderText({ "Thank you for submitting the analysis request. Please come back in at least one hour in order to download the results."})
      
      
      output$showResImages <- reactive({
        return(showResImages1() == "TRUE")
      })
      outputOptions(output, "showResImages", suspendWhenHidden = FALSE)
      
      output$ImagesOut1 <- renderUI({
        HTML(
          paste0(
            "Your results folder has the following ID: ",
            Results_index(),
            "<br>Additionally, you can download the results from: omics-oracle.com/www/Downloads/",
            Results_index(),
            ".zip"
          )
        )
      })
      output$ImagesOut2 <- renderText({ "Thank you for submitting the analysis request. Please come back in at least one hour in order to download the results."})
      
      
      output$showErrorFormat <- reactive({
        return(showErrorFormatMessage() == "TRUE")
      })
      outputOptions(output, "showErrorFormat", suspendWhenHidden = FALSE)
      
      output$showErrorFormat1 <- renderText({ "Your Input data does not have the proper format. The names must be UniProt Ids or ChEBI and the column header named UniProt_ChEBI" })
      
	  output$showWarningInput <- reactive({
        return(showWarningLowInput() == "TRUE")
      })
      outputOptions(output, "showWarningInput", suspendWhenHidden = FALSE)
      
      output$showWarningInput1 <- renderText({ "WARNING: For at least one comparison, the input datasets contain in total fewer than 150 unique identifiers. This may result in sparse networks and potentially unreliable outcomes. Please carefully evaluate the MONET decomposition performance" })
      
	  
      
      output$showErrorDataset <- reactive({
        return(showErrorDatasetMessage() == "TRUE")
      })
      outputOptions(output, "showErrorDataset", suspendWhenHidden = FALSE)
      
      output$showErrorDataset1 <- renderText({ paste0("Your input data may be incorrect. Please provide an archive with each omics measurement in an Excel table where each sheet correspond to a contrast.\n",special_error_messageDataset()) })
      
      output$showErrorSample <- reactive({
        return(showErrorSampleMessage() == "TRUE")
      })
      outputOptions(output, "showErrorSample", suspendWhenHidden = FALSE)
      
      output$showErrorSample1 <- renderText({ paste0("The condition names provided may be incorrect. Please use the exact names included in the contrasts and Excel sheets.\n",special_error_messageSample()) })
      
      
      output$showErrorContrasts <- reactive({
        return(showErrorContrastsMessage() == "TRUE")
      })
      outputOptions(output, "showErrorContrasts", suspendWhenHidden = FALSE)
      
      output$showErrorContrasts1 <- renderText({ paste0("The provided contrasts may be incorrect. Please format the contrasts (i.e. comparison groups) to have the same names as in the uploaded Excel sheets.\n",special_error_messageContrast()) })
      
      output$showErrorBioMart <- reactive({
        return(showErrorBioMartMessage() == "TRUE")
      })
      outputOptions(output, "showErrorBioMart", suspendWhenHidden = FALSE)
      
      output$showErrorBioMart1 <- renderText({ "The provided BioMart dataset may be incorrect. Please provide one in accordance with the instructions." })
      
      output$showErrorMonet <- reactive({
        return(showErrorMonetMessage() == "TRUE")
      })
      outputOptions(output, "showErrorMonet", suspendWhenHidden = FALSE)
      
      output$showErrorMonet1 <- renderText({ "The provided MONET String may be incorrect. Please format the MONET String in accordance with the instructions." })
      
      output$showErrorCPDB <- reactive({
        return(showErrorCPDBMessage() == "TRUE")
      })
      outputOptions(output, "showErrorCPDB", suspendWhenHidden = FALSE)
      
      output$showErrorCPDB1 <- renderText({ "The provided Pathway Database may be incorrect. Please provide a valid name in accordance with the instructions." })
      
      output$showErrorBiomartServer <- reactive({
        return(showErrorBiomartServerMessage() == "TRUE")
      })
      outputOptions(output, "showErrorBiomartServer", suspendWhenHidden = FALSE)
      
      output$showErrorBiomartServer1 <- renderText({ "The BioMart servers are down! We are sorry for the inconvenience!" })
      
      # #####End Output#######

}


shinyApp(mod_NOODAI_pipeline_ui, mod_NOODAI_pipeline_server)