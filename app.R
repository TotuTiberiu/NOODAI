#Main Shiny Interface of the NOODAI software, deployed on the server.

# 
#     Copyright © 2024, Empa, Tiberiu Totu.
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


setwd("/home/omics/Linux_test/")

#setwd("C:/work/Files/Algorithms/DO_server/V12/Linux_test")

working_dir <- getwd()


uniprot_ids_maps <- read.table(paste0(getwd(),"/Databases/Uniprot_FullIds.txt"),header = TRUE,fill = TRUE)
biomart_available_datasets <- read.table(paste0(getwd(),"/Databases/BioMart_datasets.txt"),header = TRUE,fill = TRUE)
CPDB_databases_list <- read.table(paste0(getwd(),"/Databases/CPDB_databases_list.txt"),header = TRUE,fill = TRUE)
Monet_exp <- "--method=[A-Z][0-9] --avgk=?[0-9]+ --linksdir=[a-z]+"


biomart_check <- function(){
  
  tryCatch(
    {
      mart <- useMart('ENSEMBL_MART_ENSEMBL')
      mart <- useDataset("hsapiens_gene_ensembl", mart)
      return(1)
    },
    error = function(e){
      tryCatch(
        {
          mart <- useMart('ENSEMBL_MART_ENSEMBL',host='https://asia.ensembl.org')
          mart <- useDataset("hsapiens_gene_ensembl", mart)
          return(1)
        },
        error = function(e){
          tryCatch(
            {
              mart <- useMart('ENSEMBL_MART_ENSEMBL',host='https://useast.ensembl.org')
              mart <- useDataset("hsapiens_gene_ensembl", mart)
              return(1)
            }, error = function(e){
              return(0)
              stop()
            }
          )
        })
    })
  
}

jobs <- list()

run_token <- function(token,working_dir,BioGRID_data_file,STRING_data_file,IntAct_data_file,file_DEA_names,
                      phenotype_names,phenotype_comparison,splicing_file_name,
                      Use_precompiled_database,LookUp_table_file,Results_index,
                      edge_file_path,monet_path,Monet_method_string,tmp_bin_folder,
                      CPDB_databases,MONET_background_file,CPDB_database_file,
                      files_edges_path,centralities_file,TF_Database,file_extension,Kinome_database,BioMart_Dataset,Client_email){
  
  
  source(paste0(working_dir,"/Integrative_analysis.R"), local = TRUE, echo=TRUE, max.deparse.length=100000)
  
  Integrative_Network_analysis(working_dir,BioGRID_data_file,STRING_data_file,IntAct_data_file,file_DEA_names,
                               phenotype_names,phenotype_comparison,splicing_file_name,
                               Use_precompiled_database,LookUp_table_file,Results_index,BioMart_Dataset)
  
  MONET_analysis(working_dir,edge_file_path,monet_path,Monet_method_string,tmp_bin_folder,Results_index)
  
  MONET_pathways(working_dir,CPDB_databases,MONET_background_file,phenotype_names,
                 phenotype_comparison,CPDB_database_file,Results_index,BioMart_Dataset)
  
  Circos_and_auxiliary(working_dir,phenotype_names,phenotype_comparison,files_edges_path,
                       centralities_file,TF_Database,file_extension,Results_index,Kinome_database,BioMart_Dataset)
  if(length(Client_email)!=0){
    system(paste0("echo ", "Your analysis associated with the results folder:",Results_index, " is done. If something did not went as expected check the log files. Otherwise, contact the authors. | mutt -s ",'"NOODAI analysis results are available" ',Client_email))
  }
}



css <- '
.tooltip {
  pointer-events: none;
}
.tooltip > .tooltip-inner {
  pointer-events: none;
  background-color: #186ea0;
  color: #FFFFFF;
  border: 1px solid green;
  padding: 5px;
  font-size: 15px;
  font-style: italic;
  text-align: justify;
  margin-left: 0;
  max-width: 1000px;
}
.tooltip > .arrow::before {
  border-right-color: #73AD21;
}
  
  
.footer {
background-attachment: scroll;
    background-position: 0% 0%;
    position: relative;
    left: 0pt;
bottom:0;
width:100%;

color: black;
padding: 0px;
background-color: gray;
z-index: 1000;
font-size: 12px;
text-align: center;
text-justify: inter-character;
font-family: Verdana
display: inline-block;
}

.tab-content {
  margin-bottom: 75px; /* Adjust the value as needed to create space for the footer */
  min-height:100vh;
}

.loading {
                display: inline-block;
                overflow: hidden;
                height: 2.3em;
                margin-top: -0.0em;
                line-height: 1.5em;
                vertical-align: text-bottom;
                box-sizing: border-box;
            }
            .loading.dots::after {
                text-rendering: geometricPrecision;
                content: "⠋\\A⠙\\A⠹\\A⠸\\A⠼\\A⠴\\A⠦\\A⠧\\A⠇\\A⠏";
                animation: spin10 1s steps(10) infinite;
                animation-duration: 1s;
                animation-timing-function: steps(10);
                animation-delay: 0s;
                animation-iteration-count: infinite;
                animation-direction: normal;
                animation-fill-mode: none;
                animation-play-state: running;
                animation-name: spin10;
            }
            .loading::after {
                display: inline-table;
                white-space: pre;
                text-align: left;
            }
            @keyframes spin10 { to { transform: translateY(-15.0em); } }

'

js <- "
$(function () {
  $('[data-toggle=tooltip]').tooltip()
})
"



DownloadButtonNoLabel <- function(outputId, label = "Download",width){
  tags$a(id = outputId, class = "btn btn-default shiny-download-link", href = "", 
         target = "_blank", download = NA, NULL, label,width)
}

ui <- shinyUI(fluidPage(useShinyjs(),
                        tags$head(
                          tags$style(HTML(css)),
                          tags$script(HTML(js))
                        ),tags$head(tags$link(rel="icon", href="https://omics-oracle.com/favicon.ico")),
                        navbarPage(title=tags$a(style="color: inherit;text-decoration: inherit;font-size:25px;",href="https://omics-oracle.com/",tags$img(src='logo.png',height="2%", width="2%"),"Network Oriented multi-Omics Data Analysis and Integration"),theme = bs_theme(version = 5, bootswatch="journal"),
                                   tags$script(
                                     HTML("var header = $('.navbar > .container-fluid');
                              header.append('<div style=\"float:right; padding-top: 8px;padding-right: 13px\"><button id=\"signin\" type=\"button\" class=\"btn btn-primary action-button\" onclick=\"RedirectSite()\">Tutorials</button></div>')")
                                   ),
                                   tags$script(
                                     HTML("
          var RedirectSite = function(tabName){window.open('https://omics-oracle.com/');}
                                          ")
                                   ),
                                   tabPanel("",
                                            tabPanel("",
                                                     fluidPage(
                                                       conditionalPanel( condition = "output.panelStatus",
                                                                         conditionalPanel( condition = "output.MemFull",
                                                                                           fluidRow(
                                                                                             column(12,align="center",span(textOutput("text_output_MemFull"),style="color:red; font-weight:bold;")),
                                                                                           )),
                                                                         
                                                                         fluidRow(radioGroupButtons(inputId = "Type_of_analysis",choices = c("Run the full pipeline", "Custom algorithms", "Results download"),status = "primary",justified = TRUE,)
                                                                         ),
                                                                         conditionalPanel( condition = "input.Type_of_analysis == 'Run the full pipeline'",
                                                                                           fluidRow(
                                                                                             column(1,div(style = "margin-top: 30px;",span(actionButton(width = "100%",inputId="Demo_values",label=tooltip(trigger=list("Demo",icon("info-circle")),"Upload the Demo data"))))),
                                                                                             column(2,span(textInput(width = "100%",inputId="phenotype_names", label=tooltip(trigger=list("Samples names",icon("info-circle")),"Provide the names of the samples (phenotypes) analyzed and written in the xlsx sheets of the input omics files (DO NOT USE '_' for names!). Format example: Cool,Hot3,mkey")))),
                                                                                             column(4,span(textInput(width = "100%",inputId="phenotype_comparison", label=tooltip(trigger=list("Samples contrasts",icon("info-circle")),"Provide the analyzed samples (phenotypes) contrasts. They should be the names of the xlsx sheets. Format example: CoolvsHot3,mkeyvsCool,Hot3vsCool")))),
                                                                                             column(4,span(fileInput(width = "100%",inputId="file_DEA_names3", label=tooltip(trigger=list("Omics files archive",icon("info-circle")),"Upload an archive containing the omics files. Each omics file needs to be an Excel table with the sheets containing a 'Protein' column with Uniprot Ids in accordance with the documentation")))),
                                                                                             column(1,div(style = "margin-top: 30px;",span(downloadButton(width = "100%",outputId="Demo_download",icon=NULL,label=tooltip(trigger=list("Demo",icon("cloud-arrow-down")),"Download the input demo data")))))
                                                                                           )
                                                                         ),
                                                                         conditionalPanel( condition = "input.Type_of_analysis == 'Custom algorithms'",   
                                                                                           fluidRow(
                                                                                             column(2,offset=5,span(textInput(inputId="Results_index", label=tooltip(trigger=list("Results directory index",icon("info-circle")),"What is the directory index where your NOODAI analysis results are found? The index is provided after the input data is submitted for the analysis")))),
                                                                                           ),
                                                                                           fluidRow(
                                                                                             column(offset=1,2,span(textInput(width = "100%",inputId="phenotype_names1", label=tooltip(trigger=list("Samples names",icon("info-circle")),"Provide the names of the samples (phenotypes) analyzed and written in the xlsx sheets of the input omics files (DO NOT USE '_' for names!). Format example: Cool,Hot3,mkey")))),
                                                                                             column(2,span(textInput(width = "100%",inputId="phenotype_comparison1", label=tooltip(trigger=list("Samples contrasts",icon("info-circle")),"Provide the analyzed samples (phenotypes) contrasts. They should be the names of the xlsx sheets. Format example: CoolvsHot3,mkeyvsCool,Hot3vsCool")))),
                                                                                             column(2,span(textInput(width = "100%",inputId="BioMart_Dataset",value = "hsapiens_gene_ensembl", label=tooltip(trigger=list("BioMart dataset",icon("info-circle")),"Provide the BioMart dataset used for the mapping of the protein Ids. For other organisms than homo sapiens, provide the pre-formated interaction table! The available database is only for humans")))),
                                                                                             column(1,span(textInput(width = "100%",inputId="splicing_file_name",value ="DTU", label=tooltip(trigger=list("DTU file",icon("info-circle")),"Provide the name of the splicing file in the archive. If you use splicing data, the contrast must be symmetric! Leave default for no DTU data")))),
                                                                                             column(3,span(selectInput(width = "100%",inputId="Use_precompiled_database",choices=c("No"="0","Yes"="1"),selected ='1', label=tooltip(trigger=list("Use Pre-compiled Interaction file",icon("info-circle")),"Choose if you would like to use a pre-compiled interaction file, formated in accordance with the instructions. A default interaction file is available. For other organisms than humans, upload a pre-compiled interaction file!"))))
                                                                                           )),
                                                                         
                                                                         conditionalPanel( condition = "input.Type_of_analysis == 'Custom algorithms'",
                                                                                           conditionalPanel( condition = "input.Use_precompiled_database == '0'",
                                                                                                             fluidRow(
                                                                                                               column(offset=1,2,span(fileInput(width = "100%",inputId="BioGRID_data_file", label=tooltip(trigger=list("BioGrid database file",icon("info-circle")),"Upload the BioGrid database file. Leave empty for default version 4.4.218")), style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                                               column(3,span(fileInput(width = "100%",inputId="STRING_data_file", label=tooltip(trigger=list("STRING database file",icon("info-circle")),"Upload the STRING database file. Leave empty for default version 11.5")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                                               column(3,span(fileInput(width = "100%",inputId="IntAct_data_file", label=tooltip(trigger=list("IntAct database file",icon("info-circle")),"Upload the IntAct database file. Leave empty for default release 245")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                                               column(2,span(fileInput(width = "100%",inputId="file_DEA_names1", label=tooltip(trigger=list("Omics files archive",icon("info-circle")),"Upload an archive containing the omics files. Each omics file needs to be an Excel table with the sheets containing a 'Protein' column with Uniprot Ids in accordance with the documentation"))))
                                                                                                             )),
                                                                                           conditionalPanel( condition = "input.Use_precompiled_database == '1'",
                                                                                                             fluidRow(
                                                                                                               column(offset=1,5,span(fileInput(width = "100%",inputId="LookUp_table_file", label=tooltip(trigger=list("Interaction table file",icon("info-circle")),"Upload the pre-formatted interaction database file. Leave empty for default")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                                               column(5,span(fileInput(width = "100%",inputId="file_DEA_names2", label=tooltip(trigger=list("Omics files archive",icon("info-circle")),"Upload an archive containing the omics files. Each omics file needs to be an Excel table with the sheets containing a 'Protein' column with Uniprot Ids in accordance with the documentation"))))
                                                                                                             ))
                                                                         ),
                                                                         fluidRow(
                                                                           conditionalPanel( condition = "input.Type_of_analysis == 'Custom algorithms'",
                                                                                             column(offset=1,10,span(actionButton(width = "100%",inputId="submit_Network_analysis",label=tooltip(trigger=list("Submit",icon("info-circle")),"Run the main network analysis algorithm extracting the top central nodes."))))
                                                                                             
                                                                           )
                                                                         ),
                                                                         
                                                                         # fluidRow(
                                                                         #   column(5,span(switchInput(inline = TRUE,inputId="Specific_analysis_run",onLabel = "Yes",offLabel = "No", label="Do you have a results folder and would like to run a specific analysis?", value = FALSE),span(`data-toggle` = "tooltip", `data-placement` = "right",title = "Provide the path where the network edges are saved.",icon("info-circle")))),
                                                                         # ),
                                                                         
                                                                         conditionalPanel( condition = "input.Type_of_analysis == 'Run the full pipeline'",
                                                                                           fluidRow(
                                                                                             column(offset=1,3,span(textInput(width = "100%",inputId="Monet_method_stringFULL", value = "--method=M1 --avgk=10 --linksdir=undirected", label=tooltip(trigger=list("MONET method",icon("info-circle")),"Provide the MONET analysis method string")))),
                                                                                             column(2,span(textInput(width = "100%",inputId="CPDB_databasesFULL",value = "Reactome", label=tooltip(trigger=list("Pathways databases",icon("info-circle")),"Provide the names of the databases against which the enrichment will be performed. Example: Reactome, Wikipathways")))),
                                                                                             column(2,span(textInput(width = "100%",inputId="BioMart_Dataset1",value = "hsapiens_gene_ensembl", label=tooltip(trigger=list("BioMart dataset",icon("info-circle")),"Provide the BioMart dataset used for the mapping of the protein Ids. For other organisms than homo sapiens, provide the pre-formated interaction table! The available database is only for humans.")))),
                                                                                             column(3,span(selectInput(width = "100%",inputId="Use_precompiled_database1",choices=c("No"="0","Yes"="1"), selected='1', label=tooltip(trigger=list("Use Pre-compiled Interaction file",icon("info-circle")),"Choose if you would like to use a pre-compiled interaction file, formated in accordance with the instructions. A default interaction file is available. For other organisms than humans, upload a pre-compiled interaction file!")))),
                                                                                           ),
                                                                                           column(offset=1,10,span(actionButton(width = "100%",inputId="submit_Full_analysis",label=tooltip(trigger=list("Submit",icon("info-circle")),"Run the full multi-omics integrative analysis."))))
                                                                                           
                                                                         ),
                                                                         
                                                                         
                                                                         conditionalPanel( condition = "input.Type_of_analysis == 'Run the full pipeline'",
                                                                                           fluidRow(),
                                                                                           fluidRow(
                                                                                             column(offset=1,1,span(width = "100%",actionButton(inputId="AdditionalDatabases",label='', icon = icon('plus'),class = 'btn-xs',style="color: #fff; background-color: #337ab7; border-color: #2e6da4;margin-top:30px;opacity:0.75; input-group{opacity:0.75;}")|>tooltip("Upload interaction databases (Default ones are avialable)"))),
                                                                                             column(9,span(width='100%',conditionalPanel( condition = " input.AdditionalDatabases%2 == 1  & input.Use_precompiled_database1 == '0'",
                                                                                                                                          fluidRow( column(4,span(fileInput(width = "100%",inputId="BioGRID_data_file1", label=tooltip(trigger=list("BioGrid database file",icon("info-circle")),"Upload the BioGrid database file. Leave empty for default version 4.4.218")), style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                                                                                    column(4,span(fileInput(width = "100%",inputId="STRING_data_file1", label=tooltip(trigger=list("STRING database file",icon("info-circle")),"Upload the STRING database file. Leave empty for default version 11.5")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                                                                                    column(4,span(fileInput(width = "100%",inputId="IntAct_data_file1", label=tooltip(trigger=list("IntAct database file",icon("info-circle")),"Upload the IntAct database file. Leave empty for default release 245")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                                                                          )
                                                                                             ),
                                                                                             conditionalPanel( condition = "input.AdditionalDatabases%2 == 1 & input.Use_precompiled_database1 == '1'",
                                                                                                               fluidRow( column(12,span(fileInput(width = "100%",inputId="LookUp_table_file1", label=tooltip(trigger=list("Interaction table file",icon("info-circle")),"Upload the pre-formatted interaction database file. Leave empty for default")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                                               )
                                                                                             )
                                                                                             )
                                                                                             )
                                                                                           ),
                                                                                           br(),
                                                                                           br(),
                                                                                           fluidRow(
                                                                                             column(offset=1,5,span(textInput(width = "100%",inputId="Client_email", label=tooltip(trigger=list("Email address (Optional)",icon("info-circle")),"If you would like to be notified when the analysis is completed, you can provide an email address."))))
                                                                                           )
                                                                                           
                                                                         ),
                                                                         
                                                                         conditionalPanel( condition = "input.Type_of_analysis == 'Custom algorithms'",
                                                                                           fluidRow(
                                                                                             column(offset=1,3,span(textInput(width = "100%",inputId="edge_file_path",value = "edge_files_PPINetworks/Symbol", label=tooltip(trigger=list("Edge file path",icon("info-circle")),"Provide the path where the network edges are saved.")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                             column(3,span(textInput(width = "100%",inputId="Monet_method_string",value = "--method=M1 --avgk=10 --linksdir=undirected", label=tooltip(trigger=list("MONET method",icon("info-circle")),"Provide the MONET analysis method string")))),
                                                                                             column(2,span(textInput(width = "100%",inputId="tmp_bin_folder", label=tooltip(trigger=list("Temporary folder",icon("info-circle")),"Provide a temporary bin server. Leave empty for default")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                             column(2,span(textInput(width = "100%",inputId="monet_path", label=tooltip(trigger=list("MONET path",icon("info-circle")),"Provide the monet executable path. Developer mode")),style="opacity:0.5; input-group{opacity:0.7;}")),
                                                                                           ),
                                                                                           fluidRow(
                                                                                             column(offset=1,10,span(actionButton(width = "100%",inputId="submit_MONET_analysis",label=tooltip(trigger=list("Submit",icon("info-circle")),"Run the MONET decomposition analysis."))))
                                                                                           ),
                                                                                           
                                                                                           fluidRow(
                                                                                             column(offset=1,3,span(textInput(width = "100%",inputId="CPDB_databases", value = "Reactome", label=tooltip(trigger=list("CPDB databases",icon("info-circle")),"Provide the names of the databases against which the enrichment is to be performed. Example: Reactome, Wikipathways")))),
                                                                                             column(4,span(fileInput(width = "100%",inputId="MONET_background_file", label=tooltip(trigger=list("MONET background file",icon("info-circle")),"Provide the MONET background file. Leave empty for default")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                             column(3,span(fileInput(width = "100%",inputId="CPDB_database_file", label=tooltip(trigger=list("CPDB database file",icon("info-circle")),"Provide the CPDB database file. Leave empty for default")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                           ),
                                                                                           fluidRow(
                                                                                             column(offset=1,10,span(actionButton(width = "100%",inputId="submit_MONET_pathways",label=tooltip(trigger=list("Submit",icon("info-circle")),"Extract the pathways associated with each network submodule."))))
                                                                                             
                                                                                           ),
                                                                                           
                                                                                           fluidRow(
                                                                                             column(offset=1,2,span(textInput(width = "100%",inputId="files_edges_path", label=tooltip(trigger=list("Edge files directory",icon("info-circle")),"Provide the edges' directory. Leave empty for default")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                             column(2,span(fileInput(width = "100%",inputId="TF_Database", label=tooltip(trigger=list("TF dataset",icon("info-circle")),"Provide the transcription factors dataset file. Leave empty for default")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                             column(2,span(fileInput(width = "100%",inputId="centralities_file", label=tooltip(trigger=list("Centralities file",icon("info-circle")),"Provide the network centralities file. Leave empty for default")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                             column(2,span(fileInput(width = "100%",inputId="Kinome_database", label=tooltip(trigger=list("Kinome Dataset",icon("info-circle")),"Provide the Kinome dataset file. Leave empty for default")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                             column(2,span(textInput(width = "100%",inputId="file_extension", label=tooltip(trigger=list("File ending",icon("info-circle")),"Provide the file ending specific for a type of analysis. Default to 'Total'. Leave empty for default")),style="opacity:0.5; input-group{opacity:0.5;}")),
                                                                                           ),
                                                                                           fluidRow(
                                                                                             column(offset=1,10,span(actionButton(width = "100%",inputId="submit_Circos_and_auxiliary",label=tooltip(trigger=list("Submit",icon("info-circle")),"Create the summary plots and report."))))
                                                                                           )
                                                                         )
                                                       ),
                                                       conditionalPanel( condition = "input.Type_of_analysis == 'Results download'",
                                                                         fluidRow(
                                                                           column(3,span(textAreaInput(width = "100%",inputId="Results_index_save", label=tooltip(trigger=list("Results directory index",icon("info-circle")),"The directory where the NOODAI analysis results are found. The index is provided after the input data is submitted for the analysis")))),
                                                                           column(4,div(width = "100%",style = "margin-top: 35px;",downloadButton(outputId="downloadData",label="Save"),span(`data-toggle` = "tooltip", `data-placement` = "right",title = "Download the results",icon("info-circle")))),
                                                                           column(offset = 3,2,div(style = "margin-top: 30px;",span(downloadButton(outputId="Demo_Reults_download",icon=NULL,label=tooltip(trigger=list("Demo Results",icon("cloud-arrow-down")),"Download the pre-compiled results for the Demo input data")))))
                                                                           
                                                                         )
                                                       )
                                                       
                                                       
                                                     ),
                                                     conditionalPanel( condition = "output.showRes",
                                                                       fluidRow(
                                                                         column(12,align="center",textOutput("text_output_string")),
                                                                         column(12,align="center",span(textOutput("verb"),style="font-weight:bold;")),
                                                                         column(12,align="center",textOutput("text_output_string2")),
                                                                         column(12,align="center",textOutput("text_output_string3"))
                                                                       )
                                                     ),
                                                     
                                                     conditionalPanel( condition = "output.showResCentralities",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("CentralitiesOut1")),
                                                                         column(12,align="center",verbatimTextOutput("CentralitiesOut2"))
                                                                       )
                                                     ),
                                                     
                                                     conditionalPanel( condition = "output.showResMONET",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("MONETOut1")),
                                                                         column(12,align="center",verbatimTextOutput("MONETOut2"))
                                                                       )
                                                     ),
                                                     
                                                     conditionalPanel( condition = "output.showResPathways",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("PathwaysOut1")),
                                                                         column(12,align="center",verbatimTextOutput("PathwaysOut2"))
                                                                       )
                                                     ),
                                                     
                                                     conditionalPanel( condition = "output.showResImages",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("ImagesOut1")),
                                                                         column(12,align="center",verbatimTextOutput("ImagesOut2"))
                                                                       )
                                                     ),
                                                     
                                                     conditionalPanel( condition = "output.showErrorFormat",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("showErrorFormat1"))
                                                                       )
                                                     ),
                                                     
                                                     conditionalPanel( condition = "output.showErrorDataset",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("showErrorDataset1"))
                                                                       )
                                                     ),
                                                     
                                                     conditionalPanel( condition = "output.showErrorSample",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("showErrorSample1"))
                                                                       )
                                                     ),
                                                     
                                                     conditionalPanel( condition = "output.showErrorContrasts",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("showErrorContrasts1"))
                                                                       )
                                                     ),
                                                     
                                                     conditionalPanel( condition = "output.showErrorBioMart",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("showErrorBioMart1"))
                                                                       )
                                                     ),
                                                     
                                                     conditionalPanel( condition = "output.showErrorMonet",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("showErrorMonet1"))
                                                                       )
                                                     ),
                                                     conditionalPanel( condition = "output.showErrorCPDB",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("showErrorCPDB1"))
                                                                       )
                                                     ),
                                                     conditionalPanel( condition = "output.showErrorBiomartServer",
                                                                       fluidRow(
                                                                         column(12,align="center",verbatimTextOutput("showErrorBiomartServer1"))
                                                                       )
                                                     )
                                                     
                                                     
                                            ),
                                            br(),
                                            fluidRow(
                                              column(offset=10,2,span(textOutput("text_output_version"),style="color:gray; float:right;"))
                                            )
                                            
                                   )),fluidPage( tags$footer("The analysis pipeline is provided for use 'as is' and 'as available' under GNU GPL v3, see <https://www.gnu.org/licenses/>. The authors make no representations or warranties of any kind, expressed or implied, regarding the performance, reliability, or suitability of the analysis.\n By using this website, you agree that the authors shall not be responsible for any data loss, corruption, or unauthorized access to your data. It is your responsibility to implement appropriate backup and security measures to safeguard your data. NOODAI Copyright (C) 2024, Empa, Tiberiu Totu.", class = "footer")
                                                 
                                   )))



server <- function(input, output, session) {
  
  output$text_output_version <- renderText({"version 0.1.0"})
  
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
  file_DEA_names5 <- reactiveVal(value=NULL)
  
  
  #####Start Network analysis#######
  
  Results_index <- reactive({
    x <- input$Results_index
    if(x==""){
      today <- as.character(Sys.Date())
      x <- paste(today, 1, sep = "_")
      dirs <- dir(working_dir)
      if(length((grep(x, dirs)[1]))>0) {
        dirs_today <- dirs[grepl(today, dirs)]
        n <- length(dirs_today)
        x <- paste(x, n + 1,paste0(sample(c(0:9, LETTERS), 6, T),collapse=''), sep = "_")
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
  
  
  BioGRID_data_file <- reactive({
    if(length(input$BioGRID_data_file)>0){
      x <- input$BioGRID_data_file$datapath
      return(x)
    }
  })
  
  BioGRID_data_file1 <- reactive({
    if(length(input$BioGRID_data_file1)>0){
      x <- input$BioGRID_data_file1$datapath
      return(x)
    }
  })
  
  STRING_data_file <- reactive({
    if(length(input$STRING_data_file)>0){
      x <- input$STRING_data_file$datapath
      return(x)
    }
  })
  
  STRING_data_file1 <- reactive({
    if(length(input$STRING_data_file1)>0){
      x <- input$STRING_data_file1$datapath
      return(x)
    }
  })
  
  IntAct_data_file <- reactive({
    if(length(input$IntAct_data_file)>0){
      x <- input$IntAct_data_file$datapath
      return(x)
    }
  })
  
  IntAct_data_file1 <- reactive({
    if(length(input$IntAct_data_file1)>0){
      x <- input$IntAct_data_file1$datapath
      return(x)
    }
  })
  
  file_DEA_names1 <- reactive({
    if(length(input$file_DEA_names1)>0){
      if(length(working_dir)!=0){
        tmp_file_loc <- paste0(working_dir,"/Results_",Results_index(),"/tmp_FilesIni")
        print(tmp_file_loc)
        unlink("tmp_file_loc", recursive=TRUE)        
        dir.create(tmp_file_loc)
        x <- archive_extract(input$file_DEA_names1$datapath,dir = tmp_file_loc)
        return(paste0(tmp_file_loc,"/",x))
      }
    }
  })
  
  file_DEA_names2 <- reactive({
    if(length(input$file_DEA_names2)>0){
      if(length(working_dir)!=0){
        tmp_file_loc <- paste0(working_dir,"/Results_",Results_index(),"/tmp_FilesIni")
        unlink("tmp_file_loc", recursive=TRUE)        
        dir.create(tmp_file_loc)
        x <- archive_extract(input$file_DEA_names2$datapath,dir = tmp_file_loc)
        return(paste0(tmp_file_loc,"/",x))
      }
    }
  })
  
  
  file_DEA_names3 <- reactive({
    if(input$Demo_values!=1){
      if(length(input$file_DEA_names3)>0){
        if(length(working_dir)!=0){
          tmp_file_loc <- paste0(working_dir,"/Results_",Results_index(),"/tmp_FilesIni")
          unlink("tmp_file_loc", recursive=TRUE)        
          dir.create(tmp_file_loc)
          x <- archive_extract(input$file_DEA_names3$datapath,dir = tmp_file_loc)
          return(paste0(tmp_file_loc,"/",x))
        }
      }
    }
  })
  
  LookUp_table_file <- reactive({
    if(length(input$LookUp_table_file)>0){
      x <- input$LookUp_table_file$datapath
      return(x)
    }
  })
  
  LookUp_table_file1 <- reactive({
    if(length(input$LookUp_table_file1)>0){
      x <- input$LookUp_table_file1$datapath
      return(x)
    }
  })
  
  
  phenotype_names <- reactive({
    if(length(input$phenotype_names)>0){
      x <- input$phenotype_names
      x <- gsub(" ", "", x, fixed = TRUE)
      x <- strsplit(x,',')
      x <- unlist(x)
      return(x)
    }
  })
  
  phenotype_comparison <- reactive({
    if(length(input$phenotype_comparison)>0){
      x <- input$phenotype_comparison
      x <- gsub(" ", "", x, fixed = TRUE)
      x <- strsplit(x,',')
      x <- unlist(x)
      return(x)
    }
  })
  
  phenotype_names1 <- reactive({
    if(length(input$phenotype_names1)>0){
      x <- input$phenotype_names1
      x <- gsub(" ", "", x, fixed = TRUE)
      x <- strsplit(x,',')
      x <- unlist(x)
      return(x)
    }
  })
  
  phenotype_comparison1 <- reactive({
    if(length(input$phenotype_comparison1)>0){
      x <- input$phenotype_comparison1
      x <- gsub(" ", "", x, fixed = TRUE)
      x <- strsplit(x,',')
      x <- unlist(x)
      return(x)
    }
  })
  
  splicing_file_name <- reactive({
    if(length(input$splicing_file_name)>0){
      x <- input$splicing_file_name
      x <- gsub(" ", "", x, fixed = TRUE)
      return(x)
    }
  })
  
  
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
    if(length(working_dir)!=0){
      dir.create(paste0(working_dir,"/Results_",Results_index()))
      tmp_file_loc <- paste0(working_dir,"/Results_",Results_index(),"/tmp_FilesIni")
      unlink("tmp_file_loc", recursive=TRUE)        
      dir.create(tmp_file_loc)
      x <- archive_extract(paste0(working_dir,"/Data_Omics_Demo.7z"),dir = tmp_file_loc)
      file_DEA_names5(paste0(tmp_file_loc,"/",x))
    }
    
  })
  
  
  
  
  
  
  observeEvent(input$submit_Network_analysis, {
    
    shinyjs::addClass(id = "submit_Network_analysis", class = "loading dots")
    shinyjs::disable("submit_Full_analysis")
    shinyjs::disable("submit_Network_analysis")
    shinyjs::disable("submit_MONET_pathways")
    shinyjs::disable("submit_Circos_and_auxiliary")
    shinyjs::disable("submit_MONET_analysis")
    
    
    if(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))<300000){
      showMemFull(TRUE)
    }else{
      
      if(length(working_dir)!=0){
        
        if(biomart_check()==0){
          showVal(FALSE)
          showErrorBiomartServerMessage(TRUE)
        }else{
          
          dir.create(paste0(working_dir,"/Results_",Results_index()))
          #con <- file(paste0(working_dir,"/Results_",Results_index(),"/log.txt"))
          #sink(paste0(working_dir,"/Results_",Results_index(),"/log.txt"), append=TRUE,type = c("output", "message"))
          
          token <- UUIDgenerate()
          message(paste0("running task for token: ", token))
          
          if(is.null(jobs[[token]])){
            if(Type_of_analysis()=='Run the full pipeline'){
              file_DEA_names <- file_DEA_names3()
              Use_precompiled_databaseVal <- Use_precompiled_database1()
              LookUp_table_fileVal <- LookUp_table_file1()
              BioGRID_data_fileVal <- BioGRID_data_file1()
              STRING_data_fileVal <- STRING_data_file1()
              IntAct_data_fileVal <- IntAct_data_file1()
              phenotype_names <- phenotype_names()
              phenotype_comparison <- phenotype_comparison()
              BioMart_Dataset <- BioMart_Dataset1()
            }else{
              Use_precompiled_databaseVal <- Use_precompiled_database()
              LookUp_table_fileVal <- LookUp_table_file()
              BioGRID_data_fileVal <- BioGRID_data_file()
              STRING_data_fileVal <- STRING_data_file()
              IntAct_data_fileVal <- IntAct_data_file()
              phenotype_names <- phenotype_names1()
              phenotype_comparison <- phenotype_comparison1()
              BioMart_Dataset <- BioMart_Dataset()
              if(Use_precompiled_database() == '0'){file_DEA_names <- file_DEA_names1()}
              if(Use_precompiled_database() == '1'){file_DEA_names <- file_DEA_names2()}
            }
            
            
            analyses_flag <- 1
            
            
            for (i in 1:length(phenotype_comparison)){
              if(length(unlist(sapply(phenotype_names,FUN=function(x){grep(x,phenotype_comparison[i])})))!=2){
                showVal(FALSE)
                showErrorContrastsMessage(TRUE)
                showErrorSampleMessage(TRUE)
                break()
              }
            }
            
            if(length(file_DEA_names)==0){
              analyses_flag<-0
              showVal(FALSE)
              showErrorDatasetMessage(TRUE)
              
            }
            if(length(phenotype_comparison)==0){
              analyses_flag<-0
              showVal(FALSE)
              showErrorContrastsMessage(TRUE)
              
            }
            if(length(phenotype_names)<2){
              analyses_flag<-0
              showVal(FALSE)
              showErrorSampleMessage(TRUE)
              
            }
            
            if(analyses_flag){
              for (j in 1:length(file_DEA_names)){
                file_Sheets1 <- readxl::excel_sheets(file_DEA_names[j])
                if(length(which(phenotype_comparison %in% file_Sheets1))!=length(phenotype_comparison)){
                  analyses_flag<-0
                  showVal(FALSE)
                  showErrorContrastsMessage(TRUE)
                  showErrorDatasetMessage(TRUE)
                  break()
                }
                for (i in 1:length(phenotype_comparison)){
                  aux <- as.data.frame(readxl::read_excel(file_DEA_names[j],sheet = phenotype_comparison[i]))
                  #Get the corresponding entrez ids from the uniprot ones#######################
                  if(length(aux)>0){
                    aux1 <- aux$Protein[which(aux$Protein %in% uniprot_ids_maps$Protein)]
                    if(length(aux1)<(length(aux$Protein)/3)){
                      showVal(FALSE)
                      showErrorFormatMessage(TRUE)
                      analyses_flag<-0
                      break()
                    }
                  }
                }
              }
            }
            
            
            if(length(which(BioMart_Dataset %in% biomart_available_datasets$Dataset))==0){
              showVal(FALSE)
              showErrorBioMartMessage(TRUE)
              analyses_flag<-0
            }
            
            
            if(analyses_flag){
              source(paste0(working_dir,"/Integrative_analysis.R"), local = TRUE,echo=TRUE, max.deparse.length=100000)
              jobs[[token]] <<- callr::r_bg(Integrative_Network_analysis, supervise = FALSE, stdout = paste0(working_dir,"/Results_",Results_index(),"/out_Centralities.txt"), package=TRUE, stderr = paste0(working_dir,"/Results_",Results_index(),"/log_Centralities.txt"),
                                            args = list(working_dir=working_dir,BioGRID_data_file=BioGRID_data_fileVal,STRING_data_file=STRING_data_fileVal,
                                                        IntAct_data_file=IntAct_data_fileVal,file_DEA_names=file_DEA_names,
                                                        phenotype_names=phenotype_names,phenotype_comparison=phenotype_comparison,
                                                        splicing_file_name=splicing_file_name(),Use_precompiled_database=Use_precompiled_databaseVal,
                                                        LookUp_table_file=LookUp_table_fileVal,Results_index=Results_index(),BioMart_Dataset=BioMart_Dataset))
              
              showVal(FALSE)
              showResCentralities1(TRUE)
            }
          }
          
        }
      }
    }
    
  })
  
  #####End Network analysis#######
  
  
  
  
  #####Start MONET analysis#######
  
  
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
    
    if(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))<300000){
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
            jobs[[token]] <<- callr::r_bg(MONET_analysis, supervise = FALSE, stdout = paste0(working_dir,"/Results_",Results_index(),"/out_MONET.txt"), package=TRUE, stderr = paste0(working_dir,"/Results_",Results_index(),"/log_MONET.txt"),
                                          args = list(working_dir=working_dir,edge_file_path=edge_file_path(),monet_path=monet_path(),
                                                      Monet_method_string=Monet_method_string(),tmp_bin_folder=tmp_bin_folder(),
                                                      Results_index=Results_index()))
            
            
            showVal(FALSE)
            showResMONET1(TRUE)
          }
        }
        
        
      }
    }
    
  })
  
  
  #####End MONET analysis#######
  
  
  #####Start MONET pathways analysis#######
  
  CPDB_databases <- reactive({
    if(length(input$CPDB_databases)>0){
      x <- input$CPDB_databases
      x <- strsplit(x,',')
      x <- unlist(x)
      return(x)
    }
  })
  
  CPDB_databasesFULL <- reactive({
    if(length(input$CPDB_databasesFULL)>0){
      x <- input$CPDB_databasesFULL
      x <- strsplit(x,',')
      x <- unlist(x)
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
    
    if(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))<300000){
      showMemFull(TRUE)
    }else{
      
      if(length(working_dir)!=0){
        
        if(biomart_check()==0){
          showVal(FALSE)
          showErrorBiomartServerMessage(TRUE)
        }else{
          
          dir.create(paste0(working_dir,"/Results_",Results_index()))
          #con <- file(paste0(working_dir,"/Results_",Results_index(),"/log.txt"))
          #sink(paste0(working_dir,"/Results_",Results_index(),"/log.txt"), append=TRUE,type = c("output", "message"))
          
          
          token <- UUIDgenerate()
          message(paste0("running task for token: ", token))
          
          if(is.null(jobs[[token]])){
            
            if(Type_of_analysis()=='Run the full pipeline'){
              phenotype_names <- phenotype_names()
              phenotype_comparison <- phenotype_comparison()
            }else{
              phenotype_names <- phenotype_names1()
              phenotype_comparison <- phenotype_comparison1()
            }
            
            analyses_flag<-1
            if(length(which(BioMart_Dataset() %in% biomart_available_datasets$Dataset))==0){
              showVal(FALSE)
              showErrorBioMartMessage(TRUE)
              analyses_flag<-0
            }
            
            if(length(which(CPDB_databases() %in% CPDB_databases_list$Database))==0){
              showVal(FALSE)
              showErrorCPDBMessage(TRUE)
              analyses_flag<-0
            }
            
            if(analyses_flag){
              
              source(paste0(working_dir,"/Integrative_analysis.R"), local = TRUE,echo=TRUE, max.deparse.length=100000)
              jobs[[token]] <<- callr::r_bg(MONET_pathways, supervise = FALSE, stdout = paste0(working_dir,"/Results_",Results_index(),"/out_Pathways.txt"), package=TRUE, stderr = paste0(working_dir,"/Results_",Results_index(),"/log_Pathways.txt"),
                                            args = list(working_dir=working_dir,CPDB_databases=CPDB_databases(),MONET_background_file=MONET_background_file(),
                                                        phenotype_names=phenotype_names,phenotype_comparison=phenotype_comparison,
                                                        CPDB_database_file=CPDB_database_file(),Results_index=Results_index(),BioMart_Dataset=BioMart_Dataset()))
              
              showVal(FALSE)
              showResPathways1(TRUE)
              
            }
          }
          
        }
      }
    }
    
  })
  
  
  #####End MONET pathways analysis#######
  
  
  #####Start MONET images generation#######
  
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
    
    if(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))<300000){
      showMemFull(TRUE)
    }else{
      
      if(length(working_dir)!=0){
        
        if(biomart_check()==0){
          showVal(FALSE)
          showErrorBiomartServerMessage(TRUE)
        }else{
          
          dir.create(paste0(working_dir,"/Results_",Results_index()))
          #con <- file(paste0(working_dir,"/Results_",Results_index(),"/log.txt"))
          # sink(paste0(working_dir,"/Results_",Results_index(),"/log.txt"), append=TRUE,type = c("output", "message"))
          
          
          token <- UUIDgenerate()
          message(paste0("running task for token: ", token))
          
          if(is.null(jobs[[token]])){
            
            if(Type_of_analysis()=='Run the full pipeline'){
              phenotype_names <- phenotype_names()
              phenotype_comparison <- phenotype_comparison()
            }else{
              phenotype_names <- phenotype_names1()
              phenotype_comparison <- phenotype_comparison1()
            }
            
            analyses_flag<-1
            if(length(which(BioMart_Dataset() %in% biomart_available_datasets$Dataset))==0){
              showVal(FALSE)
              showErrorBioMartMessage(TRUE)
              analyses_flag<-0
            }
            
            if(analyses_flag){
              source(paste0(working_dir,"/Integrative_analysis.R"), local = TRUE,echo=TRUE, max.deparse.length=100000)
              jobs[[token]] <<- callr::r_bg(Circos_and_auxiliary, supervise = FALSE, stdout = paste0(working_dir,"/Results_",Results_index(),"/out_Images.txt"), package=TRUE, stderr = paste0(working_dir,"/Results_",Results_index(),"/log_Images.txt"),
                                            args = list(working_dir=working_dir,phenotype_names=phenotype_names,phenotype_comparison=phenotype_comparison,
                                                        files_edges_path=files_edges_path(),centralities_file=centralities_file(),
                                                        TF_Database=TF_Database(),file_extension=file_extension(),Results_index=Results_index(), Kinome_database=Kinome_database(),BioMart_Dataset=BioMart_Dataset()))
              
              showVal(FALSE)
              showResImages1(TRUE)
              
            }
          }
          
        }
      }
    }
    
  })
  
  #####End MONET images generation#######  
  
  #####Run the full Network analysis####### 
  
  
  
  observeEvent(input$submit_Full_analysis, {
    
    shinyjs::addClass(id = "submit_Full_analysis", class = "loading dots")
    shinyjs::disable("submit_Full_analysis")
    shinyjs::disable("submit_Network_analysis")
    shinyjs::disable("submit_MONET_pathways")
    shinyjs::disable("submit_Circos_and_auxiliary")
    shinyjs::disable("submit_MONET_analysis")
    
    if(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))<300000){
      showMemFull(TRUE)
    }else{
      
      if(length(working_dir)!=0){
        
        if(biomart_check()==0){
          showVal(FALSE)
          showErrorBiomartServerMessage(TRUE)
        }else{
          
          dir.create(paste0(working_dir,"/Results_",Results_index()))
          token <- UUIDgenerate()
          message(paste0("running task for token: ", token))
          # the if statement is to avoid rerunning a job again
          if(is.null(jobs[[token]])){
            # call the job in the background session
            if(Type_of_analysis()=='Run the full pipeline'){
              if(input$Demo_values>=1){file_DEA_names <- file_DEA_names5()}else{file_DEA_names <- file_DEA_names3()}
              Use_precompiled_databaseVal <- Use_precompiled_database1()
              LookUp_table_fileVal <- LookUp_table_file1()
              BioGRID_data_fileVal <- BioGRID_data_file1()
              STRING_data_fileVal <- STRING_data_file1()
              IntAct_data_fileVal <- IntAct_data_file1()
              phenotype_names <- phenotype_names()
              phenotype_comparison <- phenotype_comparison()
              BioMart_Dataset <- BioMart_Dataset1()
            }else{
              Use_precompiled_databaseVal <- Use_precompiled_database()
              LookUp_table_fileVal <- LookUp_table_file()
              BioGRID_data_fileVal <- BioGRID_data_file()
              STRING_data_fileVal <- STRING_data_file()
              IntAct_data_fileVal <- IntAct_data_file()
              phenotype_names <- phenotype_names1()
              phenotype_comparison <- phenotype_comparison1()
              BioMart_Dataset <- BioMart_Dataset()
              if(Use_precompiled_database() == '0'){file_DEA_names <- file_DEA_names1()}
              if(Use_precompiled_database() == '1'){file_DEA_names <- file_DEA_names2()}
            }
            
            analyses_flag <- 1
            
            
            for (i in 1:length(phenotype_comparison)){
              if(length(unlist(sapply(phenotype_names,FUN=function(x){grep(x,phenotype_comparison[i])})))!=2){
                showVal(FALSE)
                showErrorContrastsMessage(TRUE)
                showErrorSampleMessage(TRUE)
                break()
              }
            }
            
            if(length(file_DEA_names)==0){
              analyses_flag<-0
              showVal(FALSE)
              showErrorDatasetMessage(TRUE)
              
            }
            if(length(phenotype_comparison)==0){
              analyses_flag<-0
              showVal(FALSE)
              showErrorContrastsMessage(TRUE)
              
            }
            if(length(phenotype_names)<2){
              analyses_flag<-0
              showVal(FALSE)
              showErrorSampleMessage(TRUE)
              
            }
            
            if(analyses_flag){
              for (j in 1:length(file_DEA_names)){
                file_Sheets1 <- readxl::excel_sheets(file_DEA_names[j])
                if(length(which(phenotype_comparison %in% file_Sheets1))!=length(phenotype_comparison)){
                  analyses_flag<-0
                  showVal(FALSE)
                  showErrorContrastsMessage(TRUE)
                  showErrorDatasetMessage(TRUE)
                  break()
                }
                for (i in 1:length(phenotype_comparison)){
                  aux <- as.data.frame(readxl::read_excel(file_DEA_names[j],sheet = phenotype_comparison[i]))
                  #Get the corresponding entrez ids from the uniprot ones#######################
                  if(length(aux)>0){
                    # if(length(aux)>40){aux<-aux[c(1:40),]}
                    # aux1 <- as.character(unique(get_protein_idMART(aux$Protein,mart)))
                    aux1 <- aux$Protein[which(aux$Protein %in% uniprot_ids_maps$Protein)]
                    if(length(aux1)<(length(aux$Protein)/3)){
                      showVal(FALSE)
                      showErrorFormatMessage(TRUE)
                      analyses_flag<-0
                      break()
                    }
                  }
                }
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
            
            if(analyses_flag){
              jobs[[token]] <<- callr::r_bg(run_token, supervise = FALSE, stdout = paste0(working_dir,"/Results_",Results_index(),"/out.txt"), package=TRUE, stderr = paste0(working_dir,"/Results_",Results_index(),"/log.txt"),
                                            args = list(token=token,working_dir=working_dir,BioGRID_data_file=BioGRID_data_fileVal,STRING_data_file=STRING_data_fileVal,IntAct_data_file=IntAct_data_fileVal,file_DEA_names=file_DEA_names,
                                                        phenotype_names=phenotype_names,phenotype_comparison=phenotype_comparison,splicing_file_name=splicing_file_name(),
                                                        Use_precompiled_database=Use_precompiled_databaseVal,LookUp_table_file=LookUp_table_fileVal,Results_index=Results_index(),
                                                        edge_file_path=edge_file_path(),monet_path=monet_path(),Monet_method_string=Monet_method_stringFULL(),tmp_bin_folder=tmp_bin_folder(),
                                                        CPDB_databases=CPDB_databasesFULL(),MONET_background_file=MONET_background_file(),CPDB_database_file=CPDB_database_file(),
                                                        files_edges_path=files_edges_path(),centralities_file=centralities_file(),TF_Database=TF_Database(),file_extension=file_extension(),Kinome_database=Kinome_database(),BioMart_Dataset=BioMart_Dataset,Client_email=Client_email()))
              
              showVal(FALSE)
              showRes1(TRUE)
            }
          }
        }
      }
    }
    
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
    
    
    useReactiveDownload <- function(input,output,session,df){
      
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
    
    output$downloadData <- useReactiveDownload(input,output,session,df = Results_index_save())
    
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
  
  output$verb <- renderText({ paste0("Your results folder is: ", Results_index()) })
  output$text_output_string <- renderText({ "Thank you for submitting an analysis request. Please come back in at least two hours in order to download the results. If you provided an email address you will be notified!"})
  output$text_output_string2 <- renderText({ "Please save it as it is the only way to trace back your results folder due to security. Your analysis will be automatically removed from the server after 5 days."})
  
  
  output$showResCentralities <- reactive({
    return(showResCentralities1() == "TRUE")
  })
  outputOptions(output, "showResCentralities", suspendWhenHidden = FALSE)
  
  output$CentralitiesOut1 <- renderText({ paste0("Your results folder is: ", Results_index()) })
  output$CentralitiesOut2 <- renderText({ "Thank you for submitting an analysis request. Please come back in at least one hour in order to download the results."})
  
  output$showResMONET <- reactive({
    return(showResMONET1() == "TRUE")
  })
  outputOptions(output, "showResMONET", suspendWhenHidden = FALSE)
  
  output$MONETOut1 <- renderText({ paste0("Your results folder is: ", Results_index()) })
  output$MONETOut2 <- renderText({ "Thank you for submitting an analysis request. Please come back in at least one hour in order to download the results."})
  
  output$showResPathways <- reactive({
    return(showResPathways1() == "TRUE")
  })
  outputOptions(output, "showResPathways", suspendWhenHidden = FALSE)
  
  output$PathwaysOut1 <- renderText({ paste0("Your results folder is: ", Results_index()) })
  output$PathwaysOut2 <- renderText({ "Thank you for submitting an analysis request. Please come back in at least one hour in order to download the results."})
  
  
  output$showResImages <- reactive({
    return(showResImages1() == "TRUE")
  })
  outputOptions(output, "showResImages", suspendWhenHidden = FALSE)
  
  output$ImagesOut1 <- renderText({ paste0("Your results folder is: ", Results_index()) })
  output$ImagesOut2 <- renderText({ "Thank you for submitting an analysis request. Please come back in at least one hour in order to download the results."})
  
  
  output$showErrorFormat <- reactive({
    return(showErrorFormatMessage() == "TRUE")
  })
  outputOptions(output, "showErrorFormat", suspendWhenHidden = FALSE)
  
  output$showErrorFormat1 <- renderText({ "Your Input data does not have the proper format. The names must be Uniprot Ids and the column header named Protein" })
  
  
  output$showErrorDataset <- reactive({
    return(showErrorDatasetMessage() == "TRUE")
  })
  outputOptions(output, "showErrorDataset", suspendWhenHidden = FALSE)
  
  output$showErrorDataset1 <- renderText({ "Your input data may be incorrect. Please provide an archive with each omics measurement in an Excel table where each sheet correspond to a contrast." })
  
  output$showErrorSample <- reactive({
    return(showErrorSampleMessage() == "TRUE")
  })
  outputOptions(output, "showErrorSample", suspendWhenHidden = FALSE)
  
  output$showErrorSample1 <- renderText({ "The provided samples names may be incorrect. Please use the same names that are included in the contrasts and in the Excel sheets." })
  
  
  output$showErrorContrasts <- reactive({
    return(showErrorContrastsMessage() == "TRUE")
  })
  outputOptions(output, "showErrorContrasts", suspendWhenHidden = FALSE)
  
  output$showErrorContrasts1 <- renderText({ "The provided contrasts may be incorrect. Please format the contrasts to include the samples names and be the same as the Excel sheets." })
  
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
  
  
  #####End Output#######
  
}



shinyApp(ui, server)


