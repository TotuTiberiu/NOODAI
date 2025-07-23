# NOODAI – Network-Oriented multi-Omics Data Analysis and Integration

**NOODAI** is a modular, web-based RShiny platform for integrative multi-omics analysis through network construction, analysis, decomposition, and enrichment. It enables researchers to explore, visualize, and interpret individual features as well as signaling routes supported by the combined input omics layers, such as transcriptomics, proteomics, and metabolomics.

This repository provides the local version of the application running at [https://omics-oracle.com](https://omics-oracle.com).

---

## Features

- Integrate multi-omics layers (e.g. transcriptomics, proteomics, phosphoproteomics and small molecule datasets)
- Construct regulatory networks using BioGRID, STRING, and Intact or use custom interaction files
- Perform network decomposition with MONET
- Conduct pathway enrichment on decomposed MONET modules
- Generate Cytoscape-compatible networks, circus and enrichment plots.
- Generate PDF reports

- Pre-set for 13 species:

|                             |                             |
|-----------------------------|-----------------------------|
| *Homo sapiens*              | *Rattus norvegicus*         |
| *Mus musculus*              | *Saccharomyces cerevisiae*  |
| *Bos taurus*                | *Xenopus tropicalis*        |
| *Caenorhabditis elegans*    | *Sus scrofa*                |
| *Canis familiaris*          | *Schizosaccharomyces pombe* |
| *Danio rerio*               |                             |
| *Drosophila melanogaster*   |                             |
| *Gallus gallus*             |                             |


---

## Repository Structure

```
|                                            |             |
|--------------------------------------------|-------------|
| `Databases/`                                | Reference files and Cytoscape style configuration |
| `ccRCC_inputData/`                          | Input data from the case study: [DOI:10.1101/2025.01.31.635927 (https://doi.org/10.1101/2025.01.31.635927) |
| `www/`                                      | App UI assets (CSS, JavaScript, images) |
| `app.R`                                     | Main Shiny application script |
| `install_NOODAI_packages.R`                       | List of required R packages for the full NOODAI functionality |
| `Circos_plots.R`                            | Script for generating circos plots at the level of MONET modules intersection |
| `Cytoscape_Network_Generation.R`            | Cytoscape network construction from module-specific edges |
| `Data_Omics_Demo.7z`                         | Full-size omics demo data archive |
| `Demo_data.zip`                              | Tutorial dataset from: [DOI:10.1093/narmme/ugaf013](https://doi.org/10.1093/narmme/ugaf013) |
| `Extract_TopCentralitiesInModuels.R`        | Extract top central genes or proteins per module |
| `Extract_edge_files_for_clusters.R`         | Prepares edge lists for per-module network decomposition |
| `Integrative_analysis.R`                    | Master script for cross-omics integration and module assignment |
| `LICENSE`                                   | License file (GPL-3.0) |
| `MONET_pathways_extraction_CPDB.R`          | CPDB-based pathway enrichment script using MONET modules |
| `Monet_cluster_pathways_image_joint.R`      | Pathway enrichment image export for MONET clusters |
| `Network_analysis_BioGRID_STRING_IntAct.R`  | Network construction and centrality metric computation |
| `Networks_Analyisis_results_PDF_generation.Rmd` | RMarkdown report generation for PDF summary of analysis |
| `Runtime_test_datasets.7z`                  | Comprehensive multi-omics runtime test dataset |
| `Tools-runtime.zip`                         | Runtime test for competitor tools |
| `Uplaoded_Data_Archive.zip`                 | Demo dataset for direct website upload |
```

---

## How to Deploy Locally

### 1. Prerequisites

- R version 4.3 or higher
- RStudio (recommended)
- Required R packages (see `install_NOODAI_packages.R`)
- [MONET](https://github.com/BergmannLab/MONET) installed and available via command line
- [Cytoscape](https://cytoscape.org/) installed and running in **headless mode**

### 2. Clone This Repository

Open a terminal and run:

```
git clone https://github.com/YourUsername/NOODAI.git
cd NOODAI
```

### 3. Set Working Directory

Edit the first line of `app.R` to match your local directory path:

```r
setwd("/your/full/path/to/NOODAI")
```

Replace the path with your actual folder location.

### 4. Install Required R Packages

```
Rscript install_NOODAI_packages.R
```

### 5. Start Headless Cytoscape

Cytoscape must be running in headless mode before starting NOODAI. From terminal:

```
env JAVA_HOME=$JAVA_HOME JAVA_OPTS='-Djava.awt.headless=true' ./cytoscape.sh -R -Dorg.lwjgl.util.Debug=true -Dorg.lwjgl.util.DebugLoader=true

```

(On Windows, this may be `Cytoscape.bat -R`.)

---

### 6. Launch NOODAI

From within R:

```r
shiny::runApp("app.R")
```

This will open NOODAI in your browser at a local address

---

## Configuration Notes

- MONET must be executable by the R session. If `sudo` is needed, configure securely.
- Cytoscape must be launched before running the app.
- Ensure all directories and files used by the app are readable and writable.

---

