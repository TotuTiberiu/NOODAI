# install_packages.R — NOODAI

# Check and install BiocManager if not already available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
bioc_packages <- c(
  "biomaRt", "org.Hs.eg.db", "AnnotationDbi"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Install CRAN packages
cran_packages <- c(
  "archive", "uuid", "shiny", "reshape2", "ggplot2", "ggpubr",
  "EnvStats", "moments", "bslib", "shinyWidgets", "readr",
  "RColorBrewer", "circlize", "ComplexHeatmap", "gridBase",
  "gtools", "readxl", "centiserve", "CINNA", "igraph",
  "tidygraph", "stringr", "parallel", "RCy3", "msigdbr",
  "rbioapi", "clusterProfiler", "ggchicklet", "openxlsx", "shinyjs"
)

# Install CRAN packages using install.packages()
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# CINNA must be installed via devtools due to version issues
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_version("CINNA", version = "1.2.0", repos = "http://cran.us.r-project.org")

message("All required packages for NOODAI are installed.")