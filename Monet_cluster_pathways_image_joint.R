#The function creates the representation of the top 3 pathways for the top 5 biggest subnetworks (considering the number of elements)

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

MONET_cluster_pathways_joint_image <- function(pathways_dir, file_extension, NOODAI_object) {

  library(ggplot2)
  library(ggchicklet)
  library(openxlsx)
  library(stringr)

  setwd(pathways_dir)

  theme_a <- theme(
    axis.title = element_text(
      color="black",
      size = 22,
      hjust = 0.5,
      vjust = 1
    ),
    axis.line = element_line(
      size = 1,
      colour = "black"
    ),
    axis.text.x = element_text(
      color = "black",
      size = 18
    ),
    axis.text.y = element_text(
      color = "black",
      size = 18
    ),
    plot.title = element_text(
      color = "black",
      face = 'bold.italic',
      size = 18,
      hjust = 0,
      vjust = 1
    ),
    panel.background = element_blank(),
    axis.line.y = element_blank(),
    legend.title = element_text(
      color = "black",
      size = 18,
      face = "bold"
    ),
    legend.text = element_text(
      color = "black",
      size = 18
    ),
    legend.key = element_rect(fill = "whitesmoke"),
    legend.background = element_blank(),
    legend.text.align = 0.5
  )


  files_pathways <- list.files(
    path = getwd(),
    pattern = paste0(file_extension, "_")
  )

  dir.create(paste0(pathways_dir,'/Images/'), recursive = T)

  for (ij in files_pathways) {

    xlsx_sheets <- openxlsx::getSheetNames(ij)

    if(length(xlsx_sheets) > 5) {xlsx_sheets <- xlsx_sheets[1:5]}

    data <- c()
    index <-1
    
    for (jk in xlsx_sheets) {
      
      aux <- openxlsx::read.xlsx(
        ij,
        sheet = jk,
        colNames = TRUE,
        rowNames = TRUE
      )

      if(nrow(aux) > 3) {
        aux <- aux[1:3, c(1, 4, 7)]
      } else {
        aux <- aux[1, c(1, 4, 7)]
      }

      ams <- str_split(aux[, 3], '/')
      ams <- sapply(
        ams,
        FUN = function(x) {
          as.numeric(x[[1]]) / as.numeric(x[[2]])
        }
      )
      aux$Ratio <- ams
      aux$Cluster <- paste0("Cluster", index)
      index <- index + 1
      data <- rbind(data, aux)
      
    }

    x <- data

    colnames(x) <- c('X1', 'X2', 'X3', 'X4', 'X5')

    x1 <- x
    x1$X2 <- as.numeric(x1$X2)


    x1$X2 <- -log10(as.numeric(x1$X2))

    index <- 1 
    unq_cls <- unique(x1$X5)
    for (mm in unq_cls) {
      x1$X5[which(x1$X5 %in% mm)] <- index
      index <- index + 1
    }

    x1$X1 <- as.character(x1$X1)
    x1$X1 <- factor(x1$X1, levels = unique(x1$X1))
    x1$X5 <- as.character(x1$X5)
    x1$X5 <- factor(x1$X5, levels = unique(x1$X5))

    ##

    slot_name = sapply(strsplit(ij, "_"), "[[", length(strsplit(ij, "_")[[1]]) - 1)
    enrichment_data = list(name = x1)
    names(enrichment_data) = slot_name
    NOODAI_object@pathways = c(NOODAI_object@pathways, enrichment_data)

    ##

    g <- ggplot(data = x1, aes(x = nrow(x1):1, y = X2, fill = X5)) + 
      scale_fill_manual(
        values = c(
          "#7196BE",
          "#BD3737",
          "#CD9B1D",
          "#E8768F",
          "#855C5C"
        )
      ) +
      geom_chicklet(radius = grid::unit(7.5, "mm"), width = 0.8) +
      coord_flip() + 
      ylab("-log10(q-value)")

    g$theme <- theme_a
    g$labels$x = element_blank()
    
    g <- g + 
      scale_x_discrete(labels = x1$X1, breaks = nrow(x1):1, limits = factor(1:nrow(x1))) +
      guides(fill = guide_legend(title = "Module"))

    save_name_pdf <- paste0(
      pathways_dir,
      '/Images/',
      str_remove(ij, '.xlsx'),
      '_FDR.pdf'
    )

    pdf(
      file = save_name_pdf,
      width = (trunc(max(nchar(as.character(x1$X1))) / 10) * 20 / 3) / 2,
      length(x1$X1) * 1.5 / 2
    )
    print(g, newpage = FALSE)
    dev.off()

    x1$X1 <- as.character(x1$X1)
    x1$X1 <- factor(x1$X1, levels = unique(x1$X1))
    x1$X5 <- as.character(x1$X5)
    x1$X5 <- factor(x1$X5, levels = unique(x1$X5))

    g <- ggplot(data = x1, aes(x = nrow(x1):1, y = X4, fill = X5)) + 
      scale_fill_manual(
        values = c(
          "#7196BE",
          "#BD3737",
          "#CD9B1D",
          "#E8768F",
          "#855C5C"
        )
      ) +
      geom_chicklet(radius = grid::unit(7.5, "mm"), width = 0.8) +
      coord_flip() +
      ylab("Pathway Hits/Cluster Members")

    g$theme <- theme_a
    g$labels$x = element_blank()

    g <- g + 
      scale_x_discrete(labels = x1$X1, breaks = nrow(x1):1, limits = factor(1:nrow(x1))) +
      guides(fill = guide_legend(title = "Module"))

    save_name_pdf <- paste0(
      pathways_dir,
      '/Images/',
      str_remove(ij, '.xlsx'),
      '_Ratio.pdf'
    )

    pdf(
      file = save_name_pdf,
      width = (trunc(max(nchar(as.character(x1$X1))) / 10) * 20 / 3) / 2,
      length(x1$X1) * 1.5 / 2
    )
    print(g, newpage = FALSE)
    dev.off()

  }

  return(list(1, NOODAI_object))

}
