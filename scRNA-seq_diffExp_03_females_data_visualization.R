### Visualization of Lipidomic-scRNA-seq data integration ###

library(openxlsx)
library(glue)
library(tidyr)
library(ggplot2)
library(viridus)
library(UpSetR)


cell_types = c("Astro", "L2_3_IT", "L4", "L5", "L6", "Lamp5", "Non-neuronal", "Oligo", "Pvalb", "Sncg", "Sst", "Vip")

readDEGs <- function(cellType) {
  
  files <- list.files(path=glue::glue("/Users/karineier/Documents/scRNA-seq/Females/{cellType}"), pattern = "^DEGs*")
  list <- vector(mode = "list", length = length(files))
  lipidIDs <- gsub("DEGs_", "", files)
  lipidIDs <- gsub(".xlsx", "", lipidIDs)
  names(list) = lipidIDs
  
  for (i in lipidIDs) {
    list[[i]] <- openxlsx::read.xlsx(glue::glue("/Users/karineier/Documents/scRNA-seq/Females/{cellType}/DEGs_{i}.xlsx"), rowNames = TRUE)
  }
  
  return(list)
}


diffEdata <- lapply(cell_types, readDEGs)
names(diffEdata) = cell_types

sigLipids = names(diffEdata[[1]])

plotData = tidyr::crossing(cell_types, sigLipids)

plotDataList = lapply(cell_types, function(cellType) {
  plotData = plotData %>% 
    dplyr::filter(cell_types == cellType)
  plotData$numDEGs = sapply(sigLipids, function(lipid) {
    length(which(diffEdata[[cellType]][[lipid]][,4] < 0.1))
  })
})

names(plotDataList) = cell_types

plotData = as.data.frame(do.call("rbind", plotDataList))
plotData$cell_type = rownames(plotData)

plotData_long = tidyr::gather(plotData, Lipid, numDEGs, colnames(plotData[1]):colnames(plotData[ncol(plotData)-1]))

plotData_long$numDEGs = as.numeric(plotData_long$numDEGs)

plotData_long[plotData_long == 0] <- NA

plot = ggplot(plotData_long, aes(x=cell_type, y=Lipid)) +
  geom_point(aes(size = numDEGs, color = numDEGs), shape=19) +
  scale_size_continuous(breaks=c(1:53), range=c(5,20)) +
  geom_text(aes(label=numDEGs), color="white") +
  theme_classic() +
  theme(legend.position="NULL", axis.text.x = element_text(angle=45, hjust=1)) +
  scale_color_viridis(option="plasma") +
  xlab("") +
  ylab("") 

ggsave("Summary_Lipid_DEGs_by_celltype_females.pdf", width=8.5, height=11)

### UpSet plots ###

Astro.list = lapply(sigLipids, function(x) {
  rownames(diffEdata$Astro[[x]])[which(diffEdata$Astro[[x]]$adj.P.Val<0.1)]
})

names(Astro.list) = sigLipids

Astro.upset = fromList(Astro.list)

pdf(file="Astrocytes_upset_plot_females.pdf", width=8.5, height=8.5)
UpSetR::upset(Astro.upset)
dev.off()


