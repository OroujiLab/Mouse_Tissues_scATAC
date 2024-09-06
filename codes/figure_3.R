library(ArchR)
library(Seurat)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(parallel)
library(pheatmap)
library(scran)
library(Signac)
library(sctransform)
library(SeuratObject)
library(dplyr)
library(Cairo)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Gviz)
library(Mus.musculus)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)


Mouse <- readRDS("./outputs/Archr_processed/Save-ArchR-Project.rds")
markersPeaks <- readRDS("./outputs/markerPeaks.rds")

motifsUp <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = Mouse,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.1"
)

saveRDS(motifsUp,"motifsUp.rds")

Markers_GS <- getMarkerFeatures(
  ArchRProj = Mouse, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmap <- plotMarkerHeatmap(
  seMarker = Markers_GS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  transpose = TRUE,
  returnMatrix = FALSE)
heatmap


enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = Mouse,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.1"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = TRUE, returnMatrix = TRUE)
heatmapEM <- heatmapEM[!rownames(heatmapEM) %in% c("NA"), ]

# Define desired order for y-axis labels (rows)
desired_order <- c("Endothelial Cells",
                   "Sinusoidal Endothelial Cells of Liver",
                   "Hepatocytes",
                   "Collecting Duct Cells",
                   "Proximal Tubule Cells",
                   "Distal Convoluted Tubule Cells",
                   "Pulmonary Alveolar Type I",
                   "Clara Cells",
                   "Pulmonary alveolar type II",
                   "Cardiomyocytes",
                   "Fibroblasts",
                   "Mesothelial Cells",
                   "Astrocytes",
                   "Neurons",
                   "Oligodendrocyte Precursor Cells",
                   "Oligodendrocyte",
                   "Pancreatic Beta Cells",
                   "Pancreatic Acinar Cells",
                   "Enteroendocrine Cells",
                   "Paneth Cells",
                   "Enterocytes",
                   "Goblet Cells",
                   "T Cells",
                   "B Cells",
                   "Dendritic Cells",
                   "Microglial Cells",
                   "Alveolar Macrophages",
                   "Macrophages")

heatmapEM <- heatmapEM[match(desired_order, rownames(heatmapEM)), ]

col_fun = colorRamp2(c(100, 50, 0), c("#fcad03", "white", "#036bfc"))
heatmap2 <- Heatmap(heatmapEM, 
                   name = "Norm. Enrichment\n-10log(P−adj)", 
                   col = col_fun, 
                   row_names_gp = gpar(fontsize = 10, fontface = "bold"), 
                   show_row_dend = FALSE,
                   column_names_gp = gpar(fontsize = 0),  
                   row_title_side = "left", 
                   rect_gp = gpar(col = "white", lwd = 2),
                   heatmap_legend_param = list(
                     direction = "horizontal",
                     legend_width = unit(3, "cm"),  
                     legend_height = unit(0.8, "cm"),  
                     title_position = "topcenter",
                     title_gp = gpar(fontsize = 10),  
                     labels_gp = gpar(fontsize = 8),  
                     grid_height = unit(0.4, "cm"),  
                     grid_width = unit(0.4, "cm")    
                   ),
                   row_order = 1:nrow(heatmapEM))
draw(heatmap2, heatmap_legend_side = "bottom")



##figure C

Mouse <- addGroupCoverages(ArchRProj = Mouse, groupBy = "Clusters",force=TRUE)
motifs <- c("Runx","Klf")
motifPositions <- getPositions(Mouse)
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

seFoot <- getFootprints(
  ArchRProj = Mouse, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = Mouse,
  normMethod = "None",
  addDOC = FALSE,
  smoothWindow = 5,
  height = 8,
  width = 9
)

# figure D

library(readr)

df <- read.csv("/cluster/projects/epigenomics/Aminnn/scATAC/SEACells/Mouse_9Tissues/export/gene_scores.csv")
colnames(df) <- gsub("\\.([^.]+)\\.(\\d+)$", "#\\1-\\2", colnames(df))

pdata <- read_csv("/cluster/projects/epigenomics/Aminnn/scATAC/SEACells/Mouse_9Tissues/export/SEACell_metacells.csv")

celltypesc <- c(
  "Sinusoidal Endothelial Cells of Liver", "Hepatocytes", "Collecting Duct Cells",
  "Proximal Tubule Cells", "Distal Convoluted Tubule Cells", "Pulmonary Alveolar Type I",
  "Clara Cells", "Pulmonary alveolar Type II", "Cardiomyocytes", "Fibroblasts",
  "Mesothelial Cells", "Astrocytes", "Oligodendrocyte Precursor Cells", "Oligodendrocyte",
  "Pancreatic Beta Cells", "Pancreatic Acinar Cells", "Enteroendocrine Cells",
  "Paneth Cells", "Enterocytes", "Goblet Cells", "T Cells", "Dendritic Cells",
  "Microglial Cells", "Alveolar Macrophages", "Macrophages", "Endothelial Cells",
  "B Cells", "Neurons"
)

# Function to process each cell type
process_celltype <- function(celltype, df, pdata, threshold = 2) {
  celltype_data <- df[, pdata$x[pdata$Clusters == celltype]]
  
  celltype_means <- rowMeans(celltype_data)
  
  other_data <- df[, pdata$x[pdata$Clusters != celltype]]
  other_means <- rowMeans(other_data)
  
  high_in_celltype <- celltype_means >= threshold
  low_in_others <- other_means < threshold
  
  selected_genes <- df$X[high_in_celltype & low_in_others]
  selected_scores <- celltype_means[high_in_celltype & low_in_others]
  
  return(data.frame(gene = selected_genes, score = selected_scores))
}

# Process all cell types
results_list <- list()

for (celltype in celltypesc) {
  results_list[[celltype]] <- process_celltype(celltype, df, pdata)
  
  # Print progress
  cat("Processed:", celltype, "\n")
  cat("Number of genes found:", nrow(results_list[[celltype]]), "\n\n")
}

for (celltype in celltypesc) {
  write.csv(results_list[[celltype]], 
            file = paste0(gsub(" ", "_", celltype), "_specific_genes.csv"), 
            row.names = FALSE)
}









