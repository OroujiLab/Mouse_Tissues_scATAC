library(ArchR)
library(Seurat)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(org.Mm.eg.db)
library(parallel)
library(pheatmap)
library(scran)
library(Signac)
library(sctransform)
library(SeuratObject)
library(dplyr)
library(Cairo)
library(Matrix)
library(scDataviz)
set.seed(123)


Mouse <- readRDS("./outputs/Archr_processed/Save-ArchR-Project.rds")

# Add IterativeLSI
Mouse <- addIterativeLSI(
  ArchRProj = Mouse,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 4,
  clusterParams = list(
    resolution = c(0.1, 0.2, 0.4),
    sampleCells = 60000,
    n.start = 10
  ),
  varFeatures = 65000,
  dimsToUse = 1:30
)
print(Mouse)


# Add clusters
Mouse <- addClusters(
  input = Mouse,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  maxClusters = 60,
  resolution = 1,
  force = TRUE
)
print(Mouse)

#figure C and E
Mouse <- addUMAP(
  ArchRProj = Mouse, 
  reducedDims = "IterativeLSI", 
  name = "UMAP_LowRes", 
  nNeighbors = 40, 
  minDist = 0.5, 
  metric = "cosine"
)

# Define cluster annotations (manually)
word_replacements <- c(
  "C1" = "T Cells",
  "C2" = "B Cells",
  "C3" = "Dendritic Cells",
  "C4" = "B Cells",
  "C5" = "B Cells",
  "C6" = "Alveolar Macrophages",
  "C7" = "Macrophages",
  "C8" = "Macrophages",
  "C9" = "Microglial Cells",
  "C10" = "Oligodendrocyte",
  "C11" = "NA",
  "C12" = "Astrocytes",
  "C13" = "Fibroblasts",
  "C14" = "Oligodendrocyte Precursor Cells",
  "C15" = "Enterocytes",
  "C16" = "Enterocytes",
  "C17" = "Enterocytes",
  "C18" = "Enterocytes",
  "C19" = "Paneth Cells",
  "C20" = "Goblet Cells",
  "C21" = "Goblet Cells",
  "C22" = "Enterocytes",
  "C23" = "Enterocytes",
  "C24" = "Endothelial Cells",
  "C25" = "Endothelial Cells",
  "C26" = "Endothelial Cells",
  "C27" = "Sinusoidal Endothelial Cells of Liver",
  "C28" = "Cardiomyocytes",
  "C29" = "Cardiomyocytes",
  "C30" = "Mesothelial Cells",
  "C31" = "Fibroblasts",
  "C32" = "Fibroblasts",
  "C33" = "Fibroblasts",
  "C34" = "Fibroblasts",
  "C35" = "Hepatocytes",
  "C36" = "Proximal Tubule Cells",
  "C37" = "Proximal Tubule Cells",
  "C38" = "Pancreatic Beta Cells",
  "C39" = "Astrocytes",
  "C40" = "Neurons",
  "C41" = "Oligodendrocyte Precursor Cells",
  "C42" = "Pancreatic Acinar Cells",
  "C43" = "Pancreatic Acinar Cells",
  "C44" = "Pulmonary Alveolar Type I",
  "C45" = "Clara Cells",
  "C46" = "Pulmonary Alveolar Type II",
  "C47" = "Collecting Duct Cells",
  "C48" = "Distal Convoluted Tubule Cells",
  "C49" = "Endothelial Cells",
  "C50" = "Enteroendocrine Cells"
)

# Check unique clusters
print(unique(Mouse$Clusters50))

# Apply annotations
for (key in names(word_replacements)) {
  replacement <- word_replacements[key]
  pattern <- paste0("\\b", key, "\\b")
  Mouse$Clusters50 <- gsub(pattern, replacement, Mouse$Clusters50)
}
print(head(Mouse$Clusters50))
df <- confusionMatrix(Mouse$Clusters, Mouse$Sample)
write.csv(df,"./outputs/cM.csv")
p1 <- plotEmbedding(ArchRProj = Mouse, colorBy = "cellColData", name = "Sample", embedding = "UMAP_LowRes")
p1
p2 <- plotEmbedding(ArchRProj = Mouse, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_LowRes")
p2

#figre D

umap <- getReducedDims(ArchRProj = Mouse, reducedDims = "IterativeLSI", returnMatrix = TRUE)
write.csv(umap,"./outputs/umap_values.csv")

contourPlot(umap,
            reducedDim = 'UMAP',
            bins = 50,
            subtitle = 'UMAP performed on All Tissue',
            legendLabSize = 18,
            axisLabSize = 22,
            titleLabSize = 22,
            subtitleLabSize = 18,
            captionLabSize = 18,
            gridlines.major = FALSE,
            gridlines.minor = FALSE,
            lowcol = "#036bfc",
            highcol = "#fcad03")



#figure F
# Load and preprocess data
df <- read.csv("./outputs/cM.csv", header = TRUE, row.names = 1)
df$Total <- NULL
samples <- c("Brain","Colon","Heart","Kidney","Liver","Lung","Pancreas","Small_Intestine","Spleen")
colours <- c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B", "#FEE500","#8A9FD1","#C06CAB","#F9B712")
df <- df[rownames(df) != "Unknown", ]

# Create hierarchical clustering
data_log <- log1p(df)
dist_mat <- dist(data_log)
hc <- hclust(dist_mat, method = "ward.D2")

ggtree_plot <- ggtree::ggtree(hc) + geom_tiplab(align=TRUE, linesize = 1, offset = 1) + hexpand(2)

# Prepare data for ggplot
df$Cell_Types <-rownames(df)
df_melt <- melt(df, varnames = c("Cell_Types", "Tissue_Types"))
colnames(df_melt) <- c("Cell_Types","Tissue_Types","value")
df_percent <- df_melt %>%
  group_by(Cell_Types) %>%
  mutate(Percentage = value / sum(value)) %>%
  ungroup() %>%
  dplyr::filter(Percentage > 0.05) %>%
  mutate(Tissue_Types = factor(Tissue_Types, levels = samples))

# Define desired order of cell types
matched_order <- c(
  "Mesothelial Cells", "Cardiomyocytes", "Pancreatic Beta Cells",
  "Pancreatic Acinar Cells", "Hepatocytes", "Sinusoidal Endothelial Cells of Liver",
  "Distal Convoluted Tubule Cells", "Proximal Tubule Cells", "Collecting Duct Cells",
  "Pulmonary Alveolar Type I", "Alveolar Macrophages", "Pulmonary alveolar type II",
  "Clara Cells", "Microglial Cells", "Astrocytes", "Oligodendrocyte", "Neurons",
  "Oligodendrocyte Precursor Cells", "Paneth Cells", "Enteroendocrine Cells",
  "Dendritic Cells", "Goblet Cells", "Enterocytes", "B Cells", "T Cells", "Macrophages",
  "Fibroblasts", "Endothelial Cells"
)

# Create stacked bar plot
mx <- ggplot(df_percent, aes(x = Cell_Types, y = Percentage, fill = Tissue_Types)) + 
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = colours, limits = samples, name = "Tissue") +
  labs(x = "Cell Types", y = "") +
  coord_flip() +
  scale_x_discrete(limits = rev(matched_order)) +
  theme_minimal() +
  geom_text(aes(label = paste0(sprintf("%.1f", Percentage*100), "%")), 
            position = position_fill(vjust = 0.5), color = "black", size = 2) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = 'none',
    legend.box.background = element_blank()
  )

# Combine plots
combined_plot <- plot_grid(ggtree_plot, mx, nrow = 1, rel_widths = c(0.1, 0.2), align = 'h')

# Get cell counts
cell_counts <- df_melt %>%
  group_by(Cell_Types) %>%
  summarize(Total_Cells = sum(value)) %>%
  arrange(match(Cell_Types, matched_order))

combined_plot_with_labels <- ggdraw() +
  draw_plot(combined_plot) +
  draw_text(
    text = cell_counts$Total_Cells,
    x = 0.99, y = seq(0.98, 0.07, length.out = nrow(cell_counts)),
    hjust = 0.7, vjust = 1.9,  # Adjust hjust to align to the right
    size = 11, color = "black"
  )
combined_plot_with_labels






