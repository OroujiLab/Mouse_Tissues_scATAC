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


Mouse <- readRDS("./outputs/Archr_processed/Save-ArchR-Project.rds")
rna <- readRDS("./outputs/scRNA_combined.rds")

Mouse <- addGroupCoverages(ArchRProj = Mouse, groupBy = "Clusters",force = TRUE)

# figure A
Mouse <- addGeneIntegrationMatrix(
  ArchRProj = Mouse,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = rna,
  addToArrow = FALSE,
  groupRNA = "celltype",
  nameCell = "scRNA_int_cells",
  nameGroup = "scRNA_int_groups",
  nameScore = "scRNA_score",
  force = TRUE
)

pathToMacs2 <- findMacs2()

Mouse <- addReproduciblePeakSet(
  ArchRProj = Mouse,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2
)

Mouse <- addPeakMatrix(Mouse)

markersPeaks <- getMarkerFeatures(
  ArchRProj = Mouse, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks,"./outputs/markerPeaks.rds")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
write.csv(markerList,"./outputs/PeakMarkerList_dataframe.csv")


Mouse <- addPeak2GeneLinks(
  ArchRProj = Mouse,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix"
)

p1 <- plotPeak2GeneHeatmap(ArchRProj = Mouse, groupBy = "Clusters")


######## figure B
library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(cowplot)
set.seed(123)
options(scipen=1)


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakmatrix_df <- read.csv("/cluster/projects/epigenomics/Aminnn/scATAC/ArchR/markerList_Dec18.csv", header = TRUE)

peakmatrix_df2 <- peakmatrix_df[, c("seqnames", "start", "end", "group_name","FDR")]
peakmatrix_df3 <- makeGRangesFromDataFrame(peakmatrix_df2,keep.extra.columns=TRUE)
peakmatrix_df3 <- peakmatrix_df3[peakmatrix_df3$group_name != "NA"]

celltypesc <- c(
  "Endothelial Cells", "Sinusoidal Endothelial Cells of Liver", "Hepatocytes",
  "Collecting Duct Cells", "Proximal Tubule Cells", "Distal Convoluted Tubule Cells",
  "Pulmonary Alveolar Type I", "Clara Cells", "Pulmonary alveolar type II",
  "Cardiomyocytes", "Fibroblasts", "Mesothelial Cells", "Astrocytes", "Neurons",
  "Oligodendrocyte Precursor Cells", "Oligodendrocyte", "Pancreatic Beta Cells",
  "Pancreatic Acinar Cells", "Enteroendocrine Cells", "Paneth Cells", "Enterocytes",
  "Goblet Cells", "T Cells", "B Cells", "Dendritic Cells", "Microglial Cells",
  "Alveolar Macrophages", "Macrophages"
)

generate_plot <- function(celltype) {
  tryCatch({
    Celltype2 <- peakmatrix_df3[peakmatrix_df3$group_name == celltype]
    if(length(Celltype2) == 0) {
      stop("No data available for this cell type")
    }
    Celltype2_2 = annotatePeak(Celltype2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
    Celltype2_3 = as.data.frame(Celltype2_2)
    pathway.GO2_Celltype2 <- enrichGO(as.data.frame(Celltype2_3)$geneId, org.Mm.eg.db, ont = "BP")
    
    if(nrow(pathway.GO2_Celltype2) == 0) {
      stop("No enriched GO terms found")
    }
    
    p <- barplot(pathway.GO2_Celltype2, showCategory=10) +
      ggtitle(celltype) +
      scale_fill_continuous(
        low = "#132B43", 
        high = "#56B1F7",
        name = "Adjusted p-value"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 5),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6)
      ) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 4, barheight = 0.5))
    
    return(p)
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", celltype, e$message))
    # Create a blank plot with an error message
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = paste("No data available for", celltype)) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 5)) +
      ggtitle(celltype)
    return(p)
  })
}

all_plots <- lapply(celltypesc, generate_plot)

save_plots_to_pdf <- function(plot_list, filename, plots_per_page = 9, ncol = 3) {
  total_plots <- length(plot_list)
  total_pages <- ceiling(total_plots / plots_per_page)
  
  pdf(filename, width = 11, height = 8.5)  # Standard letter size
  
  for(i in 1:total_pages) {
    start_index <- (i-1) * plots_per_page + 1
    end_index <- min(i * plots_per_page, total_plots)
    page_plots <- plot_list[start_index:end_index]
    
    while(length(page_plots) < plots_per_page) {
      page_plots <- c(page_plots, list(ggplot() + theme_void()))
    }
    
    print(plot_grid(plotlist = page_plots, ncol = ncol))
  }
  
  dev.off()
}

save_plots_to_pdf(all_plots, "celltype_GO_plots.pdf", plots_per_page = 9, ncol = 3)



# figure C

getGroupBW(
  ArchRProj = Mouse,
  groupBy = "Sample",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)


peaks_score <- readRDS("~/Downloads/Transfers/markerPeaks_Clusters50_NOV27.rds")
markerList <- getMarkers(peaks_score, cutOff = "FDR <= 0.01 & Log2FC >= 1")

tcells <- markerList@listData$`T Cells`
tcells <- head(tcells, n=10)
tcells <- makeGRangesFromDataFrame(tcells)

bw_file_tcell <- "~/Downloads/bigwig Files/Cell Types/T.Cells-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"
bw_file_tcell <- import(bw_file_tcell)

bed_tcell <- "~/Downloads/BedFiles_MarkerPeaks/T Cells.bed"
bed_tcell1 <- import(bed_tcell)

neurons <- markerList@listData$Neurons
neurons <- makeGRangesFromDataFrame(neurons)

bw_file_neuron <- "~/Downloads/bigwig Files/Cell Types/Neurons-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"
bw_file_neuron <- import(bw_file_neuron)

bed_neuron <- "~/Downloads/BedFiles_MarkerPeaks/Neurons.bed"
bed_neuron <- import(bed_neuron)

# Endothelial Cells data
endo <- markerList@listData$`Endothelial Cells`
endo <- makeGRangesFromDataFrame(endo)

bw_file_endo <- "~/Downloads/bigwig Files/Cell Types/Endothelial.Cells-TileSize-100-normMethod-ReadsInTSS.bw"
bw_file_endo <- import(bw_file_endo)

bed_endo <- "~/Downloads/BedFiles_MarkerPeaks/Endothelial Cells.bed"
bed_endo <- import(bed_endo)

create_gene_region_track <- function(chr, start, end, shape="arrow", stack_height=0.3) {
  grtrack <- GeneRegionTrack(
    TxDb.Mmusculus.UCSC.mm10.knownGene,
    chromosome = chr, start = start, end = end,
    showId = TRUE, transcriptAnnotation = "symbol", collapseTranscripts = TRUE,
    name = "",
    stackHeight = stack_height, fontsize.group = 8,
    shape = shape,
    stacking = "squish"
  )
  
  grtrack_range <- grtrack@range
  range_mapping <- select(
    Mus.musculus,
    keys = mcols(grtrack_range)$symbol,
    keytype = "TXNAME",
    columns = c("ENTREZID", "SYMBOL")
  )
  symbol(grtrack) <- range_mapping$SYMBOL
  
  return(grtrack)
}

plot_genomic_tracks <- function(chr, start, end, bed_start, bed_end, bw_start, bw_end, 
                                bw_file, bed_data, cell_type, color, shape="arrow", stack_height=0.3) {
  genome_track <- GenomeAxisTrack()
  
  bigwig_track <- DataTrack(bw_file, type = "h", name = paste(cell_type, ""), col = "black",
                            start = bw_start, end = bw_end, chromosome = chr,
                            background.panel = "white", background.title = "white")
  
  bed_track <- AnnotationTrack(bed_data, name = "",
                               start = bed_start, end = bed_end, chromosome = chr,
                               col = color, background.panel = "white")
  
  grtrack <- create_gene_region_track(chr, start, end, shape = shape, stack_height = stack_height)
  
  plotTracks(list(genome_track, bigwig_track, bed_track, grtrack),
             background.title = "gray", type = "h", from = start, to = end, chromosome = chr)
}

plot_genomic_tracks("chr6", 41000000, 42000000, 41554500, 41554999, 41554500, 41554999,
                    bw_file_tcell, tcells, "T Cells", "darkblue")

plot_genomic_tracks("chr8", 94000000, 95000000, 94413500, 94393999, 94413500, 94413999,
                    bw_file_tcell, tcells, "T Cells", "darkblue")

plot_genomic_tracks("chr19", 4900000, 5100000, 5001697, 5002197, 5001697, 5002197,
                    bw_file_neuron, neurons, "Neurons", "#9C82BA", shape="box")

plot_genomic_tracks("chr5", 146800000, 147800000, 147560196, 147727988, 147000000, 148000000,
                    bw_file_endo, endo, "EC", "#CCC759", shape="box", stack_height=0.1)



## figure D
library(vegan)
library(proxy)
library(pvclust)
set.seed(1)

file_path <- "./outputs/cM.csv"  
data <- read.csv(file_path, header = TRUE, row.names = 1)
data <- data[rownames(data) != "Unknown", ]
data <- data[, colnames(data) != "Total"]

cosine_similarity <- proxy::simil(x = as.matrix(data), method = "cosine")
cosine_matrix <- as.matrix(cosine_similarity)

heatmap(cosine_matrix, 
        main = "Cosine Similarity",
        col = colorRampPalette(c("blue", "white", "red"))(50),
        margins = c(13, 13))

transposed_data <- t(data)
pvclust_result2 <- pvclust(transposed_data, 
                           method.hclust = "average", 
                           method.dist = "correlation", 
                           nboot = 10000)

plot(pvclust_result2, 
     print.num = FALSE,     
     print.pv = FALSE,      
     print.bp = FALSE,      
     hang = -1,             
     cex = 0.8,             
     main = "",             
     sub = "",              
     xlab = "",             
     ylab = "")

cluster_results <- pvclust_result2$edges
format_columns <- c("si", "au", "bp", "se.si", "se.au", "se.bp", "v", "c", "pchi")

for (col in format_columns) {
  cluster_results[[col]] <- format(cluster_results[[col]], digits = 6)
}

print(cluster_results)
