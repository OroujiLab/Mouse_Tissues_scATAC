# Load required libraries
library(ArchR)
library(Seurat)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(org.Mm.eg.db)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(parallel)
library(pheatmap)
library(scran)
library(Signac)
library(sctransform)
library(SeuratObject)
library(dplyr)
library(Cairo)
library(Matrix)

# Set working directory and load data
setwd("/cluster/projects/epigenomics/Aminnn/scATAC/Results")
seRNA <- readRDS("/cluster/projects/epigenomics/Aminnn/scATAC/new/seRNA_Seurat.rds")

# Check data
head(colnames(seRNA))

# Define input files
inputFiles <- c(
  'Heart' = "fragments_heart.tsv.gz", 
  'Kidney' = "fragments_kidney.tsv.gz",
  'Brain' = "fragments_brain.tsv.gz", 
  'Liver' = "fragments_liver.tsv.gz", 
  'Lung' = "fragments_lung.tsv.gz", 
  'Spleen' = "fragments_spleen.tsv.gz",
  'Colon' = "fragments_colon.tsv.gz",
  'Small_Intestine' = "fragments_small_intestine.tsv.gz",
  'Pancreas' = "fragments_pancreas.tsv.gz"
)
print(inputFiles)

# Set up ArchR environment
addArchRGenome("mm10")
addArchRThreads(threads = 1)
addArchRLogging(useLogs = TRUE)

# Create gene annotation
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Mmusculus.UCSC.mm10.ensGene, OrgDb = org.Mm.eg.db)
geneAnnotation <- createGeneAnnotation(
  TSS = geneAnnotation$TSS,
  exons = geneAnnotation$exons,
  genes = geneAnnotation$genes
)

# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
print(ArrowFiles)

# Define Arrow files
ArrowFiles <- c("Kidney.arrow", "Liver.arrow", "Lung.arrow", "Heart.arrow", "Brain.arrow", "Spleen.arrow", "Small_Intestine.arrow", "Colon.arrow", "Pancreas.arrow")

# Add doublet scores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)
print(doubScores)

# Create ArchR project
Mouse <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "/cluster/projects/epigenomics/Aminnn/scATAC/ArchR/Final/Results_Final",
  copyArrows = TRUE 
)
print(Mouse)