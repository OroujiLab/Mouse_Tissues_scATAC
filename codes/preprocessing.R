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


setwd("~/Mouse_Tissues_scATAC")


inputFiles <- c(
  'Heart' = "./inputs/fragments_heart.tsv.gz", 
  'Kidney' = "./inputs/fragments_kidney.tsv.gz",
  'Brain' = "./inputs/fragments_brain.tsv.gz", 
  'Liver' = "./inputs/fragments_liver.tsv.gz", 
  'Lung' = "./inputs/fragments_lung.tsv.gz", 
  'Spleen' = "./inputs/fragments_spleen.tsv.gz",
  'Colon' = "./inputs/fragments_colon.tsv.gz",
  'Small_Intestine' = "./inputs/fragments_small_intestine.tsv.gz",
  'Pancreas' = "./inputs/fragments_pancreas.tsv.gz"
)
print(inputFiles)


addArchRGenome("mm10")
addArchRThreads(threads = 1)
addArchRLogging(useLogs = TRUE)


geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Mmusculus.UCSC.mm10.ensGene, OrgDb = org.Mm.eg.db)
geneAnnotation <- createGeneAnnotation(
  TSS = geneAnnotation$TSS,
  exons = geneAnnotation$exons,
  genes = geneAnnotation$genes
)


ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
print(ArrowFiles)


ArrowFiles <- c("Kidney.arrow", "Liver.arrow", "Lung.arrow", "Heart.arrow", "Brain.arrow", "Spleen.arrow", "Small_Intestine.arrow", "Colon.arrow", "Pancreas.arrow")


doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)
print(doubScores)


Mouse <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "./outputs/Archr_processed",
  copyArrows = TRUE 
)
print(Mouse)
