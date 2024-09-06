library(ArchR)
library(Seurat)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(org.Mm.eg.db)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome.Mmusculus.UCSC.mm10)
library(parallel)
library(pheatmap)
library(scran)
library(Signac)
library(sctransform)
library(SeuratObject)
library(dplyr)
library(Cairo)
library(circlize)


Mouse <- readRDS("./outputs/Archr_processed/Save-ArchR-Project.rds")


dir.create(sprintf("./outputs/", proj_name))
write.csv(getReducedDims(Mouse), sprintf('./outputs/svd.csv', proj_name), quote=FALSE)
write.csv(getCellColData(Mouse), sprintf('./outputs/cell_metadata.csv', proj_name), quote=FALSE)


# Gene scores
gene.scores <- getMatrixFromProject(Mouse)
scores <- assays(gene.scores)[['GeneScoreMatrix']]
scores <- as.matrix(scores)
rownames(scores) <- rowData(gene.scores)$name
write.csv(scores, './outputs/gene_scores.csv', quote=FALSE)



# Peak counts
peaks <- getPeakSet(Mouse)
peak.counts <- getMatrixFromProject(proj, 'PeakMatrix')

# Reorder peaks
# Chromosome order
chr_order <- sort(seqlevels(peaks))
reordered_features <- list()
for(chr in chr_order)
  reordered_features[[chr]] = peaks[seqnames(peaks) == chr]
reordered_features <- Reduce("c", reordered_features)

# Export counts
dir.create(sprintf("./outputs/peak_counts", proj_name))
counts <- assays(peak.counts)[['PeakMatrix']]
writeMM(counts, sprintf('./outputs/counts.mtx', proj_name))
write.csv(colnames(peak.counts), './outputs/cells.csv', quote=FALSE)
names(reordered_features) <- sprintf("Peak%d", 1:length(reordered_features))
write.csv(as.data.frame(reordered_features),'./outputs/peaks.csv', quote=FALSE)



