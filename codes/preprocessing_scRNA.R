# Download datasets
# Please copy the codes below to your terminal to download each of the datasets:
# Colon
# wget 'https://scfind.cog.sanger.ac.uk/reddim_sce/TabulaMurisFACS/TMFACS_Colon.rds'
# Kidney
# wget 'https://scfind.cog.sanger.ac.uk/reddim_sce/TabulaMurisFACS/TMFACS_Kidney.rds'
# Liver
# wget 'https://scfind.cog.sanger.ac.uk/reddim_sce/TabulaMurisFACS/TMFACS_Liver.rds'
# Heart
# wget 'https://scfind.cog.sanger.ac.uk/reddim_sce/TabulaMurisFACS/TMFACS_Heart.rds'
# Brain
# wget 'https://scfind.cog.sanger.ac.uk/reddim_sce/TabulaMurisFACS/TMFACS_BrainMicroglia.rds'
# wget 'https://scfind.cog.sanger.ac.uk/reddim_sce/TabulaMurisFACS/TMFACS_BrainNeurons.rds'
# Spleen
# wget 'https://scfind.cog.sanger.ac.uk/reddim_sce/TabulaMurisFACS/TMFACS_Spleen.rds'
# Lung
# wget 'https://scfind.cog.sanger.ac.uk/reddim_sce/TabulaMurisFACS/TMFACS_Lung.rds'
# Pancreas
# wget 'https://scfind.cog.sanger.ac.uk/reddim_sce/TabulaMurisFACS/TMFACS_Pancreas.rds'
# Small_Intestine
# wget 'https://scfind.cog.sanger.ac.uk/reddim_sce/MouseCellAtlas/MCA_SmallIntestine.rds'

# Defining each scRNA dataset
seRNA_heart <- readRDS("./outputs/TMFACS_Heart.rds")
seRNA_kidney <- readRDS("./outputs/TMFACS_Kidney.rds")
seRNA_spleen <- readRDS("./outputs/TMFACS_Spleen.rds")
seRNA_liver <- readRDS("./outputs/TMFACS_Liver.rds")
seRNA_lung <- readRDS("./outputs/TMFACS_Lung.rds")
seRNA_pancreas <- readRDS("./outputs/TMFACS_Spleen.rds")
seRNA_brain1 <- readRDS("./outputs/TMFACS_TMFACS_BrainMicroglia.rds")
seRNA_brain2 <- readRDS("./outputs/TMFACS_BrainNeurons.rds")
seRNA_colon <- readRDS("./outputs/TMFACS_Colon.rds")
seRNA_small_intestine <- readRDS("./outputs/MCA_SmallIntestine.rds'")

# Combine TMFACS files
TMFACS_Combined <- cbind(seRNA_heart, 
                         seRNA_kidney, 
                         seRNA_spleen, 
                         seRNA_liver, 
                         seRNA_lung, 
                         seRNA_pancreas, 
                         seRNA_brain1, 
                         seRNA_brain2, 
                         seRNA_colon, 
                         seRNA_small_intestine)
                         
# Load required libraries
library(Seurat)
library(SingleCellExperiment)
library(purrr)

# Define file paths
file_paths <- c(
 heart = "~/TMFACS_Heart.rds",
 kidney = "~/TMFACS_Kidney.rds",
 spleen = "~/TMFACS_Spleen.rds",
 liver = "~/TMFACS_Liver.rds",
 lung = "~/TMFACS_Lung.rds",
 pancreas = "~/TMFACS_Pancreas.rds",
 brain_microglia = "~/TMFACS_BrainMicroglia.rds",
 brain_neurons = "~/TMFACS_BrainNeurons.rds",
 colon = "~/TMFACS_Colon.rds",
 small_intestine = "~/MCA_SmallIntestine.rds"
)

# Function to read and preprocess RDS files
read_and_preprocess <- function(file_path) {
 sce <- readRDS(file_path)
 
 cols_to_remove <- c("mouse", "well", "subtissue", "FACS", "sex", 
                     "ClusterID", "Batch", "Mouse.Sex.Age", "file", "cell_type1", "Cell.Barcode")
 sce <- sce[, !colnames(colData(sce)) %in% cols_to_remove]
 
 # Rename columns
 colnames(colData(sce)) <- c("Tissue", "Cell_type")
 
 return(sce)
}

# Read and preprocess all files
sce_list <- map(file_paths, read_and_preprocess)

# Combine all SingleCellExperiment objects
combined_sce <- do.call(cbind, sce_list)

# Convert to Seurat object
scRNA <- as.Seurat(
 combined_sce,
 counts = "counts",
 data = "logcounts",
 assay = NULL,
 project = "SingleCellExperiment"
)

# Rename assay
scRNA <- RenameAssays(scRNA, originalexp = 'RNA')

# Remove rows with NA values
na_rows <- rowSums(is.na(scRNA[["RNA"]]@counts)) > 0
scRNA <- subset(scRNA, features = rownames(scRNA)[!na_rows])
saveRDS(scRNA,"./outputs/scRNA_combined.rds")