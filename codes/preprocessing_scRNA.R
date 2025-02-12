#!/bin/bash


#mkdir -p data outputs logs



datasets=(
    "TMFACS_Colon.rds"
    "TMFACS_Kidney.rds"
    "TMFACS_Liver.rds"
    "TMFACS_Heart.rds"
    "TMFACS_BrainMicroglia.rds"
    "TMFACS_BrainNeurons.rds"
    "TMFACS_Spleen.rds"
    "TMFACS_Lung.rds"
    "TMFACS_Pancreas.rds"
    "MCA_SmallIntestine.rds"
)

for dataset in "${datasets[@]}"; do
    if [ ! -f "data/$dataset" ]; then
        echo "Downloading $dataset"
        wget -P data "https://scfind.cog.sanger.ac.uk/reddim_sce/TabulaMurisFACS/$dataset" 2>>logs/download_errors.log
    else
        echo "$dataset already exists."
    fi
done

if [ ! -f "data/MCA_SmallIntestine.rds" ]; then
    wget -P data "https://scfind.cog.sanger.ac.uk/reddim_sce/MouseCellAtlas/MCA_SmallIntestine.rds" 2>>logs/download_errors.log
fi


Rscript - <<'EOF'


suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(purrr)
})


file_paths <- c(
    heart = "data/TMFACS_Heart.rds",
    kidney = "data/TMFACS_Kidney.rds",
    spleen = "data/TMFACS_Spleen.rds",
    liver = "data/TMFACS_Liver.rds",
    lung = "data/TMFACS_Lung.rds",
    pancreas = "data/TMFACS_Pancreas.rds",
    brain_microglia = "data/TMFACS_BrainMicroglia.rds",
    brain_neurons = "data/TMFACS_BrainNeurons.rds",
    colon = "data/TMFACS_Colon.rds",
    small_intestine = "data/MCA_SmallIntestine.rds"
)


read_and_preprocess <- function(file_path) {
    tryCatch({
        message(sprintf("Processing %s", file_path))
        sce <- readRDS(file_path)
        
        
        cols_to_remove <- c("mouse", "well", "subtissue", "FACS", "sex", 
                           "ClusterID", "Batch", "Mouse.Sex.Age", "file", 
                           "cell_type1", "Cell.Barcode")
        
        
        existing_cols <- intersect(colnames(colData(sce)), cols_to_remove)
        if (length(existing_cols) > 0) {
            sce <- sce[, !colnames(colData(sce)) %in% existing_cols]
        }
        
        
        colnames(colData(sce)) <- c("Tissue", "Cell_type")
        
        return(sce)
    }, error = function(e) {
        message(sprintf("Error processing %s: %s", file_path, e$message))
        return(NULL)
    })
}

sce_list <- map(file_paths, read_and_preprocess)
sce_list <- Filter(Negate(is.null), sce_list)

message("Combining SingleCellExperiment objects...")
combined_sce <- do.call(cbind, sce_list)


message("Converting to Seurat object")
scRNA <- as.Seurat(
    combined_sce,
    counts = "counts",
    data = "logcounts",
    assay = NULL,
    project = "SingleCellExperiment"
)

# Rename Assay
scRNA <- RenameAssays(scRNA, originalexp = 'RNA')

# Remove rows with NA values
message("Removing NA values...")
na_rows <- rowSums(is.na(scRNA[["RNA"]]@counts)) > 0
scRNA <- subset(scRNA, features = rownames(scRNA)[!na_rows])
saveRDS(scRNA, "outputs/scRNA_combined.rds")


EOF

echo "Script execution completed!"
