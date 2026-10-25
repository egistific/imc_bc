log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

library(yaml)
library(SpatialExperiment)
library(dplyr)
library(bluster)
library(BiocParallel)
library(ggplot2)
library(scran)
library(FlowSOM)
library(kohonen)
library(ConsensusClusterPlus)
library(data.table)

config <- snakemake@config
cofactor <- config$analysis$arcsinh_cofactor
seed <- config$project$seed
output_file <- snakemake@output[["spe"]]
exprs_file <- snakemake@output[["exprs"]]
design_file <- snakemake@output[["design"]]

cofactor <- config$analysis$arcsinh_cofactor
wq <- config$analysis$astir$winsor_quantile
batch_col <- config$analysis$astir$batch_col
spe <- readRDS(snakemake@input[["spe"]])

# Keep primary tumors and lymph node metastases
keep_tissue_types <- config$analysis$keep_tissue_types
stopifnot(all(keep_tissue_types %in% spe$TissueType))
spe <- spe[, spe$TissueType %in% keep_tissue_types]

# Keep images with at least 100 cells to exclude cores that may have missed the tumors
df_count <- as.data.frame(colData(spe)) %>% dplyr::group_by(sample_id) %>% dplyr::summarise(num_cells = dplyr::n()) %>% dplyr::arrange(num_cells)
keep_sample_ids <- df_count$sample_id[df_count$num_cells >= 100]
spe[, spe$sample_id %in% keep_sample_ids]
spe <- spe[, spe$sample_id %in% keep_sample_ids]

spe$short_px <- pmin(spe$width_px, spe$height_px)
spe$narrow_image <- spe$short_px < 200
spe$area_small <- spe$area < 10
spe$area_large <- spe$area > 500
spe$dropped <- spe$narrow_image | spe$area_small | spe$area_large

df_image <- as.data.frame(colData(spe)) %>%
  dplyr::group_by(sample_id, TissueType) %>%
  dplyr::summarise(num_cells = dplyr::n(), 
                   num_small = sum(area_small), 
                   percent_small = mean(area_small) * 100, 
                   num_large = sum(area_large), 
                   percent_large = mean(area_large) * 100, 
                   num_dropped = sum(dropped), 
                   percent_dropped = mean(dropped) * 100, 
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(percent_dropped))

# Drop cells
spe <- spe[, !spe$dropped]
dim(spe)

# Use spillover compensated values
assay(spe, "exprs") <- asinh(assay(spe, "compcounts") / cofactor)
# Filter low-quality cells
dna1_marker = rownames(spe)[rowData(spe)$marker == "DNA1"]
# Cut-off is at valley of bimodal distribution
dna1_cutoff <- 1.2
dna1_values <- assay(spe, "exprs")[dna1_marker, ]
keep <- dna1_values >= dna1_cutoff
spe <- spe[, keep]

saveRDS(spe, file = output_file, compress = FALSE)
message("Saved filtered single-cell object to: ", output_file)

astir_names <- c(Lu175 = "panCK",    Nd144 = "CK8_18",   Yb174 = "CK7",      Pr141 = "CK5",
                 Sm147 = "CK14",     Er167 = "Ecad",     Dy162 = "CD45",     Sm152 = "CD3",
                 Yb171 = "CD4",      Ho165 = "CD8a",     Dy164 = "CD20",     Nd146 = "CD68",
                 Sm154 = "CD11c",    Nd143 = "HLADR",    Nd145 = "CD15",     Yb172 = "vWF",
                 Nd148 = "SMA",      Sm149 = "Vimentin", Nd142 = "Fibronectin")

rows <- match(names(astir_names), rowData(spe)$metal)
ex <- as.matrix(assay(spe, "exprs")[rows, ]) # markers x cells

# Winsorize markers
ex <- apply(ex, 1, function(v) pmin(v, quantile(v, wq))) # cells x markers
colnames(ex) <- unname(astir_names) # drop metal tags
data.table::fwrite(data.frame(cell_id = colnames(spe), ex, check.names = FALSE), exprs_file) # use fwrite for speed

# --- design: intercept + one column per TMA except the first ---
batch <- factor(colData(spe)[[batch_col]])
table(batch)                         # check the TMA counts look right

design <- model.matrix(~ batch)
colnames(design) <- make.names(colnames(design))   # e.g. X.Intercept., batchZTMA22 ...
data.table::fwrite(data.frame(cell_id = colnames(spe), design, check.names = FALSE), design_file)
