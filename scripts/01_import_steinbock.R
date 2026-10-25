log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

library(yaml)
library(imcRtools)
library(SpatialExperiment)
library(SingleCellExperiment)


config <- snakemake@config
root_dir <- config$project$root_dir
cofactor <- config$analysis$arcsinh_cofactor
seed <- config$project$seed
panel_file <- snakemake@input[["panel"]]
patient_file <- snakemake@input[["patient"]]
img_meta_file <- snakemake@input[["img_info"]]
output_file <- snakemake@output[["spe"]]

intensities_dir <- snakemake@input[["intensities"]]
regionprops_dir <- snakemake@input[["regionprops"]]
neighbors_dir <- snakemake@input[["neighbors"]]
extract_cellid_from <- config$analysis$extract_cellid_from
extract_coords_from <- config$analysis$extract_coords_from
extract_imagemetadata_from <- config$analysis$extract_imagemetadata_from
extract_names_from <- config$analysis$extract_names_from
return_as <- config$analysis$return_as

dir.create(config$analysis$output_dir, recursive=TRUE, showWarnings=FALSE)

spe <- imcRtools::read_steinbock(
  path = root_dir,
  intensities_folder = intensities_dir,
  regionprops_folder = regionprops_dir,
  graphs_folder = neighbors_dir,
  extract_cellid_from = extract_cellid_from,
  extract_coords_from = extract_coords_from,
  image_file = img_meta_file,
  extract_imagemetadata_from = extract_imagemetadata_from,
  panel_file = panel_file,
  extract_names_from = extract_names_from,
  return_as = return_as,
  BPPARAM = BiocParallel::SerialParam()
)

short <- c(In113 = "Histone H3", La139 = "H3 (Lys28) trimethylate", Pr141 = "Cytokeratin 5", 
           Nd142 = "Fibronectin", Nd143 = "HLA-DR", Nd144 = "Cytokeratin 8/18", Nd145 = "CD15", Nd146 = "CD68", 
           Sm147 = "Keratin 14 (KRT14)", Nd148 = "SMA", Sm149 = "Vimentin", Nd150 = "c-Myc", Eu151 = "c-erbB-2-Her", 
           Sm152 = "CD3", Eu153 = "p-Histone H3", Sm154 = "CD11c", Gd155 = "Bcl-2", Gd156 = "Estrogen Receptor Alpha", 
           Gd158 = "Progesterone Receptor A/B", Tb159 = "p53", Gd160 = "CD44", Dy161 = "Androgen Receptor AR", Dy162 = "CD45", 
           Dy163 = "GATA3", Dy164 = "CD20", Ho165 = "CD8a", Er166 = "Carbonic Anhydrase IX", Er167 = "E/P-Cadherin", 
           Er168 = "Ki-67", Tm169 = "EGFR", Er170 = "p-S6", Yb171 = "CD4", Yb172 = "vWF", Yb173 = "p-mTOR",
           Yb174 = "Cytokeratin 7", Lu175 = "pan Cytokeratin", Yb176 = "Cleaved Caspase3/PARP", Ir191 = "DNA1", Ir193 = "DNA2")
# Replace any string that ends with \\)\\) with blank
rowData(spe)$metal <- sub(".*\\)\\)", "", rownames(spe))
rowData(spe)$marker <- unname(short[rowData(spe)$metal])
stopifnot(!anyNA(rowData(spe)$marker))

img_meta <- read.csv(img_meta_file)
patient <- read.csv(patient_file, fileEncoding="latin1")

# Add batch data
tok <- strsplit(img_meta$acquisition_description, "_")
img_meta$ArrayNr <- vapply(tok, function(x) x[1], character(1))
img_meta$UBTMA_location <- as.integer(vapply(tok, function(x) x[3], character(1)))
img_meta$sample_id <- sub("\\.tiff$", "", img_meta$image)

img_meta <- merge(img_meta, patient, by = c("ArrayNr", "UBTMA_location"), all.x = TRUE)
idx <- match(spe$sample_id, img_meta$sample_id)
cols <- c("acquisition_description", "ArrayNr", "UBTMA_location", "PID", "TissueType", 
          "Mol_Signature", "Grade", "pN_4gr", "OS", "Livestatus1censored")
colData(spe) <- cbind(colData(spe), img_meta[idx, cols])

# Discard images without PID
img_no_pid <- unique(spe$sample_id[is.na(spe$PID)])
spe <- spe[, !is.na(spe$PID)]
message("Dropped ", length(img_no_pid), " images without a PID: ", paste(img_no_pid, collapse = ", "))

# Name cells as <sample_id>_<ObjectNumber>, e.g. ZTMA21.1_002_1
colnames(spe) <- paste(spe$sample_id, spe$ObjectNumber, sep="_")
stopifnot(!anyDuplicated(colnames(spe)))

saveRDS(spe, file = output_file, compress = FALSE)
message("Saved single-cell object to: ", output_file)
