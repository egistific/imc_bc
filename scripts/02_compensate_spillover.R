log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

library(SpatialExperiment)
library(CATALYST)

config <- snakemake@config
cofactor <- config$analysis$arcsinh_cofactor
method = config$analysis$compensation_method
output_file <- snakemake@output[["spe"]]
spe <- readRDS(snakemake@input[["spe"]])
sm <- as.matrix(read.csv(snakemake@input[["sm"]], row.names = 1, check.names = FALSE))

# Discard ArAr80 as it is not in CATALYST's isotope list
sm <- sm[, colnames(sm) != "ArAr80Di"]

rowData(spe)$channel_name <- paste0(rowData(spe)$metal, "Di")
stopifnot(all(rowData(spe)$channel_name %in% colnames(sm)))
message("Channels that do not have a spillover row: ", 
        paste(setdiff(rowData(spe)$channel_name, rownames(sm)), collapse = ", "))

# Cofactor required by compCytof's argument check but unused with transform = FALSE
spe <- compCytof(spe, sm, transform = FALSE, cofactor = cofactor, method = method, overwrite = FALSE)

# Check that there are no NA values after spillover compensation
stopifnot(!anyNA(assay(spe, "compcounts")))
message("Number of negative values after spillover compensation: ", sum(assay(spe, "compcounts") < 0))

saveRDS(spe, file = output_file, compress = FALSE)
message("Saved single-cell object after applying spillover compensation to: ", output_file)

