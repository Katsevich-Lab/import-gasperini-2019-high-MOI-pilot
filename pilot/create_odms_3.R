###########################
# 0. Load packages; set fps
###########################
# source config file to get gasperini offsite location
gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "high-MOI-pilot/")
processed_data_dir <- paste0(gasp_offsite, "processed/")
processed_gene_dir <- paste0(processed_data_dir, "gene/")
processed_grna_expression_dir <- paste0(processed_data_dir, "grna_expression/")
processed_grna_assignment_dir <- paste0(processed_data_dir, "grna_assignment/")

dirs_to_create <- c(processed_data_dir, processed_gene_dir,
                    processed_grna_expression_dir, processed_grna_assignment_dir)
for (dir in dirs_to_create) {
  if (!dir.exists(dir)) dir.create(path = dir, recursive = TRUE)
}

# set raw directories
raw_data_dir <- paste0(gasp_offsite, "raw/")
intermediate_data_dir <- paste0(gasp_offsite, "intermediate/")

# load packages
library(magrittr)
library(ondisc)

###########################
# 1. gene expression matrix
###########################
# gene count matrix
mtx_fp <- paste0(raw_data_dir, "GSE120861_pilot_highmoi_screen.exprs.mtx")
barcodes_fp <- paste0(raw_data_dir, "GSE120861_pilot_highmoi_screen.cells.txt")
gene_ids_fp <- paste0(raw_data_dir, "GSE120861_pilot_highmoi_screen.genes.txt")

odm_fp <- paste0(processed_gene_dir, "matrix.odm")
metadata_fp <- paste0(processed_gene_dir, "metadata.rds")

# create the odm
if (!file.exists(odm_fp)) {
  gene_odm <- ondisc:::create_ondisc_matrix_from_mtx(mtx_fp = mtx_fp, barcodes_fp = barcodes_fp,
                                                     features_fp = gene_ids_fp, odm_fp = odm_fp,
                                                     metadata_fp = metadata_fp, progress = TRUE)
  
  # save the metadata (overwriting the original metadata file)
  save_odm(odm = gene_odm, metadata_fp = metadata_fp)
} else {
  gene_odm <- read_odm(odm_fp = odm_fp, metadata_fp = metadata_fp)
}

############################################
# 3. grna assignment matrix (provided in monocle object in GEO)
############################################
# next, load the grna count matrix. Write the backing .odm file and metadata file to disk.
odm_fp <- paste0(processed_grna_expression_dir, "matrix.odm")
metadata_fp <- paste0(processed_grna_expression_dir, "metadata.rds")
grna_assignment_matrix <- readRDS(paste0(intermediate_data_dir, "grna_assignment_matrix.rds"))
grna_feature_covariate_df <-  readRDS(paste0(intermediate_data_dir, "grna_feature_covariates.rds"))

# confirm that (1) cell barcodes of grna exp match those of gene odm, and (2) grna barcodes of grna exp match those of grna_feature_covariate_df
identical(ondisc:::get_cell_barcodes(gene_odm), colnames(grna_assignment_matrix))
grna_assignment_matrix <- grna_assignment_matrix[unique(grna_feature_covariate_df$barcode),]
identical(unique(grna_feature_covariate_df$barcode), unique(rownames(grna_assignment_matrix)))
cell_barcodes <- colnames(grna_assignment_matrix)
features_df <- data.frame(barcode = unique(grna_feature_covariate_df$barcode))

# skip the odm construction when it is already there
if(!file.exists(odm_fp)){
  # construct odm type of file
  odm_fp <- paste0(processed_grna_assignment_dir, "matrix.odm")
  metadata_fp <- paste0(processed_grna_assignment_dir, "metadata.rds")
  grna_odm_assign <- ondisc:::create_ondisc_matrix_from_R_matrix(r_matrix = grna_assignment_matrix,
                                                                 barcodes = cell_barcodes,
                                                                 features_df = features_df,
                                                                 odm_fp = odm_fp,
                                                                 metadata_fp = metadata_fp)
  grna_odm_assign_mod <- grna_odm_assign %>%
    mutate_feature_covariates(coef_of_variation = NULL, grna_feature_covariate_df)
  
  # save the odm file
  save_odm(odm = grna_odm_assign_mod, metadata_fp = metadata_fp)
}