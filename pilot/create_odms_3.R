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
  gene_odm <- create_ondisc_matrix_from_mtx(mtx_fp = mtx_fp, barcodes_fp = barcodes_fp,
                                            features_fp = gene_ids_fp, odm_fp = odm_fp,
                                            metadata_fp = metadata_fp, progress = TRUE)
  
  # Add p_mito and batch from cell_covariates data frame
  gasp_cell_covariates <- readRDS(paste0(intermediate_data_dir, "cell_covariates.rds"))
  gene_odm_plus_pmito_batch <- mutate_cell_covariates(gene_odm, p_mito = gasp_cell_covariates$percent.mito,
                                                      batch = factor(gasp_cell_covariates$prep_batch))
  
  # save the metadata (overwriting the original metadata file)
  save_odm(odm = gene_odm_plus_pmito_batch,
           metadata_fp = metadata_fp)
} else {
  gene_odm_plus_pmito_batch <- read_odm(odm_fp = odm_fp, metadata_fp = metadata_fp)
}

###########################
# 2. grna expression matrix
###########################
# next, load the grna count matrix. Write the backing .odm file and metadata file to disk.
odm_fp <- paste0(processed_grna_expression_dir, "matrix.odm")
metadata_fp <- paste0(processed_grna_expression_dir, "metadata.rds")
grna_count_matrix <- readRDS(paste0(intermediate_data_dir, "grna_count_matrix.rds"))
grna_feature_covariate_df <-  readRDS(paste0(intermediate_data_dir, "grna_feature_covariates.rds"))

# confirm that (1) cell barcodes of grna exp match those of gene odm, and (2) grna barcodes of grna exp match those of grna_feature_covariate_df
identical(get_cell_barcodes(gene_odm_plus_pmito_batch), colnames(grna_count_matrix))
grna_count_matrix <- grna_count_matrix[grna_feature_covariate_df$barcode,]
identical(grna_feature_covariate_df$barcode, rownames(grna_count_matrix))
cell_barcodes <- colnames(grna_count_matrix)
features_df <- data.frame(barcode = grna_feature_covariate_df$barcode)

# create the grna odm
grna_odm_exp <- create_ondisc_matrix_from_R_matrix(r_matrix = grna_count_matrix,
                                                   barcodes = cell_barcodes,
                                                   features_df = features_df,
                                                   odm_fp = odm_fp,
                                                   metadata_fp = metadata_fp)

# modify the feature covariates by adding the grna covariate information
grna_feature_covariate_df <- dplyr::mutate(grna_feature_covariate_df, barcode = NULL)
grna_odm_exp_mod <- grna_odm_exp %>%
  mutate_feature_covariates(coef_of_variation = NULL, grna_feature_covariate_df)

# overwrite metadata file
save_odm(grna_odm_exp_mod, metadata_fp)

############################################
# 3. grna assignment matrix (threshold at 5)
############################################
grna_assign_matrix <- grna_count_matrix >= 5
odm_fp <- paste0(processed_grna_assignment_dir, "matrix.odm")
metadata_fp <- paste0(processed_grna_assignment_dir, "metadata.rds")
grna_odm_assign <- create_ondisc_matrix_from_R_matrix(r_matrix = grna_assign_matrix,
                                                      barcodes = cell_barcodes,
                                                      features_df = features_df,
                                                      odm_fp = odm_fp,
                                                      metadata_fp = metadata_fp)
grna_odm_assign_mod <- grna_odm_assign %>%
  mutate_feature_covariates(coef_of_variation = NULL, grna_feature_covariate_df)

save_odm(odm = grna_odm_assign_mod, metadata_fp = metadata_fp)
