# source config file to get gasperini offsite location
gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "high-MOI-pilot/")

# create the intermediate data directory; set raw directory
intermediate_data_dir <- paste0(gasp_offsite, "intermediate/")
if (!dir.exists(intermediate_data_dir)) dir.create(path = intermediate_data_dir, recursive = TRUE)
raw_data_dir <- paste0(gasp_offsite, "raw/")

# Obtain binary grna matrix and cell metadata from monocole object
library(monocle)
library(magrittr)
monocle_obj <- readRDS(paste0(raw_data_dir, "/GSE120861_pilot_highmoi_screen.cds.rds"))
cell_metadata <- pData(monocle_obj)
rm(monocle_obj); gc()

# save the cell covariates
covariates_cols <- 1:14
cell_covariates <- cell_metadata[, covariates_cols]
saveRDS(cell_covariates, paste0(intermediate_data_dir, "cell_covariates.rds"))

# load the grna and cell barcodes information
cell_barcodes_in_use <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_pilot_highmoi_screen.cells.txt"),
                                        col_names = FALSE, col_types = "c") %>% dplyr::pull()

grna_barcodes_in_use <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_grna_groups.pilot.txt"),
                                        col_names = c("group_name", "grna_barcode"), col_types = "cc")

# obtain the grna assignment matrix
# Split the strings by "_"
split_data <- strsplit(cell_metadata$barcode, "_")
names(split_data) <- rownames(cell_metadata)

# Create an empty sparse matrix (letters as rows, original rows as columns)
grna_assignment <- Matrix(0, nrow = length(unique(grna_barcodes_in_use$grna_barcode)), 
                          ncol = length(unique(cell_barcodes_in_use)), sparse = TRUE)

# Assign row names (letters)
rownames(grna_assignment) <- unique(grna_barcodes_in_use$grna_barcode)
colnames(grna_assignment) <- cell_barcodes_in_use

# Fill the matrix: Set 1 where the letter appears in the corresponding column
for (cell_id in colnames(grna_assignment)) {
  if(any(is.na(split_data[[cell_id]]))){
    grna_assignment[, cell_id] <- 0
  }else{
    grna_assignment[split_data[[cell_id]], cell_id] <- 1
  }
}

# save the count matrix to the intermediate file directory
saveRDS(object = grna_assignment, file = paste0(intermediate_data_dir, "grna_assignment_matrix.rds"))
rm(grna_assignment)

# create the data frame of grna feature covariates
grna_id_to_group_df <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_grna_groups.pilot.txt"),
                                       col_types = "cc", col_names = c("grna_group", "barcode"))
grna_result_table <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_all_deg_results.pilot.txt"))
grna_group_to_target_df <- grna_result_table |>
  dplyr::select(grna_group = gRNA_group, target_type = site_type) |>
  dplyr::distinct() |>
  dplyr::filter(target_type != "TSS") |>
  dplyr::mutate(target_type = forcats::fct_recode(target_type, gene_tss = "selfTSS",
                                                  candidate_enhancer = "DHS",
                                                  known_enhancer = "positive_ctrl",
                                                  "non-targeting" = "NTC"))
pos_control_grna_group_to_target_gene_df <- grna_result_table |>
  dplyr::filter(site_type == "selfTSS") |>
  dplyr::pull(pairs4merge) |>
  strsplit(":") |> 
  do.call(what = rbind, args = _) |>
  as.data.frame() |>
  setNames(c("grna_group", "target_gene"))

# perform a join operation to create a data frame with columns grna ID (barcode), target (the chromosomal position that a grna targets), target_type, pos_control_gene (for positive controls, the targeted gene), and grna_group
grna_group_tbl <- dplyr::left_join(grna_group_to_target_df, pos_control_grna_group_to_target_gene_df, by = "grna_group")
grna_feature_covariates <- dplyr::left_join(grna_group_tbl, grna_id_to_group_df, by = "grna_group")
saveRDS(grna_feature_covariates, paste0(intermediate_data_dir, "grna_feature_covariates.rds"))
