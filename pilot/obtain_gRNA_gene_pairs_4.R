#######################
# 0. Load data; set fps
#######################

gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "high-MOI-pilot/")
raw_dir <- paste0(gasp_offsite, "raw/")
processed_dir <- paste0(gasp_offsite, "processed/")

# read gRNA-gene pairs and gRNA info
all_results <- readr::read_tsv(file = paste0(raw_dir, "GSE120861_all_deg_results.pilot.txt"),
                               col_names = TRUE, col_types = "c")
pairs_to_analyze <- all_results |>
  dplyr::select(gene_id = ENSG, grna_group = gRNA_group, site_type) |>
  dplyr::mutate(site_type = forcats::fct_recode(site_type, "cis" = "DHS", "neg_cntrl" = "NTC",
                                                "enhancer_pos_cntrl" = "positive_ctrl",
                                                "pos_cntrl" = "selfTSS", "other" = "TSS"))
saveRDS(pairs_to_analyze, paste0(processed_dir, "pairs_grouped.rds"))
