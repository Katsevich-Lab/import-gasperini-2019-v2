#######################
# 0. Load data; set fps
#######################

gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/")
raw_dir <- paste0(gasp_offsite, "raw/")
processed_dir <- paste0(gasp_offsite, "processed/")

# read gRNA-gene pairs and gRNA info
all_results <- readr::read_tsv(file = paste0(raw_dir, "GSE120861_all_deg_results.at_scale.txt"),
                               col_names = TRUE, col_types = "c")
pairs_to_analyze <- all_results %>%
  dplyr::select(gene_id = ENSG, gRNA_group, site_type)
saveRDS(pairs_to_analyze, paste0(processed_dir, "pairs_ungrouped.rds"))
