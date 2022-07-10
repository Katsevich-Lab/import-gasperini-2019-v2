###########################
# 0. Load packages; set fps
###########################
# source config file to get gasperini offsite location
gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/")
processed_data_dir <- paste0(gasp_offsite, "processed/")
processed_gene_dir_cp <- paste0(processed_data_dir, "gene_cp/")
if (!dir.exists(processed_gene_dir_cp)) dir.create(path = processed_gene_dir_cp, recursive = TRUE)

# set raw directories
raw_data_dir <- paste0(gasp_offsite, "raw/")

# load packages
library(ondisc)

###########################
# 1. gene expression matrix
###########################
# gene count matrix
mtx_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.exprs.mtx")
barcodes_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.cells.txt")
gene_ids_fp <- paste0(raw_data_dir, "GSE120861_at_scale_screen.genes.txt")

odm_fp <- paste0(processed_gene_dir_cp, "matrix.odm")
metadata_fp <- paste0(processed_gene_dir_cp, "metadata.rds")

# create the odm
if (!file.exists(odm_fp)) {
  gene_odm <- create_ondisc_matrix_from_mtx(mtx_fp = mtx_fp, barcodes_fp = barcodes_fp,
                                            features_fp = gene_ids_fp, odm_fp = odm_fp,
                                            metadata_fp = metadata_fp, progress = TRUE,
                                            n_lines_per_chunk = 1e8)
}
