# source config file to get gasperini offsite location
gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/")

# create the intermediate data directory; set raw directory
intermediate_data_dir <- paste0(gasp_offsite, "intermediate/")
if (!dir.exists(intermediate_data_dir)) dir.create(path = intermediate_data_dir, recursive = TRUE)
raw_data_dir <- paste0(gasp_offsite, "raw/")

# Obtain binary gRNA matrix and cell metadata from monocole object
library(monocle)
library(magrittr)
monocle_obj <- readRDS(paste0(raw_data_dir, "/GSE120861_at_scale_screen.cds.rds"))
cell_metadata <- pData(monocle_obj)
rm(monocle_obj); gc()
covariates_cols <- 1:18
cell_covariates <- cell_metadata[,covariates_cols]

# save the cell covariates
saveRDS(cell_covariates, paste0(intermediate_data_dir, "cell_covariates.rds"))

# load the gRNA "groups" at scale and gRNA count matrix
cell_barcodes_in_use_long <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_at_scale_screen.cells.txt"),
                                        col_names = FALSE, col_types = "c") %>% dplyr::pull()
# all cell barcodes have 32 characters; strip the last 9
cell_barcodes_in_use <- gsub('.{9}$', '', cell_barcodes_in_use_long)
gRNA_barcodes_in_use <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_grna_groups.at_scale.txt"),
                                        col_names = c("group_name", "gRNA_barcode"), col_types = "cc")

gRNA_counts <- readr::read_tsv(paste0(raw_data_dir, "all_libraries.gRNAcaptured.aggregated.txt"),
                               col_names = TRUE,
                               col_types = "cccc") %>% dplyr::rename(cell_barcode = cell, gRNA_barcode = barcode)
# keep rows with cells in use, gRNA barcode in use, and nonzero UMI counts
gRNA_counts_sub <- gRNA_counts %>% dplyr::filter(umi_count > 0,
                                                 cell_barcode %in% cell_barcodes_in_use,
                                                 gRNA_barcode %in% gRNA_barcodes_in_use$gRNA_barcode)
# assign integer labels to the cell_barcodes and gRNA_barcodes
cell_idxs <- match(x = gRNA_counts_sub$cell_barcode, table = cell_barcodes_in_use)
gRNA_idxs <- match(x = gRNA_counts_sub$gRNA_barcode, table = gRNA_barcodes_in_use$gRNA_barcode)

m <- Matrix::sparseMatrix(i = gRNA_idxs,
                          j = cell_idxs,
                          x = as.integer(gRNA_counts_sub$umi_count),
                          dims = c(length(unique(gRNA_barcodes_in_use$gRNA_barcode)), length(unique(cell_barcodes_in_use))),
                          repr = "T")
row.names(m) <- gRNA_barcodes_in_use$gRNA_barcode
colnames(m) <- cell_barcodes_in_use_long

# save the count matrix to the intermediate file directory
saveRDS(object = m, file = paste0(intermediate_data_dir, "gRNA_count_matrix.rds"))
rm(gRNA_counts, gRNA_counts_sub)

# create the data frame of gRNA feature covariates
gRNA_id_to_group_df <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_grna_groups.at_scale.txt"),
                                       col_types = "cc", col_names = c("gRNA_group", "barcode"))
gRNA_result_table <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_all_deg_results.at_scale.txt"))
gRNA_group_to_target_df <- gRNA_result_table |>
  dplyr::select(gRNA_group, target_site.chr, target_site.start, target_site.stop, target_type = site_type) |>
  dplyr::distinct() |>
  dplyr::mutate(target = paste0(target_site.chr, ":", target_site.start, "-", target_site.stop),
                target_site.chr  = NULL, target_site.start = NULL, target_site.stop = NULL) |>
  dplyr::filter(target_type != "TSS") |>
  dplyr::mutate(target_type = forcats::fct_recode(target_type, gene_tss = "selfTSS",
                                                  candidate_enhancer = "DHS",
                                                  known_enhancer = "positive_ctrl",
                                                  "non-targeting" = "NTC"),
                target = ifelse(target_type != "non-targeting", target, "non-targeting"))
pos_control_gRNA_group_to_target_gene_df <- gRNA_result_table |>
  dplyr::filter(site_type == "selfTSS") |>
  dplyr::pull(pairs4merge) |>
  strsplit(":") |> 
  do.call(what = rbind, args = _) |>
  as.data.frame() |>
  setNames(c("gRNA_group", "target_gene"))

# perform a join operation to create a data frame with columns gRNA ID (barcode), target (the chromosomal position that a gRNA targets), target_type, pos_control_gene (for positive controls, the targeted gene), and gRNA_group
gRNA_group_tbl <- dplyr::left_join(gRNA_group_to_target_df, pos_control_gRNA_group_to_target_gene_df, by = "gRNA_group")
gRNA_feature_covariates <- dplyr::left_join(gRNA_id_to_group_df, gRNA_group_tbl,by = "gRNA_group")
saveRDS(gRNA_feature_covariates, paste0(intermediate_data_dir, "gRNA_feature_covariates.rds"))
