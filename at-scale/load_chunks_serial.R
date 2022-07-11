gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/")
test_dir <- paste0(gasp_offsite, "raw/test")
fs <- list.files(test_dir, full.names = TRUE)

for (f in fs) {
  x <- readr::read_delim(file = f, delim = " ", skip = 2, col_types = c("iii"))
}