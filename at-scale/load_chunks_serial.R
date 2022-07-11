gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/")
test_dir <- paste0(gasp_offsite, "raw/test")
fs <- list.files(test_dir, full.names = TRUE)

for (f in fs) {
  print(paste0("Loading ", f))
  # x <- readr::read_delim(file = f, delim = " ", skip = 2, col_types = c("iii"))
  x <- data.table::fread(file = f,
                         sep = " ",
                         colClasses = c("integer", "integer", "integer"),
                         skip = 2)
  rm(x)
  gc()
}