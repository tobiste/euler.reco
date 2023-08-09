## code to prepare `DATASET` dataset goes here
tintina <- readRDS("C:/Users/tstephan/Documents/GitHub/cordillera-stress/data/faults.rds") |>
  dplyr::filter(Name %in% c("Rocky Mountain Trench", "Tintina Fault"))

usethis::use_data(tintina, overwrite = TRUE)


south_atlantic <- sf::read_sf("data-raw/south_atlantic.kml")
usethis::use_data(south_atlantic, overwrite = TRUE)
