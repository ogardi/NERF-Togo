# Rasterize RapidEye map

files <- dir("../input/RapidEye/shapefiles", pattern="\\.shp$", recursive=TRUE, full.names=TRUE)

foreach(file=files) %dopar%{
  out.file <- sub("shp$", "sdat", file)
  system(paste0("saga_cmd grid_gridding \"Shapes to Grid\" -TARGET_DEFINITION 0 ",
                "-INPUT \"", file, "\" -FIELD \"code\" -OUTPUT 2 -MULTIPLE 0 -LINE_TYPE 0 ",
                "-POLY_TYPE 0 -GRID_TYPE 2 -TARGET_USER_SIZE 30.0 -TARGET_USER_FITS 0 ",
                "-GRID \"", out.file, "\""))
}

scenes <- lapply(dir("../input/RapidEye/shapefiles", pattern="\\.sdat$", recursive=TRUE, full.names=TRUE), raster)

scenes["tolerance"] <- 0.4
RE <- do.call(merge, scenes)

AGB <- raster("../output/3_forest-biomass/2_ref-maps/TGO_2015_AGB.tif")

RE_resampled <- resample(RE, AGB, method="ngb")

writeRaster(RE_resampled, "../input/RapidEye/TGO_30m.tif", datatype="INT2S")

