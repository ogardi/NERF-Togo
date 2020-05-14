###############################################################################
# prep-ProREDD.R: rasterizer la carte d'occupation des terres RapidEye
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Définitions des variables ===================================================

IN.DIR  <- paste0(DIR.RAW.DAT, "/RapidEye/shapefiles")
OUT.DIR <- paste0(DIR.SST.DAT, "/ProREDD")
if(!dir.exists(OUT.DIR)) dir.create(OUT.DIR)

# Rasterizer Shapefiles avec l'outil grid_gridding de SAGA ====================

files <- dir(IN.DIR, pattern="\\.shp$", recursive=TRUE, full.names=TRUE)
registerDoParallel(CORES - 1)
foreach(file=files) %dopar%{
  out.file <- sub("shp$", "sdat", file)
  system(paste0("saga_cmd grid_gridding \"Shapes to Grid\" ",
                "-TARGET_DEFINITION 0 -INPUT \"", file, "\" ",
                "-FIELD \"code\" -OUTPUT 2 -MULTIPLE 0 -LINE_TYPE 0 ",
                "-POLY_TYPE 0 -GRID_TYPE 2 -TARGET_USER_SIZE 30.0 ",
                "-TARGET_USER_FITS 0 -GRID \"", out.file, "\""))
}

# Lire et fusionner les différentes scènes raster
scenes <- lapply(dir(IN.DIR, pattern="\\.sdat$", recursive=TRUE, full.names=TRUE), raster)
scenes["tolerance"] <- 0.4
RE <- do.call(merge, scenes)

# À reviser: Reprojection sur l'image Landsat
AGB <- raster("../output/3_forest-biomass/2_ref-maps/TGO_2015_AGB.tif")
RE_resampled <- resample(RE, AGB, method="ngb")
writeRaster(RE_resampled, paste0(OUT.DIR, "/TGO_30m.tif"), datatype="INT2S")