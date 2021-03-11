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
  out.file <- paste0(OUT.DIR, "/OCCSOL_5m/", sub(".*_", "", dirname(file)), "-", sub(".shp$", ".sdat", basename(file)))
  system(paste0("saga_cmd grid_gridding \"Shapes to Grid\" ",
                "-TARGET_DEFINITION 0 -INPUT \"", file, "\" ",
                "-FIELD \"code\" -OUTPUT 2 -MULTIPLE 0 -LINE_TYPE 0 ",
                "-POLY_TYPE 0 -GRID_TYPE 2 -TARGET_USER_SIZE 5.0 ",
                "-TARGET_USER_FITS 1 -GRID \"", out.file, "\""))
}

# build virtual raster mosaic
system(paste("gdalbuildvrt", paste0(OUT.DIR, "/", "OCCSOL_5m.vrt"), paste0(OUT.DIR, "/OCCSOL_5m/", "*.sdat")))

# check whether rasterization was done correctly
# r <- raster(paste0(OUT.DIR, "/", "OCCSOL_5m.vrt"))
# raster::plot(r, ext=extent(250200, 250400, 900100, 900300))
# s <- readOGR("../data/RapidEye/Class_Z_merge.shp")
# sc <- crop(s, extent(250200, 250400, 900100, 900300))
# raster::plot(sc, add=TRUE)

system(paste("gdalinfo", paste0(OUT.DIR, "/", "OCCSOL_5m.vrt")))

# export to tiff with 30m resolution and aligned with Landsat
system(paste("gdalwarp -t_srs EPSG:32631",
               "-overwrite",
               "-tr 30 30",
               "-r mode",
               "-te", TGO.EXT@xmin,  TGO.EXT@ymin, TGO.EXT@xmax, TGO.EXT@ymax,
               paste0(OUT.DIR, "/", "OCCSOL_5m.vrt"), 
               "-ot Int16",
               "-dstnodata -9999",
               paste0(OUT.DIR, "/", "OCCSOL_30m.tif")))

system(paste("gdalinfo", paste0(OUT.DIR, "/", "OCCSOL_30m.tif")))


# read raster, add raster attribute table, and plot it
r <- as.factor(raster(paste0(OUT.DIR, "/", "OCCSOL_30m.tif")))
rat <- levels(r)[[1]]
rat[["landcover"]] <- rat[["ID"]]
attr <- as.data.frame(matrix(ncol=4, byrow=TRUE, c(
  0,  "#FFFFFF", "-",                                   "-",
  1,  "#000000", "TERRES FORESTIERES",                  "-",
  11, "#184A14", "Forêts denses",                       "Forêt",
  12, "#ACFC04", "Forêts riveraines",                   "Forêt",
  13, "#38C403", "Forêts claires",                      "Forêt",
  14, "#BADD69", "Savane arborée/arbustive",            "Hors Forêt",
  15, "#000000", "Mangroves",                           "Hors Forêt",
  16, "#51741B", "Plantations",                         "Forêt",
#  17, "#000000", "Fourrés",                            "Hors Forêt",
  18, "#000000", "Savanes boisées",                     "Forêt",
  19, "#000000", "Formations marécageuses",             "Forêt",
  20, "#000000", "TERRES CULTIVEES",                    "-",
  21, "#EEEA04", "Cultures et Jachères",                "Hors Forêt",
  22, "#000000", "Cultures sans arbres",                "Hors Forêt",
  23, "#000000", "Parcs agroforestiers",                "Hors Forêt",
  30, "#000000", "FORMATIONS HERBEUSES",                "-",
#  31, "#000000", "Prairies",                           "Hors Forêt",
  32, "#1CD791", "Savanes herbeuses",                   "Hors Forêt",
  34, "#FFFFFF", "???",                                 "-",
  40, "#000000", "ETABLISSEMENTS",                      "-",
  41, "#E31A1C", "Agglomérations et infrastructures",   "Hors Forêt",
  42, "#000000", "Plantations urbaines",                "Hors Forêt",
  43, "#FFFFFF", "???",                                 "-",
  45, "#FFFFFF", "???",                                 "-",
#  5, "#000000", "TERRES HUMIDES",                      "-",
  51, "#000000", "Plans d’eau et rivières",             "Hors Forêt",
  60, "#000000", "AUTRES TERRES",                       "-",
  61, "#D1C4C4", "Sols nus",                            "Hors Forêt",
  90, "#FFFFFF", "???",                                 "-",
  92, "#FFFFFF", "???",                                 "-",
  99, "#FFFFFF", "???",                                 "-",
 100, "#FFFFFF", "???",                                 "-",
 121, "#FFFFFF", "???",                                 "-",
 149, "#FFFFFF", "???",                                 "-",
 162, "#FFFFFF", "???",                                 "-",
 216, "#FFFFFF", "???",                                 "-",
 512, "#FFFFFF", "???",                                 "-",
 612, "#FFFFFF", "???",                                 "-")))
names(attr) <- c("CODE", "COL", "OCCSOL", "TYPE")
attr$CODE <- as.integer(attr$CODE)
attr$COL <- as.character(attr$COL)
attr$OCCSOL <- as.character(attr$OCCSOL)
levels(r) <- cbind(rat, attr)
levelplot(r, att="OCCSOL", col.regions=attr$COL, xlab="", ylab="")