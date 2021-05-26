###############################################################################
# prep-SRTM.R: lire, nettoyer et empiler des données topographiques SRTM
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Définitions des variables ===================================================
IN.DIR  <- paste0(DIR.RAW.DAT, "/SRTM")
OUT.DIR <- paste0(DIR.SST.DAT, "/SRTM")
if(!dir.exists(OUT.DIR)) dir.create(OUT.DIR)


# Préparation du modèle numérique d'élévation SRTM ============================

# Lire 30m SRTM MNE (source: USGS Earthexplorer), et
# remplir vides avec 90m SRTM (source: http://srtm.csi.cgiar.org/srtmdata/)

system(paste("gdalbuildvrt", paste0(OUT.DIR, "/SRTM_1as.vrt"), 
             "-resolution highest",
             paste0(dir(paste0(IN.DIR, "/3arcsecond"), 
                        pattern=".*[.]tif$", full.names=TRUE), collapse = " "),
             paste0(dir(paste0(IN.DIR, "/1arcsecond"), 
                        pattern=".*[.]tif$", full.names=TRUE), collapse = " ")))

# Reprojection vers Landsat (résolution 30m, UTM 31, ...) 
system(paste("gdalwarp -t_srs EPSG:32631",
             paste0(OUT.DIR, "/SRTM_1as.vrt"),
             "-tr 30 30",
             paste("-te", TGO.EXT@xmin, TGO.EXT@ymin, TGO.EXT@xmax, TGO.EXT@ymax),
             paste0(OUT.DIR, "/SRTM_30m.tif"),
             "-r bilinear",             # interpolation bilinéaire
             "-ot Int16",
             "-dstnodata -9999",
             "-overwrite"))


# Créer une vignettes du MNE ==================================================
dem <- raster(paste0(OUT.DIR, "/SRTM_30m.tif"))
jpeg(paste0(OUT.DIR, "/SRTM_30m.jpeg"), width=570, height=1200)
par(mai=c(0,0,0,0))
raster::plot(TGO, col = "purple")
raster::plot(dem, add = TRUE, bgalpha = 0)
raster::plot(mask(dem, TGO, inverse=TRUE), col="#FFFFFF66", legend=FALSE, add=TRUE)
raster::plot(TGO, add=TRUE, lwd=3)
dev.off()
