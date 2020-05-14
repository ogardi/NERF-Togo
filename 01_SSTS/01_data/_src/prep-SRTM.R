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

# Lire 90m SRTM MNE sans vides (source: http://srtm.csi.cgiar.org/srtmdata/),
dem.90 <- do.call(merge, lapply(as.list(dir(paste0(IN.DIR, "/3arcsecond"), 
                                            pattern=".*[.]tif$", 
                                            full.names=TRUE)), raster))


# Lire 30m SRTM MNE (source: USGS Earthexplorer) et remplir vides avec 90m SRTM 
dem.30 <- foreach(tile=dir(paste0(IN.DIR, "/1arcsecond"), pattern=".*[.]tif$", full.names=TRUE),
                  .combine=merge, .multicombine=TRUE) %dopar% {
                    dem.30.t <- raster(tile)
                    dem.90.t <- round(projectRaster(dem.90, dem.30.t))
                    merge(dem.30.t, dem.90.t)
                  }

# Reprojection d'images MNE vers Landsat (résolution 30m, UTM 31, ...) 
writeRaster(dem.30, paste0(OUT.DIR, "/SRTM-1arcsec_raw.tif"), 
            datatype="INT2S", 
            overwrite=TRUE)
system(paste("gdalwarp",
             paste0(OUT.DIR, "/SRTM-1arcsec_raw.tif"),
             "-t_srs '+proj=utm +zone=31 +datum=WGS84'",
             "-tr 30 30",
             paste("-te", TGO.EXT@xmin, TGO.EXT@ymin, TGO.EXT@xmax, TGO.EXT@ymax),
             paste0(OUT.DIR, "/SRTM-1arcsec.tif"),
             "-ot 'Int16'",
             "-co COMPRESS='LZW'",
             "-co INTERLEAVE='BAND'",
             "-overwrite"))
file.remove(paste0(OUT.DIR, "/SRTM-1arcsec_raw.tif"))

# Calculer les statistiques GDAL (min, max, ...)
system(paste("gdalinfo -stats", paste0(OUT.DIR, "/SRTM-1arcsec.tif")))


# Créer une vignettes du MNE ==================================================
dem <- raster(paste0(OUT.DIR, "/SRTM-1arcsec.tif"))
jpeg(paste0(OUT.DIR, "/SRTM-1arcsec.jpeg"), width=1350, height=3000)
par(plt=c(0,1,0,1))
plot(dem)
plot(mask(dem, TGO, inverse=TRUE), col="#FFFFFF66", legend=FALSE, add=TRUE)
plot(TGO, add=TRUE, lwd=3)
dev.off()