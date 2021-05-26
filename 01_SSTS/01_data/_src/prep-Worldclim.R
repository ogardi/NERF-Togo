###############################################################################
# prep-Worldclim.R: lire et reprojeter les données de WorldClim
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Définitions des variables ===================================================

IN.DIR  <- paste0(DIR.RAW.DAT, "/Worldclim")
OUT.DIR <- paste0(DIR.SST.DAT, "/Worldclim")
if(!dir.exists(OUT.DIR)) dir.create(OUT.DIR)


# Reprojection images WorldClim vers Landsat (résolution 30m, UTM 31, ...) ====

registerDoParallel(CORES-1)
foreach(file=dir(IN.DIR, pattern=".*Togo[.]tif$")) %dopar% {
  
  tmp <- paste0(OUT.DIR, "/", sub("[.]tif", "_tmp.tif", file))
  
  # Multiplier avec 10 et convertir on format integer
  writeRaster(10*raster::stack(paste0(IN.DIR, "/", file)), tmp, datatype="INT2U")

  system(paste("gdalwarp -t_srs EPSG:32631",
               tmp,
               "-tr 30 30",
               paste("-te", TGO.EXT@xmin, TGO.EXT@ymin, TGO.EXT@xmax, TGO.EXT@ymax),
               paste0(OUT.DIR, "/", file),
               "-r bilinear",             # interpolation bilinéaire
               "-co COMPRESS='LZW'",
               "-co INTERLEAVE='BAND'",
               "-overwrite"))
  
  system(paste("gdalinfo -mm", paste0(OUT.DIR, "/", file)))
  
  unlink(tmp)
}


# Créer des vignettes pour les donées WorldClim ===============================

foreach(file=dir(OUT.DIR, pattern=".*[.]tif$")) %dopar% {
  image <- raster::stack(paste0(OUT.DIR, "/", file))
  type <- unlist(strsplit(file, "_"))[3]
  if      (type == "prec") { zlim <- c(0,3200);     col <- rev(topo.colors(255)) }
  else if (type == "tmin") { zlim <- c(140,278); col <- rev(heat.colors(255)) }
  else if (type == "tmax") { zlim <- c(249,375); col <- rev(heat.colors(255)) }
  else if (type == "tavg") { zlim <- c(197,327); col <- rev(heat.colors(255)) }
  else                     { zlim <- NA;           col <- rev(cm.colors(255)) }
  foreach(i=1:nlayers(image)) %dopar% {
    jpeg(paste0(OUT.DIR, "/", sub("[.]tif$", "", file), "-", str_pad(i, 2, "left", 0), ".jpeg"), 
         width=570, height=1200)
 #   par(mai=c(0,0,0,0))
    raster::plot(image[[i]], col=col, zlim=zlim)
    raster::plot(mask(image[[i]], TGO, inverse=TRUE), col="#FFFFFF66", legend=FALSE, add=TRUE)
    raster::plot(TGO, add=TRUE, lwd=3)
    dev.off()
  }
}


# Découpage en chemins WRS ====================================================

foreach(file=dir(OUT.DIR, pattern=".*Togo[.]tif$")) %:%
  foreach(path=names(WRS.EXT)) %dopar% {
    system(paste("gdalwarp",
                 paste0(OUT.DIR, "/", file),
                 "-tr 30 30",
                 paste("-te", WRS.EXT[[path]]@xmin, WRS.EXT[[path]]@ymin, WRS.EXT[[path]]@xmax, WRS.EXT[[path]]@ymax),
                 paste0(OUT.DIR, "/", sub("Togo", path, file)),
                 "-co COMPRESS='LZW'",
                 "-co INTERLEAVE='BAND'",
                 "-overwrite"))
}
