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

foreach(file=dir(IN.DIR, pattern=".*Togo[.]tif$")) %dopar% {
  system(paste("gdalwarp",
               paste0(IN.DIR, "/", file),
               "-t_srs '+proj=utm +zone=31 +datum=WGS84'",
               "-tr 30 30",
               paste("-te", TGO.EXT@xmin, TGO.EXT@ymin, TGO.EXT@xmax, TGO.EXT@ymax),
               paste0(OUT.DIR, "/", file),
               "-dstnodata -3.4e+38",
               "-co COMPRESS='LZW'",
               "-co INTERLEAVE='BAND'",
               "-overwrite"))
  
  system(paste("gdalinfo -stats",
               paste0(OUT.DIR, "/", file)))
}

# Créer des vignettes pour les donées WorldClim ===============================

foreach(file=dir(OUT.DIR, pattern=".*[.]tif$")) %dopar% {
  image <- stack(paste0(OUT.DIR, "/", file))
  type <- unlist(strsplit(file, "_"))[3]
  if      (type == "prec") { zlim <- c(0,320);     col <- rev(topo.colors(255)) }
  else if (type == "tmin") { zlim <- c(14.0,27.8); col <- rev(heat.colors(255)) }
  else if (type == "tmax") { zlim <- c(24.9,37.5); col <- rev(heat.colors(255)) }
  else if (type == "tavg") { zlim <- c(19.7,32.7); col <- rev(heat.colors(255)) }
  else                     { zlim <- NA;           col <- rev(cm.colors(255)) }
  foreach(i=1:nlayers(image)) %dopar% {
    jpeg(paste0(OUT.DIR, "/", sub("[.]tif$", "", file), "-", str_pad(i, 2, "left", 0), ".jpeg"), 
         width=1350, height=3000)
    plot(image[[i]], col=col, zlim=zlim)
    plot(mask(image[[i]], TGO, inverse=TRUE), col="#FFFFFF66", legend=FALSE, add=TRUE)
    plot(TGO, add=TRUE, lwd=3)
    dev.off()
  }
}
