####################################################################
# SSTS/data/prep-SRTM.R: cleaning and merging SRTM data
# ------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 11 March 2020


IN.DIR  <- paste0(DIR.RAW.DAT, "/SRTM")
OUT.DIR <- paste0(DIR.SST.DAT, "/SRTM")
if(!dir.exists(OUT.DIR)) dir.create(OUT.DIR)


# Prepare SRTM DEM ----------------------------

# Load 90m SRTM DEM (source: CGIAR),
dem.90 <- do.call(merge, lapply(as.list(dir(paste0(IN.DIR, "/3arcsecond"), pattern=".*[.]tif$", full.names=TRUE)), raster))


# Merge 30m SRTM DEM and fill voids with 90m SRTM (source: USGS)
dem.30 <- foreach(tile=dir(paste0(IN.DIR, "/1arcsecond"), pattern=".*[.]tif$", full.names=TRUE),
                  .combine=merge, .multicombine=TRUE) %dopar% {
                    dem.30.t <- raster(tile)
                    dem.90.t <- round(projectRaster(dem.90, dem.30.t))
                    merge(dem.30.t, dem.90.t)
                  }

# write it to the disk and reproject to reference Landsat image
writeRaster(dem.30, paste0(OUT.DIR, "/SRTM-1arcsec_raw.tif"), datatype="INT2S", overwrite=TRUE)
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

system(paste("gdalinfo -stats", paste0(OUT.DIR, "/SRTM-1arcsec.tif")))

dem <- raster(paste0(OUT.DIR, "/SRTM-1arcsec.tif"))

jpeg(paste0(OUT.DIR, "/SRTM-1arcsec.jpeg"), width=1350, height=3000)
par(plt=c(0,1,0,1))
plot(dem)
plot(mask(dem, TGO, inverse=TRUE), col="#FFFFFF66", legend=FALSE, add=TRUE)
plot(TGO, add=TRUE, lwd=3)
dev.off()