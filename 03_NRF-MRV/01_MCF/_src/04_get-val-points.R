##########################################################################
# NERF_Togo/FCC/5_create-val-points.R: validate clean forest cover maps
# ------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 20 May 2019


# merge the 1987, 2003, 2015 and 2018 maps
maps.p192 <- do.call(stack, lapply(dir(paste0(OUTPUT.DIR, "/2_forest-maps/TGO/2_raw-maps/p192"), pattern="\\.tif", full.names=TRUE)[c(2,5,10,12)], raster))
maps.p193 <- do.call(stack, lapply(dir(paste0(OUTPUT.DIR, "/2_forest-maps/TGO/2_raw-maps/p193"), pattern="\\.tif", full.names=TRUE)[c(2,6,11,13)], raster))
maps.p194 <- do.call(stack, lapply(dir(paste0(OUTPUT.DIR, "/2_forest-maps/TGO/2_raw-maps/p194"), pattern="\\.tif", full.names=TRUE)[c(2,5,10,12)], raster))
maps <- merge(maps.p193, maps.p194, maps.p192)

writeRaster(maps, filename=paste0(VALIDTN.DIR, "/1_val4k-s16trans/TGO.tif"), bylayer=TRUE, suffix=paste0(c(1987, 2003, 2015, 2018), "_FC30r"))

# maps <- stack(dir(MAPS.DIR, pattern=paste0(".*c", MAPS.SUFFIX, "[.]tif$"), full.names=TRUE))[[c(1,4,12,14)]]
change.map <- maps[[1]] + maps[[2]]*10 + maps[[3]]*100 + maps[[4]]*1000


# Allocation of validation points
# ================================

n <- 4000
alloc <- freq(change.map)[1:16,]
alloc <- cbind(alloc, prop=round(n*alloc[,"count"]/sum(alloc[,"count"])), equal=round(n/nrow(alloc)))
alloc <- cbind(alloc, balanced=round((alloc[,"prop"] + alloc[,"equal"])/ 2))
sum(alloc[,"balanced"])


# Read the frame
# ==============
frame.points        <- readOGR(paste0(INPUT.DIR, "/SSTS/TGO_frame_480m.shp"))

# extract forest transitions
frame.points$trans  <- raster::extract(change.map, frame.points)

# writeOGR(frame.points, dsn=paste0(INPUT.DIR, "/Train-Val/201908"), layer="TGO_frame_480m", driver="ESRI Shapefile", overwrite=TRUE)


# merge with training-plots data
train.plots         <- readOGR(paste0(REFMAPS.DIR, "/1_trn10k-s10ndvi/assessed_v1/COV_parcelles.shp"))
frame.points        <- merge(frame.points, train.plots@data[,c("PLOTID", "img_date", "img_src", "mod_date", "author", "ccov")], by="PLOTID", all.x=TRUE)

# initiate the sample
val.points <- frame.points[0,]

# Initialize random number generator
set.seed(1)

# for each stratum ...
for(i in 1:nrow(alloc)) {
  
  strat.n <-  alloc[i, "balanced"]   # stratum sample size
  sample.ids <- NULL                 # vectors with ids to sample
  
  # take as much samples already used as training-points
  ids.tp <- which(!is.na(frame.points$trans) & frame.points$trans==alloc[i, "value"] & !is.na(frame.points$ccov))
  n.tp   <- min(length(ids.tp), strat.n)
  if (n.tp > 0) sample.ids <- c(sample.ids, sample(ids.tp, n.tp))
  
  # and complete with other samples from the grid
  ids.r  <- which(!is.na(frame.points$trans) & frame.points$trans==alloc[i, "value"] & is.na(frame.points$ccov))
  n.r    <- min(length(ids.r), strat.n - n.tp)
  if (n.r > 0) sample.ids <- c(sample.ids, sample(ids.r, n.r)) 
  
  val.points <- rbind(val.points, frame.points[sample.ids, ])
}

# shuffle the data
val.points <- val.points[sample(1:nrow(val.points)), ]
# add ID
val.points$SAMPLEID <- paste0("val-", str_pad(string=1:nrow(val.points), width = 4, pad = "0", side = "left"))


# Convert to plots and add attributes
# ===================================

landsat.grid <- raster(change.map)
values(landsat.grid) <- 1

# convert to polygons
val.plots <- rasterToPolygons(mask(landsat.grid, val.points))

# fetch attributes
val.plots@data <- over(val.plots, val.points[, c("PLOTID", "SAMPLEID", "xcoords", "ycoords", "trans", "ccov", "img_src", "img_date", "author", "mod_date")])

# add attributes
val.plots$ccov                          <- format(round(100*val.plots$ccov, 1))
val.plots$ccov[val.plots$ccov == "  NA"]<- NA
val.plots$lc_18   <- val.plots$lc_15    <- as.character(NA)
val.plots$lc_03   <- val.plots$lc_87    <- as.character(NA)

# write plots as Shapefile and KML
# writeOGR(val.plots, dsn=paste0(INPUT.DIR, "/Train-Val/201908/val9k-s16trans/empty"), layer="UOT_parcelles", driver="ESRI Shapefile", overwrite=TRUE)
# writeKML(val.plots, kmlname="UOT_parcelles", filename=paste0(INPUT.DIR, "/Train-Val/201908/val9k-s16trans/empty/UOT_parcelles.kml"))

writeOGR(val.plots, dsn=".", layer="UOT_parcelles", driver="ESRI Shapefile", overwrite=TRUE)
writeKML(val.plots, kmlname="UOT_parcelles", filename="UOT_parcelles.kml")


# Create sample grid
# ==================

grid.size <- 7
res <- res(landsat.grid)[1]
offset <- c(res/grid.size/2 + (0:(grid.size-1))*res/grid.size)

# split the plots for parallel processing
subsets <- split(val.plots, f=1:86)

val.grids <- foreach(subset=subsets, .combine=bind, .multicombine=TRUE) %dopar% { 
  grids <- SpatialPointsDataFrame(data.frame(x = 0, y = 0), data=data.frame(PLOTID = 0, SAMPLEID = 0, GRIDPOINT = 0))[-1,]
  for(p in 1:length(subset)) {
    plot <- subset[p,]
    ext <- extent(plot)
    grids <- bind(grids, SpatialPointsDataFrame(expand.grid(ext@xmin+offset, ext@ymin+offset), data=data.frame(PLOTID = plot$PLOTID, SAMPLEID = plot$SAMPLEID,GRIDPOINT = 1:grid.size^2)))
  }
  grids
}

proj4string(val.grids) <- proj4string(val.plots)

# merge with already collected tree attributes from train.grid
train.grids <- readOGR(paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/assessed_v1/COV_parcelles_grid.shp"))
val.grids <- merge(val.grids, train.grids@data[, c("PLOTID", "GRIDPOINT", "tree")], all.x=TRUE)


writeOGR(val.grids, dsn=".", layer="UOT_parcelles_grid", driver="ESRI Shapefile", overwrite=TRUE)


# Divide into 10 subsets and export
# =================================

subsets <- split(val.plots, f=1:10)
for(i in 1:length(subsets)) {
  writeOGR(subsets[[i]], dsn=".", layer=paste0("UOT_parcelles_", i), driver="ESRI Shapefile", overwrite=TRUE)
  writeKML(subsets[[i]], kmlname=paste0("UOT_parcelles_", i) , filename=paste0("UOT_parcelles_", i, ".kml"))
  subset.grids <- val.grids[val.grids$PLOTID %in% subsets[[i]]$PLOTID,]
  writeOGR(subset.grids, dsn=".", layer=paste0("UOT_parcelles_", i, "_grid"), driver="ESRI Shapefile", overwrite=TRUE)
}




# do the assessment
# =================




# read and merge assessed validation-plots
# ========================================

subsets <- dir(paste0(VALIDTN.DIR, "/1_val4k-s16trans/3_assessed"), pattern=".*_[[:digit:]]+\\.shp$", recursive=TRUE, full.names=TRUE)
val.plots.in <- lapply(subsets, readOGR)

names(val.plots.in) <- list.dirs(paste0(VALIDTN.DIR, "/1_val4k-s16trans/3_assessed"), recursive=FALSE)

for(i in 1:length(val.plots.in)) {
  val.plots.in[[i]]$author <- as.character(val.plots.in[[i]]$author)
  val.plots.in[[i]]$author[!is.na(val.plots.in[[i]]$author)] <- names(val.plots.in[i])
}

val.plots.in <- do.call(rbind, val.plots.in)

subsets.grids <- dir(paste0(VALIDTN.DIR, "/1_val4k-s16trans/3_assessed"), pattern=".*_[[:digit:]]+_grid\\.shp$", recursive=TRUE, full.names=TRUE)
val.grids.in <- do.call(bind, lapply(subsets.grids, readOGR))


# recalculate crown-cover and merge
val.plots.in <- val.plots.in[, names(val.plots.in) != "ccov"]
tmp <- aggregate(list(ccov=val.grids.in$tree), by=list(PLOTID=val.grids.in$PLOTID), FUN=function(x) sum(!is.na(x) & x==1)/sum(!is.na(x)))
tmp$ccov[is.nan(tmp$ccov)] <- NA
val.plots.in <- merge(val.plots.in, tmp, by="PLOTID")

# clean validation plots
#-----------------------

# remove plots without ccov
val.plots.in <- val.plots.in[!is.na(val.plots.in$ccov),]

# remove plots with strange dates
year <- as.numeric(sub("-.*$", "", as.character(val.plots.in$img_date)))
val.plots.in <- val.plots.in[!is.na(year) & year >= 2000 & year <= 2019, ]

# clean grid accordingly
val.grids.in <- val.grids.in[val.grids.in$PLOTID %in% unique(val.plots.in$PLOTID), ]


# Write clean files
# -----------------
writeOGR(val.plots.in, dsn=paste0(VALIDTN.DIR, "/1_val4k-s16trans/3_assessed"), layer="UOT_parcelles", driver="ESRI Shapefile", overwrite=TRUE)
writeOGR(val.grids.in, dsn=paste0(VALIDTN.DIR, "/1_val4k-s16trans/3_assessed"), layer="UOT_parcelles_grid", driver="ESRI Shapefile", overwrite=TRUE)


# Write some descriptive information
# ----------------------------------
sink(paste0(VALIDTN.DIR, "/1_val4k-s16trans/4_descr/author-stats.txt"), split=TRUE)
table(val.plots.in$author)
sink()





