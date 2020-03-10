##################################################################################
# NERF_Togo/FCC/2_create-train-points.R: create a set of tree cover training plots
# --------------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 20 May 2019

library("spsurvey")

TRNPTS <- c(paste0(TRNPTS.DIR, "/1_trn10k-s10ndvi_rev/2_assessed"),
            paste0(TRNPTS.DIR, "/2_trn05k-prob040-060/2_assessed"),
            paste0(TRNPTS.DIR, "/3_trn1k-prob040-060/2_assessed"))

# Create a grid of observation points over whole Togo, based on Landsat images (30 x 30m) each 480 m -------------------

res <- 480

x.min <- res * extent(TGO)@xmin %/% res
x.max <- res * extent(TGO)@xmax %/% res
y.min <- res * extent(TGO)@ymin %/% res
y.max <- res * extent(TGO)@ymax %/% res


frame.points <- SpatialPoints(expand.grid(seq(x.min, x.max, by=res), seq(y.min, y.max, by=res)), proj4string=utm.31)[TGO]
# add attributes
frame.points$PLOTID  <- paste0(str_pad(frame.points@coords[,1], 7, "left", "0"), "_", str_pad(frame.points@coords[,2], 7, "left", "0"))
frame.points$xcoords <- frame.points@coords[,1]
frame.points$ycoords <- frame.points@coords[,2]

# load 2018 NDVI and mask with water, clouds and shadow 
ndvi.p192 <- mask(raster(paste0(OUTPUT.DIR, "/1_images/p192/p192_2018.tif"), band=12),
                  raster(paste0(OUTPUT.DIR, "/1_images/p192/p192_2018_qaLC08.tif")) %in% c(qa.cloud, qa.shadow, qa.water, qa.ice),
                  maskvalue=TRUE)
ndvi.p193 <- mask(raster(paste0(OUTPUT.DIR, "/1_images/p193/p193_2018.tif"), band=12),
                  raster(paste0(OUTPUT.DIR, "/1_images/p193/p193_2018_qaLC08.tif")) %in% c(qa.cloud, qa.shadow, qa.water, qa.ice),
                  maskvalue=TRUE)
ndvi.p194 <- mask(raster(paste0(OUTPUT.DIR, "/1_images/p194/p194_2018.tif"), band=12),
                  raster(paste0(OUTPUT.DIR, "/1_images/p194/p194_2018_qaLC08.tif")) %in% c(qa.cloud, qa.shadow, qa.water, qa.ice),
                  maskvalue=TRUE)
ndvi <- merge(ndvi.p192, ndvi.p193, ndvi.p194)

rm(ndvi.p192, ndvi.p193, ndvi.p194)

# read sampling frame and add NDVI for each plot
frame.points$ndvi   <- raster::extract(ndvi, frame.points)
frame.points$ndvi_c <- cut(frame.points$ndvi, 10, labels=paste0("s", 0:9))

writeOGR(frame.points, dsn=paste0(INPUT.DIR, "/Train-Val/201908"), layer="TGO_frame_480m", driver="ESRI Shapefile", overwrite=TRUE)



# Sampling frame for training-points -----------------------

# Initialize random number generator
set.seed(1)

# design for a spatially balanced sample, drawing 1500 samples of each stratum
Dsgn.grt <- list("s0"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s1"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s2"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s3"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s4"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s5"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s6"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s7"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s8"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s9"=list(panel=c(PanelOne=1500), seltype="Equal")
)


train.points <- grts(design=Dsgn.grt,         # using Dsgn design object
           DesignID='train',                  # prefix for each point name
           type.frame='finite',               # type
           src.frame='sp.object',             # sample frame is shapefile
           sp.object=frame.points,
           stratum="ndvi_c" 
)

summary(train.points)

# apply coordinate reference system 
proj4string(train.points) <- proj4string(frame.points)

# shuffle the rows
train.points <- train.points[sample(1:nrow(train.points)), ]
train.points$SAMPLEID <- paste0("trn-", str_pad(string=1:nrow(train.points), width = 4, pad = "0", side = "left"))

# Convert to plots and add attributes --------------------

landsat.grid <- raster(ndvi)
values(landsat.grid) <- 1

# convert to polygons
train.plots <- rasterToPolygons(mask(landsat.grid, train.points))

# fetch attributes
train.plots@data <- over(train.plots, train.points[, c("PLOTID", "SAMPLEID", "xcoords", "ycoords", "ndvi", "stratum")])

# add attributes
train.plots$ccov                            <- as.character(NA)
train.plots$img_src <- train.plots$img_date <- as.character(NA)
train.plots$author  <- train.plots$mod_date <- as.character(NA)

# # write plots as Shapefile and KML
writeOGR(train.plots, dsn=paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/empty"), layer="COV_parcelles", driver="ESRI Shapefile", overwrite=TRUE)
writeKML(train.plots, kmlname="COV_parcelles", filename=paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/empty/COV_parcelles.kml"))

# Create 7x7 sample grid ------------------------------------

grid.size <- 7
res <- res(landsat.grid)[1]
offset <- c(res/grid.size/2 + (0:(grid.size-1))*res/grid.size)

# split the plots for parallel processing
subsets <- split(train.plots, f=1:86)

train.grids <- foreach(subset=subsets, .combine=rbind, .multicombine=TRUE) %dopar% { 
  grids <- SpatialPointsDataFrame(data.frame(x = 0, y = 0), data=data.frame(PLOTID = 0, SAMPLEID = 0, GRIDPOINT = 0))[-1,]
  for(p in 1:length(subset)) {
    plot <- subset[p,]
    ext <- extent(plot)
    grids <- bind(grids, SpatialPointsDataFrame(expand.grid(ext@xmin+offset, ext@ymin+offset), data=data.frame(PLOTID = plot$PLOTID, SAMPLEID = plot$SAMPLEID, GRIDPOINT = 1:grid.size^2)))
  }
  grids
}

proj4string(train.grids) <- proj4string(train.plots)
train.grids$tree <- as.integer(NA) # 1 or 0

writeOGR(train.grids, dsn=paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/empty"), layer="COV_parcelles_grid", driver="ESRI Shapefile", overwrite=TRUE)


# Divide into 10 subsets and export --------------------

subsets <- split(train.plots, f=1:10)
for(i in 1:length(subsets)) {
  writeOGR(subsets[[i]], dsn=paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/empty"), layer=paste0("COV_parcelles_", i), driver="ESRI Shapefile", overwrite=TRUE)
  writeKML(subsets[[i]], kmlname=paste0("COV_parcelles_", i) , filename=paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/empty/COV_parcelles_", i, ".kml"))
  subset.grids <- train.grids[train.grids$PLOTID %in% subsets[[i]]$PLOTID,]
  writeOGR(subset.grids, dsn=paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/empty"), layer=paste0("COV_parcelles_", i, "_grid"), driver="ESRI Shapefile", overwrite=TRUE)
}




# read and merge assessed training-plots ------------------------

subsets <- dir(TRNPTS, pattern="[[:digit:]]\\.shp$", recursive=TRUE, full.names=TRUE)
train.plots.in <- lapply(subsets, readOGR)

names(train.plots.in) <- list.dirs(TRNPTS, recursive=FALSE)

for(i in 1:length(train.plots.in)) {
  train.plots.in[[i]]$author <- as.character(train.plots.in[[i]]$author)
  train.plots.in[[i]]$author[!is.na(train.plots.in[[i]]$author)] <- names(train.plots.in[i])
  train.plots.in[[i]] <- train.plots.in[[i]][, c("PLOTID", "SAMPLEID", "xcoords", "ycoords", "ccov", "img_date", "img_src", "mod_date", "author")]
  proj4string(train.plots.in[[i]]) <- utm.31
}

train.plots.in <- do.call(rbind, train.plots.in)



subsets.grids <- dir(TRNPTS, pattern=".*[[:digit:]]+_grid\\.shp$", recursive=TRUE, full.names=TRUE)
train.grids.in <- lapply(subsets.grids, readOGR)

for(i in 1:length(train.grids.in)) {
    proj4string(train.grids.in[[i]]) <- utm.31
}

train.grids.in <- do.call(rbind, train.grids.in)

# recalculate crown-cover and merge
train.plots.in <- train.plots.in[, names(train.plots.in) != "ccov"]
tmp <- aggregate(list(ccov=train.grids.in$tree), by=list(PLOTID=train.grids.in$PLOTID), FUN=function(x) sum(!is.na(x) & x==1)/sum(!is.na(x)))
tmp$ccov[is.nan(tmp$ccov)] <- NA
train.plots.in <- merge(train.plots.in, tmp, by="PLOTID")

# clean training plot data ---------------------

# remove plots without ccov
train.plots.in <- train.plots.in[!is.na(train.plots.in$ccov),]

# remove plots with strange dates
year <- as.numeric(sub("-.*$", "", as.character(train.plots.in$img_date)))
train.plots.in <- train.plots.in[!is.na(year) & year >= 2000 & year <= 2019, ]

# clean grid accordingly
train.grids.in <- train.grids.in[train.grids.in$PLOTID %in% unique(train.plots.in$PLOTID), ]

# Write clean training plots -----------------
writeOGR(train.plots.in, dsn=TRNPTS.DIR, layer="COV_parcelles", driver="ESRI Shapefile", overwrite=TRUE)
writeOGR(train.grids.in, dsn=TRNPTS.DIR, layer="COV_parcelles_grid", driver="ESRI Shapefile", overwrite=TRUE)

dir.create(paste0(TRNPTS.DIR, "/descr"), showWarnings = FALSE)

# Write some descriptive information
pdf(paste0(TRNPTS.DIR, "/descr/hist-years.pdf"))
hist(as.numeric(sub("-.*$", "", as.character(train.plots.in$img_date))), breaks=2000:2019)
dev.off()

pdf(paste0(TRNPTS.DIR, "/descr/hist-ccov.pdf"))
hist(train.plots.in$ccov)
dev.off()

sink(paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/descr/author-stats.txt"), split=TRUE)
table(train.plots.in$author)
sink()


# # Create additional 500 traininplots for regions of ambiquity (for Ayele) ------------------------
# 
# # get probability map
# prob.map <- raster(paste0(OUTPUT.DIR, "/2_forest-maps/TGO/1_ref-maps/p193/p193_2018_FC30_prob.tif"))
# prob.map[prob.map < 0.4 | prob.map > 0.6] <- NA
# prob.map <- mask(prob.map, TGO)
# 
# # sample cells with ambiguity (values 0.4 -  0.6)
# amb.plots <- rasterToPolygons(sampleRandom(prob.map, 500, asRaster=TRUE))
# 
# # fetch attributes
# names(amb.plots) <- "NDVI"
# amb.plots$PLOTID   <- paste0("add-", str_pad(string=1:nrow(amb.plots), width = 4, pad = "0", side = "left"))
# amb.plots$SAMPLEID <- amb.plots$PLOTID
# amb.plots$xcoords <- amb.plots$ycoords <- amb.plots$stratum <- as.numeric(NA)
# 
# # add attributes
# amb.plots$ccov                          <- as.character(NA)
# amb.plots$img_src <- amb.plots$img_date <- as.character(NA)
# amb.plots$author  <- amb.plots$mod_date <- as.character(NA)
# 
# # # write plots as Shapefile and KML
# writeOGR(amb.plots, dsn=paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/add"), layer="COV_parcelles_add", driver="ESRI Shapefile", overwrite=TRUE)
# writeKML(amb.plots, kmlname="COV_parcelles", filename=paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/add/COV_parcelles_add.kml"))
# 
# # Create sample grid
# 
# grid.size <- 7
# res <- res(amb.map)[1]
# offset <- c(res/grid.size/2 + (0:(grid.size-1))*res/grid.size)
# 
# # split the plots for parallel processing
# subsets <- split(amb.plots, f=1:10)
# 
# registerDoParallel(.env$numCores-1)
# amb.grids <- foreach(subset=subsets, .combine=rbind, .multicombine=TRUE) %dopar% { 
#   grids <- SpatialPointsDataFrame(data.frame(x = 0, y = 0), data=data.frame(PLOTID = 0, SAMPLEID = 0, GRIDPOINT = 0))[-1,]
#   for(p in 1:length(subset)) {
#     plot <- subset[p,]
#     ext <- extent(plot)
#     grids <- bind(grids, SpatialPointsDataFrame(expand.grid(ext@xmin+offset, ext@ymin+offset), data=data.frame(PLOTID = plot$PLOTID, SAMPLEID = plot$SAMPLEID, GRIDPOINT = 1:grid.size^2)))
#   }
#   grids
# }
# 
# proj4string(amb.grids) <- proj4string(amb.plots)
# amb.grids$tree <- as.integer(NA) # 1 or 0
# 
# writeOGR(amb.grids, dsn=paste0(INPUT.DIR, "/Train-Val/201908/trn10k-s10ndvi/add"), layer="COV_parcelles_add_grid", driver="ESRI Shapefile", overwrite=TRUE)


# # Create additional 1000 traininplots for regions of ambiquity (for Etse) ------------------------------
# 
# # get probability map
# prob.map <- merge(
#   raster(paste0(REFMAPS.DIR, "/p193_2018_FC30_R_prob.tif")),
#   raster(paste0(REFMAPS.DIR, "/p192_2018_FC30_R_prob.tif")),
#   raster(paste0(REFMAPS.DIR, "/p194_2018_FC30_R_prob.tif"))
# )
# 
# prob.map <- mask(crop(prob.map, TGO), TGO)
# prob.map[prob.map < 0.4 | prob.map > 0.6] <- NA
# 
# 
# # sample cells with ambiguity (values 0.4 -  0.6)
# amb.plots <- rasterToPolygons(sampleRandom(prob.map, 1000, asRaster=TRUE))
# 
# # fetch attributes
# names(amb.plots) <- "NDVI"
# amb.plots$PLOTID   <- paste0("add2-", str_pad(string=1:nrow(amb.plots), width = 4, pad = "0", side = "left"))
# amb.plots$SAMPLEID <- amb.plots$PLOTID
# amb.plots$xcoords <- amb.plots$ycoords <- amb.plots$stratum <- as.numeric(NA)
# 
# # add attributes
# amb.plots$ccov                          <- as.character(NA)
# amb.plots$img_src <- amb.plots$img_date <- as.character(NA)
# amb.plots$author  <- amb.plots$mod_date <- as.character(NA)
# 
# # # write plots as Shapefile and KML
# writeOGR(amb.plots, dsn=paste0(REFMAPS.DIR, "/3_trn1k-prob040-060/1_empty"), layer="COV_parcelles_add2", driver="ESRI Shapefile", overwrite=TRUE)
# writeKML(amb.plots, kmlname="COV_parcelles", filename=paste0(REFMAPS.DIR, "/3_trn1k-prob040-060/1_empty/COV_parcelles_add2.kml"))
# 
# # Create sample grid
# 
# grid.size <- 7
# res <- res(prob.map)[1]
# offset <- c(res/grid.size/2 + (0:(grid.size-1))*res/grid.size)
# 
# # split the plots for parallel processing
# subsets <- split(amb.plots, f=1:10)
# 
# registerDoParallel(.env$numCores-1)
# amb.grids <- foreach(subset=subsets, .combine=rbind, .multicombine=TRUE) %dopar% { 
#   grids <- SpatialPointsDataFrame(data.frame(x = 0, y = 0), data=data.frame(PLOTID = 0, SAMPLEID = 0, GRIDPOINT = 0))[-1,]
#   for(p in 1:length(subset)) {
#     plot <- subset[p,]
#     ext <- extent(plot)
#     grids <- bind(grids, SpatialPointsDataFrame(expand.grid(ext@xmin+offset, ext@ymin+offset), data=data.frame(PLOTID = plot$PLOTID, SAMPLEID = plot$SAMPLEID, GRIDPOINT = 1:grid.size^2)))
#   }
#   grids
# }
# 
# proj4string(amb.grids) <- proj4string(amb.plots)
# amb.grids$tree <- as.integer(NA) # 1 or 0
# 
# writeOGR(amb.grids, dsn=paste0(REFMAPS.DIR, "/3_trn1k-prob040-060/1_empty"), layer="COV_parcelles_add2_grid", driver="ESRI Shapefile", overwrite=TRUE)
















