###############################################################################
# 01_create-fc-maps.R: créer des cartes brutes du couvert forestier
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Définitions des variables ===================================================

COV.FC         <- 10             # Couverture des houppiers forêt/non-forêt

SAMPLE.DIST    <- 1              # Largeur de la lisière forêt à considérer
N.PIXELS       <- NA             # Number of non-NA cells (will be determined later)
SAMPLE.RATIO   <- 0.0025         # Share of non-NA cells to sample
CAL.RATIO      <- 0.75           # Use same amount of ref-points from cal.map as from ref.map / train.points
PREDICTORS     <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", "nbr", "ndmi", "ndvi", "evi", "BIO1", "BIO4", "BIO12", "BIO15") # "savi", "nbr2", "msavi", "x", "y" 


# Définitions des fonctions ===================================================

# Charger un image Landsat ----------------------------------------------------
#
# @param filename  Chemin du fichier landsat
#
# @return          Image Landsat avec bandes nommées 
#

load.image <- function(filename) {
  image <- brick(paste0(IMAGES.DIR, filename))
  names(image) <- SST.LSBANDS
  return(image)
}

# Tirer des points d'entraînement d'une carte autour de la lisière forêt ------
#
# @param map       Carte forêt/non-forêt à échantillonner
# @param n         Nombre d'échantillons à tirer
#
# @return          Points d'échantillon avec aatribut forêt/non-forêt
#

sample.map <- function(map, n) {
  
  tmp.src  <- tempfile(pattern = "", fileext = ".tif")             # temporary file for masked reference map
  tmp.dst1 <- tempfile(pattern = "", fileext = ".tif")             # temporary file for forest edge
  tmp.dst3 <- tempfile(pattern = "", fileext = ".tif")             # temporary file for
  
  map <- writeRaster(map, tmp.src)                               # write it to the disk
  
  # draw sample points around forest edge
  system(paste("gdal_proximity.py", tmp.src, tmp.dst1,             # forest and 3 pixel non-forest edge (distance to nearest pixel of value 1)
               "-values 1 -use_input_nodata YES -maxdist ",  SAMPLE.DIST, " -fixed-buf-val 3"))
  dst1 <- raster(tmp.dst1)
  NAvalue(dst1) <- 65535
  cat("     ")
  system(paste("gdal_proximity.py", tmp.src, tmp.dst3,             # non-forest and 3 pixel forest edge (distance to nearest pixel of value 3)
               "-values 3 -use_input_nodata YES -maxdist ",  SAMPLE.DIST, " -fixed-buf-val 1"))
  dst3 <- raster(tmp.dst3)
  NAvalue(dst3) <- 65535
  map <- mask(map, dst1)
  map <- mask(map, dst3)
  
  unlink(c(tmp.src, tmp.dst1, tmp.dst3))                           # delete temporary files
  
  n.classes <- length(unique(map))                                     # number of classes (should be two)
  cat(paste0("    -Sampling map (n=", n.classes, "*", round(n/n.classes), ") ... "))
  sample.pts <- sampleStratified(map, round(n/n.classes), sp=TRUE)[,-1]  # stratified sampling (same number of samples for each class)
  names(sample.pts) <- "CLASS"
  cat("done\n")
  return(sample.pts)
}

# Function for classifying an image ----------------------------------------

classify.image <- function(image, filename, bioclim=NULL, train.pts=NULL, ref.map=NULL, n.ref.map=NULL, cal.map=NULL, n.cal.map=NULL, mask=NULL, preds=NULL, type="classification", crossval=FALSE, prob=FALSE, n.cores=8) {
  
  txtfile <- paste0(sub("[.]tif$", "", filename), ".txt")
  cat("-- Image classification: ", basename(filename), "/", date(), " --\n", file=txtfile)
  
  if(!is.null(train.pts)) {
    cat("    -Loading training points ... ")
    train.pts <- train.pts[,1]                # only use first column in the attribute table (class)
    names(train.pts) <- "CLASS"
    set_ReplCRS_warn(FALSE)
    proj4string(train.pts) <- proj4string(image)
    cat("done\n")
    cat("Training points:", nrow(train.pts), "\n", file=txtfile, append=TRUE)
  }
  
  # add training data from ref.map if provided
  if(!is.null(ref.map)) {
    cat(paste0("    -Masking / buffering reference map ... \n"))
    ref.map  <- mask(crop(ref.map, image[[1]]), crop(image[[1]], ref.map)) # crop/mask ref.map with image
    
    if(!is.null(mask)) ref.map <- mask(ref.map, mask)                      # mask with additional mask, if provided
    
    if(!is.null(cal.map)) {
      tmp <- extend(crop(cal.map, ref.map), ref.map)                      # cut out the piece of the calibration map that overlaps ref map and extend to refmap
      ref.map <- mask(ref.map, tmp, inverse=TRUE)
    }
    
    cat("    ")
    ref.pts <- sample.map(ref.map, n.ref.map)
    cat("Ref-map points: ", nrow(ref.pts), "/", ref.map@file@name, "/", SAMPLE.DIST, "px\n", file=txtfile, append=TRUE)
    
    if(is.null(train.pts)) {
      train.pts <- ref.pts                                           # use it as training points or add to existing training points
    } else {
      train.pts <- rbind(train.pts, ref.pts)
    }
  }
  
  if(!is.null(cal.map)) {
    cat(paste0("    -Masking / buffering calibration map ... \n"))
    cal.map  <- mask(crop(cal.map, image[[1]]), crop(image[[1]], cal.map)) # crop/mask ref.map with image
    if(!is.null(mask)) cal.map <- mask(cal.map, mask)                      # mask with additional mask, if provided
    
    cat("    ")
    cal.pts <- sample.map(cal.map, n.cal.map)
    cat("Cal-map points: ", nrow(cal.pts), "from", cal.map@file@name, "/", SAMPLE.DIST, "px\n", file=txtfile, append=TRUE)
    
    if(is.null(train.pts)) {
      train.pts <- cal.pts                                           # use it as training points or add to existing training points
    } else {
      train.pts <- rbind(train.pts, cal.pts)
    }
  }
  
  cat("Total points:   ", nrow(train.pts), "\n", file=txtfile, append=TRUE)
  
  # extract spectral values
  if(is.null(preds)) {
    preds <- names(image)
    if(!is.null(bioclim)) preds <- c(preds, names(bioclim))
  }
  
  cat("    -Extracting pixel values for bands:", preds, "... ")
  train.pts <- raster::extract(image, train.pts, sp=TRUE)
  if(!is.null(bioclim)) train.pts <- raster::extract(bioclim, train.pts, sp=TRUE)
  train.dat <- na.omit(train.pts@data)[, c("CLASS", preds)]
  if(type=="classification") train.dat[,1] <- as.factor(train.dat[,1])
  cat("done\n")
  
  # calibrate RandomForest classifier 
  cat("    -Calibrating RandomForest ... ")
  sink(txtfile, append=TRUE)
  if(crossval) {
    map.model.cv <- train(y = train.dat[,1], 
                          x = train.dat[,-1],
                          method = "rf", 
                          importance = TRUE,
                          trControl = trainControl(
                            method = "repeatedcv",
                            number = 10,
                            repeats = 3))
    print(map.model.cv)
    map.model <- map.model.cv$finalModel
    print(map.model)
    cat("\n")
    print(varImp(map.model, scale=FALSE))
  } else {
    map.model <- randomForest(y=train.dat[,1], x=train.dat[,-1], importance=TRUE) # , do.trace=100)
    # Parallelization of RandomForest: confusion, err.rate, mse and rsq will be NULL
    # https://stackoverflow.com/questions/14106010/parallel-execution-of-random-forest-in-r
    # map.model <- foreach(ntree=rep(100, 5), .combine=randomForest::combine, .multicombine=TRUE, .packages='randomForest') %dopar% {
    #                 randomForest(x=ref.pts[,!(names(ref.pts) == "CLASS")], y=ref.pts$CLASS, importance=TRUE, ntree=ntree)
    #  
    print(map.model)
    cat("\n")
    print(varImp(map.model))
  }
  sink()
  if(type=="treecover") {
    cat("R2:", round(map.model$rsq[500], 2), "RMSE:", round(sqrt(map.model$mse[500]), 2), "\n")
  } else {
    cat("OOB error rate:", round(map.model$err.rate[500,1], 2), "\n")
  }
  
  # write model results
  dir.create(dirname(filename), recursive=TRUE, showWarnings=FALSE)
  # save RandomForest Model (too large)
  # save(map.model, file=paste0(sub("[.]tif$", "", filename), "r_rf.RData"))
  
  # classify image
  cat("    -Creating map ... ")
  if(!is.null(bioclim)) image <- stack(image, crop(bioclim, image))
  beginCluster(n=n.cores)
  map <- clusterR(image, predict, args=list(model=map.model))
  endCluster()
  
  # save map of classified image
  cat("writing map ... ")
  if(type=="treecover") map <- floor(map*100)
  map <- writeRaster(map, filename=filename, format="GTiff", datatype="INT2U", overwrite=TRUE)
  cat("done\n")
  
  # create probability map 
  if(prob==TRUE) {
    cat("    -Creating probability map ... ")
    beginCluster(n=n.cores)
    prob.map <- clusterR(image, predict, args=list(model=map.model, type="prob"))
    endCluster()
    cat("writing map ... ")
    writeRaster(prob.map, filename=sub("\\.tif", "_prob.tif", filename), format="GTiff", overwrite=TRUE)
    cat("done\n")
  } else {
    prob.map <- NULL
  }
  
  cat("-- Done: ", basename(filename), "/", date(), " --\n", file=txtfile, append=TRUE)
  
  invisible(list(
    "model"  = map.model,
    "map"    = map,
    "prob"   = prob.map
  ))
}


### DO THE WORK ###########################################################

# create directories if the don't exist
dir.create(FCC.REF.DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FCC.RAW.DIR, recursive=TRUE, showWarnings=FALSE)

# load 2018 images --------------------------------------------------------

ref.p192 <- brick(paste0(IMAGES.DIR, "/p192/p192_2018_m.tif"))
ref.p193 <- brick(paste0(IMAGES.DIR, "/p193/p193_2018_m.tif"))
ref.p194 <- brick(paste0(IMAGES.DIR, "/p194/p194_2018_m.tif"))
names(ref.p192) <- names(ref.p193) <- names(ref.p194) <- BANDS
ref.images <- list(p192=ref.p192, p193=ref.p193, p194=ref.p194)

N.PIXELS <- list(p192 = ncell(ref.p192[["B"]]) - summary(ref.p192)["NA's","B"],
                 p193 = ncell(ref.p193[["B"]]) - summary(ref.p193)["NA's","B"],
                 p194 = ncell(ref.p194[["B"]]) - summary(ref.p194)["NA's","B"])

bioclim.p192 <- brick(paste0(IMAGES.DIR, "/p192/p192_bioclim.tif"))
bioclim.p193 <- brick(paste0(IMAGES.DIR, "/p193/p193_bioclim.tif"))
bioclim.p194 <- brick(paste0(IMAGES.DIR, "/p194/p194_bioclim.tif"))
names(bioclim.p192) <- names(bioclim.p193) <- names(bioclim.p194) <- BIOCLIM
bioclim <- list(p192=bioclim.p192, p193=bioclim.p193, p194=bioclim.p194)


# load training plots ------------------------------------------------------

train.plots <- readOGR(paste0(TRNPTS.DIR, "/COV_parcelles.shp"))
train.plots <- train.plots[!is.na(train.plots$ccov), c("PLOTID", "ccov", "img_date", "author")]

train.plots$author <- as.factor(sub("^.*\\/\\/", "", train.plots$author))

# convert plot polygons to spatial points (centroids)
train.points <- SpatialPointsDataFrame(gCentroid(train.plots, byid=TRUE),  
                                    data.frame(author=train.plots$author, ccov=train.plots$ccov, img_date=as.Date(train.plots$img_date)))

# select only those with image date between 1.1.2017 and 31.12.2019
train.points <- train.points[!is.na(train.points$img_date) & train.points$img_date > as.Date("2017-01-01") & train.points$img_date <= as.Date("2019-12-31"), ]
pdf(paste0(FCC.REF.DIR, "/training-pts_2018.pdf"))
plot(train.points)
dev.off()

# extract image values for train.points
registerDoParallel(.env$numCores-1)
train.points <- foreach(i=1:length(ref.images), .combine=rbind) %dopar% {
  pts <- raster::extract(ref.images[[i]], train.points, sp=TRUE)
  pts <- raster::extract(bioclim[[i]], pts, sp=TRUE)
  pts$image <- names(ref.images[i])
  pts[, c("author", "image", "ccov", BANDS, BIOCLIM)]
}

# remove rows with NA's
train.points <- train.points[!is.na(rowSums(train.points@data[,-(1:2)])), ]

# discard points from authors that add confusion
# train.points <- train.points[!train.points$author %in% c("6_Eric_AGBESSI", "2_Mamalnassoh_ABIGUIME", "7_Aklasson_TOLEBA", "1_Ditorgue_BAKABIMA", "8_Yawo_KONKO"), ]

# for(author in unique(train.points$author[!train.points$author %in% discard])) {
#   print(author)
#   map.model <- randomForest(y=train.points@data[!train.points$author %in% c(author, discard),"ccov"], x=train.points@data[!train.points$author %in% c(author, discard),-(1:3)])
#   print(map.model)
# }



train.points@data <- cbind(train.points@data[,c("image", "ccov")],
                           F10=cut(train.points$ccov, breaks=c(0.0,0.1,1.0), labels=c(3, 1), right=FALSE, include.lowest=TRUE), 
                           F30=cut(train.points$ccov, breaks=c(0.0,0.3,1.0), labels=c(3, 1), right=FALSE, include.lowest=TRUE), 
                           train.points@data[,c(BANDS, BIOCLIM)])




# Parameter selection -----------------------------------------------------

cov.varsel <- rfe(y=train.points@data[train.points$image=="p193", "ccov"],
                  x=train.points@data[train.points$image=="p193", PREDICTORS],
                  sizes = c(4, 6, 8, 10),
                  rfeControl=rfeControl(
                    functions=rfFuncs,        # use RandomForest
                    method  = "repeatedcv",   # repeated cross-validation
                    number  = 10,             # 10-fold
                    repeats = 3))             # 3 repeats

print(cov.varsel)
predictors(cov.varsel)
plot(cov.varsel, type=c("g", "o"))

# 2018 reference tree cover map (TODO) -------------------------------------


# Tree cover map for 2018 for p193
set.seed(RSEED)
p193.2018.cov <- classify.image(image     = load.image("/p193/p193_2018_m.tif"), 
                                bioclim   = bioclim[["p193"]],
                                filename  = paste0(FCC.REF.DIR, "/p193_2018_COV_R.tif"),
                                train.pts = train.points[train.points$image == "p193", "ccov"],
                                preds     = PREDICTORS,
                                type      = "treecover",
                                crossval  = TRUE,
                                n.cores   = 32)

# plot observed vs. predicted
pdf(paste0(FCC.REF.DIR, "/p193_2018_COV_R.pdf"))
plot(p193.2018.cov[["model"]]$y, p193.2018.cov[["model"]]$predicted, xlim=c(0,1), ylim=c(0,1), main="Couverture houppier p193", xlab="Observation", ylab="Prédiction")
abline(0,1)
dev.off()

# 2018 reference forest cover maps -----------------------------------------

# Forest cover map for 2018 for p193
set.seed(RSEED)
classify.image(image     = load.image("/p193/p193_2018_m.tif"), 
               bioclim   = bioclim[["p193"]],
               filename  = paste0(FCC.REF.DIR, "/FC", COV.FC, "/p193_2018_FC", COV.FC, "_R.tif"),
               train.pts = train.points[train.points$image == "p193", paste0("F", COV.FC)],
               preds     = PREDICTORS,
               prob      = TRUE,
               crossval  = TRUE,
               n.cores   = 32)  

# Forest cover maps 2018 for p192 and p194, calibrating with p193
set.seed(RSEED)
registerDoParallel(.env$numCores-1)
foreach(path=c("p192", "p194")) %dopar% {
  
  train.pts <- train.points[train.points$image == path, paste0("F", COV.FC)]
  
  classify.image(image     = load.image(paste0("/", path, "/", path, "_2018_m.tif")), 
                 bioclim   = bioclim[[path]],
                 filename  = paste0(FCC.REF.DIR, "/FC", COV.FC, "/", path, "_2018_FC", COV.FC, "_R.tif"),
                 train.pts = train.pts,
                 cal.map   = raster(paste0(FCC.REF.DIR, "/FC", COV.FC, "/p193_2018_FC", COV.FC, "_R.tif")),
                 n.cal.map = max(2000, nrow(train.pts)/CAL.RATIO),
                 preds     = PREDICTORS,
                 mask      = TGO,
                 prob      = TRUE,
                 n.cores   = 32)
}

# 2003 reference forest cover maps -----------------------------------------

# Forest cover map for 2003 for p193
set.seed(RSEED)
classify.image(image     = load.image("/p193/p193_2003_m.tif"), 
               bioclim   = bioclim[["p193"]],
               filename  = paste0(FCC.REF.DIR, "/FC", COV.FC, "/p193_2003_FC", COV.FC, "_R.tif"),
               ref.map   = raster(paste0(FCC.REF.DIR, "/FC", COV.FC, "/p193_2018_FC", COV.FC, "_R.tif")),
               n.ref.map = 2 * SAMPLE.RATIO * N.PIXELS[["p193"]],
               preds     = PREDICTORS,
               mask      = TGO,
               n.cores   = 32)

# # Recalibrate forest cover map for 2018, based on 2003 for p193
#
# set.seed(RSEED)
# classify.image(image     = load.image("/p193/p193_2018.tif"), 
#                bioclim   = bioclim[["p193"]],
#                filename  = paste0(REFMAPS.DIR, "/p193_2018_FC", COV.FC, ".tif"),
#                ref.map   = raster(paste0(REFMAPS.DIR, "/p193_2003_FC", COV.FC, "_R.tif")),
#                n.ref.map = SAMPLE.RATIO * N.PIXELS[["p193"]],
#                preds     = PREDICTORS,
#                mask      = TGO,
#                n.cores   = 32)


# Forest cover maps 2003 for p192 and p194, calibrating with p193
set.seed(RSEED)
registerDoParallel(.env$numCores-1)
foreach(path=c("p192", "p194")) %dopar% {

  classify.image(image     = load.image(paste0("/", path, "/", path, "_2003_m.tif")), 
                 bioclim   = bioclim[[path]],
                 filename  = paste0(FCC.REF.DIR, "/FC", COV.FC, "/", path, "_2003_FC", COV.FC, "_R.tif"),
                 ref.map   = raster(paste0(FCC.REF.DIR, "/FC", COV.FC, "/", path, "_2018_FC", COV.FC, "_R.tif")),
                 n.ref.map = 2 * (1 - CAL.RATIO) * SAMPLE.RATIO * N.PIXELS[[path]],
                 cal.map   = raster(paste0(FCC.REF.DIR, "/FC", COV.FC, "/p193_2003_FC", COV.FC, "_R.tif")),
                 n.cal.map = 2 * CAL.RATIO * SAMPLE.RATIO * N.PIXELS[[path]],
                 preds     = PREDICTORS,
                 mask      = TGO,
                 n.cores   = 32)
  
}


# Forest cover maps for all dates ------------------------------------------

# Forest cover maps for p193
set.seed(RSEED)
registerDoParallel(.env$numCores-1)
foreach(file=dir(paste0(IMAGES.DIR, "/p193"), pattern="\\_[[:digit:]]+\\_m\\.tif")) %dopar% {
# foreach(file=c("p193_2017_m.tif", "p193_2019_m.tif")) %dopar% {     
  
  classify.image(image     = load.image(paste0("/p193/", file)), 
                 bioclim   = bioclim[["p193"]],
                 filename  = paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p193/", sub("\\_m\\.tif$", paste0("_F", COV.FC, "r.tif"), file)),
                 ref.map   = raster(paste0(FCC.REF.DIR, "/FC", COV.FC, "/p193_2003_FC", COV.FC, "_R.tif")),
                 n.ref.map = SAMPLE.RATIO * N.PIXELS[["p193"]],
                 preds     = PREDICTORS,
                 mask      = TGO,
                 n.cores   = 6)
                 # n.cores   = 32)  
}

# merge the two p193_1990 tiles
merge(raster(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p193/p193_1990_1_F", COV.FC, "r.tif")), 
      raster(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p193/p193_1990_2_F", COV.FC, "r.tif")),
      filename=paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p193/p193_1990_F", COV.FC, "r.tif"), format="GTiff", datatype="INT2U", overwrite=TRUE)

# Forest cover maps for p192 and p194 using p193 maps for calibration

set.seed(RSEED)
registerDoParallel(.env$numCores-1)
# foreach(file=c(dir(paste0(IMAGES.DIR, "/p192"), pattern="\\_[[:digit:]]+\\_m\\.tif"),
#                dir(paste0(IMAGES.DIR, "/p194"), pattern="\\_[[:digit:]]+\\_m\\.tif"))) %dopar% { 
foreach(file=c("p194_1997_m.tif", "p194_2007_m.tif", "p194_2015_m.tif", "p194_2018_m.tif")) %dopar% {

  
  path <- sub("\\_.*", "", file)
  
  if(file.exists(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p193/", sub("\\_m\\.tif$", paste0("_F", COV.FC, "r.tif"), sub(path, "p193", file))))) {
    cal.map   <- raster(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p193/", sub("\\_m\\.tif$", paste0("_F", COV.FC, "r.tif"), sub(path, "p193", file))))
    n.cal.map <- CAL.RATIO * SAMPLE.RATIO * N.PIXELS[[path]]
  } else {
    cal.map <- NULL
    n.cal.map <- NULL
  }
  
  classify.image(image     = load.image(paste0("/", path, "/", file)), 
                 bioclim   = bioclim[[path]],
                 filename  = paste0(FCC.RAW.DIR, "/FC", COV.FC, "/", path, "/", sub("\\_m\\.tif$", paste0("_F", COV.FC, "r.tif"), file)),
                 ref.map   = raster(paste0(FCC.REF.DIR, "/FC", COV.FC, "/", path, "_2003_FC", COV.FC, "_R.tif")),
                 n.ref.map = (1 - CAL.RATIO) * SAMPLE.RATIO * N.PIXELS[[path]],
                 cal.map   = cal.map,
                 n.cal.map = n.cal.map,
                 preds     = PREDICTORS,
                 mask      = TGO,
                 n.cores   = 16)
                 # n.cores   = 6)
} 


# Merging key date maps ---------------------------------------------------

for(year in VAL.YEARS) {
  merge(mask(crop(brick(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p193/p193_", year, "_F", COV.FC, "r.tif")),TGO), TGO),
        mask(crop(brick(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p192/p192_", year, "_F", COV.FC, "r.tif")),TGO), TGO),
        mask(crop(brick(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p194/p194_", year, "_F", COV.FC, "r.tif")),TGO), TGO),
        filename=paste0(FCC.RAW.DIR, "/FC", COV.FC, "/TGO/TGO_", year, "_F", COV.FC, "r.tif"), overwrite=TRUE)
}        

# and reference maps
for(map in c("2018_FC10_R", "2018_FC10_R_prob","2003_FC10_R")) {
# for(map in c("2018_FC30_R", "2018_FC30_R_prob","2003_FC30_R")) {
  merge(mask(crop(brick(paste0(FCC.REF.DIR, "/FC", COV.FC, "/p193_", map, ".tif")),TGO), TGO),
        mask(crop(brick(paste0(FCC.REF.DIR, "/FC", COV.FC, "/p192_", map, ".tif")),TGO), TGO),
        mask(crop(brick(paste0(FCC.REF.DIR, "/FC", COV.FC, "/p194_", map, ".tif")),TGO), TGO),
        filename=paste0(FCC.REF.DIR, "/FC", COV.FC, "/TGO_", map, ".tif"), overwrite=TRUE)
}   
                 