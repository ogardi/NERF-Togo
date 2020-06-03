####################################################################
# NERF_Togo/AGB/3_create-agb-maps.R: create AGB maps for different dates
# ------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 10 October 2019

### DEFINITIONS ############################################################

# Default parameters ------------------------------------------------------- 

AGB.STRATA     <- 10             # Number of strata for sampling AGB from reference maps / calibration maps
N.PIXELS       <- NA             # Number of non-NA cells (will be determined later)
SAMPLE.RATIO   <- 0.001          # Share of non-NA cells to sample
CAL.RATIO      <- 1              # Use same amount of ref-points from cal.map as from ref.map / train.points
PREDICTORS  <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", "nbr", "ndmi", "ndvi", "evi", "BIO1", "BIO4", "BIO12", "BIO15") # "savi", "nbr2", "msavi", "x", "y" 
SEED           <- 20191114

# Function for loading an image --------------------------------------------

load.image <- function(filename) {
  image <- brick(paste0(IMAGES.DIR, filename))
  names(image) <- BANDS
  return(image)
}

# function for creating biomass map based on image and training data ------------

agb.map <- function(image, filename, bioclim=NULL, train.dat=NULL, ref.map=NULL, n.ref.map=NULL, cal.map=NULL, n.cal.map=NULL, mask=NULL, preds=NULL, crossval=FALSE, bias.corr=TRUE, n.cores=8) {
  
  name <- sub("[.]tif$", "", filename)
  
  txtfile <- paste0(sub("[.]tif$", "", filename), ".txt")
  cat("-- Biomass map: ", basename(filename), "/", date(), " --\n", file=txtfile)
  
  # add training data from ref.map if provided
  if(!is.null(ref.map)) {
    cat(paste0("    -Masking / buffering reference map ... \n"))
    ref.map  <- mask(crop(ref.map, image[[1]]), crop(image[[1]], ref.map))   # crop/mask ref.map with image
     
    if(!is.null(mask)) ref.map <- mask(ref.map, mask)                      # mask with additional mask, if provided
     
    if(!is.null(cal.map)) {
       tmp <- extend(crop(cal.map, ref.map), ref.map)                      # cut out the piece of the calibration map that overlaps ref map and extend to refmap
       ref.map <- mask(ref.map, tmp, inverse=TRUE)
    }
    
    cat("    ")
    
    cat(paste0("    -Sampling map (n=", AGB.STRATA, "*", round(n.ref.map/AGB.STRATA), ") ... "))
    ref.pts <- sampleStratified(cut(ref.map, AGB.STRATA), round(n.ref.map/AGB.STRATA), sp=TRUE)[,-1]  # stratified sampling (same number of samples for each class)
    names(ref.pts) <- "AGB"
    
    cat("extracting values ... \n")
    
    ref.dat <- cbind(AGB=raster::extract(ref.map, ref.pts, df=TRUE)[,-1], 
                     raster::extract(image, ref.pts, df=TRUE)[,-1],
                     raster::extract(bioclim, ref.pts, df=TRUE)[,-1])
    
    cat("Ref-map points: ", nrow(ref.dat), "/", ref.map@file@name, file=txtfile, append=TRUE)
    
    if(is.null(train.dat)) {
      train.dat <- ref.dat                                           # use it as training points or add to existing training points
    } else {
      train.dat <- rbind(train.dat, ref.dat)
    }
  }
  
  if(!is.null(cal.map)) {
    cat(paste0("    -Masking calibration map ... \n"))
    cal.map  <- mask(crop(cal.map, image[[1]]), crop(image[[1]], cal.map))   # crop/mask ref.map with image
    if(!is.null(mask)) cal.map <- mask(cal.map, mask)                        # mask with additional mask, if provided
    
    cat(paste0("    -Sampling map (n=", AGB.STRATA, "*", round(n.cal.map/AGB.STRATA), ") ... "))
    cal.pts <- sampleStratified(cut(cal.map, AGB.STRATA), round(n.cal.map/AGB.STRATA), sp=TRUE)[,-1]  # stratified sampling (same number of samples for each class)
    names(cal.pts) <- "AGB"
    
    cat("extracting values ... \n")
    
    cal.dat <- cbind(AGB=raster::extract(cal.map, cal.pts, df=TRUE)[,-1], 
                     raster::extract(image, cal.pts, df=TRUE)[,-1],
                     raster::extract(bioclim, cal.pts, df=TRUE)[,-1])
    
    if(is.null(train.dat)) {
        train.dat <- cal.dat                                           # use it as training points or add to existing training points
    } else {
        train.dat <- rbind(train.dat, cal.dat)
    }
      
    
    cat("Cal-map points: ", nrow(cal.dat), "from", cal.map@file@name, file=txtfile, append=TRUE)
    
  }
  
  cat("Total points:   ", nrow(train.dat), "\n", file=txtfile, append=TRUE)
  
  # extract spectral values
  if(is.null(preds)) {
    preds <- names(image)
    if(!is.null(bioclim)) preds <- c(preds, names(bioclim))
  }
  
  # cat("    -Extracting pixel values for bands:", preds, "... ")
  # train.pts <- raster::extract(image, train.pts, sp=TRUE)
  # if(!is.null(bioclim)) train.pts <- raster::extract(bioclim, train.pts, sp=TRUE)
  # train.dat <- na.omit(train.pts@data)[, c("CLASS", preds)]
  # if(type=="classification") train.dat[,1] <- as.factor(train.dat[,1])
  # cat("done\n")
  
  # calibrate RandomForest classifier 
  cat("    -Calibrating RandomForest ... ")
  sink(txtfile, append=TRUE)
  if(crossval) {
    map.model.cv <- train(y = train.dat[,1], 
                          x = train.dat[,preds],
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
    map.model <- randomForest(y=train.dat[,1], x=train.dat[,preds], importance=TRUE) # , do.trace=100)
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
  cat("R2:", round(map.model$rsq[500], 2), "RMSE:", round(sqrt(map.model$mse[500]), 2), "\n")
  
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
  
  pdf(paste0(name, ".pdf"))
  plot(train.dat$AGB ~ map.model$predicted, ylab="AGB (tDM/ha)", xlab="Predicted AGB (tDM/ha)")
  abline(0,1, lty=2)
  if(bias.corr) {
    bc <- lm(train.dat$AGB ~ map.model$predicted)
    abline(bc, lty=3, col="red")
  } else {
    bc <- NULL
  }
  dev.off()
  
  if(bias.corr) {
    sink(paste0(name, "_bc.txt"), split=TRUE)
    print(summary(bc))
    sink()
    
    save(bc, file=paste0(name, "_bc.RData"))
    
    cat("Applying linear bias correction ...")
    map <- calc(map, fun=function(x){bc$coefficients[1] + bc$coefficients[2]*x})
    
  } 
  
  # save map of classified image
  cat("writing map ... ")
  map <- writeRaster(map, filename=filename, format="GTiff", datatype="INT2U", overwrite=TRUE)
  cat("done\n")
  
  
  cat("-- Done: ", basename(filename), "/", date(), " --\n", file=txtfile, append=TRUE)
  
  invisible(list(
    "rf.model" = map.model,
    "bc.model" = bc,
    "map"      = map
  ))
}
  


### DO THE WORK ############################################################
# read images 2015 and bioclim variables  ---------------------------

ref.p192 <- brick(paste0(IMAGES.DIR, "/p192/p192_2015_m.tif"))
ref.p193 <- brick(paste0(IMAGES.DIR, "/p193/p193_2015_m.tif"))
ref.p194 <- brick(paste0(IMAGES.DIR, "/p194/p194_2015_m.tif"))
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


# # Load 30m DEM and calculate slope and aspect
# dem     <- stack(raster(paste0(OUTPUT.DIR, "/1_images/SRTM/SRTM-1arcsec.tif")))
# names(dem) <- "ALT"
# dem$SLP <- terrain(dem$ALT, opt="slope")
# dem$ASP <- terrain(dem$ALT, opt="aspect")
# 
# 
# wc2 <- stack(paste0(OUTPUT.DIR, "/1_images/WCv2/wc2.0_30s_tmin_Togo.tif"),
#              paste0(OUTPUT.DIR, "/1_images/WCv2/wc2.0_30s_tmax_Togo.tif"),
#              paste0(OUTPUT.DIR, "/1_images/WCv2/wc2.0_30s_prec_Togo.tif"),
#              paste0(OUTPUT.DIR, "/1_images/WCv2/wc2.0_30s_bio_Togo.tif"))
# 
# names(wc2) <- c(paste0("tmin.", 1:12), paste0("tmax.", 1:12), paste0("prec.", 1:12), paste0("bio.", 1:19))


# load inventory data ------------------------------------------------------

plots      <- read.csv(paste0(AGB.REF.DIR, "/IFN-plots.csv")) # , fileEncoding="macintosh")
coordinates(plots) <- ~X+Y 
proj4string(plots) <- utm.31

pdf(paste0(AGB.REF.DIR, "/IFN-plots_location.pdf"),
    width=3.5, height=7)
par(mar=c(1,1,1,1))
plot(spTransform(TGO, utm.31), col="lightyellow")
plot(plots, add=TRUE, col="black", pch=16, cex=0.3)
plot(plots, add=TRUE, col="darkgreen", pch=1, cex=plots$AGB/100)
dev.off()

# convert points to polygons, for extraction of raster values
plots.poly <- SpatialPolygonsDataFrame(gBuffer(plots, byid=TRUE, width=20), plots@data)

# extract raster values
registerDoParallel(.env$numCores-1)
x <- foreach(i=1:length(ref.images), .combine=cbind) %:% 
  foreach(j=1:nlayers(ref.images[[i]]), .combine=cbind) %dopar% {
      raster::extract(ref.images[[i]][[j]], plots.poly, weights=TRUE, fun=mean, df=TRUE)[,2]  # weighted means for spectral values
  }


x2 <- foreach(i=1:length(bioclim), .combine=cbind) %:% 
  foreach(j=1:nlayers(bioclim[[i]]), .combine=cbind) %dopar% {
    raster::extract(bioclim[[i]][[j]], plots, df=TRUE)[,2]                                    # bioclimatic values
  }

dat.p192 <- na.omit(as.data.frame(cbind(plots$AGB, x[,1:13], x2[,1:19])))
dat.p193 <- na.omit(as.data.frame(cbind(plots$AGB, x[,14:26], x2[,20:38])))
dat.p194 <- na.omit(as.data.frame(cbind(plots$AGB, x[,27:39], x2[,39:57])))
names(dat.p192) <- names(dat.p193) <- names(dat.p194) <- c("AGB", BANDS, BIOCLIM)
train.data <- list(p192=dat.p192, p193=dat.p193, p194=dat.p194)


# variable selection -----------------------------

# train.data <- na.omit(x[,c(1, (2*13-13)+2:14, 41:98)])       # use '2' for selecting p193
# 
# # determine explicative variables
# 
# # TODO: include regeneration map as predictor!!
# 
# var.sel <- rfe(train.data[,2:ncol(train.data)], train.data[,1], 
#                sizes=c(1:15, 20, 50), #sizes=seq(8,26,2), #
#                rfeControl=rfeControl(
#                  functions=rfFuncs,
#                  method = "repeatedcv",
#                  number = 10,
#                  repeats = 10
#                ))
# 
# sink(paste0(OUTPUT.DIR, "/3_forest-biomass/2_ref-maps/p193_2015_AGB_varsel-fin.txt"), split=TRUE)
# predictors(var.sel) 
# print(var.sel)
# sink()
# 
# pdf(paste0(OUTPUT.DIR, "/3_forest-biomass/2_ref-maps/p193_2015_AGB_varsel-fin.pdf"))
# plot(var.sel, type=c("g", "o"))
# dev.off()
# 
# # train.data <- train.data[, c("AGB", predictors(var.sel))]
# 
# # use Landsat bands + ndvi + annual mean temp, temp seasonality, annual prec and prec seasonality
# train.data <- train.data[, c("AGB", 
#                              "B", "G", "R", "NIR", "SWIR1", "SWIR2", "ndvi",
#                              paste0("bio.", c(1, 4, 12, 15))
# )]



# create reference biomass maps 2015 ------------------------------

set.seed(SEED)
p193.2015.agb <- agb.map(image     = load.image("/p193/p193_2015_m.tif"), 
                         bioclim   = bioclim[["p193"]],
                         filename  = paste0(AGB.REF.DIR, "/p193_2015_AGB_R.tif"),
                         train.dat = train.data[["p193"]],
                         preds     = PREDICTORS,
                         crossval  = TRUE,
                         bias.corr = FALSE,
                         n.cores   = 32)

set.seed(SEED)
p192.2015.agb <- agb.map(image     = load.image("/p192/p192_2015_m.tif"), 
                         bioclim   = bioclim[["p192"]],
                         filename  = paste0(AGB.REF.DIR, "/p192_2015_AGB_R.tif"),
                         train.dat = train.data[["p192"]],
                         cal.map   = raster(paste0(AGB.REF.DIR, "/p193_2015_AGB_R.tif")),
                         n.cal.map = max(200, nrow(train.data[["p192"]])*CAL.RATIO),
                         preds     = PREDICTORS,
                         crossval  = FALSE,
                         bias.corr = FALSE,
                         n.cores   = 32,
                         mask      = TGO)

set.seed(SEED)
p194.2015.agb <- agb.map(image     = load.image("/p194/p194_2015_m.tif"), 
                         bioclim   = bioclim[["p194"]],
                         filename  = paste0(AGB.REF.DIR, "/p194_2015_AGB_R.tif"),
                         train.dat = train.data[["p194"]],
                         cal.map   = raster(paste0(AGB.REF.DIR, "/p193_2015_AGB_R.tif")),
                         n.cal.map = max(200, nrow(train.data[["p194"]])*CAL.RATIO),
                         preds     = PREDICTORS,
                         crossval  = FALSE,
                         bias.corr = FALSE,
                         n.cores   = 32,
                         mask      = TGO) 

# merge the three maps
agb.2015 <- mask(crop(mosaic(raster(paste0(AGB.REF.DIR, "/p192_2015_AGB_R.tif")),
                             raster(paste0(AGB.REF.DIR, "/p193_2015_AGB_R.tif")), 
                             raster(paste0(AGB.REF.DIR, "/p194_2015_AGB_R.tif")),
                             fun=mean),
                      TGO), 
                 TGO, 
                 filename=paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.tif"), overwrite=TRUE)

# plot biomass map
library(RColorBrewer)
pdf(paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.pdf"),
    width=3.5, height=7)
par(mar=c(1,1,1,1))
plot(agb.2015, axes=FALSE, col=brewer.pal(9, "YlGn"), zlim=c(0,220))
plot(spTransform(TGO, utm.31), add=TRUE)
dev.off()

# compare resulting biomass map with the IFN AGB values
agb.pred <- raster::extract(agb.2015, plots.poly, weights=TRUE, fun=mean, df=TRUE)[,2]  # weighted means for spectral values
pdf(paste0(AGB.REF.DIR, "/AGB-model_2015.pdf"))
plot(agb.pred ~ plots.poly$AGB, xlab="Biomasse aérienne IFN (t/ha)", ylab="Carte AGB 2015 (t/ha)", xlim=c(0,350), ylim=c(0,350))
abline(0, 1, lty="dashed")
dev.off()

plots.poly$AGB.pred <- agb.pred

# bias correction
train.control <- trainControl(method = "cv", number = 10)
cv <- train(AGB ~ AGB.pred, data=plots.poly@data[!is.na(plots.poly$AGB.pred),], method = "lm", trControl = train.control)
model <- lm(AGB ~ AGB.pred, data=plots.poly@data[!is.na(plots.poly$AGB.pred),])
sink(paste0(AGB.REF.DIR, "/TGO_2015_AGB_Rc.txt"), split=TRUE)
print(cv)
summary(model)
sink()
agb.2015 <- writeRaster(model$coefficients["(Intercept)"] + model$coefficients["AGB.pred"] * agb.2015, paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.tif"), overwrite=TRUE)

# compare bias corrected biomass map with the IFN AGB values
agb.pred.c <- raster::extract(agb.2015, plots.poly, weights=TRUE, fun=mean, df=TRUE)[,2]  # weighted means for spectral values
pdf(paste0(AGB.REF.DIR, "/AGB-model_2015_c.pdf"))
plot(agb.pred ~ plots.poly$AGB, xlab="Biomasse aérienne IFN (t/ha)", ylab="Carte AGB 2015 (t/ha)", xlim=c(0,350), ylim=c(0,350))
abline(0, 1, lty="dashed") 
dev.off()


# create reference biomass maps 2003 ------------------------------

set.seed(SEED)
p193.2003.agb <- agb.map(image     = load.image("/p193/p193_2003_m.tif"), 
                         bioclim   = bioclim[["p193"]],
                         filename  = paste0(AGB.REF.DIR, "/p193_2003_AGB_R.tif"),
                         train.dat = NULL,
                         ref.map   = raster(paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.tif")),
                         n.ref.map = SAMPLE.RATIO * N.PIXELS[["p193"]],
                         preds     = PREDICTORS,
                         crossval  = FALSE,
                         bias.corr = TRUE,
                         n.cores   = 32)

set.seed(SEED)
p192.2003.agb <- agb.map(image     = load.image("/p192/p192_2003_m.tif"), 
                         bioclim   = bioclim[["p192"]],
                         filename  = paste0(AGB.REF.DIR, "/p192_2003_AGB_R.tif"),
                         train.dat = NULL,
                         ref.map   = raster(paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.tif")),
                         n.ref.map = SAMPLE.RATIO * N.PIXELS[["p192"]] / (1 + CAL.RATIO),
                         cal.map   = raster(paste0(AGB.REF.DIR, "/p193_2003_AGB_R.tif")),
                         n.cal.map = SAMPLE.RATIO * N.PIXELS[["p192"]] / (1 + 1/CAL.RATIO),
                         preds     = PREDICTORS,
                         crossval  = FALSE,
                         bias.corr = FALSE,
                         n.cores   = 32,
                         mask      = TGO)

set.seed(SEED)
p194.2003.agb <- agb.map(image     = load.image("/p194/p194_2003_m.tif"), 
                         bioclim   = bioclim[["p194"]],
                         filename  = paste0(AGB.REF.DIR, "/p194_2003_AGB_R.tif"),
                         train.dat = NULL,
                         ref.map   = raster(paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.tif")),
                         n.ref.map = SAMPLE.RATIO * N.PIXELS[["p194"]] / (1 + CAL.RATIO),
                         cal.map   = raster(paste0(AGB.REF.DIR, "/p193_2003_AGB_R.tif")),
                         n.cal.map = SAMPLE.RATIO * N.PIXELS[["p194"]] / (1 + 1/CAL.RATIO), 
                         preds     = PREDICTORS,
                         crossval  = FALSE,
                         bias.corr = FALSE,
                         n.cores   = 32,
                         mask      = TGO)

# merge the three maps
agb.2003 <- mask(crop(mosaic(raster(paste0(AGB.REF.DIR, "/p192_2003_AGB_R.tif")),
                             raster(paste0(AGB.REF.DIR, "/p193_2003_AGB_R.tif")), 
                             raster(paste0(AGB.REF.DIR, "/p194_2003_AGB_R.tif")),
                             fun=mean),
                      TGO), 
                 TGO, 
                 filename=paste0(AGB.REF.DIR, "/TGO_2003_AGB_R.tif"), overwrite=TRUE)


# create reference biomass maps 2018 ------------------------------

set.seed(SEED)
p193.2018.agb <- agb.map(image     = load.image("/p193/p193_2018_m.tif"), 
                         bioclim   = bioclim[["p193"]],
                         filename  = paste0(AGB.REF.DIR, "/p193_2018_AGB_R.tif"),
                         train.dat = NULL,
                         ref.map   = raster(paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.tif")),
                         n.ref.map = SAMPLE.RATIO * N.PIXELS[["p193"]],
                         preds     = PREDICTORS,
                         crossval  = FALSE,
                         bias.corr = TRUE,
                         n.cores   = 32)

set.seed(SEED)
p192.2018.agb <- agb.map(image     = load.image("/p192/p192_2018_m.tif"), 
                         bioclim   = bioclim[["p192"]],
                         filename  = paste0(AGB.REF.DIR, "/p192_2018_AGB_R.tif"),
                         train.dat = NULL,
                         ref.map   = raster(paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.tif")),
                         n.ref.map = SAMPLE.RATIO * N.PIXELS[["p192"]] / (1 + CAL.RATIO),
                         cal.map   = raster(paste0(AGB.REF.DIR, "/p193_2018_AGB_R.tif")),
                         n.cal.map = SAMPLE.RATIO * N.PIXELS[["p192"]] / (1 + 1/CAL.RATIO),
                         preds     = PREDICTORS,
                         crossval  = FALSE,
                         bias.corr = FALSE,
                         n.cores   = 32,
                         mask      = TGO)

set.seed(SEED)
p194.2018.agb <- agb.map(image     = load.image("/p194/p194_2018_m.tif"), 
                         bioclim   = bioclim[["p194"]],
                         filename  = paste0(AGB.REF.DIR, "/p194_2018_AGB_R.tif"),
                         train.dat = NULL,
                         ref.map   = raster(paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.tif")),
                         n.ref.map = SAMPLE.RATIO * N.PIXELS[["p194"]] / (1 + CAL.RATIO),
                         cal.map   = raster(paste0(AGB.REF.DIR, "/p193_2018_AGB_R.tif")),
                         n.cal.map = SAMPLE.RATIO * N.PIXELS[["p194"]] / (1 + 1/CAL.RATIO), 
                         preds     = PREDICTORS,
                         crossval  = FALSE,
                         bias.corr = FALSE,
                         n.cores   = 32,
                         mask      = TGO)

# merge the three maps
agb.2018 <- mask(crop(mosaic(raster(paste0(AGB.REF.DIR, "/p192_2018_AGB_R.tif")),
                             raster(paste0(AGB.REF.DIR, "/p193_2018_AGB_R.tif")), 
                             raster(paste0(AGB.REF.DIR, "/p194_2018_AGB_R.tif")),
                             fun=mean),
                      TGO), 
                 TGO, 
                 filename=paste0(AGB.REF.DIR, "/TGO_2018_AGB_R.tif"), overwrite=TRUE)


# # biomass maps for p193 --------------------------------------
# 
# registerDoParallel(.env$numCores-1)
# foreach(file=dir(paste0(IMAGES.DIR, "/p193"), pattern="\\_[[:digit:]]+\\_m\\.tif")) %dopar% {
#   #foreach(file=c("p193_1987.tif", "p193_2003.tif", "p193_2015.tif", "p193_2018.tif")) %dopar% {     
#   # foreach(file=c("p193_1985.tif", "p193_1990_1.tif", "p193_1990_2.tif", "p193_2000.tif", "p193_2005.tif", "p193_2007.tif", "p193_2009.tif", "p193_2013.tif", "p193_2017.tif", "p193_2019.tif")) %dopar% {  
#   
#   agb.map(image     = load.image(paste0("/p193/", file)), 
#           bioclim   = bioclim[["p193"]],
#           filename  = paste0(BIOMASS.DIR, "/3_raw-maps/p193/", sub("\\.tif$", "r.tif", file)),
#           ref.map   = raster(paste0(BIOMASS.DIR, "/2_ref-maps/p193_2003_AGB_R.tif")),
#           n.ref.map = SAMPLE.RATIO * N.PIXELS[["p193"]],
#           preds     = PREDICTORS,
#           mask      = TGO,
#           n.cores   = 6)
#   # n.cores   = 32) 
# } 
# 
# # merge the two p193_1990 tiles
# merge(raster(paste0(BIOMASS.DIR, "/3_raw-maps/p193/p193_1990_1_mr.tif")), 
#       raster(paste0(BIOMASS.DIR, "/3_raw-maps/p193/p193_1990_2_mr.tif")),
#       filename=paste0(BIOMASS.DIR, "/3_raw-maps/p193/p193_1990_mr.tif"), format="GTiff", datatype="INT2U", overwrite=TRUE)
# 
# 
# # biomass maps for p192 and p194 ----------------------
# 
# registerDoParallel(.env$numCores-1)
# foreach(file=c(dir(paste0(IMAGES.DIR, "/p192"), pattern="\\_[[:digit:]]+\\_m\\.tif"),
#                dir(paste0(IMAGES.DIR, "/p194"), pattern="\\_[[:digit:]]+\\_m\\.tif"))) %dopar% {
#                  
#   path <- sub("\\_.*", "", file)
#                  
#   if(file.exists(paste0(BIOMASS.DIR, "/3_raw-maps/p193/", sub("\\.tif$", "r.tif", sub(path, "p193", file))))) {
#     cal.map   <- raster(paste0(BIOMASS.DIR, "/3_raw-maps/p193/", sub("\\.tif$", "r.tif", sub(path, "p193", file))))
#     n.cal.map <- SAMPLE.RATIO * N.PIXELS[[path]] / (1 + 1/CAL.RATIO)
#     n.ref.map <- SAMPLE.RATIO * N.PIXELS[[path]] / (1 + CAL.RATIO)
#   } else {
#     cal.map <- NULL
#     n.cal.map <- NULL
#     n.ref.map <- SAMPLE.RATIO * N.PIXELS[[path]]
#   }
#   
#   agb.map(image     = load.image(paste0("/", path, "/", file)), 
#           bioclim   = bioclim[[path]],
#           filename  = paste0(BIOMASS.DIR, "/3_raw-maps/", path, "/", sub("\\.tif$", "r.tif", file)),
#           ref.map   = raster(paste0(BIOMASS.DIR, "/2_ref-maps/", path, "_2003_AGB_R.tif")),
#           n.ref.map = n.ref.map, 
#           cal.map   = cal.map,
#           n.cal.map = n.cal.map,
#           preds     = PREDICTORS,
#           mask      = TGO,
#           n.cores   = 32)
# } 
