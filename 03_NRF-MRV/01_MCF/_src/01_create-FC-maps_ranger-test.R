###############################################################################
# 01_create-fc-maps.R: créer des cartes brutes du couvert forestier
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Définitions des variables ===================================================

source("init.R")

# Seuil de la couverture des houppiers (forêt vs. non-forêt)
COV.FC         <- 30             

# Nombre de pixels à considerer comme lisière forêt 
SAMPLE.DIST    <- 1              

# Nombre de pixels non-NA (sera déterminé plus tard)
N.PIXELS       <- NA        

# Part de pixels à prendre en compte pour la calibration des cartes
SAMPLE.RATIO   <- 0.0025        

# Part des pixels à prendre des autres chemins (calibration)
CAL.RATIO      <- 0.75 

# Bandes à utiliser pour la modélisation forêt vs. non-forêt
PREDICTORS     <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", 
                    "nbr", "ndmi", "ndvi", "evi", 
                    "BIO1", "BIO4", "BIO12", "BIO15") 

# Répertoires
LANDSAT.DIR <- DIR.SST.DAT.LST
WORLDCLIM.DIR <- DIR.SST.DAT.WC2
REF.DIR <- DIR.MRV.MCF.REF
RAW.DIR <- DIR.MRV.MCF.RAW
  

# Définitions des fonctions ===================================================

# Charger un image Landsat ----------------------------------------------------
#
# @param filename  Chemin du fichier landsat
#
# @return          Image Landsat avec bandes nommées 
#

load.image <- function(filename) {
  image <- brick(paste0(LANDSAT.DIR, filename))
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
  
  tmp.src  <- tempfile(pattern = "", fileext = ".tif")    # tmp carte forêt/non-forêt
  tmp.dst1 <- tempfile(pattern = "", fileext = ".tif")    # tmp lisière à l'extérieur
  tmp.dst3 <- tempfile(pattern = "", fileext = ".tif")    # tmp lisière à l'intérieur
  
  # Écrire la carte sur le disque
  map <- writeRaster(map, tmp.src)                        
  
  # Forêt + lisière à l'extérieur de la forêt
  system(paste("gdal_proximity.py", tmp.src, tmp.dst1,    
               "-values", FOREST, 
               "-use_input_nodata YES", 
               "-maxdist ",  SAMPLE.DIST, 
               "-fixed-buf-val", NONFOR))
  dst1 <- raster(tmp.dst1)
  NAvalue(dst1) <- 65535
  cat("     ")
  
  # Non-forêt + lisière à l'intérieur de la forêt
  system(paste("gdal_proximity.py", tmp.src, tmp.dst3,          
               "-values", NONFOR, 
               "-use_input_nodata YES",
               "-maxdist ",  SAMPLE.DIST, 
               "-fixed-buf-val", FOREST))
  dst3 <- raster(tmp.dst3)
  NAvalue(dst3) <- 65535
  
  # Masquer la carte avec les lisières
  map <- mask(map, dst1)
  map <- mask(map, dst3)
  
  # Supprimer les fichiers temporaires
  unlink(c(tmp.src, tmp.dst1, tmp.dst3))                           
  
  # Echantillonnage stratifié des lisières forêt et non-forêt
  n.classes <- length(unique(map))
  cat(paste0("    -Sampling map (n=", n.classes, "*", round(n/n.classes), ") ... "))
  sample.pts <- sampleStratified(map, round(n/n.classes), sp=TRUE)[,-1]
  names(sample.pts) <- "CLASS"
  cat("done\n")
  return(sample.pts)
}
  

# Classification d'une image -------------------------------------------
#
# @param image      Image Landsat à classifier      
# @param filename   Nom de fichier pour la sauvegarde de la carte
# @param bioclim    Raster des variables bioclimatiques à utiliser
# @param train.pts  Points d'entraînement
# @param ref.map    Carte de référence
# @param n.ref.map  Nombre de points à échantilloner
# @param cal.map    Carte de calibration (autre chemin WRS)
# @param n.cal.map  Nombre de points à échantilloner
# @param mask       Masque à utiliser pour ref.map et cal.map
# @param preds      Variables à utiliser pour le modèle
# @param type       Modèle de classification ou de regression
# @param crossval   Faire validation croisée (3 * 10-fold)
# @param prob       Produire également carte de probabilité
# @param n.cores    Nombre de processeurs à utiliser pour prédiction
#
# @return           List avec les éléments 
#                   - model, 
#                   - carte forêt/non-forêt
#                   - carte des probabilités (si disponible)
#

classify.image <- function(image, filename, bioclim=NULL, train.pts=NULL, 
                           ref.map=NULL, n.ref.map=NULL, 
                           cal.map=NULL, n.cal.map=NULL, 
                           mask=NULL, preds=NULL, type="classification", 
                           n.cores=CORES-1) {
  
  # Ouvrir le fichier journal
  txtfile <- paste0(sub("[.]tif$", "", filename), ".txt")
  cat("-- Image classification: ", basename(filename), "/", 
      date(), " --\n", file=txtfile)
  
  # Empiler les couches Landsat et bioclim
  if(!is.null(bioclim)) image <- raster::stack(image, crop(bioclim, image))
  # et utiliser toutes les variables si non-spécifiées dans les paramètres
  if(!is.null(preds)){
    image <- image[[preds]]
  } else {
    preds <- names(image)
  }
  
  # Charger des points d'entraînement ---------------------
  if(!is.null(train.pts)) {
    cat("    -Loading training points ... ")
    train.pts <- train.pts[,1]      # utiliser que la première colonne ...
    names(train.pts) <- "CLASS"     # ... et nommer "CLASS"
    set_ReplCRS_warn(FALSE)
    proj4string(train.pts) <- proj4string(image)  # Système de coordonnées CRS
    cat("done\n")
    cat("Training points:", nrow(train.pts), "\n", file=txtfile, append=TRUE)
  }
  
  # Ajouter des points d'une carte de référence -----------
  if(!is.null(ref.map)) {
    cat(paste0("    -Masking / buffering reference map ... \n"))
    # ... couper/masquer avec l'image
    ref.map  <- mask(crop(ref.map, image[[1]]), crop(image[[1]], ref.map)) 
    # ... et masque additionelle (si disponible)
    if(!is.null(mask)) ref.map <- mask(ref.map, mask)             
    # Découper la carte de calibration (si disponible)
    if(!is.null(cal.map)) {
      tmp <- extend(crop(cal.map, ref.map), ref.map)
      ref.map <- mask(ref.map, tmp, inverse=TRUE)
    }
    cat("    ")
    # Tirer des points d'échantillon ...
    ref.pts <- sample.map(ref.map, n.ref.map)
    cat("Ref-map points: ", nrow(ref.pts), "/", ref.map@file@name, "/", 
        SAMPLE.DIST, "px\n", file=txtfile, append=TRUE)
    # ... et ajouter aux points d'entraînement
    if(is.null(train.pts)) {
      train.pts <- ref.pts                                          
    } else {
      train.pts <- rbind(train.pts, ref.pts)
    }
  }
  
  # Ajouter des points d'une carte de calibration ---------
  if(!is.null(cal.map)) {
    cat(paste0("    -Masking / buffering calibration map ... \n"))
    # ... couper/masquer avec l'image
    cal.map  <- mask(crop(cal.map, image[[1]]), crop(image[[1]], cal.map)) 
    # ... et masque additionelle (si disponible)
    if(!is.null(mask)) cal.map <- mask(cal.map, mask)                     
    cat("    ")
    # Tirer des points d'échantillon ...
    cal.pts <- sample.map(cal.map, n.cal.map)
    cat("Cal-map points: ", nrow(cal.pts), "from", cal.map@file@name, "/", 
        SAMPLE.DIST, "px\n", file=txtfile, append=TRUE)
    # ... et ajouter aux points d'entraînement
    if(is.null(train.pts)) {
      train.pts <- cal.pts                                      
    } else {
      train.pts <- rbind(train.pts, cal.pts)
    }
  }
  # Nombre total des points d'entraînement
  cat("Total points:   ", nrow(train.pts), "\n", file=txtfile, append=TRUE)
  
  
  # Extraire les variables correspondantes ----------------
  
  # Extraire les variables Landsat ...
  cat("    -Extracting pixel values for bands:", preds, "... ")
  registerDoParallel(n.cores)
  train.pts <- foreach(subset=split(train.pts, f=rep_len(1:n.cores, nrow(train.pts))), .combine=rbind) %dopar% {
    raster::extract(image, subset, sp=TRUE)
  }
  stopImplicitCluster()
  
  # Ignorer des lignes avec NAs
  train.dat <- na.omit(train.pts@data)[, c("CLASS", preds)]
  
  
  # Calibration du modèle Random Forest -------------------
  
  # Variable catégorielle -> mode de classification, autrement -> régression
  if(type=="classification") train.dat[,1] <- as.factor(train.dat[,1])
  cat("done\n")
  cat("    -Calibrating RandomForest ... ")
  sink(txtfile, append=TRUE)
  # Utiliser caret::train pour validation croisée (a besoin de beaucoup de temps)
  map.model <- ranger(y=train.dat[,1], x=train.dat[,-1], 
                        importance="impurity", mtry = 2,
                        probability = is.factor(train.dat[,1]),
                        num.threads = n.cores,
                        keep.inbag = TRUE)
  print(map.model)
  cat("\n")
  imp <- importance(map.model)
  print(imp[order(imp, decreasing = TRUE)])
  sink()
  # Mesures des erreurs
  cat("Prediction error:", round(map.model$prediction.error, 2), "\n")
  
  # Classification de la carte forêt/non-forêt ------------
  dir.create(dirname(filename), recursive=TRUE, showWarnings=FALSE)
  cat("    -Creating map ... ")
  
  r_in <- readStart(image)
  val <- values(r_in)
  readStop(r_in)
  res <- rep(NA, nrow(val))
  ids <- which(complete.cases(val))
  if(type=="classification"){
    res[ids] <- predict(map.model, val[ids,], num.threads=n.cores)$predictions[,as.character(FOREST)]
    prob.map <- writeStart(raster(r_in), sub("[.]tif$", "_prob.tif", filename), datatype="INT2S", format="GTiff", overwrite=TRUE)
    prob.map <- writeValues(prob.map, 10000*res, 1)
    writeStop(prob.map)
    cat("converting to F/NF map ... ")
    map <- writeRaster(FOREST * (prob.map >= 5000), filename, datatype="INT1U", format="GTiff", overwrite=TRUE)
  } else {
    prob.map = NULL
    res[ids] <- predict(map.model, val[ids,], num.threads=n.cores)$predictions
    map <- writeStart(raster(r_in), filename, datatype="INT2S", format="GTiff", overwrite=TRUE)
    map <- writeValues(map, 10000*res, 1)
    writeStop(map)
  }
  cat("done\n")
  
  
  # Option 2: using foreach -- takes too much memory
  # b <- blockSize(image)
  # 
  # registerDoParallel(n.cores)
  # # loop over blocks
  # l.res <- foreach(i=1:length(b$row)) %dopar% {
  #   # read values for block
  #   # format is a matrix with rows the cells values and columns the layers
  #   r_in <- readStart(image)
  #   val <- getValues(r_in, row = b$row[i], nrows = b$nrows[i])
  #   readStop(r_in)
  #   res <- rep(NA, length(val))
  #   ids <- which(complete.cases(val))
  #   res[ids] <- predict(map.model, val[ids,], num.threads=1)$predictions[,as.character(FOREST)]
  #   return(res)
  # }
  # stopImplicitCluster()
  # 
  # map <- writeStart(raster(r_in), filename, datatype="INT2S", format="GTiff", overwrite=TRUE)
  # # write to output file
  # for(i in 1:length(l.res)) {
  #   writeValues(map, 10000*l.res[[i]], b$row[i])
  # }
  # map <- writeStop(map)
  
  # Option 3: using clusterR  -- does not use all nodes, small blocks uneven use of nodes
  # beginCluster(n=n.cores)
  # map <- clusterR(image, raster::predict,
  #                 args=list(model=map.model,
  #                           fun=function(model, data) {
  #                                 predict(model, data, num.threads=1)$predictions[,2]
  #                                }))
  # endCluster()
  #
  # prob.map <- raster::predict(image, map.model, sub("[.]tif$", "_prob.tif", filename), 
  # fun=function(model, data) 10000*predict(model, data, num.threads=n.cores)$predictions[,as.character(FOREST)], 
  # datatype="INT2S", format="GTiff", overwrite=TRUE)
  
  # sauvegarder la carte
  
  # convertir en %, si c'est une carte couverture houppier
  
  
  cat("-- Done: ", basename(filename), "/", date(), " --\n", file=txtfile, append=TRUE)
  
  invisible(list(
    "model"  = map.model,
    "map"    = map,
    "prob"   = prob.map
  ))
}


# COMMENCER LE TRAITEMENT #####################################################

# Charger images 2018 ---------------------------------------------------------
ref.p192 <- brick(paste0(LANDSAT.DIR, "/p192/p192_2018.tif"))
ref.p193 <- brick(paste0(LANDSAT.DIR, "/p193/p193_2018.tif"))
ref.p194 <- brick(paste0(LANDSAT.DIR, "/p194/p194_2018.tif"))
names(ref.p192) <- names(ref.p193) <- names(ref.p194) <- SST.LSBANDS
ref.images <- list(p192=ref.p192, p193=ref.p193, p194=ref.p194)

# Détérminer le nombre de pixels non-NA
N.PIXELS <- list(p192 = sum(!is.na(ref.p192[[1]][])),
                 p193 = sum(!is.na(ref.p193[[1]][])),
                 p194 = sum(!is.na(ref.p194[[1]][])))

# Charger variables bioclim
bioclim.p192 <- brick(paste0(WORLDCLIM.DIR, "/wc2.0_30s_bio_p192.tif"))
bioclim.p193 <- brick(paste0(WORLDCLIM.DIR, "/wc2.0_30s_bio_p193.tif"))
bioclim.p194 <- brick(paste0(WORLDCLIM.DIR, "/wc2.0_30s_bio_p194.tif"))
names(bioclim.p192) <- names(bioclim.p193) <- names(bioclim.p194) <- SST.BIOCLIM
bioclim <- list(p192=bioclim.p192, p193=bioclim.p193, p194=bioclim.p194)


# Charger points d'entraînement -----------------------------------------------
train.plots <- readOGR(paste0(DIR.SST.BDD.TPS, "/COV_parcelles.shp"))
train.plots <- train.plots[!is.na(train.plots$ccov),  # couverture des houppiers!
                           c("PLOTID", "ccov", "img_date", "author")]

# Convertir les polygones des parcelles en points spatiaux (centroïdes)
train.points <- SpatialPointsDataFrame(gCentroid(train.plots, byid=TRUE),  
                                       data.frame(author=train.plots$author, 
                                                  ccov=train.plots$ccov, 
                                                  img_date=as.Date(train.plots$img_date)))

# Pour la calibration de l'image 2018, sélectionner uniquement les points 
# d'entraînement dont la date de l'image GoogleEarth est entre le 1.1.2017 et le 31.12.2019
train.points <- train.points[!is.na(train.points$img_date) & 
                             train.points$img_date > as.Date("2017-01-01") & 
                             train.points$img_date <= as.Date("2019-12-31"), ]

# Illustrer la répartition des observations 
pdf(paste0(REF.DIR, "/training-pts_2018.pdf"))
raster::plot(train.points)
dev.off()

# Extraire les valeurs d'image pour les points d'entraînement (en parallèle)
registerDoParallel(length(ref.images))
train.points <- foreach(i=1:length(ref.images), .combine=rbind) %dopar% {
  pts <- raster::extract(ref.images[[i]], train.points, sp=TRUE)
  pts <- raster::extract(bioclim[[i]], pts, sp=TRUE)
  pts$image <- names(ref.images[i])
  pts[, c("author", "image", "ccov", SST.LSBANDS, SST.BIOCLIM)]
}
stopImplicitCluster()

# Supprimer les lignes avec des NA
train.points <- train.points[!is.na(rowSums(train.points@data[,-(1:2)])), ]

# Réorganiser le tableau des attributs (F10: forêt/non-forêt à 10% / F30: à 30%)
train.points@data <- cbind(train.points@data[,c("image", "ccov")],
                           F10=cut(train.points$ccov, 
                                   breaks=c(0.0,0.1,1.0), labels=c(NONFOR, FOREST), 
                                   right=FALSE, include.lowest=TRUE), 
                           F30=cut(train.points$ccov, 
                                   breaks=c(0.0,0.3,1.0), labels=c(NONFOR, FOREST), 
                                   right=FALSE, include.lowest=TRUE), 
                           train.points@data[,c(SST.LSBANDS, SST.BIOCLIM)])



# Séléction des variables explicatives ----------------------------------------
# TODO: implement with ranger https://github.com/topepo/caret/issues/555 / https://github.com/topepo/caret/issues/1062

# cov.varsel <- rfe(y = train.points@data[train.points$image=="p193", "ccov"],
#                   x = train.points@data[train.points$image=="p193", PREDICTORS],
#                   sizes = c(4, 6, 8, 10),
#                   rfeControl = rfeControl(
#                     functions = rfFuncs,      # utiliser RandomForest
#                     method  = "repeatedcv",   # validation croisée
#                     number  = 10,             # 10-fold
#                     repeats = 3))             # 3 répétitions
# print(cov.varsel)
# predictors(cov.varsel)
# plot(cov.varsel, type=c("g", "o"))


# Carte de la couverture des houppiers pour 2018 (TODO) -----------------------

set.seed(RSEED)
p193.2018.cov <- classify.image(image     = load.image("/p193/p193_2018.tif"), 
                                bioclim   = bioclim[["p193"]],
                                filename  = paste0(REF.DIR, "/p193_2018_COV_R.tif"),
                                train.pts = train.points[train.points$image == "p193", "ccov"],
                                preds     = PREDICTORS,
                                type      = "treecover")

# observé vs. prédit
pdf(paste0(REF.DIR, "/p193_2018_COV_R.pdf"))
plot(p193.2018.cov[["model"]]$y, p193.2018.cov[["model"]]$predicted, xlim=c(0,1), ylim=c(0,1), 
     main="Couverture houppier p193", xlab="Observation", ylab="Prédiction")
abline(0,1)
dev.off()

# Cartes de référence du couvert forestier 2018 -------------------------------

# Chemin WRS p193
set.seed(RSEED)
classify.image(image     = load.image("/p193/p193_2018.tif"), 
               bioclim   = bioclim[["p193"]],
               filename  = paste0(REF.DIR, "/FC", COV.FC, "/p193_2018_FC", COV.FC, "_R.tif"),
               train.pts = train.points[train.points$image == "p193", paste0("F", COV.FC)],
               preds     = PREDICTORS),
               prob      = TRUE,
               crossval  = TRUE,
               n.cores   = 32)  

# Chemins WRS p192 and p194, utilisant p193 pour calibration
set.seed(RSEED)
registerDoParallel(CORES-1)
foreach(path=c("p192", "p194")) %dopar% {   # traitement en parallèle
  train.pts <- train.points[train.points$image == path, paste0("F", COV.FC)]
  classify.image(image     = load.image(paste0("/", path, "/", path, "_2018.tif")), 
                 bioclim   = bioclim[[path]],
                 filename  = paste0(REF.DIR, "/FC", COV.FC, "/", path, "_2018_FC", COV.FC, "_R.tif"),
                 train.pts = train.pts,
                 # calibration avec carte p193 ...
                 cal.map   = raster(paste0(REF.DIR, "/FC", COV.FC, "/p193_2018_FC", COV.FC, "_R.tif")),
                 n.cal.map = max(2000, nrow(train.pts)/CAL.RATIO),  # ... avec au moin 2000 points
                 mask      = TGO,
                 preds     = PREDICTORS,
                 prob      = TRUE,
                 n.cores   = 32)
}
stopImplicitCluster()

# Cartes de référence du couvert forestier 2003 -----------------------------------------

# Chemin WRS p193
set.seed(RSEED)
classify.image(image     = load.image("/p193/p193_2003.tif"), 
               bioclim   = bioclim[["p193"]],
               filename  = paste0(REF.DIR, "/FC", COV.FC, "/p193_2003_FC", COV.FC, "_R.tif"),
               # sur base de la carte de référence 2018
               ref.map   = raster(paste0(REF.DIR, "/FC", COV.FC, "/p193_2018_FC", COV.FC, "_R.tif")),
               n.ref.map = 2 * SAMPLE.RATIO * N.PIXELS[["p193"]],
               preds     = PREDICTORS,
               mask      = TGO,
               n.cores   = 30)

# Chemins WRS p192 and p194, utilisant p193 pour calibration
set.seed(RSEED)
registerDoParallel(CORES)
foreach(path=c("p192", "p194")) %dopar% {   # traitement en parallèle
  classify.image(image     = load.image(paste0("/", path, "/", path, "_2003_m.tif")), 
                 bioclim   = bioclim[[path]],
                 filename  = paste0(REF.DIR, "/FC", COV.FC, "/", path, "_2003_FC", COV.FC, "_R.tif"),
                 # sur base de la carte de référence 2018 ...
                 ref.map   = raster(paste0(REF.DIR, "/FC", COV.FC, "/", path, "_2018_FC", COV.FC, "_R.tif")),
                 n.ref.map = 2 * (1 - CAL.RATIO) * SAMPLE.RATIO * N.PIXELS[[path]],
                 # ... et calibration avec carte p193 ...
                 cal.map   = raster(paste0(REF.DIR, "/FC", COV.FC, "/p193_2003_FC", COV.FC, "_R.tif")),
                 n.cal.map = 2 * CAL.RATIO * SAMPLE.RATIO * N.PIXELS[[path]],
                 preds     = PREDICTORS,
                 mask      = TGO,
                 n.cores   = 32)
}


# Cartes du couvert forestier brutes pour toutes les dates --------------------

# Chemin WRS p193
set.seed(RSEED)
registerDoParallel(CORES)
foreach(file=dir(paste0(LANDSAT.DIR, "/p193"), pattern="\\_[[:digit:]]+\\_m\\.tif")) %dopar% {
  classify.image(image     = load.image(paste0("/p193/", file)), 
                 bioclim   = bioclim[["p193"]],
                 filename  = paste0(RAW.DIR, "/FC", COV.FC, "/p193/", sub("\\_m\\.tif$", paste0("_F", COV.FC, "r.tif"), file)),
                 ref.map   = raster(paste0(REF.DIR, "/FC", COV.FC, "/p193_2003_FC", COV.FC, "_R.tif")),
                 n.ref.map = SAMPLE.RATIO * N.PIXELS[["p193"]],
                 preds     = PREDICTORS,
                 mask      = TGO,
                 n.cores   = 6)
}

# fusionner les deux tuiles p193_1990
merge(raster(paste0(RAW.DIR, "/FC", COV.FC, "/p193/p193_1990_1_F", COV.FC, "r.tif")), 
      raster(paste0(RAW.DIR, "/FC", COV.FC, "/p193/p193_1990_2_F", COV.FC, "r.tif")),
      filename=paste0(RAW.DIR, "/FC", COV.FC, "/p193/p193_1990_F", COV.FC, "r.tif"), 
      format="GTiff", datatype="INT2U", overwrite=TRUE)

# Chemins WRS p192 and p194, utilisant ...
set.seed(RSEED)
registerDoParallel(CORES)
foreach(file=c(dir(paste0(LANDSAT.DIR, "/p192"), pattern="\\_[[:digit:]]+\\_m\\.tif"),
               dir(paste0(LANDSAT.DIR, "/p194"), pattern="\\_[[:digit:]]+\\_m\\.tif"))) %dopar% { 
  path <- sub("\\_.*", "", file)
  # ... carte p193 pour la calibration (s'il existe, pour la même année)
  if(file.exists(paste0(RAW.DIR, "/FC", COV.FC, "/p193/", sub("\\_m\\.tif$", paste0("_F", COV.FC, "r.tif"), sub(path, "p193", file))))) {
    cal.map   <- raster(paste0(RAW.DIR, "/FC", COV.FC, "/p193/", sub("\\_m\\.tif$", paste0("_F", COV.FC, "r.tif"), sub(path, "p193", file))))
    n.cal.map <- CAL.RATIO * SAMPLE.RATIO * N.PIXELS[[path]]
  } else {
    cal.map <- NULL
    n.cal.map <- NULL
  }
  classify.image(image     = load.image(paste0("/", path, "/", file)), 
                 bioclim   = bioclim[[path]],
                 filename  = paste0(RAW.DIR, "/FC", COV.FC, "/", path, "/", sub("\\_m\\.tif$", paste0("_F", COV.FC, "r.tif"), file)),
                 ref.map   = raster(paste0(REF.DIR, "/FC", COV.FC, "/", path, "_2003_FC", COV.FC, "_R.tif")),
                 n.ref.map = (1 - CAL.RATIO) * SAMPLE.RATIO * N.PIXELS[[path]],
                 cal.map   = cal.map,
                 n.cal.map = n.cal.map,
                 preds     = PREDICTORS,
                 mask      = TGO,
                 n.cores   = 6)
} 


# Fusionner les cartes des chemins p192, p193 et p194 pour dates clés ... -----

for(year in YEARS.REF) {
  merge(mask(crop(brick(paste0(RAW.DIR, "/FC", COV.FC, "/p193/p193_", year, "_F", COV.FC, "r.tif")),TGO), TGO),
        mask(crop(brick(paste0(RAW.DIR, "/FC", COV.FC, "/p192/p192_", year, "_F", COV.FC, "r.tif")),TGO), TGO),
        mask(crop(brick(paste0(RAW.DIR, "/FC", COV.FC, "/p194/p194_", year, "_F", COV.FC, "r.tif")),TGO), TGO),
        filename=paste0(RAW.DIR, "/FC", COV.FC, "/TGO/TGO_", year, "_F", COV.FC, "r.tif"), overwrite=TRUE)
}        

# ... et le cartes de référence 2018 et 2003
for(map in c(paste0("2018_FC", COV.FC, "_R"), 
             paste0("2018_FC", COV.FC, "_R_prob"),
             paste0("2003_FC", COV.FC, "_R"))) {
  merge(mask(crop(brick(paste0(REF.DIR, "/FC", COV.FC, "/p193_", map, ".tif")),TGO), TGO),
        mask(crop(brick(paste0(REF.DIR, "/FC", COV.FC, "/p192_", map, ".tif")),TGO), TGO),
        mask(crop(brick(paste0(REF.DIR, "/FC", COV.FC, "/p194_", map, ".tif")),TGO), TGO),
        filename=paste0(REF.DIR, "/FC", COV.FC, "/TGO_", map, ".tif"), overwrite=TRUE)
}
