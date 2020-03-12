##########################################################################
# NERF_Togo/FCC/7_fc-maps-accuracy.R: validate clean forest cover maps
# ------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 20 May 2019

# based on OpenForis implementation of Olofsson et al. (2014), written by
# Antonia Ortmann, 20 October, 2014
# Source: https://github.com/openforis/accuracy-assessment/blob/master/Rscripts/error_matrix_analysis.R

VALSET <- "FC30_flex"

# Function for estimating accuracies ----------------------------------
  
accuracy.estimate <- function(areas.map, error.matrix, filename=NULL, pixelsize=30^2/10000) {
  
  # remove "x" from the category 
  colnames(error.matrix) <- rownames(error.matrix) <- gsub("x", "", colnames(error.matrix))
  
  # match category order
  maparea         <- areas.map[match(rownames(error.matrix), names(areas.map))]
  ma              <- error.matrix
  dyn             <- names(maparea)
  
  # calculate the area proportions for each map class
  aoi <- sum(maparea)
  propmaparea <- maparea/aoi
  
  # convert the absolute cross tab into a probability cross tab
  ni. <- rowSums(ma)                  # number of reference points per map class
  propma <-  as.matrix(ma/ni. * as.vector(propmaparea))
  propma[is.nan(propma)] <- 0          # for classes with ni. = 0
  
  # estimate the accuracies now
  OA <- sum(diag(propma))              # overall accuracy (Eq. 1 in Olofsson et al. 2014)
  pe <- 0                              # Agreement by chance ...
  for (i in 1:length(dyn)) {
    pe <- pe + sum(propma[i,]) * sum(propma[,i])
  }
  K  <- (OA - pe) / (1 - pe)           # ... for Cohen's Kappa
  UA <- diag(propma) / rowSums(propma) # user's accuracy (Eq. 2 in Olofsson et al. 2014)
  PA <- diag(propma) / colSums(propma) # producer's accuracy (Eq. 3 in Olofsson et al. 2014)
  
  # estimate confidence intervals for the accuracies
  V_OA <- sum(as.vector(propmaparea)^2 * UA * (1 - UA) / (ni. - 1), na.rm=T)  # variance of overall accuracy (Eq. 5 in Olofsson et al. 2014)
  V_UA <- UA * (1 - UA) / (rowSums(ma) - 1) # variance of user's accuracy (Eq. 6 in Olofsson et al. 2014)
  
  # variance of producer's accuracy (Eq. 7 in Olofsson et al. 2014)
  N.j <- array(0, dim=length(dyn))
  aftersumsign <- array(0, dim=length(dyn))
  for(cj in 1:length(dyn)) {
    N.j[cj] <- sum(maparea / ni. * ma[, cj], na.rm=T)
    aftersumsign[cj] <- sum(maparea[-cj]^2 * ma[-cj, cj] / ni.[-cj] * ( 1 - ma[-cj, cj] / ni.[-cj]) / (ni.[-cj] - 1), na.rm = T)
  }
  V_PA <- 1/N.j^2 * ( 
    maparea^2 * (1-PA)^2 * UA * (1-UA) / (ni.-1) + 
      PA^2 * aftersumsign
  ) 
  V_PA[is.nan(V_PA)] <- 0
  
  # proportional area estimation
  propAreaEst <- colSums(propma)        # proportion of area (Eq. 8 in Olofsson et al. 2014)
  AreaEst <- propAreaEst * sum(maparea) # estimated area
  
  # standard errors of the area estimation (Eq. 10 in Olofsson et al. 2014)
  V_propAreaEst <- array(0, dim=length(dyn))
  for (cj in 1:length(dyn)) {
    V_propAreaEst[cj] <- sum((as.vector(propmaparea) * propma[, cj] - propma[, cj] ^ 2) / ( rowSums(ma) - 1))
  }
  V_propAreaEst[is.na(V_propAreaEst)] <- 0
  
  # produce result tables
  
  res <- list()
  res$PREDICTED_PX      <- as.table(maparea)
  res$ERROR_MATRIX      <- ma
  res$ERROR_MATRIX_PROP <- round(propma, 3)
  res$OVERALL_ACC       <- data.frame(accuracy=round(c(OA, K), 3), 
                                      CI=round(c(1.96 * sqrt(V_OA), NA), 3),
                                      row.names=c("OA", "Kappa"))
  res$CLASSES_ACC <- data.frame(maparea=round(maparea * pixelsize, 3)) # in ha
  res$CLASSES_ACC$prop_maparea    <- round(propmaparea, 3)
  res$CLASSES_ACC$adj_proparea    <- round(propAreaEst, 3)
  res$CLASSES_ACC$CI_adj_proparea <- round(1.96 * sqrt(V_propAreaEst), 3)
  res$CLASSES_ACC$adj_area        <- round(propAreaEst * aoi * pixelsize, 3) # in ha
  res$CLASSES_ACC$CI_adj_area     <- round(1.96 * sqrt(V_propAreaEst) * aoi * pixelsize, 3) # in ha
  res$CLASSES_ACC$UA              <- round(UA, 3)
  res$CLASSES_ACC$CI_UA           <- round(1.96 * sqrt(V_UA), 3)
  res$CLASSES_ACC$PA              <- round(PA, 3)
  res$CLASSES_ACC$CI_PA           <- round(1.96 * sqrt(V_PA), 3)
  
  # write results to Excel File
  if(!is.null(filename)) {
    write.xlsx(res,
               file=filename,
               col.names=TRUE, row.names=TRUE, overwrite=TRUE)
  }
  
  return(res)
  
}


# DO THE WORK ----------------------------------

# load the error matrices (R object ct)
load(paste0(FCC.VAL.DIR, "/ConfTab_", VALSET, ".RData"))

# load the predictions and convert "potential regeneration (2)" to "non-forest (3)"
fc.2003 <- brick(paste0(FCC.CLN.DIR, "/FC30/TGO/TGO_2003_F30cf.tif")); fc.2003[fc.2003==2] <- 3
fc.2015 <- brick(paste0(FCC.CLN.DIR, "/FC30/TGO/TGO_2015_F30cf.tif")); fc.2015[fc.2015==2] <- 3
fc.2018 <- brick(paste0(FCC.CLN.DIR, "/FC30/TGO/TGO_2018_F30cf.tif")); fc.2018[fc.2018==2] <- 3

# create 3-date transition map 
fcc  <-  100 * fc.2003 + 10 * fc.2015 + 1 * fc.2018

# get pixel counts (takes time) and separate for different dates / transitions
freq <- table(fcc[])
tmp <- as.numeric(freq); names(tmp) <- names(freq); freq <- tmp

freq.03.15.18        <- c(freq["111"], freq["113"], freq["131"], freq["133"], freq["311"], freq["313"], freq["331"], freq["333"])
names(freq.03.15.18) <- c("111", "113", "131", "133", "311", "313", "331", "333"); freq.03.15.18[is.na(freq.03.15.18)] <- 0

freq.03.18        <- c(sum(freq["111"],freq["131"], na.rm=T), sum(freq["113"],freq["133"], na.rm=T),
                       sum(freq["311"],freq["331"], na.rm=T), sum(freq["313"],freq["333"], na.rm=T))
names(freq.03.18) <- c("11", "13", "31", "33"); freq.03.18[is.na(freq.03.18)] <- 0

freq.03.15        <- c(sum(freq["111"],freq["113"], na.rm=T), sum(freq["131"],freq["133"], na.rm=T),
                       sum(freq["311"],freq["313"], na.rm=T), sum(freq["331"],freq["333"], na.rm=T))
names(freq.03.15) <- c("11", "13", "31", "33"); freq.03.15[is.na(freq.03.15)] <- 0

freq.15.18        <- c(sum(freq["111"],freq["311"], na.rm=T), sum(freq["113"],freq["313"], na.rm=T),
                       sum(freq["131"],freq["331"], na.rm=T), sum(freq["133"],freq["333"], na.rm=T))
names(freq.15.18) <- c("11", "13", "31", "33"); freq.15.18[is.na(freq.15.18)] <- 0

freq.03           <- c(sum(freq["111"],freq["131"], freq["113"],freq["133"], na.rm=T),
                       sum(freq["311"],freq["331"], freq["313"],freq["333"], na.rm=T))
names(freq.03)    <- c("1", "3"); freq.03[is.na(freq.03)] <- 0
  
freq.15           <- c(sum(freq["111"],freq["311"], freq["113"],freq["313"], na.rm=T),
                       sum(freq["131"],freq["331"], freq["133"],freq["333"], na.rm=T))
names(freq.15)    <- c("1", "3"); freq.15[is.na(freq.15)] <- 0

freq.18           <- c(sum(freq["111"],freq["311"], freq["131"],freq["331"], na.rm=T),
                       sum(freq["113"],freq["313"], freq["133"],freq["333"], na.rm=T))
names(freq.18)    <- c("1", "3"); freq.18[is.na(freq.18)] <- 0

# create accuracy maps ------------------------------------

accuracy.estimate(freq.18, ct$MAP_VAL.18c$table , filename=paste0(FCC.VAL.DIR, "/", VALSET, "_Acc_18.xlsx"))
accuracy.estimate(freq.15, ct$MAP_VAL.15c$table , filename=paste0(FCC.VAL.DIR, "/", VALSET, "_Acc_15.xlsx"))
accuracy.estimate(freq.03, ct$MAP_VAL.03c$table , filename=paste0(FCC.VAL.DIR, "/", VALSET, "_Acc_03.xlsx"))

accuracy.estimate(freq.03.15, ct$MAP_VAL.03.15c$table , filename=paste0(FCC.VAL.DIR, "/", VALSET, "_Acc_03-15.xlsx"))
accuracy.estimate(freq.15.18, ct$MAP_VAL.15.18c$table , filename=paste0(FCC.VAL.DIR, "/", VALSET, "_Acc_15-18.xlsx"))
accuracy.estimate(freq.03.18, ct$MAP_VAL.03.18c$table , filename=paste0(FCC.VAL.DIR, "/", VALSET, "_Acc_03-18.xlsx"))

accuracy.estimate(freq.03.15.18, ct$MAP_VAL.03.15.18c$table , filename=paste0(FCC.VAL.DIR, "/", VALSET, "_Acc_03-15-18.xlsx"))

# calculate forest cover, loss and gains ---------------------

res <- data.frame(year     = c(2003,         2015,             2018),
                  total.ha = c(freq.03["1"], freq.15["1"],     freq.18["1"])     * 30^2/10000,
                  defor.ha = c(NA,           freq.03.15["13"], freq.15.18["13"]) * 30^2/10000,
                  regen.ha = c(NA,           freq.03.15["31"], freq.15.18["31"]) * 30^2/10000)
                  



# Example from Olofsson et al. (2014)
# -----------------------------------
# ct <- matrix(data=c(66,  0,   5,   4,
#                     0, 55,   8,  12,
#                     1,  0, 153,  11,
#                     2,  1,   9, 313),
#              nrow=4,
#              byrow=TRUE,
#              dimnames = list(Prediction=c("FN", "NF", "FF", "NN"), Reference=c("FN", "NF", "FF", "NN")))
# 
# px <- c(FN=200000, NF=150000, FF=3200000, NN=6450000)
# 
# accuracy.estimate(px, ct)

