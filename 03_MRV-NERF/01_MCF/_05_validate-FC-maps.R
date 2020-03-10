##########################################################################
# NERF_Togo/FCC/6_validate-fc-maps.R: validate clean forest cover maps
# ------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 20 May 2019

# Variables -------------------------------------------------------

COV.FC         <- 30


ct <- list()          # List of confusion tables to be written to xls file


# Load maps, training points --------------------------------

# Load raw and clean forest / non-forest maps

maps <- stack(c(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/TGO/TGO_", VAL.YEARS, "_F", COV.FC, "r.tif"),
                paste0(FCC.CLN.DIR, "/FC", COV.FC, "/TGO/TGO_", VAL.YEARS, "_F", COV.FC, "c.tif")))

names(maps) <- sub("TGO\\_", "X", sub("\\_F[[:digit:]]{2}", "", names(maps)))


# Load training points used for calibration and extract 2018 values from the maps
train.plots <- readOGR(paste0(TRNPTS.DIR, "/COV_parcelles.shp"))
train.points <- SpatialPointsDataFrame(gCentroid(train.plots, byid=TRUE),  # convert polygons to spatial points (centroids)
                                       data.frame(COV=train.plots$ccov))    # with only COV in the attribute table

train.points <- raster::extract(maps[[c("X2018r","X2018c")]], train.points, sp=TRUE)

# delete rows with NAs
train.points <- train.points[rowSums(is.na(train.points@data)) == 0, ]

# reading validation plots ------------------------

# Load validation points, clean and extract map values
val.plots  <- readOGR(paste0(VALPTS.DIR, "/UOT_parcelles.shp"))
val.points <- SpatialPointsDataFrame(gCentroid(val.plots, byid=TRUE), 
                                     data.frame(author=sub("^.*\\/", "", val.plots$author),
                                                V1987=as.numeric(substr(val.plots$lc_87, 1, 1)), V2003=as.numeric(substr(val.plots$lc_03, 1, 1)), 
                                                V2015=as.numeric(substr(val.plots$lc_15, 1, 1)), V2018=as.numeric(substr(val.plots$lc_18, 1, 1))))

val.points <- raster::extract(maps, val.points, sp=TRUE)

# delete rows with NAs
val.points <- val.points[rowSums(is.na(val.points@data)) == 0, ]

# adjusting classes of reference points ----------------------

# if(COV.FC == 30) {
#   # change the woodland (2) to non-forest
#   val.points$V1987[val.points$V1987 == 2] <- 3
#   val.points$V2003[val.points$V2003 == 2] <- 3
#   val.points$V2015[val.points$V2015 == 2] <- 3
#   val.points$V2018[val.points$V2018 == 2] <- 3
# }

if(COV.FC == 10) {
  # change the woodland (2) to forest
  val.points$V1987[val.points$V1987 == 2] <- 1
  val.points$V2003[val.points$V2003 == 2] <- 1
  val.points$V2015[val.points$V2015 == 2] <- 1
  val.points$V2018[val.points$V2018 == 2] <- 1
}

# set cloud/shadow class (4) to forest or non-forest, according to the map
val.points$V1987r <- val.points$V1987; val.points$V1987r[val.points$V1987 %in% c(2,4)] <- val.points$X1987r[val.points$V1987 %in% c(2,4)]
val.points$V2003r <- val.points$V2003; val.points$V2003r[val.points$V2003 %in% c(2,4)] <- val.points$X2003r[val.points$V2003 %in% c(2,4)]
val.points$V2015r <- val.points$V2015; val.points$V2015r[val.points$V2015 %in% c(2,4)] <- val.points$X2015r[val.points$V2015 %in% c(2,4)]
val.points$V2018r <- val.points$V2018; val.points$V2018r[val.points$V2018 %in% c(2,4)] <- val.points$X2018r[val.points$V2018 %in% c(2,4)]

# same thing for cleaned map
val.points$V1987c <- val.points$V1987; val.points$V1987c[val.points$V1987 %in% c(2,4)] <- val.points$X1987c[val.points$V1987 %in% c(2,4)]
val.points$V2003c <- val.points$V2003; val.points$V2003c[val.points$V2003 %in% c(2,4)] <- val.points$X2003c[val.points$V2003 %in% c(2,4)]
val.points$V2015c <- val.points$V2015; val.points$V2015c[val.points$V2015 %in% c(2,4)] <- val.points$X2015c[val.points$V2015 %in% c(2,4)]
val.points$V2018c <- val.points$V2018; val.points$V2018c[val.points$V2018 %in% c(2,4)] <- val.points$X2018c[val.points$V2018 %in% c(2,4)]



# adjusting classes of maps ----------------------

# for the clean maps, set regeneration (2) of the clean maps to forest/non-forest according the validation points
val.points$X1987c[val.points$X1987c == 2] <- val.points$V1987c[val.points$X1987c == 2]
val.points$X2003c[val.points$X2003c == 2] <- val.points$V2003c[val.points$X2003c == 2] 
val.points$X2015c[val.points$X2015c == 2] <- val.points$V2015c[val.points$X2015c == 2] 
val.points$X2018c[val.points$X2018c == 2] <- val.points$V2018c[val.points$X2018c == 2] 

# set to non-forest where a 2 remains (regeneration on the maps and woodland in validation points)
val.points$X1987c[val.points$X1987c == 2] <- val.points$V1987c[val.points$X1987c == 2] <- 3
val.points$X2003c[val.points$X2003c == 2] <- val.points$V2003c[val.points$X2003c == 2] <- 3
val.points$X2015c[val.points$X2015c == 2] <- val.points$V2015c[val.points$X2015c == 2] <- 3
val.points$X2018c[val.points$X2018c == 2] <- val.points$V2018c[val.points$X2018c == 2] <- 3


val.points@data <- mutate_all(val.points@data, as.factor)


## DO THE WORK -------------------------------------------------------

# Compare forest cover map 2018 with reference 2018  ---------

ref.map <- raster(paste0(FCC.REF.DIR, "/FC", COV.FC, "/TGO_2018_FC", COV.FC, "_R.tif"))
ct[["MAP_REF.18r"]] <- confusionMatrix(as.factor(maps$X2018r[]), as.factor(ref.map[]))

pdf(paste0(FCC.VAL.DIR, "/FC2018_vs_COV2018_FC", COV.FC, ".pdf"))
plot(factor(train.points$X2018r, labels=c("Forest", "Non-Forest")) ~ train.points$COV, xlab="Tree Cover 2018", ylab="Forest Cover Map 2018")
dev.off()

pdf(paste0(FCC.VAL.DIR, "/COV2018_vs_FC2018_FC", COV.FC, ".pdf"))
plot(train.points$COV ~ factor(train.points$X2018r, labels=c("Forest", "Non-Forest")), xlab="Forest Cover Map 2018", ylab="Tree Cover 2018")
dev.off()


# Map validation / confusion matrices ----------------------------------

confMat <- function(xtrans, vtrans) {
  levels <- unique(c(xtrans, vtrans))
  return(confusionMatrix(factor(xtrans, levels=levels[order(levels)]), factor(vtrans, levels=levels[order(levels)])))
}

# val.points.bu <- val.points        # backup
val.points <- val.points[!val.points$author %in% c("7_Aklasson_TOLEBA", "8_Yawo_KONKO", "2_Mamalnassoh_ABIGUIME"), ]

# check for authors whose validation points reduce precision
# print(confMat(paste0("x", tmp$X2003c, tmp$X2015c, tmp$X2018c), 
#              paste0("x", tmp$VX2003, tmp$VX2015, tmp$VX2018))$overall)
# 
# for(author in unique(tmp$author)) {
#   print(author)
#   print(confMat(paste0("x", tmp$X2003c[tmp$author != author], tmp$X2015c[tmp$author != author], tmp$X2018c[tmp$author != author]), 
#                 paste0("x", tmp$VX2003[tmp$author != author], tmp$VX2015[tmp$author != author], tmp$VX2018[tmp$author != author]))$overall)
# }


for(t in c("r", "c")) {
  
  # confusion matrices for individual dates
  ct[[paste0("MAP_VAL.87", t)]]          <- confMat(paste0(val.points[[paste0("X1987", t)]], "x",         "x",         "x"        ), paste0(val.points[[paste0("V1987", t)]], "x",         "x",         "x"        ))
  ct[[paste0("MAP_VAL.03", t)]]          <- confMat(paste0("x",         val.points[[paste0("X2003", t)]], "x",         "x"        ), paste0("x",         val.points[[paste0("V2003", t)]], "x",         "x"        ))
  ct[[paste0("MAP_VAL.15", t)]]          <- confMat(paste0("x",         "x",         val.points[[paste0("X2015", t)]], "x"        ), paste0("x",         "x",         val.points[[paste0("V2015", t)]], "x"        ))
  ct[[paste0("MAP_VAL.18", t)]]          <- confMat(paste0("x",         "x",         "x",         val.points[[paste0("X2018", t)]]), paste0("x",         "x",         "x",         val.points[[paste0("V2018", t)]]))
  
  # confusion matrices for 2-date transitions
  ct[[paste0("MAP_VAL.87.03", t)]]       <- confMat(paste0(val.points[[paste0("X1987", t)]], val.points[[paste0("X2003", t)]], "x",         "x"        ), paste0(val.points[[paste0("V1987", t)]], val.points[[paste0("V2003", t)]], "x",         "x"        ))
  ct[[paste0("MAP_VAL.03.15", t)]]       <- confMat(paste0("x",         val.points[[paste0("X2003", t)]], val.points[[paste0("X2015", t)]], "x"        ), paste0("x",         val.points[[paste0("V2003", t)]], val.points[[paste0("V2015", t)]], "x"        ))
  ct[[paste0("MAP_VAL.15.18", t)]]       <- confMat(paste0("x",         "x",         val.points[[paste0("X2015", t)]], val.points[[paste0("X2018", t)]]), paste0("x",         "x",         val.points[[paste0("V2015", t)]], val.points[[paste0("V2018", t)]]))
  ct[[paste0("MAP_VAL.87.18", t)]]       <- confMat(paste0(val.points[[paste0("X1987", t)]], "x",         "x",         val.points[[paste0("X2018", t)]]), paste0(val.points[[paste0("V1987", t)]], "x",         "x",         val.points[[paste0("V2018", t)]]))
  ct[[paste0("MAP_VAL.03.18", t)]]       <- confMat(paste0("x",         val.points[[paste0("X2003", t)]], "x",         val.points[[paste0("X2018", t)]]), paste0("x",         val.points[[paste0("V2003", t)]], "x",         val.points[[paste0("V2018", t)]]))
  
  # confusion matrices for 3-date transition
  ct[[paste0("MAP_VAL.03.15.18", t)]]    <- confMat(paste0("x",val.points[[paste0("X2003", t)]], val.points[[paste0("X2015", t)]], val.points[[paste0("X2018", t)]]), paste0("x", val.points[[paste0("V2003", t)]], val.points[[paste0("V2015", t)]], val.points[[paste0("V2018", t)]]))
  
  # confusion matrices for 4-date transition
  ct[[paste0("MAP_VAL.87.03.15.18", t)]] <- confMat(paste0(val.points[[paste0("X1987", t)]], val.points[[paste0("X2003", t)]], val.points[[paste0("X2015", t)]], val.points[[paste0("X2018", t)]]), paste0(val.points[[paste0("V1987", t)]], val.points[[paste0("V2003", t)]], val.points[[paste0("V2015", t)]], val.points[[paste0("V2018", t)]]))
  
}

# Validation RapidEye map 2015 --------------------------------

rey <- raster(paste0(DATA.DIR, "/RapidEye/TGO_30m.tif"))
names(rey) <- "R2015"
val.points <- raster::extract(rey, val.points, sp=TRUE)

val.points$Rr2015 <- 3
val.points$Rr2015[val.points$R2015 %in% c(11, 12, 16, 18, 19)] <- 1
val.points$VR2015 <- val.points$V2015
val.points$VR2015[val.points$V2015 == 2] <- val.points$Rr2015[val.points$V2015 == 2]

ct[["REY_VAL.15"]]         <- confMat(paste0("x",        "x",        val.points$Rr2015, "x"       ), paste0("x",         "x",         val.points$VR2015, "x"        ))

# Write confusion tables -------------------------------------
save(ct, file=paste0(FCC.VAL.DIR, "/FC", COV.FC, "_flex_ConfTab.RData"))

# write error matrices to Excel File
write.xlsx(lapply(ct, "[[", "table"),
           file=paste0(FCC.VAL.DIR, "/FC", COV.FC, "_flex_ConfTab.xlsx"),
           colnames=TRUE, overwrite=TRUE) 

