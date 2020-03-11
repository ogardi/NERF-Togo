################################################################
# .Rprofile: oading libraries and defining functions
# --------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 11 March 2019

# Libraries  -----------------
library("sp")           # Classes and methods for spatial data
library("rgdal")        # Bindings for the Geospatial Data Abstraction Library
library("rgeos")        # Interface to Geometry Engine - Open Source
library("raster")       # Geographic Data Analysis and Modeling
library("randomForest") # Breiman and Cutler's random forests for classification and regression
library("RStoolbox")    # for terrain correction (topCor)
library("caret")        # for confusion matrix
library("openxlsx")     # pour directement lire des fichiers excel (xlsx)
library("dplyr")        # pour cr??er des tableaux crois??s
library("tidyr")
library("ggplot2")
library("foreach")
library("doParallel")
library("gdalUtils")
library("stringr")
library("maptools")     # kmlPolygon
library("spsurvey")

.snsf = new.env()

# Years ----------------------------------

.snsf$YEARS.ALL   <-   1985:2019
.snsf$YEARS.JNT   <- c(1987,             2003, 2005, 2007,             2015, 2017, 2018)
.snsf$YEARS.REF   <- c(1987,             2003,                         2015,       2018)
.snsf$YEARS.NRF   <- c(                  2003,                         2015,       2018)

# Directories ----------------------------------

.snsf$DIR.RAW.DAT   <- "../data"
.snsf$DIR.SST       <- "./01_SSTS"
.snsf$DIR.IFN       <- "./02_IFN"
.snsf$DIR.MRV       <- "./03_NRF-MRV"

.snsf$DIR.SST.DAT   <- paste0(.snsf$DIR.SST, "/01_data")
.snsf$DIR.SST.BDD   <- paste0(.snsf$DIR.SST, "/02_BdD")
.snsf$DIR.IFN.DAT   <- paste0(.snsf$DIR.IFN, "/01_field-data")

# CRS, AOI, Extents -------------------
.snsf$UTM.30 <- crs("+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
.snsf$UTM.31 <- crs("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Layers of SSTS images
.snsf$SST.LSBANDS     <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", "evi", "msavi", "nbr", "nbr2", "ndmi", "ndvi", "savi")
.snsf$SST.BIOCLIM     <- paste0("BIO", 1:19)

# Boundary of Togo
.snsf$TGO     <- spTransform(readOGR(paste0(.snsf$DIR.RAW.DAT, "/GADM/gadm36_TGO_0.shp")), .snsf$UTM.31)
.snsf$TGO.REG <- spTransform(readOGR(paste0(.snsf$DIR.RAW.DAT, "/GADM/gadm36_TGO_1.shp")), .snsf$UTM.31)
.snsf$TGO.EXT <- extent(151155, 373005, 670665, 1238175)  # define extent by xmin, xmax, ymin and ymax


# Misc -------------------
# prepare for parallel computing
.snsf$CORES <- detectCores()

# Seed for random number generator
.snsf$RSEED <- 20191114

attach(.snsf)

message("
=> chargé librairies et variables définit en .Rprofile
  
###############################
# BIENVENUE DANS LE SNSF TOGO #
###############################
\n",
date(),
"\n"
)




