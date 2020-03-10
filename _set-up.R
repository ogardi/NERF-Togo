################################################################
# NERF_Togo/0_set-up.R: Loading libraries and defining functions
# --------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 November 2019

Sys.setlocale("LC_CTYPE", "de_CH.UTF-8")

# set working directory
setwd("~/NERF_Togo/NERF_v1/R")

# Hidden environment for variables and functions 
.env = new.env()

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


# Years ----------------------------------
# .env$YEARS       <- c(1987, 1991, 2000, 2003, 2005, 2007, 2010, 2013, 2015, 2017, 2018)

.env$PERIOD      <- 1985:2019
.env$JNT.YEARS   <- c(1987,             2003, 2005, 2007,             2015, 2017, 2018)
.env$VAL.YEARS   <- c(1987,             2003,                         2015,       2018)
.env$REF.YEARS   <- c(                  2003,                         2015,       2018)

# Directories ----------------------------------

.env$DATA.DIR   <- "../../data"

.env$INPUT.DIR   <- "../input"
.env$IMAGES.DIR  <- paste0(.env$INPUT.DIR, "/1_images")
.env$TRNPTS.DIR  <- paste0(.env$INPUT.DIR, "/2_train-plots")
.env$VALPTS.DIR  <- paste0(.env$INPUT.DIR, "/3_val-plots")
.env$INVENT.DIR  <- paste0(.env$INPUT.DIR, "/4_IFN")

.env$OUTPUT.DIR  <- "../output"
.env$FCC.DIR     <- paste0(.env$OUTPUT.DIR, "/1_forest-cover")
.env$FCC.REF.DIR <- paste0(.env$FCC.DIR, "/1_ref-maps")
.env$FCC.RAW.DIR <- paste0(.env$FCC.DIR, "/2_raw-maps")
.env$FCC.CLN.DIR <- paste0(.env$FCC.DIR, "/3_clean-maps")
.env$FCC.VAL.DIR <- paste0(.env$FCC.DIR, "/4_validation")
.env$FCC.RES.DIR <- paste0(.env$FCC.DIR, "/5_results")

.env$AGB.DIR     <- paste0(.env$OUTPUT.DIR, "/2_biomass")
.env$AGB.REF.DIR <- paste0(.env$AGB.DIR, "/1_ref-maps")
.env$AGB.RAW.DIR <- paste0(.env$AGB.DIR, "/2_raw-maps")
.env$AGB.CLN.DIR <- paste0(.env$AGB.DIR, "/3_clean-maps")
.env$AGB.VAL.DIR <- paste0(.env$AGB.DIR, "/4_validation")
.env$AGB.RES.DIR <- paste0(.env$AGB.DIR, "/5_results")

.env$BANDS       <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", "evi", "msavi", "nbr", "nbr2", "ndmi", "ndvi", "savi")
.env$BIOCLIM     <- paste0("BIO", 1:19)

# CRS, AOI, Extents -------------------
.env$utm.30 <- crs("+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
.env$utm.31 <- crs("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Boundary of Togo
.env$TGO     <- spTransform(readOGR(paste0(.env$DATA.DIR, "/GADM/gadm36_TGO_0.shp")), .env$utm.31)
.env$TGO.reg <- spTransform(readOGR(paste0(.env$DATA.DIR, "/GADM/gadm36_TGO_1.shp")), .env$utm.31)
.env$TGO.ext <- extent(151155, 373005, 670665, 1238175)  # define extent by xmin, xmax, ymin and ymax

# Cloud and water masks for Landsat 4-7 and Landsat 8   ------------------------
# see Landsat 4-7 and Landsat 8 Surface Reflectance Product Guides

.env$qa.water <- c(
  68, 132,                      # Landsat 4-7
  324, 388, 836, 900, 1348      # Landsat 8
)

.env$qa.shadow <- c(
  72, 136,
  328, 392, 840, 904, 1350
)

.env$qa.ice <- c(
  80, 112, 144, 176,
  336, 368, 400, 432, 848, 880, 912, 944, 1352
)

.env$qa.cloud <- c(
  96, 112, 160, 176, 224,
  352, 368, 416, 432, 480, 864, 880, 928, 944, 992
)

# Misc -------------------
# prepare for parallel computing
.env$numCores <- detectCores()

attach(.env)


