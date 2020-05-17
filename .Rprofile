###############################################################################
# .Rprofile: Préparer l'environnement (libraries, variables)
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020


# Charger libraries  ==========================================================
library("sp")           # Classes et méthodes pour les données spatiales
library("rgdal")        # Geospatial Data Abstraction Library
library("raster")       # Analyse et modélisation des données géographiques
library("randomForest") # Algorithme de classification et régression
library("caret")        # Outils pour classification et régression
library("openxlsx")     # Lire et écrire des fichiers Excel (xlsx)
library("dplyr")        # Fonctions pour manipuler des données
library("tidyr")        # Fonctions pour reorganiser des données
library("ggplot2")      # Production des figures
library("foreach")      # Faire des calcules en parallèle ...
library("doParallel")   # ... sur plusieurs processeurs


# Créer un environnement caché  ===============================================
.snsf = new.env()

# Années / Périodes -----------------------------

.snsf$YEARS.ALL <-   1985:2019                                  # - tous
.snsf$YEARS.JNT <- c(1987, 2003, 2005, 2007, 2015, 2017, 2018)  # - conjointes
.snsf$YEARS.REF <- c(1987, 2003,             2015,       2018)  # - référence
.snsf$YEARS.NRF <- c(      2003,             2015,       2018)  # - NRF

# Répertoires -----------------------------------

.snsf$DIR.RAW.DAT   <- "../data"

.snsf$DIR.SST       <- "./01_SSTS"
.snsf$DIR.SST.DAT   <- paste0(.snsf$DIR.SST, "/01_data")
.snsf$DIR.SST.BDD   <- paste0(.snsf$DIR.SST, "/02_BdD")

.snsf$DIR.IFN       <- "./02_IFN"
.snsf$DIR.IFN.DAT   <- paste0(.snsf$DIR.IFN, "/01_field-data")

.snsf$DIR.MRV       <- "./03_NRF-MRV"
.snsf$DIR.MRV.MCF   <- paste0(.snsf$DIR.MRV, "/01_MCF")
.snsf$DIR.MRV.AGB   <- paste0(.snsf$DIR.MRV, "/02_AGB")

# Système de référence des coordonnées ----------

.snsf$UTM.30 <- crs("+proj=utm +zone=30 +datum=WGS84
                     +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

.snsf$UTM.31 <- crs("+proj=utm +zone=31 +datum=WGS84 
                     +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Domaine d'intérêt  ----------------------------

.snsf$TGO     <- spTransform(
                   readOGR(paste0(.snsf$DIR.RAW.DAT, "/GADM/gadm36_TGO_0.shp")),
                     .snsf$UTM.31
                 )

.snsf$TGO.EXT <- extent(151155, 373005, 670665, 1238175)  # xmin, xmax, ymin, ymax

# Noms des couches de données -------------------

.snsf$SST.LSBANDS     <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2", 
                           "evi", "msavi", "nbr", "nbr2", "ndmi", "ndvi", "savi")

.snsf$SST.BIOCLIM     <- paste0("BIO", 1:19)

# Codes pour les classes ------------------------

.snsf$NONFOR <- 0   # non-forêt
.snsf$PREGEN <- 1   # régénération potentielle
.snsf$REGEN  <- 2   # régénération
.snsf$FOREST <- 3   # forêt / forêt initiale

# Divers ----------------------------------------

# Processeurs disponibles pour le calcul parallèle
.snsf$CORES <- detectCores()

# Semence pour le générateur de nombres aléatoires
.snsf$RSEED <- 20191114

# Attacher l'environnement
attach(.snsf)


# Message de bienvenue  =======================================================

message("
=> chargé librairies et variables définit en .Rprofile

.oPYo. o    o .oPYo.  ooooo   ooooo                      
8      8b   8 8       8         8                        
`Yooo. 8`b  8 `Yooo. o8oo       8   .oPYo. .oPYo. .oPYo. 
    `8 8 `b 8     `8  8         8   8    8 8    8 8    8 
     8 8  `b8      8  8         8   8    8 8    8 8    8 
`YooP' 8   `8 `YooP'  8         8   `YooP' `YooP8 `YooP' 
:.....:..:::..:.....::..::::::::..:::.....::....8 :.....:
:::::::::::::::::::::::::::::::::::::::::::::ooP'.:::::::
:::::::::::::::::::::::::::::::::::::::::::::...:::::::::
\n",
date(),
"\n"
)




