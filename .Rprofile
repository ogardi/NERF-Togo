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
library("rgeos")        # Opérations spatiales
library("rasterVis")    # illuster raster catégoriel
library("randomForest") # Algorithme de classification et régression
library("ranger")       # Algorithme de classification et régression efficace
library("caret")        # Outils pour classification et régression
library("openxlsx")     # Lire et écrire des fichiers Excel (xlsx)
library("dplyr")        # Fonctions pour manipuler des données
library("tidyr")        # Fonctions pour reorganiser des données
library("stringr")      # Pour fonction str_pad
library("ggplot2")      # Production des figures
library("RColorBrewer") # palettes de couleurs
library("foreach")      # Faire des calcules en parallèle ...
library("doParallel")   # ... sur plusieurs processeurs
library("knitr")        # pour la documentation html


# Créer un environnement caché  ===============================================
.snsf = new.env()

# Années / Périodes -----------------------------

.snsf$YEARS.ALL <-   1985:2021                                  # - tous
.snsf$YEARS.JNT <- c(1987, 2003, 2005, 2007, 2015, 2017, 2018, 2019, 2020)  # - conjointes
.snsf$YEARS.REF <- c(1987, 2003,             2015,       2018, 2019, 2020)  # - référence
.snsf$YEARS.NRF <- c(      2003,             2015,       2018, 2019, 2020)  # - NRF

# Répertoires -----------------------------------

.snsf$DIR.RAW.DAT   <- "../data"

.snsf$DIR.SST       <- "./01_SSTS"
.snsf$DIR.SST.DAT   <- paste0(.snsf$DIR.SST, "/01_data")
.snsf$DIR.SST.DAT.LST   <- paste0(.snsf$DIR.SST.DAT, "/Landsat")
.snsf$DIR.SST.DAT.WC2   <- paste0(.snsf$DIR.SST.DAT, "/Worldclim")

.snsf$DIR.SST.BDD   <- paste0(.snsf$DIR.SST, "/02_BdD")
.snsf$DIR.SST.BDD.GRD   <- paste0(.snsf$DIR.SST.BDD, "/01_reseau-SSTS")
.snsf$DIR.SST.BDD.TPS   <- paste0(.snsf$DIR.SST.BDD, "/02_train-plots")
.snsf$DIR.SST.BDD.VPS   <- paste0(.snsf$DIR.SST.BDD, "/03_val-plots")

.snsf$DIR.IFN       <- "./02_IFN"
.snsf$DIR.IFN.DAT   <- paste0(.snsf$DIR.IFN, "/01_field-data")

.snsf$DIR.MRV       <- "./03_NRF-MRV"
.snsf$DIR.MRV.MCF   <- paste0(.snsf$DIR.MRV, "/01_MCF")
.snsf$DIR.MRV.MCF.REF   <- paste0(.snsf$DIR.MRV.MCF, "/01_ref-maps")
.snsf$DIR.MRV.MCF.RAW   <- paste0(.snsf$DIR.MRV.MCF, "/02_raw-maps")
.snsf$DIR.MRV.MCF.CLN   <- paste0(.snsf$DIR.MRV.MCF, "/03_cln-maps")
.snsf$DIR.MRV.MCF.VAL   <- paste0(.snsf$DIR.MRV.MCF, "/04_validation")
.snsf$DIR.MRV.MCF.RES   <- paste0(.snsf$DIR.MRV.MCF, "/05_results")

.snsf$DIR.MRV.AGB   <- paste0(.snsf$DIR.MRV, "/02_AGB")
.snsf$DIR.MRV.AGB.REF   <- paste0(.snsf$DIR.MRV.AGB, "/01_ref-maps")
.snsf$DIR.MRV.AGB.RES   <- paste0(.snsf$DIR.MRV.AGB, "/02_results")

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

.snsf$TGO.REG     <- spTransform(
                      readOGR(paste0(.snsf$DIR.RAW.DAT, "/GADM/gadm36_TGO_1.shp")),
                        .snsf$UTM.31
                    )

.snsf$TGO.EXT <-  extent(151155, 373005,  670665, 1238175)  # xmin, xmax, ymin, ymax

.snsf$WRS.EXT <- list(
  p192 = extent(216885, 373005,  670665, 1075215),
  p193 = extent(151155, 373005,  683085, 1238175),
  p194 = extent(151155, 222165, 1017495, 1238175))

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
.snsf$RSEED <- 20210316

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




