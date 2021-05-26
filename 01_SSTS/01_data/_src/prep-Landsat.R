###############################################################################
# prep-Landsat.R: lire, nettoyer et empiler des images Landsat
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Définitions des variables ===================================================
IN.DIR  <- paste0(DIR.RAW.DAT, "/Landsat")

OUT.DIR <- paste0(DIR.SST.DAT, "/Landsat")
if(!dir.exists(OUT.DIR)) dir.create(OUT.DIR)

# WRS     <- readOGR(paste0(IN.DIR, "/WRS2/WRS2_descending.shp"))    # Scènes WRS

# Masques de nuages et d'eau pour Landsat 4-7 et Landsat 8   ------------------
# voir les guides des produits de réflexion de surface Landsat 4-7 et Landsat 8
# https://www.usgs.gov/land-resources/nli/landsat/landsat-surface-reflectance

QA.WATER <- c(
  68, 132,                      # Landsat 4-7
  324, 388, 836, 900, 1348      # Landsat 8
)

QA.SHADOW <- c(
  72, 136,
  328, 392, 840, 904, 1350
)

QA.ICE <- c(
  80, 112, 144, 176,
  336, 368, 400, 432, 848, 880, 912, 944, 1352
)

QA.CLOUD <- c(
  96, 112, 160, 176, 224,
  352, 368, 416, 432, 480, 864, 880, 928, 944, 992
)

# Liste des images à préparer et à fusionner ==================================

in.image.list <- list(

  # Chemin WRS p192 -----------------------------

  p192.1986 = list(paste0(IN.DIR, "/192_054/1986/LT051920541986011301T1-SC20190405164223/"),
                   paste0(IN.DIR, "/192_055/1986/LT051920551986011301T1-SC20190405164227/"),
                   paste0(IN.DIR, "/192_056/1986/LT051920561986011301T1-SC20190405164153/")),

  p192.1987 = list(paste0(IN.DIR, "/192_054/1987/LT051920541986123101T1-SC20190405164150/"),
                   paste0(IN.DIR, "/192_055/1987/LT051920551986123101T1-SC20190405163521/"),
                   paste0(IN.DIR, "/192_056/1987/LT051920561986123101T1-SC20190405164444/")),

  p192.1991 = list(paste0(IN.DIR, "/192_054/1991/LT041920541991010301T1-SC20190405164201/"),
                   paste0(IN.DIR, "/192_055/1991/LT041920551991010301T1-SC20190405165911/"),
                   paste0(IN.DIR, "/192_056/1991/LT041920561991010301T1-SC20190405163911/")),

  p192.2001 = list(paste0(IN.DIR, "/192_054/2001/LE071920542000121301T1-SC20190405165521/"),
                   paste0(IN.DIR, "/192_055/2001/LE071920552000121301T1-SC20190405165645/"),
                   paste0(IN.DIR, "/192_056/2001/LE071920562000121301T1-SC20190405164029/")),

  p192.2003 = list(paste0(IN.DIR, "/192_054/2003/LE071920542003010401T1-SC20190520111322/"),
                   paste0(IN.DIR, "/192_055/2003/LE071920552003010401T1-SC20190520100402/"),
                   paste0(IN.DIR, "/192_056/2003/LE071920562003010401T1-SC20190520100206/")),

  p192.2005 = list(paste0(IN.DIR, "/192_054/2005/LE071920542004122401T1-SC20190405165520/"),
                   paste0(IN.DIR, "/192_055/2005/LE071920552004122401T1-SC20190405164050/"),
                   paste0(IN.DIR, "/192_056/2005/LE071920562004122401T1-SC20190405164030/")),

  p192.2007 = list(paste0(IN.DIR, "/192_054/2007/LE071920542006123001T1-SC20190406034211/"),
                   paste0(IN.DIR, "/192_055/2007/LE071920552006123001T1-SC20190406034231/"),
                   paste0(IN.DIR, "/192_056/2007/LE071920562006123001T1-SC20190406034202/")),

  p192.2011 = list(paste0(IN.DIR, "/192_054/2011/LE071920542011011001T1-SC20190406034214/"),
                   paste0(IN.DIR, "/192_055/2011/LE071920552011011001T1-SC20190406034114/"),
                   paste0(IN.DIR, "/192_056/2011/LE071920562011011001T1-SC20190406034155/")),

  p192.2013 = list(paste0(IN.DIR, "/192_054/2013/LE071920542013013101T1-SC20190406034224/"),
                   paste0(IN.DIR, "/192_055/2013/LE071920552013013101T1-SC20190406034046/"),
                   paste0(IN.DIR, "/192_056/2013/LE071920562013013101T1-SC20190406034057/")),

  p192.2015 = list(paste0(IN.DIR, "/192_054/2015/LC081920542015011301T1-SC20190405163446/"),
                   paste0(IN.DIR, "/192_055/2015/LC081920552015011301T1-SC20190405163723/"),
                   paste0(IN.DIR, "/192_056/2015/LC081920562015011301T1-SC20190405164231/")),

  p192.2017 = list(paste0(IN.DIR, "/192_054/2017/LC081920542017021901T1-SC20190405163339/"),
                   paste0(IN.DIR, "/192_055/2017/LC081920552017021901T1-SC20190405163342/"),
                   paste0(IN.DIR, "/192_056/2017/LC081920562017021901T1-SC20190405163222/")),

  p192.2018 = list(paste0(IN.DIR, "/192_054/2018/LC081920542018010501T1-SC20190405164304/"),
                   paste0(IN.DIR, "/192_055/2018/LC081920552018010501T1-SC20190405163402/"),
                   paste0(IN.DIR, "/192_056/2018/LC081920562018010501T1-SC20190405163250/")),

  p192.2019 = list(paste0(IN.DIR, "/192_054/2019/LC081920542018122301T1-SC20190405164258/"),
                   paste0(IN.DIR, "/192_055/2019/LC081920552018122301T1-SC20190405163359/"),
                   paste0(IN.DIR, "/192_056/2019/LC081920562018122301T1-SC20190405163342/")),
  
  p192.2020 = list(paste0(IN.DIR, "/192_054/2020/LC081920542020012701T1-SC20200224112049/"),
                   paste0(IN.DIR, "/192_055/2020/LC081920552020012701T1-SC20200224112213/"),
                   paste0(IN.DIR, "/192_056/2020/LC081920562020012701T1-SC20200224112158/")),
  
  p192.2021 = list(paste0(IN.DIR, "/192_054/2021/LC081920542021011301T1-SC20210311141418/"),
                   paste0(IN.DIR, "/192_055/2021/LC081920552021011301T1-SC20210311141448/"),
                   paste0(IN.DIR, "/192_056/2021/LC081920562021011301T1-SC20210311141403/")),

  # Chemin WRS p193 -----------------------------

  p193.1985 = list(paste0(IN.DIR, "/193_052/1985/LT051930521985030601T1-SC20190520100259/"),
                   paste0(IN.DIR, "/193_053/1985/LT051930531985030601T1-SC20190520100324/"),
                   paste0(IN.DIR, "/193_054/1985/LT051930541985030601T1-SC20190520100340/"),
                   paste0(IN.DIR, "/193_055/1985/LT051930551985030601T1-SC20190520100140/")),

  p193.1987 = list(paste0(IN.DIR, "/193_052/1987/LT051930521987012301T1-SC20190405182322/"),
                   paste0(IN.DIR, "/193_053/1987/LT051930531987012301T1-SC20190405182335/"),
                   paste0(IN.DIR, "/193_054/1987/LT051930541987012301T1-SC20190405182331/"),
                   paste0(IN.DIR, "/193_055/1987/LT051930551987012301T1-SC20190405182328/")),

  # Attention: là nous avons des images avec des dates différentes !
  p193.1990.1 = list(paste0(IN.DIR, "/193_052/1990/LT051930521989112801T1-SC20190520100201/"),
                     paste0(IN.DIR, "/193_053/1990/LT051930531989112801T1-SC20190520100233/")),
  p193.1990.2 = list(paste0(IN.DIR, "/193_054/1991/LT041930541991011001T1-SC20190402043117/"),
                     paste0(IN.DIR, "/193_055/1991/LT041930551991011001T1-SC20190402042453/")),

  p193.2000 = list(paste0(IN.DIR, "/193_052/2000/LE071930522000020401T1-SC20190520100729/"),
                   paste0(IN.DIR, "/193_053/2000/LE071930532000020401T1-SC20190520100345/"),
                   paste0(IN.DIR, "/193_054/2000/LE071930542000020401T1-SC20190402045232/"),
                   paste0(IN.DIR, "/193_055/2000/LE071930552000020401T1-SC20190402043121/")),

  p193.2003 = list(paste0(IN.DIR, "/193_052/2003/LE071930522002122601T1-SC20190405182352/"),
                   paste0(IN.DIR, "/193_053/2003/LE071930532002122601T1-SC20190405182309/"),
                   paste0(IN.DIR, "/193_054/2003/LE071930542002122601T1-SC20190405182226/"),
                   paste0(IN.DIR, "/193_055/2003/LE071930552002122601T1-SC20190405190255/")),

  p193.2005 = list(paste0(IN.DIR, "/193_052/2005/LE071930522005021701T1-SC20190405190117/"),
                   paste0(IN.DIR, "/193_053/2005/LE071930532005021701T1-SC20190405190003/"),
                   paste0(IN.DIR, "/193_054/2005/LE071930542005021701T1-SC20190405182210/"),
                   paste0(IN.DIR, "/193_055/2005/LE071930552005021701T1-SC20190405190021/")),

  p193.2007 = list(paste0(IN.DIR, "/193_052/2007/LE071930522007012201T1-SC20190405182221/"),
                   paste0(IN.DIR, "/193_053/2007/LE071930532007012201T1-SC20190405182607/"),
                   paste0(IN.DIR, "/193_054/2007/LE071930542007012201T1-SC20190405182139/"),
                   paste0(IN.DIR, "/193_055/2007/LE071930552007012201T1-SC20190405182418/")),

  p193.2009 = list(paste0(IN.DIR, "/193_052/2009/LE071930522009012701T1-SC20190405182143/"),
                   paste0(IN.DIR, "/193_053/2009/LE071930532009012701T1-SC20190405182301/"),
                   paste0(IN.DIR, "/193_054/2009/LE071930542009012701T1-SC20190405182103/"),
                   paste0(IN.DIR, "/193_055/2009/LE071930552009012701T1-SC20190405182754/")),

  p193.2013 = list(paste0(IN.DIR, "/193_052/2013/LE071930522013022301T1-SC20190405182200/"),
                   paste0(IN.DIR, "/193_053/2013/LE071930532013022301T1-SC20190405182213/"),
                   paste0(IN.DIR, "/193_054/2013/LE071930542013022301T1-SC20190405182152/"),
                   paste0(IN.DIR, "/193_055/2013/LE071930552013022301T1-SC20190405182331/")),

  p193.2015 = list(paste0(IN.DIR, "/193_052/2015/LC081930522015010401T1-SC20190405181512/"),
                   paste0(IN.DIR, "/193_053/2015/LC081930532015010401T1-SC20190405181751/"),
                   paste0(IN.DIR, "/193_054/2015/LC081930542015010401T1-SC20190402042510/"),
                   paste0(IN.DIR, "/193_055/2015/LC081930552015010401T1-SC20190402042446/")),

  p193.2017 = list(paste0(IN.DIR, "/193_052/2017/LC081930522017012501T1-SC20190405181511/"),
                   paste0(IN.DIR, "/193_053/2017/LC081930532017012501T1-SC20190405181440/"),
                   paste0(IN.DIR, "/193_054/2017/LC081930542017012501T1-SC20190405181458/"),
                   paste0(IN.DIR, "/193_055/2017/LC081930552017012501T1-SC20190405181444/")),

  p193.2018 = list(paste0(IN.DIR, "/193_052/2018/LC081930522018011201T1-SC20190405181524/"),
                   paste0(IN.DIR, "/193_053/2018/LC081930532018011201T1-SC20190405181459/"),
                   paste0(IN.DIR, "/193_054/2018/LC081930542018011201T1-SC20190405181510/"),
                   paste0(IN.DIR, "/193_055/2018/LC081930552018011201T1-SC20190405181442/")),

  p193.2019 = list(paste0(IN.DIR, "/193_052/2019/LC081930522019021601T1-SC20190405183839/"),
                   paste0(IN.DIR, "/193_053/2019/LC081930532019021601T1-SC20190405181518/"),
                   paste0(IN.DIR, "/193_054/2019/LC081930542019021601T1-SC20190405183609/"),
                   paste0(IN.DIR, "/193_055/2019/LC081930552019021601T1-SC20190405181507/")),
  
  p193.2020 = list(paste0(IN.DIR, "/193_052/2020/LC081930522020010201T1-SC20210121071838/"),
                   paste0(IN.DIR, "/193_053/2020/LC081930532020010201T1-SC20200224112134/"),
                   paste0(IN.DIR, "/193_054/2020/LC081930542020010201T1-SC20200224112140/"),
                   paste0(IN.DIR, "/193_055/2020/LC081930552020010201T1-SC20200224112235/")),
  
  p193.2021 = list(paste0(IN.DIR, "/193_052/2021/LC081930522021012001T1-SC20210311141507/"),
                   paste0(IN.DIR, "/193_053/2021/LC081930532021012001T1-SC20210311141519/"),
                   paste0(IN.DIR, "/193_054/2021/LC081930542021012001T1-SC20210311141506/"),
                   paste0(IN.DIR, "/193_055/2021/LC081930552021012001T1-SC20210311141446/")),

  # Chemin WRS p194 -----------------------------

  p194.1986 = list(paste0(IN.DIR, "/194_052/1986/LT051940521986011101T1-SC20190405172804/"),
                   paste0(IN.DIR, "/194_053/1986/LT051940531986011101T1-SC20190405172758/")),

  p194.1987 = list(paste0(IN.DIR, "/194_052/1987/LT051940521986122901T1-SC20190405172903/"),
                   paste0(IN.DIR, "/194_053/1987/LT051940531986122901T1-SC20190405174433/")),

  p194.1997 = list(paste0(IN.DIR, "/194_052/1997/LT051940521997021001T1-SC20190405181746/"),
                   paste0(IN.DIR, "/194_053/1997/LT051940531997021001T1-SC20190405173130/")),

  p194.2000 = list(paste0(IN.DIR, "/194_052/2000/LE071940522000012601T1-SC20190405172721/"),
                   paste0(IN.DIR, "/194_053/2000/LE071940532000012601T1-SC20190405172733/")),

  p194.2003 = list(paste0(IN.DIR, "/194_052/2003/LE071940522002121701T1-SC20190405172823/"),
                   paste0(IN.DIR, "/194_053/2003/LE071940532002121701T1-SC20190405172739/")),

  p194.2005 = list(paste0(IN.DIR, "/194_052/2005/LE071940522004122201T1-SC20190405172700/"),
                   paste0(IN.DIR, "/194_053/2005/LE071940532004122201T1-SC20190405172612/")),

  p194.2007 = list(paste0(IN.DIR, "/194_052/2007/LT051940522007010501T1-SC20190405172919/"),
                   paste0(IN.DIR, "/194_053/2007/LT051940532007010501T1-SC20190405172216/")),

  p194.2010 = list(paste0(IN.DIR, "/194_052/2010/LE071940522010012101T1-SC20190405172745/"),
                   paste0(IN.DIR, "/194_053/2010/LE071940532010012101T1-SC20190405173304/")),
  
  p194.2012 = list(paste0(IN.DIR, "/194_052/2012/LE071940522012011101T1-SC20190405173146/"),
                   paste0(IN.DIR, "/194_053/2012/LE071940532012011101T1-SC20190405172236/")),

  p194.2015 = list(paste0(IN.DIR, "/194_052/2015/LC081940522015012701T1-SC20190405172055/"),
                   paste0(IN.DIR, "/194_053/2015/LC081940532015012701T1-SC20190405172042/")),

  p194.2017 = list(paste0(IN.DIR, "/194_052/2017/LC081940522016123101T1-SC20190405172058/"),
                   paste0(IN.DIR, "/194_053/2017/LC081940532016123101T1-SC20190405172040/")),

  p194.2018 = list(paste0(IN.DIR, "/194_052/2018/LC081940522017121801T1-SC20190405172038/"),
                   paste0(IN.DIR, "/194_053/2018/LC081940532017121801T1-SC20190405174114/")),

  p194.2019 = list(paste0(IN.DIR, "/194_052/2019/LC081940522019012201T1-SC20190405172019/"),
                   paste0(IN.DIR, "/194_053/2019/LC081940532019012201T1-SC20190405172055/")),
  
  p194.2020 = list(paste0(IN.DIR, "/194_052/2020/LC081940522020010901T1-SC20200224112249/"),
                   paste0(IN.DIR, "/194_053/2020/LC081940532020010901T1-SC20200224112242/")),
  
  p194.2021 = list(paste0(IN.DIR, "/194_052/2021/LC081940522021012701T1-SC20210311141413/"),
                   paste0(IN.DIR, "/194_053/2021/LC081940532021012701T1-SC20210311141450/"))

)

# Fonction pour traiter et fusionner un ensemble d'images Landsat =============
#
# @param in.image.dirs liste des répertoires des images à traiter
# @param ext           l'étendue à utiliser pour le recadrage des images
# @param filename      nom du fichier pour la sauvegarde de l'image traitée
# @param overwrite     retraiter et écraser des images déjà existantes
# @param log           écrire les informations sur le processus dans un fichier
#
# @return              objet raster de l'image traitée (invisible)
#
# prepare.landsat <- function(in.image.dirs, 
#                             ext=NULL, 
#                             filename=NULL, 
#                             overwrite=FALSE, 
#                             log=TRUE) {
#   
#   if(!file.exist(filename) | overwrite==TRUE) {
#   
#     # Ouvrir le fichier journal si un nom de fichier et log==true
#     if(!is.null(filename) & log==TRUE) {
#       dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
#       logfile <- file(sub("\\.[[:alnum:]]+$", ".log", filename), open="wt")
#       sink(logfile, type="output")
#       sink(logfile, type="message")
#       message(date())
#     }
#     
#     # Pour chaque image ... 
#     for(image.dir in in.image.dirs) {
#       
#       tmp1 <- paste0(image.dir, "tmp1.tif")
#       tmp2 <- paste0(image.dir, "tmp2.tif")
#       unlink(c(tmp1, tmp2))
#       
#       image.name   <- substr(dir(image.dir)[1], 1, 40)
#       image.sensor <- substr(image.name, 0,4)
#       image.scene  <- paste0(substr(image.name, 11, 13), "_", substr(image.name, 14, 16))
#       image.date   <- substr(image.name, 18, 21)
#       image.path   <- as.numeric(substr(image.scene, 1, 3))
#       image.row    <- as.numeric(substr(image.scene, 5, 7))
#        
#       if(image.sensor=="LC08") {
#         bands <- "^.*_(pixel_qa|band2|band3|band4|band5|band6|band7|evi|msavi|nbr|nbr2|ndmi|ndvi|savi).tif$"
#       } else {
#         bands <- "^.*_(pixel_qa|band1|band2|band3|band4|band5|band7|evi|msavi|nbr|nbr2|ndmi|ndvi|savi).tif$"
#       }
#       
#       message("- ", image.name, ": ", appendLF = FALSE)
#       
#       # ... empiler les couches
#       message("stacking ... ", appendLF = FALSE)
#       system(paste("gdal_merge.py -q -n -9999 -separate -ot Int16 -of GTiff -o", tmp1, paste(grep(bands, dir(image.dir, full.names=TRUE), value=TRUE), collapse=" ")))
#       
#       # ... couper et masquer l'image avec l'étendue et le WRS
#       message("crop/mask ... ", appendLF = FALSE)
#       wrs <- spTransform(WRS[WRS$PATH==image.path & WRS$ROW==image.row, ], CRS=UTM.31)
#       if(!is.null(ext)) wrs <- crop(wrs, ext)
#       wrs.file <- tempfile()
#       writeOGR(wrs, dirname(wrs.file), basename(wrs.file), driver="ESRI Shapefile")
#       system(paste("gdalwarp -q -t_srs EPSG:32631 -cutline", paste0(wrs.file, ".shp"), "-crop_to_cutline", tmp1, tmp2))
#       message("done")
#     }
#     
#     # Fusionner les images (scènes) dans les listes
#     message("- merging to ", basename(filename), " ... ", appendLF = FALSE)
#     system(paste("gdal_merge.py -q -n -9999 -a_nodata -9999 -ot Int16 -of GTiff -o", filename, paste0(unlist(in.image.dirs), "tmp2.tif", collapse=" ")))
#   
#     # Delete temporary files
#     message("unlink tmpfiles ... ", appendLF = FALSE)
#     unlink(c(paste0(unlist(in.image.dirs), "tmp1.tif"), paste0(unlist(in.image.dirs), "tmp2.tif")))
#     message("done")
#     
#     # Fermer le fichier journal
#     if(!is.null(filename) & log==TRUE) {
#       sink(type="output")
#       sink(type="message")
#     }
#   }
# }


prepare.landsat <- function(in.image.dirs, 
                            filename,
                            ext, 
                            overwrite = FALSE,
                            log = TRUE) {
  
  if(!file.exists(filename) | overwrite==TRUE) {
    
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    
    # Ouvrir le fichier journal si log==true
    if(log==TRUE) {
      logfile <- file(sub(".tif", ".log", filename), open="wt")
      sink(logfile, type="output")
      sink(logfile, type="message")
      message(filename, " -- ", date())
    }
    
    # répértoire pour les différents couches
    vrt.dir <- sub(".tif", "", filename)
    dir.create(vrt.dir, showWarnings = FALSE)
    
    # tous les fichiers 
    files <- unlist(lapply(in.image.dirs, dir, full.names = TRUE))
    
    # étendue pour couper
    tmp <- tempfile()
    system(paste("gdalbuildvrt", tmp, paste(grep("_band7.tif$", files, value=TRUE), collapse=" ")))
    img.tmp <- raster(tmp)
    img.ext <- extent(img.tmp)
    p194 <- compareCRS(img.tmp, UTM.30)
    
    if(p194) {
      ext.orig <- ext
      ext <- as(ext, "SpatialPolygons")
      proj4string(ext) <- UTM.31
      ext <- extent(spTransform(ext, UTM.30))
    }
    
    # xmin <- max(ext@xmin, img.ext@xmin)
    # ymin <- max(ext@ymin, img.ext@ymin)
    # xmax <- min(ext@xmax, img.ext@xmax)
    # ymax <- min(ext@ymax, img.ext@ymax)
    
    stripe <- basename(files[1])
    stripe <- paste0(substr(stripe, 1, 13), "xxx", substr(stripe, 17, 40))
    
    sensor <- substr(stripe, 1, 4)
     
    # cominber images, bandes 1-5 (L8 2-6)
    for(b in c(paste0("sr_band", 1:5+(sensor == "LC08")), "sr_band7", "sr_evi", "sr_msavi", "sr_nbr", "sr_nbr2", "sr_ndmi", "sr_ndvi", "sr_savi", "pixel_qa")) {
      system(paste("gdalbuildvrt", paste0(vrt.dir, "/", stripe, "_", b, ".vrt"), 
                   "-te", ext@xmin, ext@ymin, ext@xmax, ext@ymax,
                   paste(grep(paste0("^.*_", b, ".tif$"), files, value=TRUE), collapse=" ")))
    }
    
    # empiler les bandes en virtual raster vrt
    system(paste("gdalbuildvrt -separate", sub(".tif", ".vrt", filename),
                 paste0(vrt.dir, "/*_sr_*.vrt")))
    
    
    # Using single band as mask for vrt is not working!
    # system(paste("gdal_translate -b 1 -b 2 -b 3 -b 4 -b 5 -b 6 -mask 14", sub(".tif", ".vrt", filename), sub(".tif", "_mask.tif", filename)))
    
    # masquer et sauvegarder comme tif 
    maskvalues <- c(QA.CLOUD, QA.ICE, QA.SHADOW, QA.WATER)
    system(paste("pksetmask -i", sub(".tif", ".vrt", filename), 
                           "-m", paste0(vrt.dir, "/*_pixel_qa.vrt"), 
                           paste("-msknodata",  maskvalues, collapse=" "),  
                           paste(rep("-nodata -9999", length(maskvalues)), collapse = " "),
                           "-o", filename))
    
    # transformer en UTM31 (pour p194)
    if(p194) {
      system(paste("gdalwarp -t_srs EPSG:32631",
                   "-tr 30 30",
                   "-te", ext.orig@xmin, ext.orig@ymin, ext.orig@xmax, ext.orig@ymax,
                   filename, sub(".tif", "_repr.tif", filename)))
      file.rename(sub(".tif", "_repr.tif", filename), filename)
    }
    
    # créer vignette
    image <- brick(filename)
    jpeg(sub("[.]tif$", ".jpeg", filename), width=1350, height=3000)
    par(plt=c(0,1,0,1))
    raster::plot(TGO, col="yellow")
    plotRGB(image, r=6, g=5, b=3, stretch="lin", bgalpha = 0, add=TRUE)
    raster::plot(mask(image[[1]], TGO, inverse=TRUE), col="#FFFFFF66", legend=FALSE, add=TRUE)
    raster::plot(TGO, add=TRUE, lwd=3)
    dev.off()
    
    # Fermer le fichier journal
    if(log==TRUE) {
      sink(type="output")
      sink(type="message")
    }
  }
}


# COMMENCER LE TRAITEMENT #####################################################

# Pour chaque ensembles d'images (année/chemin), en parallèle, ...
registerDoParallel(CORES-1)
foreach(i=1:length(in.image.list)) %dopar% {
  
  # ... extraire la liste des images à traiter
  in.image.dirs <- in.image.list[[i]]
  
  name <- unlist(strsplit(names(in.image.list[i]), "[.]"))
  path <- name[1]   # p.ex. "p192"
  year <- name[2]   # p.ex. "1986"
  tile <- name[3]   # p.ex. "NA" (ou 1,2, ...)
  
  # ... répertoire pour sauvegarder l'image
  out.image.dir <- paste0(OUT.DIR, "/", path)
  
  # ... nom du fichier
  if(is.na(tile)) {
    filename <- paste0(out.image.dir, "/", path, "_", year, ".tif")
  } else {
    filename <- paste0(out.image.dir, "/", path, "_", year, "_", tile, ".tif")
  }
  
  # ... et faire le travail
  message("Processing ", path, "_", year)
  prepare.landsat(in.image.dirs, 
                  ext = WRS.EXT[[path]], 
                  filename = filename, 
                  overwrite=FALSE, 
                  log=TRUE)
  
}


# Fusionner les fichiers journaux 
for(d in list.dirs(paste0(OUT.DIR), recursive=FALSE, full.names=TRUE)) {
  path <- basename(d)
  system(paste0("tail -n +1 ", d, "/*.log > ", d, "/", path, ".tmp"))
  system(paste0("rm ", d, "/*.log"))
  system(paste0("mv ", d, "/", path, ".tmp ", d, "/", path, ".log"))
}