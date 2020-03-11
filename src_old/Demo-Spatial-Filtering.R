###########################################################################
# Filtern Waldkarte Togo mit gdal_sieve.py und C-code 
###########################################################################
# R-Nutzer-Gruppe HAFL
# 22.10.2019
# oliver.gardi@bfh.ch
#

# Rausfiltern der Waldflächen < 0.5ha mit gdal_sieve.py -------------------

# Wald/Nicht-Wald Karten der Jahre 2003, 2015 und 2018 laden (Wald = 1 / Nicht-Wald = 3)
maps <- stack("../output/2_forest-maps/TGO/2_raw-maps/TGO/TGO_2003_F30r.tif",
              "../output/2_forest-maps/TGO/2_raw-maps/TGO/TGO_2018_F30r.tif")

# Erstellen temporärer Dateien für Waldmaske
tmp1 <- tempfile(pattern = "", fileext = ".tif")
tmp2 <- tempfile(pattern = "", fileext = ".tif")

# wo einst Wald war ...
onceforest.map <- raster(maps)
onceforest.map[] <- NA
onceforest.map[sum(maps) != nlayers(maps) * 3] <- 1
writeRaster(onceforest.map, tmp1)                                  # rohe Waldmaske rausschreiben

# Waldinseln < 0.5 ha (6 Pixels) mit gdal_sieve rausfiltern
system(paste0("gdal_sieve.py -st 6 -8 -nomask ", tmp1, " ", tmp2))

onceforest.map.clean <- raster(tmp2)                               # wieder einlesen
onceforest.map.clean[onceforest.map.clean == -2147483648] <- NA    # kleine Bereinigung Rasterdatei (NA value)

# Maskieren der Wald/Nicht-Wald Karten mit der bereinigten Waldmaske
maps.clean <- mask(maps, onceforest.map.clean, updatevalue=3)
writeRaster(maps.clean, filename=paste0("./", names(maps.clean), "_05ha.tif"), bylayer=TRUE, format="GTiff", datatype="INT2U", overwrite=TRUE)



# Rausfiltern von Kleinst-Entwaldungsflächen mit eigenem C-Code -----------

maps <- stack("./TGO_2003_F30r_05ha.tif", "./TGO_2018_F30r_05ha.tif")

fcc <- raster(maps)
fcc[] <- 0                                       # Mit 0 initialisieren 
fcc[maps[[1]] == 1 & maps[[2]] == 3] <- 1        # Pixels mit Entwaldung -> 1
fcc[maps[[1]] == 1 & maps[[2]] == 1] <- 2        # Pixels mit stabilem Wald -> 2


# mit C-Algorithums filtern
# pixels with value 1 will only remain if they are in squares of 4 or adjacent to that
# otherwise, if they are adjacent to a pixel with value '2' they will be converted to '2'
# otherwise, they will be converted to zero

# C code kompilieren und laden
system("R CMD SHLIB ./phcf_filter/phcf_filter.c")
dyn.load("./phcf_filter/phcf_filter.so")

# Wrapper-Funktion zur einfacheren Bedienung der C-Schnittstelle
my_c_filter <- function(defor) {
  defor[] <- .C("phcf_filter", nrow(defor), ncol(defor), as.integer(defor[]))[[3]]
  return(defor)
}

# et voilà, Bereinigung der Maske
fcc.clean <- my_c_filter(fcc)

# Maske auf unsere Karten anwenden
defor <- maps$TGO_2003_F30r_05ha == 1 & maps$TGO_2018_F30r_05ha == 3    # Entwaldungspixels
maps.clean <- maps
maps.clean$TGO_2018_F30r_05ha[defor & fcc.clean == 2] <- 1              # Wald wo Wald auf der bereinigten Karten
maps.clean$TGO_2003_F30r_05ha[defor & fcc.clean == 0] <- 3              # Nichtwald wo Nichtwald auf der bereinigten Karten

writeRaster(maps.clean, filename=paste0("./", names(maps.clean), "_c.tif"), bylayer=TRUE, format="GTiff", datatype="INT2U", overwrite=TRUE)






# # mask for removing "deforestation noise"
# defor.phcf <- Y
# defor.phcf[defor.phcf==13] <- 1                 # setting deforested pixels to 1
# defor.phcf[defor.phcf==11] <- 2                 # forest pixels to 2
# defor.phcf[(defor.phcf!=1 & defor.phcf!=2) | is.na(defor.phcf)] <- 0  # all the rest to 0
# defor.phcf[] <- phcf_pt(defor.phcf)
# 
# # apply filter to the image and write it to the disk
# Y.phcf <- Y
# Y.phcf[Y.phcf==13 & defor.phcf==2] <- 11
# Y.phcf[Y.phcf==13 & defor.phcf==0] <- 33
# Y.phcf <- writeRaster(Y.phcf, filename="./results/class_phcf.tif", format="GTiff", datatype="INT2U", overwrite=TRUE)
# 
# # add: isolated forest pixels in deforested patches -> deforested area
# 
# # Step 4: Filter forest area and forest gaps (not defor!) < 0.5ha (6 pixels = 0.54 hectares), 8-connectedness 
# # TODO: don't filter forest gaps.
# forest <- calc(Y.phcf, fun=function(Y){Y==11 | Y==13 | Y==99}, filename="./results/forest.tif", format="GTiff", datatype="INT2U", overwrite=TRUE)
# system("/Library/Frameworks/GDAL.framework/Programs/gdal_sieve.py -st 6 -8 -nomask ./results/forest.tif ./results/forest_ha.tif")
# forest.ha <- raster("./results/forest_ha.tif")
# Y.ha <- Y.phcf
# Y.ha[forest==1 & forest.ha==0] <- 33            # Prairie where forest smaller than 1 ha
# Y.ha[forest==0 & forest.ha==1] <- 11            # Forest, where prairie smaller than 1 ha
# Y.ha <- writeRaster(Y.ha, filename="./results/class_ha.tif", format="GTiff", datatype="INT2U", overwrite=TRUE)
# unlink(c("./results/forest.tif", "./results/forest_ha.tif"))


