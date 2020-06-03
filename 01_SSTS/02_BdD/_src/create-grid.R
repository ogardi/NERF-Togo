###############################################################################
# create-grid.R: créer une grille de points d'observation SSTS
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Définitions des variables ===================================================

OUT.DIR <- paste0(DIR.SST.BDD, "/01_reseau-SSTS")
if(!dir.exists(OUT.DIR)) dir.create(OUT.DIR)

RES <- 480   # Résolution de la grille (mètres)

# Grille d'observation ========================================================

# Coordonnées min/max en ligne avec images Landsat
x.min <- RES * extent(TGO)@xmin %/% RES
x.max <- RES * extent(TGO)@xmax %/% RES
y.min <- RES * extent(TGO)@ymin %/% RES
y.max <- RES * extent(TGO)@ymax %/% RES

# Creation de la grille
frame.points <- SpatialPoints(expand.grid(seq(x.min, x.max, by=RES), 
                                          seq(y.min, y.max, by=RES)), 
                              proj4string=UTM.31)[TGO]

# Ajouter attribues PLOTID, xcoords et ycoords
frame.points$xcoords <- frame.points@coords[,1]
frame.points$ycoords <- frame.points@coords[,2]
frame.points$PLOTID  <- paste0(str_pad(frame.points@xcoords, 7, "left", "0"), "_", 
                               str_pad(frame.points@ycoords, 7, "left", "0"))

# Sauveguarder comme fichier Shapefile
writeOGR(frame.points, dsn=OUT.DIR, layer="TGO_frame_480m", 
         driver="ESRI Shapefile", overwrite=TRUE)
