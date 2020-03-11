##################################################################################
# SSTS/BdD/create-grid.R: create a grid of SSTS plots
# --------------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 11 March 2020

RES <- 480
OUT.DIR <- paste0(DIR.SST.BDD, "/01_reseau-SSTS")
if(!dir.exists(OUT.DIR)) dir.create(OUT.DIR)

# Create a grid of observation points over whole Togo, based on Landsat images (30 x 30m) each 480 m -------------------

x.min <- RES * extent(TGO)@xmin %/% RES
x.max <- RES * extent(TGO)@xmax %/% RES
y.min <- RES * extent(TGO)@ymin %/% RES
y.max <- RES * extent(TGO)@ymax %/% RES

frame.points <- SpatialPoints(expand.grid(seq(x.min, x.max, by=RES), seq(y.min, y.max, by=RES)), proj4string=UTM.31)[TGO]

# add attributes
frame.points$PLOTID  <- paste0(str_pad(frame.points@coords[,1], 7, "left", "0"), "_", str_pad(frame.points@coords[,2], 7, "left", "0"))
frame.points$xcoords <- frame.points@coords[,1]
frame.points$ycoords <- frame.points@coords[,2]

writeOGR(frame.points, dsn=OUT.DIR, layer="TGO_frame_480m", driver="ESRI Shapefile", overwrite=TRUE)