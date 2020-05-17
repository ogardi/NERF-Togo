###############################################################################
# create-train-plots.R: Créer un ensemble de parcelles d'entraînement
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Préparation des variables ===================================================

OUT.DIR <- paste0(DIR.SST.BDD, "/02_train-plots/empty")
if(!dir.exists(OUT.DIR)) dir.create(OUT.DIR, recursive=TRUE)

# Charger la grille d'échantillonnage du SSTS
frame.points <- readOGR(paste0(DIR.SST.BDD, "/01_reseau-SSTS/TGO_frame_480m.shp"))

# Charger NDVI 2018 masquée (sans nuages, ombre, ...)
ndvi.band <- which(SST.LSBANDS == "ndvi")
ndvi <- merge(raster(paste0(DIR.SST.DAT, "/Landsat/p192/p192_2018_m.tif"), band=ndvi.band), 
              raster(paste0(DIR.SST.DAT, "/Landsat/p193/p193_2018_m.tif"), band=ndvi.band),
              raster(paste0(DIR.SST.DAT, "/Landsat/p194/p194_2018_m.tif"), band=ndvi.band))
  
# Extraire le NDVI pour chaque parcelle et le découper en 10 strates
frame.points$ndvi   <- raster::extract(ndvi, frame.points)
frame.points$ndvi_c <- cut(frame.points$ndvi, 10, labels=paste0("s", 0:9))


# Echantillonnage des parcelles d'entraînement ================================

# Initialiser le générateur de nombres aléatoires
set.seed(RSEED)

# Définition d'échantillonnage (1500 échantillons par strate NDVI s0 - s9)
Dsgn.grt <- list("s0"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s1"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s2"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s3"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s4"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s5"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s6"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s7"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s8"=list(panel=c(PanelOne=1500), seltype="Equal"),
                 "s9"=list(panel=c(PanelOne=1500), seltype="Equal")
)

# Échantillonnage spatialement équilibrée (Generalized Random Tesselation)
train.points <- grts(sp.object=frame.points,  # grille d'échantillonnage
                     src.frame="sp.object",   # type de grille (objet spatiale)
                     design=Dsgn.grt,         # définition d'échantillonnage
                     stratum="ndvi_c",        # colonne avec strate
                     type.frame="finite",     # type d'échantillonnage
                     DesignID="train"         # préfixe pour chaque point
)

# Appliquer le système de référence des coordonnées 
proj4string(train.points) <- proj4string(frame.points)

# Mélanger les points d'entraînement et ajouter un identifiant
train.points <- train.points[sample(1:nrow(train.points)), ]
train.points$SAMPLEID <- paste0("trn-", str_pad(string=1:nrow(train.points), 
                                                width = 4, 
                                                pad = "0", 
                                                side = "left"))


# Conversion des points en parcelles et ajouter des attributs =================

# Grille Landsat
landsat.grid <- raster(ndvi)
values(landsat.grid) <- 1

# Selectionner les pixels Landsat correspondantes et convertir en polygone
train.plots <- rasterToPolygons(mask(landsat.grid, train.points))

# Extraire les attributs des points d'entraînement ...
train.plots@data <- over(train.plots, 
                         train.points[, c("PLOTID", "SAMPLEID", 
                                          "xcoords", "ycoords", 
                                          "ndvi", "stratum")])

# ... et compléter avec autres attributs
train.plots$ccov                            <- as.character(NA)  # couverture houppiers
train.plots$img_src <- train.plots$img_date <- as.character(NA)  # source et date d'image 
train.plots$author  <- train.plots$mod_date <- as.character(NA)  # auteur et date

# Sauvegarder parcelles sous format Shapefile et KML
writeOGR(train.plots, dsn=OUT.DIR, layer="COV_parcelles", driver="ESRI Shapefile", overwrite=TRUE)
writeKML(train.plots, kmlname="COV_parcelles", filename=paste0(OUT.DIR, "/COV_parcelles.kml"))


# Créer une grille d'échantillon 7x7 dans chaque parcelle =====================

# Détérminer la grille
grid.size <- 7
res <- res(landsat.grid)[1]
offset <- c(res/grid.size/2 + (0:(grid.size-1))*res/grid.size)

# Diviser les parcelles pour un traitement parallèle
subsets <- split(train.plots, f=1:(CORES-1))
registerDoParallel(CORES-1)
train.grids <- foreach(subset=subsets, .combine=rbind, .multicombine=TRUE) %dopar% { 
  # Créer un couche de points vides ...
  grids <- SpatialPointsDataFrame(data.frame(x = 0, y = 0), 
                                  data=data.frame(PLOTID = 0, 
                                                  SAMPLEID = 0, 
                                                  GRIDPOINT = 0))[-1,]
  # ... et ajoute la grille d'échantillon pour chaque parcelle
  for(p in 1:length(subset)) {
    plot <- subset[p,]
    ext <- extent(plot)
    grids <- bind(grids, SpatialPointsDataFrame(expand.grid(ext@xmin+offset, ext@ymin+offset), 
                                                data=data.frame(PLOTID = plot$PLOTID, 
                                                                SAMPLEID = plot$SAMPLEID, 
                                                                GRIDPOINT = 1:grid.size^2)))
  }
  grids
}

# Appliquer le système de référence des coordonnées
proj4string(train.grids) <- proj4string(train.plots)

# Ajouter un attribut "arbre" (1: oui, 0: non)
train.grids$tree <- as.integer(NA)

# Sauveguarder les grille d'échantillon sous format Shapefile
writeOGR(train.grids, dsn=OUT.DIR, layer="COV_parcelles_grid", 
         driver="ESRI Shapefile", overwrite=TRUE)


# Diviser parcelles et grilles d'échantillon en 10 sous-ensembles =============
# pour le traitement par différents photo-interprètes

subsets <- split(train.plots, f=1:10)
for(i in 1:length(subsets)) {
  # Sauvegarder parcelles sous format Shapefile et KML
  writeOGR(subsets[[i]], dsn=OUT.DIR, layer=paste0("COV_parcelles_", i), 
           driver="ESRI Shapefile", overwrite=TRUE)
  writeKML(subsets[[i]], kmlname=paste0("COV_parcelles_", i) , 
           filename=paste0(OUT.DIR, "/COV_parcelles_", i, ".kml"))
  # Sauvegarder les grilles correspondantes sous format Shapefile
  subset.grids <- train.grids[train.grids$PLOTID %in% subsets[[i]]$PLOTID,]
  writeOGR(subset.grids, dsn=OUT.DIR, layer=paste0("COV_parcelles_", i, "_grid"), 
           driver="ESRI Shapefile", overwrite=TRUE)
}