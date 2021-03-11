###############################################################################
# prep-SoilGrids.R: lire et reprojeter des données SoilGrids
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Source des données: https://files.isric.org/soilgrids/data/recent/
# Documentation: https://www.isric.org/explore/soilgrids/faq-soilgrids-2017

# Définitions des variables ===================================================

IN.DIR  <- paste0(DIR.RAW.DAT, "/SoilGrids")
OUT.DIR <- paste0(DIR.SST.DAT, "/SoilGrids")
if(!dir.exists(OUT.DIR)) dir.create(OUT.DIR)

# Reprojection oilGrids vers Landsat (résolution 30m, UTM 31, ...) ============

foreach(file=dir(IN.DIR, pattern=".*[.]tif$")) %do% {
  system(paste("gdalwarp -t_srs EPSG:32631",
               paste0(IN.DIR, "/", file),
               "-tr 30 30",
               paste("-te", TGO.EXT@xmin, TGO.EXT@ymin, TGO.EXT@xmax, TGO.EXT@ymax),
               paste0(OUT.DIR, "/", file),
               "-overwrite"))
}


# Créer des vignettes pour les donées SoilGrids ===============================

# Types de sol WRB ------------------------------
jpeg(paste0(OUT.DIR, "/TAXNWRB_250m_ll.jpeg"), width=675, height=1200)
soil.types <- as.factor(raster(paste0(OUT.DIR, "/TAXNWRB_250m_ll.tif")))
cats <- levels(soil.types)[[1]]
attr <- read.csv(paste0(IN.DIR, "/TAXNWRB_250m_ll.tif.csv"))
cats <- merge(cats, attr[,c("Number", "Group", "WRB_group", "R", "G", "B")], 
              by.x = "ID", by.y = "Number")
levels(soil.types) <- cats
levelplot(soil.types, att="Group",
          col.regions=rgb(cats$R, cats$G, cats$B, maxColorValue = 255), 
          xlab="", ylab="") +
  levelplot(mask(soil.types, TGO, inverse=TRUE), col.regions="#FFFFFF66") +
  latticeExtra::layer(sp.polygons(TGO, lwd=3))
dev.off()

# Réservoire carbone à 30 cm --------------------
jpeg(paste0(OUT.DIR, "/OCSTHA_M_30cm_250m_ll.jpeg"), width=570, height=1200)
image <- raster(paste0(OUT.DIR, "/OCSTHA_M_30cm_250m_ll.tif"))
raster::plot(image, zlim=c(30,100))
raster::plot(mask(image, TGO, inverse=TRUE), col="#FFFFFF66", legend=FALSE, add=TRUE)
sp::plot(TGO, add=TRUE, lwd=3)
dev.off()


# Tableau des surfaces des catégories de sol selon GIEC  ======================
tmp <- as.data.frame(table(mask(soil.types, TGO)[]) * 30^2/10000)
tmp <- merge(tmp, cats, by.x="Var1", by.y="ID")
tmp$IPCC <- NA
tmp$IPCC[tmp$WRB_group %in% c("Leptosols","Vertisols","Kastanozems","Chernozems",
                              "Phaeozems","Luvisols","Alisols","Albeluvisols",
                              "Solonetz","Calcisols","Gypsisols","Umbrisols",
                              "Cambisols","Regosols")] <- "HAC"
tmp$IPCC[tmp$WRB_group %in% c("Acrisols","Lixisols","Nitisols","Ferralsols",
                              "Durisols")] <- "LAC"
tmp$IPCC[tmp$WRB_group %in% c("Arenosols")] <- "SAN"
tmp$IPCC[tmp$WRB_group %in% c("Podzols")] <- "POD"
tmp$IPCC[tmp$WRB_group %in% c("Andosols")] <- "VOL"
tmp$IPCC[tmp$WRB_group %in% c("Gleysols")] <- "WET"
tmp$IPCC[is.na(tmp$WRB_group)] <- "Other"
