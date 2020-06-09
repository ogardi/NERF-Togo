###############################################################################
# 04_analyze-ER.R: Analyse des émissions et séquestrations CO2
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Définitions des variables ===================================================

# Rapports racines-tige forêts tropicales sèches selon Mokany et al. (2006)
# https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2486.2005.001043.x
RSR.Y    <- 0.563     # moyenne pour la jeune forêt      
RSR.Y.SE <- 0.086     # écart type pour la jeune forêt  
RSR.O    <- 0.275     # moyenne pour les forêts matures 
RSR.O.SE <- 0.003     # écart type pour les forêts matures
RSR.AGB  <- 20        # Seuil biomasse aérienne jeunes forêts <-> forêts matures

# Teneur en carbone de la biomasse (valeur défault IPCC 2006)
# https://www.ipcc-nggip.iges.or.jp/public/2006gl/pdf/4_Volume4/V4_04_Ch4_Forest_Land.pdf#page=48
C.RATIO  <- 0.47

# Répétoires des cartes biomasse et des résultats
AGB.REF.DIR <- DIR.MRV.AGB.REF
AGB.RES.DIR <- DIR.MRV.AGB.RES



# Analyse des cartes de biomasse aérienne =====================================

# Biomasse aérienne par strate IFN (carte RapidEye) ---------------------------

agb.2015 <- raster(paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.tif"))
strata  <- raster(paste0(DIR.RAW.DAT, "/RapidEye/TGO_30m.tif"))

zonal(agb.2015, strata, fun="mean")
zonal(agb.2015, strata, fun="sd")

# Distribution des biomasses aériennes (densité) ------------------------------

# Charger cartes biomasse aérienne 2003, 2015 et 2018
files <- dir(AGB.REF.DIR, pattern="^TGO.*\\.tif$", full.names=TRUE)
maps <- stack(files)
names(maps) <- substr(names(maps), 5, 8)

# Sélectionner les années 2003 et 2018
maps.dat <- as.data.frame(na.omit(maps[[c("X2003", "X2018")]][]))

# production de la figure
dens.plot <- maps.dat %>% gather(variable, value, X2003, X2018) %>%
  ggplot(aes(x=value, color=variable, linetype=variable)) +
  geom_density() +
  theme_light() 

pdf(paste0(AGB.RES.DIR, "/AGB-densities_03-18.pdf"))
print(dens.plot)
dev.off() 


# Analyze des émissions et séquestrations CO2 2003-2018 =======================

# charger cartes de biomasse 2003 et 2018 et calcule de la biomasse racinaire
agb.2003 <- raster(paste0(AGB.REF.DIR, "/TGO_2003_AGB_R.tif"))
bgb.2003 <- agb.2003 * RSR.O
bgb.2003[agb.2003 <= RSR.AGB] <- agb.2003[agb.2003 <= RSR.AGB] * RSR.Y

agb.2018 <- raster(paste0(AGB.REF.DIR, "/TGO_2018_AGB_R.tif"))
bgb.2018 <- agb.2018 * RSR.O
bgb.2018[agb.2018 <= RSR.AGB] <- agb.2018[agb.2018 <= RSR.AGB] * RSR.Y

# charger cartes des surfaces forestiers 2003 et 2018
fc.2003 <- raster(paste0(DIR.MRV.MCF.CLN, "/FC30/TGO/TGO_2003_F30cf.tif"))
fc.2018 <- raster(paste0(DIR.MRV.MCF.CLN, "/FC30/TGO/TGO_2018_F30cf.tif"))


# Carte des pertes de biomasse ------------------------------------------------

# Masque forêt en 2003 et non-forêt en 2018
defor <- fc.2003 == FOREST & fc.2018 != FOREST       

# moyenne biomasse non-forêt (aérienne et racinaire)
nf.agb <- mean(agb.2018[fc.2018 == 3], na.rm=TRUE)   
nf.bgb <- mean(bgb.2018[fc.2018 == 3], na.rm=TRUE)

# Perte de biomasse due à la déforestation = biomasse 2003 - moyenne biomasse non-forêt
defor.agb <- mask(agb.2003, defor, maskvalue=0) - nf.agb
defor.bgb <- mask(bgb.2003, defor, maskvalue=0) - nf.bgb

# Mettre à 0, là ou la perte est négative
defor.agb[defor.agb < 0] <- 0  
defor.bgb[defor.bgb < 0] <- 0 


# Carte des gains de biomasse -------------------------------------------------

# Masque non-forêt en 2003 et forêt en 2018
regen <- fc.2003 != FOREST & fc.2018 == FOREST

# Gain de biomasse due à la régénération = biomasse 2018 - biomasse 2003
regen.agb.2003 <- mask(agb.2003, regen, maskvalue=0)
regen.bgb.2003 <- mask(bgb.2003, regen, maskvalue=0)
regen.agb.2018 <- mask(agb.2018, regen, maskvalue=0)
regen.bgb.2018 <- mask(bgb.2018, regen, maskvalue=0)
regen.agb <- regen.agb.2018 - regen.agb.2003
regen.bgb <- regen.bgb.2018 - regen.bgb.2003

# Mettre à 0, là ou le gain est négative
regen.agb[regen.agb < 0] <- 0
regen.bgb[regen.bgb < 0] <- 0


# Calculer émissions ets séquestrations pour une région -----------------------
#
# @param defor      Carte des pixels déforestés     
# @param defor.agb  Carte des pertes de biomasse aérienne
# @param defor.bgb  Carte des pertes de biomasse racinaire
# @param regen      Carte des pixels régénérés
# @param regen.agb  Carte des gains de biomasse aérienne
# @param regen.bgb  Carte des gains de biomasse racinaire
# @param aoi        Zone d'intérêt (région)
# @param years      Nombre d'années entre les observations
#
# @return           List avec les éléments 
#                   - defor.area     Surface déforestée
#                   - defor.area.a   Surface déforestée par année
#                   - defor.agb.ha   Perte de biomasse aérienne par hectare déforestée
#                   - defor.agb.a    Perte de biomasse aérienne par année
#                   - defor.bgb.ha   ... même chose pour la biomasse racinaire ...
#                   - defor.bgb.a
#                   - defor.co2.ha   ... et converti en CO2eq ...
#                   - defor.co2.a
#                   - regen.area     Surface régénérée
#                   - regen.area.a   Surface régénérée par année
#                   - regen.agb.ha   Gain de biomasse aérienne par hectare régénérée
#                   - regen.agb.a    Gain de biomasse aérienne par année
#                   - regen.bgb.ha   ... même chose pour la biomasse racinaire ...
#                   - regen.bgb.a
#                   - regen.co2.ha   ... et converti en CO2eq ...
#                   - regen.co2.a
#
#

emissions <- function(defor, defor.agb, defor.bgb, 
                      regen, regen.agb, regen.bgb, aoi=NULL, years=15) {
  
  # couper la zone d'intérêt
  if(!is.null(aoi)) {
    defor     <- mask(crop(defor, aoi), aoi)
    defor.agb <- mask(crop(defor.agb, aoi), aoi)
    defor.bgb <- mask(crop(defor.agb, aoi), aoi)
    regen     <- mask(crop(regen, aoi), aoi)
    regen.agb <- mask(crop(regen.agb, aoi), aoi)
    regen.bgb <- mask(crop(regen.bgb, aoi), aoi)
  }
  
  # convertir biomasse aérienne et racinaire en CO2eq
  defor.co2    <- (defor.agb + defor.bgb) * C.RATIO * 44/12
  # déterminer la surface (taille d'un pixel Landsat = 30mx30m)
  defor.area   <- sum(defor[], na.rm=TRUE) * 30^2 / 10000
  defor.area.a <- defor.area / years
  # moyenne des pertes (par hectare)
  defor.agb.ha <- mean(defor.agb[], na.rm=TRUE)
  defor.bgb.ha <- mean(defor.bgb[], na.rm=TRUE)
  defor.co2.ha <- mean(defor.co2[], na.rm=TRUE)
  
  # la même chose pour la régénération
  regen.co2    <- (regen.agb + regen.bgb) * C.RATIO * 44/12
  regen.area   <- sum(regen[], na.rm=TRUE) * 30^2 / 10000
  regen.area.a <- regen.area / years
  regen.agb.ha <- mean(regen.agb[], na.rm=TRUE)
  regen.bgb.ha <- mean(regen.bgb[], na.rm=TRUE)
  regen.co2.ha <- mean(regen.co2[], na.rm=TRUE)
  
  return(list(
    defor.area   = defor.area, 
    defor.area.a = defor.area.a,
    defor.agb.ha = defor.agb.ha,
    defor.agb.a  = defor.agb.ha * defor.area.a,
    defor.bgb.ha = defor.bgb.ha,
    defor.bgb.a  = defor.bgb.ha * defor.area.a,
    defor.co2.ha = defor.co2.ha,
    defor.co2.a  = defor.co2.ha * defor.area.a,
    regen.area   = regen.area, 
    regen.area.a = regen.area.a,
    regen.agb.ha = regen.agb.ha,
    regen.agb.a  = regen.agb.ha * regen.area.a,
    regen.bgb.ha = regen.bgb.ha,
    regen.bgb.a  = regen.bgb.ha * regen.area.a,
    regen.co2.ha = regen.co2.ha,
    regen.co2.a  = regen.co2.ha * regen.area.a
  ))
}


# COMMENCER LE TRAITEMENT #####################################################

# Analyse des pertes de biomasse et CO2 dà cause de la déforestation
# et les gains à cause de la régénération pour Togo et les régions ------------

# traitement en parallèle
registerDoParallel(CORES)
res <- foreach(i=0:length(TGO.REG), .combine=rbind) %dopar% {
  
  if(i==0) {
    # Analyse pour l'ensemble du Togo
    data.frame(reg = "TGO", 
               emissions(defor, defor.agb, defor.bgb, regen, regen.agb, regen.bgb))
  } else {
    # Analyse pour une région
    region <- TGO.reg[i,]
    data.frame(reg = region$NAME_1, 
               emissions(defor, defor.agb, defor.bgb, regen, regen.agb, regen.bgb, 
                         aoi=region))
  }
}

# sauvegarder le tableau des résultats
write.csv(res, paste0(AGB.RES.DIR, "/NERF-Results.csv"), row.names=FALSE)







