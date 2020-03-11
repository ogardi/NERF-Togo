####################################################################
# NERF_Togo/AGB/5_analyze-AGB.R: analyze AGB maps
# ------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 10 October 2019

RSR.Y    <- 0.563
RSR.Y.SE <- 0.086
RSR.O    <- 0.275
RSR.O.SE <- 0.003
RSR.AGB  <- 20
C.RATIO  <- 0.47


# evaluate 2015 biomass map with IFN strata ------------------

agb.2015 <- raster(paste0(AGB.REF.DIR, "/TGO_2015_AGB_R.tif"))
strata  <- raster(paste0(DATA.DIR, "/RapidEye/TGO_30m.tif"))

zonal(agb.2015, strata, fun="mean")
zonal(agb.2015, strata, fun="sd")

# AGB density plots ---------

files <- dir(AGB.REF.DIR, pattern="^TGO.*\\.tif$", full.names=TRUE)
maps <- stack(files)
names(maps) <- substr(names(maps), 5, 8)

maps.dat <- as.data.frame(na.omit(maps[[c("X2003", "X2018")]][]))

dens.plot <- maps.dat %>% gather(variable, value, X2003, X2018) %>%
  ggplot(aes(x=value, color=variable, linetype=variable)) +
  geom_density() +
  theme_light() 

# print(dens.plot)

pdf(paste0(AGB.RES.DIR, "/AGB-densities_03-18.pdf"))
print(dens.plot)
dev.off() 


# Analyze emissions and enhancements 2003/2018 ----------------

agb.2003 <- raster(paste0(AGB.REF.DIR, "/TGO_2003_AGB_R.tif"))
bgb.2003 <- agb.2003 * RSR.O
bgb.2003[agb.2003 <= RSR.AGB] <- agb.2003[agb.2003 <= RSR.AGB] * RSR.Y

agb.2018 <- raster(paste0(AGB.REF.DIR, "/TGO_2018_AGB_R.tif"))
bgb.2018 <- agb.2018 * RSR.O
bgb.2018[agb.2018 <= RSR.AGB] <- agb.2018[agb.2018 <= RSR.AGB] * RSR.Y

fc.2003 <- raster(paste0(FCC.CLN.DIR, "/FC30/TGO/TGO_2003_F30cf.tif"))
fc.2018 <- raster(paste0(FCC.CLN.DIR, "/FC30/TGO/TGO_2018_F30cf.tif"))





# emissions from deforestation

defor <- fc.2003 == 1 & fc.2018 != 1

nf.agb <- mean(agb.2018[fc.2018 == 3], na.rm=TRUE)   # average non-forest AGB
pd.agb <- mean(agb.2018[defor], na.rm = TRUE)        # average post-defor AGB
defor.agb <- mask(agb.2003, defor, maskvalue=0) - nf.agb
defor.agb[defor.agb < 0] <- 0

nf.bgb <- mean(bgb.2018[fc.2018 == 3], na.rm=TRUE)   # average non-forest BGB
pd.bgb <- mean(bgb.2018[defor], na.rm = TRUE)        # average post-defor BGB
defor.bgb <- mask(bgb.2003, defor, maskvalue=0) - nf.bgb
defor.bgb[defor.bgb < 0] <- 0 

# enhancements from reforestation

regen <- fc.2003 != 1 & fc.2018 == 1

regen.agb.2003 <- mask(agb.2003, regen, maskvalue=0)
regen.agb.2018 <- mask(agb.2018, regen, maskvalue=0)
regen.agb <- regen.agb.2018 - regen.agb.2003
regen.agb[regen.agb < 0] <- 0

regen.agb.2003 <- mask(bgb.2003, regen, maskvalue=0)
regen.agb.2018 <- mask(bgb.2018, regen, maskvalue=0)
regen.bgb <- regen.bgb.2018 - regen.bgb.2003
regen.bgb[regen.bgb < 0] <- 0



# statistics for TGO and regions

emissions <- function(defor, defor.agb, defor.bgb, 
                      regen, regen.agb, regen.bgb, aoi=NULL, years=15) {
  
  if(!is.null(aoi)) {
    defor     <- mask(crop(defor, aoi), aoi)
    defor.agb <- mask(crop(defor.agb, aoi), aoi)
    defor.bgb <- mask(crop(defor.agb, aoi), aoi)
    regen     <- mask(crop(regen, aoi), aoi)
    regen.agb <- mask(crop(regen.agb, aoi), aoi)
    regen.bgb <- mask(crop(regen.bgb, aoi), aoi)
  }
  
  defor.co2    <- (defor.agb + defor.bgb) * C.RATIO * 44/12
  defor.area   <- sum(defor[], na.rm=TRUE) * 30^2 / 10000
  defor.area.a <- defor.area / years
  defor.agb.ha <- mean(defor.agb[], na.rm=TRUE)
  defor.bgb.ha <- mean(defor.bgb[], na.rm=TRUE)
  defor.co2.ha <- mean(defor.co2[], na.rm=TRUE)
  
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



# DO THE WORK

registerDoParallel(.env$numCores-1)
res <- foreach(i=0:length(TGO.reg), .combine=rbind) %dopar% {
  
  if(i==0) {
    data.frame(reg = "TGO", emissions(defor, defor.agb, defor.bgb, regen, regen.agb, regen.bgb))
  } else {
    region <- TGO.reg[i,]
    data.frame(reg = region$NAME_1, emissions(defor, defor.agb, defor.bgb, regen, regen.agb, regen.bgb, aoi=region))
  }
}

write.csv(res, paste0(AGB.RES.DIR, "/NERF-Results.csv"), row.names=FALSE)








# Aggregated uncertainty (Monte Carlo Simulation)

# activity data (uncertainty)

mc.defor  <- rnorm(1000, mean.defor, sd)
mc.regen  <- rnorm(1000, mean.regen, sd)


# emission factor
mean.agb.2003 <- 
mean.agb.2018 <- mean(agb.2018[defor])

mc.agb.2003   <- rnorm(1000, mean(agb.2003[defor]), sqrt((RMSE^2)/n) ) 
mc.bgb.2003   <- mc.agb.2003 * rnorm(1000, mean(bgb.2003), RSR.O.SE)

mc.agb.2018   <- rnorm(1000, mean(agb.2018[defor]), sqrt((RMSE^2)/n) ) 
mc.bgb.2018   <- mc.agb.2003 * rnorm(1000, mean(bgb.2003), RSR.O.SE)






