####################################################################
# NERF_Togo/AGB/4_clean-agb-maps.R: clean AGB raw maps
# ------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 10 October 2019

### DEFINITIONS ############################################################

PERIOD    <- 1985:2019

path <- "p193"


# LOESS smoothing of time series

loess.filter <- function (agb, years, pred, span) {
  ifelse( all(is.na(agb)), 
          return(rep(NA, length(pred))), 
          return(predict(loess(agb ~ years,
                               degree = 1,
                               span = span), 
                         pred)
          ))
}

clean.agb <- function(path) {
  
  # Preparation -------------------------------------------------------------
  
  maps <- stack(dir(paste0(BIOMASS.DIR, "/3_raw-maps/", path), pattern=".*[[:digit:]]{4}\\_mr\\.tif$", full.names=TRUE))
  map.names <- sub("r$", "", names(maps))
  map.cols  <- sub(paste0(path, "\\_"), "X", sub("\\_[[:alnum:]]+$", "M", map.names))
  map.years <- as.numeric(substr(map.cols, 2, 5))
  
  maps.values <- values(maps)                              # extract vector with cell values (huge matrix, takes some time)
  colnames(maps.values) <- map.cols
  
  nsubsets <- numCores - 1                                 # define subsets for parallel processing
  subsets  <- c(0, floor((1:nsubsets)*(nrow(maps.values)/nsubsets)))    
  
  
  
  # Parallel cleaning of pixel trajectories ###############################
  
  registerDoParallel(.env$numCores-1)
  loess.dat <- foreach(i=1:nsubsets, .combine=rbind) %dopar% {
    
    val <- maps.values[(subsets[i]+1):subsets[i+1], ]       # get one tile of the matrix
    
    t(apply(val, 1, loess.filter, years=map.years, pred=map.years, span=2))
  }
  
  
  maps.loess <- maps
  maps.loess[] <- loess.dat
  writeRaster(maps.loess, paste0(BIOMASS.DIR, "/4_clean-maps/", path, "_loess.tif"), overwrite=TRUE)
  
}
    
    
    





################################

dates <- c(1987, 1991, 2000, 2003, 2005, 2007, 2008, 2009, 2010, 2012, 2013, 2015, 2017, 2018, 2019)


for(year in dates) {
  print(paste("Creating AGB map", year, "..."))
  image <- dir(path="../results/images/", pattern=paste0("^", year, ".*_L2.tif$"), full.names=TRUE)
  base  <- paste0("../results/carbone/AGB_", year, "_2015")
  biomass.map       <- create.biomass.map(ref =      raster("../results/carbone/AGB_2015_REF_corr.tif"), 
                                          x =        brick(image))
  writeRaster(biomass.map$map,            filename = paste0(base, ".tif"), overwrite=TRUE)
  writeRaster(biomass.map$map.bias.corr,  filename = paste0(base, "_corr.tif"), overwrite=TRUE)
  sink(                                       file = paste0(base, "_models.txt"))
  print(biomass.map$rf); cat("\n")
  print(importance(biomass.map$rf)); cat("\n")
  print(summary(biomass.map$lm.bias.corr))
  sink()
}



# Show the effect of bias correction

ggplot(biomass.map$sample.pts@data, aes(y=ref.agb, x=agb)) +
  geom_point(alpha=0.01) +
  geom_abline(slope=1, linetype=3, alpha=0.9) +
  geom_abline(intercept = biomass.map$lm.bias.corr$coefficients[1], slope = biomass.map$lm.bias.corr$coefficients[2], linetype=2, alpha=0.9) +
  theme_light()



# load the raw carbon maps and do some analysis
# ≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠

files <- dir(path="../results/carbone", pattern="^.*2015.tif$", full.names=TRUE)
maps <- stack(files)
names(maps) <- dates

maps.dat <- as.data.frame(na.omit(maps[]))

dens.plot <- maps.dat %>% gather(variable, value, X1987, X1991, X2000, X2003, X2008, X2013, X2015, X2018) %>%
  ggplot(aes(x=value, color=variable, linetype=variable)) +
  geom_density() +
  theme_light() 

print(dens.plot)

pdf("../results/carbone/densities_all.pdf")
print(dens.plot)
dev.off()






# Validation of loess-smoothed 2015
map.2015 <- stack("../results/carbone/level2/AGBmaps_loess.tif")[[12]]
map.2003 <- stack("../results/carbone/level2/AGBmaps_loess.tif")[[4]]
plots <- raster::extract(map.2015, plots, sp=TRUE, fun=mean, buffer=20, normalizeWeights=TRUE)
plots <- raster::extract(map.2003, plots, sp=TRUE, fun=mean, buffer=20, normalizeWeights=TRUE)
names(plots)

tmp <- plots@data[,c("BaHa", "AGBmaps_loess.12", "AGBmaps_loess.4")]
tmp$cex <- 1+abs(tmp$AGBmaps_loess.12-tmp$AGBmaps_loess.4)/30
tmp$col[tmp$AGBmaps_loess.12-tmp$AGBmaps_loess.4 <= 0] <- "red"
tmp$col[tmp$AGBmaps_loess.12-tmp$AGBmaps_loess.4 > 0] <- "blue"

pdf("../results/carbone/figures/Validation_LoessMap2015.pdf")
plot(tmp$BaHa, tmp$AGBmaps_loess.12, main="AGB Carbon Stocks on Permanent Inventory Plots",
     xlab="National Forest Inventory 2015 (tC/ha)", ylab="LOESS-smoothed Carbon Map 2015")
points(tmp$BaHa[tmp$cex > 1.7], tmp$AGBmaps_loess.12[tmp$cex > 1.7], cex=tmp$cex[tmp$cex > 1.7], col=tmp$col[tmp$cex > 1.7])
abline(0,1, lty=2)
abline(0,0.5, lty=3)
abline(0,1/0.5, lty=3)
legend("bottomright", col=c("red", "blue"), pch=1, legend=c("Loss of >= 20 tC since 2003", "Gain of >= 20 tC since 2003"))
dev.off()

pdf("../results/carbone/figures/Validation_LoessMap2015_Residuals.pdf")
plot(tmp$AGBmaps_loess.12-tmp$BaHa~tmp$BaHa, main="AGB Carbon Stocks on Permanent Inventory Plots",
     xlab="National Forest Inventory 2015 (tC/ha)", ylab="Residuals of LOESS-smoothed Carbon Map 2015")
points(tmp$BaHa[tmp$cex > 1.7], tmp$AGBmaps_loess.12[tmp$cex > 1.7]-tmp$BaHa[tmp$cex > 1.7], cex=tmp$cex[tmp$cex > 1.7], col=tmp$col[tmp$cex > 1.7])
bias.corr <- lm(tmp$AGBmaps_loess.12-tmp$BaHa~tmp$BaHa)
abline(0,0)
abline(bias.corr, lty=2)
legend("topright", col=c("red", "blue"), pch=1, legend=c("Loss of >= 20 tC since 2003", "Gain of >= 20 tC since 2003"))
dev.off()




# Load the cleaned maps and create tables
# =======================================

biomass.tabs <- function(map, aoi=NULL) {
  
  if(!is.null(aoi)) map <- crop(map, aoi)
  
  agb.dat        <- map[]
  agb.change.dat <- agb.dat[,2:ncol(agb.dat)] - agb.dat[,1:ncol(agb.dat)-1]
  agb.loss.dat   <- agb.change.dat; agb.loss.dat[agb.loss.dat > 0] <- 0
  agb.gain.dat   <- agb.change.dat; agb.gain.dat[agb.gain.dat < 0] <- 0
  
  dates <- as.numeric(substr(names(map), 2, 5))
  
  agb.tab <- data.frame(year      = dates, 
                        total.agb = colSums(agb.dat, na.rm=TRUE) * 30^2/10000, 
                        period    = c(NA, dates[-1] - dates[-length(dates)]), 
                        center    = c(NA, (dates[-1] + dates[-length(dates)])/2))
  
  agb.tab$loss.t.yr              <- c(NA, colSums(agb.loss.dat, na.rm=TRUE) * 30^2/10000)/agb.tab$period
  agb.tab$gain.t.yr              <- c(NA, colSums(agb.gain.dat, na.rm=TRUE) * 30^2/10000)/agb.tab$period
  agb.tab$netchange.t.yr         <- agb.tab$loss.t.yr + agb.tab$gain.t.yr 
  
  agb.tab$loss.pc.yr             <- round(100 * c(NA, 1/agb.tab$period[2:(nrow(agb.tab))] * log((agb.tab$total.agb[1:nrow(agb.tab)-1] + agb.tab$period[2:(nrow(agb.tab))] * agb.tab$loss.t.yr[2:(nrow(agb.tab))]) / agb.tab$total.agb[1:nrow(agb.tab)-1])), 3)
  agb.tab$gain.pc.yr             <- round(100 * c(NA, 1/agb.tab$period[2:(nrow(agb.tab))] * log((agb.tab$total.agb[1:nrow(agb.tab)-1] + agb.tab$period[2:(nrow(agb.tab))] * agb.tab$gain.t.yr[2:(nrow(agb.tab))]) / agb.tab$total.agb[1:nrow(agb.tab)-1])), 3)
  agb.tab$netchange.pc.yr        <- round(100 * c(NA, 1/agb.tab$period[2:(nrow(agb.tab))] * log( agb.tab$total.agb[2:nrow(agb.tab)]                                                                               / agb.tab$total.agb[1:nrow(agb.tab)-1])), 3)
  
  return(agb.tab)
}



plot.biomass <- function(fc, zone, filename=NULL) {
  
  if(zone=="Zone IV") {
    title    <- "Net and Gross Changes in Forest Biomass 1987 - 2019 in Ecological Zone IV"
    ylim     <- c(-750000, 750000)
    ybreaks  <- seq(-750000, 750000, 100000)
    ymbreaks <- seq(-750000, 750000, 50000)
  } 
  
  if (zone=="AOI.1") {
    title    <- "Net and Gross Changes in Forest Biomass 1987 - 2019 in Plaine de Litimé (AOI.1)"
    ylim     <- c(-18000, 15000)
    ybreaks  <- seq(-18000, 15000, 2000)
    ymbreaks <- seq(-18000, 15000, 1000)
  } 
  
  if (zone=="AOI.2") {
    title    <- "Net and Gross Changes in Forest Biomass 1987 - 2019 in Kpélé (AOI.2)"
    ylim     <- c(-18000, 15000)
    ybreaks  <- seq(-18000, 15000, 2000)
    ymbreaks <- seq(-18000, 15000, 1000)
  } 
  
  if(!is.null(filename)) pdf(filename)
  
  fc$netloss.t.yr <- fc$netchange.t.yr 
  fc$netloss.t.yr[fc$netloss.t.yr > 0] <- 0
  fc$netgain.t.yr <- fc$netchange.t.yr 
  fc$netgain.t.yr[fc$netgain.t.yr < 0] <- 0
  
  print(
    fc[-1,] %>% gather(variable, value, loss.t.yr, gain.t.yr, netloss.t.yr, netgain.t.yr) %>%
      ggplot(aes(y=value, x = center, width=period-0.7)) +
      geom_bar(data=. %>% filter(variable %in% c("loss.t.yr", "gain.t.yr")), aes(fill=variable, alpha=0.9), stat = "identity", show.legend=FALSE) +
      guides(alpha = FALSE) +
      geom_bar(data=. %>% filter(variable %in% c("netloss.t.yr", "netgain.t.yr")), aes(fill=variable), stat = "identity") +
      scale_fill_manual(name = NULL, breaks = c("gain.t.yr", "loss.t.yr"), values=c("#00BFC4", "#F8766D", "#00BFC4", "#F8766D"), labels=c("Biomass gain    ", "Biomass loss    ")) +
      xlab("Year") + ylab("Tonnes AGB per Year") + ggtitle(title) +
      scale_x_continuous(minor_breaks = NULL, breaks=c(1987, 2000, 2005, 2010, 2015, 2019)) + 
      scale_y_continuous(breaks=ybreaks, minor_breaks = ymbreaks) + coord_cartesian(ylim=ylim) + 
      theme_light() + theme(legend.position=c(0,1), legend.box = "horizontal", legend.justification=c(-0.2,1.2))
  )
  
  if(!is.null(filename)) dev.off()
}


maps <- brick("../results/carbone/level2/AGBmaps_loess.tif")
names(maps) <- dates

maps.dat <- as.data.frame(na.omit(maps[]))

dens.plot <- maps.dat %>% gather(variable, value, X1991, X2000, X2005, X2010, X2015, X2018) %>%
  ggplot(aes(x=value, color=variable, linetype=variable)) +
  geom_density() +
  theme_light() 

print(dens.plot)

pdf("../results/carbone/figures/densities.pdf")
print(dens.plot)
dev.off()

res.all      <- biomass.tabs(map=maps)
res          <- biomass.tabs(map=maps[[c("X1991", "X2000", "X2005", "X2010", "X2015", "X2018")]])
res.aoi1.all <- biomass.tabs(map=maps, AOI.1)
res.aoi1     <- biomass.tabs(map=maps[[c("X1991", "X2000", "X2005", "X2010", "X2015", "X2018")]], AOI.1)
res.aoi2.all <- biomass.tabs(map=maps, AOI.2)
res.aoi2     <- biomass.tabs(map=maps[[c("X1991", "X2000", "X2005", "X2010", "X2015", "X2018")]], AOI.2)

res.tmp          <- biomass.tabs(map=maps[[c("X1991", "X2003", "X2018")]])

plot.biomass(res,      "Zone IV", "../results/carbone/figures/biomass.pdf")
plot.biomass(res.aoi1, "AOI.1",    "../results/carbone/figures/biomass_aoi1.pdf")
plot.biomass(res.aoi2, "AOI.2",    "../results/carbone/figures/biomass_aoi2.pdf")


write.xlsx(list("Total" = res,
                "AOI.1" = res.aoi1,
                "AOI.2" = res.aoi2),
           file = "../results/carbone/figures/biomass-change.xlsx")




Y.diff <- raster("../results/carbone/AGB_2018_2015_corr.tif") - Y.corr
Y.gain <- Y.diff; Y.gain[Y.gain<=0] <- 0
Y.loss <- Y.diff; Y.loss[Y.loss>=0] <- 0
Y.diff <- writeRaster(round(Y.diff), "../results/carbone/AGB_2018_1991_diff.tif", datatype="INT2S", overwrite=TRUE)

-sum(Y.loss[], na.rm=TRUE)/sum(raster("../results/carbone/AGB_2003_2015_corr.tif")[], na.rm=TRUE)/15
sum(Y.gain[], na.rm=TRUE)/sum(raster("../results/carbone/AGB_2003_2015_corr.tif")[], na.rm=TRUE)/15


dat.c <- data.frame(X1991=raster("../results/carbone/AGB_1991_2015_corr.tif")[],
                    X2003=raster("../results/carbone/AGB_2003_2015_corr.tif")[],
                    X2018=raster("../results/carbone/AGB_2018_2015_corr.tif")[])

dat.c <- na.omit(dat.c)

dat.c %>% gather(variable, value, X1991:X2018) %>% 
  ggplot(aes(y=value, x=variable, color=variable)) +
  geom_boxplot() +
  coord_flip() +
  theme_light() 

pdf("../results/carbone/densities_91-03-18.pdf")

dat.c %>% gather(variable, value, X1991:X2018) %>%
  ggplot(aes(x=value, color=variable)) +
  geom_density() +
  theme_light() 

dev.off()





# verify the map

val.points <- raster::extract(Y.corr, dat.placette, sp=TRUE, fun=mean, buffer=20, normalizeWeights=TRUE)
val.points <- raster::extract(Y.2003.corr, val.points, sp=TRUE, fun=mean, buffer=20, normalizeWeights=TRUE)
val.points <- raster::extract(Y.2015, val.points, sp=TRUE, fun=mean, buffer=20, normalizeWeights=TRUE)
plot(val.points$BaHa, val.points$AGB_2015_2003_corr)
abline(0,1)
mod <- lm(AGB_2015_2003_corr ~ BaHa, data=val.points)
abline(mod, lty=2)




############ compare 2015 carbon map with the forest change map #######################

map.carb <- raster("../results/carbone/AGB_2015_REF.tif")
map.gain <- raster("../results/fcc/3_alldates-direct2003/level3/forest-gain.tif")
map.loss <- raster("../results/fcc/3_alldates-direct2003/level3/forest-loss.tif")

map.loss.a2015 <- map.loss; map.loss.a2015[map.loss.a2015<2015] <- NA
map.carb.loss.a2015 <- mask(map.carb, map.loss.a2015)
hist(map.carb.loss.a2015[]); summary(map.carb.loss.a2015[]); sum(map.carb.loss.a2015[], na.rm=TRUE)*(30*30/10000)

map.gain.b2015 <- map.gain; map.gain.b2015[map.gain.b2015>=2015] <- NA
map.gain.b2015[] <- 2015 - map.gain.b2015[]
map.carb.gain.b2015 <- mask(map.carb, map.gain.b2015)
hist(map.carb.gain.b2015[]); summary(map.carb.gain.b2015[]); sum(map.carb.gain.b2015[], na.rm=TRUE)*(30*30/10000)

boxplot(map.carb.gain.b2015[] ~ map.gain.b2015[])

map.carb.loss.a2015.2 <- mask(Y.loss, map.loss.a2015)






# ######### Try with linear regression ############
# X$evi <- evi(nir=X$NIR, r=X$R, b=X$B)
# X$ndvi <- nd(X$NIR, X$R)
# X$ndmi <- nd(X$NIR, X$SWIR1)
# X$nbr <- nd(X$NIR, X$SWIR2)
# cor(-X$SWIR1[], Y.2015[], use="pairwise.complete.obs")
# s <- sample(1:ncell(Y.2015), 10000)
# plot(log(Y.2015[s]) ~ log(X$SWIR1[s]))
# 
# 
# val.points <- raster::extract(X, val.points, sp=TRUE, fun=mean, buffer=20, normalizeWeights=TRUE)
# plot(log(val.points$BaHa) ~ log(val.points$SWIR1))
# val.points$BaHa.log <- log(val.points$BaHa)
# val.points$SWIR1.log <- log(val.points$SWIR1)
# lm.swir <- lm(BaHa.log ~ SWIR1.log, data=val.points@data)
# abline(lm.swir)
# 
# plot(val.points$BaHa ~ val.points$SWIR1)
# lines(seq(1000, 2500, 50), exp(lm.swir$coeff[1] + lm.swir$coeff[2]*log(seq(1000, 2500, 50))), lty=4)
# 
# Y.2015.swir <- exp(lm.swir$coeff[1] + lm.swir$coeff[2]*log(X$SWIR1))
# val.points <- raster::extract(Y.2015.swir, val.points, sp=TRUE, fun=mean, buffer=20, normalizeWeights=TRUE)
# 
# plot(val.points$BaHa, val.points@data[,25])
# abline(0,1)
# 
# 
# 
# Y.2015.log <- log(Y.2015)
# SWIR1.log <- log(X$SWIR1)
# lm.swir.map <- lm(Y.2015.log[] ~ SWIR1.log[])
# abline(lm.swir.map, lty=2)
# 
# plot(Y.2015[s] ~ X$SWIR1[s])
# lines(seq(1000, 4000, 50), exp(lm.swir.map$coeff[1] + lm.swir.map$coeff[2]*log(seq(1000, 4000, 50))), lty=4, col="red")
# Y.2015.swir <- exp(lm.swir.map$coeff[1] + lm.swir.map$coeff[2]*log(X$SWIR1))
# Y.2015.swir <- writeRaster(Y.2015.swir, "../results/carbone/AGB_2015_SWIR1.tif")
# 
# 
# Y <- Y.2015.swir
# X <- stack("../results/images/20030127_L7_B123457_L2.tif")
# names(X) <- c("B", "G", "R", "NIR", "SWIR1", "SWIR2")
# 
# Y.log <- log(Y)
# SWIR.log <- log(X$SWIR1)
# 
# lm.swir <- lm(Y.log[] ~ SWIR.log[])
# Y.swir <- exp(lm.swir$coeff[1] + lm.swir$coeff[2]*SWIR.log)
# Y.swir <- writeRaster(Y.swir, "../results/carbone/AGB_2003_SWIR1.tif")
# 
# Y.diff <- Y.2015.swir - Y.swir
# Y.diff <- writeRaster(Y.diff, "../results/carbone/AGB_diff_SWIR1.tif")









ggplot(plots@data[], aes(x=cov18, y=BaHa)) + 
  geom_point(aes(colour=OccSol),               # colour depends on cond2
             size=3)    +
  facet_wrap( ~ OccSol)





