###############################################################################
# 01_compile-AGB.R: Évaluation de la biomasse des données de l'IFN
# -----------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 13 Mai 2020

# Définitions des variables ===================================================

PLOT.SIZE <- 20^2*pi  # Taille de la parcelle en mètres carrés

# Rapports racines-tige forêts tropicales sèches selon Mokany et al. (2006)
# https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2486.2005.001043.x
RSR.Y    <- 0.563     # moyenne pour la jeune forêt      
RSR.Y.SE <- 0.086     # écart type pour la jeune forêt  
RSR.O    <- 0.275     # moyenne pour les forêts matures 
RSR.O.SE <- 0.003     # écart type pour les forêts matures
RSR.AGB  <- 20        # Seuil biomasse aérienne jeunes forêts <-> forêts matures


# Lire et préparer les données IFN-1 ==========================================

plots <- read.xlsx(paste0(DIR.RAW.DAT, "/IFN/IFN-Togo-2015.xlsx"), "placettes")[,1:4]
names(plots) <- c("PlotID", "X", "Y", "LULC")

trees <- read.xlsx("~/Downloads/IFN-Togo-2015.xlsx", "arbres")[,c(1,4:6,8)]
names(trees) <- c("PlotID", "Species", "Status", "DBH", "H")

# Tableau avec les espèces et le nombre d'observations
species <- as.data.frame(table(trees$Species))
names(species) <- c("Species", "Count")

# fusionner avec les densités de bois basales ---------------------------------
# Densité basale: rapport entre mass anhydre et volume vert (à l'état saturé)
# Source des données: étude biomasse de Fonton et al. 2018
fonton <- read.csv2(paste0(DIR.RAW.DAT, "/IFN/donnees_Tg_Fonton.csv"), encoding="latin1")[,c(4,7)]
fonton <- aggregate(list(D=fonton$wsg), by=list(Species=fonton$NOM), FUN=modal)
species <- merge(species[,c("Species", "Count")], 
                     fonton[,c("Species", "D")], by="Species")
species$Source <- "Fonton"

# fusionner avec les arbres dans l'IFN-1
trees <- merge(trees, species, by="Species", all.x=TRUE)

# estimer la biomasse aérienne des arbres  ------------------------------------
# en utilisant la fonction allométrique de Chave et al. (2014)
# https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.12629
trees$AGB <- 0.0673 * (trees$D * trees$DBH^2 * trees$H)^0.976

# distinguer la biomasse des arbres vivants et celle des arbres morts
trees$AGBm <- trees$AGBv <- trees$AGB    # Dupliquer les valeurs de la biomasse ...
trees$AGBv[trees$Status != "V"] <- 0     # ... et mettre à 0 aux endroits appropriés
trees$AGBm[trees$Status == "V"] <- 0       

# biomasse aérienne et bois mort par parcelle (somme des arbres) --------------
plots <- merge(plots, 
               aggregate(trees[,c("AGBv", "AGBm")], 
                         by=list(PlotID=trees$PlotID), 
                         # somme dbiomasse es arbres -> à l'héctare -> en tonne
                         FUN=function(x) sum(x) * 10000 / PLOT.SIZE / 1000),  
               by="PlotID", all.x=TRUE)  
# mettre à 0 la biomasse et le bois mort pour les parcelles sans valeurs (NA)
plots$AGBv[is.na(plots$AGBv)] <- 0                                            
plots$AGBm[is.na(plots$AGBm)] <- 0
plots$AGB <- plots$AGBv + plots$AGBm


# estimer la biomasse racinaire par parcelle ----------------------------------
# avec les facteurs root-shoot de Mokany et al. (2006) 

# jeunes forêts
plots$BGB[plots$AGBv <= RSR.AGB] <- plots$AGBv[plots$AGBv <= RSR.AGB] * RSR.Y
# forêts matures
plots$BGB[plots$AGBv  > RSR.AGB] <- plots$AGBv[plots$AGBv  > RSR.AGB] * RSR.O 

# Biomasse totale = biomasse aérienne (vivant et mort) + biomasse racinaire
plots$BM <- plots$AGB + plots$BGB


# Sauvegarder le résultat et production des figures et tableaux ===============
# Note: sur Mac utiliser fileEncoding = "macintosh" 

write.csv(plots, paste0(DIR.MRV.AGB.REF, "/IFN-plots.csv"), 
          row.names = FALSE) #, fileEncoding = "macintosh")


# distribution des biomasses par strate IFN ------------------------------------

# boxplot biomasse aérienne
pdf(paste0(DIR.MRV.AGB.REF, "/AGB-vs-LULC.pdf"))
par(mar=c(11,5,1,1), cex.axis=0.7)
boxplot(plots$AGBv~plots$LULC, las=2, ylab="AGB (t/ha)", xlab=NULL)
dev.off()

# boxplot bois mort
pdf(paste0(AGB.REF.DIR, "/AGBm-vs-LULC.pdf"))
par(mar=c(11,5,1,1), cex.axis=0.7)
boxplot(plots$AGBm~plots$LULC, las=2, ylab="Bmort (t/ha)", xlab=NULL)
dev.off()

# tableau biomass per LU/LC category (moyenne et écart type)
bm.lulc.tab <- plots %>%
  group_by(LULC) %>%            # grouper par strate
  summarise(n=length(AGB),      # definir colonnes et calcul des valeurs
            AGBv.mean=mean(AGBv), AGBv.sd=sd(AGBv),     
            BGB.mean=mean(BGB), BGB.sd=sd(BGB),
            AGBm.mean=mean(AGBm), AGBm.sd=sd(AGBm), 
            BM.mean =mean(BM),  BM.sd =sd(BM))
write.csv(bm.lulc.tab, paste0(AGB.REF.DIR, "/AGB_LULC.csv"), 
          row.names = FALSE) #, fileEncoding = "macintosh")

# différences biomasse par type de conversion LU/LC ---------------------------
dbm.lulc.tab <- bm.lulc.tab %>% 
  expand(FROM = nesting(LULC, AGB.mean, AGB.sd, BGB.mean, BGB.sd, BM.mean, BM.sd), 
         TO   = nesting(LULC, AGB.mean, AGB.sd, BGB.mean, BGB.sd, BM.mean, BM.sd)) %>%
  mutate(from = FROM$LULC, to = TO$LULC,
         dAGB.mean = TO$AGB.mean - FROM$AGB.mean, dAGB.sd=sqrt(FROM$AGB.sd^2 + TO$AGB.sd^2),
         dBGB.mean = TO$BGB.mean - FROM$BGB.mean, dBGB.sd=sqrt(FROM$BGB.sd^2 + TO$BGB.sd^2),
         dBM.mean  = TO$BM.mean  - FROM$BM.mean,  dBM.sd =sqrt(FROM$BM.sd^2  + TO$BM.sd^2)) %>%
  select(c(from, to, dAGB.mean, dAGB.sd, dBGB.mean, dBGB.sd, dBM.mean, dBM.sd))
write.csv(dbm.lulc.tab, paste0(AGB.REF.DIR, "/AGB_LULC-diff.csv"), 
          row.names = FALSE) #, fileEncoding = "macintosh")
