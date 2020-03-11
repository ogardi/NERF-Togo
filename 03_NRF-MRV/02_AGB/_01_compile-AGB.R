####################################################################
# NERF_Togo/AGB/2_compile-IFN.R: evaluating AGB of IFN plots
# ------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 3 December 2019

PLOT.SIZE <- 20^2*pi          # Plot size in square meters
RSR.Y    <- 0.563             # Root-shoot ratios according Mokany
RSR.Y.SE <- 0.086
RSR.O    <- 0.275
RSR.O.SE <- 0.003
RSR.AGB  <- 20 

# read inventory data
# -------------------
plots <- read.xlsx(paste0(DATA.DIR, "/IFN/IFN-Togo-2015.xlsx"), "placettes")[,1:4]
names(plots) <- c("PlotID", "X", "Y", "LULC")

trees <- read.xlsx(paste0(DATA.DIR, "/IFN/IFN-Togo-2015.xlsx"), "arbres")[,c(1,4:6,8)]
names(trees) <- c("PlotID", "Species", "Status", "DBH", "H")

# merge with specific wood densities used by from Fonton et al. 2018
species <- as.data.frame(table(trees$Species))
names(species) <- c("Species", "Count")
fonton <- read.csv2(paste0(DATA.DIR, "/IFN/donnees_Tg_Fonton.csv"), encoding="latin1")[,c(4,7)]
fonton <- aggregate(list(D=fonton$wsg), by=list(Species=fonton$NOM), FUN=modal)
species <- merge(species[,c("Species", "Count")], 
                     fonton[,c("Species", "D")], by="Species")
species$Source <- "Fonton"


# merge tree table with species densities
# -------------------------------------------------------------------------
trees <- merge(trees, species, by="Species", all.x=TRUE)

# estimer la biomasse des arbres avec la fonction de Chave et al. (2014)
# ----------------------------------------------------------------------
trees$AGB <- 0.0673 * (trees$D * trees$DBH^2 * trees$H)^0.976
trees$AGBm <- trees$AGBv <- trees$AGB                # copier les valeurs dans une colonne bois mort
trees$AGBv[trees$Status != "V"] <- 0       # mettre ?? zero la biomass l?? o?? l'arbre n'est pas vivant
trees$AGBm[trees$Status == "V"] <- 0       # mettre ?? zero le bois mort l?? o?? l'arbre est vivant

# aggregation de la biomasse est le bois mort par plots (somme)
# ----------------------------------------------------------------
plots <- merge(plots, 
               aggregate(trees[,c("AGBv", "AGBm")], by=list(PlotID=trees$PlotID), 
                         FUN=function(x) sum(x) * 10000 / PLOT.SIZE / 1000),  
               by="PlotID", all.x=TRUE)                        # joindre tmp avec notre tableau plotss
plots$AGBv[is.na(plots$AGBv)] <- 0                                            # mettre ?? zero la biomasse et le bois mort pour les plotss sans valeurs (NA)
plots$AGBm[is.na(plots$AGBm)] <- 0
plots$AGB <- plots$AGBv + plots$AGBm


# estimer la biomasse racinaire par plot avec les facteurs root-shoot de Mokany et al. (2006) pour les for??ts tropicales s??ches
# -----------------------------------------------------------------------------------------------------------------------------
plots$BGB[plots$AGBv <= RSR.AGB] <- plots$AGBv[plots$AGBv <= RSR.AGB] * RSR.Y
plots$BGB[plots$AGBv  > RSR.AGB] <- plots$AGBv[plots$AGBv  > RSR.AGB] * RSR.O 

plots$BM <- plots$AGB + plots$BGB

# write plot coordinates and biomass to file
# ------------------------------------------

write.csv(plots, paste0(AGB.REF.DIR, "/IFN-plots.csv"),
          row.names = FALSE) #, fileEncoding = "macintosh")


# produire le tableau crois?? avec les biomasses par strate (fonctionalit?? de l'extension dplyr)
# ---------------------------------------------------------------------------------------------
pdf(paste0(AGB.REF.DIR, "/AGB-vs-LULC.pdf"))
par(mar=c(11,5,1,1), cex.axis=0.7)
boxplot(plots$AGBv~plots$LULC, las=2, ylab="AGB (t/ha)", xlab=NULL)
dev.off()

pdf(paste0(AGB.REF.DIR, "/AGBm-vs-LULC.pdf"))
par(mar=c(11,5,1,1), cex.axis=0.7)
boxplot(plots$AGBm~plots$LULC, las=2, ylab="Bmort (t/ha)", xlab=NULL)
dev.off()

# biomass per LU/LC category
bm.lulc.tab <- plots %>%
  group_by(LULC) %>%                                                       # grouper par strate
  summarise(n=length(AGB),
            AGBv.mean=mean(AGBv), AGBv.sd=sd(AGBv),     # definir les colonnes et le calcul des valeurs
            BGB.mean=mean(BGB), BGB.sd=sd(BGB),
            AGBm.mean=mean(AGBm), AGBm.sd=sd(AGBm),     # definir les colonnes et le calcul des valeurs
            BM.mean =mean(BM),  BM.sd =sd(BM))

write.csv(bm.lulc.tab, paste0(AGB.REF.DIR, "/AGB_LULC.csv"), 
          row.names = FALSE) #, fileEncoding = "macintosh")

# differences biomass per LU/LC conversion
dbm.lulc.tab <- bm.lulc.tab %>% 
  expand(FROM=nesting(LULC, AGB.mean, AGB.sd, BGB.mean, BGB.sd, BM.mean, BM.sd), 
         TO=nesting(LULC, AGB.mean, AGB.sd, BGB.mean, BGB.sd, BM.mean, BM.sd)) %>%
  mutate(from = FROM$LULC, to = TO$LULC,
         dAGB.mean=TO$AGB.mean - FROM$AGB.mean, dAGB.sd=sqrt(FROM$AGB.sd^2 + TO$AGB.sd^2),
         dBGB.mean=TO$BGB.mean - FROM$BGB.mean, dBGB.sd=sqrt(FROM$BGB.sd^2 + TO$BGB.sd^2),
         dBM.mean =TO$BM.mean  - FROM$BM.mean,  dBM.sd =sqrt(FROM$BM.sd^2  + TO$BM.sd^2)) %>%
  select(c(from, to, dAGB.mean, dAGB.sd, dBGB.mean, dBGB.sd, dBM.mean, dBM.sd))

write.csv(dbm.lulc.tab, paste0(AGB.REF.DIR, "/AGB_LULC-diff.csv"), 
          row.names = FALSE) #, fileEncoding = "macintosh")


