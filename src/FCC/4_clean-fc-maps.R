##########################################################################
# NERF_Togo/FCC/4_clean-fc-maps_orig.R: clean time-serie of raw forest cover maps
# ------------------------------------------------------------------------
# Bern University of Applied Sciences
# Oliver Gardi, <oliver.gardi@bfh.ch>
# 20 May 2019


# Default parameters -------------------------------------------------------

COV.FC         <- 30

# Function for temporal cleaning of path ------------------------------------

clean.temporal <- function(path) {

  # Preparation -------------------------------------------------------------
  
  maps <- stack(dir(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/", path), pattern=".*[[:digit:]]{4}\\_F.*\\.tif$", full.names=TRUE))
  map.names <- sub("r$", "", names(maps))
  map.cols  <- sub(paste0(path, "\\_"), "X", sub("\\_[[:alnum:]]+$", "M", map.names))
  
  no.map.cols <- paste0("X", PERIOD[!PERIOD %in% gsub("[[:alpha:]]", "", map.cols)], "_")
  col.order <- c(map.cols, no.map.cols)[order(c(map.cols, no.map.cols))]
  
  maps.values <- values(maps)                              # extract vector with cell values (huge matrix, takes some time)
  colnames(maps.values) <- map.cols
  
  nsubsets <- numCores - 1                                 # define subsets for parallel processing
  subsets  <- c(0, floor((1:nsubsets)*(nrow(maps.values)/nsubsets)))    
  
  
  
  # Parallel cleaning of pixel trajectories ###############################
  
  registerDoParallel(.env$numCores-1)
  maps.values.clean <- foreach(i=1:nsubsets, .combine=rbind) %dopar% {
    
    # 0. get subset from maps.values
    val <- maps.values[(subsets[i]+1):subsets[i+1], ]       # get one tile of the matrix
    val[!(is.na(val) | val %in% c(1,3))] <- NA              # set everything to NA that is not NA already or 1,3
    
    
    # 1. remove isolated NAs
    str   <- apply(val, 1, paste, collapse="")              # convert to strings
    str.c <- gsub("NA", "9", str)                           # replace NA with 9
    while(!identical(str, str.c)) {                         # clean until convergence
      str <- str.c
      str.c <- gsub("^(.*3)9(9*3.*)$", "\\13\\2", str.c)    # set NA to the corrsponding class
      str.c <- gsub("^(.*1)9(9*1.*)$", "\\11\\2", str.c)                 
    }                                                       # convert back to numeric vector
    val <- matrix(as.numeric(unlist(strsplit(str.c, ""))), ncol=ncol(val), byrow=TRUE)   
    val[val==9] <- NA                                       # and set NAs again
    colnames(val) <- map.cols
    
    # 2. clean trajectories with sliding window
    val.c   <- val
    val.o <- val[] <- 0
    iter    <- 0                                            # clean until convergence
    while(!identical(val, val.c) & !identical(val.o, val.c)) {   
      iter <- iter+1
      message("    -Clean modal: iteration ", iter, " ... ", appendLF = FALSE)
      val.o <- val
      val <- val.c
      for(l in 3:(ncol(val)-2)) {                             # 3rd - 3rd last year: modal window size 5 
        val.c[,l] <- apply(val[,(l-2):(l+2)], 1, modal, na.rm=TRUE, ties='NA')
      }
      for(l in c(2, ncol(val)-1)) {                           # 2nd and 2nd last year: modal window size 3 
        val.c[,l] <- apply(val.c[,(l-1):(l+1)], 1, modal, na.rm=TRUE, ties='NA')
      }
      message("done")
    }
    val <- val.c
    
    # 3. final regexp cleaning
    val   <- cbind(val, matrix(nrow=nrow(val),              # add years with no observation               
                               ncol=length(no.map.cols), dimnames=list(NULL, no.map.cols)))
    val   <- val[, col.order]                               # order correctly           
    str   <- apply(val, 1, paste, collapse="")              # convert to strings
    str.c <- gsub("NA", "9", str)                           # replace NA with 9
    
    # replace remaining NAs with preceeding land-cover
    while(!identical(str, str.c)) {                         # clean until convergence
      str <- str.c
      str.c <- gsub("39", "33", str.c)                      # set NA to preceeding class
      str.c <- gsub("19", "11", str.c)                 
    }
    str <- ""
    
    # replace trailing NAs with following land-cover
    while(!identical(str, str.c)) {                         # clean until convergence
      str <- str.c
      str.c <- gsub("93", "33", str.c)                      # set NA to following class
      str.c <- gsub("91", "11", str.c)                 
    }
    str <- ""
    
    # removing isolated 1s (regeneration that is lost again, not)
    while(!identical(str, str.c)) {
      str <- str.c
      str.c <-gsub("^(.*3)1(1{0,9}3.*)$", "\\13\\2", str.c)  # removing series of 1s up to length 10
      #str.c <-gsub("^(.*33)1(1*33.*)$", "\\13\\2", str.c)
    }
    str <- ""
    
    # consider regeneration only as forest after 10 years (before it is considered as potential regeneration 2)
    while(!identical(str, str.c)) {
      str <- str.c
      # irgendetwas -> nichtwald -> (potentielle regen, max 8) -> Wald -> irgendetwas
      str.c <-gsub("^(.*32{0,8})1(.*)$", "\\12\\2", str.c)
    }
    str <- ""
    
    
    val <- matrix(as.numeric(unlist(strsplit(str.c, ""))), ncol=ncol(val), byrow=TRUE) # convert back to numeric vector
    val[val==9] <- NA                                                                  # and set NAs again
    colnames(val) <- col.order
    val[,grepl("M$", colnames(val))]                                                   # and get data for map layers
    
    # in order to be considered as forest in 2003, a pixels needs to be forest as well in 2000 and 1991
    # in order to be considered as "regeneration" in 2018, a pixel needs to be forest from 2008 onwards
    
  }
  
  # Write the result to the disk --------------------------------------------
  # without the "uncleaned" first and the last layer
  values(maps) <- maps.values.clean
  
  writeRaster(dropLayer(maps, c(1,nlayers(maps))), 
              filename=paste0(FCC.CLN.DIR, "/FC", COV.FC, "/", path, "/", map.names[2:(nlayers(maps)-1)], "c.tif"), 
              bylayer=TRUE, format="GTiff", datatype="INT2U", overwrite=TRUE)
  
  
  # write trajectories to textfile
  maps.strings <- apply(maps.values.clean, 1, paste, collapse="")
  sink(paste0(FCC.CLN.DIR, "/FC", COV.FC, "/", path, "/Trajectories.txt"))
  print(table(maps.strings))
  sink()

}











# Functions for spatial cleaning of a set of images -----------------------

# Remove forest/defor/regen pixels that NEVER has been part of "forest > 0.5ha" since 2003
clean.spatial.forest <- function(maps, exclude = NULL, size=6, connectedness=8){
  # default: at least 6 connected pixels (0.54 ha) with 8-connectedness
  
  tmp1 <- tempfile(pattern = "", fileext = ".tif")
  tmp2 <- tempfile(pattern = "", fileext = ".tif")
  
  fcc.map <- apply(maps[[-exclude]][], 1, paste, collapse="")
  onceforest.map <- raster(maps)
  onceforest.map[] <- NA 
  
  # onceforest.map <- fcc.map                             # create "once-forest" map
  onceforest.map[grepl("^3*$",    fcc.map)] <- 3
  onceforest.map[grepl("^.*1.*$", fcc.map)] <- 1
  onceforest.map[grepl("^.*2.*$", fcc.map)] <- 1
  writeRaster(onceforest.map, tmp1)
                                                        # remove isolated forest patches < xy ha
  system(paste0("gdal_sieve.py -st ", size, " -", connectedness, " -nomask ", tmp1, " ", tmp2))
  onceforest.map.clean <- raster(tmp2)
  onceforest.map.clean[onceforest.map.clean == -2147483648] <- NA
                                                        # mask the maps with this cleaned forest map
  maps.clean <- mask(maps, onceforest.map.clean, maskvalue=3, updatevalue=3)
  writeRaster(maps.clean, filename=paste0(FCC.CLN.DIR, "/FC", COV.FC, "/TGO/", names(maps.clean), "f.tif"), 
              bylayer=TRUE, format="GTiff", datatype="INT2U", overwrite=TRUE)
  
}

clean.spatial.fcc <- function(maps, size=3, connectedness=8){
  
  # at least 3 connected pixels (8-connectedness)
  
  # combine the maps 2003, 2005 and 2018 
  # change:            0    "321" "112" "122" "322" "132" "332" "113" "213" "313"  
  # stable non-forest: 3    "333"
  # stable forest:     1    "111"
  
  # filter the map (forest islands and non-forest islands are also filtered out)
  system("/Library/Frameworks/GDAL.framework/Programs/gdal_sieve.py -st 3 -8 -nomask ./results/forest.tif ./results/forest_ha.tif")
  
}




### DO THE WORK ###########################################################

# change extent of p192_2019 to the extent of other p192 images (TODO: to be done already in 1_prepare-images.R)
extend(brick(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p192/p192_2019_F", COV.FC, "r.tif")), raster(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p192/p192_2018_F", COV.FC, "r.tif")), 
       filename=paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p192/p192_2019_F", COV.FC, "rt.tif"), format="GTiff", datatype="INT2U")
file.rename(paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p192/p192_2019_F", COV.FC, "rt.tif"), paste0(FCC.RAW.DIR, "/FC", COV.FC, "/p192/p192_2019_F", COV.FC, "r.tif"))

# Temporal cleaning of paths ----------------------------------------------

for(path in c("p192", "p193", "p194")) {
  clean.temporal(path)
} 



# Merging paths -----------------------------------------------------------

for(year in JNT.YEARS) {
  merge(mask(crop(brick(paste0(FCC.CLN.DIR, "/FC", COV.FC, "/p193/p193_", year, "_F", COV.FC, "c.tif")), TGO), TGO),
        mask(crop(brick(paste0(FCC.CLN.DIR, "/FC", COV.FC, "/p192/p192_", year, "_F", COV.FC, "c.tif")), TGO), TGO),
        mask(crop(brick(paste0(FCC.CLN.DIR, "/FC", COV.FC, "/p194/p194_", year, "_F", COV.FC, "c.tif")), TGO), TGO),
        filename=paste0(FCC.CLN.DIR, "/FC", COV.FC, "/TGO/TGO_", year, "_F", COV.FC, "c.tif"), overwrite=TRUE)
}

# merge(mask(crop(brick(paste0(FCC.CLN.DIR, "/p193/p193_1990_F30c.tif")), TGO), TGO),
#       mask(crop(brick(paste0(FCC.CLN.DIR, "/p192/p192_1991_F30c.tif")), TGO), TGO),
#       mask(crop(brick(paste0(FCC.CLN.DIR, "/p194/p194_1997_F30c.tif")), TGO), TGO),
#       filename=paste0(FCC.CLN.DIR, "/TGO/1_clean/TGO_1991_F30c.tif"), overwrite=TRUE)
# 
# merge(mask(crop(brick(paste0(FCC.CLN.DIR, "/p193/p193_2000_F30c.tif")), TGO), TGO),
#       mask(crop(brick(paste0(FCC.CLN.DIR, "/p192/p192_2001_F30c.tif")), TGO), TGO),
#       mask(crop(brick(paste0(FCC.CLN.DIR, "/p194/p194_2000_F30c.tif")), TGO), TGO),
#       filename=paste0(FCC.CLN.DIR, "/TGO/1_clean/TGO_2000_F30c.tif"), overwrite=TRUE)
# 
# merge(mask(crop(brick(paste0(FCC.CLN.DIR, "/p193/p193_2009_F30c.tif")), TGO), TGO),
#       mask(crop(brick(paste0(FCC.CLN.DIR, "/p192/p192_2011_F30c.tif")), TGO), TGO),
#       mask(crop(brick(paste0(FCC.CLN.DIR, "/p194/p194_2010_F30c.tif")), TGO), TGO),
#       filename=paste0(FCC.CLN.DIR, "/TGO/1_clean/TGO_2010_F30c.tif"), overwrite=TRUE)
# 
# merge(mask(crop(brick(paste0(FCC.CLN.DIR, "/p193/p193_2013_F30c.tif")), TGO), TGO),
#       mask(crop(brick(paste0(FCC.CLN.DIR, "/p192/p192_2013_F30c.tif")), TGO), TGO),
#       mask(crop(brick(paste0(FCC.CLN.DIR, "/p194/p194_2012_F30c.tif")), TGO), TGO),
#       filename=paste0(FCC.CLN.DIR, "/TGO/1_clean/TGO_2013_F30c.tif"), overwrite=TRUE)


# Spatial cleaning of results (only from 2003 onwards) -----------------------------------

clean.spatial.forest(maps = stack(dir(paste0(FCC.CLN.DIR, "/FC", COV.FC, "/TGO"), pattern = "c\\.tif", full.names = TRUE)),
                     exclude = c(1), size = 6, connectedness = 8) # exclude the first layer (1987) for creating the forest / non-forest mask  

# # and put together in one stack
# tmp <- stack(dir(paste0(FCC.CLN.DIR, "/TGO/2_clean_05ha"), pattern=".*[[:digit:]]{4}.*", full.names=TRUE))
# writeRaster(tmp, paste0(FCC.CLN.DIR, "/TGO/2_clean_05ha/TGO_stack_F30c_05ha.tif"), overwrite=TRUE) 
# 
# # Loading PHCF-Filter -----------------------------------------------------
# # pixels with value 1 will only remain if they are in squares of 4 or adjacent to that
# # otherwise, if they are adjacent to a pixel with value '2' they will be converted to '2'
# # otherwise, they will be converted to zero
# 
# # compiling and loading C code that implements the PHCF filter
# system("R CMD SHLIB ./phcf_filter/phcf_filter.c")
# dyn.load("../phcf_filter/phcf_filter.so")
# 
# # defining function for facilitating the use of the C function
# phcf_pt <- function(defor) {
#   res <- .C("phcf_filter", nrow(defor), ncol(defor), as.integer(defor[]))
#   return(res[[3]])
# }
# 
# clean.spatial.fcc(maps = stack(raster(paste0(FCC.CLN.DIR, "/TGO/2_clean_05ha/TGO_2003_F30c_05ha.tif")), 
#                                raster(paste0(FCC.CLN.DIR, "/TGO/2_clean_05ha/TGO_2015_F30c_05ha.tif")), 
#                                raster(paste0(FCC.CLN.DIR, "/TGO/2_clean_05ha/TGO_2018_F30c_05ha.tif"))),
#                   size = 3, connectedness = 8)










