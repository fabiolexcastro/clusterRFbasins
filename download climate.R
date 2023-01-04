

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(sf, tidyverse, rgeos, gtools, cclust, fs, Hmisc, dismo, multcomp, glue, exactextractr, geodata, ggspatial, randomForest)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

source('functionsRF.R')

# Load data ---------------------------------------------------------------
bsns <- st_read('gpkg/basins_col.gpkg')
bsns <- st_transform(bsns, crs = st_crs(4326))

# Download ----------------------------------------------------------------
prec <- worldclim_country(country = 'COL', var = 'prec', path = 'tmpr')
tmax <- worldclim_country(country = 'COL', var = 'tmax', path = 'tmpr')
tavg <- worldclim_country(country = 'COL', var = 'tavg', path = 'tmpr')
tmin <- worldclim_country(country = 'COL', var = 'tmin', path = 'tmpr')

# Masking and cropping ----------------------------------------------------
prec <- terra::crop(prec, vect(bsns)) %>% terra::mask(., vect(bsns))
tmax <- terra::crop(tmax, vect(bsns)) %>% terra::mask(., vect(bsns))
tavg <- terra::crop(tavg, vect(bsns)) %>% terra::mask(., vect(bsns))
tmin <- terra::crop(tmin, vect(bsns)) %>% terra::mask(., vect(bsns))

names(prec) <- glue('prec_{1:12}')
names(tmax) <- glue('tmax_{1:12}')
names(tavg) <- glue('tavg_{1:12}')
names(tmin) <- glue('tmin_{1:12}')

# To calculate the bioclim ------------------------------------------------
bioc <- dismo::biovars(prec = stack(prec), tmin = stack(tmin), tmax = stack(tmax))
colnames(bioc) <- glue('bioc_{1:19}')

# To calculate the zonal --------------------------------------------------
bioc_znal <- exact_extract(bioc, bsns, 'mean') %>% as_tibble()
colnames(bioc_znal) <- gsub('mean.', '', colnames(bioc_znal))


# Cluster using Random Forest ---------------------------------------------

rf.clust <- function(occ, nforest, ntrees, nVars, nclasses){
  
  datRF_presences <- occ[,3:ncol(occ)]
  print(nrow(datRF))
  
  attach(datRF_presences)
  no.forests <- nforest
  no.trees <- ntrees
  distRF_presences <- RFdist(datRF_presences, mtry1 = nVars, no.trees, no.forests, addcl1 = T, addcl2 = F, imp = T, oob.prox1 = T)
  no.presencesclasses <- nclasses
  labelRF <- pamNew(distRF_presences$cl1, no.presencesclasses)
  print(table(labelRF))
  clusterdata <- hclust(as.dist(distRF_presences$cl1), method = 'ward.D2')
  
  return(list(labelRF, clusterdata))
  
}



