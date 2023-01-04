

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(sf, tidyverse, rgeos, gtools, cclust, fs, extrafont, Hmisc, showtext, dismo, multcomp, glue, exactextractr, geodata, ggspatial, randomForest)

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
names(bioc) <- glue('bioc_{1:19}')

# To calculate the zonal --------------------------------------------------
bioc_znal <- exact_extract(bioc, bsns, 'mean') %>% as_tibble()
colnames(bioc_znal) <- gsub('mean.', '', colnames(bioc_znal))

# Join the two tables -----------------------------------------------------
tble <- bsns %>% st_drop_geometry %>% dplyr::select(NOMAH, NOMZH, NOMSZH) %>% as_tibble() %>% cbind(., bioc_znal) %>% as_tibble()

# Cluster using Random Forest ---------------------------------------------
mtrx <- tble[,4:ncol(tble)]
nfrs <- 50
ntrs <- 500
rfds <- RFdist(mtrx, mtry1 = 9, ntrs, nfrs, addcl1 = T, addcl2 = F, imp = T, oob.prox1 = T)
lbrf <- pamNew(rfds$cl1, 5)

bsns <- mutate(bsns, cluster = lbrf)
plot(dplyr::select(bsns, cluster))

tble <- mutate(tble, cluster = lbrf)
tble <- mutate(tble, cluster = factor(cluster, levels = as.character(1:5)))

# To make a boxplot -------------------------------------------------------
mtrx <- dplyr::select(tble, cluster, bioc_1:bioc_19)
mtrx <- gather(mtrx, variable, value, -cluster)
mtrx <- filter(mtrx, variable %in% c('bioc_1', 'bioc_10', 'bioc_11', 'bioc_12', 'bioc_16', 'bioc_17', 'bioc_18', 'bioc_19'))

lbls <- tibble(variable = c('bioc_1', 'bioc_10',  'bioc_11', 'bioc_12', 'bioc_16', 'bioc_17', 'bioc_18', 'bioc_19'), 
               nombre = c('Temp prom.', 'Temp. mes más cálido', 'Temp. mes más frío', 'Prec acum', 'Prec trim más húmedo', 'Prec trim más seco', 'Prec trim más cálido', 'Prec trim más frío'))

mtrx <- inner_join(mtrx, lbls, by = 'variable')
mtrx <- mutate(mtrx, nombre = factor(nombre, levels = c('Temp prom.', 'Temp. mes más cálido', 'Temp. mes más frío', 'Prec acum', 'Prec trim más húmedo', 'Prec trim más seco', 'Prec trim más cálido', 'Prec trim más frío')))


windowsFonts(Abadi = windowsFont('Abadi'))
font_add_google("Montserrat", "Montserrat")

gbox <- ggplot(data = mtrx, aes(x = cluster, y = value)) + 
  geom_boxplot() + 
  facet_wrap(~nombre, scales = 'free') + 
  labs(x = 'Clúster', y = 'Valor variable (C) - (mm)') +
  theme_minimal() + 
  theme(text = element_text(family = 'Montserrat'),
        strip.text.x = element_text(face = 'bold'), 
        axis.text.y = element_text(angle = 90, hjust = 0.5))
gbox

ggsave(plot = gbox, filename = 'png/graphs/boxplot.png', units = 'in', width = 9, height = 7, dpi = 300)


