# working off of this tutorial:
# https://mgimond.github.io/simple_moransI_example/

library(sf)
#library(sp)
library(spdep)
#library(tmap)
library(tidyverse)
library(data.table)
#library(broom)


## loading shape file
loc_shp_MSOA <- "~/scratch3/key_data/infuse_MSOA_2011/infuse_msoa_lyr_2011_clipped.shp"
loc_shp_LA <- "~/scratch3/key_data/infuse_LA_2011/infuse_dist_lyr_2011_clipped.shp"
s <- sf::st_read(loc_shp_LA)

## loading coord data
dir_script <- "~/jobs/PXS_pipeline/code/"
dir_scratch <- "~/scratch3/PXS_pipeline/"
#loc_pheno_full1 <- "/n/groups/patel/uk_biobank/main_data_34521/ukb34521.tab"
loc_pheno_full2 <- "/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"

cols_to_keep2 <- c("eid",paste0(c(129,130,22702,22704),"-0.0"))
#pheno2 <- as_tibble(fread(loc_pheno_full2, select = cols_to_keep2))
#pheno2 <- pheno2 %>% select("IID"="eid","birth_x"="130-0.0","birth_y"="129-0.0","home_x"="22702-0.0","home_y"="22704-0.0")
#write.table(pheno2, "~/scratch3/12-21_sf/UKB_coords.txt",sep="\t",quote=FALSE,row.names=FALSE)
pheno2 <- as_tibble(fread("~/scratch3/12-21_sf/UKB_coords.txt"))

subset <- pheno2 %>% drop_na() %>% filter(row_number() %in% sample(1:nrow(pheno2),5000))
#st_coords <- st_transform(st_sfc(st_point(coords),crs = 4326), 2163)

## plotting
ggplot(s) +
  #geom_sf() +
  geom_sf(data = s[380,], fill="blue") +
  geom_point(data=subset[1,],aes(x=home_x,home_y),color="red",alpha=1)

##
t1 <- Sys.time()
.GlobalEnv$i <- 1
subset$region_index <- apply(subset[1:40,], 1,function(row) {
  print(.GlobalEnv$i)
  .GlobalEnv$i <- .GlobalEnv$i + 1
  point <- st_point(row[c("home_x","home_y")] %>% unlist())
  #index <- (1:nrow(s))[st_intersects(point,s, sparse=FALSE)]
  index <- (1:nrow(s))[st_within(point,s, sparse=FALSE)]
  #index <- (1:nrow(s))[st_contains(s,point, sparse=FALSE)]
  return(index)
})
Sys.time() - t1
# st_intersects LA,20: 42.91s
# st_within LA,20: 15.13s
# st_within LA,40: 29.74s

point <- st_point(subset[1,c("home_x","home_y")] %>% unlist())
out1 <- st_intersects(s,point, sparse=FALSE)

#

simple <- as.data.frame(subset[1:20,] %>% select(x=home_x,y=home_y))
points <- st_as_sf(simple, coords=c("x","y"), crs = st_crs(s))
grid <- st_make_grid(points)
points_with_regions <- st_join(points, grid)

###
points <- st_as_sf(subset, coords=c("home_x","home_y"), crs = st_crs(s))

t1 <- Sys.time()
out1 <- st_contains(s,points, sparse=FALSE)
#out1 <- st_within(points, s, sparse=FALSE)
Sys.time() - t1

#str(out1)
for (i in 1:nrow(subset)) {
  n_regions <- nrow(s)
  region_index <- (1:n_regions)[out1[,i]]
  if (length(region_index)==0) {
    #region_index <- as.numeric(NA)
  }
  subset$region_index[i] <- region_index
}
# water dwellers
ggplot(s) +
  geom_sf() +
  geom_sf(data=s[out2,], fill="blue") +
  xlim(500000,600000) +
  ylim(173000,181000) +
  #geom_sf(data = s[380,], fill="blue") +
  geom_point(data=subset %>% filter(is.na(region_index)),aes(x=home_x,home_y),color="red",alpha=1)

water_dwellers <- st_as_sf(subset[is.na(subset$region_index),],
                           coords=c("home_x","home_y"), crs = st_crs(s))
out2 <- st_nearest_feature(water_dwellers, s, sparse=FALSE)
subset[is.na(subset$region_index),"region_index"] <- out2

#####
t1 <- Sys.time()
nb <- poly2nb(s, queen=TRUE)
Sys.time() - t1
# 10: 1.27s
# 40: 5.31s
# 100: 18.34s
# 200: 51.22s
# 404: 975.48s