library(tidyverse)
library(data.table)
library(sf)
library(spdep)

dir_script <- "~/jobs/PXS_pipeline/code/"
dir_scratch <- "~/scratch3/PXS_pipeline/"
loc_pheno_EC <- paste0(dir_scratch,"pheno_EC.txt")
loc_fields <- paste0(dir_scratch,"fields_tbl.txt")
# phenotype file containing home and birthplace coordinates
loc_pheno_full2 <- "/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"
## loading shape file
loc_shp_MSOA <- "~/scratch3/key_data/infuse_MSOA_2011/infuse_msoa_lyr_2011_clipped.shp"
loc_shp_LA <- "~/scratch3/key_data/infuse_LA_2011/infuse_dist_lyr_2011_clipped.shp"
s <- sf::st_read(loc_shp_LA)

# extracts home and birthplace coordinates
cols_to_keep2 <- c("eid",paste0(c(129,130,22702,22704),"-0.0"))
pheno2 <- as_tibble(fread(loc_pheno_full2, select = cols_to_keep2)) %>%
  select("IID"="eid","birth_x"="130-0.0","birth_y"="129-0.0","home_x"="22702-0.0","home_y"="22704-0.0") %>%
  drop_na()

# joins phenotypes/exposures file with coordinates
pheno <- as_tibble(fread(loc_pheno_EC)) %>%
  inner_join(pheno2, by="IID")
# 379,785 individuals remain

# places each individual in their region
temp <- matrix(rep(1:nrow(s),nrow(pheno)), nrow=nrow(s), byrow=FALSE)

points_home <- st_as_sf(pheno, coords=c("home_x","home_y"), crs = st_crs(s))
out_home <- st_contains(s,points_home, sparse=FALSE)
region_home_indices <- colSums(temp * out_home)
home_water_dwellers <- st_as_sf(pheno[region_home_indices==0,], coords=c("home_x","home_y"), crs = st_crs(s))
region_home_indices[region_home_indices==0] <- st_nearest_feature(home_water_dwellers, s, sparse=FALSE)
pheno$region_home_index <- region_home_indices

points_birth <- st_as_sf(pheno, coords=c("birth_x","birth_y"), crs = st_crs(s))
out_birth <- st_contains(s,points_birth, sparse=FALSE)
region_birth_indices <- colSums(temp * out_birth)
birth_water_dwellers <- st_as_sf(pheno[region_birth_indices==0,], coords=c("birth_x","birth_y"), crs = st_crs(s))
region_birth_indices[region_birth_indices==0] <- st_nearest_feature(birth_water_dwellers, s, sparse=FALSE)
pheno$region_birth_index <- region_birth_indices

# removes N.Ireland regions and indviduals
s$region_index <- 1:nrow(s)
s$n_home_pop <- rowSums(out_home)
s$n_birth_pop <- rowSums(out_birth)
NI_indices <- s[!(substring(s$geo_code,1,1) %in% c("E","W","S")),]$region_index
#s <- s[-NI_indices,]
pheno <- pheno %>% filter(!(region_home_index %in% NI_indices),
                          !(region_birth_index %in% NI_indices))

# gets mean exposure value for each region
fields_tbl <- as_tibble(fread(loc_fields))
expos <- fields_tbl$term[fields_tbl$use_type=="exposure"]
#expos <- fields_tbl$term

for (i in s$region_index) {
  colmeans_home <- colMeans(pheno[pheno$region_home_index == i,expos], na.rm=TRUE)
  colmeans_birth <- colMeans(pheno[pheno$region_birth_index == i,expos], na.rm=TRUE)
  
  if (i == s$region_index[1]) {
    region_home_means <- data.frame(t(colmeans_home))
    region_birth_means <- data.frame(t(colmeans_birth))
  } else {
    region_home_means <- rbind(region_home_means,t(colmeans_home))
    region_birth_means <- rbind(region_birth_means,t(colmeans_birth))
  }
  if (i %% 50 == 0) {print(i)}
}
s_home <- cbind(s,region_home_means)
s_birth <- cbind(s,region_birth_means)

# figures out neighboring regions (TAKES A WHILE, LIKE 16 minutes)
#nb <- poly2nb(s, queen=TRUE)
#saveRDS(nb, paste0(dir_scratch,"UK_LA_neighbors.rds"))
nb <- readRDS(file=paste0(dir_scratch,"UK_LA_neighbors.rds"))
#nb <- nb[-NI_indices]
# adds weights to neighbors, use style "B" like in Abdellaoui et al. 2019
spdep::set.ZeroPolicyOption(TRUE)
lw <- nb2listw(nb, style="B", zero.policy=TRUE)

# runs through each exposure and computes I
home_Is <- c()
home_ps <- c()
birth_Is <- c()
birth_ps <- c()
for (expo in expos) {
  moran_home <- moran.test(as.data.frame(s_home)[,expo], lw, na.action=na.exclude)
  home_I <- moran_home$estimate[1] %>% unname()
  home_p <- moran_home$p.value
  home_Is <- c(home_Is, home_I)
  home_ps <- c(home_ps, home_p)
  
  moran_birth <- moran.test(as.data.frame(s_birth)[,expo], lw, na.action=na.exclude)
  birth_I <- moran_birth$estimate[1] %>% unname()
  birth_p <- moran_birth$p.value
  birth_Is <- c(birth_Is, birth_I)
  birth_ps <- c(birth_ps, birth_p)
  
  print(expo)
}

fields_tbl$home_I[fields_tbl$use_type=="exposure"] <- home_Is
#fields_tbl$home_I <- home_Is
fields_tbl$home_p[fields_tbl$use_type=="exposure"] <- p.adjust(home_ps,method="fdr")
#fields_tbl$home_p <- home_ps

#fields_tbl$birth_I <- birth_Is
fields_tbl$birth_I[fields_tbl$use_type=="exposure"] <- birth_Is
#fields_tbl$birth_p <- birth_ps
fields_tbl$birth_p[fields_tbl$use_type=="exposure"] <- p.adjust(birth_ps,method="fdr")

fields_tbl$delta_I <- fields_tbl$home_I - fields_tbl$birth_I

fields_tbl %>% arrange(birth_p)
fields_tbl %>% arrange(-abs(delta_I))

library(ggrepel)
ggplot(fields_tbl, aes(x=birth_I,y=home_I)) +
  geom_abline(slope=1) +
  geom_smooth(method="lm") +
  geom_point() +
  geom_label_repel(data=fields_tbl %>% filter(home_I>birth_I), aes(label=fieldname))

# ggplot(s_birth[-NI_indices]) +
#   geom_sf(aes(fill=pc11), size=0.2) +
#   scale_fill_gradient(low="red",high="white") #scale_fill_brewer(palette = "YlOrRd")
ggplot(s) +
  geom_sf(aes(fill=f728), size=0.2) +
  labs(title="Spatial distribution of 'Number of vehicles in household'",
       substitle=paste0("Values are inverse-ranked transformed. I = ",fields_tbl$home_I[fields_tbl$term == "f728"]))


### f728 - number of vehicles in the household ###
s_f728 <- tibble(region_index = 1:nrow(s), f728 = s_home$f728)
#f728_home = s_home$f728,
#f728_birth = s_birth$f728)
movers <- pheno %>% filter(region_birth_index != region_home_index) %>%
  select(FID, IID, region_home_index, region_birth_index) %>%
  left_join(s_f728 %>% rename("home_f728"="f728"), by=c("region_home_index"="region_index")) %>%
  left_join(s_f728 %>% rename("birth_f728"="f728"), by=c("region_birth_index"="region_index")) %>%
  mutate(delta_f728 = home_f728 - birth_f728)

t1 <- t.test(movers$delta_f728)
ggplot(movers, aes(x=delta_f728)) +
  geom_histogram(binwidth = 0.05) +
  xlim(-1,1) +
  geom_vline(xintercept = 0) +
  #geom_rect(ymin=-Inf,ymax=Inf, xmin=t1$conf.int[1], xmax=t1$conf.int[2], fill="red", alpha=0.5) +
  geom_vline(xintercept = t1$estimate, color="red")


### f1160 - sleep duration ###
s_f1160 <- tibble(region_index = 1:nrow(s), f1160 = s_home$f1160)
#f1160_home = s_home$f1160,
#f1160_birth = s_birth$f1160)
movers <- pheno %>% filter(region_birth_index != region_home_index) %>%
  select(FID, IID, region_home_index, region_birth_index) %>%
  left_join(s_f1160 %>% rename("home_f1160"="f1160"), by=c("region_home_index"="region_index")) %>%
  left_join(s_f1160 %>% rename("birth_f1160"="f1160"), by=c("region_birth_index"="region_index")) %>%
  mutate(delta_f1160 = home_f1160 - birth_f1160)

t1 <- t.test(movers$delta_f1160)
ggplot(movers, aes(x=delta_f1160)) +
  geom_histogram(binwidth = 0.05) +
  xlim(-1,1) +
  geom_vline(xintercept = 0) +
  #geom_rect(ymin=-Inf,ymax=Inf, xmin=t1$conf.int[1], xmax=t1$conf.int[2], fill="red", alpha=0.5) +
  geom_vline(xintercept = t1$estimate, color="red")
