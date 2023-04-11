library(tidyverse)
library(data.table)

loc_pheno3 <- "/n/groups/patel/uk_biobank/project_22881_671028/ukb671028.csv"
cols_keep <- c("eid","26285-0.0")
pheno3 <- as_tibble(fread(loc_pheno3, select=cols_keep)) %>%
  rename(PRS_T2D = `26285-0.0`)

dir_scratch <- "~/scratch3/PXS_pipeline/"
loc_pheno <- paste0(dir_scratch, "pheno_EC.txt")
pheno <- as_tibble(fread(loc_pheno)) %>%
  left_join(pheno3, by=c("IID"="eid"))


pheno_T2D <- pheno %>% select("IID",ends_with("T2D")) %>% drop_na()

ggplot(pheno_T2D, aes(x=PRS_T2D,y=PXS_T2D)) +
  geom_point(alpha=0.01) +
  geom_smooth(method="lm")


ggplot(pheno_T2D, aes(x=PRS_T2D,y=PXS_T2D)) +
  geom_hex()

cor.test(pheno_T2D$PXS_T2D, pheno_T2D$PRS_T2D)


####
library(sf)
library(spdep)

loc_pheno_full2 <- "/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"
## loading shape file
#loc_shp_MSOA <- "~/scratch3/key_data/infuse_MSOA_2011/infuse_msoa_lyr_2011_clipped.shp"
loc_shp_LA <- "~/scratch3/key_data/infuse_LA_2011/infuse_dist_lyr_2011_clipped.shp"
s <- sf::st_read(loc_shp_LA)

# extracts home and birthplace coordinates
cols_to_keep2 <- c("eid",paste0(c(129,130,22702,22704),"-0.0"))
pheno2 <- as_tibble(fread(loc_pheno_full2, select = cols_to_keep2)) %>%
  select("IID"="eid","birth_x"="130-0.0","birth_y"="129-0.0","home_x"="22702-0.0","home_y"="22704-0.0") %>%
  drop_na() %>%
  filter(birth_x > 0, birth_y > 0, home_x > 0, home_y > 0)

# joins phenotypes/exposures file with coordinates
pheno <- pheno %>% inner_join(pheno2, by="IID")

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
pheno <- pheno %>% filter(!(region_home_index %in% NI_indices),
                          !(region_birth_index %in% NI_indices))

# gets mean PRS and mean PXS per region
cols <- c("PXS_T2D","PRS_T2D")
for (i in s$region_index) {
  colmeans_home <- colMeans(pheno[pheno$region_home_index == i,cols], na.rm=TRUE)
  colmeans_birth <- colMeans(pheno[pheno$region_birth_index == i,cols], na.rm=TRUE)
  
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

#nb <- poly2nb(s, queen=TRUE)
#saveRDS(nb, paste0(dir_scratch,"UK_LA_neighbors.rds"))
nb <- readRDS(file=paste0(dir_scratch,"UK_LA_neighbors.rds"))
spdep::set.ZeroPolicyOption(TRUE)
lw <- nb2listw(nb, style="B", zero.policy=TRUE)

home_Is <- c()
home_ps <- c()
birth_Is <- c()
birth_ps <- c()

for (col in cols) {
  moran_home <- moran.test(as.data.frame(s_home)[,col], lw, na.action=na.exclude)
  home_I <- moran_home$estimate[1] %>% unname()
  home_p <- moran_home$p.value
  home_Is <- c(home_Is, home_I)
  home_ps <- c(home_ps, home_p)
  
  moran_birth <- moran.test(as.data.frame(s_birth)[,col], lw, na.action=na.exclude)
  birth_I <- moran_birth$estimate[1] %>% unname()
  birth_p <- moran_birth$p.value
  birth_Is <- c(birth_Is, birth_I)
  birth_ps <- c(birth_ps, birth_p)
  
  print(col)
}

morans_table <- tibble(term = cols) %>%
  mutate(
    home_I = home_Is, home_p = home_ps,
    birth_I = birth_Is, birth_p = birth_ps)

# function for plotting a trait
plot_geospatial <- function(shapefile, trait_col, trait_name, I, p, coords) {
  title <- paste0("Geospatial clustering of '",trait_name,"'\nbased on ", coords, " coordinates")
  subtitle <- paste0("Moran's I = ", round(I,3), " :: adjusted p-value = ", formatC(p, format="E", digits=2))
  
  
  ggplot(shapefile) +
    geom_sf(aes(fill=!!sym(trait_col)), size=0.2) +
    scale_fill_gradient2(low="green",mid="white", high="red",
                         midpoint = median(shapefile[[trait_col]], na.rm=TRUE)) +
    labs(title = title, subtitle = subtitle, fill = "Trait Value")
}

plot_geospatial(s_home[-NI_indices], "PRS_T2D", "PRS for T2D", morans_table$home_I[2], morans_table$home_p[2], "home")
plot_geospatial(s_birth[-NI_indices], "PRS_T2D", "PRS for T2D", morans_table$birth_I[2], morans_table$birth_p[2], "birth")


cor1 <- cor.test(s_home$PRS_T2D, s_home$PXS_T2D)
temp <- as_tibble(s_home) %>% select(n_home_pop,PRS_T2D, PXS_T2D) %>% drop_na()
cor2 <- cov.wt(temp[,-1], wt = temp$n_home_pop, cor=TRUE)
ggplot(s_home, aes(x=PRS_T2D, y=PXS_T2D)) +
  geom_point(aes(size = n_home_pop), alpha = 0.25) +
  #geom_smooth(method="lm") +
  labs(title = "Correlation between mean PXS and PRS for T2D in each HOME region",
       subtitle = paste0("r = ", round(cor1$estimate,3),
                         " :: p = ", formatC(cor1$p.value, 3, format="E"),
                         " :: population-weighted r = ", round(cor2$cor[1,2], 3)))

cor1 <- cor.test(s_birth$PRS_T2D, s_birth$PXS_T2D)
temp <- as_tibble(s_birth) %>% select(n_birth_pop,PRS_T2D, PXS_T2D) %>% drop_na()
cor2 <- cov.wt(temp[,-1], wt = temp$n_birth_pop, cor=TRUE)
ggplot(s_birth, aes(x=PRS_T2D, y=PXS_T2D)) +
  geom_point(aes(size = n_home_pop), alpha = 0.25) +
  #geom_smooth(method="lm") +
  labs(title = "Correlation between mean PXS and PRS for T2D in each BIRTH region",
       subtitle = paste0("r = ", round(cor1$estimate,3),
                         " :: p = ", formatC(cor1$p.value, 3, format="E"),
                         " :: population-weighted r = ", round(cor2$cor[1,2], 3)))

####

pheno <- pheno %>%
  left_join(s_home %>% select(region_home_index = region_index,
                              home_PXS_T2D = PXS_T2D,
                              home_PRS_T2D = PRS_T2D),
            by="region_home_index") %>%
  left_join(s_birth %>% select(region_birth_index = region_index,
                              birth_PXS_T2D = PXS_T2D,
                              birth_PRS_T2D = PRS_T2D),
            by="region_birth_index")

pheno <- pheno %>%
  mutate(
    resid_home_PXS_T2D = PXS_T2D - home_PXS_T2D,
    resid_home_PRS_T2D = PRS_T2D - home_PRS_T2D,
    resid_birth_PXS_T2D = PXS_T2D - birth_PXS_T2D,
    resid_birth_PRS_T2D = PRS_T2D - birth_PRS_T2D
  )
##
cor1 <- cor.test(pheno$PRS_T2D, pheno$PXS_T2D)
ggplot(pheno, aes(x=PRS_T2D,y=PXS_T2D)) +
  geom_point(alpha=0.01) +
  geom_smooth(method="lm") +
  labs(title = "Correlation between PXS and PRS for T2D",
       subtitle = paste0("r = ", round(cor1$estimate,4))) +
  xlab("PRS") +
  ylab("PXS")

##
cor1 <- cor.test(pheno$resid_home_PRS_T2D, pheno$resid_home_PXS_T2D)
ggplot(pheno, aes(x=resid_home_PRS_T2D,y=resid_home_PXS_T2D)) +
  geom_point(alpha=0.01) +
  geom_smooth(method="lm") +
  labs(title = "Correlation between PXS and PRS for T2D, adjsuted for home region",
       subtitle = paste0("r = ", round(cor1$estimate,4))) +
  xlab("home-adjusted PRS") +
  ylab("home-adjusted PXS")

##
cor1 <- cor.test(pheno$resid_birth_PRS_T2D, pheno$resid_birth_PXS_T2D)
ggplot(pheno, aes(x=resid_birth_PRS_T2D,y=resid_birth_PXS_T2D)) +
  geom_point(alpha=0.01) +
  geom_smooth(method="lm") +
  labs(title = "Correlation between PXS and PRS for T2D, adjsuted for birth region",
       subtitle = paste0("r = ", round(cor1$estimate,4))) +
  xlab("birth-adjusted PRS") +
  ylab("birth-adjusted PXS")
