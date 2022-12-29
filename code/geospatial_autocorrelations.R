library(tidyverse)
library(data.table)
library(sf)
library(spdep)
library(ggrepel)

dir_script <- "~/jobs/PXS_pipeline/code/"
dir_scratch <- "~/scratch3/PXS_pipeline/"
loc_pheno_EC <- paste0(dir_scratch,"pheno_EC.txt")
loc_fields <- paste0(dir_scratch,"fields_tbl.txt")
# phenotype file containing home and birthplace coordinates
loc_pheno_full2 <- "/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"
## loading shape file
#loc_shp_MSOA <- "~/scratch3/key_data/infuse_MSOA_2011/infuse_msoa_lyr_2011_clipped.shp"
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
pheno <- pheno %>% filter(!(region_home_index %in% NI_indices),
                          !(region_birth_index %in% NI_indices))

# gets mean exposure value for each region
fields_tbl <- as_tibble(fread(loc_fields))
expos <- c(fields_tbl$term[fields_tbl$use_type=="exposure"],"PXS_T2D")

for (i in s$region_index) {
  colmeans_home <- colMeans(pheno[pheno$region_home_index == i,expos], na.rm=TRUE)
  colmeans_birth <- colMeans(pheno[pheno$region_birth_index == i,expos], na.rm=TRUE)
  colmeans_stayers <- colMeans(pheno[(pheno$region_birth_index == i) & (pheno$region_home_index == i),expos], na.rm=TRUE)
  
  if (i == s$region_index[1]) {
    region_home_means <- data.frame(t(colmeans_home))
    region_birth_means <- data.frame(t(colmeans_birth))
    region_stayers_means <- data.frame(t(colmeans_stayers))
  } else {
    region_home_means <- rbind(region_home_means,t(colmeans_home))
    region_birth_means <- rbind(region_birth_means,t(colmeans_birth))
    region_stayers_means <- rbind(region_stayers_means,t(colmeans_stayers))
  }
  if (i %% 50 == 0) {print(i)}
}
s_home <- cbind(s,region_home_means)
s_birth <- cbind(s,region_birth_means)
s_stayers <- cbind(s,region_stayers_means)

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
morans_table <- tibble(term = expos) %>%
  left_join(fields_tbl %>% select(term, fieldname, Meaning), by="term") %>%
  mutate(
    home_I = home_Is, home_p = p.adjust(home_ps,method="fdr"),
    birth_I = birth_Is, birth_p = p.adjust(birth_ps,method="fdr"))
morans_table[morans_table$term=="PXS_T2D","fieldname"] <- "PXS for Type II Diabetes"
# function for plotting a trait
plot_geospatial <- function(shapefile, trait_col, trait_name, I, p, coords) {
  title <- paste0("Geospatial clustering of '",trait_name,"'\nbased on ", coords, " coordinates")
  subtitle <- paste0("Moran's I = ", round(I,3), " :: adjusted p-value = ", formatC(p, format="E", digits=2))
  
  
  ggplot(shapefile) +
    geom_sf(aes(fill=!!sym(trait_col)), size=0.2) +
    scale_fill_gradient2(low="green",mid="white", high="red",
                         midpoint = median(shapefile[[trait_col]])) +
    labs(title = title, subtitle = subtitle, fill = "Trait Value")
}
dir.create(paste0(dir_scratch,"figures/birth_I"), recursive=TRUE)
dir.create(paste0(dir_scratch,"figures/home_I"), recursive=TRUE)
# loops through each trait and birth/home and saves geospatial plot of exposure values
for (coords in c("birth", "home")) {
  for (i in 1:nrow(morans_table)) {
    if (coords=="birth") {
      shapefile <- s_birth[-NI_indices,]
      I <- morans_table$birth_I[i]
      p <- morans_table$birth_p[i]
      dir_out <- paste0(dir_scratch,"figures/birth_I/")
    } else if (coords=="home") {
      shapefile <- s_home[-NI_indices,]
      I <- morans_table$home_I[i]
      p <- morans_table$home_p[i]
      dir_out <- paste0(dir_scratch,"figures/home_I/")
    }
    
    trait_col <- morans_table$term[i]
    
    if (is.na(morans_table$Meaning[i])) {
      trait_name <- morans_table$fieldname[i]
    } else {
      trait_name <- paste0(morans_table$fieldname[i],": ", morans_table$Meaning[i])
    }
    
    gg <- plot_geospatial(shapefile = shapefile,
                    trait_col = trait_col, trait_name = trait_name,
                    I = I, p = p, coords = coords)
    
    loc_out <- paste0(dir_scratch,"figures/",coords,"_I/",trait_col,"_",coords,"_geospatial.png")
    ggsave(loc_out, plot=gg, width = 2000, height = 3000, units = "px")
    print(paste("Saved plot for", coords,"Moran's I for", trait_col))
  }
}


# plots home_I vs birth_I
cor1 <- cor(morans_table$birth_I, morans_table$home_I)
ggplot(morans_table, aes(x=birth_I,y=home_I)) +
  geom_abline(slope=1) +
  geom_smooth(method="lm") +
  geom_point(aes(color = home_I>birth_I)) +
  geom_label_repel(data=morans_table %>% filter(home_I>birth_I), aes(label=fieldname)) +
  xlim(-0.05,0.7) + ylim(-0.05,0.7) +
  xlab("Moran's I using birthplace coordinates") +
  ylab("Moran's I using current home coordinates") +
  labs(title = "Comparison of geospatial clustering of exposures based on birthplace vs current home",
       subtitle = paste0("r = ", round(cor1,3))) +
  theme(legend.position = "none")
loc_out <- paste0(dir_scratch,"figures/home_vs_birth_I.png")
ggsave(loc_out, width = 4000, height = 3000, units = "px")


# Working with movers and stayers
movers <- pheno %>% filter(region_birth_index != region_home_index)
stayers <- pheno %>% filter(region_birth_index == region_home_index)

dir.create(paste0(dir_scratch,"figures/migration_shifts"), recursive=TRUE)
for (i in 1:nrow(morans_table)) {
  trait_col <- morans_table$term[i]
  
  if (is.na(morans_table$Meaning[i])) {trait_name <- morans_table$fieldname[i]
  } else {trait_name <- paste0(morans_table$fieldname[i],": ", morans_table$Meaning[i])}
  
  s_trait <- tibble(region_index = 1:nrow(s), trait = s_stayers[[trait_col]])
  
  movers_trait <- pheno %>% filter(region_birth_index != region_home_index) %>%
    select(FID, IID, region_home_index, region_birth_index) %>%
    left_join(s_trait %>% rename(trait_home = trait), by=c("region_home_index"="region_index")) %>%
    left_join(s_trait %>% rename(trait_birth = trait), by=c("region_birth_index"="region_index")) %>%
    mutate(trait_delta = trait_home - trait_birth)
  
  t1 <- t.test(movers_trait$trait_delta)
  movers_trait2 <- movers_trait %>% select(-trait_delta) %>% pivot_longer(
    cols = starts_with("trait_"),
    names_to = "coords",
    names_prefix = "trait_",
    values_to = "mean_value"
  )
  # plots birth and home distributions
  ggplot(movers_trait2) +
    geom_histogram(aes(x=mean_value, fill = coords), alpha = 0.5) +
    geom_vline(xintercept = mean(movers_trait$trait_birth, na.rm=TRUE), color="#F8766D") +
    geom_vline(xintercept = mean(movers_trait$trait_home, na.rm=TRUE), color="#00BFC4") +
    xlab("Mean trait value in individuals' region (may be inverse rank normalized)") +
    ylab("Number of indviduals (movers)") +
    labs(title = paste0("Distribution of mean '",trait_name,"' value at individuals' birthplace/current home region among movers"),
         subtitle = paste0("Mean difference = ", formatC(t1$estimate,digits=3, format = "fg"),
                           " :: p-value (t-test) = ", formatC(t1$p.value,digits=2, format = "E"),
                           " :: Mean region value determined by current inhabitants who were also born in same region (stayers)"),
         fill = "Region")
  
  loc_out <- paste0(dir_scratch,"figures/migration_shifts/",trait_col,"_migration_shift.png")
  ggsave(loc_out, width = 3000, height = 2000, units = "px")
  
  # plots delta
  ggplot(movers_trait) +
    geom_histogram(aes(x=trait_delta), alpha = 0.5, fill="#00BA38") +
    geom_vline(xintercept = 0, color="black") +
    geom_vline(xintercept = mean(movers_trait$trait_delta, na.rm=TRUE), color="darkgreen") +
    xlab("Mean region shift (home - birth) (may be inverse rank normalized)") +
    ylab("Number of indviduals (movers)") +
    labs(title = paste0("Distribution of birth-->home mean region shift for '",trait_name,"' among movers"),
         subtitle = paste0("Mean difference = ", formatC(t1$estimate,digits=3, format = "fg"),
                           " :: p-value (t-test) = ", formatC(t1$p.value,digits=2, format = "E"),
                           " :: Mean region value determined by current inhabitants who were also born in same region (stayers)"))
  
  loc_out <- paste0(dir_scratch,"figures/migration_shifts/",trait_col,"_migration_shift2.png")
  ggsave(loc_out, width = 3000, height = 2000, units = "px")
  
  print(trait_col)
}
