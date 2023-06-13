loc_40PCs <- "~/scratch3/key_data/UKB_40PCs_500k.txt"
PCs40 <- as_tibble(fread(loc_40PCs)) %>% select(-FID)
PXS_T2D <- as_tibble(fread("../scratch/T2D/PXS_T2D.txt"))

PCs40_filtered <- PCs40 %>% filter(IID %in% PXS_T2D$IID)

ggplot(PCs40, aes(x=pc1, y=pc2)) +
  geom_point(alpha=0.02)

ggplot(PCs40_filtered, aes(x=pc1, y=pc2)) +
  geom_point(data=PCs40 %>% filter(!IID %in% PXS_T2D$IID), alpha=0.02, color="black") +
  geom_point(data=PCs40_filtered, alpha=0.02, color="red")

###

sf <- as_tibble(fread("../scratch/T2D/LMM_PXS_T2D_bgen.txt"))
median(sf$CHISQ_BOLT_LMM_INF) / qchisq(0.5,1)
sf_mini <- as_tibble(fread("../scratch/T2D/LMM_PXS_T2D.txt"))
median(sf_mini$CHISQ_BOLT_LMM_INF) / qchisq(0.5,1)
