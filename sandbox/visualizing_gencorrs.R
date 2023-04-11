#saveRDS(genCorr_CRF_tbl, "/home/nur479/scratch3/PXS_pipeline/genCorr_CRF_tbl.rds")
genCorr_CRF_tbl <- readRDS("/home/nur479/scratch3/PXS_pipeline/genCorr_CRF_tbl.rds")

ggplot(genCorr_CRF_tbl, aes(x=CRF_fieldname, y=expo_fieldname, fill=gencorr)) +
  geom_tile() +
  geom_text(aes(label = round(gencorr,4))) +
  scale_fill_gradient2(low="red",mid="white", high="green", midpoint = 0) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=30)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=50)) +
  xlab("Clinical Risk Factor (CRF)") +
  ylab("Exposure") +
  labs(title = "Genetic Correlations between exposures and CRFs")
