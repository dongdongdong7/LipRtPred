library(magrittr)
set.seed(1)
lipidmapsCmpTb <- pubcmpR::load_lipidmapsCmpTb()
idx <- sort(sample(1:nrow(lipidmapsCmpTb), size = 100, replace = FALSE))
cmpDf_demo <- lipidmapsCmpTb[idx, ]
cmpDf_demo <- cmpDf_demo %>%
  dplyr::select(LM_ID, SMILES)
colnames(cmpDf_demo) <- c("id", "smiles")
#usethis::use_data(cmpDf_demo)
