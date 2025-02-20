library(magrittr)
lipidmapsCmpTb <- pubcmpR::load_lipidmapsCmpTb()
set.seed(2)
idx <- sort(sample(1:nrow(lipidmapsCmpTb), size = 100, replace = FALSE))
cmpDf_demo <- lipidmapsCmpTb[idx, ]
cmpDf_demo <- cmpDf_demo %>%
  dplyr::select(LM_ID, SMILES)
colnames(cmpDf_demo) <- c("id", "smiles")
if(length(which(is.na(cmpDf_demo$smiles))) == 0){
  # usethis::use_data(cmpDf_demo)
}
