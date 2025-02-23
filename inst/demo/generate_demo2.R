# This is the script to generate cmpDf_demo2 object.
# The cmpDf_demo2.xlsx is subset of CQS lipid dataset.
# Please see G:\spd\Projects\2025\Lipid RT Prediction\Progress\Experiment3 LipRtPred predicts lipid retention time\CQS241215_7600_EC C18_Biomatrix sample\raw
cmpDf_demo2 <- openxlsx::read.xlsx(system.file("demo", "cmpDf_demo2.xlsx", package = "LipRtPred"), sheet = 1)
cmpDf_demo2 <- dplyr::as_tibble(cmpDf_demo2)
if(length(which(is.na(cmpDf_demo2$smiles))) == 0){
  # usethis::use_data(cmpDf_demo2)
}
