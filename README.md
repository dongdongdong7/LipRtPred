# LipRtPred

This is an R package to predict lipid retention time. User can use this package to calculate three types of molecule descriptors and six types of molecule fingerprints.

## Download

Please install miniconda before use this package.

```R
# install.packages("rcdk")
# install.packages("reticulate")
# options(timeout = 600)
devtools::install_github("dongdongdong7/LipRtPred")
```

## Load demo data

This demo data is subset of LipidMaps database.

```R
data("cmpDf_demo2", package = "LipRtPred")
```

## Calculate molecule descriptors

```R
# CDK
descsDf <- GetCDK_MD(cmpDf = cmpDf_demo2, thread = 3)
# RDKit
descsDf <- GetRDKit_MD(cmpDf = cmpDf_demo2, thread = 3)
# Mordred
descsDf <- GetMordred_MD(cmpDf = cmpDf_demo2, thread = 3)
```

## Calculate molecule fingerprints

```R
# Morgan
fpsDf <- GetMorgan_FP(cmpDf = cmpDf_demo2, thread = 1)
# Feature Morgan
fpsDf <- GetFMorgan_FP(cmpDf = cmpDf_demo2, thread = 1)
# RDKit
fpsDf <- GetRDKit_FP(cmpDf = cmpDf_demo2, thread = 1)
# Atom Pairs
fpsDf <- GetAp_FP(cmpDf = cmpDf_demo2, thread = 1)
# Topological Torsion
fpsDf <- GetTt_FP(cmpDf = cmpDf_demo2, thread = 1)
# MACCS Keys
fpsDf <- GetMaccs_FP(cmpDf = cmpDf_demo2, thread = 1)
# CQS
fpsDf <- GetCQS_FP(cmpDf = cmpDf_demp2, thread = 3)
```

## Build machine learning model

### Preprocess raw data

```R
inputDf <- filterColumns(inputDf = descsDf)
tmp <- split_data(inputDf = inputDf, percentage = 0.8)
trainingDf <- tmp$training
testingDf <- tmp$testing
```

### Training

```R
model_rf <- build_rf(trainingDf = trainingDf)
model_xgb <- build_xgb(trainingDf = trainingDf)
model_brnn <- build_brnn(trainingDf = trainingDf)
```

### Evaluate

```R
evaluate_model(testingDf = testingDf, model = model_rf)
evaluate_model(testingDf = testingDf, model = model_xgb)
evaluate_model(testingDf = testingDf, model = model_brnn)
```

