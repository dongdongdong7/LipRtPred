# LipRtPred

This is an R package to predict lipid retention time. User can use this package to calculate three types of molecule descriptors and six types of molecule fingerprints.

## Load demo data

This demo data is subset of LipidMaps database.

```R
data("cmpDf_demo", package = "LipRtPred")
```

## Calculate molecule descriptors

```R
# CDK
descsDf <- GetCDK_MD(cmpDf = cmpDf_demo, thread = 1)
# RDKit
descsDf <- GetRDKit_MD(cmpDf = cmpDf_demo, thread = 1)
# Mordred
descsDf <- GetMordred_MD(cmpDf = cmpDf_demo, thread = 1)
```

## Calculate molecule fingerprints

```R
# Morgan
fpsDf <- GetMorgan_FP(cmpDf = cmpDf_demo, thread = 1)
# Feature Morgan
fpsDf <- GetFMorgan_FP(cmpDf = cmpDf_demo, thread = 1)
# RDKit
fpsDf <- GetRDKit_FP(cmpDf = cmpDf_demo, thread = 1)
# Atom Pairs
fpsDf <- GetAp_FP(cmpDf = cmpDf_demo, thread = 1)
# Topological Torsion
fpsDf <- GetTt_FP(cmpDf = cmpDf_demo, thread = 1)
# MACCS Keys
fpsDf <- GetMaccs_FP(cmpDf = cmpDf_demo, thread = 1)
```
