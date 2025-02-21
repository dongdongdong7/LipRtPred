# LipRtPred

This is an R package to predict lipid retention time. User can use this package to calculate three types of molecule descriptors and six types of molecule fingerprints.

## Load demo data

This demo data is subset of LipidMaps database.

```R
data("cmpDf_demo", package = "LipRtPred")
```

## Calculate molecule descriptors

### CDK

```R
descsDf <- GetCDK_MD(cmpDf = cmpDf_demo, thread = 1)
```

### RDKit

```R
descsDf <- GetRDKit_MD(cmpDf = cmpDf_demo, thread = 1)
```

### Mordred

```R
descsDf <- GetMordred_MD(cmpDf = cmpDf_demo, thread = 1)
```

## Calculate molecule fingerprints

### Morgan

```R
fpsDf <- GetMorgan_FP(cmpDf = cmpDf_demo, thread = 1)
```

### Feature Morgan

```R
fpsDf <- GetFMorgan_FP(cmpDf = cmpDf_demo, thread = 1)
```

### RDKit

```R
fpsDf <- GetRDKit_FP(cmpDf = cmpDf_demo, thread = 1)
```

### Atom Pairs

```R
fpsDf <- GetAp_FP(cmpDf = cmpDf_demo, thread = 1)
```

### Topological Torsion

```R
fpsDf <- GetTt_FP(cmpDf = cmpDf_demo, thread = 1)
```

### MACCS Keys

```R
fpsDf <- GetMaccs_FP(cmpDf = cmpDf_demo, thread = 1)
```

