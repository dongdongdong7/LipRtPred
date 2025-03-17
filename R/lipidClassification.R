# Classify lipid
# Barry Song
# 250314

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# .lipidClassification(smi = "C1C2[C@@]3([H])CC[C@]4(C)C(=O)CC[C@@]4([H])[C@]3([H])CCC=2C=C(O)C=1",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .lipidClassification(smi = "[C@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)(COP(CC[N+](C)(C)C)(=O)[O-])COC(CCCCCCCCCCCCCCC)=O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .lipidClassification(smi = "C(=O)(COP(=O)(O)O)COCCCCCCCCCCCCCCCC",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .lipidClassification(smi = "O(C(CCCCCCCCCCCCCCCCCC)CO[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1)C(=O)CCCCCCCCCCCCCCC",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.lipidClassification <- function(smi, scriptPath){

  # Whether Sterol
  steroidSkeleton_position <- .searchSteroidSkeleton(smi = smi, scriptPath = scriptPath)
  secosteroidSkeleton_position <- .searchSecosteroidSkeleton(smi = smi, scriptPath = scriptPath)

  if(length(steroidSkeleton_position) != 0 | length(secosteroidSkeleton_position) != 0) return("Sterol")

  # Whether Sphingolipid
  sphingosie_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
  spisulosine_position <- .searchSpisulosine(smi = smi, scriptPath = scriptPath)

  if(length(sphingosie_position) != 0 | length(spisulosine_position) != 0) return("Sphingolipid")

  # Whether Glycerolipid or Glycerophospholipid
  glycerolEster_position <- .searchGlycerolEster(smi = smi, scriptPath = scriptPath)
  glycerolEther_position <- .searchGlycerolEther(smi = smi, scriptPath = scriptPath)
  dihydroxyacetoneEster_position <- .searchDihydroxyacetoneEster(smi = smi, scriptPath = scriptPath)
  dihydroxyacetoneEther_position <- .searchDihydroxyacetoneEther(smi = smi, scriptPath = scriptPath)
  glycerolPhosphate_position <- .searchGlycerolPhosphate(smi = smi, scriptPath = scriptPath)

  if(length(glycerolPhosphate_position) != 0) return("Glycerophospholipid")
  if(length(glycerolEster_position) != 0 | length(glycerolEther_position) != 0 | length(dihydroxyacetoneEster_position) != 0 | length(dihydroxyacetoneEther_position) != 0) return("Glycerolipid")

  # Whether Fatty Acyls
  acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  amide_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
  thioester_position <- .searchThioester(smi = smi, scriptPath = scriptPath)

  if(length(acyloxy_position) != 0 | length(amide_position) != 0 | length(thioester_position) != 0) return("Fatty Acyls")

  fattyAlcohol_position <- .searchFattyAlcohols(smi = smi, scriptPath = scriptPath)
  if(length(fattyAlcohol_position) != 0) return("Fatty Alcohol")

  fattyAldehyde_position <- .searchFattyAldehydes(smi = smi, scriptPath = scriptPath)
  if(length(fattyAldehyde_position) != 0) return("Fatty Aldehyde")

  fattyNitrile <- .searchFattyNitriles(smi = smi, scriptPath = scriptPath)
  if(length(fattyNitrile) != 0) return("Fatty nitrile")

  fattyEther_position <- .searchFattyEthers(smi = smi, scriptPath = scriptPath)
  if(length(fattyEther_position) != 0) return("Fatty Ether")

  hydrocarbons_position <- .searchHydrocarbons(smi = smi, scriptPath = scriptPath)
  if(length(hydrocarbons_position) != 0) return("Hydrocarbons")

  return("No Lipid")
}
