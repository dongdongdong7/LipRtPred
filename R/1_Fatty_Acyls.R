# Search substructure of Fatty Acyls (FA)
# Barry Song
# 250307

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py
# minium_C: the minimum C number of hydrocarbons

# 1. Acyloxy R1-C(=O)O-R2
# (1) Search acyloxy group
# .searchAcyloxy(smi = "O=C(C(C(C)C(O)C)N)O",
#                scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchAcyloxy <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CX3;$(C=O);$(C-O)](=O)[OX2]",
                       scriptPath = scriptPath)[[1]]
}

# (2) Search acyl chain of acyloxy group
# .searchAcyloxy_AcylChain(smi = "C(CC(C)CCCCCCCCCC(=O)O)(=O)O",
#                          scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchAcyloxy_AcylChain <- function(smi, scriptPath){
  acyloxy_position <- .searchAcyloxy(smi = smi,
                                     scriptPath = scriptPath)
  lapply(acyloxy_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# (3) Search alkoxy chain of acyloxy group
# .searchAcyloxy_AlkoxyChain(smi = "CCCCCCCCCCCCCCCC(OCCCCCCCCCCCCCCCC)=O",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchAcyloxy_AlkoxyChain <- function(smi, scriptPath){
  acyloxy_position <- .searchAcyloxy(smi = smi,
                                     scriptPath = scriptPath)
  if(length(acyloxy_position) != 0){
    posi <- lapply(acyloxy_position, function(x) {
      .TraverseMolecule(smi = smi,
                        start_atom_idx = x[3],
                        non_traversable_atom_idx = x[-3],
                        scriptPath = scriptPath)[-1]
    })
    posi <- posi[sapply(posi, function(x) {
      if(length(x) == 0) return(FALSE)
      else return(TRUE)
    })]
  }else posi <- list()
  return(posi)
}

# 2. Amide R1-C(=O)N-R2
# (1) Search amide group
# .searchAmide(smi = "C(=C/C=C/C=C/C(=O)NC1=C(O)CCC1=O)\\[C@@]1(C=C(C(=O)C[C@H]1O)NC(=O)/C=C/C(C)CC(C)CCCC)O",
#              scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchAmide <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CX3;$(C=O);$(C-N)](=O)N",
                       scriptPath = scriptPath)[[1]]
}

# (2) Search acyl chain of amide group
# .searchAmide_AcylChain(smi = "C(=C/C=C/C=C/C(=O)NC1=C(O)CCC1=O)\\[C@@]1(C=C(C(=O)C[C@H]1O)NC(=O)/C=C/C(C)CC(C)CCCC)O",
#                        scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchAmide_AcylChain <- function(smi, scriptPath){
  amide_position <- .searchAmide(smi = smi,
                                 scriptPath = scriptPath)

  lapply(amide_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# (3) Search alkylamino chain of amide group
# .searchAmide_AlkylaminoChain(smi = "C(=C/C=C/C=C/C(=O)NC1=C(O)CCC1=O)\\[C@@]1(C=C(C(=O)C[C@H]1O)NC(=O)/C=C/C(C)CC(C)CCCC)O",
#                              scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchAmide_AlkylaminoChain <- function(smi, scriptPath){
  amide_position <- .searchAmide(smi = smi,
                                 scriptPath = scriptPath)
  if(length(amide_position) != 0){
    posi <- lapply(amide_position, function(x) {
      .TraverseMolecule(smi = smi,
                        start_atom_idx = x[3],
                        non_traversable_atom_idx = x[-3],
                        scriptPath = scriptPath)[-1]
    })
    posi <- posi[sapply(posi, function(x) {
      if(length(x) == 0) return(FALSE)
      else return(TRUE)
    })]
  }else posi <- list()
  return(posi)
}

# 3. Thioester R1-C(=O)S-R2
# (1) Search thioester group
# .searchThioester(smi = "C(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)CCC(=O)SCCNC(=O)CCNC([C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@@H]1[C@@H](OP(=O)(O)O)[C@@H](O)[C@@H](O1)N1C=NC2C(N)=NC=NC1=2)=O",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchThioester <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CX3;$(C=O);$(C-S)](=O)S",
                       scriptPath = scriptPath)[[1]]
}

# (2) Search acyl chain of thioester group
# .searchThioester_AcylChain(smi = "C(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)CCC(=O)SCCNC(=O)CCNC([C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@@H]1[C@@H](OP(=O)(O)O)[C@@H](O)[C@@H](O1)N1C=NC2C(N)=NC=NC1=2)=O",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchThioester_AcylChain <- function(smi, scriptPath){
  thioester_position <- .searchThioester(smi = smi,
                                         scriptPath = scriptPath)
  lapply(thioester_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# (3) Search thioalkyl chain of thioester group
# .searchThioester_ThioalkylChain(smi = "C(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)CCC(=O)SCCNC(=O)CCNC([C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@@H]1[C@@H](OP(=O)(O)O)[C@@H](O)[C@@H](O1)N1C=NC2C(N)=NC=NC1=2)=O",
#                                 scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchThioester_ThioalkylChain <- function(smi, scriptPath){
  thioester_position <- .searchThioester(smi = smi,
                                         scriptPath = scriptPath)
  if(length(thioester_position) != 0){
    posi <- lapply(thioester_position, function(x) {
      .TraverseMolecule(smi = smi,
                        start_atom_idx = x[3],
                        non_traversable_atom_idx = x[-3],
                        scriptPath = scriptPath)[-1]
    })
    posi <- posi[sapply(posi, function(x) {
      if(length(x) == 0) return(FALSE)
      else return(TRUE)
    })]
  }else posi <- list()
  return(posi)
}

# 4. Fatty alcohols R-C-OH
# (1) Search fatty alcohols group
# .searchFattyAlcohols(smi = "OC=CCC(/C#CC#CC(O)/C=C/CCCCCCC)=C\\CO",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchFattyAlcohols(smi = "OC=CCC(/C#CC#CC(O)/C=C/CCCCCCC)=C\\C(=O)O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchFattyAlcohols <- function(smi, scriptPath){
  acyloxy_position <- .searchAcyloxy(smi = smi,
                                     scriptPath = scriptPath)
  amide_position <- .searchAmide(smi = smi,
                                 scriptPath = scriptPath)
  thioester_position <- .searchThioester(smi = smi,
                                         scriptPath = scriptPath)
  if(length(acyloxy_position) == 0 & length(amide_position) == 0 & length(thioester_position) == 0){
    matchRes <- .GetSubstructMatches(smis = smi,
                                     SMARTS = "[C;$([CH2,$([CH;$(C=C)])]);!$(C=O);$(C-O)]-[OH]",
                                     scriptPath = scriptPath)[[1]]
  }else{
    matchRes <- list()
  }
  return(matchRes)
}

# (2) Search fatty alcohol chain
# .searchFattyAlcohols_Chain(smi = "OC=CCC(/C#CC#CC(O)/C=C/CCCCCCC)=C\\CO",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchFattyAlcohols_Chain(smi = "C(CC(C)CCCCCCCCCC(=O)O)(=O)O",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchFattyAlcohols_Chain <- function(smi, scriptPath){
  fattyAlcohols_position <- .searchFattyAlcohols(smi = smi,
                                                 scriptPath = scriptPath)
  lapply(fattyAlcohols_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# 5. Fatty aldehydes R-C(=O)H
# (1) Search fatty aldehydes group
# .searchFattyAldehydes(smi = "C(CCC/C=C\\C=C/CCCC)(=O)[H]",
#                       scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchFattyAldehydes <- function(smi, scriptPath){
  acyloxy_position <- .searchAcyloxy(smi = smi,
                                     scriptPath = scriptPath)
  amide_position <- .searchAmide(smi = smi,
                                 scriptPath = scriptPath)
  thioester_position <- .searchThioester(smi = smi,
                                         scriptPath = scriptPath)
  if(length(acyloxy_position) == 0 & length(amide_position) == 0 & length(thioester_position) == 0){
    matchRes <- .GetSubstructMatches(smis = smi,
                                     SMARTS = "[CX3;CH;$(C=O);!$(C-O)]=O",
                                     scriptPath = scriptPath)[[1]]
  }else{
    matchRes <- list()
  }
  return(matchRes)
}

# (2) Search fatty aldehydes chain
# .searchFattyAldehydes_Chain(smi = "C(CCC/C=C\\C=C/CCCC)(=O)[H]",
#                             scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchFattyAldehydes_Chain <- function(smi, scriptPath){
  fattyAldehydes_position <- .searchFattyAldehydes(smi = smi,
                                                   scriptPath = scriptPath)
  lapply(fattyAldehydes_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# 6. Fatty nitriles R-C#N
# (1) Search fatty nitriles group
# .searchFattyNitriles(smi = "N#CCCCCC#CC#CC#CCCCCC#N",
#                     scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchFattyNitriles <- function(smi, scriptPath){
  acyloxy_position <- .searchAcyloxy(smi = smi,
                                     scriptPath = scriptPath)
  amide_position <- .searchAmide(smi = smi,
                                 scriptPath = scriptPath)
  thioester_position <- .searchThioester(smi = smi,
                                         scriptPath = scriptPath)
  if(length(acyloxy_position) == 0 & length(amide_position) == 0 & length(thioester_position) == 0){
    matchRes <- .GetSubstructMatches(smis = smi,
                                     SMARTS = "[CX2;$(C#N)]#N",
                                     scriptPath = scriptPath)[[1]]
  }else{
    matchRes <- list()
  }
  return(matchRes)
}

# (2) Search fatty nitriles chain
# .searchFattyNitriles_Chain(smi = "N#CCCCCC#CC#CC#CCCCCC#N",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchFattyNitriles_Chain <- function(smi, scriptPath){
  fattyNitries_position <- .searchFattyNitriles(smi = smi,
                                               scriptPath = scriptPath)
  lapply(fattyNitries_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# 7. Fatty ethers R1-O-R2
# (1) Search fatty ethers group
# .searchFattyEthers(smi = "C(CCCCCCC/C=C\\C=C\\O/C=C/C=C\\CC)",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchFattyEthers <- function(smi, scriptPath){
  acyloxy_position <- .searchAcyloxy(smi = smi,
                                     scriptPath = scriptPath)
  amide_position <- .searchAmide(smi = smi,
                                 scriptPath = scriptPath)
  thioester_position <- .searchThioester(smi = smi,
                                         scriptPath = scriptPath)
  if(length(acyloxy_position) == 0 & length(amide_position) == 0 & length(thioester_position) == 0){
    matchRes <- .GetSubstructMatches(smis = smi,
                                     SMARTS = "C-O-C",
                                     scriptPath = scriptPath)[[1]]
  }else{
    matchRes <- list()
  }
  return(matchRes)
}

# (2) Search fatty ethers chain
# .searchFattyEthers_chain(smi = "C(CCCCCCC/C=C\\C=C\\O/C=C/C=C\\CC)",
#                          scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchFattyEthers_chain <- function(smi, scriptPath){
  fattyEthers_position <- .searchFattyEthers(smi = smi,
                                             scriptPath = scriptPath)
  if(length(fattyEthers_position) != 0){
    posi <- lapply(fattyEthers_position, function(x) {
      .TraverseMolecule(smi = smi,
                        start_atom_idx = x[2],
                        non_traversable_atom_idx = c(),
                        scriptPath = scriptPath)[-1]
    })
    posi <- posi[sapply(posi, function(x) {
      if(length(x) == 0) return(FALSE)
      else return(TRUE)
    })]
  }else posi <- list()
  return(posi)
}

# 8. Hydrocarbons
# (1) Search hydrocarbons
# .searchHydrocarbons(smi = "CCCC(C)CCCC(C)CCCC(C)CCCC(C)CCCCCCCCCCCCCCCC",
#                     minium_C = 6,
#                     scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchHydrocarbons(smi = "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)O",
#                     minium_C = 6,
#                     scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchHydrocarbons <- function(smi, minium_C = 6, scriptPath){
  acyloxy_position <- .searchAcyloxy(smi = smi,
                                     scriptPath = scriptPath)
  amide_position <- .searchAmide(smi = smi,
                                 scriptPath = scriptPath)
  thioester_position <- .searchThioester(smi = smi,
                                         scriptPath = scriptPath)
  fattyAlcohols_position <- .searchFattyAlcohols(smi = smi,
                                                 scriptPath = scriptPath)
  fattyAldehydes_position <- .searchFattyAldehydes(smi = smi,
                                                   scriptPath = scriptPath)
  fattyNitries_position <- .searchFattyAldehydes(smi = smi,
                                                 scriptPath = scriptPath)
  fattyEthers_position <- .searchFattyEthers(smi = smi,
                                             scriptPath = scriptPath)
  groupNum <- sum(length(amide_position), length(thioester_position),
                  length(acyloxy_position), length(fattyAlcohols_position),
                  length(fattyAldehydes_position), length(fattyNitries_position),
                  length(fattyEthers_position))
  pattern <- rep("C", minium_C)
  pattern <- paste0(pattern, collapse = "~")
  if(groupNum == 0){
    matchRes <- .GetSubstructMatches(smis = smi,
                                     SMARTS = pattern,
                                     scriptPath = scriptPath)[[1]]
    if(length(matchRes) != 0){
      hydrocarbons_position <- .TraverseMolecule(smi, start_atom_idx = 0,
                                                 non_traversable_atom_idx = c(),
                                                 scriptPath = scriptPath)
    }
    else{
      hydrocarbons_position <- list()
    }
  }else{
    hydrocarbons_position <- list()
  }
  return(hydrocarbons_position)
}

# Search FA Main Chains
# .searchFA_MainChains(smi = "[C@@H]1([C@H](O)[C@H](OP(=O)(O)O)[C@@H](COP(O)(=O)OP(O)(=O)OCC(C)([C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C/C=C\\CC(O)=O)C)O1)N1C=NC2C(N)=NC=NC1=2",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchFA_MainChains <- function(smi, minium_C = 6, scriptPath){
  # Acyloxy
  acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  names(acyloxy_position) <- rep("acyloxy", length(acyloxy_position))
  acyloxy_acylChain_position <- .searchAcyloxy_AcylChain(smi = smi, scriptPath = scriptPath)
  names(acyloxy_acylChain_position) <- rep("acyloxy_acylChain", length(acyloxy_acylChain_position))
  acyloxy_alkoxyChain_position <- .searchAcyloxy_AlkoxyChain(smi = smi, scriptPath = scriptPath)
  names(acyloxy_alkoxyChain_position) <- rep("acyloxy_alkoxyChain", length(acyloxy_alkoxyChain_position))
  # Amide
  amide_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
  names(amide_position) <- rep("amide", length(amide_position))
  amide_acylChain_position <- .searchAmide_AcylChain(smi = smi, scriptPath = scriptPath)
  names(amide_acylChain_position) <- rep("amide_acylChain", length(amide_acylChain_position))
  amide_alkylaminoChain_position <- .searchAmide_AlkylaminoChain(smi = smi, scriptPath = scriptPath)
  names(amide_alkylaminoChain_position) <- rep("amide_alkylaminoChain", length(amide_alkylaminoChain_position))
  # Thioester
  thioester_position <- .searchThioester(smi = smi, scriptPath = scriptPath)
  names(thioester_position) <- rep("thioester", length(thioester_position))
  thioester_acylChain_position <- .searchThioester_AcylChain(smi = smi, scriptPath = scriptPath)
  names(thioester_acylChain_position) <- rep("thioester_acylChain", length(thioester_acylChain_position))
  thioester_thioalkylChain_position <- .searchThioester_ThioalkylChain(smi = smi, scriptPath = scriptPath)
  names(thioester_thioalkylChain_position) <- rep("thioester_thioalkylChain", length(thioester_thioalkylChain_position))
  # Fatty Alcohol
  fattyAlcohol_position <- .searchFattyAlcohols(smi = smi, scriptPath = scriptPath)
  names(fattyAlcohol_position) <- rep("fattyAlcohol", length(fattyAlcohol_position))
  fattyAlcoholChain_position <- .searchFattyAlcohols_Chain(smi = smi, scriptPath = scriptPath)
  names(fattyAlcoholChain_position) <- rep("fattyAlcoholChain", length(fattyAlcoholChain_position))
  # Fatty Aldehyde
  fattyAldehyde_position <- .searchFattyAldehydes(smi = smi, scriptPath = scriptPath)
  names(fattyAldehyde_position) <- rep("fattyAldehyde", length(fattyAldehyde_position))
  fattyAldehydeChain_position <- .searchFattyAldehydes_Chain(smi = smi, scriptPath = scriptPath)
  names(fattyAldehydeChain_position) <- rep("fattyAldehydeChain", length(fattyAldehydeChain_position))
  # Fatty Nitrile
  fattyNitrile_position <- .searchFattyNitriles(smi = smi, scriptPath = scriptPath)
  names(fattyNitrile_position) <- rep("fattyNitrile", length(fattyNitrile_position))
  fattyNitrileChain_position <- .searchFattyNitriles_Chain(smi = smi, scriptPath = scriptPath)
  names(fattyNitrileChain_position) <- rep("fattyNitrileChain", length(fattyNitrileChain_position))
  # Fatty Ether
  fattyEther_position <- .searchFattyEthers(smi = smi, scriptPath = scriptPath)
  names(fattyEther_position) <- rep("fattyEther", length(fattyEther_position))
  fattyEtherChain_position <- .searchFattyEthers_chain(smi = smi, scriptPath = scriptPath)
  names(fattyEtherChain_position) <- rep("fattyEtherChain", length(fattyEtherChain_position))
  # hydrocarbons
  hydrocarbons_position <- .searchHydrocarbons(smi = smi, scriptPath = scriptPath)
  names(hydrocarbons_position) <- rep("hydrocarbons", length(hydrocarbons_position))

  c(acyloxy_position, acyloxy_acylChain_position, acyloxy_alkoxyChain_position,
    amide_position, amide_acylChain_position, amide_alkylaminoChain_position,
    thioester_position, thioester_acylChain_position, thioester_thioalkylChain_position,
    fattyAlcohol_position, fattyAlcoholChain_position,
    fattyAldehyde_position, fattyAldehydeChain_position,
    fattyNitrile_position, fattyNitrileChain_position,
    fattyEther_position, fattyEtherChain_position,
    hydrocarbons_position)
}
