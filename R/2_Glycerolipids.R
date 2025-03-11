# Search substructure of Glycerolipids (GL)
# Barry Song
# 250307

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py

# 1. Glycero C(O)C(O)C(O)
# (1) Search glycerol
# .searchGlycerol(smi = "OC[C@]([H])(O)COC(CCCCCCCCCCC)=O",
#                 scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CH2;$(C-O)](O)-[CH;$(C-O)](O)-[CH2;$(C-O)](O)",
                       scriptPath = scriptPath)[[1]]
}

# (2) Search glycerides
# .searchGlycerides(smi = "OC[C@]([H])(O)COC(CCCCCCCCCCC)=O",
#                   scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerides <- function(smi, scriptPath){
  glycerol_position <- .searchGlycerol(smi = smi, scriptPath = scriptPath)
  O_sn1 <- sapply(glycerol_position, function(x) {x[2]})
  O_sn2 <- sapply(glycerol_position, function(x) {x[4]})
  O_sn3 <- sapply(glycerol_position, function(x) {x[6]})
  acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  sn_info <- sapply(acyloxy_position, function(x) {
    if(x[3] %in% O_sn1) return("sn1")
    else if(x[3] %in% O_sn2) return("sn2")
    else if(x[3] %in% O_sn3) return("sn3")
    else return("none")
  })
  acyloxy_position <- acyloxy_position[sn_info != "none"]
  sn_info <- sn_info[sn_info != "none"]
  names(acyloxy_position) <- sn_info
  return(acyloxy_position)
}

# (3) Search glycerol ethers
# .searchGlycerolEthers(smi = "OC[C@]([H])(O)CO/C=C\\C#C/C=C\\CCCCCCC(C)C",
#                       scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerolEthers <- function(smi, scriptPath){
  glycerol_position <- .searchGlycerol(smi = smi, scriptPath = scriptPath)
  O_sn1 <- sapply(glycerol_position, function(x) {x[2]})
  O_sn2 <- sapply(glycerol_position, function(x) {x[4]})
  O_sn3 <- sapply(glycerol_position, function(x) {x[6]})
  ethers_position <- .GetSubstructMatches(smis = smi,
                                          SMARTS = "C-O-C",
                                          scriptPath = scriptPath)[[1]]
  sn_info <- sapply(ethers_position, function(x) {
    if(x[2] %in% O_sn1) return("sn1")
    else if(x[2] %in% O_sn2) return("sn2")
    else if(x[2] %in% O_sn3) return("sn3")
    else return("none")
  })
  ethers_position <- ethers_position[sn_info != "none"]
  sn_info <- sn_info[sn_info != "none"]
  names(ethers_position) <- sn_info
  return(ethers_position)
}

# 2. Glycosyl
# (1) Search glycosyl
# .searchGlycosyl(smi = "C(O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1)[C@]([H])(OC(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)=O)COC(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)=O",
#                 scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycosyl <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)C1C(O)C(O)C(O)C(O)O1",
                       scriptPath = scriptPath)[[1]]
}

# 3. Betaine
# (1) Search betaine
# .searchBetaine(smi = "[C@](COCCC(C(=O)[O-])[N+](C)(C)C)([H])(O)COC(CCCCCCCCCCCCCCCCC)=O",
#                scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchBetaine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[N+](C)(C)(C)CC(=O)[O-]",
                       scriptPath = scriptPath)[[1]]
}
