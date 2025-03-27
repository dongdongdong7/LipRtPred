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
# .searchGlycerol(smi = "C(OC(CCCCCCCCCCCCCCC)=O)C(=O)COP(=O)(O)O",
#                 scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CH2;$(C-O)](O)-[CH&$(C-O),CH0&$(C=O)](~O)-[CH2;$(C-O)](O)",
                       scriptPath = scriptPath)[[1]]
}

# (2) Search glycerol ester
# .searchGlycerolEster(smi = "OC[C@]([H])(O)COC(CCCCCCCCCCC)=O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerolEster <- function(smi, scriptPath){
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

# (3) Search glycerol ester chain
# .searchGlycerolEster_Chain(smi = "OC[C@]([H])(O)COC(CCCCCCCCCCC)=O",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGlycerolEster_Chain(smi = "C(OC(CCCCCCCCCCCCCCC)=O)C(=O)COP(=O)(O)O",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerolEster_Chain <- function(smi, scriptPath){
  glycerolEster_position <- .searchGlycerolEster(smi = smi,
                                                 scriptPath = scriptPath)
  lapply(glycerolEster_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# (4) Search glycerol ether
# .searchGlycerolEther(smi = "OC[C@]([H])(O)CO/C=C\\C#C/C=C\\CCCCCCC(C)C",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGlycerolEther(smi = "[C@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)(COP(CC[N+](C)(C)C)(=O)[O-])COC(CCCCCCCCCCCCCCC)=O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGlycerolEther(smi = "C(=O)(COP(=O)(O)O)COCCCCCCCCCCCCCCCC",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerolEther <- function(smi, scriptPath){
  glycerol_position <- .searchGlycerol(smi = smi, scriptPath = scriptPath)
  O_sn1 <- sapply(glycerol_position, function(x) {x[2]})
  O_sn2 <- sapply(glycerol_position, function(x) {x[4]})
  O_sn3 <- sapply(glycerol_position, function(x) {x[6]})
  ethers_position <- .GetSubstructMatches(smis = smi,
                                          SMARTS = "[CH2,CH]-O-[C;!$(C=O)]",
                                          scriptPath = scriptPath)[[1]]
  ethers_position <- lapply(ethers_position, function(x) {
    if(all(x[c(1,2)] %in% unlist(glycerol_position))) return(x)
    else return(rev(x))
  })
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

# (5) Search glycerol ether chain
# .searchGlycerolEther_Chain(smi = "OC[C@]([H])(O)CO/C=C\\C#C/C=C\\CCCCCCC(C)C",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGlycerolEther_Chain(smi = "OC[C@]1(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H](CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@@]([H])(CO)COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC1)C)[H]",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGlycerolEther_Chain(smi = "OC[C@@](OCCCCCCCCC1CC2C3C4CCC4C3C2CC1)([H])COCCCCCCC1CC2C3C=CCCC3C2C=C1",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGlycerolEther_Chain(smi = "C(=O)(COP(=O)(O)O)COCCCCCCCCCCCCCCCC",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerolEther_Chain <- function(smi, scriptPath){
  glycerolEther_position <- .searchGlycerolEther(smi = smi,
                                                 scriptPath = scriptPath)
  ringsInfo <- .GetRingAtom(smi = smi, scriptPath = scriptPath)
  rings <- lapply(ringsInfo, function(ring) {
    lapply(glycerolEther_position, function(x) {
      if(all(x %in% ring)) return(x)
      else return(NULL)
    })
  })
  rings <- lapply(rings, function(ring) {
    purrr::compact(ring)
  })
  rings <- rings[!sapply(rings, function(ring) {length(ring) == 0})]
  if(length(rings) == 0){
    lapply(glycerolEther_position, function(x) {
      .TraverseMolecule(smi = smi,
                        start_atom_idx = x[3],
                        non_traversable_atom_idx = x[-3],
                        scriptPath = scriptPath)
    })
  }else{
    lapply(glycerolEther_position, function(x) {
      ring <- rings[[which(sapply(ringsInfo, function(y) {
        if(all(x %in% y)) return(TRUE)
        else return(FALSE)
      }))]]
      ring <- ring[!sapply(ring, function(y) {
        if(all(x %in% y)) return(TRUE)
        else return(FALSE)
      })]
      ring <- lapply(ring, function(y) {y[1:2]})
      .TraverseMolecule(smi = smi,
                        start_atom_idx = x[3],
                        non_traversable_atom_idx = c(x[-3], unlist(ring)),
                        scriptPath = scriptPath)
    })
  }
}

# Search GL Main Chains
# .searchGL_MainChains(smi = "OC[C@]1(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H](CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@@]([H])(CO)COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC1)C)[H]",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGL_MainChains <- function(smi, scriptPath){
  glycerol_position <- .searchGlycerol(smi = smi, scriptPath = scriptPath)
  names(glycerol_position) <- rep("glycerol", length(glycerol_position))
  glycerolEster_position <- .searchGlycerolEster(smi = smi, scriptPath = scriptPath)
  names(glycerolEster_position) <- rep("glycerolEster", length(glycerolEster_position))
  glycerolEther_position <- .searchGlycerolEther(smi = smi, scriptPath = scriptPath)
  names(glycerolEther_position) <- rep("glycerolEther", length(glycerolEther_position))
  glycerolEsterChain_position <- .searchGlycerolEster_Chain(smi = smi, scriptPath = scriptPath)
  glycerolEtherChain_position <- .searchGlycerolEther_Chain(smi = smi, scriptPath = scriptPath)

  main_chains <- c(glycerolEsterChain_position, glycerolEtherChain_position)
  if(length(main_chains) != 0){
    for(i in 1:length(main_chains)){
      x <- main_chains[[i]]
      if(all(x == "none")) next
      for(j in (1:length(main_chains))[-i]){
        y <- main_chains[[j]]
        if(all(y == "none")) next
        if(any(y %in% x)) main_chains[[j]] <- "none"
      }
    }
    main_chains <- lapply(main_chains, function(x){
      if(all(x == "none")) return(NULL)
      else return(x)
    })
    main_chains <- main_chains[!sapply(main_chains, is.null)]
  }

  c(main_chains, glycerol_position, glycerolEster_position, glycerolEther_position)
}
