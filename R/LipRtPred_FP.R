# Calculate LipRtPred Fingerprints
# 250319
# Barry Song

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py
# minimumohNum: minimum number of OH to be taken into account when calculating the OH distance

.calLipRtPred_FP(smi = "OC[C@]([H])(NC(=O)CCCCCCCCCCCCCCC)C(=O)CCCCCCCCCCCCCCC",
                 scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.calLipRtPred_FP <- function(smi, scriptPath, minimumohNum = 6){
  # Main Chains
  FA_MainChains <- .searchFA_MainChains(smi = smi, scriptPath = scriptPath)
  GL_MainChains <- .searchGL_MainChains(smi = smi, scriptPath = scriptPath)
  GP_MainChains <- .searchGP_MainChains(smi = smi, scriptPath = scriptPath)
  SP_MainChains <- .searchSP_MainChains(smi = smi, scriptPath = scriptPath)
  ST_MainChains <- .searchST_MainChains(smi = smi, scriptPath = scriptPath)

  main_C_chains <- c(FA_MainChains, GL_MainChains, GP_MainChains, ST_MainChains[names(ST_MainChains) == "sterylEster_chain"], SP_MainChains[names(SP_MainChains) == "sphingosineN_chain"], SP_MainChains[names(SP_MainChains) == "SP_other_chain"])
  sphingosineC_chain <- as.integer(unlist(SP_MainChains[names(SP_MainChains) == "sphingosineC_chain"]))
  spisulosineC_chain <- as.integer(unlist(SP_MainChains[names(SP_MainChains) == "spisulosineC_chain"]))
  steroidSkeleton_chain <- as.integer(unlist(ST_MainChains[names(ST_MainChains) == "steroidSkeleton_chain"]))
  steroidSkeleton_derivative <- as.integer(unlist(ST_MainChains[names(ST_MainChains) == "steroidSkeleton_derivative"]))
  secosteroidSkeleton_chain <- as.integer(unlist(ST_MainChains[names(ST_MainChains) == "secosteroidSkeleton_chain"]))
  secosteroidSkeleton_derivative <- as.integer(unlist(ST_MainChains[names(ST_MainChains) == "secosteroidSkeleton_derivative"]))

  # 一般基团
  # 双键
  db_position <- .GetSubstructMatches(smis = smi,
                                      SMARTS = "C=C",
                                      scriptPath = scriptPath)[[1]]
  if(length(db_position) != 0){
    dbNum <- length(which(sapply(db_position, function(x) {
      if(all(x %in% unlist(c(main_C_chains)))) return(TRUE)
      else return(FALSE)
    })))
  }else dbNum <- 0
  # 羟基
  oh_position <- .GetSubstructMatches(smis = smi,
                                      SMARTS = "[C;!$(C=O);$(C-[OH])]-[OH]",
                                      scriptPath = scriptPath)[[1]]
  if(length(oh_position) != 0){
    ohNum <- length(which(sapply(oh_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    })))
  }else ohNum <- 0
  # 酮
  dbo_position <- .GetSubstructMatches(smis = smi,
                                       SMARTS = "[C;$(C=O);!$(C-[N,O,S,P])]=[OX1]",
                                       scriptPath = scriptPath)[[1]]
  if(length(dbo_position) != 0){
    dboNum <- length(which(sapply(dbo_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    })))
  }else dboNum <- 0
  # 酰氧键
  coo_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  if(length(coo_position) != 0){
    cooNum <- length(which(sapply(coo_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    })))
  }else cooNum <- 0
  # 酰胺键
  con_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
  if(length(con_position) != 0){
    conNum <- which(sapply(con_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    }))
  }else conNum <- 0
  # 硫代酰
  cos_position <- .searchThioester(smi = smi, scriptPath = scriptPath)
  if(length(cos_position) != 0){
    cosNum <- which(sapply(cos_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    }))
  }else cosNum <- 0
  # 醚键
  coc_position <- .GetSubstructMatches(smis = smi,
                                       SMARTS = "[C;!$(C=O);$(C-O)]-[OX2;OH0]-[C;!$(C=O);$(C-O)]",
                                       scriptPath = scriptPath)[[1]]
  if(length(coc_position) != 0){
    cocNum <- which(sapply(coc_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    }))
  }else cocNum <- 0
  # 过氧化键
  oo_position <- .GetSubstructMatches(smis = smi,
                                      SMARTS = "[O]-[O]",
                                      scriptPath = scriptPath)[[1]]
  if(length(oo_position) != 0){
    ooNum <- which(sapply(oo_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    }))
  }else ooNum <- 0

  # 鞘氨醇骨架及其变体
  if(length(sphingosineC_chain) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = sphingosineC_chain, scriptPath = scriptPath)
    sphC <- length(which(symbols == "C"))
    # 羟基与位置
    if(length(oh_position) != 0){
      sphoh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% sphingosineC_chain)) return(TRUE)
        else return(FALSE)
      })]
      sphohNum <- length(sphoh_position)
      if(sphohNum != 0){
        sphohPos <- sapply(sphoh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = sphingosineC_chain[1], scriptPath = scriptPath))
        })
        sphohPos <- sphohPos[1:minimumohNum]
        names(sphohPos) <- paste0("sphohPos", 1:length(sphohPos))
        sphohPos[is.na(sphohPos)] <- 0
      }else{
        sphohPos <- rep(0, minimumohNum)
        names(sphohPos) <- paste0("sphohPos", 1:length(sphohPos))
      }
    }else{
      sphohNum <- 0
      sphohPos <- rep(0, minimumohNum)
      names(sphohPos) <- paste0("sphohPos", 1:length(sphohPos))
    }
    # 双键
    if(length(db_position) != 0){
      sphdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% sphingosineC_chain)) return(TRUE)
        else return(FALSE)
      })]
      sphdbNum <- length(sphdb_position)
    }else {
      sphdbNum <- 0
    }
    # 酮
    if(length(dbo_position) != 0){
      sphdbo_position <- dbo_position[sapply(dbo_position, function(x) {
        if(all(x %in% sphingosineC_chain)) return(TRUE)
        else return(FALSE)
      })]
      sphdboNum <- length(sphdbo_position)
    }else{
      sphdboNum <- 0
    }
  }
  browser()

  # 水苏碱骨架及其变体


  FA_c <- sum(sapply(FA_MainChains, function(x) {
    x_symbol <- .GetAtomSymbol(smi = smi, atom_idx_vector = x, scriptPath = scriptPath)
    length(which(x_symbol == "C"))
  }))
}
