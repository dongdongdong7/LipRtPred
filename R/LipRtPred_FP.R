# Calculate LipRtPred Fingerprints
# 250319
# Barry Song

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py
# minimumohNum: minimum number of OH to be taken into account when calculating the OH distance

.calLipRtPred_FP(smi = "C1C[C@H](O)C(C)(C)[C@]2([H])CCC3[C@]4(C)CC[C@]([H])([C@H](C)CCC(=C)C(C)C)[C@@]4(C)CCC=3[C@@]12C",
                 scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.calLipRtPred_FP <- function(smi, scriptPath, minimumohNum = 6){
  # Main Chains
  FA_MainChains <- .searchFA_MainChains(smi = smi, scriptPath = scriptPath)
  GL_MainChains <- .searchGL_MainChains(smi = smi, scriptPath = scriptPath)
  GP_MainChains <- .searchGP_MainChains(smi = smi, scriptPath = scriptPath)
  SP_MainChains <- .searchSP_MainChains(smi = smi, scriptPath = scriptPath)
  ST_MainChains <- .searchST_MainChains(smi = smi, scriptPath = scriptPath)

  main_C_chains <- c(FA_MainChains, GL_MainChains, GP_MainChains, ST_MainChains[names(ST_MainChains) == "sterylEster_chain"], SP_MainChains[names(SP_MainChains) == "sphingosineN_chain"], SP_MainChains[names(SP_MainChains) == "spisulosineN_chain"], SP_MainChains[names(SP_MainChains) == "SP_other_chain"])
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

  # 水苏碱骨架及其变体
  if(length(spisulosineC_chain) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = spisulosineC_chain, scriptPath = scriptPath)
    spiC <- length(which(symbols == "C"))
    # 羟基与位置
    if(length(oh_position) != 0){
      spioh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% spisulosineC_chain)) return(TRUE)
        else return(FALSE)
      })]
      spiohNum <- length(spioh_position)
      if(length(spiohNum) != 0){
        spiohPos <- sapply(spioh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = spisulosineC_chain[1], scriptPath = scriptPath))
        })
        spiohPos <- spiohPos[1:minimumohNum]
        names(spiohPos) <- paste0("spiohPos", 1:length(spiohPos))
        spiohPos[is.na(spiohPos)] <- 0
      }else{
        spiohPos <- rep(0, minimumohNum)
        names(spiohPos) <- paste0("spiohPos", 1:length(spiohPos))
      }
    }else{
      spiohNum <- 0
      spiohPos <- rep(0, minimumohNum)
      names(spiohPos) <- paste0("spiohPos", 1:length(spiohPos))
    }
    # 双键
    if(length(db_position) != 0){
      spidb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% spisulosineC_chain)) return(TRUE)
        else return(FALSE)
      })]
      spidbNum <- length(spidb_position)
    }else{
      spidbNum <- 0
    }
    # 酮
    if(length(dbo_position) != 0){
      spidbo_position <- dbo_position[sapply(dbo_position, function(x) {
        if(all(x %in% spisulosineC_chain)) return(TRUE)
        else return(FALSE)
      })]
      spidboNum <- length(spidbo_position)
    }else{
      spidboNum <- 0
    }
  }

  # 类固醇骨架变体
  if(length(steroidSkeleton_derivative) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = steroidSkeleton_derivative, scriptPath = scriptPath)
    steroidSkDC <- length(symbols == "C")
    # 羟基与位置
    if(length(oh_position)){
      steroidSkDoh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% steroidSkeleton_derivative)) return(TRUE)
        else return(FALSE)
      })]
      steroidSkDohNum <- length(steroidSkDoh_position)
      if(steroidSkDohNum != 0){
        steroidSkDohPos <- sapply(steroidSkDoh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = steroidSkeleton_derivative[1], scriptPath = scriptPath))
        })
        steroidSkDohPos <- steroidSkDohPos[1:minimumohNum]
        names(steroidSkDohPos) <- paste0("steroidSkDohPos", 1:length(steroidSkDohPos))
        steroidSkDohPos[is.na(steroidSkDohPos)] <- 0
      }else{
        steroidSkDohPos <- rep(0, minimumohNum)
        names(steroidSkDohPos) <- paste0("steroidSkDohPos", 1:length(steroidSkDohPos))
      }
    }else{
      steroidSkDohNum <- 0
      steroidSkDohPos <- rep(0, minimumohNum)
      names(steroidSkDohPos) <- paste0("steroidSkDohPos", 1:length(steroidSkDohPos))
    }
    # 双键
    if(length(db_position) != 0){
      steroidSkDdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% steroidSkeleton_derivative)) return(TRUE)
        else return(FALSE)
      })]
      steroidSkDdbNum <- length(steroidSkDdb_position)
    }else{
      steroidSkDdbNum <- 0
    }
    # 酮
    if(length(dbo_position) != 0){
      steroidSkDdbo_position <- dbo_position[sapply(dbo_position, function(x) {
        if(all(x %in% steroidSkeleton_derivative)) return(TRUE)
        else return(FALSE)
      })]
      steroidSkDdboNum <- length(steroidSkDdbo_position)
    }else{
      steroidSkDdboNum <- 0
    }
  }

  # 类固醇碳链
  if(length(steroidSkeleton_chain) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = steroidSkeleton_chain, scriptPath = scriptPath)
    steroidSkCC <- length(which(symbols == "C"))
    # 羟基与位置
    if(length(oh_position)){
      steroidSkCoh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% steroidSkeleton_chain)) return(TRUE)
        else return(FALSE)
      })]
      steroidSkCohNum <- length(steroidSkCoh_position)
      if(steroidSkCohNum != 0){
        steroidSkCohPos <- sapply(steroidSkCoh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = steroidSkeleton_chain[1], scriptPath = scriptPath))
        })
        steroidSkCohPos <- steroidSkCohPos[1:minimumohNum]
        names(steroidSkCohPos) <- paste0("steroidSkCohPos", 1:length(steroidSkCohPos))
        steroidSkCohPos[is.na(steroidSkCohPos)] <- 0
      }else{
        steroidSkCohPos <- rep(0, minimumohNum)
        names(steroidSkCohPos) <- paste0("steroidSkCohPos", 1:length(steroidSkCohPos))
      }
    }else{
      steroidSkCohNum <- 0
      steroidSkCohPos <- rep(0, minimumohNum)
      names(steroidSkCohPos) <- paste0("steroidSkCohPos", 1:length(steroidSkCohPos))
    }
    # 双键
    if(length(db_position) != 0){
      steroidSkCdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% steroidSkeleton_chain)) return(TRUE)
        else return(FALSE)
      })]
      steroidSkCdbNum <- length(steroidSkCdb_position)
    }else{
      steroidSkCdbNum <- 0
    }
    # 酮
    if(length(dbo_position) != 0){
      steroidSkCdbo_position <- dbo_position[sapply(dbo_position, function(x) {
        if(all(x %in% steroidSkeleton_chain)) return(TRUE)
        else return(FALSE)
      })]
      steroidSkCdboNum <- length(steroidSkCdbo_position)
    }else{
      steroidSkCdboNum <- 0
    }
  }


  browser()

  FA_c <- sum(sapply(FA_MainChains, function(x) {
    x_symbol <- .GetAtomSymbol(smi = smi, atom_idx_vector = x, scriptPath = scriptPath)
    length(which(x_symbol == "C"))
  }))
}
