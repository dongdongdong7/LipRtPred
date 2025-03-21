# Calculate LipRtPred Fingerprints
# 250319
# Barry Song

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py
# minimumohNum: minimum number of OH to be taken into account when calculating the OH distance

# .calLipRtPred_FP(smi = "C(OCCCCCCCCCCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C)=O)COC(=O)CCCCCCCCCCCC(C)C",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .calLipRtPred_FP(smi = "[C@](COP(=O)(O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1)([H])(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .calLipRtPred_FP(smi = "O(P(OC[C@@]1([H])OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@]([H])(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H](CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC1)C)CO)(O)=O)CCN",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .calLipRtPred_FP(smi = "C(OC(=O)CCCCCCCCCCCCCCC)[C@]([H])(NC(CCCCCCCCCCCCCCC)=O)[C@H](O)/C=C/CCCCCCCCCCCCC",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
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

  # 特征基团
  # 甘油
  glycerolNum <- length(.searchGlycerol(smi = smi, scriptPath = scriptPath))
  names(glycerolNum) <- "glycerolNum"
  # 二羟基丙酮
  oh2dboNum <- length(.searchDihydroxyacetone(smi = smi, scriptPath = scriptPath))
  names(oh2dboNum) <- "oh2dboNum"
  # 磷酸
  phosphateNum <- length(.searchPhosphate(smi = smi, scriptPath = scriptPath))
  names(phosphateNum) <- "phosphateNum"
  # 甜菜碱
  betaineNum <- length(.searchBetaine(smi = smi, scriptPath = scriptPath))
  names(betaineNum) <- "betaineNum"
  # 五碳糖
  pentose_position <- .searchPentose(smi = smi, scriptPath = scriptPath)
  pentoseNum <- length(pentose_position)
  # 六碳糖
  hexose_position <- .searchHexose(smi = smi, scriptPath = scriptPath)
  hexoseNum <- length(hexose_position)
  carbohydrate_position <- unlist(c(pentose_position, hexose_position))
  main_C_chains <- lapply(main_C_chains, function(x) {
    x[!x %in% carbohydrate_position]
  })
  carbohydrateNum <- sum(c(pentoseNum, hexoseNum))
  names(carbohydrateNum) <- "carbohydrateNum"
  # 胆碱
  cholineNum <- length(.searchCholine(smi = smi, scriptPath = scriptPath))
  names(cholineNum) <- "cholineNum"
  # 乙醇胺
  ethanolamineNum <- length(.searchEthanolamine(smi = smi, scriptPath = scriptPath))
  names(ethanolamineNum) <- "ethanolamineNum"
  # 丝氨酸
  serineNum <- length(.searchSerine(smi = smi, scriptPath = scriptPath))
  names(serineNum) <- "serineNum"
  # 肌醇
  inositolNum <- length(.searchInositol(smi = smi, scriptPath = scriptPath))
  names(inositolNum) <- "inositolNum"
  # 乙醇
  ethanolNum <- length(.searchEthanol(smi = smi, scriptPath = scriptPath))
  names(ethanolNum) <- "ethanolNum"
  # 苏氨酸
  threonineNum <- length(.searchThreonine(smi = smi, scriptPath = scriptPath))
  names(threonineNum) <- "threonineNum"
  # 肉碱
  carnitineNum <- length(.searchCarnitine(smi = smi, scriptPath = scriptPath))
  names(carnitineNum) <- "carnitineNum"
  # 磺酰基
  sulfonylNum <- length(.searchSulfonyl(smi = smi, scriptPath = scriptPath))
  names(sulfonylNum) <- "sulfonylNum"


  # 一般基团
  # 羟基
  ohNum <- 0
  ohPos <- rep(0, minimumohNum)
  oh_position <- .GetSubstructMatches(smis = smi,
                                      SMARTS = "[C;!$(C=O);$(C-[OH])]-[OH]",
                                      scriptPath = scriptPath)[[1]]
  if(length(oh_position) != 0){
    ohNum <- length(which(sapply(oh_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    })))
    if(length(main_C_chains) != 0 & ohNum != 0){
      ohPos <- sapply(oh_position, function(x) {
        idx <- which(sapply(main_C_chains, function(y) {
          if(all(x %in% y)) return(TRUE)
          else return(FALSE)
        }))
        if(length(idx) == 0) return(NA)
        nmc <- main_C_chains[[idx]]
        length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = nmc[1], scriptPath = scriptPath))
      })
      ohPos <- sort(ohPos)
      ohPos <- ohPos[1:minimumohNum]
      ohPos[is.na(ohPos)] <- 0
    }
  }
  names(ohNum) <- "ohNum"
  names(ohPos) <- paste0("ohPos", 1:length(ohPos))
  # 双键
  dbNum <- 0
  db_position <- .GetSubstructMatches(smis = smi,
                                      SMARTS = "C=C",
                                      scriptPath = scriptPath)[[1]]
  if(length(db_position) != 0){
    dbNum <- length(which(sapply(db_position, function(x) {
      if(all(x %in% unlist(c(main_C_chains)))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(dbNum) <- "dbNum"
  # 酮
  dboNum <- 0
  dbo_position <- .GetSubstructMatches(smis = smi,
                                       SMARTS = "[C;$(C=O);!$(C-[N,O,S,P])]=[OX1]",
                                       scriptPath = scriptPath)[[1]]
  if(length(dbo_position) != 0){
    dboNum <- length(which(sapply(dbo_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(dboNum) <- "dboNum"
  # 酰氧键
  cooNum <- 0
  coo_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  if(length(coo_position) != 0){
    cooNum <- length(which(sapply(coo_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(cooNum) <- "cooNum"
  # 酰胺键
  conNum <- 0
  con_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
  if(length(con_position) != 0){
    conNum <- length(which(sapply(con_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(conNum) <- "conNum"
  # 硫代酰
  cosNum <- 0
  cos_position <- .searchThioester(smi = smi, scriptPath = scriptPath)
  if(length(cos_position) != 0){
    cosNum <- length(which(sapply(cos_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(cosNum) <- "cosNum"
  # 醚键
  cocNum <- 0
  coc_position <- .GetSubstructMatches(smis = smi,
                                       SMARTS = "[C;!$(C=O);$(C-O)]-[OX2;OH0]-[C;!$(C=O);$(C-O)]",
                                       scriptPath = scriptPath)[[1]]
  if(length(coc_position) != 0){
    cocNum <- length(which(sapply(coc_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(cocNum) <- "cocNum"
  # 过氧化键
  ooNum <- 0
  oo_position <- .GetSubstructMatches(smis = smi,
                                      SMARTS = "[O]-[O]",
                                      scriptPath = scriptPath)[[1]]
  if(length(oo_position) != 0){
    ooNum <- length(which(sapply(oo_position, function(x) {
      if(all(x %in% unlist(main_C_chains))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(ooNum) <- "ooNum"

  # 鞘氨醇骨架及其变体
  sphC <- 0
  sphohNum <- 0
  sphohPos <- rep(0, minimumohNum)
  sphdbNum <- 0
  sphdboNum <- 0
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
        sphohPos <- sort(sphohPos)
        sphohPos <- sphohPos[1:minimumohNum]
        sphohPos[is.na(sphohPos)] <- 0
      }
    }
    # 双键
    if(length(db_position) != 0){
      sphdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% sphingosineC_chain)) return(TRUE)
        else return(FALSE)
      })]
      sphdbNum <- length(sphdb_position)
    }
    # 酮
    if(length(dbo_position) != 0){
      sphdbo_position <- dbo_position[sapply(dbo_position, function(x) {
        if(all(x %in% sphingosineC_chain)) return(TRUE)
        else return(FALSE)
      })]
      sphdboNum <- length(sphdbo_position)
    }
  }
  names(sphC) <- "sphC"
  names(sphohNum) <- "sphohNum"
  names(sphohPos) <- paste0("sphohPos", 1:length(sphohPos))
  names(sphdbNum) <- "sphdbNum"
  names(sphdboNum) <- "sphdboNum"

  # 水苏碱骨架及其变体
  spiC <- 0
  spiohNum <- 0
  spiohPos <- rep(0, minimumohNum)
  spidbNum <- 0
  spidboNum <- 0
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
        spiohPos <- sort(spiohPos)
        spiohPos <- spiohPos[1:minimumohNum]
        spiohPos[is.na(spiohPos)] <- 0
      }
    }
    # 双键
    if(length(db_position) != 0){
      spidb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% spisulosineC_chain)) return(TRUE)
        else return(FALSE)
      })]
      spidbNum <- length(spidb_position)
    }
    # 酮
    if(length(dbo_position) != 0){
      spidbo_position <- dbo_position[sapply(dbo_position, function(x) {
        if(all(x %in% spisulosineC_chain)) return(TRUE)
        else return(FALSE)
      })]
      spidboNum <- length(spidbo_position)
    }
  }
  names(spiC) <- "spiC"
  names(spiohNum) <- "spiohNum"
  names(spiohPos) <- paste0("spiohPos", 1:length(spiohPos))
  names(spidbNum) <- "spidbNum"
  names(spidboNum) <- "spidboNum"

  # 类固醇骨架变体
  steroidSkDC <- 0
  steroidSkDohNum <- 0
  steroidSkDohPos <- rep(0, minimumohNum)
  steroidSkDdbNum <- 0
  steroidSkDdboNum <- 0
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
        steroidSkDohPos <- sort(steroidSkDohPos)
        steroidSkDohPos <- steroidSkDohPos[1:minimumohNum]
        steroidSkDohPos[is.na(steroidSkDohPos)] <- 0
      }
    }
    # 双键
    if(length(db_position) != 0){
      steroidSkDdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% steroidSkeleton_derivative)) return(TRUE)
        else return(FALSE)
      })]
      steroidSkDdbNum <- length(steroidSkDdb_position)
    }
    # 酮
    if(length(dbo_position) != 0){
      steroidSkDdbo_position <- dbo_position[sapply(dbo_position, function(x) {
        if(all(x %in% steroidSkeleton_derivative)) return(TRUE)
        else return(FALSE)
      })]
      steroidSkDdboNum <- length(steroidSkDdbo_position)
    }
  }
  names(steroidSkDC) <- "steroidSkDC"
  names(steroidSkDohNum) <- "steroidSkDohNum"
  names(steroidSkDohPos) <- paste0("steroidSkDohPos", 1:length(steroidSkDohPos))
  names(steroidSkDdbNum) <- "steroidSkDdbNum"
  names(steroidSkDdboNum) <- "steroidSkDdboNum"

  # 类固醇骨架碳链
  steroidSkCC <- 0
  steroidSkCohNum <- 0
  steroidSkCohPos <- rep(0, minimumohNum)
  steroidSkCdbNum <- 0
  steroidSkCdboNum <- 0
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
        steroidSkCohPos <- sort(steroidSkCohPos)
        steroidSkCohPos <- steroidSkCohPos[1:minimumohNum]
        steroidSkCohPos[is.na(steroidSkCohPos)] <- 0
      }
    }
    # 双键
    if(length(db_position) != 0){
      steroidSkCdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% steroidSkeleton_chain)) return(TRUE)
        else return(FALSE)
      })]
      steroidSkCdbNum <- length(steroidSkCdb_position)
    }
    # 酮
    if(length(dbo_position) != 0){
      steroidSkCdbo_position <- dbo_position[sapply(dbo_position, function(x) {
        if(all(x %in% steroidSkeleton_chain)) return(TRUE)
        else return(FALSE)
      })]
      steroidSkCdboNum <- length(steroidSkCdbo_position)
    }
  }
  names(steroidSkCC) <- "steroidSkCC"
  names(steroidSkCohNum) <- "steroidSkCohNum"
  names(steroidSkCohPos) <- paste0("steroidSkCohPos", 1:length(steroidSkCohPos))
  names(steroidSkCdbNum) <- "steroidSkCdbNum"
  names(steroidSkCdboNum) <- "steroidSkCdboNum"

  # 断链甾醇骨架及其变体
  secosteroidSkDC <- 0
  secosteroidSkDohNum <- 0
  secosteroidSkDohPos <- rep(0, minimumohNum)
  secosteroidSkDdbNum <- 0
  secosteroidSkDdboNum <- 0
  if(length(secosteroidSkeleton_derivative) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = secosteroidSkeleton_derivative, scriptPath = scriptPath)
    secosteroidSkDC <- length(symbols == "C")
    # 羟基与位置
    if(length(oh_position)){
      secosteroidSkDoh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% secosteroidSkeleton_derivative)) return(TRUE)
        else return(FALSE)
      })]
      secosteroidSkDohNum <- length(secosteroidSkDoh_position)
      if(secosteroidSkDohNum != 0){
        secosteroidSkDohPos <- sapply(secosteroidSkDoh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = secosteroidSkeleton_derivative[1], scriptPath = scriptPath))
        })
        secosteroidSkDohPos <- sort(secosteroidSkDohPos)
        secosteroidSkDohPos <- secosteroidSkDohPos[1:minimumohNum]
        secosteroidSkDohPos[is.na(secosteroidSkDohPos)] <- 0
      }
    }
    # 双键
    if(length(db_position) != 0){
      secosteroidSkDdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% secosteroidSkeleton_derivative)) return(TRUE)
        else return(FALSE)
      })]
      secosteroidSkDdbNum <- length(secosteroidSkDdb_position)
    }
    # 酮
    if(length(dbo_position) != 0){
      secsteroidSkDdbo_position <- dbo_position[sapply(dbo_position, function(x) {
        if(all(x %in% secosteroidSkeleton_derivative)) return(TRUE)
        else return(FALSE)
      })]
      secosteroidSkDdboNum <- length(secsteroidSkDdbo_position)
    }
  }
  names(secosteroidSkDC) <- "secosteroidSkDC"
  names(secosteroidSkDohNum) <- "secosteroidSkDohNum"
  names(secosteroidSkDohPos) <- paste0("secosteroidSkDohPos", 1:length(secosteroidSkDohPos))
  names(secosteroidSkDdbNum) <- "secsteroidSkDdbNum"
  names(secosteroidSkDdboNum) <- "secsteroidSkDdboNum"

  secosteroidSkCC <- 0
  secosteroidSkCohNum <- 0
  secosteroidSkCohPos <- rep(0, minimumohNum)
  secosteroidSkCdbNum <- 0
  secosteroidSkCdboNum <- 0
  # 断链甾醇骨架碳链
  if(length(secosteroidSkeleton_chain) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = secosteroidSkeleton_chain, scriptPath = scriptPath)
    secosteroidSkCC <- length(symbols == "C")
    # 羟基与位置
    if(length(oh_position)){
      secosteroidSkCoh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% secosteroidSkeleton_chain)) return(TRUE)
        else return(FALSE)
      })]
      secosteroidSkCohNum <- length(secosteroidSkCoh_position)
      if(secosteroidSkCohNum != 0){
        secosteroidSkCohPos <- sapply(secosteroidSkCoh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = secosteroidSkeleton_chain[1], scriptPath = scriptPath))
        })
        secosteroidSkCohPos <- sort(secosteroidSkCohPos)
        secosteroidSkCohPos <- secosteroidSkCohPos[1:minimumohNum]
        secosteroidSkCohPos[is.na(secosteroidSkCohPos)] <- 0
      }
    }
    # 双键
    if(length(db_position) != 0){
      secosteroidSkCdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% secosteroidSkeleton_chain)) return(TRUE)
        else return(FALSE)
      })]
      secosteroidSkCdbNum <- length(secosteroidSkCdb_position)
    }
    # 酮
    if(length(dbo_position) != 0){
      secsteroidSkCdbo_position <- dbo_position[sapply(dbo_position, function(x) {
        if(all(x %in% secosteroidSkeleton_chain)) return(TRUE)
        else return(FALSE)
      })]
      secosteroidSkCdboNum <- length(secsteroidSkCdbo_position)
    }
  }
  names(secosteroidSkCC) <- "secosteroidSkCC"
  names(secosteroidSkCohNum) <- "secosteroidSkCohNum"
  names(secosteroidSkCohPos) <- paste0("secosteroidSkCohPos", 1:length(secosteroidSkCohPos))
  names(secosteroidSkCdbNum) <- "secosteroidSkCdbNum"
  names(secosteroidSkCdboNum) <- "secosteroidSkCdboNum"

  mainChainC <- 0
  if(length(main_C_chains) != 0){
    mainChainC <- sum(sapply(main_C_chains, function(x) {
      length(which(.GetAtomSymbol(smi = smi, atom_idx_vector = x, scriptPath = scriptPath) == "C"))
    }))
  }
  names(mainChainC) <- "C"

  snInfo <- names(main_C_chains)[names(main_C_chains) %in% c("sn1", "sn2", "sn3")]
  sn1Num <- length(which(snInfo == "sn1" | snInfo == "sn3"))
  sn2Num <- length(which(snInfo == "sn2"))
  names(sn1Num) <- "sn1Num"
  names(sn2Num) <- "sn2Num"

  LipRtPredFP <- c(mainChainC, ohNum, ohPos, dbNum, dboNum, cooNum, conNum, cosNum, cocNum, ooNum, sn1Num, sn2Num,
                   sphC, sphohNum, sphohPos, sphdbNum, sphdboNum,
                   spiC, spiohNum, spiohPos, spidbNum, spidboNum,
                   steroidSkDC, steroidSkDohNum, steroidSkDohPos, steroidSkDdbNum, steroidSkDdboNum,
                   steroidSkCC, steroidSkCohNum, steroidSkCohPos, steroidSkCdbNum, steroidSkCdboNum,
                   secosteroidSkDC, secosteroidSkDohNum, secosteroidSkDohPos, secosteroidSkDdbNum, secosteroidSkDdboNum,
                   secosteroidSkCC, secosteroidSkCohNum, secosteroidSkCohPos, secosteroidSkCdbNum, secosteroidSkCdboNum,
                   glycerolNum, oh2dboNum, phosphateNum, betaineNum, carbohydrateNum, cholineNum, ethanolamineNum, serineNum, inositolNum, ethanolNum, threonineNum, carnitineNum, sulfonylNum)
  return(LipRtPredFP)
}


# test <- lapply(1:nrow(cmpDf_demo2), function(i) { #cmpDf_demo2 需要Check TODO
#   print(i)
#   x <- cmpDf_demo2$smiles[i]
#   .calLipRtPred_FP(smi = x, scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# })
