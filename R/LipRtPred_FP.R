# Calculate LipRtPred Fingerprints
# 250319
# Barry Song

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py
# minimumohNum: minimum number of OH to be taken into account when calculating the OH distance
# 脂肪酰
# .calLipRtPred_FP(smi = "N(C(=O)[C@H](CC/C=C\\CCCCC(CCCCC(CCCCCCCCCCO[C@@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(O)=O)=O)=O)OC)CCC1=CC(O)=C(O)C=C1",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# 甘油
# .calLipRtPred_FP(smi = "OC[C@]1(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H](CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@@]([H])(CO)COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC1)C)[H]",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# 甘油磷酸
# .calLipRtPred_FP(smi = "O(CC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C[C@]([H])(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)COP(O)(=O)OC[C@]([H])(O)COP(OC)(O)=O",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# # 鞘脂
# .calLipRtPred_FP(smi = "C(OC(=O)CCCCCCCCCCCCCCC)[C@]([H])(NC(CCCCCCCCCCCCCCC)=O)[C@H](O)/C=C/CCCCCCCCCCCCC",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# # 固醇
# .calLipRtPred_FP(smi = "C1[C@H](OC(=O)CCCCCCCCCCC)CC2=CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@]([H])(C)CCCC(C)C)[C@@]4(C)CC[C@]3([H])[C@@]2(C)C1",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .calLipRtPred_FP(smi = "C1C(=C)/C(=C\\C=C2\\[C@]3([H])CC[C@@]([H])([C@@]3(C)CCC\\2)[C@]([H])(C)CCCC(C)C)/C[C@@H](O)C1",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.calLipRtPred_FP <- function(smi, scriptPath, minimumohNum = 6){

  C_Chains_position <- list()
  sphingosineChain_position <- list()
  steroidSkD_position <- list()
  steroidSkC_position <- list()
  secoSkD_position <- list()
  secoSkC_position <- list()

  category <- .lipidClassification(smi = smi, scriptPath = scriptPath)
  if(category == "Sterol"){
    ST_MainChains <- .searchST_MainChains(smi = smi, scriptPath = scriptPath)
    C_Chains_position <- ST_MainChains[names(ST_MainChains) == "sterylEsterChain"]
    steroidSkD_position <- ST_MainChains[names(ST_MainChains) == "steroidSkD"]
    steroidSkC_position <- ST_MainChains[names(ST_MainChains) == "steroidSkC"]
    secoSkD_position <- ST_MainChains[names(ST_MainChains) == "secoSkD"]
    secoSkC_position <- ST_MainChains[names(ST_MainChains) == "secoSkC"]
  }else if(category == "Sphingolipid"){
    SP_MainChains <- .searchSP_MainChains(smi = smi, scriptPath = scriptPath)
    C_Chains_position <- SP_MainChains[names(SP_MainChains) %in% c("SphEsterChain", "SphAmide_AcylChain")]
    sphingosineChain_position <- SP_MainChains[names(SP_MainChains) == "sphingosine_chain"]
  }else if(category == "Glycerophospholipid" | category == "Glycerolipid"){
    GL_MainChains <- .searchGL_MainChains(smi = smi, scriptPath = scriptPath)
    C_Chains_position <- GL_MainChains[names(GL_MainChains) %in% c("sn1", "sn2", "sn3")]
  }else if(category == "Fatty Acyls"){
    FA_MainChains <- .searchFA_MainChains(smi = smi, scriptPath = scriptPath)
    C_Chains_position <- FA_MainChains[names(FA_MainChains) %in% c("acyloxy_acylChain", "amide_acylChain", "thioester_acylChain", "acyloxy_alkoxyChain", "amide_alkylaminoChain", "thioester_thioalkylChain",
                                                                   "fattyAlcoholChain", "fattyAldehydeChain", "fattyNitrileChain", "fattyEtherChain", "hydrocarbons")]
  }

  # Special Group
  # Glycerol
  glcerol_position <- .searchGlycerol(smi = smi, scriptPath = scriptPath)
  glycerolNum <- length(glcerol_position)
  names(glycerolNum) <- "glycerolNum"
  # Phosphate
  phosphate_position <- .searchPhosphate(smi = smi, scriptPath = scriptPath)
  phosphateNum <- length(phosphate_position)
  names(phosphateNum) <- "phosphateNum"
  # Betaine
  betaine_position <- .searchBetaine(smi = smi, scriptPath = scriptPath)
  betaineNum <- length(betaine_position)
  names(betaineNum) <- "betaineNum"
  # Pentose
  pentose_position <- .searchPentose(smi = smi, scriptPath = scriptPath)
  pentoseNum <- length(pentose_position)
  names(pentoseNum) <- "pentoseNum"
  # Hexose
  hexose_position <- .searchHexose(smi = smi, scriptPath = scriptPath)
  hexoseNum <- length(hexose_position)
  names(hexoseNum) <- "hexoseNum"
  # Choline
  choline_position <- .searchCholine(smi = smi, scriptPath = scriptPath)
  cholineNum <- length(choline_position)
  names(cholineNum) <- "cholineNum"
  # Ethanolamine
  ethanolamine_position <- .searchEthanolamine(smi = smi, scriptPath = scriptPath)
  ethanolamineNum <- length(ethanolamine_position)
  names(ethanolamineNum) <- "ethanolamineNum"
  # Serine
  serine_position <- .searchSerine(smi = smi, scriptPath = scriptPath)
  serineNum <- length(serine_position)
  names(serineNum) <- "serineNum"
  # Inositol
  inositol_position <- .searchInositol(smi = smi, scriptPath = scriptPath)
  inositolNum <- length(inositol_position)
  names(inositolNum) <- "inositolNum"
  # Ethanol
  ethanol_position <- .searchEthanol(smi = smi, scriptPath = scriptPath)
  ethanolNum <- length(ethanol_position)
  names(ethanolNum) <- "ethanolNum"
  # Threonine
  threonine_position <- .searchThreonine(smi = smi, scriptPath = scriptPath)
  threonineNum <- length(threonine_position)
  names(threonineNum) <- "threonineNum"
  # Carnitine
  carnitine_position <- .searchCarnitine(smi = smi, scriptPath = scriptPath)
  carnitineNum <- length(carnitine_position)
  names(carnitineNum) <- "carnitineNum"
  if(carnitineNum > 0 & cholineNum > 0) cholineNum <- cholineNum - carnitineNum # Delete choline in carnitine
  # Sulfonyl
  sulfonyl_position <- .searchSulfonyl(smi = smi, scriptPath = scriptPath)
  sulfonylNum <- length(sulfonyl_position)
  names(sulfonylNum) <- "sulfonylNum"
  # Glucuronic acid
  gluacid_position <- .searchGluAcid(smi = smi, scriptPath = scriptPath)
  gluacidNum <- length(gluacid_position)
  names(gluacidNum) <- "gluacidNum"

  # Remove special group
  remove_position <- unique(unlist(c(glcerol_position, phosphate_position, betaine_position, pentose_position, hexose_position,
                                     choline_position, ethanolamine_position, serine_position, inositol_position, ethanol_position,
                                     threonine_position, carnitine_position, sulfonyl_position, gluacid_position)))
  if(length(remove_position) != 0) remove_position <- remove_position[.GetAtomSymbol(smi = smi, atom_idx_vector = remove_position, scriptPath = scriptPath) == "C"]
  C_Chains_position <- lapply(C_Chains_position, function(x) {
    x[!x %in% remove_position]
  })

  # Common group
  # hydroxyl group
  ohNum <- 0
  ohPos <- rep(0, minimumohNum)
  oh_position <- .GetSubstructMatches(smis = smi,
                                      SMARTS = "[C,c;!$(C=O);$([C,c]-[OH])]-[OH]",
                                      scriptPath = scriptPath)[[1]]
  if(length(oh_position) != 0){
    ohNum <- length(which(sapply(oh_position, function(x) {
      if(all(x %in% unlist(C_Chains_position))) return(TRUE)
      else return(FALSE)
    })))
    if(length(C_Chains_position) != 0 & ohNum != 0){
      ohPos <- sapply(oh_position, function(x) {
        idx <- which(sapply(C_Chains_position, function(y) {
          if(all(x %in% y)) return(TRUE)
          else return(FALSE)
        }))
        if(length(idx) == 0) return(NA)
        nmc <- C_Chains_position[[idx]]
        length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = nmc[1], scriptPath = scriptPath))
      })
      ohPos <- sort(ohPos)
      ohPos <- ohPos[1:minimumohNum]
      ohPos[is.na(ohPos)] <- 0
    }
  }
  names(ohNum) <- "ohNum"
  names(ohPos) <- paste0("ohPos", 1:length(ohPos))
  # Double bond
  dbNum <- 0
  db_position <- .GetSubstructMatches(smis = smi,
                                      SMARTS = "C=C",
                                      scriptPath = scriptPath)[[1]]
  if(length(db_position) != 0){
    dbNum <- length(which(sapply(db_position, function(x) {
      if(all(x %in% unlist(C_Chains_position))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(dbNum) <- "dbNum"
  # Triple bond
  tpNum <- 0
  tp_position <- .GetSubstructMatches(smis = smi,
                                      SMARTS = "C#C",
                                      scriptPath = scriptPath)[[1]]
  if(length(tp_position) != 0){
    tpNum <- length(which(sapply(tp_position, function(x) {
      if(all(x %in% unlist(C_Chains_position))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(tpNum) <- "tpNum"
  # Ketone
  ktNum <- 0
  kt_position <- .GetSubstructMatches(smis = smi,
                                          SMARTS = "[C;$(C=O);!$(C-[N,O,S,P])]=[OX1]",
                                          scriptPath = scriptPath)[[1]]
  if(length(kt_position) != 0){
    ktNum <- length(which(sapply(kt_position, function(x) {
      if(all(x %in% unlist(C_Chains_position))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(ktNum) <- "ktNum"
  # Acyloxy
  cooNum <- 0
  coo_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  if(length(coo_position) != 0){
    cooNum <- length(which(sapply(coo_position, function(x) {
      if(all(x %in% unlist(C_Chains_position))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(cooNum) <- "cooNum"
  # Amide
  conNum <- 0
  con_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
  if(length(con_position) != 0){
    conNum <- length(which(sapply(con_position, function(x) {
      if(all(x %in% unlist(C_Chains_position))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(conNum) <- "conNum"
  # Thioester
  cosNum <- 0
  cos_position <- .searchThioester(smi = smi, scriptPath = scriptPath)
  if(length(cos_position) != 0){
    cosNum <- length(which(sapply(cos_position, function(x) {
      if(all(x %in% unlist(C_Chains_position))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(cosNum) <- "cosNum"
  # Ether
  cocNum <- 0
  coc_position <- .GetSubstructMatches(smis = smi,
                                       SMARTS = "[C;!$(C=O);$(C-O)]-[OX2;OH0]-[C;!$(C=O);$(C-O)]",
                                       scriptPath = scriptPath)[[1]]
  if(length(coc_position) != 0){
    cocNum <- length(which(sapply(coc_position, function(x) {
      if(all(x %in% unlist(C_Chains_position))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(cocNum) <- "cocNum"
  # Peroxide bond
  ooNum <- 0
  oo_position <- .GetSubstructMatches(smis = smi,
                                      SMARTS = "[O]-[O]",
                                      scriptPath = scriptPath)[[1]]
  if(length(oo_position) != 0){
    ooNum <- length(which(sapply(oo_position, function(x) {
      if(all(x %in% unlist(C_Chains_position))) return(TRUE)
      else return(FALSE)
    })))
  }
  names(ooNum) <- "ooNum"

  C <- 0
  if(length(C_Chains_position) != 0){
    C <- sum(sapply(C_Chains_position, function(x) {
      length(which(.GetAtomSymbol(smi = smi, atom_idx_vector = x, scriptPath = scriptPath) == "C"))
    }))
  }
  names(C) <- "C"

  snInfo <- names(C_Chains_position)[names(C_Chains_position) %in% c("sn1", "sn2", "sn3")]
  sn1Num <- length(which(snInfo == "sn1" | snInfo == "sn3"))
  sn2Num <- length(which(snInfo == "sn2"))
  names(sn1Num) <- "sn1Num"
  names(sn2Num) <- "sn2Num"

  # Skeleton
  # Sphingosine
  sphC <- 0
  sphohNum <- 0
  sphohPos <- rep(0,minimumohNum)
  sphdbNum <- 0
  sphktNum <- 0
  if(length(sphingosineChain_position) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = unlist(sphingosineChain_position), scriptPath = scriptPath)
    sphC <- length(which(symbols == "C"))
    if(length(oh_position) != 0){
      sphoh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% unlist(sphingosineChain_position))) return(TRUE)
        else return(FALSE)
      })]
      sphohNum <- length(sphoh_position)
      if(sphohNum != 0){
        sphohPos <- sapply(sphoh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = unlist(sphingosineChain_position)[1], scriptPath = scriptPath))
        })
        sphohPos <- sort(sphohPos)
        sphohPos <- sphohPos[1:minimumohNum]
        sphohPos[is.na(sphohPos)] <- 0
      }
    }
    if(length(db_position) != 0){
      sphdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% unlist(sphingosineChain_position))) return(TRUE)
        else return(FALSE)
      })]
      sphdbNum <- length(sphdb_position)
    }
    if(length(kt_position) != 0){
      sphkt_position <- kt_position[sapply(kt_position, function(x) {
        if(all(x %in% unlist(sphingosineChain_position))) return(TRUE)
        else return(FALSE)
      })]
      sphktNum <- length(sphkt_position)
    }
  }
  names(sphC) <- "sphC"
  names(sphohNum) <- "sphohNum"
  names(sphohPos) <- paste0("sphohPos", 1:length(sphohPos))
  names(sphdbNum) <- "sphdbNum"
  names(sphktNum) <- "sphktNum"

  # Steroid Skeleton Derivative
  stdSkDC <- 0
  stdSkDohNum <- 0
  stdSkDohPos <- rep(0,minimumohNum)
  stdSkDdbNum <- 0
  stdSkDktNum <- 0
  if(length(steroidSkD_position) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = unlist(steroidSkD_position), scriptPath = scriptPath)
    stdSkDC <- length(which(symbols == "C"))
    if(length(oh_position) != 0){
      stdSkDoh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% unlist(steroidSkD_position))) return(TRUE)
        else return(FALSE)
      })]
      stdSkDohNum <- length(stdSkDoh_position)
      if(stdSkDohNum != 0){
        stdSkDohPos <- sapply(stdSkDoh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = unlist(steroidSkD_position)[1], scriptPath = scriptPath))
        })
        stdSkDohPos <- sort(stdSkDohPos)
        stdSkDohPos <- stdSkDohPos[1:minimumohNum]
        stdSkDohPos[is.na(stdSkDohPos)] <- 0
      }
    }
    if(length(db_position) != 0){
      stdSkDdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% unlist(steroidSkD_position))) return(TRUE)
        else return(FALSE)
      })]
      stdSkDdbNum <- length(stdSkDdb_position)
    }
    if(length(kt_position) != 0){
      stdSkDkt_position <- kt_position[sapply(kt_position, function(x) {
        if(all(x %in% unlist(steroidSkD_position))) return(TRUE)
        else return(FALSE)
      })]
      stdSkDktNum <- length(stdSkDkt_position)
    }
  }
  names(stdSkDC) <- "stdSkDC"
  names(stdSkDohNum) <- "stdSkDohNum"
  names(stdSkDohPos) <- paste0("stdSkDohPos", 1:length(stdSkDohPos))
  names(stdSkDdbNum) <- "stdSkDdbNum"
  names(stdSkDktNum) <- "stdSkDktNum"

  # Steroid Skeleton Chain
  stdSkCC <- 0
  stdSkCohNum <- 0
  stdSkCohPos <- rep(0,minimumohNum)
  stdSkCdbNum <- 0
  stdSkCktNum <- 0
  if(length(steroidSkC_position) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = unlist(steroidSkC_position), scriptPath = scriptPath)
    stdSkCC <- length(which(symbols == "C"))
    if(length(oh_position) != 0){
      stdSkCoh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% unlist(steroidSkC_position))) return(TRUE)
        else return(FALSE)
      })]
      stdSkCohNum <- length(stdSkCoh_position)
      if(stdSkCohNum != 0){
        stdSkCohPos <- sapply(stdSkCoh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = unlist(steroidSkC_position)[1], scriptPath = scriptPath))
        })
        stdSkCohPos <- sort(stdSkCohPos)
        stdSkCohPos <- stdSkCohPos[1:minimumohNum]
        stdSkCohPos[is.na(stdSkCohPos)] <- 0
      }
    }
    if(length(db_position) != 0){
      stdSkCdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% unlist(steroidSkC_position))) return(TRUE)
        else return(FALSE)
      })]
      stdSkCdbNum <- length(stdSkCdb_position)
    }
    if(length(kt_position) != 0){
      stdSkCkt_position <- kt_position[sapply(kt_position, function(x) {
        if(all(x %in% unlist(steroidSkC_position))) return(TRUE)
        else return(FALSE)
      })]
      stdSkCktNum <- length(stdSkCkt_position)
    }
  }
  names(stdSkCC) <- "stdSkCC"
  names(stdSkCohNum) <- "stdSkCohNum"
  names(stdSkCohPos) <- paste0("stdSkCohPos", 1:length(stdSkCohPos))
  names(stdSkCdbNum) <- "stdSkCdbNum"
  names(stdSkCktNum) <- "stdSkCktNum"

  # Secosteroid Skeleton Derivative
  secSkDC <- 0
  secSkDohNum <- 0
  secSkDohPos <- rep(0, minimumohNum)
  secSkDdbNum <- 0
  secSkDktNum <- 0
  if(length(secoSkD_position) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = unlist(secoSkD_position), scriptPath = scriptPath)
    secSkDC <- length(symbols == "C")
    # 羟基与位置
    if(length(oh_position)){
      secSkDoh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% unlist(secoSkD_position))) return(TRUE)
        else return(FALSE)
      })]
      secSkDohNum <- length(secSkDoh_position)
      if(secSkDohNum != 0){
        secSkDohPos <- sapply(secSkDoh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = unlist(secoSkD_position)[1], scriptPath = scriptPath))
        })
        secSkDohPos <- sort(secSkDohPos)
        secSkDohPos <- secSkDohPos[1:minimumohNum]
        secSkDohPos[is.na(secSkDohPos)] <- 0
      }
    }
    # 双键
    if(length(db_position) != 0){
      secSkDdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% unlist(secoSkD_position))) return(TRUE)
        else return(FALSE)
      })]
      secSkDdbNum <- length(secSkDdb_position)
    }
    # 酮
    if(length(kt_position) != 0){
      secSkDkt_position <- kt_position[sapply(kt_position, function(x) {
        if(all(x %in% unlist(secoSkD_position))) return(TRUE)
        else return(FALSE)
      })]
      secSkDktNum <- length(secSkDkt_position)
    }
  }
  names(secSkDC) <- "secSkDC"
  names(secSkDohNum) <- "secSkDohNum"
  names(secSkDohPos) <- paste0("secSkDohPos", 1:length(secSkDohPos))
  names(secSkDdbNum) <- "secSkDdbNum"
  names(secSkDktNum) <- "secSkDktNum"

  # Secosteroid Skeleton Chain
  secSkCC <- 0
  secSkCohNum <- 0
  secSkCohPos <- rep(0, minimumohNum)
  secSkCdbNum <- 0
  secSkCktNum <- 0
  if(length(secoSkC_position) != 0){
    symbols <- .GetAtomSymbol(smi = smi, atom_idx_vector = unlist(secoSkC_position), scriptPath = scriptPath)
    secSkCC <- length(symbols == "C")
    # 羟基与位置
    if(length(oh_position)){
      secSkCoh_position <- oh_position[sapply(oh_position, function(x) {
        if(all(x %in% unlist(secoSkC_position))) return(TRUE)
        else return(FALSE)
      })]
      secSkCohNum <- length(secSkCoh_position)
      if(secSkCohNum != 0){
        secSkCohPos <- sapply(secSkCoh_position, function(x) {
          length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = unlist(secoSkC_position)[1], scriptPath = scriptPath))
        })
        secSkCohPos <- sort(secSkCohPos)
        secSkCohPos <- secSkCohPos[1:minimumohNum]
        secSkCohPos[is.na(secSkCohPos)] <- 0
      }
    }
    # 双键
    if(length(db_position) != 0){
      secSkCdb_position <- db_position[sapply(db_position, function(x) {
        if(all(x %in% unlist(secoSkC_position))) return(TRUE)
        else return(FALSE)
      })]
      secSkCdbNum <- length(secSkCdb_position)
    }
    # 酮
    if(length(kt_position) != 0){
      secSkCkt_position <- kt_position[sapply(kt_position, function(x) {
        if(all(x %in% unlist(secoSkC_position))) return(TRUE)
        else return(FALSE)
      })]
      secSkCktNum <- length(secSkCkt_position)
    }
  }
  names(secSkCC) <- "secSkCC"
  names(secSkCohNum) <- "secSkCohNum"
  names(secSkCohPos) <- paste0("secSkCohPos", 1:length(secSkCohPos))
  names(secSkCdbNum) <- "secSkCdbNum"
  names(secSkCktNum) <- "secSkCktNum"

  LipRtPredFP <- c(C, ohNum, ohPos, dbNum, tpNum, ktNum, cooNum, conNum, cosNum, cocNum, ooNum, sn1Num, sn2Num,
                   sphC, sphohNum, sphohPos, sphdbNum, sphktNum,
                   stdSkDC, stdSkDohNum, stdSkDohPos, stdSkDdbNum, stdSkDktNum,
                   stdSkCC, stdSkCohNum, stdSkCohPos, stdSkCdbNum, stdSkCktNum,
                   secSkDC, secSkDohNum, secSkDohPos, secSkDdbNum, secSkDktNum,
                   secSkCC, secSkCohNum, secSkCohPos, secSkCdbNum, secSkCktNum,
                   glycerolNum, phosphateNum, betaineNum, pentoseNum, hexoseNum, cholineNum, ethanolamineNum, serineNum, inositolNum, ethanolNum, threonineNum, carnitineNum, sulfonylNum, gluacidNum)

  return(LipRtPredFP)
}

# lg <- as.logical(sapply(cmpDf_test$smiles, function(x) {
#   .CheckSMILES(smi = x, scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# }))
# cmpDf_test <- cmpDf_test[lg, ]
# test <- lapply(1:nrow(cmpDf_test), function(i) {
#   print(i)
#   x <- cmpDf_test$smiles[i]
#   .calLipRtPred_FP(smi = x, scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# })
# i <- 149
# draw_smi(smi = cmpDf_test$smiles[i], width = 1600, height = 800)
# test[[i]]
