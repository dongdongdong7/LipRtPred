# Calculate CQS Fingerprints
# Barry Song
# 250327

# system.file("python", "molecule_operation.py", package = "LipRtPred")
# TODO EC_raw4ML 看到了 1253 磷酸后加甲基的怎么办 1515 胆固醇多了个甲基 1536 一个不认识的结构
.calCQS_FP <- function(smi, minimumohNum = 1, scriptPath){

  category <- .lipidClassification(smi = smi, scriptPath = scriptPath)
  if(category == "Sterol"){
    ST_MainChains <- .searchST_MainChains(smi = smi, scriptPath = scriptPath)
    C_Chains_position <- ST_MainChains[names(ST_MainChains) == "sterylEsterChain"]
  }else if(category == "Sphingolipid"){
    SP_MainChains <- .searchSP_MainChains(smi = smi, scriptPath = scriptPath)
    C_Chains_position <- SP_MainChains[names(SP_MainChains) %in% c("SphEsterChain", "sphingosine_chain", "SphAmide_AcylChain")]
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
  # DEMA
  dema_position <- .searchDEMA(smi = smi, scriptPath = scriptPath)
  demaNum <- length(dema_position)
  names(demaNum) <- "demaNum"
  # MMEA
  mmea_position <- .searchMMEA(smi = smi, scriptPath = scriptPath)
  mmeaNum <- length(mmea_position)
  names(mmeaNum) <- "mmeaNum"
  # Sialic acid
  sialicAcid_position <- .searchSialicAcid(smi = smi, scriptPath = scriptPath)
  sialicAcidNum <- length(sialicAcid_position)
  names(sialicAcidNum) <- "sialicAcidNum"
  # sulfatedgalactosyl
  sulfatedgalactosyl_position <- .searchSulfatedGalactosyl(smi = smi, scriptPath = scriptPath)
  sulfatedgalactosylNum <- length(sulfatedgalactosyl_position)
  names(sulfatedgalactosylNum) <- "sulfatedgalactosylNum"
  if(sulfatedgalactosylNum > 0 & sulfonylNum > 0 & hexoseNum > 0){
    sulfonylNum <- sulfonylNum - sulfatedgalactosylNum
    hexoseNum <- hexoseNum - sulfatedgalactosylNum
  }
  # Trimethylhomoserine
  trimethylhomoserine_position <- .searchTrimethylhomoserine(smi = smi, scriptPath = scriptPath)
  trimethylhomoserineNum <- length(trimethylhomoserine_position)
  names(trimethylhomoserineNum) <- "trimethylhomoserineNum"

  # Remove special group
  remove_position <- unique(unlist(c(glcerol_position, phosphate_position, betaine_position, pentose_position, hexose_position,
                                     choline_position, ethanolamine_position, serine_position, inositol_position, ethanol_position,
                                     threonine_position, carnitine_position, sulfonyl_position, gluacid_position, dema_position,
                                     mmea_position, sialicAcid_position, sulfatedgalactosyl_position, trimethylhomoserine_position)))
  #if(length(remove_position) != 0) remove_position <- remove_position[.GetAtomSymbol(smi = smi, atom_idx_vector = remove_position, scriptPath = scriptPath) == "C"]
  C_Chains_position <- lapply(C_Chains_position, function(x) {
    x[!x %in% remove_position]
  })
  if(length(C_Chains_position) > 0){
    C_Chains_position <- C_Chains_position[sapply(C_Chains_position, function(x) {
      length(which(.GetAtomSymbol(smi = smi, atom_idx_vector = x, scriptPath = scriptPath) == "C"))
    }) > 0]
  }

  phytosphingosine_position <- .searchPhytosphingosine(smi = smi, scriptPath = scriptPath)
  phytosphingosineNum <- length(phytosphingosine_position)
  names(phytosphingosineNum) <- "phytosphingosineNum"

  dihydrosphingosine_position <- .searchDihydrosphingosine(smi = smi, scriptPath = scriptPath)
  dihydrosphingosineNum <- length(dihydrosphingosine_position)
  names(dihydrosphingosineNum) <- "dihydrosphingosineNum"

  sphingosine_position <- .searchSphingosine0(smi = smi, scriptPath = scriptPath)
  sphingosineNum <- length(sphingosine_position)
  names(sphingosineNum) <- "sphingosineNum"

  sphPosi <- unique(unlist(c(phytosphingosine_position, dihydrosphingosine_position, sphingosine_position)))
  if(length(sphPosi) != 0){
    sphPosi <- sphPosi[.GetAtomSymbol(smi = smi, atom_idx_vector = sphPosi, scriptPath = scriptPath) != "C"]
    C_Chains_position <- lapply(C_Chains_position, function(x) {
      x[!x %in% sphPosi]
    })
    if(length(C_Chains_position) > 0){
      C_Chains_position <- C_Chains_position[sapply(C_Chains_position, function(x) {
        length(which(.GetAtomSymbol(smi = smi, atom_idx_vector = x, scriptPath = scriptPath) == "C"))
      }) > 0]
    }
  }

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
        length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = nmc[1], scriptPath = scriptPath)) - 1
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
  # Acyloxy
  cooNum <- 0
  coo_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  if(length(coo_position) != 0){
    coo_position <- coo_position[sapply(coo_position, function(x) {
      if(all(x %in% unlist(C_Chains_position))) return(TRUE)
      else return(FALSE)
    })]
    cooNum <- length(coo_position)
  }
  names(cooNum) <- "cooNum"
  # COOH
  coohNum <- 0
  cooh_position <- .GetSubstructMatches(smis = smi,
                                        SMARTS = "[CX3;$(C=O);$(C-O)](=O)[OX2;OH]",
                                        scriptPath = scriptPath)[[1]]
  if(length(cooh_position) != 0){
    cooh_position <- cooh_position[sapply(cooh_position, function(x) {
      if(all(x %in% unlist(coo_position))) return(TRUE)
      else return(FALSE)
    })]
    coohNum <- length(cooh_position)
  }
  names(coohNum) <- "coohNum"
  if(coohNum > 0 & cooNum > 0) cooNum <- cooNum - coohNum
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
  # allyl ether
  allylEtherNum <- 0
  allylEther_position <- .GetSubstructMatches(smis = smi,
                                              SMARTS = "C=C-O-C",
                                              scriptPath = scriptPath)[[1]]
  if(length(allylEther_position) != 0){
    allylEther_position <- allylEther_position[sapply(allylEther_position, function(x) {
      idx <- which(sapply(C_Chains_position, function(y) {
        if(all(x %in% y)) return(TRUE)
        else return(FALSE)
      }))
      if(length(idx) == 0) return(NA)
      nmc <- C_Chains_position[[idx]]
      len1 <- length(.GetShortestPath(smi = smi, start_atom_idx = x[2], end_atom_idx = nmc[1], scriptPath = scriptPath))
      len2 <- length(.GetShortestPath(smi = smi, start_atom_idx = x[3], end_atom_idx = nmc[1], scriptPath = scriptPath))
      if(len1 > len2) return(TRUE)
      else return(FALSE)
    })]
  }
  allylEtherNum <- length(allylEther_position)
  names(allylEtherNum) <- "allylEtherNum"
  if(allylEtherNum > 0 & cocNum > 0) cocNum <- cocNum - allylEtherNum
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
  sn1Num <- length(which(snInfo == "sn1"))
  sn2Num <- length(which(snInfo == "sn2"))
  sn3Num <- length(which(snInfo == "sn3"))
  names(sn1Num) <- "sn1Num"
  names(sn2Num) <- "sn2Num"
  names(sn3Num) <- "sn3Num"

  cholesterol_position <- .searchCholesterol(smi = smi, scriptPath = scriptPath)
  cholesterolNum <- length(cholesterol_position)
  names(cholesterolNum) <- "cholesterolNum"

  ergosterol_position <- .searchErgosterol(smi = smi, scriptPath = scriptPath)
  ergosterolNum <- length(ergosterol_position)
  names(ergosterolNum) <- "ergosterolNum"

  brassicasterol_position <- .searchBrassicasterol(smi = smi, scriptPath = scriptPath)
  brassicasterolNum <- length(brassicasterol_position)
  names(brassicasterolNum) <- "brassicasterolNum"

  ethylcholesterol_position <- .searchEthylcholesterol(smi = smi, scriptPath = scriptPath)
  alpha_ethylcholesterolNum <- 0
  beta_ethylcholesterolNum <- 0
  ethylcholesterolNum <- length(ethylcholesterol_position)
  if(ethylcholesterolNum > 0){
    if(names(ethylcholesterol_position) == "alpha"){
      alpha_ethylcholesterolNum <- ethylcholesterolNum
    }else if(names(ethylcholesterol_position) == "beta"){
      beta_ethylcholesterolNum <- ethylcholesterolNum
    }
  }
  names(alpha_ethylcholesterolNum) <- "alpha_ethylcholesterolNum"
  names(beta_ethylcholesterolNum) <- "beta_ethylcholesterolNum"

  stigmasterol_position <- .searchStigmasterol(smi = smi, scriptPath = scriptPath)
  stigmasterolNum <- length(stigmasterol_position)
  names(stigmasterolNum) <- "stigmasterolNum"

  cholestanol_position <- .searchCholestanol(smi = smi, scriptPath = scriptPath)
  cholestanolNum <- length(cholestanol_position)
  names(cholestanolNum) <- "cholestanolNum"

  CQS_FP <- c(C, dbNum, ohNum, ohPos, coohNum, cooNum, sn1Num, sn2Num, sn3Num, carnitineNum,
              allylEtherNum, cocNum, cholesterolNum, cholineNum, demaNum, ethanolNum,
              ethanolamineNum, gluacidNum, glycerolNum, hexoseNum, inositolNum, mmeaNum,
              sialicAcidNum, ooNum, phosphateNum, phytosphingosineNum, ergosterolNum, brassicasterolNum,
              alpha_ethylcholesterolNum, stigmasterolNum, serineNum, dihydrosphingosineNum, sphingosineNum,
              cholestanolNum, beta_ethylcholesterolNum, sulfatedgalactosylNum, sulfonylNum, trimethylhomoserineNum)
  names(CQS_FP) <- c("c", "d", "h", "n", "o", "a", "alpha", "beta", "gama", "m", "f", "j",
                     "k", "q", "r", "u", "e", "g", "v", "w", "i", "x", "y", "z", "p", "delta",
                     "epsilon", "theta", "xi", "pi", "s", "eta", "miu", "omega", "lambda",
                     "psi", "phi", "t")
  return(CQS_FP)
}

# Search Phytosphingosine
# .searchPhytosphingosine(smi = "CCCCCCCCCCCCCC[C@H]([C@H]([C@H](CO)N)O)O",
#                         scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhytosphingosine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CH2](O)[CX4;CH](N)[CX4;CH](O)[CX4;CH](O)", # [CH2](O)[CX4;CH](N)[CX4;CH](O)[CX4;CH](O)CCCCCCCCCCCCC[CH3]
                       scriptPath = scriptPath)[[1]]
}

# Search Dihydrosphingosine
# .searchDihydrosphingosine(smi = "CCCCCCCCCCCCCCC[C@H]([C@H](CO)N)O",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchDihydrosphingosine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CH2](O)[CX4;CH](N)[CX4;CH](O)[CX4;CH2]",
                       scriptPath = scriptPath)[[1]]
}

# Search Sphingosine
# .searchSphingosine0(smi = "CCCCCCCCCCCCC/C=C/[C@H]([C@H](CO)N)O",
#                     scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSphingosine0 <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CH2](O)[CX4;CH](N)[CX4;CH](O)C=C",
                       scriptPath = scriptPath)[[1]]
}

# Search Cholesterol
# .searchCholesterol(smi = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchCholesterol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)CCCC(C)C)CCC4C3CC=C2C1",
                       scriptPath = scriptPath)[[1]]
}

# Search Ergosterol
# .searchErgosterol(smi = "C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
#                   scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchErgosterol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)C=CC(C)C(C)C)CCC4C3=CC=C2C1",
                       scriptPath = scriptPath)[[1]]
}

# Search Brassicasterol
# .searchBrassicasterol(smi = "C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
#                       scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchBrassicasterol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)C=CC(C)C(C)C)CCC4C3CC=C2C1",
                       scriptPath = scriptPath)[[1]]
}

# Search Ethylcholesterol
# alpha
# .searchEthylcholesterol(smi = "CC[C@@H](C(C)C)CC[C@@H](C)[C@@]1([H])CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@@]21C",
#                         scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# beta
# .searchEthylcholesterol(smi = "CC[C@H](C(C)C)CC[C@@H](C)[C@@]1([H])CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@@]21C",
#                         scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchEthylcholesterol <- function(smi, scriptPath){
  posi <- .GetSubstructMatches(smis = smi,
                               SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)CCC(CC)C(C)C)CCC4C3CC=C2C1",
                               scriptPath = scriptPath)[[1]]
  if(length(posi) != 0){
    ch <- .GetAtomCip(smi = smi, atom_idx_vector = unlist(posi)[c(1, 17)], scriptPath = scriptPath)
    if(ch[1] == ch[2]) names(posi) <- "beta"
    else names(posi) <- "alpha"
  }
  return(posi)
}

# Search Stigmasterol
# .searchStigmasterol(smi = "CC[C@H](/C=C/[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C",
#                     scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchStigmasterol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)C=CC(CC)C(C)C)CCC4C3CC=C2C1",
                       scriptPath = scriptPath)[[1]]
}

# Search Cholestanol
# .searchCholestanol(smi = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CC[C@@H](C4)O)C)C",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchCholestanol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)CCCC(C)C)CCC4C3CCC2C1",
                       scriptPath = scriptPath)[[1]]
}
# .getCQSFP(cmpDf = cmpDf_demo2)
.getCQSFP <- function(cmpDf, minimumohNum = 1, thread = 1){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  wrong_idx <- .check_smiles(smiles = cmpDf$smiles)
  if(length(wrong_idx) != 0){
    message(paste0("Wrong Idx: ", paste0(wrong_idx, collapse = " ")))
    cmpDf <- cmpDf[-wrong_idx, ]
  }
  message("Calculate CQS FP...")
  smis <- cmpDf$smiles
  scriptPath <- system.file("python", "molecule_operation.py", package = "LipRtPred")
  pb <- utils::txtProgressBar(max = length(smis), style = 3)
  progress <- function(n){utils::setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  cl <- snow::makeCluster(thread)
  doSNOW::registerDoSNOW(cl)
  #parallel::clusterExport(cl, c(".calCQS_FP"), envir = environment())
  fpsList <- foreach::`%dopar%`(foreach::foreach(smi = smis, n = 1:length(smis),
                                                 .packages = c(),
                                                 .export = c(".calCQS_FP"),
                                                 .options.snow = opts),
                                {
                                  CQS <- .calCQS_FP(smi = smi, minimumohNum = minimumohNum, scriptPath = scriptPath)
                                  CQSDF <- data.frame(matrix(CQS, nrow = 1))
                                  colnames(CQSDF) <- names(CQS)
                                  return(CQSDF)
                                })
  snow::stopCluster(cl)
  gc()
  fpsDf <- purrr::list_rbind(fpsList)
  return(dplyr::as_tibble(cbind(cmpDf, fpsDf)))
}

#' @rdname getMD
#' @param minimumohNum minimum number of OH to be taken into account when calculating the OH distance
#' @export
#'
#' @examples
#' fpsDf <- GetCQS_FP(cmpDf = cmpDf_demo2, thread = 1)
GetCQS_FP <- function(cmpDf, flavor = "CxSmiles", minimumohNum = 1, thread = 1){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  fpsDf <- .getCQSFP(cmpDf = cmpDf, minimumohNum = minimumohNum, thread = thread)
  return(fpsDf)
}
