# Calculate CQS molecular fingerprints
# Barry Song
# 250317

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py

# Search FA main chains
# .searchFA_MainChains(smi = "[C@@H]1([C@H](O)[C@H](OP(=O)(O)O)[C@@H](COP(O)(=O)OP(O)(=O)OCC(C)([C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C/C=C\\CC(O)=O)C)O1)N1C=NC2C(N)=NC=NC1=2",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchFA_MainChains(smi = "OC(CCCCCCC(=O)O)=O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchFA_MainChains(smi = "C(/C#CC#CC(O)/C=C/CCCCCCC)=C\\CO",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchFA_MainChains(smi = "C(CCC/C=C\\C=C/CCCC)(=O)[H]",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchFA_MainChains <- function(smi, scriptPath){
  main_class <- .lipidClassification(smi = smi, scriptPath = scriptPath)
  if(main_class == "Fatty Acyls"){
    acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
    amide_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
    thioester_position <- .searchThioester(smi = smi, scriptPath = scriptPath)
    acyloxy_acylChain_position <- .searchAcyloxy_AcylChain(smi = smi, scriptPath = scriptPath)
    names(acyloxy_acylChain_position) <- rep("acyloxy_acylChain", length(acyloxy_acylChain_position))
    amide_acylChain_position <- .searchAmide_AcylChain(smi = smi, scriptPath = scriptPath)
    names(amide_acylChain_position) <- rep("amide_acylChain", length(amide_acylChain_position))
    thioester_acylChain_position <- .searchThioester_AcylChain(smi = smi, scriptPath = scriptPath)
    names(thioester_acylChain_position) <- rep("thioester_acylChain", length(thioester_acylChain_position))
    acyloxy_alkoxyChain_position <- .searchAcyloxy_AlkoxyChain(smi = smi, scriptPath = scriptPath)
    names(acyloxy_alkoxyChain_position) <- rep("acyloxy_alkoxyChain", length(acyloxy_alkoxyChain_position))
    amide_alkylaminoChain_position <- .searchAmide_AlkylaminoChain(smi = smi, scriptPath = scriptPath)
    names(amide_alkylaminoChain_position) <- rep("amide_alkylaminoChain", length(amide_alkylaminoChain_position))
    thioester_thioalkylChain_position <- .searchThioester_ThioalkylChain(smi = smi, scriptPath = scriptPath)
    names(thioester_thioalkylChain_position) <- rep("thioester_thioalkylChain", length(thioester_thioalkylChain_position))
    # 酰氧基链加上酰氧基键的位置
    if(length(acyloxy_acylChain_position) != 0){
      acyloxy_acylChain_position <- lapply(1:length(acyloxy_acylChain_position), function(i) {
        unique(c(acyloxy_position[[i]], acyloxy_acylChain_position[[i]]))
      })
    }
    # 酰胺基链加上酰胺基键的位置
    if(length(amide_acylChain_position) != 0){
      amide_acylChain_position <- lapply(1:length(amide_acylChain_position), function(i) {
        unique(c(amide_position[[i]], amide_acylChain_position[[i]]))
      })
    }
    # 硫代酰碳链加上硫代酰键的位置
    if(length(thioester_acylChain_position) != 0){
      thioester_acylChain_position <- lapply(1:length(thioester_acylChain_position), function(i) {
        unique(c(thioester_position[[i]], thioester_acylChain_position[[i]]))
      })
    }
    # 酰氧基链加上酰氧基键的位置
    if(length(acyloxy_alkoxyChain_position) != 0){
      acyloxy_alkoxyChain_position <- lapply(1:length(acyloxy_alkoxyChain_position), function(i) {
        unique(c(acyloxy_position[[i]], acyloxy_alkoxyChain_position[[i]]))
      })
    }
    # 酰胺基链加上酰胺基键的位置
    if(length(amide_alkylaminoChain_position) != 0){
      amide_alkylaminoChain_position <- lapply(1:length(amide_alkylaminoChain_position), function(i) {
        unique(c(amide_position[[i]], amide_alkylaminoChain_position[[i]]))
      })
    }
    # 硫代酰碳链加上硫代酰键的位置
    if(length(thioester_thioalkylChain_position) != 0){
      thioester_thioalkylChain_position <- lapply(1:length(thioester_thioalkylChain_position), function(i) {
        unique(c(thioester_position[[i]], thioester_thioalkylChain_position[[i]]))
      })
    }
    main_Chains <- c(acyloxy_acylChain_position, amide_acylChain_position, thioester_acylChain_position,
                     acyloxy_alkoxyChain_position, amide_alkylaminoChain_position, thioester_thioalkylChain_position)
    namesofmainChains <- c(rep("acyloxy_acylChain", length(acyloxy_acylChain_position)),
                           rep("amide_acylChain", length(amide_acylChain_position)),
                           rep("thioester_acylChain", length(thioester_acylChain_position)),
                           rep("acyloxy_alkoxyChain", length(acyloxy_alkoxyChain_position)),
                           rep("amide_alkylaminoChain", length(amide_alkylaminoChain_position)),
                           rep("thioester_thioalkylChain", length(thioester_thioalkylChain_position)))
    for(head in c("acyloxy_acylChain", "amide_acylChain", "thioester_acylChain", "acyloxy_alkoxyChain", "amide_alkylaminoChain", "thioester_thioalkylChain")){
      is <- which(namesofmainChains == head)
      if(length(is) == 0) next
      for(i in is){
        x <- main_Chains[[i]]
        if(all(x == "none")) next
        for(j in (1:length(main_Chains))[-i]){
          y <- main_Chains[[j]]
          if(all(y == "none")) next
          if(all(y[1:3] %in% x[-c(1:3)])) main_Chains[[j]] <- "none"
        }
      }
    }
    main_Chains <- lapply(1:length(main_Chains), function(i) {
      x <- main_Chains[[i]]
      if(all(x == "none")) return(NULL)
      else{
        if(namesofmainChains[i] %in% c("acyloxy_acylChain", "amide_acylChain", "thioester_acylChain")){
          return(x[-c(2,3)])
        }else{
          return(x[-c(1,2,3)])
        }
      }
    })
    names(main_Chains) <- namesofmainChains
    main_Chains <- main_Chains[sapply(main_Chains, function(x) {length(x) != 0})]
    #main_Chains <- main_Chains[!sapply(main_Chains, is.null)]
    return(main_Chains)
  }
  else if(main_class == "Fatty Alcohol"){
    main_Chains <- .searchFattyAlcohols_Chain(smi = smi, scriptPath = scriptPath)
    names(main_Chains) <- "Fatty Alcohol"
  }else if(main_class == "Fatty Aldehyde"){
    main_Chains <- .searchFattyAldehydes_Chain(smi = smi, scriptPath = scriptPath)
    names(main_Chains) <- "Fatty Aldehyde"
  }else if(main_class == "Fatty Nitrile"){
    main_Chains <- .searchFattyNitriles_Chain(smi = smi, scriptPath = scriptPath)
    names(main_Chains) <- "Fatty Nitrile"
  }else if(main_class == "Fatty Ether"){
    main_Chains <- .searchFattyEthers_chain(smi = smi, scriptPath = scriptPath)
    names(main_Chains) <- "Fatty Ether"
  }else if(main_class == "Hydrocarbons"){
    main_Chains <- .searchHydrocarbons(smi = smi, scriptPath = scriptPath)
    names(main_Chains) <- "Hydrocarbons"
  }
  else{
    main_Chains <- list()
  }
  return(main_Chains)
}

# Search GL main chains
# .searchGL_MainChains(smi = "OC[C@]1(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H](CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@@]([H])(CO)COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC1)C)[H]",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGL_MainChains(smi = "OC[C@@](OCCCCCCCCC1CC2C3C4CCC4C3C2CC1)([H])COCCCCCCC1CC2C3C=CCCC3C2C=C1",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGL_MainChains(smi = "C(O[C@@H]1O[C@@H]([C@@H]([C@@H]([C@H]1O)O)O)CO)[C@]([H])(O)COC(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGL_MainChains(smi = "O[C@@H]1[C@@H](O)[C@@H](O)[C@@H](CO[C@@H]2[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O2)O[C@H]1OC[C@]([H])(O)COC(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGL_MainChains(smi = "[C@](CO[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CS(O)(=O)=O)O1)([H])(O)COCCCCCCCCCCCCCCCC",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGL_MainChains <- function(smi, scriptPath){
  main_class <- .lipidClassification(smi = smi, scriptPath = scriptPath)
  if(main_class != "Glycerolipid"){
    return(list())
  }
  glycerol_position <- .searchGlycerol(smi = smi, scriptPath = scriptPath)
  dha_position <- .searchDihydroxyacetone(smi = smi, scriptPath = scriptPath)
  glycerolEster_position <- .searchGlycerolEster(smi = smi, scriptPath = scriptPath)
  glycerolEther_position <- .searchGlycerolEther(smi = smi, scriptPath = scriptPath)
  dhaEster_position <- .searchDihydroxyacetoneEster(smi = smi, scriptPath = scriptPath)
  dhaEther_position <- .searchDihydroxyacetoneEther(smi = smi, scriptPath = scriptPath)
  glycerolEsterChain_position <- .searchGlycerolEster_Chain(smi = smi, scriptPath = scriptPath)
  glycerolEtherChain_position <- .searchGlycerolEther_Chain(smi = smi, scriptPath = scriptPath)
  dhaEsterChain_position <- .searchDihydroxyacetoneEster_Chain(smi = smi, scriptPath = scriptPath)
  dhaEtherChain_position <- .searchDihydroxyacetoneEther_Chain(smi = smi, scriptPath = scriptPath)
  pentose_position <- .searchPentose(smi = smi, scriptPath = scriptPath)
  hexose_position <- .searchHexose(smi = smi, scriptPath = scriptPath)
  glycerolEtherChain_position <- lapply(glycerolEtherChain_position, function(x) {
    if(x[1] %in% unlist(c(pentose_position, hexose_position))) return(NULL)
    else return(x)
  })
  glycerolEtherChain_position <- glycerolEtherChain_position[!sapply(glycerolEtherChain_position, is.null)]
  dhaEtherChain_position <- lapply(dhaEtherChain_position, function(x) {
    if(x[1] %in% unlist(c(pentose_position, hexose_position))) return(NULL)
    else return(x)
  })
  dhaEtherChain_position <- dhaEtherChain_position[!sapply(dhaEtherChain_position, is.null)]
  main_chains <- c(glycerolEsterChain_position, glycerolEtherChain_position, dhaEsterChain_position, dhaEtherChain_position)
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
  return(main_chains)
}

# Search GP main chains
# .searchGP_MainChains(smi = "O(CC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C[C@]([H])(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)COP(O)(=O)OC[C@]([H])(O)COP(OC)(O)=O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGP_MainChains(smi = "O(P(OC[C@@]1([H])OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@]([H])(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H](CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC1)C)CO)(O)=O)CCN",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGP_MainChains(smi = "C([C@@]1([H])OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H](CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@]([H])(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC1)COP(O)(=O)OCCN)C)(C(O)(CO)C(O)C(O)C(O)CO)O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGP_MainChains(smi = "[C@](COP(=O)(O)OCCNC[C@@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O)([H])(OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)CO/C=C\\CCCCCCCCCCCCCC",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGP_MainChains <- function(smi, scriptPath){
  main_class <- .lipidClassification(smi = smi, scriptPath = scriptPath)
  if(main_class != "Glycerophospholipid") return(list())
  glycerolPhosphate_position <- .searchGlycerolPhosphate(smi = smi, scriptPath = scriptPath)
  glycerolEsterChain_position <- .searchGlycerolEster_Chain(smi = smi, scriptPath = scriptPath)
  glycerolEtherChain_position <- .searchGlycerolEther_Chain(smi = smi, scriptPath = scriptPath)
  dhaEsterChain_position <- .searchDihydroxyacetoneEster_Chain(smi = smi, scriptPath = scriptPath)
  dhaEtherChain_position <- .searchDihydroxyacetoneEther_Chain(smi = smi, scriptPath = scriptPath)
  main_chains <- c(glycerolEsterChain_position, glycerolEtherChain_position, dhaEsterChain_position, dhaEtherChain_position)
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
  return(main_chains)
}

# Search SP main chains
# .searchSP_MainChains(smi = "C(OC(=O)CCCCCCCCCCCCCCC)[C@]([H])(NC(CCCCCCCCCCCCCCC)=O)[C@H](O)/C=C/CCCCCCCCCCCCC",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSP_MainChains(smi = "[C@](C)([H])(NC(CCCCCCCCCCCCCCCCCCCCCCC)=O)[C@]([H])(O)CCCCCCCCCCCCCCC",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSP_MainChains <- function(smi, scriptPath){
  main_class <- .lipidClassification(smi = smi, scriptPath = scriptPath)
  if(main_class != "Sphingolipid") return(list())
  sphingosine_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
  spisulosine_position <- .searchSpisulosine(smi = smi, scriptPath = scriptPath)
  sphingosineC_chain_position <- .searchSphingosineC_Chain(smi = smi, scriptPath = scriptPath)
  sphingosineN_chain_position <- .searchSphingosineN_Chain(smi = smi, scriptPath = scriptPath)
  spisulosineC_chain_position <- .searchSpisulosineC_Chain(smi = smi, scriptPath = scriptPath)
  spisulosineN_chain_position <- .searchSpisulosineN_Chain(smi = smi, scriptPath = scriptPath)
  acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  acyloxy_acylChain_position <- .searchAcyloxy_AcylChain(smi = smi, scriptPath = scriptPath)
  if(length(acyloxy_acylChain_position) != 0){
    acyloxy_acylChain_position <- lapply(1:length(acyloxy_acylChain_position), function(i) {
      unique(c(acyloxy_position[[i]], acyloxy_acylChain_position[[i]]))
    })
  }
  SP_other_chain_position <- lapply(acyloxy_acylChain_position, function(x){
    if(x[3] %in% unlist(c(sphingosine_position, spisulosine_position))){
      return(x[-c(2,3)])
    }else{
      return(NULL)
    }
  })
  names(sphingosineC_chain_position) <- rep("sphingosineC_chain", length(sphingosineC_chain_position))
  names(sphingosineN_chain_position) <- rep("sphingosineN_chain", length(sphingosineN_chain_position))
  names(spisulosineC_chain_position) <- rep("spisulosineC_chain", length(spisulosineC_chain_position))
  names(spisulosineN_chain_position) <- rep("spisulosineN_chain", length(spisulosineN_chain_position))
  names(SP_other_chain_position) <- rep("SP_other_chain", length(SP_other_chain_position))
  main_chains <- c(sphingosineC_chain_position, sphingosineN_chain_position, spisulosineC_chain_position, spisulosineN_chain_position, SP_other_chain_position)
  return(main_chains)
}

# Search ST main chains
# .searchST_MainChains(smi = "C1[C@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)CC2=C(O)C(=O)C3=C(CC[C@]4(C3=CC[C@]4([H])[C@](O)(C)[C@H](O)CCC(C)C)C)[C@]21C",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchST_MainChains(smi = "[C@]12(CC=C3C[C@@H](OC(CCCCC/C=C\\CCCCCCCCC)=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@@]([H])([C@@](C)([H])CCCC(C)C)CC[C@@]21[H])[H]",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchST_MainChains(smi = "C1C(=C)/C(=C\\C=C2\\[C@]3([H])CC[C@@]([H])([C@@]3(C)CCC\\2)[C@]([H])(C)/C=C/[C@H](C)C(C)C)/C[C@@H](O)C1",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchST_MainChains(smi = "C1[C@]2(C)[C@@]3([H])CC[C@]4(C)[C@@]([H])([C@]([H])(C)CCC(O)=O)CC[C@@]4([H])[C@]3([H])C(=O)C[C@@]2([H])CCC1",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchST_MainChains <- function(smi, scriptPath){
  main_class <- .lipidClassification(smi = smi, scriptPath = scriptPath)
  if(main_class != "Sterol") return(list())
  steroidSkeleton_position <- .searchSteroidSkeleton(smi = smi, scriptPath = scriptPath)
  steroidSkeleton_chain <- .searchSteroidSkeleton_Chain(smi = smi, scriptPath = scriptPath)
  steroidSkeleton_derivative <- .searchSteroidSkeleton_Derivative(smi = smi, scriptPath = scriptPath)
  secosteroidSkeleton_chain <- .searchSecosteroidSkeleton_Chain(smi = smi, scriptPath = scriptPath)
  secosteroidSkeleton_derivative <- .searchSecosteroidSkeleton_Derivative(smi = smi, scriptPath = scriptPath)
  sterylEsterChain_position <- .searchSterylEster_Chain(smi = smi, scriptPath = scriptPath)
  acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  if(length(sterylEsterChain_position) != 0){
    sterylEsterChain_position_new <- lapply(sterylEsterChain_position, function(x) {
      l <- sapply(acyloxy_position, function(y) {
        y[1] %in% x
      })
      acy <- acyloxy_position[[which(l)]]
      unique(c(acy, x))
    })
    steroidSkeleton_derivative <- lapply(steroidSkeleton_derivative, function(x) {
      x[!x %in% unlist(sterylEsterChain_position_new)]
    })
  }

  names(steroidSkeleton_chain) <- rep("steroidSkeleton_chain", length(steroidSkeleton_chain))
  names(sterylEsterChain_position) <- rep("sterylEster_chain", length(sterylEsterChain_position))
  names(steroidSkeleton_derivative) <- rep("steroidSkeleton_derivative", length(steroidSkeleton_derivative))
  names(secosteroidSkeleton_chain) <- rep("secosteroidSkeleton_chain", length(secosteroidSkeleton_chain))
  names(secosteroidSkeleton_derivative) <- rep("secosteroidSkeleton_derivative", length(secosteroidSkeleton_derivative))

  main_chains <- c(steroidSkeleton_chain, sterylEsterChain_position, steroidSkeleton_derivative,
                   secosteroidSkeleton_chain, secosteroidSkeleton_derivative)
  return(main_chains)
}

# .calCQS_FP_old <- function(smi, scriptPath){
#   browser()
#   acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
#   amide_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
#   thioester_position <- .searchThioester(smi = smi, scriptPath = scriptPath)
#   acyloxy_aclyChain_position <- .searchAcyloxy_AcylChain(smi = smi, scriptPath = scriptPath)
#   amide_acylChain_position <- .searchAmide_AcylChain(smi = smi, scriptPath = scriptPath)
#   thioester_acylChain_position <- .searchThioester_AcylChain(smi = smi, scriptPath = scriptPath)
#   acyloxy_alkoxyChain_position <- .searchAcyloxy_AlkoxyChain(smi = smi, scriptPath = scriptPath)
#   amide_alkylaminoChain_position <- .searchAmide_AlkylaminoChain(smi = smi, scriptPath = scriptPath)
#   thioester_thioalkylChain_position <- .searchThioester_ThioalkylChain(smi = smi, scriptPath = scriptPath)
#
  # # 酰氧基链加上酰氧基键的位置
  # acyloxy_aclyChain_position <- lapply(1:length(acyloxy_aclyChain_position), function(i) {
  #   unique(c(acyloxy_position[[i]], acyloxy_aclyChain_position[[i]]))
  # })
  # # 酰胺基链加上酰胺基键的位置
  # amide_acylChain_position <- lapply(1:length(amide_acylChain_position), function(i) {
  #   unique(c(amide_position[[i]], amide_acylChain_position[[i]]))
  # })
  # # 硫代酰碳链加上硫代酰键的位置
  # thioester_acylChain_position <- lapply(1:length(thioester_acylChain_position), function(i) {
  #   unique(c(thioester_position[[i]], thioester_acylChain_position[[i]]))
  # })
#   # 删去与酰氧碳链逆位的酰胺碳链
#   l <- lapply(acyloxy_aclyChain_position, function(x) {
#     sapply(amide_position, function(y){
#       if(all(y %in% x)){
#         z <- x[x %in% y]
#         z <- z[z != y[2]]
#         if(all(z == y[-2])) return(TRUE)
#         else return(FALSE)
#       }else return(FALSE)
#     })
#   })
#   amide_acylChain_position <- amide_acylChain_position[!sapply(1:length(thioester_position), function(i) {
#     any(sapply(l, function(x) {l[[i]]}))
#   })]
#   # 删去与酰氧碳链逆位的硫代碳链
#   l <- lapply(acyloxy_aclyChain_position, function(x) {
#     sapply(thioester_position, function(y){
#       if(all(y %in% x)){
#         z <- x[x %in% y]
#         z <- z[z != y[2]]
#         if(all(z == y[-2])) return(TRUE)
#         else return(FALSE)
#       }else return(FALSE)
#     })
#   })
#   thioester_acylChain_position <- thioester_acylChain_position[!sapply(1:length(thioester_position), function(i) {
#     any(sapply(l, function(x) {l[[i]]}))
#   })]
#
#   names(acyloxy_aclyChain_position) <- rep("acyloxy", length(acyloxy_aclyChain_position))
#   names(amide_acylChain_position) <- rep("amide", length(amide_acylChain_position))
#   names(thioester_acylChain_position) <- rep("thioester", length(thioester_acylChain_position))
#   acylChains <- c(acyloxy_aclyChain_position, amide_acylChain_position, thioester_acylChain_position)
#   namesofacylChains <- names(acylChains)
#   acylChains <- lapply(1:length(acylChains), function(i) {
#     x <- acylChains[[i]]
#     Del_logical <- sapply((1:length(acylChains))[-i], function(j) {
#       y <- acylChains[[j]]
#       if(all(y %in% x)) return(TRUE)
#       else return(FALSE)
#     })
#     setdiff(x, unique(unlist(acylChains[-i][Del_logical])))[-c(2, 3)] # 去除 OH 和 =O
#   })
#   names(acylChains) <- namesofacylChains
#
#   acyloxy_alkoxyChain_position <- lapply(acyloxy_alkoxyChain_position, function(x) {
#     if(length(x) == 1) return(NULL)
#     if(any(x[-1] %in% unlist(acylChains))) return(NULL)
#     else return(x)
#   })
#   acyloxy_alkoxyChain_position <- acyloxy_alkoxyChain_position[!sapply(acyloxy_alkoxyChain_position, is.null)]
#   amide_alkylaminoChain_position <- lapply(amide_alkylaminoChain_position, function(x) {
#     if(length(x) == 1) return(NULL)
#     if(any(x[-1] %in% unlist(acylChains))) return(NULL)
#     else return(x)
#   })
#   amide_alkylaminoChain_position <- amide_alkylaminoChain_position[!sapply(amide_alkylaminoChain_position, is.null)]
#   thioester_thioalkylChain_position <- lapply(thioester_thioalkylChain_position, function(x) {
#     if(length(x) == 1) return(NULL)
#     if(any(x[-1] %in% unlist(acylChains))) return(NULL)
#     else return(x)
#   })
#   thioester_thioalkylChain_position <- thioester_thioalkylChain_position[!sapply(thioester_thioalkylChain_position, is.null)]
# }
