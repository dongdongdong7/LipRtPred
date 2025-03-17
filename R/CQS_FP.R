# Calculate CQS molecular fingerprints
# Barry Song
# 250317

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py

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
          return(x[-c(1,2)])
        }
      }
    })
    names(main_Chains) <- namesofmainChains
    main_Chains <- main_Chains[!sapply(main_Chains, is.null)]
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



# .calCQS_FP(smi = "[C@@H]1([C@H](O)[C@H](OP(=O)(O)O)[C@@H](COP(O)(=O)OP(O)(=O)OCC(C)([C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C/C=C\\CC(O)=O)C)O1)N1C=NC2C(N)=NC=NC1=2",
#            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.calCQS_FP <- function(smi, scriptPath){
  .lipidClassification(smi = smi, scriptPath = scriptPath)
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
  acyloxy_acylChain_position <- lapply(1:length(acyloxy_acylChain_position), function(i) {
    unique(c(acyloxy_position[[i]], acyloxy_acylChain_position[[i]]))
  })
  # 酰胺基链加上酰胺基键的位置
  amide_acylChain_position <- lapply(1:length(amide_acylChain_position), function(i) {
    unique(c(amide_position[[i]], amide_acylChain_position[[i]]))
  })
  # 硫代酰碳链加上硫代酰键的位置
  thioester_acylChain_position <- lapply(1:length(thioester_acylChain_position), function(i) {
    unique(c(thioester_position[[i]], thioester_acylChain_position[[i]]))
  })
  # 酰氧基链加上酰氧基键的位置
  acyloxy_alkoxyChain_position <- lapply(1:length(acyloxy_alkoxyChain_position), function(i) {
    unique(c(acyloxy_position[[i]], acyloxy_alkoxyChain_position[[i]]))
  })
  # 酰胺基链加上酰胺基键的位置
  amide_alkylaminoChain_position <- lapply(1:length(amide_alkylaminoChain_position), function(i) {
    unique(c(amide_position[[i]], amide_alkylaminoChain_position[[i]]))
  })
  # 硫代酰碳链加上硫代酰键的位置
  thioester_thioalkylChain_position <- lapply(1:length(thioester_thioalkylChain_position), function(i) {
    unique(c(thioester_position[[i]], thioester_thioalkylChain_position[[i]]))
  })
  main_Chains <- c(acyloxy_acylChain_position, amide_acylChain_position, thioester_acylChain_position,
                   acyloxy_alkoxyChain_position, amide_alkylaminoChain_position, thioester_thioalkylChain_position)
  namesofmainChains <- c(rep("acyloxy_acylChain", length(acyloxy_acylChain_position)),
                         rep("amide_acylChain", length(amide_acylChain_position)),
                         rep("thioester_acylChain", length(thioester_acylChain_position)),
                         rep("acyloxy_alkoxyChain", length(acyloxy_alkoxyChain_position)),
                         rep("amide_alkylaminoChain", length(amide_alkylaminoChain_position)),
                         rep("thioester_thioalkylChain", length(thioester_thioalkylChain_position)))
  for(head in c("acyloxy_acylChain", "amide_acylChain", "thioester_acylChain")){
    is <- which(namesofmainChains == head)
    for(i in is){
      x <- main_Chains[[i]]
      if(all(x == "none")) next
      for(j in (1:length(main_Chains))[-i]){
        y <- main_Chains[[j]]
        if(all(y[1:3] %in% x[-c(1:3)])) main_Chains[[j]] <- "none"
      }
    }
  }
  main_Chains <- lapply(main_Chains, function(x) {
    if(all(x == "none")) return(NULL)
    else return(x[-c(2,3)])
  })
  names(main_Chains) <- namesofmainChains
  main_Chains <- main_Chains[!sapply(main_Chains, is.null)]
  browser()
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
