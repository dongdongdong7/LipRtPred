# Search substructure of Sphingolipids (SP)
# Barry Song
# 250312

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py

# 1. Sphingosine
# (1) Search Sphingosine
# .searchSphingosine(smi = "OC[C@]([H])(N)[C@]([H])(O)/C=C/CCCCCCCCCCCCC",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSphingosine(smi = "[C@](CO)([H])(N)C(=O)/C=C/CCCCCCCCCCCCC",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSphingosine(smi = "[C@](C)([H])(NC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)[C@]([H])(O)/C=C/CCCCCCCCCCCCC",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSphingosine <- function(smi, scriptPath){
  # sphingosine
  position <- .GetSubstructMatches(smis = smi,
                                   SMARTS = "C(O)C(N)[C;$(C-O),$(C=O)](~[O;$(O-C),$(O=C)])[C;$(C-C),$(C=C),$(C#C)]~C",
                                   scriptPath = scriptPath)[[1]]
  names(position) <- rep("sphingosine", length(position))
  if(length(position) == 0){
    # spisulosine
    position <- .GetSubstructMatches(smis = smi,
                                     SMARTS = "[CH3]C(N)[C;$(C-O),$(C=O)](~[O;$(O-C),$(O=C)])[C;$(C-C),$(C=C),$(C#C)]~C",
                                     scriptPath = scriptPath)[[1]]
    names(position) <- rep("spisulosine", length(position))
  }
  return(position)
}

# (2) Search sphingosine C chain
# .searchSphingosine_Chain(smi = "[C@](CO)([H])(N)C(=O)/C=C/CCCCCCCCCCCCC",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSphingosine_Chain(smi = "[C@](C)([H])(NC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)[C@]([H])(O)/C=C/CCCCCCCCCCCCC",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSphingosine_Chain <- function(smi, scriptPath){
  sphingosine_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
  type <- names(sphingosine_position)[1]
  if(is.na(type)) return(list())
  if(type == "sphingosine"){
   position <- lapply(sphingosine_position, function(x) {
     c(.TraverseMolecule(smi = smi,
                         start_atom_idx = x[1],
                         non_traversable_atom_idx = x[c(2, 4, 6)],
                         scriptPath = scriptPath),
       x[2],x[6])
   })
   names(position) <- rep("sphingosine_chain", length(position))
  }else if(type == "spisulosine"){
    position <- lapply(sphingosine_position, function(x) {
      c(.TraverseMolecule(smi = smi,
                          start_atom_idx = x[1],
                          non_traversable_atom_idx = x[c(3, 5)],
                          scriptPath = scriptPath),
        x[5])
    })
    names(position) <- rep("spisulosine_chain", length(position))
  }else stop("Unkonwn sphingosine type")
  return(position)
}

# (3) Search sphingosine amide acylChain
# .searchSphAmide_AcylChain(smi = "[C@](C)([H])(NC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)[C@]([H])(O)/C=C/CCCCCCCCCCCCC",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSphAmide_AcylChain <- function(smi, scriptPath){
  sph_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
  if(length(sph_position) > 1) warnings("More than one sphingosine skeleton!")
  if(names(sph_position)[1] == "sphingosine") sph_N_position <- unlist(sph_position)[4]
  else if(names(sph_position)[1] == "spisulosine") sph_N_position <- unlist(sph_position)[3]
  amide_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
  if(length(amide_position) == 0) amide_acylChain_position <- list()
  else{
    l <- sapply(amide_position, function(x) {
      if(x[3] == sph_N_position) return(TRUE)
      else return(FALSE)
    })
    amide_acylChain_position <- .searchAmide_AcylChain(smi = smi, scriptPath = scriptPath)
    amide_acylChain_position <- amide_acylChain_position[l]
  }
  names(amide_acylChain_position) <- "SphAmide_AcylChain"
  return(amide_acylChain_position)
}

# (4) Search sphingosine ester
# .searchSphEster(smi = "C(OC(=O)CCCCCCCCCCCCCCC)[C@]([H])(NC(CCCCCCCCCCCCCCC)=O)[C@H](O)/C=C/CCCCCCCCCCCCC",
#                 scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSphEster <- function(smi, scriptPath){
  sph_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
  if(length(sph_position) > 1) warnings("More than one sphingosine skeleton!")
  if(names(sph_position)[1] == "sphingosine") sph_O_position <- unlist(sph_position)[c(2, 6)]
  else if(names(sph_position)[1] == "spisulosine") sph_O_position <- unlist(sph_position)[5]
  acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  if(length(acyloxy_position) == 0) sphEster_chain <- list()
  else{
    l <- sapply(acyloxy_position, function(x) {
      if(x[3] %in% sph_O_position) return(TRUE)
      else return(FALSE)
    })
    sphEster_position <- acyloxy_position[l]
  }
  names(sphEster_position) <- "sphEster"
  return(sphEster_position)
}

# (5) Search sphingosine ester chain
# .searchSphEster_Chain(smi = "C(OC(=O)CCCCCCCCCCCCCCC)[C@]([H])(NC(CCCCCCCCCCCCCCC)=O)[C@H](O)/C=C/CCCCCCCCCCCCC",
#                       scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSphEster_Chain <- function(smi, scriptPath){
  sph_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
  if(length(sph_position) > 1) warnings("More than one sphingosine skeleton!")
  if(names(sph_position)[1] == "sphingosine") sph_O_position <- unlist(sph_position)[c(2, 6)]
  else if(names(sph_position)[1] == "spisulosine") sph_O_position <- unlist(sph_position)[5]
  acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  if(length(acyloxy_position) == 0) sphEster_chain <- list()
  else{
    l <- sapply(acyloxy_position, function(x) {
      if(x[3] %in% sph_O_position) return(TRUE)
      else return(FALSE)
    })
    sphEster_chain <- .searchAcyloxy_AcylChain(smi = smi, scriptPath = scriptPath)
    sphEster_chain <- sphEster_chain[l]
  }
  names(sphEster_chain) <- "SphEsterChain"
  return(sphEster_chain)
}

# Search SP Main Chains
# .searchSP_MainChains(smi = "C(OC(=O)CCCCCCCCCCCCCCC)[C@]([H])(NC(CCCCCCCCCCCCCCC)=O)[C@H](O)/C=C/CCCCCCCCCCCCC",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSP_MainChains <- function(smi, scriptPath){
  sphingosine_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
  sphingosineChain_position <- .searchSphingosine_Chain(smi = smi, scriptPath = scriptPath)
  sphAmide_acylChain_position <- .searchSphAmide_AcylChain(smi = smi, scriptPath = scriptPath)
  sphEster_position <- .searchSphEster(smi = smi, scriptPath = scriptPath)
  sphEsterChain_position <- .searchSphEster_Chain(smi = smi, scriptPath = scriptPath)

  c(sphingosine_position, sphingosineChain_position, sphAmide_acylChain_position, sphEster_position, sphEsterChain_position)
}

# # 1. Fatty Chain Head
# # (1) Search sphingosine(homologs and variants)
# # .searchSphingosine(smi = "OC[C@]([H])(N)[C@]([H])(O)/C=C/CCCCCCCCCCCCC",
# #                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# # .searchSphingosine(smi = "[C@](CO)([H])(N)C(=O)/C=C/CCCCCCCCCCCCC",
# #                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSphingosine <- function(smi, scriptPath){
#   .GetSubstructMatches(smis = smi,
#                        SMARTS = "C(O)C(N)[C;$(C-O),$(C=O)](~[O;$(O-C),$(O=C)])[C;$(C-C),$(C=C)]~C",
#                        scriptPath = scriptPath)[[1]]
# }
#
# # (2) Search sphingosine C chain
# # .searchSphingosineC_Chain(smi = "[C@](CO)([H])(N)C(=O)/C=C/CCCCCCCCCCCCC",
# #                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSphingosineC_Chain <- function(smi, scriptPath){
#   sphingosine_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
#   lapply(sphingosine_position, function(x) {
#     .TraverseMolecule(smi = smi,
#                       start_atom_idx = x[3],
#                       non_traversable_atom_idx = x[c(1:2, 4)],
#                       scriptPath = scriptPath)
#   })
# }
# # (3) Search sphingosie N chain
# # .searchSphingosineN_Chain(smi = "P(=O)(OC[C@]([H])(NC(=O)CCCCCCCCCCC)[C@]([H])(O)/C=C/CCCCCCCCCCCCC)(O)O",
# #                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSphingosineN_Chain <- function(smi, scriptPath){
#   sphingosine_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
#   sphingosieAmide_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
#   sphingosieAmide_logical <- sapply(sphingosieAmide_position, function(x) {
#     if(x[3] %in% unlist(sphingosine_position)) return(TRUE)
#     else return(FALSE)
#   })
#   if(length(sphingosieAmide_logical) == 0) return(list())
#   sphingosieAmide_position <- sphingosieAmide_position[sphingosieAmide_logical]
#   lapply(sphingosieAmide_position, function(x) {
#     .TraverseMolecule(smi = smi,
#                       start_atom_idx = x[1],
#                       non_traversable_atom_idx = x[-1],
#                       scriptPath = scriptPath)
#   })
# }
#
# # (4) Search spisulosine
# # .searchSpisulosine(smi = "C(C)(N)C(O)CCCCCCCCCCC(C)CC(C)CC",
# #                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSpisulosine <- function(smi, scriptPath){
#   .GetSubstructMatches(smis = smi,
#                        SMARTS = "[CH3]C(N)[C;$(C-O),$(C=O)](~[O;$(O-C),$(O=C)])[C;$(C-C),$(C=C)]~C",
#                        scriptPath = scriptPath)[[1]]
# }
#
# # (5) Search spisulosine C Chain
# # .searchSpisulosineC_Chain(smi = "C(C)(N)C(O)CCCCCCCCCCC(C)CC(C)CC",
# #                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSpisulosineC_Chain <- function(smi, scriptPath){
#   spisulosine_position <- .searchSpisulosine(smi = smi, scriptPath = scriptPath)
#   lapply(spisulosine_position, function(x) {
#     .TraverseMolecule(smi = smi,
#                       start_atom_idx = x[2],
#                       non_traversable_atom_idx = x[c(1, 3)],
#                       scriptPath = scriptPath)
#   })
# }
# # (6) Search spisulosine N Chain
# # .searchSpisulosineN_Chain(smi = "[C@](C)([H])(NC(CCCCCCCCCCCCCCCCCCCCCCC)=O)[C@]([H])(O)CCCCCCCCCCCCCCC",
# #                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSpisulosineN_Chain <- function(smi, scriptPath){
#   spisulosine_position <- .searchSpisulosine(smi = smi, scriptPath = scriptPath)
#   spisulosineAmide_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
#   spisulosineAmide_logical <- sapply(spisulosineAmide_position, function(x) {
#     if(x[3] %in% unlist(spisulosine_position)) return(TRUE)
#     else return(FALSE)
#   })
#   if(length(spisulosineAmide_logical) == 0) return(list())
#   spisulosineAmide_position <- spisulosineAmide_position[spisulosineAmide_logical]
#   lapply(spisulosineAmide_position, function(x) {
#     .TraverseMolecule(smi = smi,
#                       start_atom_idx = x[1],
#                       non_traversable_atom_idx = x[-1],
#                       scriptPath = scriptPath)
#   })
# }
