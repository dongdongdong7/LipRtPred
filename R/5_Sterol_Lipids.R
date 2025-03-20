# Search sterol lipids
# Barry Song
# 250313

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py

# 1. Fatty Chain Head
# (1) Search Steroid skeleton
# .searchSteroidSkeleton(smi = "C1C[C@H](O)C(C)(C)[C@]2([H])CCC3[C@]4(C)CC[C@]([H])([C@H](C)CCC(=C)C(C)C)[C@@]4(C)CCC=3[C@@]12C",
#                        scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSteroidSkeleton(smi = "C1C2[C@@]3([H])CC[C@]4(C)C(=O)CC[C@@]4([H])[C@]3([H])CCC=2C=C(O)C=1",
#                        scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSteroidSkeleton <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[#6]1~[#6]~[#6]~[#6]2~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4[#6]~3~[#6]~[#6]~[#6]~2[#6]~1",
                       scriptPath = scriptPath)[[1]]
}

# (2) Search Steroid skeleton derivative
# .searchSteroidSkeleton_Derivative(smi = "C1C[C@H](O)C(C)(C)[C@]2([H])CCC3[C@]4(C)CC[C@]([H])([C@H](C)CCC(=C)C(C)C)[C@@]4(C)CCC=3[C@@]12C",
#                                   scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSteroidSkeleton_Derivative <- function(smi, scriptPath){
  steroidSkeleton_position <- .searchSteroidSkeleton(smi = smi, scriptPath = scriptPath)
  lapply(steroidSkeleton_position, function(x) {
    c(.TraverseMolecule(smi = smi,
                        start_atom_idx = x[1],
                        non_traversable_atom_idx = x[c(9, 10)],
                        scriptPath = scriptPath),
      x[c(9, 10)])
  })
}

# (3) Search steryl ester
# .searchSterylEster_Chain(smi = "[C@]12(CC=C3C[C@@H](OC(CCCCC/C=C\\CCCCCCCCC)=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@@]([H])([C@@](C)([H])CCCC(C)C)CC[C@@]21[H])[H]",
#                          scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSterylEster_Chain(smi = "[C@]12(CC=C3C[C@@H](O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@@]([H])([C@@](C)([H])CCCC(C)C)CC[C@@]21[H])[H]",
#                          scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSterylEster_Chain <- function(smi, scriptPath){
  posi <- .GetSubstructMatches(smis = smi,
                               SMARTS = "[#6]1(-[OX2]-[CX3;$(C=O);$(C-O)](=O))~[#6]~[#6]~[#6]2~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4[#6]~3~[#6]~[#6]~[#6]~2[#6]~1",
                               scriptPath = scriptPath)[[1]]
  posi <- lapply(posi, function(x) {x[3]})
  acyloxy_acylChain_position <- .searchAcyloxy_AcylChain(smi = smi, scriptPath = scriptPath)
  acyloxy_acylChain_position <- lapply(acyloxy_acylChain_position, function(x) {
    if(x[1] %in% unlist(posi)) return(x)
    else return(NULL)
  })
  acyloxy_acylChain_position <- acyloxy_acylChain_position[!sapply(acyloxy_acylChain_position, is.null)]
  return(acyloxy_acylChain_position)
}

# (4) Search steroid skeleton Chain
# .searchSteroidSkeleton_Chain(smi = "C1[C@]2(C)[C@@]3([H])CC[C@]4(C)[C@@]([H])([C@]([H])(C)CC/C(=C/CC)/C(C)C)CC[C@@]4([H])[C@]3([H])CC=C2C[C@@H](O)C1",
#                              scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSteroidSkeleton_Chain <- function(smi, scriptPath){
  steroidSkeleton_position <- .searchSteroidSkeleton(smi = smi, scriptPath = scriptPath)
  lapply(steroidSkeleton_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[9],
                      non_traversable_atom_idx = x[-9],
                      scriptPath = scriptPath)
  })
}

# (5) Search secosteroid skeleton
# .searchSecosteroidSkeleton(smi = "C1=C(C)C(CC[C@@]2(O)[C@]3([H])CC[C@@]([H])([C@@]3(C)CC[C@H]2O)[C@]([H])(C)CC/C(=C/C)/C(C)C)=CC(O)=C1",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSecosteroidSkeleton(smi = "C1[C@@H](O)C[C@]2([H])CC[C@@]3([H])[C@]4([H])CCC(=O)[C@@]4(C)CC[C@]3([H])[C@@]2(C)C1",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSecosteroidSkeleton <- function(smi, scriptPath){
  steroidSkeletion_position <- .searchSteroidSkeleton(smi = smi, scriptPath = scriptPath)
  secosteroidSkeletion_position <- .GetSubstructMatches(smis = smi,
                                                        SMARTS = "[#6](O)(~[#6]~[#6]~[#6]1)~[#6]~[#6]~1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]~[#6]~2~3", # 15
                                                        scriptPath = scriptPath)[[1]]
  if(length(secosteroidSkeletion_position) == 0) return(secosteroidSkeletion_position)
  secosteroidSkeletion_position <- secosteroidSkeletion_position[sapply(secosteroidSkeletion_position, function(x) {
    if(all(x[-2] %in% unlist(steroidSkeletion_position))) return(FALSE)
    else return(TRUE)
  })]
  return(secosteroidSkeletion_position)
}

# (6) Search secosteroid skeleton derivative
# .searchSecosteroidSkeleton_Derivative(smi = "C1C(=C)/C(=C/C=C2\\[C@]3([H])CC[C@@]([H])([C@@]3(C)CCC\\2)[C@@](C)([H])CCCC(O)(C)C)/C[C@H](O)C1",
#                                       scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSecosteroidSkeleton_Derivative <- function(smi, scriptPath){
  secosteroidSkeleton_position <- .searchSecosteroidSkeleton(smi = smi, scriptPath = scriptPath)
  lapply(secosteroidSkeleton_position, function(x) {
    c(.TraverseMolecule(smi = smi,
                        start_atom_idx = x[1],
                        non_traversable_atom_idx = x[15],
                        scriptPath = scriptPath),
      x[15])
  })
}

# (7) Search secosteroid skeleton chain
# .searchSecosteroidSkeleton_Chain(smi = "C1=C(C)C(CC[C@@]2(O)[C@]3([H])CC[C@@]([H])([C@@]3(C)CC[C@H]2O)[C@]([H])(C)CC/C(=C/C)/C(C)C)=CC(O)=C1",
#                                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSecosteroidSkeleton_Chain <- function(smi, scriptPath){
  secosteroidSkeleton_position <- .searchSecosteroidSkeleton(smi = smi, scriptPath = scriptPath)
  lapply(secosteroidSkeleton_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[15],
                      non_traversable_atom_idx = x[-15],
                      scriptPath = scriptPath)
  })
}
