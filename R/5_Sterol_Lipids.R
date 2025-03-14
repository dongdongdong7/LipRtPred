# Search sterol lipids
# Barry Song
# 250313

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py

# 1. Fatty Chain Head
# (1) Search sterol skeleton
# .searchSterolSkeleton(smi = "C1C[C@H](O)C(C)(C)[C@]2([H])CCC3[C@]4(C)CC[C@]([H])([C@H](C)CCC(=C)C(C)C)[C@@]4(C)CCC=3[C@@]12C",
#                       scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSterolSkeleton(smi = "C1C2[C@@]3([H])CC[C@]4(C)C(=O)CC[C@@]4([H])[C@]3([H])CCC=2C=C(O)C=1",
#                       scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSterolSkeleton <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[#6]1~[#6]~[#6]~[#6]2~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4[#6]~3~[#6]~[#6]~[#6]~2[#6]~1",
                       scriptPath = scriptPath)[[1]]
}

# (2) Search sterol skeleton derivative
# .searchSterolSkeleton_Derivative(smi = "C1C[C@H](O)C(C)(C)[C@]2([H])CCC3[C@]4(C)CC[C@]([H])([C@H](C)CCC(=C)C(C)C)[C@@]4(C)CCC=3[C@@]12C",
#                                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSterolSkeleton_Derivative <- function(smi, scriptPath){
  sterolSkeleton_position <- .searchSterolSkeleton(smi = smi, scriptPath = scriptPath)
  lapply(sterolSkeleton_position, function(x) {
    c(x[1],
      .TraverseMolecule(smi = smi,
                        start_atom_idx = x[2],
                        non_traversable_atom_idx = c(x[1], x[9]),
                        scriptPath = scriptPath),
      x[9])
  })
}

# (3) Search sterol skeleton Chain 1
# .searchSterolSkeleton_Chain1(smi = "C1[C@H](OC(=O)CCCCCCCCCCC)CC2=CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@]([H])(C)CCCC(C)C)[C@@]4(C)CC[C@]3([H])[C@@]2(C)C1",
#                              scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSterolSkeleton_Chain1 <- function(smi, scriptPath){
  sterolSkeleton_position <- .searchSterolSkeleton(smi = smi, scriptPath = scriptPath)
  lapply(sterolSkeleton_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# (4) Search sterol skeleton Chain 2
# .searchSterolSkeleton_Chain2(smi = "C1[C@]2(C)[C@@]3([H])CC[C@]4(C)[C@@]([H])([C@]([H])(C)CC/C(=C/CC)/C(C)C)CC[C@@]4([H])[C@]3([H])CC=C2C[C@@H](O)C1",
#                              scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSterolSkeleton_Chain2 <- function(smi, scriptPath){
  sterolSkeleton_position <- .searchSterolSkeleton(smi = smi, scriptPath = scriptPath)
  lapply(sterolSkeleton_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[9],
                      non_traversable_atom_idx = x[-9],
                      scriptPath = scriptPath)
  })
}
