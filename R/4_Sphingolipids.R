# Search substructure of Sphingolipids (SP)
# Barry Song
# 250312

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py

# 1. Fatty Chain Head
# (1) Search sphingosine(homologs and variants)
# .searchSphingosine(smi = "OC[C@]([H])(N)[C@]([H])(O)/C=C/CCCCCCCCCCCCC",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchSphingosine(smi = "[C@](CO)([H])(N)C(=O)/C=C/CCCCCCCCCCCCC",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSphingosine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)C(N)[C;$(C-O),$(C=O)](~[O;$(O-C),$(O=C)])[C;$(C-C),$(C=C)]~C",
                       scriptPath = scriptPath)[[1]]
}

# (2) Search sphingosine C chain
# .searchSphingosineC_Chain(smi = "[C@](CO)([H])(N)C(=O)/C=C/CCCCCCCCCCCCC",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSphingosineC_Chain <- function(smi, scriptPath){
  sphingosine_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
  lapply(sphingosine_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[3],
                      non_traversable_atom_idx = x[c(1:2, 4)],
                      scriptPath = scriptPath)
  })
}
# (3) Search sphingosie N chain
# .searchSphingosineN_Chain(smi = "P(=O)(OC[C@]([H])(NC(=O)CCCCCCCCCCC)[C@]([H])(O)/C=C/CCCCCCCCCCCCC)(O)O",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSphingosineN_Chain <- function(smi, scriptPath){
  sphingosine_position <- .searchSphingosine(smi = smi, scriptPath = scriptPath)
  sphingosieAmide_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
  sphingosieAmide_logical <- sapply(sphingosieAmide_position, function(x) {
    if(x[3] %in% unlist(sphingosine_position)) return(TRUE)
    else return(FALSE)
  })
  sphingosieAmide_position <- sphingosieAmide_position[sphingosieAmide_logical]
  lapply(sphingosieAmide_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# (4) Search spisulosine
# .searchSpisulosine(smi = "C(C)(N)C(O)CCCCCCCCCCC(C)CC(C)CC",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSpisulosine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CH3]C(N)[C;$(C-O),$(C=O)](~[O;$(O-C),$(O=C)])[C;$(C-C),$(C=C)]~C",
                       scriptPath = scriptPath)[[1]]
}

# (5) Search spisulosine C Chain
# .searchSpisulosineC_Chain(smi = "C(C)(N)C(O)CCCCCCCCCCC(C)CC(C)CC",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSpisulosineC_Chain <- function(smi, scriptPath){
  spisulosine_position <- .searchSpisulosine(smi = smi, scriptPath = scriptPath)
  lapply(spisulosine_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[2],
                      non_traversable_atom_idx = x[c(1, 3)],
                      scriptPath = scriptPath)
  })
}
# (6) Search spisulosine N Chain
# .searchSpisulosineN_Chain(smi = "[C@](C)([H])(NC(CCCCCCCCCCCCCCCCCCCCCCC)=O)[C@]([H])(O)CCCCCCCCCCCCCCC",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSpisulosineN_Chain <- function(smi, scriptPath){
  spisulosine_position <- .searchSpisulosine(smi = smi, scriptPath = scriptPath)
  spisulosineAmide_position <- .searchAmide(smi = smi, scriptPath = scriptPath)
  spisulosineAmide_logical <- sapply(spisulosineAmide_position, function(x) {
    if(x[3] %in% unlist(spisulosine_position)) return(TRUE)
    else return(FALSE)
  })
  spisulosineAmide_position <- spisulosineAmide_position[spisulosineAmide_logical]
  lapply(spisulosineAmide_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}
