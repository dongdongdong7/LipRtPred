# Wrapper function of molecule_operation.py
# Barry Song
# 250307

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# addHs: whether add Hs
# scriptPath: path of molecule_operation.py
# smis: SMILES string vector
# SMARTS: SMARTS pattern
# start_atom_idx: start atom index
# end_atom_idx: end atom index
# atom_idx_vector: atom index vector
# non_traversable_atom_idx: vector of atom index that are not allowed to be traversed

# Get the number of atom
# .GetNumAtoms(smi = "CCCCC(=O)O",
#              scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.GetNumAtoms <- function(smi, addHs = FALSE, scriptPath){
  reticulate::source_python(scriptPath)
  GetNumAtoms_py(smi, addHs = addHs)
}

# Get SMILES with atom index
# .Get_SMILES_with_index(smi = "CCCCC(=O)O",
#                        scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.Get_SMILES_with_index <- function(smi, scriptPath){
  reticulate::source_python(scriptPath)
  Get_SMILES_with_index_py(smi)
}

# Searching for substructures of molecules using SMARTS
# .GetSubstructMatches(smis = c("CC(=O)OC", "CCCO"),
#                      SMARTS = "CO",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.GetSubstructMatches <- function(smis, SMARTS, scriptPath){
  reticulate::source_python(scriptPath)
  matchRes <- GetSubstructMatches_py(smis = as.list(smis),
                                     SMARTS = SMARTS)
  lapply(matchRes, function(x) {
    lapply(x, function(y) unlist(y))
  })
}

# Get the shortest path between two atoms
# .GetShortestPath(smi = "CC(CC(C))COCC",
#                  start_atom_idx = 0, end_atom_idx = 6,
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.GetShortestPath <- function(smi, start_atom_idx, end_atom_idx, scriptPath){
  reticulate::source_python(scriptPath)
  path <- GetShortestPath_py(smi = smi,
                             start_atom_idx = as.integer(start_atom_idx), end_atom_idx = as.integer(end_atom_idx))
  unlist(path)
}

# Get symbol of atoms
# .GetAtomSymbol(smi = "CCOCC(O)N",
#                atom_idx_vector = c(0,1,2,3,4,5,6),
#                scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.GetAtomSymbol <- function(smi, atom_idx_vector, scriptPath){
  reticulate::source_python(scriptPath)
  GetAtomSymbol_py(smi = smi, atom_idx_list = as.list(as.integer(atom_idx_vector)))
}

# Get atom index of rings
# .GetRingAtom(smi = "OC[C@]1(OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H](CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@@]([H])(CO)COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC1)C)[H]",
#              scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.GetRingAtom <- function(smi, scriptPath){
  reticulate::source_python(scriptPath)
  RingRes <- GetRingAtom_py(smi = smi)

  lapply(RingRes, function(x) {
    unlist(x)
  })
}

# Traverse the entire molecule starting from the start atom, excluding the specified atoms
# .TraverseMolecule(smi = "CC(CC)C(=O)O",
#                   start_atom_idx = 4, non_traversable_atom_idx = c(5, 6),
#                   scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.TraverseMolecule <- function(smi, start_atom_idx, non_traversable_atom_idx, scriptPath){
  reticulate::source_python(scriptPath)
  TraverseMolecule_py(smi = smi,
                      start_atom_idx = as.integer(start_atom_idx),
                      non_traversable_atom_idx = as.list(as.integer(non_traversable_atom_idx)))
}
