#' @title Match substructure using SMARTS
#'
#' @param smiles smiles vector.
#' @param SMARTS SMARTS pattern.
#'
#' @return list.
#' @export
#'
#' @examples
#' smartsMatch(smiles = c("CC(O)CC(O)C", "C(O)CC(O)C"), SMARTS = "CO")
smartsMatch <- function(smiles, SMARTS){
  reticulate::source_python(system.file("python", "SMARTS.py", package = "LipRtPred"))
  matchRes <- GetSubstructMatches_py(smis = as.list(smiles), SMARTS = SMARTS)
  lapply(matchRes, function(x) {
    lapply(x, function(y) unlist(y))
  })
}

#' @title C_C_count
#' @description
#' Atom number between two C atom.
#'
#' @param smiles A smiles string.
#' @param start_atom_idx A number.
#' @param end_atom_idx A number.
#'
#' @return Atom number.
#' @export
#'
#' @examples
#' draw_smiles(smiles = "C(O)(=O)CC(O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC", SMARTS = "[CX3;$(C=O);$(C-O)](=O)O", subImgSize = c(800, 800))
#' smartsMatch(smiles = "C(O)(=O)CC(O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC", SMARTS = "[CX3;$(C=O);$(C-O)](=O)O")
#' C_C_count("C(O)(=O)CC(O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC", start_atom_idx = 0, end_atom_idx = 7)
C_C_count <- function(smiles, start_atom_idx, end_atom_idx){
  reticulate::source_python(system.file("python", "SMARTS.py", package = "LipRtPred"))
  atom_count <- C_C_count_py(smi = smiles, start_atom_idx = as.integer(start_atom_idx), end_atom_idx = as.integer(end_atom_idx))
  return(atom_count)
}
