#' @title get_atom_distance
#' @description
#' Calculate 3D distance between two atom.
#'
#' @param smiles A smiles string.
#' @param start_atom_idx Start atom index.
#' @param end_atom_idx End atom index.
#'
#' @return A distance value.
#' @export
#'
#' @examples
#' smiles <- "C(OCC)(/C=C/C(O)CCCCC)C(O)CCCCCCCC(=O)O"
#' get_atom_distance(smiles, start_atom_idx = 6, end_atom_idx = 22)
#' get_atom_distance(smiles, start_atom_idx = 13, end_atom_idx = 22)
get_atom_distance <- function(smiles, start_atom_idx, end_atom_idx){
  reticulate::source_python(system.file("python", "3D.py", package = "LipRtPred"))
  get_atom_distance_py(smi = smiles, atom_idx1 = as.integer(start_atom_idx), atom_idx2 = as.integer(end_atom_idx))
}
