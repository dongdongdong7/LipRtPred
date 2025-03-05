#' @rdname LipRtPred_MD
#' @title LipRtPred Molecule Descriptors
#' @description
#' Calculate LipRtPred MD using RDKit package with SMARTS rules.
#' @details
#' \code{smartsMatch()} function is used to find the SMARTS pattern.\cr
#' \code{C_C_count()} function can calculate the atom number between two C atom.
#'

#' @param smiles smiles vector
#' @param SMARTS SMARTS pattern
#'
#' @return A list.
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
# This function is smartsMatch's parallel version
#' @rdname LipRtPred_MD
#' @examples
#' .smartsMatch(smiles = c("CC(O)CC(O)C", "C(O)CC(O)C"), SMARTS = "CO",
#'              scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.smartsMatch <- function(smiles, SMARTS, scriptPath){
  reticulate::source_python(scriptPath)
  matchRes <- GetSubstructMatches_py(smis = as.list(smiles), SMARTS = SMARTS)
  lapply(matchRes, function(x) {
    lapply(x, function(y) unlist(y))
  })
}

#' @rdname LipRtPred_MD
#' @param smi A SMILES string
#' @param start_atom_idx The index of start atom
#' @param end_atom_idx The index of end atom
#'
#' @return Atom number.
#' @export
#'
#' @examples
#' draw_smiles(smiles = "C(O)(=O)CC(O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC", SMARTS = "[CX3;$(C=O);$(C-O)](=O)O", subImgSize = c(800, 800))
#' smartsMatch(smiles = "C(O)(=O)CC(O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC", SMARTS = "[CX3;$(C=O);$(C-O)](=O)O")
#' C_C_count("C(O)(=O)CC(O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC", start_atom_idx = 0, end_atom_idx = 7)
C_C_count <- function(smi, start_atom_idx, end_atom_idx){
  reticulate::source_python(system.file("python", "SMARTS.py", package = "LipRtPred"))
  atom_count <- C_C_count_py(smi = smi, start_atom_idx = as.integer(start_atom_idx), end_atom_idx = as.integer(end_atom_idx))
  return(atom_count)
}
#' @rdname LipRtPred_MD
#' @examples
#' .C_C_count("C(O)(=O)CC(O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC", start_atom_idx = 0, end_atom_idx = 7,
#'            scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.C_C_count <- function(smi, start_atom_idx, end_atom_idx, scriptPath){
  reticulate::source_python(scriptPath)
  atom_count <- C_C_count_py(smi = smi, start_atom_idx = as.integer(start_atom_idx), end_atom_idx = as.integer(end_atom_idx))
  return(atom_count)
}

# .getAtomSymbol(smi = "CCC=CCC(O)OC=CC=CCCCCCCCC(=O)O", atom_idx = c(6, 19, 7),
#                scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.getAtomSymbol <- function(smi, atom_idx, scriptPath){
  reticulate::source_python(scriptPath)
  getAtomSymbol_py(smi = smi, atom_idx_list = as.list(as.integer(atom_idx)))
}

# Find atoms position of fatty acyl's main chain
# min_C: Minimum value of main chain length to be considered
# max_C: Maximum value of main chain length to be considered
# .searchCOO(smi = "C(OC(=O)CCCCCCCCCCCCCCCCC)[C@]([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O",
#            scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
# .searchCOO(smi = "CCC=CCC(O)OC=CC=CCCCCCCCC(=O)O",
#            scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
# .searchCOO(smi = "C(OCCCCCCCC(CC)(CCC)CCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C)=O)COC(=O)CCCCCCC(CC)CCCCC(C)C",
#            scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.searchCOO <- function(smi, min_C = 1, max_C = 24, scriptPath){
  if(min_C <= 0 | max_C > 50 | min_C >= max_C){
    stop("Please reset min_C or max_C!")
  }
  pattern <- sapply((min_C - 1):(max_C - 1), function(x) paste(rep("[C,OH0;!$(O=C);!$(O-[CX3;$(C=O);$(C-O)])]~", x), collapse = ""))
  pattern <- paste0(pattern, "[CX3;$(C=O);$(C-O)]")
  chains_lenList <- lapply(pattern, function(x) {
    searchList <- .smartsMatch(smi, x, scriptPath = scriptPath)[[1]]
    if(length(searchList) == 0) chains_len <- 0
    else chains_len <- sapply(searchList, length)
    return(chains_len)
  })
  chain_num <- max(sapply(chains_lenList, length))
  chains_len <- sapply(chain_num:1, function(i) {
    max(purrr::map_int(chains_lenList, ~ .x[i]), na.rm = TRUE)
  })
  if(all(sapply(chains_len, function(len) len == 0))) return(list())
  chains_len <- sort(chains_len, decreasing = TRUE)
  pattern_target <- sapply(chains_len - 1, function(x) paste(rep("[C,OH0;!$(O=C);!$(O-[CX3;$(C=O);$(C-O)])]~", x), collapse = ""))
  pattern_target <- paste0(pattern_target, "[CX3;$(C=O);$(C-O)]")
  idx_target <- c()
  idx_target_list <- list()
  for(pattern_tmp in pattern_target){
    idxList <- .smartsMatch(smi, pattern_tmp, scriptPath = scriptPath)[[1]]
    for(idx in idxList){
      if(any(idx %in% idx_target)) next
      else{
        idx_target <- c(idx_target, idx)
        idx_target_list <- append(idx_target_list, list(idx))
      }
    }
  }
  return(idx_target_list)
}

# Find atom position of start carbon of branch carbon chain
# .searchBranch(smi = "C(OCCCCCCCC(CC)(CCC)CCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C)=O)COC(=O)CCCCCCC(CC)CCCCC(C)C",
#               scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
# .searchBranch(smi = "C(OCCCCCCCC(CC)(CCC)CCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C)=O)COC(=O)CCCCCCC(OCCC)CCCCC(C)C",
#               scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.searchBranch <- function(smi, scriptPath){
  branch_C_position <- .smartsMatch(smiles = smi,
                                  SMARTS = "[C;$([CH0,CH1]);$([$(C(C)(C)(C)),$(C(C)(C)(C)(C))])]",
                                  scriptPath = scriptPath)[[1]]
  branch_O_position <- .smartsMatch(smiles = smi,
                                    SMARTS = "[O;$(O(C)C(C)(C));!$(O-[CX3;$(C=O);$(C-O)])]",
                                    scriptPath = scriptPath)[[1]]
  branch_O_position <- .smartsMatch(smiles = smi,
                                    SMARTS = "[C;$([CH0,CH1]);$(C(C)(C)[O;$(O(C)C(C)(C));!$(O-[CX3;$(C=O);$(C-O)])]C)]",
                                    scriptPath = scriptPath)[[1]]
  branch_position <- append(branch_C_position, branch_O_position)
  return(branch_position)
}

# Traverse all carbon atoms starting from a specific atom and count them."
# start_atom_idx: The idx of start atom
# non_traversable_atom_ids: Vector of idx of atoms that are forbidden to be traversed
# FA_main_C_position <- .searchCOO(smi = "C(OCCCCCCCC(CC)(CCC)CCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C)=O)COC(=O)CCCCCCC(OCCC)CCCCC(C)C",
#                                  scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
# FA_branch_start_C_position <- .searchBranch(smi = "C(OCCCCCCCC(CC)(CCC)CCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C)=O)COC(=O)CCCCCCC(OCCC)CCCCC(C)C",
#                                             scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
# Select an atom on FA main chain
# .walk_away(smi = "C(OCCCCCCCC(CC)(CCC)CCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C)=O)COC(=O)CCCCCCC(OCCC)CCCCC(C)C",
#            start_atom_idx = FA_branch_start_C_position[[5]],
#            non_traversable_atom_ids = unlist(FA_main_C_position)[unlist(FA_main_C_position) != FA_branch_start_C_position[[5]]],
#            scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.walk_away <- function(smi, start_atom_idx, non_traversable_atom_ids, scriptPath){
  reticulate::source_python(scriptPath)
  count <- walk_away_py(smi = smi, start_atom_idx = start_atom_idx, non_traversable_atom_ids = non_traversable_atom_ids)
  return(count)
}
# Calculate C number of FA chains.
# .cal_c(smi = "C(OC(=O)CCCCCCC(CC)(CCC)CCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C))COC(=O)CCCCCCC(OCCC)CCCCC(C)C",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
# .cal_c(smi = "C(OCCCCCCCCCCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C)=O)COC(=O)CCCCCCCCCCCC(C)C",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_c <- function(smi, min_C = 1, max_C = 24, scriptPath){
  FA_position <- .searchCOO(smi, min_C = min_C, max_C = max_C, scriptPath = scriptPath)
  if(length(FA_position) == 0) return(0)
  branch_C_position <- .searchBranch(smi = smi, scriptPath = scriptPath)
  branch_C_position <- lapply(branch_C_position, function(x){
    if(x %in% unlist(FA_position)) return(x)
    else return(NULL)
  })
  branch_C_position <- branch_C_position[!sapply(branch_C_position, is.null)]
  if(length(branch_C_position) == 0) strand_C_num <- 0
  else{
    branch_C_num <- sum(sapply(branch_C_position, function(i) {
      .walk_away(smi = smi, start_atom_idx = i, non_traversable_atom_ids = unlist(FA_position)[unlist(FA_position) != i], scriptPath = scriptPath)
    }))
  }
  FA_position <- unlist(FA_position)
  FA_symbols <- .getAtomSymbol(smi = smi, atom_idx = FA_position, scriptPath = scriptPath)
  FA_C_num <- sum(sapply(FA_symbols, function(x){
    if(x == "C") return(1)
    else return(0)
  }))
  return(as.integer(FA_C_num + branch_C_num))
}
# Calculate = number on FA's C-Chains.
# .cal_d(smi = "C(OC(=O)CCCCCCC/C=C\\CCCCCC)[C@]([H])(OC(CCCCCCC/C=C\\CCCCCC)=O)COC(CCCCCCC/C=C\\CCCCCC)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_d <- function(smi, min_C = 1, max_C = 24, scriptPath){
  FA_position <- .searchCOO(smi = smi, min_C = min_C, max_C = max_C, scriptPath = scriptPath)
  if(length(FA_position) == 0) return(0)
  FA_position <- unlist(FA_position)
  matchList <- .smartsMatch(smiles = smi, SMARTS = "C=C", scriptPath = scriptPath)[[1]]
  if(length(matchList) == 0) return(0)
  matchNum <- sapply(matchList, function(x) {
    if(all(x %in% FA_position)) return(1)
    else return(0)
  })
  return(sum(matchNum))
}
# Calculate OH number on FA's C-Chains
# .cal_h(smi = "C(CC/C=C\\C/C=C\\CC(O)C(O)C/C=C\\C/C=C\\C/C=C\\CC)(=O)O",
#        min_C = 1, max_C = 24,
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
# .cal_h(smi = "CCC=CCC(O)OC=CC=CCCCCCCCC(=O)O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_h <- function(smi, min_C = 1, max_C = 24, scriptPath){
  FA_position <- .searchCOO(smi = smi, min_C = min_C, max_C = max_C, scriptPath = scriptPath)
  if(length(FA_position) == 0) return(0)
  FA_position <- unlist(FA_position)
  matchList <- .smartsMatch(smiles = smi, SMARTS = "[CX4;$(C-[OH])]", scriptPath = scriptPath)[[1]]
  if(length(matchList) == 0) return(0)
  matchNum <- sapply(matchList, function(x) {
    if(all(x %in% FA_position)) return(1)
    else return(0)
  })
  return(sum(matchNum))
}

# Calculate OH position on FA's C-Chains
# .cal_n(smi = "C(OC(=O)CCCCCCCCCCCCCC(O)CCC)[C@]([H])(OC(CCCCCC(O)CC(O)CCCCCCC)=O)COC(CCCCCCCCCCC)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
# .cal_n(smi = "CCC=CCC(O)OC=CC=CCCCCCCCC(=O)O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_n <- function(smi, min_C = 1, max_C = 24, max_OH = 5,scriptPath){
  FA_position <- .searchCOO(smi = smi, min_C = min_C, max_C = max_C, scriptPath = scriptPath)
  if(length(FA_position) == 0) return(0)
  FAC_position <- .smartsMatch(smiles = smi, SMARTS = "[CX3;$(C=O);$(C-O)]", scriptPath = scriptPath)[[1]]
  #FA_position <- unlist(FA_position)
  FAOH_position <- .smartsMatch(smiles = smi, SMARTS = "[CX4;$(C-[OH])]", scriptPath = scriptPath)[[1]]
  if(length(FAOH_position) == 0) return(rep(0, max_OH))
  # .C_C_count("C(O)C(=O)OC") will return 0
  # distance can not be 0 so + 1
  distance <- sapply(FAOH_position, function(x) {
    chain_rightNow <- FA_position[sapply(FA_position, function(y) {
      x %in% y
    })][[1]]
    C_chain_rightNow <- FAC_position[sapply(FAC_position, function(y) {
      y %in% chain_rightNow
    })][[1]]
    .C_C_count(smi = smi, start_atom_idx = x, end_atom_idx = C_chain_rightNow, scriptPath = scriptPath)
  }) + 1
  distance <- sort(distance)[1:max_OH]
  names(distance) <- paste0("n", 1:max_OH)
  return(distance)
}

# .cal_o(smi = "CCC(C(=O)OCC)CC(=O)OC",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_o <- function(smi, scriptPath){
  o_position <- .smartsMatch(smiles = smi, SMARTS = "C(C)(C)[CX3;$(C=O);$(C-O)](=O)O",
                             scriptPath = scriptPath)[[1]]
  if(length(o_position) == 0) return(0)
  else return(length(o_position))
}

# .cal_a(smi = "CCC(OC(=O)CC)CC(=O)OC",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_a <- function(smi, scriptPath){
  a_position <- .smartsMatch(smiles = smi, SMARTS = "C(C)(C)O[CX3;$(C=O);$(C-O)](=O)",
                             scriptPath = scriptPath)[[1]]
  if(length(a_position) == 0) return(0)
  else return(length(a_position))
}

# Search glycerol position
# .searchGlycerol(smi = "C(OC(=O)CCCCCCCCCCCCCCCCC)[C@]([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O",
#                 scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.searchGlycerol <- function(smi, scriptPath){
  return(.smartsMatch(smiles = smi, SMARTS = "[CH2;$(C-O)]-[CH;$(C-O)]-[CH2;$(C-O)]", scriptPath = scriptPath)[[1]])
}
# Search glycerol FA
# .cal_alpha(smi = "C(OC(=O)CCCCCCCCCCCCCCCCC)[C@]([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O",
#            scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_alpha <- function(smi, min_C = 1, max_C = 24, scriptPath){
  alpha_gama <- .smartsMatch(smiles = smi, SMARTS = "[CX3;$(C=O);$(C-O)](=O)O[CH2;$(C-O)]-[CH;$(C-O)]-[CH2;$(C-O)]", scriptPath = scriptPath)[[1]]
  if(length(alpha_gama) == 0){
    alpha <- 0;gama <- 0
  }else if(length(alpha_gama) == 1){
    alpha <- 1;gama <- 0
  }else if(length(alpha_gama) == 2){
    alpha <- 1;gama <- 1
  }
  beta_position <- .smartsMatch(smiles = smi, SMARTS = "[CH2;$(C-O)]-[CH;$(C-O)](OC(=O))-[CH2;$(C-O)]", scriptPath = scriptPath)[[1]]
  if(length(beta_position) == 0) beta <- 0
  else if(length(beta_position)) beta <- 1
  return(c(alpha = alpha, beta = beta, gama = gama))
}
# .cal_m(smi = "O=C(O[C@](CC([O-])=O)(C[N+](C)(C)C)[H])CCCCCCC",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_m <- function(smi, scriptPath){
  car_position <- .smartsMatch(smiles = smi, SMARTS = "C(C[N+](C)(C)(C))CC(=O)[O-]", scriptPath = scriptPath)[[1]]
  if(length(car_position) == 0) return(0)
  else return(length(car_position))
}
# .cal_f(smi = "CCC=CCCOC=CC=CCCCCCCCC(=O)O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_f <- function(smi, scriptPath){
  alkether_position <- .smartsMatch(smiles = smi, SMARTS = "C=C-O-C", scriptPath = scriptPath)[[1]]
  if(length(alkether_position) == 0) return(0)
  else return(length(alkether_position))
}
# .cal_j(smi = "CCC=CCCOCCC=CCCCCCCCC(=O)O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_j <- function(smi, scriptPath){
  ether_position <- .smartsMatch(smiles = smi, SMARTS = "[CX4]-O-[CX4]", scriptPath = scriptPath)[[1]]
  if(length(ether_position) == 0) return(0)
  else return(length(ether_position))
}
# .cal_k(smi = "C1[C@H](OC(=O)CCCCCCCCCCCCC)CC2=CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@]([H])(C)CCCC(C)C)[C@@]4(C)CC[C@]3([H])[C@@]2(C)C1",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_k <- function(smi, scriptPath){
  cholesterol_position <- .smartsMatch(smiles = smi,
                                       SMARTS = "C[CH](CCCC(C)C)[CH]1CC[CH]2[C]1(CC[CH]3[CH]2CC=C4[C]3(CC[CH](C4)O)C)C",
                                       scriptPath = scriptPath)[[1]]
  if(length(cholesterol_position) == 0) return(0)
  else return(length(cholesterol_position))
}
# .cal_q(smi = "[C@](COP(=O)([O])OCC[N+](C)(C)C)([H])(O)CO/C=C\\CCCCCCCCCCCC",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_q <- function(smi, scriptPath){
  choline_position <- .smartsMatch(smiles = smi,
                                   SMARTS = "OCC[N+](C)(C)(C)",
                                   scriptPath = scriptPath)[[1]]
  if(length(choline_position) == 0) return(0)
  else return(length(choline_position))
}
# .cal_r(smi = "[C@](COP(=O)(O)OCCN(C)C)([H])(OC(CCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_r <- function(smi, scriptPath){
  r_position <- .smartsMatch(smiles = smi,
                             SMARTS = "OCCN(C)(C)",
                             scriptPath = scriptPath)[[1]]
  if(length(r_position) == 0) return(0)
  else return(length(r_position))
}
# .cal_u(smi = "[C@](COP(=O)(O)OCC)([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_u <- function(smi, scriptPath){
  u_position <- .smartsMatch(smiles = smi,
                             SMARTS = "O[CH2][CH3]",
                             scriptPath = scriptPath)[[1]]
  if(length(u_position) == 0) return(0)
  else return(length(u_position))
}
# .cal_e(smi = "[C@](COP(=O)(O)OCCN)([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_e <- function(smi, scriptPath){
  e_position <- .smartsMatch(smiles = smi,
                             SMARTS = "O[CH2][CH2][NH2]",
                             scriptPath = scriptPath)[[1]]
  if(length(e_position) == 0) return(0)
  else return(length(e_position))
}
# .cal_g(smi = "C1CCC(C)=C(/C=C/C(=C/C=C/C(=C/C(=O)O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](C(=O)O)O2)/C)/C)C1(C)C",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_g <- function(smi, scriptPath){
  g_position <- .smartsMatch(smiles = smi,
                             SMARTS = "O=C(O)[CH]1O[CH]([CH](O)[CH](O)[CH]1O)O",
                             scriptPath = scriptPath)[[1]]
  if(length(g_position) == 0) return(0)
  else return(length(g_position))
}
# .cal_v(smi = "C(OC(=O)CCCCCCCCCCCCCCCCC)[C@]([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_v <- function(smi,scriptPath){
  v_position <- .smartsMatch(smiles = smi,
                             SMARTS = "[CH2;$(C-O)]-[CH;$(C-O)]-[CH2;$(C-O)]",
                             scriptPath = scriptPath)[[1]]
  if(length(v_position) == 0) return(0)
  else return(length(v_position))
}
# .cal_w(smi = "C[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O)C)C",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_w <- function(smi, scriptPath){
  w_position <- .smartsMatch(smiles = smi,
                             SMARTS = "C(O)[CH]1O[CH]([CH](O)[CH](O)[CH]1O)O",
                             scriptPath = scriptPath)[[1]]
  if(length(w_position) == 0) return(0)
  else return(length(w_position))
}
# .cal_i(smi = "OC(=O)C/C(/C)=C/C=C/CC(C)C(O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_i <- function(smi, scriptPath){
  i_position <- .smartsMatch(smiles = smi,
                             SMARTS = "C1(C(C(C(C(C1O)O)O)O)O)O",
                             scriptPath = scriptPath)[[1]]
  if(length(i_position) == 0) return(0)
  else return(length(i_position))
}
# .cal_x(smi = "[C@](COP(=O)(O)OCCNC)([H])(OC(CCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCC)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_x <- function(smi, scriptPath){
  x_position <- .smartsMatch(smiles = smi,
                             SMARTS = "O[CH2][CH2][NH][CH3]",
                             scriptPath = scriptPath)[[1]]
  if(length(x_position) == 0) return(0)
  else return(length(x_position))
}
# .cal_y(smi = "OC[C@H]1O[C@@H](O[C@H]2[C@@H](O)[C@H](O[C@@H]2CO)O[C@@H]3[C@@H](O)[C@H](O[C@@H]3CO)O)[C@H](O)[C@@H](NC(C)=O)C1=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_y <- function(smi, scriptPath){
  y_position <- .smartsMatch(smiles = smi,
                             SMARTS = "OC[C@H]1O[C@@H](O[C@H]2[C@@H](O)[C@H](O[C@@H]2CO)O[C@@H]3[C@@H](O)[C@H](O[C@@H]3CO)O)[C@H](O)[C@@H](NC(C)=O)C1=O",
                             scriptPath = scriptPath)[[1]]
  if(length(y_position) == 0) return(0)
  else return(length(y_position))
}
# .cal_z(smi = "CCCCC=CC=CC(CCCCCCCC(=O)O)OO",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_z <- function(smi, scriptPath){
  z_position <- .smartsMatch(smiles = smi,
                             SMARTS = "O-[OH]",
                             scriptPath = scriptPath)[[1]]
  if(length(z_position) == 0) return(0)
  else return(length(z_position))
}
# .cal_p(smi = "[C@](COP(=O)(O)OCCNC)([H])(OC(CCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCC)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_p <- function(smi, scriptPath){
  p_position <- .smartsMatch(smiles = smi,
                             SMARTS = "[PX4;$(P=O);$(P(-O)-O)]",
                             scriptPath = scriptPath)[[1]]
  if(length(p_position) == 0) return(0)
  else return(length(p_position))
}
# .cal_delta(smi = "[C@](CO)([H])(NC(CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCOC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)=O)[C@]([H])(O)[C@H](O)CCCCCCCCCCCCCCCCCC",
#            scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_delta <- function(smi, scriptPath){
  delta_position <- .smartsMatch(smiles = smi,
                                 SMARTS = "CCCCCCCCCCCCCC[C@H]([C@H]([C@H](CO)N)O)O",
                                 scriptPath = scriptPath)[[1]]
  if(length(delta_position) == 0) return(0)
  else return(length(delta_position))
}
# .cal_epsilon(smi = "C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
#              scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_epsilon <- function(smi, scriptPath){
  epsilon_position <- .smartsMatch(smiles = smi,
                                   SMARTS = "C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
                                   scriptPath = scriptPath)[[1]]
  if(length(epsilon_position) == 0) return(0)
  else return(length(epsilon_position))
}
# .cal_theta(smi = "C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
#            scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_theta <- function(smi, scriptPath){
  theta_position <- .smartsMatch(smiles = smi,
                                 SMARTS = "C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
                                 scriptPath = scriptPath)[[1]]
  if(length(theta_position) == 0) return(0)
  else return(length(theta_position))
}
# .cal_xi(smi = "CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)C(C)C)C)C",
#         scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_xi <- function(smi, scriptPath){
  xi_position <- .smartsMatch(smiles = smi,
                              SMARTS = "CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)C(C)C)C)C",
                              scriptPath = scriptPath)[[1]]
  if(length(xi_position) == 0) return(0)
  else return(length(xi_position))
}
# .cal_pi(smi = "C1[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)O2)CC2=CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@@]([H])(/C=C/[C@@H](CC)C(C)C)C)[C@@]4(C)CC[C@]3([H])[C@@]2(C)C1",
#         scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_pi <- function(smi, scriptPath){
  pi_position <- .smartsMatch(smiles = smi,
                              SMARTS = "CC[C@H](/C=C/[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C",
                              scriptPath = scriptPath)[[1]]
  if(length(pi_position) == 0) return(0)
  else return(length(pi_position))
}
# .cal_s(smi = "C(O)[C@](C(=O)O)([H])NC(C[C@@H](O)CCCCCCC)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_s <- function(smi, scriptPath){
  s_position <- .smartsMatch(smiles = smi,
                             SMARTS = "C([CH](C(=O)O)N)O",
                             scriptPath = scriptPath)[[1]]
  if(length(s_position) == 0) return(0)
  else return(length(s_position))
}

# .cal_eta(smi = "[C@](CO)([H])(NC(CCCCCCCCCCCCCCCCC)=O)[C@]([H])(O)CCCCCCCCCCCCCCC",
#          scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_eta <- function(smi, scriptPath){
  eta_position <- .smartsMatch(smiles = smi,
                               SMARTS = "CCCCCCCCCCCCCCC[C@H]([C@H](CO)N)O",
                               scriptPath = scriptPath)[[1]]
  if(length(eta_position) == 0) return(0)
  else return(length(eta_position))
}
# .cal_miu(smi = "[C@](CO)([H])(NC(C(O)CCCCCCCCCCCCCCCC)=O)[C@]([H])(O)/C=C/CCCCCCCCCCCCC",
#          scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_miu <- function(smi, scriptPath){
  miu_position <- .smartsMatch(smiles = smi,
                               SMARTS = "CCCCCCCCCCCCCC=C[CH]([CH](CO)N)O",
                               scriptPath = scriptPath)[[1]]
  if(length(miu_position) == 0) return(0)
  else return(length(miu_position))
}
# .cal_omega(smi = "C1[C@]2(C)[C@@]3([H])CC[C@]4(C)[C@@]([H])([C@]([H])(C)CCCC(C)C)CC[C@@]4([H])[C@]3([H])CC[C@]2([H])C[C@H](O)C1",
#            scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_omega <- function(smi, scriptPath){
  omega_position <- .smartsMatch(smiles = smi,
                                 SMARTS = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CC[C@@H](C4)O)C)C",
                                 scriptPath = scriptPath)[[1]]
  if(length(omega_position) == 0) return(0)
  else return(length(omega_position))
}
# .cal_lambda(smi = "CC[C@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)C(C)C)C)C",
#             scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_lambda <- function(smi, scriptPath){
  lambda_position <- .smartsMatch(smiles = smi,
                                  SMARTS = "CC[C@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)C(C)C)C)C",
                                  scriptPath = scriptPath)[[1]]
  if(length(lambda_position) == 0) return(0)
  else return(length(lambda_position))
}
# .cal_psi(smi = "C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)OS(=O)(=O)c1ccccc1",
#          scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_psi <- function(smi, scriptPath){
  psi_position <- .smartsMatch(smiles = smi,
                               SMARTS = "C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)OS(=O)(=O)c1ccccc1",
                               scriptPath = scriptPath)[[1]]
  if(length(psi_position) == 0) return(0)
  else return(length(psi_position))
}
# .cal_phi(smi = "CC(C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)OS(=O)(=O)C5=CC=CC=C5)C)C",
#          scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_phi <- function(smi, scriptPath){
  phi_position <- .smartsMatch(smiles = smi,
                               SMARTS = "c1ccccc1S(=O)(=O)O",
                               scriptPath = scriptPath)[[1]]
  if(length(phi_position) == 0) return(0)
  else return(length(phi_position))
}
# .cal_t(smi = "[C@](COCCC(C(=O)O)[N+](C)(C)C)([H])(O)COC(CCCCCCCCCCCCC)=O",
#        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_t <- function(smi, scriptPath){
  t_position <- .smartsMatch(smiles = smi,
                             SMARTS = "[CH3][N+]([CH3])([CH3])C(CCO)C(=O)O",
                             scriptPath = scriptPath)[[1]]
  if(length(t_position) == 0) return(0)
  else return(length(t_position))
}

# .getLipRtPredMD(smi = "CC(O)CC(=O)OC", scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.getLipRtPredMD <- function(smi, scriptPath){
  browser()
  reticulate::source_python(scriptPath)
}
#GetLipRtPred_MD(cmpDf = cmpDf_demo2)
GetLipRtPred_MD <- function(cmpDf, flavor = "CxSmiles", min_C = 1, max_C = 24, thread = 3){
  browser()
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  wrong_idx <- .check_smiles(smiles = cmpDf$smiles)
  if(length(wrong_idx) != 0){
    message(paste0("Wrong Idx: ", paste0(wrong_idx, collapse = " ")))
    cmpDf <- cmpDf[-wrong_idx, ]
  }
}
