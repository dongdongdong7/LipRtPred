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

# This function is used to search FA's C position.
#' @rdname LipRtPred_MD
#' @param min_C Minimum value of C-chain length to be considered
#' @param max_C Maximum value of C-chain length to be considered
#' @examples
#' .searchCOO(smi = "C(OC(=O)CCCCCCCCCCCCCCCCC)[C@]([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O",
#'            scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.searchCOO <- function(smi, min_C = 1, max_C = 24, scriptPath){
  if(min_C <= 0 | max_C > 50 | min_C >= max_C){
    stop("Please reset min_C or max_C!")
  }
  pattern <- sapply((min_C - 1):(max_C - 1), function(x) paste(rep("C~", x), collapse = ""))
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
  pattern_target <- sapply(chains_len - 1, function(x) paste(rep("C~", x), collapse = ""))
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

# Calculate strand C number on FA C-Chains.
.walk_away <- function(smi, start_atom_idx, main_chains_atom_idx, scriptPath){
  reticulate::source_python(scriptPath)
  count <- walk_away_py(smi = smi, start_atom_idx = start_atom_idx, main_chains_atom_idx = main_chains_atom_idx)
  return(count)
}
# Calculate C number of FA C-Chains.
#' @rdname LipRtPred_MD
#' @examples
#' .cal_c(smi = "C(OC(=O)CCCCCCCCCCCCCCCCC)[C@]([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O",
#'        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
#' .cal_c(smi = "C(OCCCCCCCCCCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C)=O)COC(=O)CCCCCCCCCCCC(C)C",
#'        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
.cal_c <- function(smi, min_C = 1, max_C = 24, scriptPath){
  FA_position <- .searchCOO(smi, min_C = min_C, max_C = max_C, scriptPath = scriptPath)
  if(length(FA_position) == 0) return(0)
  strand_C_position <- .smartsMatch(smiles = smi, SMARTS = "[C;$([CH0,CH1]);$(C(C)C);!$(C-O);!$(C-N);!$(C-P);!$(C-S)]", scriptPath = scriptPath)[[1]]
  strand_C_position <- lapply(strand_C_position, function(x){
    if(x %in% unlist(FA_position)) return(x)
    else return(NULL)
  })
  strand_C_position <- strand_C_position[!sapply(strand_C_position, is.null)]
  if(length(strand_C_position) == 0) strand_C_num <- 0
  else{
    strand_C_num <- sum(sapply(strand_C_position, function(i) {
      .walk_away(smi = smi, start_atom_idx = i, main_chains_atom_idx = unlist(FA_position)[unlist(FA_position) != i], scriptPath = scriptPath)
    }))
  }
  return(as.integer(sum(sapply(FA_position, length))) + strand_C_num)
}
# Calculate = number on FA's C-Chains.
#' @rdname LipRtPred_MD
#' @examples
#' .cal_d(smi = "C(OC(=O)CCCCCCC/C=C\\CCCCCC)[C@]([H])(OC(CCCCCCC/C=C\\CCCCCC)=O)COC(CCCCCCC/C=C\\CCCCCC)=O",
#'        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
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
#' @rdname LipRtPred_MD
#' @examples
#' .cal_h(smi = "C(CC/C=C\\C/C=C\\CC(O)C(O)C/C=C\\C/C=C\\C/C=C\\CC)(=O)O",
#'        min_C = 1, max_C = 24,
#'        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
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
#' @rdname LipRtPred_MD
#' @examples
#' .cal_n(smi = "C(OC(=O)CCCCCCCCCCCCCC(O)CCC)[C@]([H])(OC(CCCCCC(O)CC(O)CCCCCCC)=O)COC(CCCCCCCCCCC)=O",
#'        scriptPath = system.file("python", "SMARTS.py", package = "LipRtPred"))
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
  FA_position <- .searchCOO(smi = smi, min_C = min_C, max_C = max_C, scriptPath = scriptPath)
  if(length(FA_position) == 0) return(c(alpha = 0, beta = 0, gama = 0))
  Clycerol_position <- .searchGlycerol(smi = smi, scriptPath = scriptPath)
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
