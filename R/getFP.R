# This script is used to calculate fingerprint
# 250221
# Barry Song

# data("cmpDf_demo", package = "LipRtPred")
# .getMorganFP(cmpDf = cmpDf_demo)
.getMorganFP <- function(cmpDf, radius = 2, fpSize = 2048, count = FALSE, thread = 1){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  message("Calculate Morgan FP...")
  reticulate::source_python(system.file("python", "RDKit_FP.py", package = "LipRtPred"))
  fp <- getMorganFP_py(smis = as.list(cmpDf$smiles), radius = as.integer(radius), fpSize = as.integer(fpSize), count = count, thread = as.integer(thread))
  colnames(fp) <- paste0("Morgan", as.integer(colnames(fp)) + 1)
  return(dplyr::as_tibble(cbind(cmpDf, fp)))
}
# data("cmpDf_demo", package = "LipRtPred")
# .getFMorganFP(cmpDf = cmpDf_demo)
.getFMorganFP <- function(cmpDf, radius = 2, fpSize = 2048, count = FALSE, thread = 1){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  message("Calculate Feature Morgan FP...")
  reticulate::source_python(system.file("python", "RDKit_FP.py", package = "LipRtPred"))
  fp <- getFMorganFP_py(smis = as.list(cmpDf$smiles), radius = as.integer(radius), fpSize = as.integer(fpSize), count = count, thread = as.integer(thread))
  colnames(fp) <- paste0("FMorgan", as.integer(colnames(fp)) + 1)
  return(dplyr::as_tibble(cbind(cmpDf, fp)))
}
# data("cmpDf_demo", package = "LipRtPred")
# .getRDKitFP(cmpDf = cmpDf_demo)
.getRDKitFP <- function(cmpDf, minPath = 1, maxPath = 7, useHs = TRUE, fpSize = 2048, count = FALSE, thread = 1){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  message("Calculate RDKit FP...")
  reticulate::source_python(system.file("python", "RDKit_FP.py", package = "LipRtPred"))
  fp <- getRDKitFP_py(smis = as.list(cmpDf$smiles), minPath = as.integer(minPath), maxPath = as.integer(maxPath), useHs = useHs, fpSize = as.integer(fpSize), count = count, thread = as.integer(thread))
  colnames(fp) <- paste0("RDKit", as.integer(colnames(fp)) + 1)
  return(dplyr::as_tibble(cbind(cmpDf, fp)))
}
# data("cmpDf_demo", package = "LipRtPred")
# .getApFP(cmpDf = cmpDf_demo)
.getApFP <- function(cmpDf, minDistance = 1, maxDistance = 30, fpSize = 2048, count = FALSE, thread = 1){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  message("Calculate Atom Pairs FP...")
  reticulate::source_python(system.file("python", "RDKit_FP.py", package = "LipRtPred"))
  fp <- getApFP_py(smis = as.list(cmpDf$smiles), minDistance = as.integer(minDistance), maxDistance = as.integer(maxDistance), fpSize = as.integer(fpSize), count = count, thread = as.integer(thread))
  colnames(fp) <- paste0("Ap", as.integer(colnames(fp)) + 1)
  return(dplyr::as_tibble(cbind(cmpDf, fp)))
}
# data("cmpDf_demo", package = "LipRtPred")
# .getTtFP(cmpDf = cmpDf_demo)
.getTtFP <- function(cmpDf, fpSize = 2048, count = FALSE, thread = 1){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  message("Calculate Topological Torsion FP...")
  reticulate::source_python(system.file("python", "RDKit_FP.py", package = "LipRtPred"))
  fp <- getTtFP_py(smis = as.list(cmpDf$smiles), fpSize = as.integer(fpSize), count = count, thread = as.integer(thread))
  colnames(fp) <- paste0("Tt", as.integer(colnames(fp)) + 1)
  return(dplyr::as_tibble(cbind(cmpDf, fp)))
}
# data("cmpDf_demo", package = "LipRtPred")
# .getMaccsFP(cmpDf = cmpDf_demo)
.getMaccsFP <- function(cmpDf, thread = 1){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  message("Calculate MACCS FP...")
  reticulate::source_python(system.file("python", "RDKit_FP.py", package = "LipRtPred"))
  fp <- getMaccsFP_py(smis = as.list(cmpDf$smiles), thread = as.integer(thread))
  colnames(fp) <- paste0("Maccs", as.integer(colnames(fp)) + 1)
  return(dplyr::as_tibble(cbind(cmpDf, fp)))
}

#' @rdname getFP
#' @title Calculate molecule fingerprint
#' @description
#' Use python package RDKit to calculate molecule fingerprint.
#'
#' @param cmpDf A tibble with two column, id and smiles.
#' @param flavor Default is CxSmiles. More flavor please see \code{\link[rcdk]{smiles.flavors}}.
#' @param radius The radius defines the range of atomic neighbourhoods considered for th fingerprint calculation.
#' @param fpSize Size of the generated fingerprint.
#' @param count Whether it is a counted fingerprint.
#' @param thread Number of threads in parallel.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' data("cmpDf_demo", package = "LipRtPred")
#' fpsDf <- GetMorgan_FP(cmpDf = cmpDf_demo)
GetMorgan_FP <- function(cmpDf, flavor = "CxSmiles", radius = 2, fpSize = 2048, count = FALSE, thread = 1){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  fpsDf <- .getMorganFP(cmpDf = cmpDf, radius = radius, fpSize = fpSize, count = count, thread = thread)
  return(fpsDf)
}
#' @rdname getFP
#' @export
#'
#' @examples
#' fpsDf <- GetFMorgan_FP(cmpDf = cmpDf_demo)
GetFMorgan_FP <- function(cmpDf, flavor = "CxSmiles", radius = 2, fpSize = 2048, count = FALSE, thread = 1){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  fpsDf <- .getFMorganFP(cmpDf = cmpDf, radius = radius, fpSize = fpSize, count = count, thread = thread)
  return(fpsDf)
}
#' @rdname getFP
#'
#' @param minPath The minimum path length (in bonds) to be included.
#' @param maxPath The maximum path length (in bonds) to be included.
#' @param useHs Toggles inclusion of Hs in paths (if the molecule has explicit Hs)
#'
#' @export
#'
#' @examples
#' fpsDf <- GetRDKit_FP(cmpDf = cmpDf_demo)
GetRDKit_FP <- function(cmpDf, flavor = "CxSmiles", minPath = 1, maxPath = 7, useHs = TRUE, fpSize = 2048, count = FALSE, thread = 1){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  fpsDf <- .getRDKitFP(cmpDf = cmpDf, minPath = minPath, maxPath = maxPath, useHs = useHs, fpSize = fpSize, count = count, thread = thread)
  return(fpsDf)
}
#' @rdname getFP
#'
#' @param minDistance minimum distance between atoms to be considered in a pair, default is 1 bond.
#' @param maxDistance maximum distance between atoms to be considered in a pair, default is 7 bond.
#'
#' @export
#'
#' @examples
#' fpsDf <- GetAp_FP(cmpDf = cmpDf_demo)
GetAp_FP <- function(cmpDf, flavor = "CxSmiles", minDistance = 1, maxDistance = 30, fpSize = 2048, count = FALSE, thread = 1){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  fpsDf <- .getApFP(cmpDf = cmpDf, minDistance = minDistance, maxDistance = maxDistance, fpSize = fpSize, count = count, thread = thread)
  return(fpsDf)
}
#' @rdname getFP
#' @export
#'
#' @examples
#' fpsDf <- GetTt_FP(cmpDf = cmpDf_demo)
GetTt_FP <- function(cmpDf, flavor = "CxSmiles", fpSize = 2048, count = FALSE, thread = 1){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  fpsDf <- .getTtFP(cmpDf = cmpDf, fpSize = fpSize, count = count, thread = thread)
  return(fpsDf)
}
#' @rdname getFP
#' @export
#'
#' @examples
#' fpsDf <- GetMaccs_FP(cmpDf = cmpDf_demo)
GetMaccs_FP <- function(cmpDf, flavor = "CxSmiles", thread = 1){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  fpsDf <- .getMaccsFP(cmpDf = cmpDf, thread = thread)
  return(fpsDf)
}
