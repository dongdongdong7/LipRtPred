# This script is used to calculate descriptors
# 250221
# Barry Song

# data("cmpDf_demo", package = "LipRtPred")
# ?rcdk::smiles.flavors
# cmpDf_demo$smiles <- .convertSMILES(smiles = cmpDf_demo$smiles)
.convertSMILES <- function(smiles, flavor = "CxSmiles"){
  if(any(is.na(smiles))){
    na_idx <- which(is.na(smiles))
    stop(paste0("smiles[", na_idx, ", ] is NA."))
  }
  message("Converting SMILES...")
  mols <- rcdk::parse.smiles(smiles)
  pb <- utils::txtProgressBar(max = length(mols), style = 3)
  smiles <- sapply(1:length(mols), function(i) {
    utils::setTxtProgressBar(pb, i)
    rcdk::get.smiles(mols[[i]], rcdk::smiles.flavors(flavors = flavor))
  })
  message("\n")
  return(smiles)
}
# .check_smiles(smiles = cmpDf_demo2$smiles)
.check_smiles <- function(smiles){
  message("Check SMILES...")
  reticulate::source_python(system.file("python", "RDKit_MD.py", package = "LipRtPred"))
  resList <- check_smiles_py(smis = as.list(smiles))
  wrong_idx <- which(sapply(resList, is.null))
  warning(paste0("These SMILES can not be generated to a Mol object!\n", paste(smiles[wrong_idx], collapse = "\n")))
  return(wrong_idx)
}
# data("cmpDf_demo", package = "LipRtPred")
# category <- c("all", "protein", "hybrid", "constitutional", "topological", "electronic", "geometrical")[1]
# .getCDKMD(cmpDf = cmpDf_demo)
.getCDKMD <- function(cmpDf, category = "all", thread = 1){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  message("Calculate CDK MD...")
  smis <- cmpDf$smiles
  dc <- rcdk::get.desc.categories()
  if(category == "all") dn <- unique(unlist(sapply(dc, rcdk::get.desc.names)))
  else dn <- rcdk::get.desc.names(category)
  pb <- utils::txtProgressBar(max = length(smis), style = 3)
  progress <- function(n){utils::setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  cl <- snow::makeCluster(thread)
  doSNOW::registerDoSNOW(cl)
  descsList <- foreach::`%dopar%`(foreach::foreach(smi = smis, n = 1:length(smis),
                                                   .packages = c("rcdk"),
                                                   .options.snow = opts),
                                  {
                                    mol <- rcdk::parse.smiles(smi)
                                    rcdk::eval.desc(mol, dn)
                                  })
  snow::stopCluster(cl)
  gc()
  descs <- purrr::list_rbind(descsList)
  return(dplyr::as_tibble(cbind(cmpDf, descs)))
}
# data("cmpDf_demo", package = "LipRtPred")
# .getRDKitMD(cmpDf = cmpDf_demo)
.getRDKitMD <- function(cmpDf, thread = 1){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  wrong_idx <- .check_smiles(smiles = cmpDf$smiles)
  if(length(wrong_idx) != 0){
    message(paste0("Wrong Idx: ", paste0(wrong_idx, collapse = " ")))
    cmpDf <- cmpDf[-wrong_idx, ]
  }
  message("Calculate RDKit MD...")
  reticulate::source_python(system.file("python", "RDKit_MD.py", package = "LipRtPred"))
  descs <- getRDKitMD_py(smis = as.list(cmpDf$smiles), thread = thread)
  return(dplyr::as_tibble(cbind(cmpDf, descs)))
}
# data("cmpDf_demo", package = "LipRtPred")
# .getMordredMD(cmpDf = cmpDf_demo)
.getMordredMD <- function(cmpDf, ignore_3D = TRUE, thread = 1){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  wrong_idx <- .check_smiles(smiles = cmpDf$smiles)
  if(length(wrong_idx) != 0){
    message(paste0("Wrong Idx: ", paste0(wrong_idx, collapse = " ")))
    cmpDf <- cmpDf[-wrong_idx, ]
  }
  message("Calculate Mordred MD...")
  reticulate::source_python(system.file("python", "Mordred_MD.py", package = "LipRtPred"))
  descs <- getMordredMD_py(smis = as.list(cmpDf$smiles), ignore_3D = ignore_3D, thread = thread)
  return(dplyr::as_tibble(cbind(cmpDf, descs)))
}
#' @rdname getMD
#' @title Calculate molecule descriptors.
#' @description
#' Three molecule descriptors can be computed, CDK, RDKit and Mordred.
#'
#' @inherit GetMorgan_FP
#' @param category CDK molecule descriptors category include "all", "protein", "hybrid", "constitutional", "topological", "electronic" and "geometrical"
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' data("cmpDf_demo", package = "LipRtPred")
#' descsDf <- GetCDK_MD(cmpDf = cmpDf_demo)
GetCDK_MD <- function(cmpDf, flavor = "CxSmiles", category = "all", thread = 1){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  descsDf <- .getCDKMD(cmpDf = cmpDf, category = category, thread = thread)
  return(descsDf)
}
#' @rdname getMD
#' @export
#'
#' @examples
#' descsDf <- GetRDKit_MD(cmpDf = cmpDf_demo)
GetRDKit_MD <- function(cmpDf, flavor = "CxSmiles", thread = 1){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  descsDf <- .getRDKitMD(cmpDf = cmpDf, thread = thread)
  return(descsDf)
}
#' @rdname getMD
#' @param ignore_3D Whether to ignore 3D descriptors.
#' @export
#'
#' @examples
#' descsDf <- GetMordred_MD(cmpDf = cmpDf_demo)
GetMordred_MD <- function(cmpDf, flavor = "CxSmiles", ignore_3D = TRUE, thread = 1){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles, flavor = flavor)
  descsDf <- .getMordredMD(cmpDf = cmpDf, ignore_3D = ignore_3D, thread = as.integer(thread))
  return(descsDf)
}
