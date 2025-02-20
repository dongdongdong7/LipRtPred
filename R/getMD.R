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
# data("cmpDf_demo", package = "LipRtPred")
# category <- c("all", "protein", "hybrid", "constitutional", "topological", "electronic", "geometrical")[1]
# .getCDKMD(cmpDf = cmpDf_demo)
.getCDKMD <- function(cmpDf, category = "all"){
  if(any(is.na(cmpDf$smiles))){
    na_idx <- which(is.na(cmpDf$smiles))
    stop(paste0("cmpDf[", na_idx, ", ] has NA smiles!"))
  }
  message("Calculate CDK MD...")
  mols <- rcdk::parse.smiles(cmpDf$smiles)
  dc <- rcdk::get.desc.categories()
  if(category == "all") dn <- unique(unlist(sapply(dc, rcdk::get.desc.names)))
  else dn <- rcdk::get.desc.names(category)
  pb <- utils::txtProgressBar(max = length(mols), style = 3)
  # The parallel method doesn't work
  # Because th Java-Object jobjRef can't be propagated and the result is null
  descsList <- lapply(1:length(mols), function(i) {
    utils::setTxtProgressBar(pb, i)
    mol <- mols[[i]]
    rcdk::eval.desc(mol, dn)
  })
  descs <- purrr::list_rbind(descsList)
  return(dplyr::as_tibble(cbind(cmpDf, descs)))
}

#' @title GetCDK_MD
#' @description
#' Calculate cdk molecule descriptors for a cmpDf.
#'
#' @param cmpDf A tibble or data.frame with two column, id and smiles.
#' @param flavor SMILES flavor.
#' @param category cdk molecule descriptors category.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' data("cmpDf_demo", package = "LipRtPred")
#' descsDf <- GetCDK_MD(cmpDf = cmpDf_demo)
GetCDK_MD <- function(cmpDf, flavor = "CxSmiles", category = "all"){
  cmpDf$smiles <- .convertSMILES(smiles = cmpDf$smiles)
  descsDf <- .getCDKMD(cmpDf = cmpDf_demo)
  return(descsDf)
}
