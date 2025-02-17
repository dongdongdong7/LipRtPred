#' @title Draw molecules using smiles
#'
#' @param smiles A vector with smiles.
#' @param SMARTS SMARTS pattern.
#' @param molsPerRow mols number per row.
#' @param subImgSize Size of sub image. The higher.
#' The larger the molecule the larger the subImgSize should be assigned to make the image clearer
#'
#' @return None.
#' @export
#'
#' @examples
#' draw_smiles(smiles = c("CC(O)CC(O)C", "C(O)CC(O)C"), SMARTS = "CO")
draw_smiles <- function(smiles, SMARTS = NULL, molsPerRow = 1L, subImgSize = c(400L, 400L)){
  reticulate::source_python(system.file("python", "draw_smis.py", package = "LipRtPred"))
  tempfile_name <- paste0(tempfile(), ".png")
  draw_smis(smis = as.list(smiles),
            filePath = tempfile_name,
            SMARTS = SMARTS,
            molsPerRow = as.integer(molsPerRow),
            subImgSize = as.integer(subImgSize))
  img <- png::readPNG(tempfile_name)
  grid::grid.raster(img)
  file.remove(tempfile_name)
  invisible(NULL)
}
