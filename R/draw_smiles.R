#' @title Draw molecules using smiles
#' @param smiles A vector with smiles.
#'
#' @return None.
#' @export
#'
#' @examples
#' draw_smiles(smiles = c("c1ccccc1O", "CCCCC"))
draw_smiles <- function(smiles){
  reticulate::source_python(system.file("python", "draw_smis.py", package = "LipRtPred"))
  tempfile_name <- paste0(tempfile(), ".png")
  draw_smis(smis = as.list(smiles),
            filePath = tempfile_name)
  img <- png::readPNG(tempfile_name)
  grid::grid.raster(img)
}
