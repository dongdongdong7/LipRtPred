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

#' @title Draw a molecule using a SMILES
#'
#' @param smi A SMILES string
#' @param width Width of pic
#' @param height Height of pic
#' @param addAtomIndices Whether to add atom indices
#' @param addBondIndices Whether to add bond indices
#' @param addStereoAnnotation Whether to add stereoisomorphic annotations
#' @param explicitMethyl Whether to display the terminal methyl group
#'
#' @return NULL
#' @export
#'
#' @examples
#' draw_smi(smi = "C(OCCCCCCCC(CC)(CCC)CCCCC(C)C)C(OC(CCCCCCCCCCCC(C)C)=O)COC(=O)CCCCCCC(OCCC)CCCCC(C)C",
#'          width = 1600, height = 800)
draw_smi <- function(smi, width = 450, height = 150,
                     addAtomIndices = TRUE,
                     addBondIndices = FALSE,
                     addStereoAnnotation = FALSE,
                     explicitMethyl = FALSE){
  reticulate::source_python(system.file("python", "draw_smis.py", package = "LipRtPred"))
  tempfile_name <- paste0(tempfile(), ".png")
  draw_smi_py(smi = smi,
              filePath = tempfile_name,
              width = as.integer(width), height = as.integer(height),
              addAtomIndices = addAtomIndices,
              addBondIndices = addBondIndices,
              addStereoAnnotation = addStereoAnnotation,
              explicitMethyl = explicitMethyl)
  img <- png::readPNG(tempfile_name)
  grid::grid.raster(img, width = grid::unit(1, "npc"))
  file.remove(tempfile_name)
  invisible(NULL)
}
