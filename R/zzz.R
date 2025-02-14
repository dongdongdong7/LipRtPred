# " _      _       _____  _   _____              _"
# "| |    (_)     |  __ \| | |  __ \            | |"
# "| |     _ _ __ | |__) | |_| |__) | __ ___  __| |"
# "| |    | | '_ \|  _  /| __|  ___/ '__/ _ \/ _` |"
# "| |____| | |_) | | \ \| |_| |   | | |  __/ (_| |"
# "|______|_| .__/|_|  \_\\__|_|   |_|  \___|\__,_|"
# "         | |                                     "
# "         |_|                                     "
# Ascii art Big
.onLoad <- function(libname, pkgname){
  # Install conda
  conda_version <- tryCatch({
    reticulate::conda_version()
  }, error = function(e) {
    NULL
  })
  if(is.null(conda_version)){
    cat("Conda is not found!\n",
        "Please install miniconda first!\n",
        "https://www.anaconda.com/download/success",
        sep = "")
  }else{
    cat("Conda is installed. Version: ", conda_version, "\n")
    # LipRtPred environment
    if(!reticulate::condaenv_exists("LipRtPred")){
      cat("Create LipRtPred environment...")
      reticulate::conda_create(envname = "LipRtPred", python_version = "3.11")
      reticulate::py_install("rdkit", envname = "LipRtPred", pip = TRUE)
    }
    reticulate::use_condaenv(condaenv = "LipRtPred")
    #rdkit <<- reticulate::import("rdkit", delay_load = TRUE)
  }
}
.onAttach <- function(libname, pkgname){
  packageStartupMessage(
    paste("Welcome to LipRtPred!",
          " _      _       _____  _   _____              _",
          "| |    (_)     |  __ \\| | |  __ \\            | |",
          "| |     _ _ __ | |__) | |_| |__) | __ ___  __| |",
          "| |    | | '_ \\|  _  /| __|  ___/ '__/ _ \\/ _` |",
          "| |____| | |_) | | \\ \\| |_| |   | | |  __/ (_| |",
          "|______|_| .__/|_|  \\_\\\\__|_|   |_|  \\___|\\__,_|",
          "         | |                                     ",
          "         |_|                                     ",
          "This is a package to predict lipid retention time.",
          sep = "\n")
  )
}
