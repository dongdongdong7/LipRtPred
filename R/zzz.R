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
  if(!reticulate::virtualenv_exists(envname = "LipRtPred")){
    message("You are building python virtual environment LipRtPred...")
    reticulate::install_python(version = "3.11")
    reticulate::virtualenv_create(
      envname = "LipRtPred",
      version = "3.11",
    )
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
