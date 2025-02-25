# Generate parallel environment
# Barry Song
# 250225

#' @title Generate parallel environment
#'
#' @param thread Number of parallel threads
#'
#' @return cl
#' @export
#'
#' @examples
#' parallelEnv()
parallelEnv <- function(thread = parallel::detectCores()){
  cl <- parallel::makeCluster(thread)
  doParallel::registerDoParallel(cores = cl)
  return(cl)
}
