#' Generate log normal variables from mean and sd inputs
#'
#' @param nsim number of simulations
#' @param X mean value from user input for single contaminant
#' @param SD SD value from user input for single contaminant
#' 
#' @export
#' 
#' @details http://yasai.rutgers.edu/yasai-guide-27.html, but see macros code for YASAI_sccwrp
genlognorm_fun <- function(nsim, X, SD){
  
  # genlognormal, see link for doc
  mu <- log(X) - 0.5 * log(1 + SD ^ 2 / X ^ 2)
  sigma <- (log(SD ^ 2 / X ^ 2 + 1)) ^ 0.5
  sims <- suppressWarnings({exp(rnorm(nsim, mu, sigma))})
  
  simi <- seq(1:nsim)
  out <- data.frame(i = simi, sims = sims)
  return(out)
  
}