#' Fit all multiplicative and additive models
#'
#' Wrapper function to fit all 2K+1 models (multiplicative, additive, and nonperiodic).
#'
#' @param K the number of frequency components in the largest model to fit.
#' @param spikes a list of spike trains.
#' @param f.common.table a table whose names contain the high-amplitude frequency components as computed by \code{\link{FindTopFrequencies}}.
#' @param setup.pars a list of additional parameters for the likelihood function, computed by \code{\link{SetupLikelihoods}}.
#' @param terminal.points a numeric vector containing the time points at which \eqn{c(t)} changes.
#' @param ct a numeric vector containing the estimated piecewise constant intensity function \eqn{c(t)}. The length of \eqn{c(t)} should be a whole number power of 2.
#' @param user.select whether to allow the user to select the frequencies for the model; must be FALSE for this function to run effectively.
#'
#' @return A list of length 3 is returned.
#' The first item in the list is a list of frequency estimates for each model.
#' The second item in the list is a list of phase estimates for each model.
#' The third item in the list is a list of eta/gamma estimates and fit criteria for each model.
#'
#' @export

FitAllMultAddModels <- function(K, spikes, f.common.table, setup.pars, terminal.points, ct, user.select = FALSE) {
## This abstracts Step 6 and 7 of the original script to run for each neuron
## For some value K, it fits Multiplicative 1-K, Additive 1-K, and the no periodicity model (i.e. c(t))


 f.hat.list <- vector("list", length = 2*K+1)  # length 2K+1 list of NULLs
 w0.hat.list <- vector("list", length = 2*K+1)
 K.list <- vector("list", length = 2*K+1)

for(k in 1:K){
   f.hat.list[[2*k]] <- f.hat.list[[2*(k-1)+1]] <- SelectTopFrequencies(f.common.table, k, user.select = user.select)
   w0.hat.list[[2*(k-1)+1]] <-  w0.hat.list[[2*k]] <- EstimatePhase(spikes,f.hat.list[[2*k]])
   K.list[[2*(k-1)+1]] <- FitMultiplicativeModel(spikes, f.hat.list[[2*(k-1)+1]], w0.hat.list[[2*(k-1)+1]], setup.pars, terminal.points, ct)
   K.list[[2*k]] <- FitAdditiveModel(spikes, f.hat.list[[2*k]], w0.hat.list[[2*k]], setup.pars, terminal.points, ct)
   names(f.hat.list)[c(2*k-1, 2*k)] <- paste(c("Multiplicative", "Additive"), k)
}

   f.hat.list[[2*K+1]] <- 0
   w0.hat.list[[2*K+1]] <-  EstimatePhase(spikes,f.hat.list[[2*K+1]])
   K.list[[2*K+1]] <- FitNonperiodicModel(spikes, setup.pars, terminal.points, ct)

   names(f.hat.list)[2*K+1] <- "Nonperiodic"
   names(w0.hat.list) <- names(f.hat.list)
   names(K.list) <- names(f.hat.list)

return(list(f.hat.list = f.hat.list, w0.hat.list = w0.hat.list, K.list = K.list))
}
