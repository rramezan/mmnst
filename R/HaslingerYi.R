#' Calculate y(i)
#'
#' Calculates Yi in equation (39) of Haslinger \emph{et al.} (2010)
#'
#' @importFrom stats runif
#'
#' @param spike.train a list of numeric vectors, each of which contains the spike times for one trial of the experiment.
#' @param endpoints.bins a numeric vector containing the endpoints of the bins for the binned spike train.
#' @param q a numeric vector containing the q values calculated by the \code{\link{HaslingerQ}} function.
#'
#' @return A numeric vector of Yi values as defined in equation (39) of Haslinger \emph{et al.} (2010).
#'
#' @references Haslinger, R., Pipa G., and Brown, E. (2010). Discrete time rescaling theorem: determining goodness of fit for discrete time statistical models of neural spiking. \emph{Neural Computation}. \strong{22}(10):2477-506. doi: 10.1162/NECO_a_00015.


HaslingerYi<-function(spike.train,endpoints.bins,q){
##spike.train is the spike train
##endpoints.bins is the vector of start/endpoints of the bins

N <- length(spike.train)
xi <- numeric(N-1)

for(i in 1:(N-1)){
	##fix to get around "no non-missing arguments to min/max" bug
	indx1 <- min(length(endpoints.bins),which(endpoints.bins>spike.train[i]))
	indx2 <- max(1,which(endpoints.bins<spike.train[i+1]))-1
	sum.term <- (indx2>indx1)*sum(q[indx1:indx2])

	r <- runif(1)
	xi[i] <- sum.term-log(1-r*(1-exp(-q[(indx2+1)])))
}##end for loop

y <- 1-exp(-xi)

return(y)
}
