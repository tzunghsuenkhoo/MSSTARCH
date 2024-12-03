#' This functions simulate the data of n locations in a period of T1 which have a Markovian structure of 2 regimes. 
#' @importFrom spdep nb2mat cell2nb
#' @importFrom markovchain markovchainSequence
#' @importFrom methods new
#' @importFrom FastGP rcppeigen_invert_matrix
#' @importFrom stats rnorm
#' @param n  An integer. The number of spatial locations 
#' @param T1 An integer. The length of the period considered.
#' @param param1 A vector of length 4. param1 should contains lambda1, gamma1, rho1, alpha1. Details is described in the function "MSSTARCH".
#' @param param2 A vector of length 4. param2 should contains lambda2, gamma2, rho2, alpha2. Details is described in the function "MSSTARCH".
#' @param sigma2 A positive real number. The variance of the mean-zero Gaussian error variable.
#' @param Wn  A \eqn{n \times n} matrix. Spatial weight matrix.
#' @param trans_probs A vector of length 2. This vector contains the staying probability of the 2 states.
#'
#' @return \item{Y}{A \eqn{n \times T1} matrix. This matrix consists of simulated data of n locations in a period of T1.}
#'         \item{out}{A vector of length T1. This vector consists of realized states at each t in 1,2,...,T1.}
#' 
#' @export
#' 
#' @examples
#' # Define variables
#' 
#' n      = 36
#' T1     = 200
#' lamb1  = 0.2
#' gamma1 = 0.2
#' rho1   = 0.2
#' alpha1 = 0.1
#' lamb2  = 0.2
#' gamma2 = 0.8
#' rho2   = -0.2
#' alpha2 = 0.1
#' sigma  = 1
#' trans_probs = c(0.97,0.93)
#' 
#' # Simulate a n times n queen-contiguity weight matrix
#' 
#' Wn = spdep::nb2mat(spdep::cell2nb(sqrt(n), sqrt(n), type = "queen"))
#' 
#' # Generate the simulated data.
#' Y  =  sim_tworegdat(n,T1,Wn,c(lamb1,gamma1,rho1,alpha1), c(lamb2,gamma2,rho2,alpha2), sigma, trans_probs) 
sim_tworegdat <- function(n, T1, Wn, param1, param2, sigma2, trans_probs){
  # Assigns parameters
  lambW1 <- param1[1]*Wn
  gamma1 <- param1[2]
  rhoW1  <- param1[3]*Wn 
  alpha1 <- param1[4]
  
  lambW2 <- param2[1]*Wn
  gamma2 <- param2[2]
  rhoW2  <- param2[3]*Wn
  alpha2 <- param1[4]
  
  p <- trans_probs[1]
  q <- trans_probs[2]
  
  # To create 1st observation Y_00 (Burn-in period)
  eps        <- rnorm(n,mean=0,sd=(sigma2)^0.5)
  mu_i0_1    <- rep(alpha1,n)
  mu_i0_2    <- rep(alpha2,n)
  lagged_var <- mu_i0_1
  tau        <- FastGP::rcppeigen_invert_matrix(diag(n) - lambW1)%*%(lagged_var+log(eps^2))
  X          <- lambW1%*%tau + lagged_var 
  Y_00       <- diag(as.vector(exp(X)))^(0.5)%*%eps
  
  for(b in 1:20){
    eps        <- rnorm(n,mean=0,sd=(sigma2)^0.5)
    lagged_var <- mu_i0_1  + gamma1*log(Y_00^2) + rhoW1%*%log(Y_00^2)
    tau        <- FastGP::rcppeigen_invert_matrix(diag(n) - lambW1)%*%(lagged_var+log(eps^2))
    X          <- lambW1%*%tau + lagged_var 
    Y          <- diag(as.vector(exp(X)))^(0.5)%*%eps
    Y_00       <- Y
  }
  
  # Create actual simulated Y_i. First pre-generates the states at each time t in 1,2,...,T1
  t2          <- T1
  Y           <- matrix(,nrow = n, ncol = t2)
  statesNames <- c("l", "h")
  mcB         <- methods::new("markovchain", states = statesNames, 
                     transitionMatrix = matrix(c(p,1-p,1-q,q), 
                                       nrow = 2, byrow = TRUE, dimnames = list(statesNames, statesNames)))
  # Show the sequence
  outs <- markovchain::markovchainSequence(n = T1, markovchain = mcB, t0 = "l")
  for (b in 1:t2){
    if (outs[b] == "l"){
      eps        <- rnorm(n,mean=0,sd=(sigma2)^0.5)
      lagged_var <- mu_i0_1 + gamma1*log(Y_00^2) + rhoW1%*%log(Y_00^2)
      tau        <- FastGP::rcppeigen_invert_matrix(diag(n) - lambW1)%*%(lagged_var+log(eps^2))
      X          <- lambW1%*%tau + lagged_var
      Y[,b]      <- diag(as.vector(exp(X)))^(0.5)%*%eps
      Y_00       <- Y[,b]           
    }
    if (outs[b] == "h"){
      eps        <- rnorm(n,mean=0,sd=(sigma2)^0.5)
      lagged_var <- mu_i0_2 + gamma2*log(Y_00^2) + rhoW2%*%log(Y_00^2)
      tau        <- FastGP::rcppeigen_invert_matrix(diag(n) - lambW2)%*%(lagged_var+log(eps^2))
      X          <- lambW2%*%tau + lagged_var
      Y[,b]      <- diag(as.vector(exp(X)))^(0.5)%*%eps
      Y_00       <- Y[,b]          
    } 
    
  }
  return(list(Y,outs))
}