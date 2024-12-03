#' This function estimates the parameters of n locations with T observations in 2 regimes. The smoothed probabilities of two regimes are also generated using the estimated parameters. 
#' @importFrom spdep nb2mat cell2nb
#' @importFrom FastGP rcppeigen_invert_matrix rcppeigen_get_det
#' @importFrom Rsolnp solnp
#' @importFrom TSclust diss
#' @importFrom ggplot2 ggplot geom_point geom_line scale_x_date theme aes ylab element_blank element_line
#' @importFrom msm deltamethod
#' @param Y_i  A \eqn{T \times n} matrix consisting \eqn{T} observations of \eqn{n} spatial locations.
#' @param Wn   A \eqn{n \times n} spatial weight matrix. 
#' @return    \item{df.est}{A data frame containing the estimated parameters, associated standard errors, Z-values, and p-values.
#'            See details for the descriptions for each parameter.}
#'            \item{smoothed.prob}{A matrix containing the estimated smoothed probabilities for 2 regimes.}
#'            
#' @details   df.est contains the following parameters: 
#'   1. lambda1/2: Estimated spatial parameter of the 1st/2nd regime.
#'   2. gamma1/2: Estimated temporal parameter of the 1st/2nd regime.
#'   3. rho1/2: Estimated spatio-temporal parameter of the 1st/2nd regime.
#'   4. p/q: Estimated staying probability of the 1st/2nd regime.
#'   5. a1/a2: Estimated homogeneous fixed effect for regime 1 for \eqn{n} spatial locations.
#'   6. s1: Estimated of the variance in the assumed i.i.d. Gaussian variable in the quasi-MLE.
#' @export
#'
#' @examples
#' # Example 1:
#' # Define variables
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
#' Res = MSSTARCH(t(Y[[1]]), Wn)
#' 
#' # Example 2:
#' Y              <- stock_exchanges_28  # the dataset in this package
#' Date_ymd       <- rownames(stock_exchanges_28)
#' 
#' # Calculate the weight matrix using the Piccolo distance
#' # k-nearest neighbors (k=5) along with AR.PIC
#' k              <- 5 
#' dist_mat       <- as.matrix(TSclust::diss(t(log(stock_exchanges_28^2)), "AR.PIC"))
#' Wmat           <- t(sapply(1:ncol(dist_mat), function(i) ifelse(dist_mat[i,] < sort(dist_mat[i,])[k+2] & dist_mat[i,] > 0, 1/k, 0))) 
#' # k+2, because of the diagonal zero entry and the strict inequality
#' diag(Wmat)     <- 0
#' 
#' Res            <- MSSTARCH(stock_exchanges_28, Wmat)
#' 
#' # Plot the smoothed probability of regime 1. 
#' 
#' smoothed       <- Res[[2]]
#' Date           <- as.Date(Date_ymd)
#' smoothed.reg1  <- smoothed[1:dim(stock_exchanges_28)[1],2]
#' df             <- data.frame(Date,smoothed.reg1)
#' ggplot2::ggplot(df, ggplot2::aes(Date,smoothed.reg1)) +       
#'     ggplot2::geom_point(col="black") +    
#'     ggplot2::scale_x_date(date_breaks = "year", date_labels = "%b\n%Y") +
#'     ggplot2::geom_line(col="black") + ggplot2::ylab("Smoothed Probability of Regime 1")+
#'     ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
#'                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))      
MSSTARCH <- function(Y_i, Wn)
{
  T1  <- dim(Y_i)[1]
  n   <- dim(Y_i)[2]
  Y_i <- log(Y_i^2)
  objfun <- function(x) 
  { 
    #Parameters setting
    lamb1   =  -1 + 2/(1+exp(-x[1]))
    gamma1  =  -1 + 2/(1+exp(-x[2]))
    rho1    =  -1 + 2/(1+exp(-x[3]))
    lamb2   =  -1 + 2/(1+exp(-x[4]))
    gamma2  =  -1 + 2/(1+exp(-x[5]))
    rho2    =  -1 + 2/(1+exp(-x[6]))
    p       = 1/(1+exp(-x[7]))
    q       = 1/(1+exp(-x[8]))
    a1      = x[9]
    a2      = x[10]
    s1      = exp(x[11])
    P       = matrix(c(p,1-q,1-p,q),nrow=2,ncol=2,byrow=TRUE)
    pi_inf  = c( (1-q)/(2-p-q) , (1-p)/(2-p-q))
    ksi     = pi_inf
    L       = c(0,0)
    LL      = rep(0,T1-1)
    Sn1     = (diag(n)-lamb1*Wn)
    Sn2     = (diag(n)-lamb2*Wn)
    Det_Sn1 <- FastGP::rcppeigen_get_det(Sn1)
    Det_Sn2 <- FastGP::rcppeigen_get_det(Sn2)
    # Build LL
    for (tt in 2:T1){
      e1     = Sn1%*%Y_i[tt,] - gamma1*Y_i[tt-1,] - rho1*Wn%*%Y_i[tt-1,] - rep(a1,n) + rep(1.27,n)
      e2     = Sn2%*%Y_i[tt,] - gamma2*Y_i[tt-1,] - rho2*Wn%*%Y_i[tt-1,] - rep(a2,n) + rep(1.27,n)
      f_yt_1 = Det_Sn1*((s1*2*pi)^(-n/2))*exp((-0.5/(s1*1))*t(e1)%*%e1)
      f_yt_2 = Det_Sn2*((s1*2*pi)^(-n/2))*exp((-0.5/(s1*1))*t(e2)%*%e2)
      L      = c(ksi[1]*f_yt_1, ksi[2]*f_yt_2)
      xsi    = L/sum(L)
      ksi    = P%*%xsi
      LL[tt] = log(L[1]+L[2])
    }
    return(-sum(LL))
  }
  
  init = c(0.1,0.1,0.1,
           0.4,0.1,0.1,
           0.8,0.8,
           -1,-1, 
           2)
  
  initmap = c( log((1+init[1])/(1-init[1])),
               log((1+init[2])/(1-init[2])),
               log((1+init[3])/(1-init[3])),
               log((1+init[4])/(1-init[4])),
               log((1+init[5])/(1-init[5])),
               log((1+init[6])/(1-init[6])),
               log((init[7])/(1-init[7])),
               log((init[8])/(1-init[8])),
               init[9],
               init[10],
               log(init[11]))
  
  
  confun1 <- function(x){
    
    lamb1  =  -1 + 2/(1+exp(-x[1]))
    gamma1 =  -1 + 2/(1+exp(-x[2]))
    rho1   =  -1 + 2/(1+exp(-x[3]))
    
    lamb2  =  -1 + 2/(1+exp(-x[4]))
    gamma2 =  -1 + 2/(1+exp(-x[5]))
    rho2   =  -1 + 2/(1+exp(-x[6]))
    
    p      = 1/(1+exp(-x[7]))
    q      = 1/(1+exp(-x[8]))
    
    Sn1    = (diag(n)-lamb1*Wn)
    Sn2    = (diag(n)-lamb2*Wn)
    A1     = FastGP::rcppeigen_invert_matrix(Sn1)%*%(gamma1*diag(n) + rho1*Wn)
    A2     = FastGP::rcppeigen_invert_matrix(Sn2)%*%(gamma2*diag(n) + rho2*Wn)
    M1     = p*norm(A1,type="1")^2 + (1-p)*norm(A2,type="1")^2
    M2     = (1-q)*norm(A1,type="1")^2 +  q*norm(A2,type="1")^2
    
    return(c(M1,M2))
  }
  
  Res         <- Rsolnp::solnp( pars = initmap, fun = objfun, ineqfun = confun1, ineqLB = c(0,0), ineqUB = c(0.9,0.9))
  # Calculate transformed parameters
  x           <- Res$pars
  lamb1  =  -1 + 2/(1+exp(-x[1]))
  gamma1 =  -1 + 2/(1+exp(-x[2]))
  rho1   =  -1 + 2/(1+exp(-x[3]))
  lamb2  =  -1 + 2/(1+exp(-x[4]))
  gamma2 =  -1 + 2/(1+exp(-x[5]))
  rho2   =  -1 + 2/(1+exp(-x[6]))
  p      = 1/(1+exp(-x[7]))
  q      = 1/(1+exp(-x[8]))
  a1     = x[9]
  a2     = x[10]
  s1     = exp(x[11])
  Parameters  <- c(lamb1 = lamb1, gamma1 = gamma1,
                   rho1 = rho1, lamb2= lamb2,
                   gamma2 = gamma2, rho2 = rho2,
                   p = p, q = q,
                   a1 = a1, a2 = a2, s1 = s1)
  Std.Errors  <- msm::deltamethod(list(~-1 + 2/(1+exp(-x1)), 
                                       ~-1 + 2/(1+exp(-x2)), 
                                       ~-1 + 2/(1+exp(-x3)),
                                       ~-1 + 2/(1+exp(-x4)),
                                       ~-1 + 2/(1+exp(-x5)),
                                       ~-1 + 2/(1+exp(-x6)),
                                       ~ 1/(1+exp(-x7)),
                                       ~ 1/(1+exp(-x8)),
                                       ~ x9,
                                       ~ x10,
                                       ~ exp(x11)), mean=Res$pars, cov = solve(Res$hessian[-c(1,2),-c(1,2)]))
  z_values    <- Parameters/Std.Errors
  Pr          <- exp(-0.717*z_values-0.416*z_values^2)
  # Construct output data frame for the estimated parameters
  df.est      <- data.frame(Parameters, Std.Errors, z_values, Pr)
  # 2nd step: Computing smoothed probabilities
  P       = matrix(c(p,1-q,1-p,q),nrow=2,ncol=2,byrow=TRUE)
  pi_inf  = c( (1-q)/(2-p-q) , (1-p)/(2-p-q))
  ksi     = matrix(rep(0,2*T1+2),nrow=2,ncol=T1+1)
  ksi[,2] = pi_inf  
  xsi     = matrix(rep(0,2*T1),nrow=2,ncol=T1)
  L       = c(0,0)
  LL      = rep(0,T1-1)
  Sn1     = (diag(n)-lamb1*Wn)
  Sn2     = (diag(n)-lamb2*Wn)
  Det_Sn1 = (rcppeigen_get_det(Sn1))
  Det_Sn2 = (rcppeigen_get_det(Sn2))
  
  for (t in 2:T1){
    e1        = Sn1%*%Y_i[t,] - gamma1*Y_i[t-1,] - rho1*Wn%*%Y_i[t-1,] - rep(a1,n) + rep(1.27,n)
    e2        = Sn2%*%Y_i[t,] - gamma2*Y_i[t-1,] - rho2*Wn%*%Y_i[t-1,] - rep(a2,n) + rep(1.27,n)
    f_yt_1    = Det_Sn1*((s1*2*pi)^(-n/2))*exp((-0.5/s1)*t(e1)%*%e1)
    f_yt_2    = Det_Sn2*((s1*2*pi)^(-n/2))*exp((-0.5/s1)*t(e2)%*%e2)
    L         = c(ksi[1,t]*f_yt_1, ksi[2,t]*f_yt_2)
    LL[t]     = log(L[1]+L[2])
    xsi[,t]   = L/sum(L)
    ksi[,t+1] = P%*%xsi[,t]
  }
  
  eps_prob = matrix(rep(0,2*T1),nrow=2,ncol=T1)
  eps_prob[,T1] = xsi[,T1]
  for (t in 1:(T1-1)){
    eps_prob[,T1-t] = xsi[,T1-t]*(t(P)%*%(eps_prob[,T1-t+1]/ksi[,T1-t+1]))
  }
  
  smoothed.prob <- t(eps_prob)
  return(list(df.est, smoothed.prob))
}






