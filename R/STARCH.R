#' QMLE estimation of spatio-temporal data with 1 regime.
#' @importFrom bbmle mle2
#' @importFrom msm deltamethod
#' @param Y_i  A \eqn{T \times n} matrix consisting \eqn{T} observations of \eqn{n} spatial locations.
#' @param Wn   A \eqn{n \times n} spatial weight matrix. 
#'
#' @return    \item{df.est}{A data frame containing the estimated parameters, associated standard errors, z-values, and p-values.
#'            See details for the descriptions for each parameter.}
#'            
#' @details   df.est contains the following parameters: 
#'   1. lambda: Estimated spatial parameter.
#'   2. gamma: Estimated temporal parameter.
#'   3. rho: Estimated spatio-temporal parameter.
#'   5. a1: Estimated homogeneous fixed effect for \eqn{n} spatial locations.
#'   6. s1: Estimated of the variance in the assumed i.i.d. Gaussian variable in the quasi-MLE.
#'               
#' @export
#'
#' @examples
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
#' Res            <- STARCH(stock_exchanges_28, Wmat)
STARCH <- function(Y_i,Wn)
{
  T1  = dim(Y_i)[1]
  n   = dim(Y_i)[2]
  Y_i = log(Y_i^2)
  LL  = rep(0,T1-1)
  objfun <- function(x1,x2,x3,x4,x5) 
  { 
    lamb1    = (exp(x1) - 1)/(exp(x1)+exp(x3)+1)
    gamma1   = -1 + 2/(1+exp(-x2))
    rho1     = (exp(x3) - 1)/(exp(x1)+exp(x3)+1)
    a1       = x4
    s1       = exp(x5)
    Sn1      = (diag(n)-lamb1*Wn)
    Det_Sn1  = (rcppeigen_get_det(Sn1))
    for (tt in 2:T1){
      e1     = Sn1%*%Y_i[tt,] - gamma1*Y_i[tt-1,] - rho1*Wn%*%Y_i[tt-1,] - rep(a1,n) + rep(1.27,n)
      f_yt_1 = Det_Sn1*((s1*2*pi)^(-n/2))*exp((-0.5/s1)*t(e1)%*%e1)
      LL[tt] = log(f_yt_1)
    }
    return(-sum(LL))
  }
  
  
  init = c(0.2,0.2,0.2,-2,4)
  initmap = list( 
    x1 = log((init[1]+1)/(1-init[1]) + init[1]*(1-init[1]+2*init[3])/(1-init[1]-init[3])), 
    x2 = log((1+init[2])/(1-init[2])),
    x3 = log((1-init[1]+2*init[3])/(1-init[1]-init[3])),
    x4  = init[4],
    x5  = log(init[5])
  )
  
  Res         <- bbmle::mle2( start = initmap , minuslogl = objfun, method="BFGS")
  x           <- Res@coef
  lamb1       <- (exp(x[1]) - 1)/(exp(x[1])+exp(x[3])+1)
  gamma1      <- -1 + 2/(1+exp(-x[2]))
  rho1        <- (exp(x[3]) - 1)/(exp(x[1])+exp(x[3])+1)
  a1          <- x[4]
  s1          <- exp(x[5])
  Parameters  <- c(lamb1 = lamb1,gamma1 = gamma1,rho1 = rho1,a1 = a1,s1 = s1)
  Std.Errors  <- msm::deltamethod(list(~ (exp(x1) - 1)/(exp(x1)+exp(x3)+1), 
                                       ~ -1 + 2/(1+exp(-x2)), 
                                       ~ (exp(x3) - 1)/(exp(x1)+exp(x3)+1),
                                       ~ x4,
                                       ~ exp(x5)), mean=Res@coef, cov = Res@vcov)
  z_values      <- Parameters/Std.Errors
  Pr            <- exp(-0.717*z_values-0.416*z_values^2)
  df.est        <- data.frame(Parameters, Std.Errors, z_values, Pr)
  return(df.est)
}
