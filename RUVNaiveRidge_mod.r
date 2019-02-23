#' Removal of unwanted variation for gene correlations.
#' 
#' \code{RUVNaiveRidge} applies the ridged version of global removal of unwanted variation 
#' to simulated or real gene expression data. 
#'
#' @param Y A matrix of gene expression values or an object of 
#' class \code{simulateGEdata}. 
#' @param center A logical scalar; if \code{TRUE} the data is centered, 
#' if \code{FALSE} data is assumed to be already centered.
#' @param nc_index A vector of indices of negative controls.
#' @param nu A numeric scalar value of \code{nu} \eqn{\geq 0}.
#' @param kW An integer setting the number of dimensions for the estimated noise.
#' @param check.input A logical scalar; if \code{TRUE} all input is 
#' checked (not advisable for large simulations).
#' @return \code{RUVNaiveRidge} returns a list with a matrix of the cleaned 
#' (RUV-treated) centered gene expression values called Yhat and a matrix
#' systematic noise Walphahat.
#' @details 
#' The parameter \code{kW} controls how much noise is cleaned, whereas the 
#' parameter \code{nu} controls the amount of ridging to deal with possible dependence of 
#' the noise and the factor of interest.
RUVNaiveRidge<-function(
      Y, 
      center=TRUE, ##set equal to FALSE in case of centered data
      nc_index, ##column index for negative controls 
      nu, ## Ridge factor
      kW, ## number of noise dimensions
      check.input=FALSE) UseMethod("RUVNaiveRidge")

#' \code{RUVNaiveRidge.default} applies the ridged version of global 
#' removal of unwanted variation to matrices.
#'
#' @rdname RUVNaiveRidge
#' @export
RUVNaiveRidge.default<-function(
      Y, ##matrix of gene expression data 
      center=TRUE, ##set equal to FALSE in case of centered data
      nc_index, ##column index for negative controls 
      nu, ## Ridge factor
      kW, ## number of noise dimensions
      check.input=FALSE)
{ 
  
  if(check.input){
    if(is.matrix(Y)==FALSE){stop("Y needs to be a matrix.")}
    if(nu<0){stop("nu has to be positive or 0.")}
    if(kW>dim(Y)[1]){ stop("kW is too big.") }
   }
  
  if(center){
    Y<-scale(Y, center=TRUE, scale=FALSE)
    ## center data
  }
  
  Yc<-Y[, nc_index]
  ## subset negative controls
  
  tmp<-svd(Yc, nu=kW, nv=kW)
  S.d<-diag(tmp$d[1:kW], nrow=kW, ncol=kW)
  ## SVD of negative controls
  
  W.hat<-tmp$u%*%S.d
  ## estimate W.hat
  
  alpha.hat<-solve(t(W.hat)%*%W.hat+nu*diag(dim(W.hat)[2]))%*%t(W.hat)%*%Y
  ## estimate alpha.hat
  
  return(list(Yhat=Y-W.hat%*%alpha.hat, Walphahat=W.hat%*%alpha.hat, What=data.frame(W.hat))) 
  ## calculate Y.hat (note that this is the mean centered Y)
}  