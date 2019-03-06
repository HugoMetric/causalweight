#' Causal mediation analysis with instruments for treatment and mediator based on weighting
#' @description Causal mediation analysis (evaluation of natural direct and indirect effects) with instruments for a binary treatment and a continuous mediator based on weighting as suggested in Frölich and Huber (2017),  Theorem 1.
#' @param y Dependent variable, must not contain missings.
#' @param d Treatment, must be binary (either 1 or 0), must not contain missings.
#' @param m Mediator(s),must be a continuous scalar, must not contain missings.
#' @param zd Instrument for the treatment,  must be binary (either 1 or 0), must not contain missings.
#' @param zm Instrument for the mediator, must contain at least one continuous element, may be a scalar or a vector, must not contain missings. If no user-specified bandwidth is provided for the regressors when estimating the conditional cumulative distribution function F(M|Z2,X), i.e. if \code{bwreg=NULL}, then \code{zm} must be exclusively numeric.
#' @param x Pre-treatment confounders, may be a scalar or a vector, must not contain missings. If no user-specified bandwidth is provided for the regressors when estimating the conditional cumulative distribution function F(M|Z2,X), i.e. if \code{bwreg=NULL}, then \code{x} must be exclusively numeric.
#' @param trim Trimming rule for discarding observations with extreme weights. Discards observations whose relative weight would exceed the value in \code{trim} in the estimation of any of the potential outcomes. Default is 0.1 (i.e. a maximum weight of 10\% per observation).
#' @param csquared If TRUE, then not only the control function C, but also its square is used as regressor in any estimated function that conditions on C. Default is FALSE.
#' @param boot Number of bootstrap replications for estimating standard errors. Default is 1999.
#' @param cminobs Minimum number of observations to compute the control function C, see the numerator of equation (7) in Frölich and Huber (2017). A larger value increases boundary bias when estimating the control function for lower values of M, but reduces the variance. Default is 40, but should be adapted to sample size and the number of variables in Z2 and X.
#' @param bwreg Bandwidths for \code{zm} and \code{x} in the estimation of the conditional cumulative distribution function F(M|Z2,X) based on the np package by Hayfield and Racine (2008). The length of the numeric vector must correspond to the joint number of elements in \code{zm} and \code{x} and will be used both in the original sample for effect estimation and in bootstrap samples to compute standard errors. If set to \code{NULL}, then the rule of thumb is used for bandwidth calculation, see the np package for details. In the latter case, all elements in the regressors must be numeric. Default is \code{NULL}.
#' @param bwm Bandwidth for \code{m} in the estimation of the conditional cumulative distribution function F(M|Z2,X) based on the np package by Hayfield and Racine (2008). Must be scalar and will be used both in the original sample for effect estimation and in bootstrap samples to compute standard errors. If set to \code{NULL}, then the rule of thumb is used for bandwidth calculation, see the np package for details. Default is \code{NULL}.
#' @param logit If FALSE, probit regression is used for any propensity score estimation. If TRUE, logit regression is used. Default is FALSE.
#' @param cluster A cluster ID for block or cluster bootstrapping when units are clustered rather than iid. Must be numerical. Default is NULL (standard bootstrap without clustering).
#' @details Estimation of causal mechanisms (natural direct and indirect effects) of a binary treatment among treatment compliers based on distinct instruments for the treatment and the mediator. The treatment and its instrument are assumed to be binary, while the mediator and its instrument are assumed to be continuous, see Theorem 1 in Frölich and Huber (2017). The instruments are assumed to be conditionally valid given a set of observed confounders. A control function is used to tackle mediator endogeneity. Standard errors are obtained by bootstrapping the effects.
#' @return A medlateweight object contains two components, \code{results} and \code{ntrimmed}:
#' @return \code{results}: a 3x7 matrix containing the effect estimates in the first row ("effects"), standard errors in the second row ("se"), and p-values in the third row ("p-value").
#' The first column provides the total effect, namely the local average treatment effect (LATE) on the compliers.
#' The second and third columns provide the direct effects under treatment and control, respectively ("dir.treat", "dir.control").
#' The fourth and fifth columns provide the indirect effects under treatment and control, respectively ("indir.treat", "indir.control").
#' The sixth and seventh columns provide the parametric direct and indirect effect estimates ("dir.para", "indir.para") without intercation terms, respectively. For the parametric estimates, probit or logit specifications are used for the treatment model and OLS specifications for the mediator and outcome models.
#' @return \code{ntrimmed}: number of discarded (trimmed) observations due to large weights.
#' @references Frölich, M. and Huber, M. (2017): "Direct and indirect treatment effects: Causal chains and mediation analysis with instrumental variables", Journal of the Royal Statistical Society Series B, 79, 1645–1666.
#' @examples # A little example with simulated data (3000 observations)
#' \dontrun{
#' n=3000; sigma=matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),3,3)
#' e=(rmvnorm(n,rep(0,3),sigma))
#' x=rnorm(n)
#' zd=(0.5*x+rnorm(n)>0)*1
#' d=(-1+0.5*x+2*zd+e[,3]>0)
#' zm=0.5*x+rnorm(n)
#' m=(0.5*x+2*zm+0.5*d+e[,2])
#' y=0.5*x+d+m+e[,1]
#' # The true direct and indirect effects on compliers are equal to 1 and 0.5, respectively
#' medlateweight(y,d,m,zd,zm,x,trim=0.1,csquared=FALSE,boot=19,cminobs=40,
#'               bwreg=NULL,bwm=NULL,logit=FALSE)}
#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm quantile
#' @importFrom np npcdensbw npcdist
#' @import mvtnorm
#' @export
medlateweight=function(y,d,m,zd, zm, x, trim=0.1,  csquared=FALSE, boot=1999, cminobs=40, bwreg=NULL, bwm=NULL, logit=FALSE, cluster=NULL){
  temp<-effects.late.x(y=y,d=d,m=m,zd=zd, zm=zm, x=x, trim=trim, csquared=csquared, bwreg=bwreg, bwm=bwm, cminobs=cminobs, logit=logit)
  temp2<-bootstrap.mediation.late.x(y=y,d=d,m=m,zd=zd,zm=zm, x=x, boot=boot,trim=trim, csquared=csquared, bwreg=bwreg, bwm=bwm, cminobs=cminobs, logit=logit, cluster=cluster)
  ntrimmed=temp[length(temp)];
  temp=temp[1:(length(temp)-1)]
  se=apply(temp2[,1:(ncol(temp2)-1)], 2, sd)
  temp3=2*pnorm(-abs(temp/se))
  results=rbind(temp, se, temp3)
  colnames(results)=c("LATE", "dir.treat", "dir.control", "indir.treat", "indir.control", "dir.para", "indir.para")
  rownames(results)=c("effect", "se", "p-value")
  list(results=results, ntrimmed=ntrimmed)
}
