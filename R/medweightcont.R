#' Causal mediation analysis with a continuous treatment based on weighting by the inverse of generalized propensity scores
#' @description Causal mediation analysis (evaluation of natural direct and indirect effects) of a continuous treatment based on weighting by the inverse of generalized propensity scores as suggested in Hsu, Huber, Lee, and Pipoz (2018).
#' @param y Dependent variable, must not contain missings.
#' @param d Continuous treatment, must not contain missings.
#' @param m Mediator(s), may be a scalar or a vector, must not contain missings.
#' @param x Pre-treatment confounders of the treatment, mediator, and/or outcome, must not contain missings.
#' @param d0 Value of \code{d} under non-treatment. Effects are based on pairwise comparisons, i.e. differences in potential outcomes evaluated at \code{d1} and \code{d0}.
#' @param d1 Value of \code{d} under treatment. Effects are based on pairwise comparisons, i.e. differences in potential outcomes evaluated at \code{d1} and \code{d0}.
#' @param ATET If FALSE, the average treatment effect (ATE) and the corresponding direct and indirect effects are estimated. If TRUE, the average treatment effect on the treated (ATET)  and the corresponding direct and indirect effects are estimated. Default is FALSE.
#' @param trim Trimming rule for discarding observations with extreme generalized propensity scores. If \code{lognorm=FALSE}, observations with f(D=\code{d1}|M,X)<\code{trim} or f(D=\code{d0}|M,X)<\code{trim} are dropped, with f denoting the generalized propensity score (or conditional density of treatment). If \code{lognorm=TRUE}, then \code{trim} corresponds to the share of lowest f(D=\code{d1}|M,X) or f(D=\code{d0}|M,X), respectively, that are dropped.
#' @param lognorm If FALSE, a linear model with normally distributed errors is assumed for generalized propensity score estimation. If TRUE, a lognormal model is assumed. Default is FALSE.
#' @param bw Bandwith for the second order Epanechnikov kernel functions of the treatment. If set to NULL, the rule of thumb for Epanechnikov kernels is used for bandwidth computation. Default is NULL.
#' @param boot Number of bootstrap replications for estimating standard errors. Default is 1999.
#' @param cluster A cluster ID for block or cluster bootstrapping when units are clustered rather than iid. Must be numerical. Default is NULL (standard bootstrap without clustering).
#' @details Estimation of causal mechanisms (natural direct and indirect effects) of a continuous treatment under a selection on observables assumption assuming that all confounders of the treatment and the mediator, the treatment and the outcome, or the mediator and the outcome are observed. Units are weighted by the inverse of their conditional treatment densities (known as generalized propensity scores) given the mediator and/or observed confounders, which are estimated by linear or loglinear regression.
#' Standard errors are obtained by bootstrapping the effects.
#' @return A medweightcont object contains two components, \code{results} and \code{ntrimmed}:
#' @return \code{results}: a 3X5 matrix containing the effect estimates in the first row ("effects"), standard errors in the second row ("se"), and p-values in the third row ("p-value").
#' The first column provides the total effect, namely the average treatment effect (ATE) if \code{ATET=FALSE} or the average treatment effect on the treated (ATET), i.e. those with D=\code{d1}, if \code{ATET=TRUE}.
#' The second and third columns provide the direct effects under treatment and control, respectively ("dir.treat", "dir.control"). The fourth and fifth columns provide the indirect effects under treatment and control, respectively ("indir.treat", "indir.control").
#' @return \code{ntrimmed}: number of discarded (trimmed) observations due to extreme propensity score values.
#' @references Hsu, Y.-C., Huber, M., Lee, Y.-Y., Pipoz, L.: (2018): "Direct and indirect effects of continuous treatments based on generalized propensity score weighting",  SES working paper 495, University of Fribourg.
#' @examples # A little example with simulated data (10000 observations)
#' n=10000
#' x=runif(n=n,min=-1,max=1)
#' d=0.25*x+runif(n=n,min=-2,max=2)
#' d=d-min(d)
#' m=0.5*d+0.25*x+runif(n=n,min=-2,max=2)
#' y=0.5*d+m+0.25*x+runif(n=n,min=-2,max=2)
#' # The true direct and indirect effects are all equal to 0.5
#' output=medweightcont(y,d,m,x, d0=2, d1=3, ATET=FALSE, trim=0.05, lognorm=FALSE, bw=NULL, boot=19)
#' round(output$results,3)
#' output$ntrimmed
#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm dnorm quantile
#' @importFrom np npksum
#' @export

medweightcont<-function(y,d,m,x, d0, d1, ATET=FALSE, trim=0.05, lognorm=FALSE, bw=NULL, boot=1999, cluster=NULL){
  if(is.null(bw)) bw=sd(d)*2.34/(length(d)^0.2)
  temp=mediation.cont(y=y,d=d,m=m,x=x, d0=d0, d1=d1,  ATET=ATET, trim=trim, lognorm=lognorm, bw=bw)
  ntrimmed=temp[length(temp)]
  temp=temp[1:(length(temp)-1)]
  temp2=bootstrap.mediation.cont(y=y,d=d,m=m,x=x, d0=d0, d1=d1,  ATET=ATET, trim=trim, lognorm=lognorm, bw=bw, boot=boot, cluster=cluster)
  se=apply(temp2[,1:(ncol(temp2)-1)], 2, sd)
  temp3=2*pnorm(-abs(temp/se))
  results=rbind(temp, se, temp3)
  if (ATET==FALSE) colnames(results)=c("ATE", "dir.treat", "dir.control", "indir.treat", "indir.control")
  if (ATET==TRUE) colnames(results)=c("ATET", "dir.treat", "dir.control", "indir.treat", "indir.control")
  rownames(results)=c("effect", "se", "p-value")
  list(results=results, ntrimmed=ntrimmed)
}
