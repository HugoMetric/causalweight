#' Treatment evaluation based on inverse probability weighting with optional sample selection correction.
#' @description Treatment evaluation based on inverse probability weighting with optional sample selection correction.
#' @param y Dependent variable.
#' @param d Treatment, must be binary (either 1 or 0), must not contain missings.
#' @param x Confounders of the treatment and outcome, must not contain missings.
#' @param s Selection indicator. Must be one if \code{y} is observed (non-missing) and zero if \code{y} is not observed (missing). Default is \code{NULL}, implying that \code{y} does not contain any missings.
#' @param z Optional instrumental variable(s) for selection \code{s}. If \code{NULL}, outcome selection based on observables (\code{x},\code{d}) - known as "missing at random" - is assumed.
#' If \code{z} is defined, outcome selection based on unobservables - known as "non-ignorable missingness" - is assumed. Default is \code{NULL}. If \code{s} is \code{NULL}, \code{z} is ignored.
#' @param selpop Only to be used if both \code{s} and \code{z} are defined. If \code{TRUE}, the effect is estimated for the selected subpopulation with \code{s}=1 only. If \code{FALSE}, the effect is estimated for the total population.
#' (note that this relies on somewhat stronger statistical assumptions). Default is \code{FALSE}. If \code{s} or  \code{z} is \code{NULL}, \code{selpop} is ignored.
#' @param ATET If \code{FALSE}, the average treatment effect (ATE) is estimated. If \code{TRUE}, the average treatment effect on the treated (ATET)  is estimated. Default is \code{FALSE}.
#' @param trim Trimming rule for discarding observations with extreme propensity scores. If \code{ATET=FALSE}, observations with Pr(D=1|X)<\code{trim} or Pr(D=1|X)>(1-\code{trim}) are dropped.
#' If \code{ATET=TRUE}, observations with Pr(D=1|X)>(1-\code{trim}) are dropped. If \code{s} is defined and \code{z} is \code{NULL}, observations with extremely low selection propensity scores, Pr(S=1|D,X)<\code{trim}, are discarded, too. If \code{s} and \code{z} are defined, the treatment propensity scores to be trimmed change to Pr(D=1|X,Pr(S=1|D,X,Z)). If in addition \code{selpop} is \code{FALSE}, observation with Pr(S=1|D,X,Z)<\code{trim} are discarded, too. Default for \code{trim} is 0.05.
#' @param logit If \code{FALSE}, probit regression is used for propensity score estimation. If \code{TRUE}, logit regression is used. Default is \code{FALSE}.
#' @param boot Number of bootstrap replications for estimating standard errors. Default is 1999.
#' @param cluster A cluster ID for block or cluster bootstrapping when units are clustered rather than iid. Must be numerical. Default is NULL (standard bootstrap without clustering).
#' @details Estimation of treatment effects of a binary treatment under a selection on observables assumption assuming that all confounders of the treatment and the outcome are observed. Units are weighted by the inverse of their conditional treatment propensities given the observed confounders, which are estimated by probit or logit regression. Standard errors are obtained by bootstrapping the effect.
#' If \code{s} is defined, the procedure allows correcting for sample selectiondue to missing outcomes based on the inverse of the conditional selection probability. The latter might either be related to observables, which implies a missing at random assumption, or in addition also to unobservables, if an instrument for sample selection is available. See Huber (2012, 2014) for further details.
#' @return A treatweight object contains six components: \code{effect}, \code{se}, \code{pval}, \code{y1}, \code{y0}, and \code{ntrimmed}.
#' @return \code{effect}: average treatment effect (ATE) if \code{ATET=FALSE} or the average treatment effect on the treated (ATET) if \code{ATET=TRUE}.
#' @return \code{se}: bootstrap-based standard error of the effect.
#' @return \code{pval}: p-value of the effect.
#' @return \code{y1}: mean potential outcome under treatment.
#' @return \code{y0}: mean potential outcome under control.
#' @return \code{ntrimmed}: number of discarded (trimmed) observations due to extreme propensity score values.
#' @references Horvitz, D. G., and Thompson, D. J. (1952): "A generalization of sampling without replacement from a finite universe", Journal of the American Statistical Association, 47, 663–685.
#' @references Huber, M. (2012): "Identification of average treatment effects in social experiments under alternative forms of attrition", Journal of Educational and Behavioral Statistics, 37 , 443-474.
#' @references Huber, M. (2014): "Treatment evaluation in the presence of sample selection", Econometric Reviews, 33, 869-905.
#' @examples # A little example with simulated data (10000 observations)
#' n=10000
#' x=rnorm(n); d=(0.25*x+rnorm(n)>0)*1
#' y=0.5*d+0.25*x+rnorm(n)
#' # The true ATE is equal to 0.5
#' output=treatweight(y=y,d=d,x=x, trim=0.05, ATET=FALSE, logit=TRUE, boot=19)
#' cat("ATE: ",round(c(output$effect),3),", standard error: ",
#'     round(c(output$se),3), ", p-value: ",round(c(output$pval),3))
#' output$ntrimmed
#' @examples # An example with non-random outcome selection and an instrument for selection
#' n=10000
#' sigma=matrix(c(1,0.6,0.6,1),2,2)
#' e=(2*rmvnorm(n,rep(0,2),sigma))
#' x=rnorm(n)
#' d=(0.5*x+rnorm(n)>0)*1
#' z=rnorm(n)
#' s=(0.25*x+0.25*d+0.5*z+e[,1]>0)*1
#' y=d+x+e[,2]; y[s==0]=0
#' # The true ATE is equal to 1
#' output=treatweight(y=y,d=d,x=x, s=s, z=z, selpop=FALSE, trim=0.05, ATET=FALSE, logit=TRUE, boot=19)
#' cat("ATE: ",round(c(output$effect),3),", standard error: ",
#'     round(c(output$se),3), ", p-value: ",round(c(output$pval),3))
#' output$ntrimmed
#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm quantile
#' @import mvtnorm
#' @export
treatweight<-function(y,d,x, s=NULL, z=NULL, selpop=FALSE, ATET=FALSE, trim=0.05, logit=FALSE, boot=1999, cluster=NULL){
  if (is.null(s)==FALSE) y[s==0]=0
  temp=ipw(y=y,d=d,x=x, s=s, z=z, selpop=selpop, trim=trim, ATET=ATET, logit=logit)
  temp2=bootstrap.ipw(y=y,d=d,x=x,s=s, z=z, selpop=selpop, boot=boot,trim=trim, ATET=ATET, logit=logit, cluster=cluster)
  se=sd(temp2[,1])
  list(effect=temp[1], se=se, pval=2*pnorm(-abs(temp[1]/se)), y1=temp[2], y0=temp[3], ntrimmed=temp[4])
}



