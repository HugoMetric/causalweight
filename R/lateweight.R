#' Local average treatment effect estimation based on inverse probability weighting
#' @description Instrumental variable-based evaluation of local average treatment effects using weighting by the inverse of the instrument propensity score.
#' @param y Dependent variable, must not contain missings.
#' @param d Treatment, must be binary (either 1 or 0), must not contain missings.
#' @param z Instrument for the endogenous treatment, must be binary (either 1 or 0), must not contain missings.
#' @param x Confounders of the instrument and outcome, must not contain missings.
#' @param LATT If FALSE, the local average treatment effect (LATE) among compliers (whose treatment reacts to the instrument) is estimated. If TRUE, the local average treatment effect on the treated compliers (LATT)  is estimated. Default is FALSE.
#' @param trim Trimming rule for discarding observations with extreme propensity scores. If \code{LATT=FALSE}, observations with Pr(Z=1|X)<\code{trim} or Pr(Z=1|X)>(1-\code{trim}) are dropped.
#' If \code{LATT=TRUE}, observations with Pr(Z=1|X)>(1-\code{trim}) are dropped. Default is 0.05.
#' @param cluster A cluster ID for block or cluster bootstrapping when units are clustered rather than iid. Must be numerical. Default is NULL (standard bootstrap without clustering).
#' @param logit If FALSE, probit regression is used for propensity score estimation. If TRUE, logit regression is used. Default is FALSE.
#' @param boot Number of bootstrap replications for estimating standard errors. Default is 1999.
#' @details Estimation of local average treatment effects of a binary endogenous treatment based on a binary instrument that is conditionally valid, implying that all confounders of the instrument and the outcome are observed. Units are weighted by the inverse of their conditional instrument propensities given the observed confounders, which are estimated by probit or logit regression. Standard errors are obtained by bootstrapping the effect.
#' @return A lateweight object contains 10 components, \code{effect}, \code{se.effect}, \code{pval.effect}, \code{first}, \code{se.first}, \code{pval.first}, \code{ITT}, \code{se.ITT}, \code{pval.ITT}, and \code{ntrimmed}:
#' @return \code{effect}: local average treatment effect (LATE) among compliers if \code{LATT=FALSE} or the local average treatment effect on treated compliers (LATT) if \code{LATT=TRUE}.
#' @return \code{se.effect}: bootstrap-based standard error of the effect.
#' @return \code{pval.effect}: p-value of the effect.
#' @return \code{first}: first stage estimate of the complier share if \code{LATT=FALSE} or the first stage estimate among treated if \code{LATT=TRUE}.
#' @return \code{se.first}: bootstrap-based standard error of the first stage effect.
#' @return \code{pval.first}: p-value of the first stage effect.
#' @return \code{ITT}: intention to treat effect (ITT) of \code{z} on \code{y} if \code{LATT=FALSE} or the ITT among treated if \code{LATT=TRUE}.
#' @return \code{se.ITT}: bootstrap-based standard error of the ITT.
#' @return \code{pval.ITT}: p-value of the ITT.
#' @return \code{ntrimmed}: number of discarded (trimmed) observations due to extreme propensity score values.
#' @references Frölich, M. (2007): "Nonparametric IV estimation of local average treatment effects with covariates", Journal of Econometrics, 139, 35-75.
#' @examples # A little example with simulated data (10000 observations)
#' n=10000
#' u=rnorm(n)
#' x=rnorm(n)
#' z=(0.25*x+rnorm(n)>0)*1
#' d=(z+0.25*x+0.25*u+rnorm(n)>0.5)*1
#' y=0.5*d+0.25*x+u
#' # The true LATE is equal to 0.5
#' output=lateweight(y=y,d=d,z=z, x=x, trim=0.05, LATT=FALSE, logit=TRUE, boot=19)
#' cat("LATE: ",round(c(output$effect),3),", standard error: ",
#'              round(c(output$se.effect),3), ", p-value: ",
#'              round(c(output$pval.effect),3))
#' output$ntrimmed
#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm quantile
#' @import mvtnorm
#' @export
lateweight<-function(y,d,z,x, LATT=FALSE, trim=0.05, logit=FALSE, boot=1999, cluster=NULL){
  temp=late(y=y,d=d,z=z, x=x,trim=trim, LATT=LATT, logit=logit)
  ntrimmed=temp[length(temp)]
  temp=temp[1:(length(temp)-1)]
  temp2=bootstrap.late(y=y,d=d,z=z, x=x,boot=boot,trim=trim, LATT=LATT, logit=logit, cluster=cluster)
  se=apply(temp2[,1:(ncol(temp2)-1)], 2, sd)
  temp3=2*pnorm(-abs(temp/se))
  list(effect=temp[1], se.effect=se[1], pval.effect=temp3[1],  first=temp[2], se.first=se[2], pval.first=temp3[2], ITT=temp[3], se.ITT=se[3], pval.ITT=temp3[3], ntrimmed=ntrimmed)
}



