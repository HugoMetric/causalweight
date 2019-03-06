### R code from vignette source 'bodory-huber.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: bodory-huber.Rnw:273-281 (eval = FALSE)
###################################################
## n=10000
## set.seed(100); x=rnorm(n)
## set.seed(101); d=(0.25*x+rnorm(n)>0)*1
## set.seed(102); y=0.5*d+0.25*x+rnorm(n)
## output=treatweight(y=y,d=d,x=x,trim=0.05,ATET=FALSE,logit=TRUE, boot=19)
## cat("ATE: ",round(c(output$effect),3),", standard error: ",
##     round(c(output$se),3), ", p-value: ",round(c(output$pval),3))
## output$ntrimmed


###################################################
### code chunk number 2: bodory-huber.Rnw:299-312 (eval = FALSE)
###################################################
## n=10000
## sigma=matrix(c(1,0.6,0.6,1),2,2)
## set.seed(100); e=(2*rmvnorm(n,rep(0,2),sigma))
## set.seed(101); x=rnorm(n)
## set.seed(102); d=(0.5*x+rnorm(n)>0)*1
## set.seed(103); z=rnorm(n)
## s=(0.25*x+0.25*d+0.5*z+e[,1]>0)*1
## y=d+x+e[,2]; y[s==0]=0
## output=treatweight(y=y,d=d,x=x, s=s,z=z,selpop=FALSE,trim=0.05,ATET=FALSE,
##                    logit=TRUE,boot=19)
## cat("ATE: ",round(c(output$effect),3),", standard error: ",
##     round(c(output$se),3), ", p-value: ",round(c(output$pval),3))
## output$ntrimmed


###################################################
### code chunk number 3: bodory-huber.Rnw:484-494 (eval = FALSE)
###################################################
## n=10000
## set.seed(100); x=rnorm(n);
## set.seed(101); d=(0.25*x+rnorm(n)>0)*1
## set.seed(102); w=0.2*d+0.25*x+rnorm(n);
## set.seed(103); m=0.5*w+0.5*d+0.25*x+rnorm(n)
## set.seed(104); y=0.5*d+m+w+0.25*x+rnorm(n)
## output=medweight(y=y,d=d,m=m,x=x,w=w,trim=0.05,ATET=FALSE,logit=TRUE,
## boot=19)
## round(output$results,3)
## output$ntrimmed


###################################################
### code chunk number 4: bodory-huber.Rnw:633-644 (eval = FALSE)
###################################################
## n=10000
## set.seed(100); u=rnorm(n)
## set.seed(101); x=rnorm(n)
## set.seed(102); z=(0.25*x+rnorm(n)>0)*1
## set.seed(103); d=(z+0.25*x+0.25*u+rnorm(n)>0.5)*1
## y=0.5*d+0.25*x+u
## output=lateweight(y=y,d=d,z=z,x=x,trim=0.05,LATT=FALSE,logit=TRUE,boot=19)
## cat("LATE: ",round(c(output$effect),3),", standard error: ",
##              round(c(output$se.effect),3), ", p-value: ",
##              round(c(output$pval.effect),3))
## output$ntrimmed


###################################################
### code chunk number 5: bodory-huber.Rnw:800-812 (eval = FALSE)
###################################################
## n=3000
## sigma=matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),3,3)
## set.seed(100); e=(rmvnorm(n,rep(0,3),sigma))
## set.seed(101); x=rnorm(n)
## set.seed(102); zd=(0.5*x+rnorm(n)>0)*1
## d=(-1+0.5*x+2*zd+e[,3]>0)
## set.seed(103); zm=0.5*x+rnorm(n)
## m=(0.5*x+2*zm+0.5*d+e[,2])
## y=0.5*x+d+m+e[,1]
## options(digits=3)
## medlateweight(y,d,m,zd,zm,x,trim=0.1,csquared=FALSE,boot=19,cminobs=40,
##               bwreg=NULL,bwm=NULL,logit=FALSE)


