mediation<-function(y,d,m,x,w=NULL,s=NULL,z=NULL, selpop=FALSE, trim=0.05, ATET=FALSE, logit=FALSE){
  if (is.null(w)==TRUE){
    if (logit==FALSE){
      if (is.null(z)==TRUE){
        pscore.mx=glm(d~cbind(m,x),family=binomial(probit))$fitted
        pscore.x=glm(d~x,family=binomial(probit))$fitted
      }
      if (is.null(s)==FALSE & is.null(z)==TRUE) pscore.s=glm(s~cbind(d,m,x),family=binomial(probit))$fitted
      if (is.null(s)==FALSE & is.null(z)==FALSE) {
        pscore.s=glm(s~cbind(d,m,x,z),family=binomial(probit))$fitted
        pscore.mx=glm(d~cbind(m,x,pscore.s),family=binomial(probit))$fitted
        pscore.x=glm(d~cbind(x,pscore.s),family=binomial(probit))$fitted
      }
    }
    if (logit==TRUE){
      if (is.null(z)==TRUE){
        pscore.mx=glm(d~cbind(m,x),family=binomial(logit))$fitted
        pscore.x=glm(d~x,family=binomial(logit))$fitted
      }
      if (is.null(s)==FALSE & is.null(z)==TRUE) pscore.s=glm(s~cbind(d,m,x),family=binomial(logit))$fitted
      if (is.null(s)==FALSE & is.null(z)==FALSE) {
        pscore.s=glm(s~cbind(d,m,x,z),family=binomial(logit))$fitted
        pscore.mx=glm(d~cbind(m,x,pscore.s),family=binomial(logit))$fitted
        pscore.x=glm(d~cbind(x,pscore.s),family=binomial(logit))$fitted
      }
    }
    if (is.null(s)==TRUE | selpop==TRUE){
      ind=((pscore.mx<trim) | (pscore.mx>(1-trim)) )
      y=y[ind==0]; d=d[ind==0]; pscore.mx=pscore.mx[ind==0]; pscore.x=pscore.x[ind==0]
    }
    if (is.null(s)==FALSE & selpop==FALSE) {
      ind=((pscore.mx<trim) | (pscore.mx>(1-trim)) | pscore.s<trim)
      y=y[ind==0]; d=d[ind==0]; s=s[ind==0]; pscore.mx=pscore.mx[ind==0]; pscore.x=pscore.x[ind==0]; pscore.s=pscore.s[ind==0]
    }
    if ((is.null(s)==TRUE) | (is.null(s)==FALSE & selpop==TRUE)){
      if (is.null(s)==FALSE & selpop==TRUE) {y=y[s==1]; d=d[s==1]; pscore.mx=pscore.mx[s==1]; pscore.x=pscore.x[s==1]}
      if (ATET==FALSE){
        y1m1=sum(y*d/pscore.x)/sum(d/pscore.x)
        y1m0=sum(y*d*(1-pscore.mx)/(pscore.mx*(1-pscore.x)))/sum(d*(1-pscore.mx)/(pscore.mx*(1-pscore.x)))
        y0m0=sum(y*(1-d)/(1-pscore.x))/sum((1-d)/(1-pscore.x))
        y0m1=sum(y*(1-d)*pscore.mx/((1-pscore.mx)*pscore.x))/sum((1-d)*pscore.mx/((1-pscore.mx)*pscore.x))
      }
      if (ATET==TRUE){
        y1m1=sum(y*d)/sum(d)
        y1m0=sum(y*d*(1-pscore.mx)*pscore.x/(pscore.mx*(1-pscore.x)))/sum(d*(1-pscore.mx)*pscore.x/(pscore.mx*(1-pscore.x)))
        y0m0=sum(y*(1-d)*pscore.x/(1-pscore.x))/sum((1-d)*pscore.x/(1-pscore.x))
        y0m1=sum(y*(1-d)*pscore.mx/((1-pscore.mx)))/sum((1-d)*pscore.mx/((1-pscore.mx)))
      }
    }

   if (is.null(s)==FALSE & selpop==FALSE){
      if (ATET==FALSE){
        y1m1=sum(y*d*s/(pscore.s*pscore.x))/sum(d*s/(pscore.s*pscore.x))
        y1m0=sum(y*d*s*(1-pscore.mx)/(pscore.s*pscore.mx*(1-pscore.x)))/sum(d*s*(1-pscore.mx)/(pscore.s*pscore.mx*(1-pscore.x)))
        y0m0=sum(y*(1-d)*s/(pscore.s*(1-pscore.x)))/sum((1-d)*s/(pscore.s*(1-pscore.x)))
        y0m1=sum(y*(1-d)*s*pscore.mx/(pscore.s*(1-pscore.mx)*pscore.x))/sum((1-d)*s*pscore.mx/(pscore.s*(1-pscore.mx)*pscore.x))
      }
      if (ATET==TRUE){
        y1m1=sum(y*d*s/pscore.s)/sum(d*s/pscore.s)
        y1m0=sum(y*d*s*(1-pscore.mx)*pscore.x/(pscore.s*pscore.mx*(1-pscore.x)))/sum(d*s*(1-pscore.mx)*pscore.x/(pscore.s*pscore.mx*(1-pscore.x)))
        y0m0=sum(y*(1-d)*s*pscore.x/(pscore.s*(1-pscore.x)))/sum((1-d)*s*pscore.x/(pscore.s*(1-pscore.x)))
        y0m1=sum(y*(1-d)*s*pscore.mx/(pscore.s*(1-pscore.mx)))/sum((1-d)*s*pscore.mx/(pscore.s*(1-pscore.mx)))
      }
    }
    results=c(y1m1-y0m0, y1m1 - y0m1, y1m0-y0m0, y1m1 - y1m0, y0m1 - y0m0, sum(ind))
  }

  if (is.null(w)==FALSE){
    if (logit==FALSE){
      pscore.mwx=glm(d~cbind(m,w,x),family=binomial(probit))$fitted
      pscore.x=glm(d~x,family=binomial(probit))$fitted
      pscore.wx=glm(d~cbind(w,x),family=binomial(probit))$fitted
      pscore.mx=glm(d~cbind(m,x),family=binomial(probit))$fitted
    }
    if (logit==TRUE){
      pscore.mwx=glm(d~cbind(m,w,x),family=binomial(logit))$fitted
      pscore.x=glm(d~x,family=binomial(logit))$fitted
      pscore.wx=glm(d~cbind(w,x),family=binomial(logit))$fitted
      pscore.mx=glm(d~cbind(m,x),family=binomial(logit))$fitted
    }
    ind=((pscore.mwx<trim) | (pscore.mwx>(1-trim)) )
    y=y[ind==0]; d=d[ind==0]; pscore.mx=pscore.mx[ind==0]; pscore.x=pscore.x[ind==0]; pscore.mwx=pscore.mwx[ind==0]; pscore.wx=pscore.wx[ind==0]

    if (ATET==FALSE){
      y1m1=sum(y*d/pscore.x)/sum(d/pscore.x)
      y1m0<-(sum(y*d*(1-pscore.mwx)/((1-pscore.x)*pscore.mwx))/sum(d*(1-pscore.mwx)/((1-pscore.x)*pscore.mwx)))
      y0m0=sum(y*(1-d)/(1-pscore.x))/sum((1-d)/(1-pscore.x))
      y0m1<-(sum(y*(1-d)* pscore.mwx/(pscore.x*(1-pscore.mwx)))/sum((1-d)* pscore.mwx/(pscore.x*(1-pscore.mwx))))
      y1m0p=(sum( y*d/pscore.mwx * (1-pscore.mwx)/(1-pscore.wx)* (pscore.wx)/pscore.x )/sum(d/pscore.mwx * (1-pscore.mwx)/(1-pscore.wx)* (pscore.wx)/pscore.x ))
      y0m1p=sum( y*(1-d)/(1-pscore.mwx) * (pscore.mwx)/(pscore.wx)* (1-pscore.wx)/(1-pscore.x) )/sum((1-d)/(1-pscore.mwx) * (pscore.mwx)/(pscore.wx)* (1-pscore.wx)/(1-pscore.x) )
    }
    if (ATET==TRUE){
      y1m1=sum(y*d)/sum(d)
      y1m0<-(sum(y*d*pscore.x*(1-pscore.mwx)/((1-pscore.x)*pscore.mwx))/sum(d*pscore.x*(1-pscore.mwx)/((1-pscore.x)*pscore.mwx)))
      y0m0=(sum(y*(1-d)*pscore.x/(1-pscore.x))/sum((1-d)*pscore.x/(1-pscore.x)))
      y0m1<-(sum(y*(1-d)* pscore.mwx/((1-pscore.mwx)))/sum((1-d)* pscore.mwx/((1-pscore.mwx))))
      y1m0p=(sum( y*d*pscore.x/pscore.mwx * (1-pscore.mwx)/(1-pscore.wx)* (pscore.wx)/pscore.x )/sum(d*pscore.x/pscore.mwx * (1-pscore.mwx)/(1-pscore.wx)* (pscore.wx)/pscore.x ))
      y0m1p=sum( y*(1-d)*pscore.x/(1-pscore.mwx) * (pscore.mwx)/(pscore.wx)* (1-pscore.wx)/(1-pscore.x) )/sum((1-d)*pscore.x/(1-pscore.mwx) * (pscore.mwx)/(pscore.wx)* (1-pscore.wx)/(1-pscore.x) )
    }
    results=c(y1m1-y0m0, y1m1 - y0m1, y1m0-y0m0, y1m1 - y1m0p, y0m1p - y0m0, sum(ind))

  }
  results
}

bootstrap.mediation<-function(y,d,m,x,w=NULL,s=NULL,z=NULL,boot=1999, selpop=FALSE, trim=0.05, ATET=FALSE, logit=FALSE, cluster=NULL){
  if (is.null(cluster)){
    obs<-length(y)
    bsamples=matrix( ,boot,6)
    for(i in 1:boot){
      sboot<-sample(1:obs,obs,TRUE)
      yb=y[sboot]
      db=d[sboot]
      if (is.null(s)==FALSE) sb=s[sboot]
      if (is.null(s)==TRUE) sb=NULL
      if (is.null(ncol(m))) mb<-m[sboot]
      if (is.null(ncol(m))==0) mb<-m[sboot,]
      if (is.null(ncol(x))) xb<-x[sboot]
      if (is.null(ncol(x))==0) xb<-x[sboot,]
      if ( (is.null(w)==FALSE) & (length(w)==length(y))) wb=w[sboot]
      if ( (is.null(w)==FALSE) & (length(w)!=length(y))) wb=w[sboot,]
      if (is.null(w)==TRUE) wb=NULL
      if ( (is.null(z)==FALSE) & (length(z)==length(y))) zb=z[sboot]
      if ( (is.null(z)==FALSE) & (length(z)!=length(y))) zb=z[sboot,]
      if (is.null(z)==TRUE) zb=NULL

      bsamples[i,]=c(mediation(y=yb,d=db,m=mb,x=xb,w=wb, s=sb, z=zb, selpop=selpop, trim=trim, ATET=ATET, logit=logit))
    }
  }
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; mb=c(); wb=c(); sb=c(); zb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]])
        if (is.null(s)==FALSE) sb<-c(sb,s[key==sboot[k]])
        if (is.null(ncol(m))) mb<-c(mb,m[key==sboot[k]])
        if (is.null(ncol(m))==0) mb=rbind(mb,m[key==sboot[k],])
        if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
        if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
        if ((is.null(w)==FALSE) & is.null(ncol(w))) wb<-c(wb,w[key==sboot[k]])
        if ((is.null(w)==FALSE) & is.null(ncol(w))==0) wb=rbind(wb,w[key==sboot[k],])
        if ((is.null(z)==FALSE) & is.null(ncol(z))) zb<-c(zb,z[key==sboot[k]])
        if ((is.null(z)==FALSE) & is.null(ncol(z))==0) zb=rbind(zb,z[key==sboot[k],])
      }
      if (is.null(w)==TRUE) wb=NULL
      if (is.null(s)==TRUE) sb=NULL
      if (is.null(z)==TRUE) zb=NULL
      est=c(mediation(y=yb,d=db,m=mb,x=xb,w=wb, s=sb, z=zb, trim=trim, ATET=ATET, logit=logit))
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0) cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  bsamples
}


ipw<-function(y,d,x,s,z, selpop=FALSE, trim=0.05, ATET=FALSE, logit=FALSE){
   if (logit==FALSE) {
      if ( (is.null(s))  | (is.null(s)==0 & is.null(z)) )  pscore.x=glm(d~x,family=binomial(probit))$fitted
      if (is.null(s)==0 & is.null(z)) selscore=glm(s~cbind(d,x),family=binomial(probit))$fitted
      if (is.null(s)==0 & is.null(z)==0) {
        selscore=glm(s~cbind(d,x,z),family=binomial(probit))$fitted
        pscore.x=glm(d~cbind(x,selscore),family=binomial(probit))$fitted
      }
    }
    if (logit==TRUE)  {
      if ( (is.null(s))  | (is.null(s)==0 & is.null(z)) )  pscore.x=glm(d~x,family=binomial(logit))$fitted
      if (is.null(s)==0 & is.null(z)) selscore=glm(s~cbind(d,x),family=binomial(logit))$fitted
      if (is.null(s)==0 & is.null(z)==0) {
        selscore=glm(s~cbind(d,x,z),family=binomial(logit))$fitted
        pscore.x=glm(d~cbind(x,selscore),family=binomial(logit))$fitted
      }
    }
    if (ATET==FALSE){
      if (is.null(s)){
        ind=((pscore.x<trim) | (pscore.x>(1-trim)) )
        y=y[ind==0]; d=d[ind==0];  pscore.x=pscore.x[ind==0]
        y1=sum(y*d/pscore.x)/sum(d/pscore.x)
        y0=(sum(y*(1-d)/(1-pscore.x))/sum((1-d)/(1-pscore.x)))
      }
      if ((is.null(s)==0 & is.null(z)) | (is.null(s)==0 & is.null(z)==0 & selpop==FALSE))  {
        ind=((pscore.x<trim) | (pscore.x>(1-trim)) | selscore<trim )
        y=y[ind==0]; d=d[ind==0]; s=s[ind==0];  pscore.x=pscore.x[ind==0]; selscore=selscore[ind==0]
        y1=sum(y*d*s/(pscore.x*selscore))/sum(d*s/(pscore.x*selscore))
        y0=(sum(y*(1-d)*s/((1-pscore.x)*selscore))/sum((1-d)*s/((1-pscore.x)*selscore)))
      }
      if  (is.null(s)==0 & is.null(z)==0 & selpop==TRUE)  {
        ind=((pscore.x<trim) | (pscore.x>(1-trim)) )
        y=y[ind==0 & s==1]; d=d[ind==0 & s==1];  pscore.x=pscore.x[ind==0 & s==1]
        y1=sum(y*d/pscore.x)/sum(d/pscore.x)
        y0=(sum(y*(1-d)/(1-pscore.x))/sum((1-d)/(1-pscore.x)))
      }
    }
    if (ATET==TRUE){
      if (is.null(s)){
        ind= (pscore.x>(1-trim))
        y=y[ind==0]; d=d[ind==0];  pscore.x=pscore.x[ind==0]
        y1=sum(y*d)/sum(d)
        y0=(sum(y*(1-d)*pscore.x/(1-pscore.x))/sum((1-d)*pscore.x/(1-pscore.x)))
      }
      if ((is.null(s)==0 & is.null(z)) | (is.null(s)==0 & is.null(z)==0 & selpop==FALSE))  {
        ind= (pscore.x>(1-trim) | selscore<trim )
        y=y[ind==0]; d=d[ind==0]; s=s[ind==0]; pscore.x=pscore.x[ind==0]; selscore=selscore[ind==0]
        y1=sum(y*d*s/selscore)/sum(d*s/selscore)
        y0=(sum(y*(1-d)*s*pscore.x/((1-pscore.x)*selscore))/sum((1-d)*s*pscore.x/((1-pscore.x)*selscore)))
      }
      if  (is.null(s)==0 & is.null(z)==0 & selpop==TRUE)  {
        ind= (pscore.x>(1-trim))
        y=y[ind==0 & s==1]; d=d[ind==0 & s==1];  pscore.x=pscore.x[ind==0 & s==1]
        y1=sum(y*d)/sum(d)
        y0=(sum(y*(1-d)*pscore.x/(1-pscore.x))/sum((1-d)*pscore.x/(1-pscore.x)))
      }
    }
  results=c(y1-y0, y1, y0, sum(ind))
  results
}

bootstrap.ipw<-function(y,d,x,s=NULL,z=NULL, selpop=FALSE, boot=1999,trim=0.05, ATET=FALSE, logit=FALSE, cluster=NULL){
  if (is.null(cluster)){
   obs<-length(y)
    bsamples=matrix( ,boot,4)
    for(i in 1:boot){
      sboot<-sample(1:obs,obs,TRUE)
      yb=y[sboot]; db<-d[sboot]
      if (is.null(s)==0) sb=s[sboot]
      if (is.null(s)) sb=NULL
      if (is.null(ncol(x))) xb<-x[sboot]
      if (is.null(ncol(x))==0) xb<-x[sboot,]
      if (is.null(z)==0){
        if (is.null(ncol(z))) zb<-z[sboot]
        if (is.null(ncol(z))==0) zb<-z[sboot,]
      }
      if (is.null(z)) zb=NULL
      bsamples[i,]=c(ipw(y=yb,d=db,x=xb, s=sb, z=zb, selpop=selpop, trim=trim, ATET=ATET, logit=logit))
    }
  }
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; sb=c(); zb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]])
        if (is.null(s)==0) sb<-c(sb,s[key==sboot[k]])
        if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
        if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
        if ((is.null(z)==FALSE) & is.null(ncol(z))) zb<-c(zb,z[key==sboot[k]])
        if ((is.null(z)==FALSE) & is.null(ncol(z))==0) zb=rbind(zb,z[key==sboot[k],])
      }
      if (is.null(s)) sb=NULL
      if (is.null(z)) zb=NULL
      est=c(ipw(y=yb,d=db,x=xb, s=sb, z=zb, selpop=selpop, trim=trim, ATET=ATET, logit=logit))
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0){
    cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  }
  bsamples
}



late<-function(y,d,z, x,trim=0.05, LATT=FALSE, logit=FALSE){
  if (logit==FALSE) pscore.x=glm(z~x,family=binomial(probit))$fitted
  if (logit==TRUE)  pscore.x=glm(z~x,family=binomial(logit))$fitted
  if (LATT==FALSE){
    ind=((pscore.x<trim) | (pscore.x>(1-trim)) )
    y=y[ind==0]; d=d[ind==0]; z=z[ind==0]; pscore.x=pscore.x[ind==0]
    firststage=sum(d*z/pscore.x)/(sum(z/pscore.x))-sum(d*(1-z)/(1-pscore.x))/(sum((1-z)/(1-pscore.x)))
    ITT=sum(y*z/pscore.x)/(sum(z/pscore.x))-sum(y*(1-z)/(1-pscore.x))/(sum((1-z)/(1-pscore.x)))
  }
  if (LATT==TRUE){
    ind= (pscore.x>(1-trim))
    y=y[ind==0]; d=d[ind==0]; z=z[ind==0]; pscore.x=pscore.x[ind==0]
    firststage=sum(d*z)/(sum(z))-sum(d*(1-z)*pscore.x/(1-pscore.x))/(sum((1-z)*pscore.x/(1-pscore.x)))
    ITT=sum(y*z)/(sum(z))-sum(y*(1-z)*pscore.x/(1-pscore.x))/(sum((1-z)*pscore.x/(1-pscore.x)))
  }
  results=c(ITT/firststage,  firststage, ITT, sum(ind))
  results
}


bootstrap.late<-function(y,d,z,x,boot=1999,trim=0.05, LATT=FALSE, logit=FALSE, cluster=NULL){
  if (is.null(cluster)){
    obs<-length(y)
    bsamples=matrix( ,boot,4)
    for(i in 1:boot){
      sboot=sample(1:obs,obs,TRUE)
      yb=y[sboot]; db=d[sboot]; zb=z[sboot];
      if (is.null(ncol(x))) xb<-x[sboot]
      if (is.null(ncol(x))==0) xb<-x[sboot,]
      bsamples[i,]=c(late(y=yb,d=db, z=zb, x=xb, trim=trim, LATT=LATT, logit=logit))
    }
  }
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; zb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]]); zb<-c(zb,z[key==sboot[k]])
        if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
        if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
      }
      est=c(late(y=yb,d=db, z=zb, x=xb, trim=trim, LATT=LATT, logit=logit))
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0){
    cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  }
  bsamples
}


effects.late.x<-function(y,d,m,zd,  x, zm, trim=0.05, csquared=FALSE, bwreg=bwreg, bwm=bwm, cminobs=40, logit=FALSE){
  if (is.null(bwreg) | is.null(bwm)) temp<-npcdensbw(ydat=m, xdat=data.frame(zm,x), ckertype="gaussian", bwmethod="normal-reference")
  if (is.null(bwreg))  bwreg<-temp$xbw
  if (is.null(bwm))  bwm<-temp$ybw
  m.dist=npcdist(bws=c(bwreg,bwm), tydat=m, txdat=data.frame(zm,x), ckertype="gaussian")$condist
  x<-as.matrix(cbind(x))
  zm<-as.matrix(cbind(zm))
  nobs=length(y)
  dzd=d*zd
  one_dzd=(1-d)*(1-zd)
  if (logit==FALSE) {
    pscore2=glm(zd~x,family=binomial(probit))$fitted
    pred.d<-fitted.values(glm(d~zd+x,family=binomial(probit)))
  }
  if (logit==TRUE) {
    pscore2=glm(zd~x,family=binomial(logit))$fitted
    pred.d<-fitted.values(glm(d~zd+x,family=binomial(logit)))
  }
  pred.m<-fitted.values(lm(m~cbind(zm,pred.d,x)))
  templm<-lm(y~cbind(pred.d,pred.m,x))
  templm2<-lm(m~cbind(pred.d,x))$coef[2]

  c.d=d*(zd-pscore2)
  c.d_1=(d-1)*(zd-pscore2)
  c.den.d=lm(c.d~zm+x)$fitted

    mm<-sort(m)
     c.num<-c()
    for (i in 1:nrow(x)){
      ind= ((m<=m[i]) | (m<=mm[cminobs]))
      c.num.d=lm(c.d[ind==1]~zm[ind==1,]+x[ind==1,])
      c.num.d_1=lm(c.d_1[ind==1]~zm[ind==1,]+x[ind==1,])
      c.num<-c(c.num, (d[i]*(c(1,zm[i,],x[i,])%*%c.num.d$coef)+(1-d[i])*(c(1,zm[i,],x[i,])%*%c.num.d_1$coef)))
    }
  c=c.num/c.den.d*m.dist

  ind= (is.infinite(c)==0) & (is.na(c)==0)
  y=y[ind==1]; d=d[ind==1]; m=m[ind==1]; x=x[ind==1,]; zm=zm[ind==1,]; zd=zd[ind==1]; c=c[ind==1]; dzd=dzd[ind==1]; one_dzd=one_dzd[ind==1]; pscore2=pscore2[ind==1]


  regs=cbind(c,m,x)
  if (csquared==TRUE)  regs=cbind(regs, c^2)
    if (logit==FALSE){
      pscore1<-glm(zd~regs,family=binomial(probit))$fitted.values
      pscored<-glm(d~regs,family=binomial(probit))$fitted.values
      pscored[pscored<0]=0; pscored[pscored>1]=1
      pscoredz<-glm(dzd~regs,family=binomial(probit))$fitted.values
      pscoreone_dz<-glm(one_dzd~regs,family=binomial(probit))$fitted.values
      pscoreone_dz[pscoreone_dz<0]=0; pscoreone_dz[pscoreone_dz>1]=1
    }
    if (logit==TRUE){
      pscore1<-glm(zd~regs,family=binomial(logit))$fitted.values
      pscored<-glm(d~regs,family=binomial(logit))$fitted.values
      pscored[pscored<0]=0; pscored[pscored>1]=1
      pscoredz<-glm(dzd~regs,family=binomial(logit))$fitted.values
      pscoreone_dz<-glm(one_dzd~regs,family=binomial(logit))$fitted.values
      pscoreone_dz[pscoreone_dz<0]=0; pscoreone_dz[pscoreone_dz>1]=1
    }


  wgt=zd/pscore2/sum(zd/pscore2)-(1-zd)/(1-pscore2)/sum((1-zd)/(1-pscore2))
  firststage=(d*wgt)
  omega=1-(pscore1-pscore2)/(pscoredz-pscored*pscore2)
  one_omega=1/omega
  ind= (is.infinite(firststage)==0) & (is.infinite(wgt)==0) & (is.infinite(one_omega)==0) & (is.infinite(omega)==0)  &  (is.na(firststage)==0) &  (is.na(wgt)==0) & (is.na(one_omega)==0) & (is.na(omega)==0) & (d*omega*wgt/sum(d*omega*wgt)<=trim) & (d*wgt/sum(d*wgt)<=trim) & ((d-1)*one_omega*wgt/sum((d-1)*one_omega*wgt)<=trim) & ((d-1)*wgt/sum((d-1)*wgt)<=trim)
  y=y[ind==1]; d=d[ind==1]; omega=omega[ind==1]; wgt=wgt[ind==1]; firststage=firststage[ind==1]; one_omega=one_omega[ind==1]
  firststage=sum(firststage)
  y1m0=sum(y*d*omega*wgt)/firststage
  y1m1=sum( y*d*wgt )/firststage
  y0m1=sum(y*(d-1)*one_omega*wgt )/firststage
  y0m0= sum(y*(d-1)*wgt)/firststage
  results=c(y1m1-y0m0, y1m1-y0m1, y1m0-y0m0,  y1m1-y1m0, y0m1-y0m0, templm$coef[2], templm$coef[3]*templm2, nobs-length(y))
  results
}

bootstrap.mediation.late.x<-function(y,d,m,zd,zm,x, boot=1999,trim=0.05, csquared=FALSE, bwreg=bwreg, bwm= bwm, cminobs=40, logit=FALSE, cluster=NULL){
  if (is.null(cluster)){
    obs<-length(y)
    bsamples=matrix( ,boot,8)
    for(i in 1:boot){
      sboot<-sample(1:obs,obs,TRUE)
      yb=y[sboot]; db<-d[sboot]; zdb<-zd[sboot]; mb=m[sboot]
      if (is.null(ncol(zm))) zmb<-zm[sboot]
      if (is.null(ncol(zm))==0) zmb<-zm[sboot,]
      if (is.null(ncol(x))) xb<-x[sboot]
      if (is.null(ncol(x))==0) xb<-x[sboot,]
      bsamples[i,]=effects.late.x(y=yb,d=db,m=mb,zd=zdb, zm=zmb, x=xb, trim=trim, csquared=csquared, bwreg=bwreg, bwm=bwm, cminobs=cminobs, logit=logit)
    }
  }

  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; zdb=c(); zmb=c(); mb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]]);
        zdb<-c(zdb,zd[key==sboot[k]]); mb<-c(mb,m[key==sboot[k]])
        if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
        if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
        if (is.null(ncol(zm))) zmb<-c(zmb,zm[key==sboot[k]])
        if (is.null(ncol(zm))==0) zmb=rbind(zmb,zm[key==sboot[k],])
      }
      est=effects.late.x(y=yb,d=db,m=mb,zd=zdb, zm=zmb, x=xb, trim=trim, csquared=csquared, bwreg=bwreg, bwm=bwm, cminobs=cminobs, logit=logit)
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0){
    cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  }
  bsamples
}


mediation.cont<-function(y,d,m,x, d0, d1, ATET=FALSE, trim=0.05, lognorm=FALSE, bw){
  if(lognorm==TRUE){
    dd=d;dd[d==0]=0.00001
    ggg=glm(log(dd)~x)
    if (d0==0) d0=0.00001; if (d1==0) d1=0.00001
    pscore1d0=(dnorm( (log(d0)-cbind(1,x)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/d0)
    pscore1d1=(dnorm( (log(d1)-cbind(1,x)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/d1)
    ggg=glm(log(dd)~cbind(x,m))
    pscore2d0=(dnorm( (log(d0)-cbind(1,x,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/d0)
    pscore2d1=(dnorm( (log(d1)-cbind(1,x,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/d1)
  }
  if(lognorm==FALSE){
    ggg=glm(d~x)
    pscore1d0=(dnorm( (d0-cbind(1,x)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
    pscore1d1=(dnorm( (d1-cbind(1,x)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
    ggg=glm(d~cbind(x,m))
    pscore2d0=(dnorm( (d0-cbind(1,x,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
    pscore2d1=(dnorm( (d1-cbind(1,x,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
  }
  kernwgtd0=npksum(bws=bw, txdat = d, tydat = y, exdat = d0, return.kernel.weights=TRUE, ckertype="epanechnikov", ckerorder=2)$kw
  kernwgtd1=npksum(bws=bw, txdat = d, tydat = y, exdat = d1, return.kernel.weights=TRUE, ckertype="epanechnikov", ckerorder=2)$kw
  if(lognorm==FALSE) ind=((pscore2d0>=trim) & (pscore2d1>=trim) )
  if(lognorm==TRUE)  ind=((pscore2d0>=quantile(pscore2d0, trim)) & (pscore2d1>=quantile(pscore2d1, trim)) )
  y=y[ind];  pscore1d0=pscore1d0[ind];pscore1d1=pscore1d1[ind];pscore2d0=pscore2d0[ind]; pscore2d1=pscore2d1[ind];
  kernwgtd0=kernwgtd0[ind]; kernwgtd1= kernwgtd1[ind]
  if (ATET==FALSE){
    yd1m1=sum(y*kernwgtd1/pscore1d1)/sum(kernwgtd1/pscore1d1)
    yd1m0=sum(y*kernwgtd1*pscore2d0/(pscore2d1*pscore1d0))/sum(kernwgtd1*pscore2d0/(pscore2d1*pscore1d0))
    yd0m0=sum(y*kernwgtd0/pscore1d0)/sum(kernwgtd0/pscore1d0)
    yd0m1=sum(y*kernwgtd0*pscore2d1/(pscore2d0*pscore1d1))/sum(kernwgtd0*pscore2d1/(pscore2d0*pscore1d1))
  }
  if (ATET==TRUE){
    yd1m1=sum(y*kernwgtd1)/sum(kernwgtd1)
    yd1m0=sum(y*kernwgtd1*pscore2d0*pscore1d1/(pscore2d1*pscore1d0))/sum(kernwgtd1*pscore2d0*pscore1d1/(pscore2d1*pscore1d0))
    yd0m0=sum(y*kernwgtd0*pscore1d1/pscore1d0)/sum(kernwgtd0*pscore1d1/pscore1d0)
    yd0m1=sum(y*kernwgtd0*pscore2d1/pscore2d0)/sum(kernwgtd0*pscore2d1/pscore2d0)
  }
  results=c(yd1m1-yd0m0, yd1m1 - yd0m1, yd1m0-yd0m0, yd1m1 - yd1m0, yd0m1 - yd0m0, sum(1-ind))
}

bootstrap.mediation.cont<-function(y,d,m,x,d0,d1,ATET=FALSE, trim=0.05, lognorm=FALSE, bw, boot=1999, cluster=NULL){
  if (is.null(cluster)){
    obs<-length(y)
    bsamples=matrix( ,boot,6)
    for(i in 1:boot){
      sboot<-sample(1:obs,obs,TRUE)
      yb=y[sboot]
      db<-d[sboot]
      if (is.null(ncol(m))) mb<-m[sboot]
      if (is.null(ncol(m))==0) mb<-m[sboot,]
      if (is.null(ncol(x))) xb<-x[sboot]
      if (is.null(ncol(x))==0) xb<-x[sboot,]
      bsamples[i,]=c(mediation.cont(y=yb,d=db,m=mb,x=xb, d0=d0, d1=d1,  ATET=ATET, trim=trim, lognorm=lognorm, bw=bw))
    }
  }
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      db<-c(); yb<-c(); xb<-c() ; mb=c()
      for (k in 1:length(sboot)) {
        db<-c(db,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]])
        if (is.null(ncol(m))) mb<-c(mb,m[key==sboot[k]])
        if (is.null(ncol(m))==0) mb=rbind(mb,m[key==sboot[k],])
        if (is.null(ncol(x))) xb<-c(xb,x[key==sboot[k]])
        if (is.null(ncol(x))==0) xb=rbind(xb,x[key==sboot[k],])
      }
      est=c(mediation.cont(y=yb,d=db,m=mb,x=xb, d0=d0, d1=d1,  ATET=ATET, trim=trim, lognorm=lognorm, bw=bw))
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  bna=apply(bsamples, 1, sum)
  bsamples=bsamples[is.na(bna)==0,]
  if (sum(is.na(bna))>0) cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  bsamples
}

