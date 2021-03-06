ynegbinompower<-function(nsize=200,r0=1.0,r1=0.5,shape0=1,shape1=shape0,pi1=0.5,alpha=0.05,twosided=1,fixedfu=1,type=1,u=c(0.5,0.5,1),ut=c(0.5,1.0,1.5),tfix=ut[length(ut)]+0.5,maxfu=10.0,tchange=c(0,0.5,1),ratec1=c(0.15,0.15,0.15),ratec0=ratec1,eps=1.0e-03){
  ####nsize: total number of subjects from both groups
  ####r0,r1: event rates for the control and treatment
  ####shape0, shape1: dispersion parameters for the control and treatment
  ####pi1: allocation probability to the treatment
  ####alpha: type-1 error
  ####two-sided: two-sided(1) vs one-sided(other values)
  ####fixedfu: fixed fu time (0.5,1 year etc)
  ####u: recruitment probabilities within each interval
  ####ut: recruitment intervals
  ####tfix: fixed study time (FPI to LPLV)
  ####maxfu: maximum follow-up time (1yr,1.5 yr etc)
  ####tchange: cut-point times when the drop-out rate changes
  ####ratec1,ratec0: same length as tchange, the drop-out rates fro treatment and control
  #type: follow-up time type,
  ##     type=1: fixed fu with fu time fixedfu
  ##     type=2: fixed fu with fu time fixedfu but subject to censoring
  ##     type=3: depending on entry time, minimum fu is fixedfu and maximum fu is maxfu
  ##     type=4: depending on entry time, minimum fu is fixedfu and maximum fu is maxfu but subject to censoring
  
  
  pi0<-1-pi1
  eff2<-(log(r1/r0))^2
  alpha1<-alpha
  if (twosided==1)alpha1<-alpha/2
  qalpha<--qnorm(alpha1)
  
  tilder0<-tilder1<-rep(0,3)
  tilder0[1]<-tilder1[1]<-r0
  tilder0[2]<-r0;tilder1[2]<-r1
  tilder0[3]<-tilder1[3]<-pi0*r0+pi1*r1
  
  exposure<-SDexp<-rep(0,3)
  tildeXPWR<-XPWR<-tildemineffsize<-mineffsize<-rep(0,3)
  
  if (type==1){
    exposure<-fixedfu
    SDexp<-0
    if (shape0==0){a<-pi0*r0*exposure[1]}
    else {a<-shape0*r0*exposure[1];a<-a/(1+a);a<-pi0/shape0*a}
    if (shape1==0){b<-pi1*r1*exposure[2]}
    else {b<-shape1*r1*exposure[2];b<-b/(1+b);b<-pi1/shape1*b}
    
    mta<-a
    mtb<-b
    V1<-1/a+1/b
    tildeV1<-1/mta+1/mtb
    for (j in 1:3){
      ###j=1,estimate V0 under H0: r0=r1=r
      ###j=2, estimate V0 under alternative
      ###j=3, estimate V0 using MLE
      if (shape0==0){anull<-pi0*tilder0[j]*exposure[1]}
      else {anull<-shape0*tilder0[j]*exposure[1];anull<-anull/(1+anull);anull<-pi0/shape0*anull}
      if (shape1==0){bnull<-pi1*tilder1[j]*exposure[2]}
      else {bnull<-shape1*tilder1[j]*exposure[2];bnull<-bnull/(1+bnull);bnull<-pi1/shape1*bnull}
      
      mtanull=anull
      mtbnull=bnull
      V0<-1/anull+1/bnull
      tildeV0<-1/mtanull+1/mtbnull
      tildeXPWR[j]<-pnorm(-qalpha*sqrt(tildeV0/tildeV1)+sqrt(nsize)*sqrt(eff2/tildeV1))
      XPWR[j]<-pnorm(-qalpha*sqrt(V0/V1)+sqrt(nsize)*sqrt(eff2/V1))
      tildemineffsize[j]<-1-exp(-qalpha*sqrt(tildeV0/nsize))
      mineffsize[j]<-1-exp(-qalpha*sqrt(V0/nsize))
    }
  }
  else {
    abc0=negint2(ux=r0*shape0,fixedfu=fixedfu,type=type,u=u,ut=ut,tfix=tfix,maxfu=maxfu,tchange=tchange,ratec=ratec0,eps=eps)
    abc1=negint2(ux=r1*shape1,fixedfu=fixedfu,type=type,u=u,ut=ut,tfix=tfix,maxfu=maxfu,tchange=tchange,ratec=ratec1,eps=eps)
    
    exposure[1]=abc0$tt;SDexp[1]=sqrt(abc0$vt)
    exposure[2]=abc1$tt;SDexp[2]=sqrt(abc1$vt)
    exposure[3]=pi0*exposure[1]+pi1*exposure[2]
    SDexp[3]=pi0*SDexp[1]+pi1*SDexp[2]
    
    if (shape0==0) {a<-pi0*r0*exposure[1];mta=a}
    else {a<-shape0*r0*exposure[1];a<-a/(1+a);a<-pi0/shape0*a;mta<-pi0/shape0*abc0$mt}
    if (shape1==0) {b<-pi1*r1*exposure[2];mtb=b}
    else {b<-shape1*r1*exposure[2];b<-b/(1+b);b<-pi1/shape1*b;mtb<-pi1/shape1*abc1$mt}
    
    V1<-1/a+1/b
    tildeV1<-1/mta+1/mtb
    for (j in 1:3){
      ###j=1,estimate V0 under H0: r0=r1=r
      ###j=2, estimate V0 under alternative
      ###j=3, estimate V0 using MLE
      def0=negint2(ux=tilder0[j]*shape0,fixedfu=fixedfu,type=type,u=u,ut=ut,tfix=tfix,maxfu=maxfu,tchange=tchange,ratec=ratec0,eps=eps)
      def1=negint2(ux=tilder1[j]*shape1,fixedfu=fixedfu,type=type,u=u,ut=ut,tfix=tfix,maxfu=maxfu,tchange=tchange,ratec=ratec1,eps=eps)
      
      if (shape0==0){anull<-pi0*tilder0[j]*exposure[1];mtanull=anull}
      else {anull<-shape0*tilder0[j]*exposure[1];anull<-anull/(1+anull);anull<-pi0/shape0*anull;mtanull<-pi0/shape0*def0$mt}
      if (shape1==0){bnull<-pi1*tilder1[j]*exposure[2];mtbnull=bnull}
      else {bnull<-shape1*tilder1[j]*exposure[2];bnull<-bnull/(1+bnull);bnull<-pi1/shape1*bnull;mtbnull<-pi1/shape1*def1$mt}
      
      V0<-1/anull+1/bnull
      tildeV0<-1/mtanull+1/mtbnull
      tildeXPWR[j]<-pnorm(-qalpha*sqrt(tildeV0/tildeV1)+sqrt(nsize)*sqrt(eff2/tildeV1))
      XPWR[j]<-pnorm(-qalpha*sqrt(V0/V1)+sqrt(nsize)*sqrt(eff2/V1))
      tildemineffsize[j]<-1-exp(-qalpha*sqrt(tildeV0/nsize))
      mineffsize[j]<-1-exp(-qalpha*sqrt(V0/nsize))
    }
  }
  list(tildeXPWR=tildeXPWR*100,XPWR=XPWR*100,tildemineffsize=tildemineffsize,mineffsize=mineffsize,Exposure=exposure,SDexp=SDexp)
}
