ynegbinompowersim<-function(nsize=200,r0=1.0,r1=0.5,shape0=1,shape1=shape0,pi1=0.5,alpha=0.05,twosided=1,fixedfu=1,type=1,u=c(0.5,0.5,1),ut=c(0.5,1.0,1.5),tfix=ut[length(ut)]+0.5,maxfu=10.0,tchange=c(0,0.5,1),ratec1=c(0.15,0.15,0.15),ratec0=ratec1,rn=10000){
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
  ##rn: number of simulations
  
  wtest=rep(0,rn)
  
  for (r in 1:rn){
    zz<-rbinom(nsize,size=1,prob=pi1)
    yy=rep(0,nsize)
    n1=sum(zz==1);n0=sum(zz==0)
    rt1=rep(fixedfu,n1);rt0=rep(fixedfu,n0)
    tneg=rep(0,nsize)
    if (type==2){
      rt1=pmin(rpwe(n1,rate=ratec1,tchange=tchange)$r,fixedfu)
      rt0=pmin(rpwe(n0,rate=ratec0,tchange=tchange)$r,fixedfu)
    } 
    else if (type==3){
      rt1<-pmin(tfix-rpwu(nr=n1,u=u,ut=ut)$r,maxfu)
      rt0<-pmin(tfix-rpwu(nr=n0,u=u,ut=ut)$r,maxfu)
    }
    else if (type==4){
      eu1<-pmin(tfix-rpwu(nr=n1,u=u,ut=ut)$r,maxfu)
      eu0<-pmin(tfix-rpwu(nr=n0,u=u,ut=ut)$r,maxfu)
      rt1<-pmin(rpwe(n1,rate=ratec1,tchange=tchange)$r,eu1)
      rt0<-pmin(rpwe(n0,rate=ratec0,tchange=tchange)$r,eu0)
    }
    a0=1/(1+r0*shape0*rt0);a1=1/(1+r1*shape1*rt1)
    rneg0=rnbinom(n0,size=1/shape0,prob=a0)
    rneg1=rnbinom(n1,size=1/shape1,prob=a1)
    yy[zz==1]=rneg1;yy[zz==0]=rneg0
    tneg[zz==1]=rt1;tneg[zz==0]=rt0
    #abc0=negml(yy,zz,tt=tneg)
    #abc <- glm.nb(yy ~ zz+offset(log(tneg)),start=c(abc0$beta,abc0$gamma),init.theta=abc0$k,control=glm.control(maxit=0))
    abc <- glm.nb(yy ~ zz+offset(log(tneg)),start=c(log(r0),log(r1/r0)),init.theta=1/shape0,control=glm.control(maxit=30)) 
    #if (abc$theta.warn=="iteration limit reached"){
    #   abc<-glm(yy ~ zz+offset(log(tneg)),family=poisson)
    #}
    wtest[r]=summary(abc)$coefficients[2,1]/summary(abc)$coefficients[2,2]
  }
  ptest=(1-pnorm(abs(wtest)))
  if (twosided==1)ptest=2*(1-pnorm(abs(wtest)))
  power=sum(ptest<alpha)/rn
  list(power=power*100)         
}
