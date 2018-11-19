negint2<-function(ux=0.5,fixedfu=1,type=2,u=c(0.5,0.5,1),ut=c(0.5,1.0,1.5),tfix=ut[length(ut)]+0.5,maxfu=10.0,tchange=c(0,0.5,1),ratec=c(0.15,0.15,0.15),eps=1.0e-03){
  #ux: the parameter a in (a*t)/(1+a*t)
  #fixedfu: the minimum fu
  #type: follow-up time type,  
  ##     type=2: fixed fu with fu time fixedfu but subject to censoring
  ##     type=3: depending on entry time, minimum fu is fixedfu and maximum fu is maxfu
  ##     type=4: depending on entry time, minimum fu is fixedfu and maximum fu is maxfu but subject to censoring
  ##     u: piecewise constant recruitment rate
  ##     ut: timepoints where recruitment rates change, the same length as u
  ##     tfix: study length, should be the last value of ut plus the minimum fu time fixedfu
  ##     maxfu: maximum fu time
  ##     tchange: the changing points of the piecewise exponential drop-out rates
  ##     ratec: drop-out rates
  ##     eps: error tolerance in the numerical integration
  if (type==2){
    amax=max(ux,ratec,fixedfu)+1
    aseq=seq(0,fixedfu,by=1/amax*eps)
    bseq=c(aseq,0)
    cseq=unique(bseq[bseq<=fixedfu])
    nc=length(cseq)
    cseq0=c(0,cseq[-nc])
    dseq=c((cseq+cseq0)/2,fixedfu)
    nd=length(dseq)
    ss=pwe(t=dseq,rate=ratec,tchange=tchange)$surv
    temp=ux/(1+ux*dseq)^2
    dtemp=cseq-cseq0
    mt=sum(temp[-nd]*ss[-nd]*dtemp)
    tt=sum(ss[-nd]*dtemp)
    vt=2*sum(dseq[-nd]*ss[-nd]*dtemp)-tt^2
  }
  else if (type==3){
    amax=max(ux,u,ut,tfix,maxfu)+1
    aseq=seq(0,maxfu,by=1/amax*eps)
    bseq=c(aseq,0)
    cseq=unique(bseq[bseq<=maxfu])
    nc=length(cseq)
    cseq0=c(0,cseq[-nc])
    dseq=c((cseq+cseq0)/2,maxfu)
    nd=length(dseq)
    ss=pwu(t=tfix-dseq,u=u,ut=ut)$dist
    temp=ux/(1+ux*dseq)^2
    dtemp=cseq-cseq0
    mt=sum(temp[-nd]*ss[-nd]*dtemp)
    tt=sum(ss[-nd]*dtemp)
    vt=2*sum(dseq[-nd]*ss[-nd]*dtemp)-tt^2
  }
  else if (type==4){
    amax=max(ux,u,ut,tfix,maxfu,ratec)+1
    aseq=seq(0,maxfu,by=1/amax*eps)
    bseq=c(aseq,0)
    cseq=unique(bseq[bseq<=maxfu])
    nc=length(cseq)
    cseq0=c(0,cseq[-nc])
    dseq=c((cseq+cseq0)/2,maxfu)
    nd=length(dseq)
    ss=pwu(t=tfix-dseq,u=u,ut=ut)$dist*pwe(t=dseq,rate=ratec,tchange=tchange)$surv
    temp=ux/(1+ux*dseq)^2
    dtemp=cseq-cseq0
    mt=sum(temp[-nd]*ss[-nd]*dtemp)
    tt=sum(ss[-nd]*dtemp)
    vt=2*sum(dseq[-nd]*ss[-nd]*dtemp)-tt^2
  }
  list(mt=mt,tt=tt,vt=pmax(vt,0))
  ##mt: mean of (a*t)/(1+a*t)
  ##tt: mean of t
  ##vt: variance of t
}
