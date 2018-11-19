\name{negint2}
\alias{negint2}
\title{A utility functon to calculate the mean exposure under different scenarios}
\description{This will calculate the mean exposure under different scenarios: 
2: fixed follow-up with drop-out, 3: variable follow-up with a maximum (maxfu), 4: variable follow-up with a maximum and drop-out
}
\usage{
negint2(ux=0.5,fixedfu=1,type=2,u=c(0.5,0.5,1),ut=c(0.5,1.0,1.5),
  tfix=ut[length(ut)]+0.5,maxfu=10.0,tchange=c(0,0.5,1),
  ratec=c(0.15,0.15,0.15),eps=1.0e-03)
}
\arguments{
  \item{ux}{the parameter a in (a*t)/(1+a*t)}
  \item{fixedfu}{the minimum follow-up time}
  \item{type}{follow-up type, type=2: fixed fu with fu time \code{fixedfu} but subject to censoring; type=3: depending on entry time, minimum fu is \code{fixedfu} and maximum fu is \code{maxfu}; type=4: same as 3 but subject to censoring}
    \item{u}{recruitment rate}
  \item{ut}{recruitment interval, must have the same length as \code{u}}
    \item{tfix}{fixed study duration, often equals to recruitment time plus minimum follow-up}
  \item{maxfu}{maximum follow-up time, should not be greater than \code{tfix}}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which the drop-out rate changes. The first element of tchange must be zero. The above rates and \code{tchange} must have the same length.}
  \item{ratec}{piecewise constant drop-out rate}
  \item{eps}{error tolerance for the numerical intergration}
}
\details{Let \eqn{\tau_{min}} and \eqn{\tau_{max}} correspond to the minimum follow-up time \code{taumin} and the maximum follow-up time \code{taumax}. Let \eqn{T_f}, \eqn{C}, \eqn{E} and \eqn{R} be the follow-up time, the drop-out time, the study entry time and the total recruitment period(\eqn{R} is the last element of \code{ut}). For type 2 follow-up \eqn{T_f=min(C,\tau_{min})}. For type 3 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max})}. For type 4 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max},C)}. Let \eqn{f} be the density of \eqn{T_f}. We calculate
\deqn{\int_0^{\infty} t f(t)dt}
and
\deqn{\int_0^{\infty} \frac{a t}{1+a t} f(t)dt}
where \eqn{a} is the \code{ux}.
}
\value{\item{mt}{mean of (a*t)/(1+a*t)}
       \item{tt}{mean of t}
       \item{vt}{variance of t}
}
\author{Xiaodong Luo}
%\note{This provides the calculation for the mean and variance of the exposure. }
%\seealso{\code{\link{rpwe}}}
\examples{
##calculating the exposure for type 4 follow-up
exp4=negint2(ux=0.5,fixedfu=1,type=2,u=c(0.5,0.5,1),ut=c(0.5,1.0,1.5),
  tfix=2.0,maxfu=1.0,tchange=c(0,0.5,1),
  ratec=c(0.15,0.15,0.15),eps=1.0e-03)
#mean exposure            
meanexp=exp4$tt
#var exposure
varexp=exp4$vt
c(meanexp,sqrt(varexp))
#mean of (ux*t)/(1+ux*t)
meanuxt=exp4$mt
}
\keyword{piecewise exponential}
\keyword{piecewise uniform}
\keyword{variable follow-up}
