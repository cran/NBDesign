\name{ynegbinompowersim}
\alias{ynegbinompowersim}
\title{Two-sample sample size calculation for negative binomial distribution with variable follow-up}
\description{This will calculate the power for the negative binomial distribution for the 2-sample case under different follow-up scenarios: 1: fixed follow-up, 2: fixed follow-up with drop-out, 3: variable follow-up with a minimum fu and a maximum fu, 4: variable follow-up with a minimum fu and a maximum fu and drop-out.
}
\usage{
ynegbinompowersim(nsize=200,r0=1.0,r1=0.5,shape0=1,shape1=shape0,pi1=0.5,
   alpha=0.05,twosided=1,fixedfu=1,type=1,u=c(0.5,0.5,1),ut=c(0.5,1.0,1.5),
   tfix=ut[length(ut)]+0.5,maxfu=10.0,tchange=c(0,0.5,1),
   ratec1=c(0.15,0.15,0.15),ratec0=ratec1,rn=10000)
}
\arguments{
  \item{nsize}{total number of subjects in two groups}
  \item{r0}{event rate for the control}
  \item{r1}{event rate for the treatment}
  \item{shape0}{dispersion parameter for the control}
  \item{shape1}{dispersion parameter for the treatment}
  \item{pi1}{allocation prob for the treatment}
  \item{alpha}{type-1 error}
  \item{twosided}{1: two-side, others: one-sided}
  \item{fixedfu}{fixed follow-up time for each patient}
  \item{type}{follow-up time type, type=1: fixed fu with fu time \code{fixedfu};  type=2: same as 1 but subject to censoring; type=3: depending on entry time, minimum fu is \code{fixedfu} and maximum fu is \code{maxfu}; type=4: same as 3 but subject to censoring}
  \item{u}{recruitment rate}
  \item{ut}{recruitment interval, must have the same length as \code{u}}
  \item{tfix}{fixed study duration, often equals to recruitment time plus minimum follow-up}
  \item{maxfu}{maximum follow-up time, should not be greater than \code{tfix}}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which the drop-out rate changes. The first element of tchange must be zero. The above rates and \code{tchange} must have the same length.}
  \item{ratec1}{piecewise constant drop-out rate for the treatment}
  \item{ratec0}{piecewise constant drop-out rate for the control}
  \item{rn}{Number of repetitions}
}
\details{Let \eqn{\tau_{min}} and \eqn{\tau_{max}} correspond to the minimum follow-up time \code{fixedfu} and the maximum follow-up time \code{maxfu}. Let \eqn{T_f}, \eqn{C}, \eqn{E} and \eqn{R} be the follow-up time, the drop-out time, the study entry time and the total recruitment period(\eqn{R} is the last element of \code{ut}). For type 1 follow-up, \eqn{T_f=\tau_{min}}. For type 2 follow-up \eqn{T_f=min(C,\tau_{min})}. For type 3 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max})}. For type 4 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max},C)}. Let \eqn{f} be the density of \eqn{T_f}. 
Suppose that \eqn{Y_i} is the number of event obsevred in follow-up time \eqn{t_i} for patient \eqn{i} with treatment assignment \eqn{Z_i}, \eqn{i=1,\ldots,n}. Suppose that \eqn{Y_i} follows a negative binomial distribution such that 
\deqn{P(Y_i=y\mid Z_i=j)=\frac{\Gamma(y+1/k_j)}{\Gamma(y+1)\Gamma(1/k_j)}\Bigg(\frac{k_ju_i}{1+k_ju_i}\Bigg)^y\Bigg(\frac{1}{1+k_ju_i}\Bigg)^{1/k_j},}
where \eqn{k_j, j=0,1} are the dispersion parameters for control \eqn{j=0} and treatment \eqn{j=1} and
\deqn{\log(u_i)=\log(t_i)+\beta_0+\beta_1 Z_i.}
The data will be gnerated according to the above model. Note that the piecewise exponential distribution and the piecewise uniform distribution are genrated using R package PWEALL functions "rpwe" and "rpwu", respectively. 

The parameters in the model are estimated by MLE using the R package MASS function "glm.nb". 
}
\value{
\item{power}{simulation power (in percentage)}
}
\author{Xiaodong Luo}
\examples{
##calculating the sample sizes
abc=ynegbinompowersim(nsize=200,r0=1.0,r1=0.5,shape0=1,
        pi1=0.5,alpha=0.05,twosided=1,fixedfu=1,
        type=4,u=c(0.5,0.5,1),ut=c(0.5,1.0,1.5),
        tchange=c(0,0.5,1),
        ratec1=c(0.15,0.15,0.15),rn=10)
###Power
abc$power
}
\keyword{negative binomial}
\keyword{power}
\keyword{piecewise exponential}
\keyword{piecewise uniform}

