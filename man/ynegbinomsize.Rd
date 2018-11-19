\name{ynegbinomsize}
\alias{ynegbinomsize}
\title{Two-sample sample size calculation for negative binomial distribution with variable follow-up}
\description{This will calculate the sample size for the negative binomial distribution for the 2-sample case under different follow-up scenarios: 1: fixed follow-up, 2: fixed follow-up with drop-out, 3: variable follow-up with a minimum fu and a maximum fu, 4: variable follow-up with a minimum fu and a maximum fu and drop-out.
}
\usage{
ynegbinomsize(r0=1.0,r1=0.5,shape0=1,shape1=shape0,pi1=0.5,
      alpha=0.05,twosided=1,beta=0.2,fixedfu=1,
      type=1,u=c(0.5,0.5,1),ut=c(0.5,1.0,1.5),tfix=ut[length(ut)]+0.5,maxfu=10.0,
      tchange=c(0,0.5,1),ratec1=c(0.15,0.15,0.15),ratec0=ratec1,eps=1.0e-03)
}
\arguments{
  \item{r0}{event rate for the control}
  \item{r1}{event rate for the treatment}
  \item{shape0}{dispersion parameter for the control}
  \item{shape1}{dispersion parameter for the treatment}
  \item{pi1}{allocation prob for the treatment}
  \item{alpha}{type-1 error}
  \item{twosided}{1: two-side, others: one-sided}
  \item{beta}{tyep-2 error}
  \item{fixedfu}{fixed follow-up time for each patient}
  \item{type}{follow-up time type, type=1: fixed fu with fu time \code{fixedfu};  type=2: same as 1 but subject to censoring; type=3: depending on entry time, minimum fu is \code{fixedfu} and maximum fu is \code{maxfu}; type=4: same as 3 but subject to censoring}
  \item{u}{recruitment rate}
  \item{ut}{recruitment interval, must have the same length as \code{u}}
  \item{tfix}{fixed study duration, often equals to recruitment time plus minimum follow-up \code{fixedfu}}
  \item{maxfu}{maximum follow-up time, should not be greater than \code{tfix}}
  \item{tchange}{a strictly increasing sequence of time points starting from zero at which the drop-out rate changes. The first element of tchange must be zero. }
  \item{ratec1}{piecewise constant drop-out rate for the treatment. The rate and \code{tchange} must have the same length.}
  \item{ratec0}{piecewise constant drop-out rate for the control. The rate and \code{tchange} must have the same length.}
  \item{eps}{error tolerance for the numerical intergration}
}
\details{Let \eqn{\tau_{min}} and \eqn{\tau_{max}} correspond to the minimum follow-up time \code{fixedfu} and the maximum follow-up time \code{maxfu}. Let \eqn{T_f}, \eqn{C}, \eqn{E} and \eqn{R} be the follow-up time, the drop-out time, the study entry time and the total recruitment period(\eqn{R} is the last element of \code{ut}). For type 1 follow-up, \eqn{T_f=\tau_{min}}. For type 2 follow-up \eqn{T_f=min(C,\tau_{min})}. For type 3 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max})}. For type 4 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max},C)}. Let \eqn{f} be the density of \eqn{T_f}. 
Suppose that \eqn{Y_i} is the number of event obsevred in follow-up time \eqn{t_i} for patient \eqn{i} with treatment assignment \eqn{Z_i}, \eqn{i=1,\ldots,n}. Suppose that \eqn{Y_i} follows a negative binomial distribution such that 
\deqn{P(Y_i=y\mid Z_i=j)=\frac{\Gamma(y+1/k_j)}{\Gamma(y+1)\Gamma(1/k_j)}\Bigg(\frac{k_ju_i}{1+k_ju_i}\Bigg)^y\Bigg(\frac{1}{1+k_ju_i}\Bigg)^{1/k_j},}
where 
\deqn{\log(u_i)=\log(t_i)+\beta_0+\beta_1 Z_i.}
Let \eqn{\hat{\beta}_0} and \eqn{\hat{\beta}_1} be the MLE of \eqn{\beta_0} and \eqn{\beta_1}. 
The varaince of \eqn{\hat{\beta}_1} is 
\deqn{\mbox{var}(\hat{\beta}_1)=1/\tilde{a}_0(r_0)+1/\tilde{a}_1(r_1)}
where 
\deqn{\tilde{a}_j(r)=\sum_{i=1}^n I(Z_i=j)k_jrt_i/(1+k_jrt_i), \hspace{0.5cm}j=0,1,}
and \eqn{k_j, j=0,1} are the dispersion parameters for control \eqn{j=0} and treatment \eqn{j=1}. Note that Zhu and Lakkis (2014) use 
\deqn{a_j(r)=\sum_{i=1}^n I(Z_i=j)k_jrE(t_i)/\{1+k_jrE(t_i)\}, }
to replace \eqn{\tilde{a}_j(r)}, \eqn{j=0,1}. Using Jensen's inequality, we can show \eqn{a_j(r)\ge \tilde{a}_j(r)}, which means 
Zhu and Lakkis's method will underestimate variance of \eqn{\hat{\beta}_1}, which leads to either smaller than required sample size or inflated power. For comparison, I provide sample sizes under both \eqn{\tilde{a}_j(r)} and \eqn{a_j(r)}. 

Zhu and Lakkis (2014) discuss three types of the variance under the null. The first way is to set \eqn{\tilde{r}_0=\tilde{r}_1=r_0}, using event rate from the control group. The second way is to set \eqn{\tilde{r}_0=r_0, \tilde{r}_1=r_1}, using true event rates. The third way is to set \eqn{\tilde{r}_0=\tilde{r}_1=\tilde{r}}, where \eqn{\tilde{r}=\pi_1 r_1+\pi_0 r_0}, using maximum likelihood estimation. 

Therefore, for each type of follow-up, there are 3 sample sizes calculated (because there are 3 varainces under the null) for with and without approximation of Zhu and Lakkis (2014). 

Note that PASS14.0 provides 3 ways of null varaince with the default being the MLE. PASS does not allow different dispersion parameters between treatmetn and control. EAST only provides the second way of null varaince but allows for different dispersion parameters. Both of these softwares base on the approximatin method of Zhu and Lakkis (2014), which underestimate the required sample sizes.    
}
\value{
\item{tildeXN}{sample sizes based on current approach, i.e. not based on the Zhu and Lakkis's approximation}
\item{XN}{sample sizes based on the Zhu and Lakkis's approximation}
\item{Exposure}{mean exposure under different follow-up types with element 1 for control, element 2 for treatment and element 3 for overall.}
\item{SDExp}{Sd of the exposure under different follow-up types with element 1 for control, element 2 for treatment and column 3 for overall.}
}
\author{Xiaodong Luo}
\examples{
##calculating the sample sizes
abc=ynegbinomsize(r0=1.0,r1=0.5,shape0=1,pi1=0.5,alpha=0.05,twosided=1,
    beta=0.2,fixedfu=1,type=4,u=c(0.5,0.5,1),ut=c(0.5,1.0,1.5),
    tfix=1.5,maxfu=1,tchange=c(0,0.5,1),ratec1=c(0.15,0.15,0.15),
    eps=1.0e-03)
###Zhu and Lakkis's sample sizes (i.e. with approximation) 
abc$XN
###Our sample sizes (i.e. without approximation)
abc$tildeXN
}
\references{Zhu~H and Lakkis~H. Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 2014, 33: 376-387.}
\keyword{negative binomial}
\keyword{sample size}
\keyword{piecewise exponential}
\keyword{piecewise uniform}

