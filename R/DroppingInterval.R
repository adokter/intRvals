#' Analysing time-ordered event data with missed observations
#'
#' \pkg{intRvals} calculates means and variances of arrival intervals (and arrival rates)
#' corrected for missed arrival observations, and compares means
#'  and variances of groups of interval data.
#' @details
#' \subsection{General}{
#' The central function of package \pkg{intRvals} is
#' \code{\link{estinterval}}, which is used to estimate the
#' mean arrival interval (and its standard deviation) from interval
#' data with missed arrivals. This is
#' achieved by fitting the theoretical probability density
#' \code{\link{intervalpdf}} to the interval data
#'
#' The package can be used to analyse general interval data where
#' intervals are derived from distinct arrival observations.
#' For example, the authors have used it to analyze dropping intervals
#' of grazing geese for estimating their faecal output.
#'
#' Intervals are defined as the time between observed arrival events (e.g. the time between one excreted droppings to the next)
#' The package provides a way of taking into account missed observations
#' (e.g. defecations), which lead to occasional observed intervals at integer multiples of the
#' true arrival interval.
#' }
#' \subsection{Typical workflow}{
#' \enumerate{
#'   \item{Fit interval model \code{m} to an interval dataset \code{d} using \code{\link{estinterval}}, as in \code{m=estinterval(d)}}.
#'   \item{Visually inspect model fits using \code{\link{plot.intRvals}}, as in \code{plot(m)}}.
#'   \item{Use \code{\link{anova.intRvals}} to check whether the missed event probability was signficantly different from zero, as in \code{anova(m)}.}
#'   \item{Also use \code{\link{anova.intRvals}} to perform model selection between competing models \code{m1},\code{m2} for the same interval dataset \code{d}, as in \code{anova(m1,m2)}.}
#'   \item{Compare means and variances between different interval datasets \code{d1},\code{d2} using \code{\link{ttest}} and \code{\link{vartest}}.}
#' }
#' }
#' \subsection{Other useful functionality}{
#' \code{\link{fold}} provides functionality to fold observed intervals back to their fundamental interval
#'
#' \code{\link{fundamental}} tests which intervals are fundamental, i.e. intervals not containing a missed arrival observation
#'
#' \code{\link{interval2rate}} converts interval estimates to rates
#'
#' \code{\link{partition}} estimates and tests for the presence of within-subject variation
#'
#' \code{\link{intervalsim}} simulates a set of observed intervals
#'
#' The package comes with a example interval dataset \code{\link{goosedrop}}
#'
#' Please cite this package using the publication "Analysing time-ordered event data with missed observations, Ecology and Evolution, 2017" by Dokter et al.
#' }
#'
#' @references
#' Dokter, A.M., van Loon, E.E., Fokkema, W., Lameris, T.K, , Nolet, B.A. and van der Jeugd, H.P. 2017. Analysing time-ordered event data with missed observations, Ecology and Evolution, 2017, in press.
#'
#' B\'{e}dard, J. & Gauthier, G. 1986. Assessment of faecal output in geese. Journal of Applied Ecology, 23, 77-90.
#'
#' Owen, M. 1971. The Selection of Feeding Site by White-Fronted Geese in Winter. Journal of Applied Ecology 8: 905-917.
#'
#' @import plyr
#' @import lme4
#' @importFrom graphics curve hist
#' @importFrom stats anova dgamma dnorm integrate lm median optim pchisq pf pgamma plogis pnorm pt qf qt rexp rgamma rnorm runif sd setNames

"_PACKAGE"
#> [1] "_PACKAGE"

## to generate pdf of the manual do:
# R CMD Rd2pdf '/Library/Frameworks/R.framework/Versions/3.3/Resources/library/intRvals'

## Note: convergence is still conditional on the starting values provided by users:
# this works:
# dr.schiermonnikoog=estinterval(goosedrop[goosedrop$site=="schiermonnikoog",]$interval,fun="normal",mu=247,sigma=83)
# this not:
# dr.schiermonnikoog=estinterval(goosedrop[goosedrop$site=="schiermonnikoog",]$interval,fun="normal")

#terschelling=read.csv("~/Dropbox/metawad/veldwerk terschelling/dropping_intervals_terschelling_elaborate.csv",colClasses=c("POSIXct","numeric","numeric","numeric","integer"))
#schiermonnikoog=read.csv("~/Dropbox/metawad/veldwerk terschelling/dropping_intervals_schier_elaborate.csv",colClasses=c("POSIXct","numeric","numeric","numeric","integer"))
#schiermonnikoog$bout_id=schiermonnikoog$bout_id+max(terschelling$bout_id)
#terschelling$site="terschelling"
#schiermonnikoog$site="schiermonnikoog"
#goosedrop=rbind(terschelling,schiermonnikoog)
# break up data in periods according to date
MidPoints=as.POSIXct(c("2013-3-19", "2013-4-2", "2013-4-15", "2013-4-29", "2013-5-14"))
#goosedrop$period=NA
#goosedrop[goosedrop$date>=MidPoints[1] & goosedrop$date<MidPoints[2],]$period=1
#goosedrop[goosedrop$date>=MidPoints[2] & goosedrop$date<MidPoints[3],]$period=2
#goosedrop[goosedrop$date>=MidPoints[3] & goosedrop$date<MidPoints[4],]$period=3
#goosedrop[goosedrop$date>=MidPoints[4] & goosedrop$date<MidPoints[5],]$period=4
#goosedrop[goosedrop$date>=MidPoints[5],]$period=5
#save(goosedrop,file="~/git/R/intRvals/data/goosedrop.RData")
#load("~/git/R/intRvals/data/goosedrop.RData")
#' Dataset with dropping intervals observed for foraging Brent Geese (Branta bernicla bernicla)
#'
#' The dataset contains observations from two sites: the island of Schiermonnikoog (saltmarsh) and Terschelling (agricultural grassland).
#' Brent geese were observed continuously with spotting scopes, and the time when geese excreted a dropping was written down.
#' The time in seconds between wo subsequently observed dropping arrivals of a single individual refers to one dropping interval.
#' The variables are as follows:
#' \describe{
#'   \item{\code{date}}{observation start time of the interval}
#'   \item{\code{interval}}{length of the interval in seconds}
#'   \item{\code{bout_length}}{total observation time of individual}
#'   \item{\code{prop_abdomen_seen}}{proportion of total observation time when abdomen could be observed}
#'   \item{\code{bout_id}}{intervals belonging to the same observation bout of an individual have the same \code{bout_id}}
#'   \item{\code{site}}{observation site. One of 'terschelling' (agricultural grassland) or 'schiermonnikoog' (salt marsh)}
#'   \item{\code{period}}{two-week observation period (1-5)}
#' }
#'
#' @author Adriaan Dokter \email{a.m.dokter@uva.nl}
"goosedrop"

# small number, used to maintain numerical stability
eps=.Machine$double.xmin

# normal distribution for (i-1) missed arrivals component of the
# probability density function (PDF)
normi=  function(x,mu,sigma,p,i,sigmap=sigma) (p^(i-1)-p^i)*dnorm(x,i*mu,sqrt((i-1)*sigma^2+sigmap^2))
normpi= function(x,mu,sigma,p,i,sigmap=sigma) (p^(i-1)-p^i)*pnorm(x,i*mu,sqrt((i-1)*sigma^2+sigmap^2))
# gamma distribution for (i-1) missed arrivals component of the PDF; sigmap not implemented, dummy only
#gammai= function(x,mu,sigma,p,i,sigmap=sigma) (p^(i-1)-p^i)*dgamma(x,shape=i*mu^2/(sigma^2),scale=(sigma^2)/abs(mu))
#gammapi=function(x,mu,sigma,p,i,sigmap=sigma) (p^(i-1)-p^i)*pgamma(x,shape=i*mu^2/(sigma^2),scale=(sigma^2)/abs(mu))
gammai= function(x,mu,sigma,p,i,sigmap=sigma) (p^(i-1)-p^i)*dgamma(x,shape=i^2*mu^2/((i-1)*sigma^2+sigmap^2),scale=((i-1)*sigma^2+sigmap^2)/abs(i*mu))
gammapi=function(x,mu,sigma,p,i,sigmap=sigma) (p^(i-1)-p^i)*pgamma(x,shape=i^2*mu^2/((i-1)*sigma^2+sigmap^2),scale=((i-1)*sigma^2+sigmap^2)/abs(i*mu))


# normalized list of all components i of the PDF
pdfcomponents=function(x,mu,sigma,p,N=5L,fun=normi,sigmap=sigma) sapply(1:N,FUN=fun,x=x,p=p,mu=mu,sigma=sigma,sigmap=sigmap)

# the cumulative density function
cdf=function(x,mu,sigma,p,N=5L,fun=normpi,fpp=0,sigmap=sigma) {
  if(length(x)>1){
    if(fpp==0){
      return(rowSums(pdfcomponents(x,mu,sigma,p,N,fun=fun,sigmap=sigmap)))
    } else return((1-fpp)*rowSums(pdfcomponents(x,mu,sigma,p,N,fun=fun,sigmap=sigmap)) + fpp*rowSums(pdfcomponents(x,mu,mu,p,N,fun=gammapi)))
  } else{
    if(fpp==0){
      return(sum(pdfcomponents(x,mu,sigma,p,N,fun=fun,sigmap=sigmap)))
    }
    else return((1-fpp)*sum(pdfcomponents(x,mu,sigma,p,N,fun=fun,sigmap=sigmap))+fpp*sum(pdfcomponents(x,mu,mu,p,N,fun=gammapi)))
  }
}

# the PDF
probdens=function(x,mu,sigma,p,N=5L,fun=normi,fpp=0,funcdf=normpi,trunc=c(0,Inf),sigmap=sigma) {
  if(is.na(sigmap)) sigmap=sigma
  if(is.numeric(trunc)){
    normfactorA=cdf(trunc[2],mu,sigma,p,N,fun=funcdf,fpp=fpp,sigmap=sigmap)-cdf(trunc[1],mu,sigma,p,N,fun=funcdf,fpp=fpp,sigmap=sigmap)
    normfactorB=cdf(trunc[2],mu,mu,p,N,fun=funcdf,fpp=fpp,sigmap=sigmap)-cdf(trunc[1],mu,mu,p,N,fun=funcdf,fpp=fpp,sigmap=sigmap)
    normfactor=(1-fpp)*normfactorA+fpp*normfactorB
  } else normfactor=1
  if(length(x)>1){
    if(fpp==0){
      return(rowSums(pdfcomponents(x,mu,sigma,p,N,fun=fun,sigmap=sigmap))/normfactor)
    } else return(((1-fpp)*rowSums(pdfcomponents(x,mu,sigma,p,N,fun=fun,sigmap=sigmap)) + fpp*rowSums(pdfcomponents(x,mu,mu,p,N,fun=gammai)))/normfactor)
  } else{
    if(fpp==0){
      return(sum(pdfcomponents(x,mu,sigma,p,N,fun=fun,sigmap=sigmap))/normfactor)
    }
    else return(((1-fpp)*sum(pdfcomponents(x,mu,sigma,p,N,fun=fun,sigmap=sigmap))+fpp*sum(pdfcomponents(x,mu,mu,p,N,fun=gammai)))/normfactor)
  }
}

# check of common input arguments
checkargs = function(data=c(1,2),mu=1,sigma=1,p=0.2,N=5L,n=1L,fun="gamma",trunc=c(0,Inf),fpp=0,fpp.method="fixed",p.method="auto", conf.level = 0.95,sigma.within=NA,n.ind=NA,group=NA){
  if(!is.numeric(data)) stop("'data' should be numeric vector with intervals")
  if(min(data)<0) stop("data contains one or more negative intervals, only positive intervals allowed.")
  if(!is.na(sigma.within) & sigma.within!="auto"){
    if(!missing(group) & length(group)!=length(data)) stop("'group' should be a vector of equal length as 'data'")
    # need this only in intervalsim function ... if(is.na(n.ind)) stop("n.ind should be a numeric integer")
    if(sigma.within>sigma) stop(paste("within-subject standard deviation 'sigma.within' (",sigma.within,") should be smaller than total standard deviation 'sigma' (",sigma,")",sep=""))
    if(sigma.within<=0) stop("within-subject standard deviation 'sigma.within' should be a positive numeric value")
  }
  if (!missing(mu) & (length(mu) != 1 || mu<=0 || is.na(mu))) stop("'mu' must be a single positive number")
  if (!missing(sigma) & (length(sigma) != 1 || sigma<=0 || is.na(sigma))) stop("'sigma' must be a single positive number")
  if (!missing(p) & (length(p) != 1 || !is.finite(p) ||
                      p < 0 || p > 1)) stop("'p' must be a single number between 0 and 1")
  if (!missing(fpp) & (length(fpp) != 1 || !is.finite(fpp) ||
                        fpp < 0 || fpp > 1)) stop("'fpp' must be a single number between 0 and 1")
  if(!(fpp.method=="auto" | fpp.method=="fixed")) stop("'fpp.method' expected to be one of 'fixed' or 'auto'")
  if(!(p.method=="auto" | p.method=="fixed")) stop("'p.method' expected to be one of 'fixed' or 'auto'")
  if (!(fun=="normal" || fun=="gamma")) stop("fun needs to be either 'normal' or 'gamma'")
  if(!(N%%1==0 & N>=0)) stop("'N' expected to be an integer larger than 0")
  if(!(n%%1==0 & n>=0)) stop("'n' expected to be an integer larger than 0")
  if(!is.na(n.ind) && !(n.ind%%1==0 & n>=0)) stop("'n.ind' expected to be an integer larger than 0")
  if(!(is.numeric(trunc) & length(trunc) == 2)) stop("'trunc' expected to be a numeric vector of length 2")
  if(min(data)<trunc[1] | max(data)>trunc[2]) warning("'data' contains values outside truncation range 'trunc'")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1)) stop("'conf.level' must be a single number between 0 and 1")
}


#' Probability density function of an observed interval distribution
#'
#' Observed intervals are assumed to be sampled through observation of continuous
#' distinct arrivals in time. Two subsequently observed arrivals mark the start and end
#' of an interval. The probability that an arrival is not observed can be nonzero, leading to
#' observed intervals at integer multiples of the true interval.
#' @param data A list of intervals for which to calculate the probability density
#' @param mu The mean of the true interval distribution
#' @param sigma The standard deviation of the true interval distribution
#' @param p The probability that an arrival that marks the start or end of an interval is not observed
#' @param N The maximum number of consecutive missed arrivals to take into consideration
#' @param fun assumed distribution family of the true interval distribution, one of
#' "\code{normal}" or "\code{gamma}", corresponding
#' to the \link[stats]{Normal} and \link[stats]{GammaDist} distributions.
#' @param trunc Use a truncated probability density function with range \code{trunc}
#' @param fpp Baseline proportion of intervals distributed as a random poisson process with mean arrival rate \code{mu}
#' @param sigma.within within-subject standard deviation, only available when \code{fun} is "normal"
#' @export
#' @return This function returns a list of points describing the interval distribution
#' @details
#' \subsection{General}{
#' intervals x are assumed to follow a standard distribution (either a normal
#' or gamma distribution) with probability density function \eqn{\phi(x|\mu,\sigma)}
#' with \eqn{\mu} the mean arrival interval and \eqn{\sigma} its associated standard deviation.
#' The probability density function \eqn{\phi_{obs}} of observed arrival intervals
#' in a scenario where the probability to not observe an arrival is nonzero,
#' will be a superposition of several standard distributions, at multiples of the fundamental mean
#' arrival interval. Standard distribution \eqn{i} will correspond to those intervals where \eqn{i} arrivals have been
#' missed consecutively. If \eqn{p} equals this probability of not observing an arrival, then the
#' probability \eqn{P(i)} to miss \eqn{i} consecutive arrivals equals
#' \deqn{P(i)=p^i-p^{i+1}}{P(i)=p^i-p^{i+1}}
#' The width of standard distribution i will be broadened relative to the fundamental, according to
#' standard uncertainty propagation in the case of addition. Both in the case
#' of normal and gamma-distributed intervals (see next subsections) we may write for the observed
#' probability density function, \eqn{\phi_{obs}}:
#' \deqn{\phi_{obs}(x | \mu, \sigma, p)=\sum_{i=1}^\infty \phi_{obs}(x,i | \mu,\sigma,p)}{\phi_{obs}(x | \mu, \sigma, p)=\sum_{i=1}^\infty \phi_{obs}(x,i | \mu,\sigma,p)}
#' with
#' \deqn{\phi_{obs}(x,i | \mu, \sigma, p)= P(i-1) \phi(x | i \mu,\sqrt i \sigma)}{\phi_{obs}(x,i | \mu, \sigma, p)= P(i-1) \phi(x | i \mu,\sqrt i \sigma)}
#' In practice, this probability density function is well approximate when the infinite sum is capped at a finite integer N.
#' Be default the sum is ran up to N=5.
#'
#' }
#' \subsection{Gamma-distributed intervals}{
#' By default intervals x are assumed to follow a Gamma (\link[stats]{GammaDist}) distribution \eqn{Gamma(\mu,\sigma)}~\code{\link[stats]{dgamma}(shape=}\eqn{\mu^2/\sigma^2}\code{, scale=}\eqn{\sigma^2/\mu)}
#' with a probability density function \eqn{\phi(x)}:
#' \deqn{\phi(x|\mu,\sigma)~Gamma(\mu,\sigma)}{\phi(x|\mu,\sigma)~Gamma(\mu,\sigma)}
#' which has a mean \eqn{\mu} and standard deviation \eqn{\sigma}.
#' }
#' \subsection{Normal-distributed intervals}{
#' intervals x may also be assumed to follow a \link[stats]{Normal} distribution \eqn{N(\mu,\sigma)}~\code{\link[stats]{dnorm}(mean=}\eqn{\mu}\code{,sd=}\eqn{\sigma)},
#' with a probability density function \eqn{\phi(x)}:
#' \deqn{\phi(x|\mu,\sigma)~N(\mu,\sigma)}{\phi(x|\mu,\sigma)~N(\mu,\sigma)}
#' which also has a mean \eqn{\mu} and standard deviation \eqn{\sigma}. Because intervals
#' are by definition non-negative, the Normal distribution is always truncated at zero.
#' In the limit that \eqn{\mu>\sigma} the gamma distribution tends to the normal distribution.
#' }
#' \subsection{Within and between-subject variation}{
#' To account for witin-subject and between-subject differences in mean interval length we define
#' \eqn{\sigma_w} as within-subject standard deviation in interval length,
#' and \eqn{\sigma_b} as between-subject standard deviation in interval length,
#' with \eqn{\sigma^2 = \sigma^2_b + \sigma^2_w}.
#' In the normal limit (\eqn{\mu>\sigma}) the population pdf will be a convolution between \eqn{\phi(x|\mu,\sigma_b)}
#' and \eqn{\phi(x|\mu,\sigma_w)} equal to:
#' \deqn{\phi_{obs}(x | \mu, \sigma,\sigma_w,p)=\sum_{i=1}^\infty P(i-1) \phi(x | i \mu,\sqrt i \sigma)}{\phi_{obs}(x | \mu,\sigma,p)=\sum_{i=1}^\infty P(i-1) \phi(x | i \mu,\sqrt{i {\sigma_w}^2 + \sigma^2})}
#' }
#' @examples
#' # a low probability of not observing an arrival
#' # results in an observed PDF with primarily
#' # a single peak, with a mean and standard
#' # deviation almost identical to the true interval
#' # distribution:
#' plot(intervalpdf(mu=200,sigma=40,p=0.01),type='l',col='red')
#'
#' # a higher probability to miss an arrival
#' # results in an observed PDF with multiple
#' # peaks at integer multiples of the mean of the true
#' # interval distribution
#' plot(intervalpdf(mu=200,sigma=40,p=0.4),type='l',col='red')
intervalpdf=function(data=seq(0,1000),mu=200,sigma=40,p=0.3,N=5L,fun="gamma",trunc=c(0,Inf),fpp=0,sigma.within=NA){
  checkargs(data=data,mu=mu,sigma=sigma,p=p,N=N,fun=fun,trunc=trunc,fpp=fpp,sigma.within=sigma.within)
  if(fun=="normal"){
    funpdf=normi
    funcdf=normpi
  }
  else{
    funpdf=gammai
    funcdf=gammapi
  }
  if(!is.na(sigma.within)){
    sigmap=sigma
    sigma=sigma.within
  }
  else sigmap=sigma
  data.frame(interval=data,density=probdens(data,mu,sigma,p,N,fpp=fpp,fun=funpdf,funcdf=funcdf,trunc=trunc,sigmap=sigmap))
}

# log-likelihood of an observed interval distribution
loglik=function(data,mu,sigma,p,N=5L,fun=normi,funcdf=normpi,fpp=0,trunc=c(0,Inf)) sum(log(probdens(data,mu,sigma,p,N,fun=fun,fpp=fpp,trunc=trunc)))

# log-likelihood, with parameters to be optimised in list params
loglikfn1=function(params,data,N=5L,fun=gammai,funcdf=gammapi,fpp,p,trunc=c(0,Inf)) sum(log(eps+probdens(data,params[1],params[2],p=p,N=N,fun=fun,funcdf=funcdf,fpp=fpp,trunc=trunc)))
loglikfn2=function(params,data,N=5L,fun=gammai,funcdf=gammapi,fpp,p,trunc=c(0,Inf)) sum(log(eps+probdens(data,params[1],params[2],plogis(params[3]),N,fun=fun,funcdf=funcdf,fpp=fpp,trunc=trunc)))
loglikfn3=function(params,data,N=5L,fun=gammai,funcdf=gammapi,fpp,p,trunc=c(0,Inf)) sum(log(eps+probdens(data,params[1],params[2],p=p,N=N,fun=fun,funcdf=funcdf,trunc=trunc,fpp=plogis(params[3]))))
loglikfn4=function(params,data,N=5L,fun=gammai,funcdf=gammapi,fpp,p,trunc=c(0,Inf)) sum(log(eps+probdens(data,params[1],params[2],plogis(params[3]),N,fun=fun,funcdf=funcdf,trunc=trunc,fpp=plogis(params[4]))))
# log-likelihood of a null model without missed event probability (i.e. N=1, p=0)
logliknull1=function(params,data,N=5L,fun=gammai,funcdf=gammapi,fpp=0,trunc=c(0,Inf),...) sum(log(probdens(data,params[1],params[2],0,N,fun=fun,funcdf=funcdf,trunc=trunc,fpp=fpp)))
logliknull2=function(params,data,N=5L,fun=gammai,funcdf=gammapi,trunc=c(0,Inf),...) sum(log(probdens(data,params[1],params[2],0,N,fun=fun,funcdf=funcdf,trunc=trunc,fpp=plogis(params[3]))))

# log-likelihood fns, within/between individual variation separate:
loglikfn1WB=function(params,data,N=5L,fun=gammai,funcdf=gammapi,fpp,p,trunc=c(0,Inf),sigma.within) sum(log(eps+probdens(data,params[1],sigma.within,p=p,N=N,fun=fun,funcdf=funcdf,fpp=fpp,trunc=trunc,sigmap=params[2])))
loglikfn2WB=function(params,data,N=5L,fun=gammai,funcdf=gammapi,fpp,p,trunc=c(0,Inf),sigma.within) sum(log(eps+probdens(data,params[1],sigma.within,plogis(params[3]),N,fun=fun,funcdf=funcdf,fpp=fpp,trunc=trunc,sigmap=params[2])))
loglikfn3WB=function(params,data,N=5L,fun=gammai,funcdf=gammapi,fpp,p,trunc=c(0,Inf),sigma.within) sum(log(eps+probdens(data,params[1],sigma.within,p=p,N=N,fun=fun,funcdf=funcdf,trunc=trunc,fpp=plogis(params[3]),sigmap=params[2])))
loglikfn4WB=function(params,data,N=5L,fun=gammai,funcdf=gammapi,fpp,p,trunc=c(0,Inf),sigma.within) sum(log(eps+probdens(data,params[1],sigma.within,plogis(params[3]),N,fun=fun,funcdf=funcdf,trunc=trunc,fpp=plogis(params[4]),sigmap=params[2])))
# log-likelihood of a null model without missed event probability (i.e. N=1, p=0)
logliknull1WB=function(params,data,N=5L,fun=gammai,funcdf=gammapi,fpp=0,trunc=c(0,Inf),sigma.within,...) sum(log(probdens(data,params[1],sigma.within,0,N,fun=fun,funcdf=funcdf,trunc=trunc,fpp=fpp,sigmap=params[2])))
logliknull2WB=function(params,data,N=5L,fun=gammai,funcdf=gammapi,trunc=c(0,Inf),sigma.within,...) sum(log(probdens(data,params[1],sigma.within,0,N,fun=fun,funcdf=funcdf,trunc=trunc,fpp=plogis(params[3]),sigmap=params[2])))


#' log-likelihood of an observed interval distribution
#'
#' @param data A numeric list of intervals.
#' @param mu mean arrival interval.
#' @param sigma standard deviation of the arrival interval.
#' @param p chance to not observe an arrival.
#' @param N Maximum number of missed observations to be taken into account (default N=5).
#' @param fun Assumed distribution for the intervals, one of "\code{normal}" or "\code{gamma}", corresponding
#' to the \link[stats]{Normal} and \link[stats]{GammaDist} distributions
#' @param trunc Use a truncated probability density function with range \code{trunc}
#' @param fpp Baseline proportion of intervals distributed as a random poisson process with mean arrival interval \code{mu}
#' @details
#' Refer to \link[intRvals]{intervalpdf} for details on the functional form of
#' the probability density function of an observed interval distribution \eqn{\phi_{obs}}.
#' The log-likelihood \eqn{L} given a set of intervals {\eqn{x_j}} in \code{data} is given by
#' \deqn{L(\mu,\sigma,p)=\log \sum_j \phi_{obs}(x_j | \mu,\sigma,p)}
#' The function is provided to allow likelihood maximisation by other optimization
#' tools than the default by \link[stats]{optim}.
#' @export
#' @return returns the value of the loglikelihood
#' @examples
#' data(goosedrop)
#' loglikinterval(goosedrop$interval,mu=200,sigma=50,p=.3)
loglikinterval=function(data,mu,sigma,p,N=5L,fun="gamma",trunc=c(0,Inf),fpp=0){
  checkargs(data=data,mu=mu,sigma=sigma,p=p,N=N,fun=fun,trunc=trunc,fpp=fpp)
  if(!missing(data) & length(data)<=2) stop("no or insufficient data points")
  if(fun=="normal"){
    funpdf=normi
    funcdf=normpi
  }
  else{
    funpdf=gammai
    funcdf=gammapi
  }
  loglik(data=data,mu=mu,sigma=sigma,p=p,N=N,fun=funpdf,funcdf=funcdf,trunc=trunc,fpp=fpp)
}

#' Plot an interval histogram and fit of intRvals object
#'
#' @param x An intRvals class object
#' @param binsize Width of the histogram bins
#' @param line.col Color of the plotted curve for the model fit
#' @param line.lwd Line width of the plotted curve for the model fit
#' @param main an overall title for the plot
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param ... Additional arguments to be passed to the low level \link[graphics]{hist} plotting function
#' @export
#' @return This function returns a list with data, corresponding to the model fit
#' @examples
#' data(goosedrop)
#' dr=estinterval(goosedrop$interval)
#' plot(dr)
#' plot(dr,binsize=10,line.col='blue')
plot.intRvals=function(x,binsize=20,xlab="Interval",ylab="Density",main="Interval histogram and fit", line.col='red', line.lwd=1, ...){
  object=x
  stopifnot(inherits(object, "intRvals"))
  if(object$distribution=="normal"){
    funpdf=normi
    funcdf=normpi
  }
  else{
    funpdf=gammai
    funcdf=gammapi
  }
  hist(object$data,freq=F,breaks=seq(0,max(object$data)+binsize,binsize),xlab=xlab,ylab=ylab,main=main,...)
  if(!is.na(object$sigma.within)) plotfunc=function(x) probdens(x,object$mu,object$sigma.within,object$p,object$N,fun=funpdf,funcdf=funcdf,trunc=object$trunc,fpp=object$fpp,sigmap=object$sigma)
  else plotfunc=function(x) probdens(x,object$mu,object$sigma,object$p,object$N,fun=funpdf,funcdf=funcdf,trunc=object$trunc,fpp=object$fpp)

  curve(plotfunc,0,1500,col=line.col,lwd=line.lwd,add=T)
}

# helper function to set up a call to optim()
prepare.optim=function(mu=200,sigma=50,p=0.2,N=5L,fun="gamma",trunc=c(0,Inf),fpp=0,fpp.method="fixed",p.method="auto",sigma.within=NA){
  fpp.logit=log(fpp/(1-fpp))
  p.logit=log(p/(1-p))
  if(fun=="normal"){
    funpdf=normi
    funcdf=normpi
  }
  else{
    funpdf=gammai
    funcdf=gammapi
  }
  if(fpp.method=="fixed" && p.method=="fixed"){
    par=c(mu,sigma)
    par.null=c(mu,sigma)
    if(is.na(sigma.within)){
      L=loglikfn1
      L.null=logliknull1
    }
    else{
      L=loglikfn1WB
      L.null=logliknull1WB
    }
  }
  if(fpp.method=="fixed" && p.method=="auto"){
    par=c(mu,sigma,p.logit)
    par.null=c(mu,sigma)
    if(is.na(sigma.within)){
      L=loglikfn2
      L.null=logliknull1
    }
    else{
      L=loglikfn2WB
      L.null=logliknull1WB
    }
  }
  if(fpp.method=="auto" && p.method=="fixed"){
    par=c(mu,sigma,fpp.logit)
    par.null=c(mu,sigma,fpp.logit)
    if(is.na(sigma.within)){
      L=loglikfn3
      L.null=logliknull2
    }
    else{
      L=loglikfn3WB
      L.null=logliknull2WB
    }
  }
  if(fpp.method=="auto" && p.method=="auto"){
    par=c(mu,sigma,p.logit,fpp.logit)
    par.null=c(mu,sigma,fpp.logit)
    if(is.na(sigma.within)){
      L=loglikfn4
      L.null=logliknull2
    }
    else{
      L=loglikfn4WB
      L.null=logliknull2WB
    }
  }
  list(par=par,L=L,par.null=par.null,L.null=L.null,pdf=funpdf,cdf=funcdf,fpp.logit=fpp.logit,p.logit=p.logit)
}

# prepares optimization-specific output
prepare.output=function(opt,fpp,p,fpp.method="fixed",p.method="auto"){
  if(fpp.method=="fixed" && p.method=="fixed"){
    fpp.out=fpp
    p.out=p
    # manual adjustment of missed event probability treated as optimization, i.e. loss of one df
    if(p==0) n.opt=2 else n.opt=3
  }
  if(fpp.method=="fixed" && p.method=="auto"){
    fpp.out=fpp
    p.out=plogis(opt$par[3])
    n.opt=3
  }
  if(fpp.method=="auto" && p.method=="fixed"){
    p.out=p
    fpp.out=plogis(opt$par[3])
    # manual adjustment of missed event probability treated as optimization, i.e. loss of one df
    if(p==0) n.opt=3 else n.opt=4
  }
  if(fpp.method=="auto" && p.method=="auto"){
    p.out=plogis(opt$par[3])
    fpp.out=plogis(opt$par[4])
    n.opt=4
  }
  list(mu=opt$par[1],sigma=opt$par[2],loglik=opt$value,p=p.out,fpp=fpp.out,n.opt=n.opt)
}


#' Estimate interval model accounting for missed arrival observations
#'
#' Estimate interval mean and variance accounting for missed arrival observations,
#' by fitting the probability density function \link[intRvals]{intervalpdf} to the interval data.
#' @param data A numeric list of intervals.
#' @param mu Start value for the numeric optimization for the mean arrival interval.
#' @param sigma Start value for the numeric optimization for the standard deviation of the arrival interval.
#' @param p Start value for the numeric optimization for the probability to not observe an arrival.
#' @param N Maximum number of missed observations to be taken into account (default N=5).
#' @param fun Assumed distribution for the intervals, one of "\code{normal}" or "\code{gamma}", corresponding
#' to the \link[stats]{Normal} and \link[stats]{GammaDist} distributions
#' @param trunc Use a truncated probability density function with range \code{trunc}
#' @param fpp Baseline proportion of intervals distributed as a random poisson process with mean arrival interval \code{mu}
#' @param fpp.method A string equal to 'fixed' or 'auto'. When 'auto' \code{fpp} is optimized as a free model parameter,
#' in which case \code{fpp} is taken as start value in the optimisation
#' @param p.method A string equal to 'fixed' or 'auto'. When 'auto' \code{p} is optimized as a free model parameter,
#' in which case \code{p} is taken as start value in the optimisation
#' @param conf.level Confidence level for deviance test that checks whether model with nonzero missed event probability
#' \code{p} significantly outperforms a model without a missed event probability (\code{p=0}).
#' @param group optional vector of equal length as data, indicating the group or subject in which the interval was observed
#' @param sigma.within optional within-subject standard deviation. When equal to default 'NA', assumes
#' no additional between-subject effect, with \code{sigma.within} equal to \code{sigma}. When equal to 'auto'
#' an estimate is provided by iteratively calling \link[intRvals]{partition}
#' @param iter maximum number of iterations in numerical iteration for \code{sigma.within}
#' @param tol tolerance in the iteration, when \code{sigma.within} changes less than this value in one iteration step, the optimization is considered converged.
#' @param silent logical. When \code{TRUE} print no information to console
#' @param ... Additional arguments to be passed to \link[stats]{optim}
#' @details
#' The probability density function for observed intervals \link[intRvals]{intervalpdf}
#' is fit to \code{data} by maximization of the
#' associated log-likelihood using \link[stats]{optim}.
#'
#' Within-group variation \code{sigma.within} may be separated from the total variation \code{sigma} in an iterative fit of \link[intRvals]{intervalpdf} on the interval data.
#' In the iteration \link[intRvals]{partition} is used to (1) determine which intervals according to the fit are a fundamental interval at a confidence level \code{conf.level},
#' and (2) to partition the within-group variation from the total variation in interval length.
#'
#' Within- and between-group variation is estimated on the subset of fundamental intervals with repeated measures only.
#' As the set of fundamental interval depends on the precise value of \code{sigma.within}, the fit of \link[intRvals]{intervalpdf} and the subsequent estimation of
#' \code{sigma.within} using \link[intRvals]{partition} is iterated until both converge to a stable solution. Parameters \code{tol}
#' and \code{iter} set the threshold for convergence and the maximum number of iterations.
#'
#' We note that an exponential interval model can be fitted by setting \code{fpp=1} and \code{fpp.method=fixed}.
#' @export
#' @return This function returns an object of class \code{intRvals}, which is a list containing the following:
#' \describe{
#'   \item{\code{data}}{the interval data}
#'   \item{\code{mu}}{the modelled mean interval}
#'   \item{\code{mu.se}}{the modelled mean interval standard error}
#'   \item{\code{sigma}}{the modelled interval standard deviation}
#'   \item{\code{p}}{the modelled probability to not observe an arrival}
#'   \item{\code{fpp}}{the modelled fraction of arrivals following a random poisson process, see \link[intRvals]{intervalpdf}}
#'   \item{\code{N}}{the highest number of consecutive missed arrivals taken into account, see \link[intRvals]{intervalpdf}}
#'   \item{\code{convergence}}{convergence field of \link[stats]{optim}}
#'   \item{\code{counts}}{counts field of \link[stats]{optim}}
#'   \item{\code{loglik}}{vector of length 2, with first element the log-likelihood of the fitted model, and second element the log-likelihood of the model without a missed event probability (i.e. \code{p}=0)}
#'   \item{\code{df.residual}}{degrees of freedom, a 2-vector (1, number of intervals - \code{n.param})}
#'   \item{\code{n.param}}{number of optimized model parameters}
#'   \item{\code{p.chisq}}{p value for a likelihood-ratio test of a model including a miss probability relative against a model without a miss probability}
#'   \item{\code{distribution}}{assumed interval distribution, one of 'gamma' or 'normal'}
#'   \item{\code{trunc}}{interval range over which the interval pdf was truncated and normalized}
#'   \item{\code{fpp.method}}{A string equal to 'fixed' or 'auto'. When 'auto' \code{fpp} has been optimized as a free model parameter}
#'   \item{\code{p.method}}{A string equal to 'fixed' or 'auto'. When 'auto' \code{p} has been optimized as a free model parameter}
#' }
#' @examples
#' data(goosedrop)
#' # calculate mean and standard deviation of arrival intervals, accounting for missed observations:
#' dr=estinterval(goosedrop$interval)
#' # plot some summary information
#' summary(dr)
#' # plot a histogram of the intervals and fit:
#' plot(dr)
#' # test whether the mean arrival interval is greater than 200 seconds:
#' ttest(dr,mu=200,alternative="greater")
#'
#' # let's estimate mean and variance of dropping intervals by site
#' # (schiermonnikoog vs terschelling) for time period 5.
#' # first prepare the two datasets:
#' set1=goosedrop[goosedrop$site=="schiermonnikoog" & goosedrop$period == 5,]
#' set2=goosedrop[goosedrop$site=="terschelling"  & goosedrop$period == 5,]
#' # allowing a fraction of intervals to be distributed randomly (fpp='auto')
#' dr1=estinterval(set1$interval,fpp.method='auto')
#' dr2=estinterval(set2$interval,fpp.method='auto')
#' # plot the fits:
#' plot(dr1,xlim=c(0,1000))
#' plot(dr2,xlim=c(0,1000))
#' # mean dropping interval are not significantly different
#' # at the two sites (on a 0.95 confidence level):
#' ttest(dr1,dr2)
#' # now compare this test with a t-test not accounting for unobserved intervals:
#' t.test(set1$interval,set2$interval)
#' # not accounting for missed observations leads to a (spurious)
#' # larger difference in means, which also increases
#' # the apparent statistical significance of the difference between means
estinterval=function(data,mu=median(data),sigma=sd(data)/2,p=0.2,N=5L,fun="gamma",trunc=c(0,Inf),fpp=(if(fpp.method=="fixed") 0 else 0.1),fpp.method="auto",p.method="auto", conf.level = 0.9, group = NA, sigma.within=NA, iter=10, tol=0.001,silent=F, ...){
  checkargs(data=data,mu=mu,sigma=sigma,p=p,N=N,fun=fun,trunc=trunc,fpp=fpp,fpp.method=fpp.method,p.method=p.method,conf.level=conf.level,sigma.within=sigma.within)

  #discard data outside truncation interval
  data=data[which(data>trunc[1] & data<trunc[2])]
  call=sys.call()

  if(is.na(sigma.within) || sigma.within=='auto'){
    est=estintervalHelper(data,mu,sigma,p,N,fun,trunc,fpp,fpp.method,p.method,conf.level,group,NA,call)
    est$p.within=NA
  }
  if(!is.na(sigma.within) && sigma.within=='auto'){
    iter=0
    found=F
    # initialize sigma.within:
    est$sigma.within=sigma/2
    while(!found && iter<10){
      sigmaw.old=est$sigma.within
      part=partition(est,silent=silent,conf.level=conf.level)
      if(is.na(part$sigma.within)) break;
      sigmaw.new=part$sigma.within
      est=estintervalHelper(data,mu,sigma,p,N,fun,trunc,fpp,fpp.method,p.method,conf.level,group,sigmaw.new,call)
      iter=iter+1
      if(abs(sigmaw.old-sigmaw.new)<tol) found=T
    }
    if(found) est$p.within=part$p.within else est$p.within=NA
    if(!silent & found) cat("sigma.within converged in",iter,"iterations...\n")
    if(!silent & !found & !is.na(est$p.within)) cat("no convergence of sigma.within in",iter,"iterations...\n")
    if(!silent & !found & is.na(est$p.within)) cat("no convergence of sigma.within, insufficient fundamental intervals...\n")
  }
  if(is.numeric(sigma.within)){
    est=estintervalHelper(data,mu,sigma,p,N,fun,trunc,fpp,fpp.method,p.method,conf.level,group,sigma.within,call)
    return(est)
    part=partition(est,silent=silent,conf.level=conf.level)
    est$p.within=part$p.within
  }

  est
}

estintervalHelper=function(data,mu=median(data),sigma=sd(data)/2,p=0.2,N=5L,fun="gamma",trunc=c(0,Inf),fpp=(if(fpp.method=="fixed") 0 else 0.1),fpp.method="fixed",p.method="auto", conf.level = 0.95, group = NA, sigma.within=NA,call,...){

  # set up the optimization
  optpars=prepare.optim(mu=mu,sigma=sigma,p=p,N=N,fun=fun,trunc=trunc,fpp=fpp,fpp.method=fpp.method,p.method=p.method,sigma.within=sigma.within)
  # run optimization for model and null model
  if(!is.na(sigma.within)){
    opt=optim(par=optpars$par,optpars$L,data=data,p=p,fpp=fpp,N=N,fun=optpars$pdf,funcdf=optpars$cdf,trunc=trunc,sigma.within=sigma.within,control=list(fnscale=-1),method="BFGS")
    optnull=optim(par=optpars$par.null,fn=optpars$L.null,data=data,N=1,fun=optpars$pdf,funcdf=optpars$cdf,trunc=trunc,fpp=fpp,sigma.within=sigma.within,control=list(fnscale=-1),method="BFGS")
  }
  else{
    opt=optim(par=optpars$par,optpars$L,data=data,p=p,fpp=fpp,N=N,fun=optpars$pdf,funcdf=optpars$cdf,trunc=trunc,control=list(fnscale=-1),method="BFGS")
    optnull=optim(par=optpars$par.null,fn=optpars$L.null,data=data,N=1,fun=optpars$pdf,funcdf=optpars$cdf,trunc=trunc,fpp=fpp,control=list(fnscale=-1),method="BFGS")
  }
  # prepare optimization-specific output
  out.prep=prepare.output(opt,fpp,p,fpp.method,p.method)
  # prepare output
  out=list(call=call,data=data,mu=opt$par[1],sigma=opt$par[2],p=out.prep$p,fpp=out.prep$fpp,N=N,convergence=opt$convergence,counts=opt$counts,loglik=c(opt$value,optnull$value),df.residual=c(1,length(data)-out.prep$n.opt),n.param=out.prep$n.opt,p.chisq=NA,distribution=fun,trunc=trunc,fpp.method=fpp.method,p.method=p.method,group=group,sigma.within=sigma.within,p.within=NA)
  # add standard error of mu
  out$mu.se=out$mu/sqrt(out$df.residual[2])

  if(!is.na(sigma.within) && out$sigma<sigma.within) warning("'sigma.within' too large, fit did not converge ...")
  if(out$p>0.99) warning("'p' too large, fit did not converge ...")

  # test whether model model outperforms null model
  out$p.chisq=pchisq(2*(out$loglik[1]-out$loglik[2]), df=out$df.residual[1], lower.tail=FALSE)
  if(out$p.chisq>1-conf.level){
    warning("no support for a model including a miss probability relative to a model without a miss probability at specified confidence level.\n\nCheck convergence using different starting values (arguments: mu,sigma and p).\nConsider including a random poisson process background (argument: fpp),\nor changing the interval distribution (argument: fun).")
  }

  class(out)="intRvals"
  out
}

fundamentalProb = function(x,mu,sigma,p,fpp,N,fun){
  # get the probabilities of each component of the PDF
  probs=(1-fpp)*pdfcomponents(x,mu,sigma,p,N,fun)
  # get the probability that x is part of fpp component
  probs.fpp=fpp*pdfcomponents(x,mu,mu,p,N,fun=gammai)
  probs=c(probs,probs.fpp)
  # normalize
  probs=(probs/sum(probs))
  # return relative probably belonging to fundamental
  probs[1]
}

#' Estimate which intervals are fundamental
#'
#' Estimates which intervals in a dataset are fundamental intervals, i.e. an
#' interval not containing a missed arrival observation
#' @param x object inheriting from class \code{intRvals}, usually a result of a call to \code{\link[intRvals]{estinterval}}
#' @param conf.level confidence level for identifying intervals as fundamental
#' @return logical atomic vector of the same length as \code{x$data}
#' @details
#' This functions thus determines for each interval \code{x$data} whether it has a probabiliy > \code{conf.level} to be
#' a fundamental interval, given the model fit generated by \link[intRvals]{estinterval} for object \code{x}.
#'
#' The fit of an \code{intRvals} object gives the decomposition of the likelihood of an interval observation
#' into partial likelihoods \eqn{\phi_{obs}(x,i | \mu, \sigma, p)} (see \link[intRvals]{intervalpdf}).
#' If the amplitude of the partial likelihood with i=0 (i.e. the likelihood component without missed observations)
#' is at least a proportion \code{conf.level} of the sum of all terms i=0..N,
#' an interval is considered to be fundamental (not containing a missed event observation).
#' @export
fundamental = function(x, conf.level=0.9){
  stopifnot(inherits(x,"intRvals"))
  if(x$distribution=="gamma") fun = gammai
  else fun = normi
  output=sapply(x$data,function(ival) fundamentalProb(ival,x$mu,x$sigma,x$p,x$fpp,x$N,fun)>conf.level)
  output
}


#' Simulate a set of observed intervals
#'
#' @param n Number of simulated interval observations.
#' @param mu Mean arrival interval.
#' @param sigma Standard deviation of the arrival interval.
#' @param p Probability to not observe an arrival.
#' @param fun Assumed distribution for the intervals, one of "\code{normal}" or "\code{gamma}", corresponding
#' to the \link[stats]{Normal} and \link[stats]{GammaDist} distributions
#' @param trunc Observational range of intervals (intervals outside this range won't be observed)
#' @param fpp Baseline proportion of intervals distributed as a random poisson process with mean arrival interval \code{mu}
#' @param n.ind Number of intervals per group. Ignored without a numeric value for \code{sigma.within}.
#' @param sigma.within The within-group standard-deviation. When a numeric value is given for \code{sigma.within}, \code{sigma} denotes the total (within+between subject) standard deviation
#' @details
#' Simulates the observations process of arrival intervals.
#'
#' The default is to not differentiate between
#' within- and between-group variance.
#'
#' If both \code{n.ind} and \code{sigma.within} have numeric values, intervals are simulated
#' with separate within-group variation (\code{sigma.within}) and between-group variation,
#' for groups of size \code{n.ind}. Intervals belonging to the same group have:
#' \itemize{
#'  \item a within-group mean interval length that has been randomly drawn from a distribution with mean \code{mu}
#'  and between-group standard deviation \eqn{\sqrt{sigma^2 - sigma.within^2}}
#'  \item a within-group standard deviation in interval length equal to \code{sigma.within}
#' }
#'
#' @export
#' @return This function returns a dataframe containing the following:
#' \describe{
#'   \item{\code{interval}}{the simulated interval data}
#'   \item{\code{group_id}}{a group identifier}
#' }
#' @examples
#' # simulate observed intervals:
#' intervals=intervalsim(n=50,mu=200,sigma=40,trunc=c(0,600),fpp=0.1)
#' # check whether we retrieve the simulation parameters:
#' estinterval(goosedrop$interval)
intervalsim=function(n=500,mu=200,sigma=40,p=0.3,fun="gamma",trunc=c(0,600),fpp=0,n.ind=NA,sigma.within=NA){
  checkargs(n=n,mu=mu,sigma=sigma,p=p,fun=fun,trunc=trunc,fpp=fpp,sigma.within=sigma.within,n.ind=n.ind)
  if(!is.na(sigma.within)){
    sigma.between=sqrt(sigma^2-sigma.within^2)
    sigma=sigma.within
  } else sigma.between=NA
  if(fun=="normal") pdfsample = function(x,mu,sigma) rnorm(x,mu,sigma)
  else pdfsample = function(x,mu,sigma) rgamma(x,shape=mu^2/(sigma^2),scale=(sigma^2)/mu)
  siminterval = function(mu) {
    found=F
    while(!found){
      ival=pdfsample(1,mu,sigma)
      if(fpp>runif(1)){
        ival = rexp(1,1/mu)
      }
      while(p>runif(1)){
        if(fpp>runif(1)){
          ival = ival + rexp(1,1/mu)
        }
        else{
          ival = ival + pdfsample(1,mu,sigma)
        }
      }
      if(ival<trunc[2]) found=T
    }
    return(ival)
  }
  if (is.na(n.ind) || is.na(sigma.between)){
    data = data.frame(interval=replicate(n,siminterval(mu)),group_id=1)
  }
  else{
    mu.ind=pdfsample(1,mu,sigma.between)
    data=data.frame(interval=replicate(n.ind,siminterval(mu.ind)),group_id=rep(1,n.ind))
    group=1
    while(group*n.ind<n){
      mu.ind=pdfsample(1,mu,sigma.between)
      repdata=replicate(n.ind,siminterval(mu.ind))
      data=rbind(data,data.frame(interval=repdata,group_id=rep(group,n.ind)))
      group=group+1
    }
  }
  data[1:n,]
}

#' Conversion of interval estimates to rates
#'
#' @param data An object of class \code{intRvals}, usually a result of a call to \link[intRvals]{estinterval}
#' @param minint the minimum interval value from which numerical integrations converting to rates are started
#' @param maxint the maximum interval value up to which numerical integrations converting to rates are continued
#' @param digits the number of digits for printing to screen
#' @param method A string equal to 'exact' or 'taylor'. When 'exact' exact formula or numeric integration
#' is used. When 'taylor' a Taylor approximation is used as in standard propagation of uncertainty in the case of division.
#' @details
#' \subsection{Gamma-distributed intervals}{
#' When inter-arrival times (intervals) follow a gamma distribution with mean \eqn{\mu} and
#' standard deviation \eqn{\sigma}, i.e. follow the probability density function
#' \code{\link[stats]{GammaDist}(shape=}\eqn{\alpha=\mu^2/\sigma^2}\code{, scale=}\eqn{\beta=\sigma^2/\mu)},
#'  then the associated distribution of rates is given by an inverse gamma distribution
#'  with shape parameter \eqn{\alpha} and scale parameter \eqn{1/\beta}.
#'
#' The mean of this inverse gamma distribution is given by the formula
#' \deqn{\mu_{rate}=\mu/(\mu^2 - \sigma^2)}
#' provided that \eqn{\alpha > 1}, i.e. \eqn{\mu > \sigma}.
#'
#' The variance of this inverse gamma distribution is given by the formula
#' \deqn{\sigma^2_{rate}=\mu^2\sigma^2/((\mu^2 - \sigma^2)(\mu^2 - 2\sigma^2)^2}
#' provided that \eqn{\alpha > 2}, i.e. \eqn{\mu > sqrt(2) * \sigma}.
#'
#' Values \eqn{\mu} and \eqn{\sigma} are estimated on the interval data, and
#' above formula are used to calculate the estimated mean and variance of the arrival rate.
#'
#' If these formula cannot be used (because the provisions on the value
#' of \eqn{\alpha} are not met), numerical integration is used instead,
#' analagous to the procedure for normal-distributed intervals, see below.
#' }
#' \subsection{Normal-distributed intervals}{
#' When inter-arrival times (intervals) \eqn{x} follow a normal distribution with mean \eqn{\mu} and
#' standard deviation \eqn{\sigma}, i.e. follow the probability density function
#' \code{\link[stats]{Normal}(mean=}\eqn{\mu}\code{, sd=}\eqn{\sigma)},
#'  then the mean rate (\eqn{\mu_{rate}}) can be calculated numerically by:
#'  \deqn{\mu_{rate}=\int_0^\infty (1/x) * \phi(x | \mu,\sigma)}
#'  and the variance of the rate (\eqn{\sigma^2_{rate}}) by:
#'  \deqn{\sigma^2_{rate}=\int_0^\infty (1/x^2) * \phi(x | \mu,\sigma) -\mu_{rate}^2}
#' This approximation is only valid for distributions that have a negligable
#' density near \eqn{x=0}, such that the distribution can be effectively
#'  truncated before \eqn{x} approaches zero, where the integral is not defined.
#' For interval data with intervals \eqn{x}
#' near zero, use of a gamma distribution is recommended instead.
#' }
#' @export
#' @return The function \code{interval2rate} computes and returns a named vector with the rate mean and standard deviation
#' @examples
#' data(goosedrop)
#' dr=estinterval(goosedrop$interval)
#' interval2rate(dr)
interval2rate=function(data,minint=data$mu/100,maxint=data$mu+3*data$sigma,digits = max(3L, getOption("digits") - 3L), method="exact"){
  if (!(method=="exact" || method=="taylor")) stop("method needs to be either 'exact' or 'taylor'")
  stopifnot(inherits(data, "intRvals"))
  if(method=="exact"){
    if(data$distribution=="normal"){
      cat(paste("Numerically calculating rate mean and standard deviation\n  truncating normal distribution of intervals over range",format(signif(minint,digits)),"to",format(signif(maxint,digits)),"\n"))
      intinv=integrate(function(y) (1/y)*dnorm(y,mean=data$mu,sd=data$sigma),minint,maxint)
      intinvcubed=integrate(function(y) (1/y^2)*dnorm(y,mean=data$mu,sd=data$sigma),minint,maxint)
      rate.mean=intinv$value
      rate.stdev=sqrt(intinvcubed$value-intinv$value^2)

    }
    if(data$distribution=="gamma"){
      if(data$mu>sqrt(2)*data$sigma){
        rate.mean=data$mu/(data$mu^2-data$sigma^2)
        rate.stdev=sqrt(data$mu^2*data$sigma^2/((data$mu^2-2*data$sigma^2)*(data$mu^2-data$sigma^2)^2))
      }
      else if(data$mu>data$sigma){
        rate.mean=data$mu/(data$mu^2-data$sigma^2)
        cat(paste("Numerically calculating rate standard deviation\n  truncating gamma distribution of intervals over range",format(signif(minint,digits)),"to",format(signif(maxint,digits)),"\n"))
        intinvcubed=integrate(function(y) (1/y^2)*dgamma(y,shape=data$mu^2/(data$sigma^2),scale=(data$sigma^2)/data$mu),minint,maxint)
        rate.stdev=sqrt(intinvcubed$value-rate.mean^2)
      }
      else{
        cat(paste("Numerically calculating rate mean and standard deviation\n  truncating gamma distribution of intervals over range",format(signif(minint,digits)),"to",format(signif(maxint,digits)),"\n"))
        intinv=integrate(function(y) (1/y)*dgamma(y,shape=data$mu^2/(data$sigma^2),scale=(data$sigma^2)/data$mu),minint,maxint)
        intinvcubed=integrate(function(y) (1/y^2)*dgamma(y,shape=data$mu^2/(data$sigma^2),scale=(data$sigma^2)/data$mu),minint,maxint)
        rate.mean=intinv$value
        rate.stdev=sqrt(intinvcubed$value-intinv$value^2)
      }
    }
  }
  if(method=="taylor"){
    rate.mean=1/data$mu
    rate.stdev=data$sigma/(data$mu)^2
  }
  output=c(mean=rate.mean,stdev=rate.stdev)
  output
}

foldHelper=function(x,mu,sigma,p,fpp,N=5L,fun=normi,take.sample=F){
  components=c((1-fpp)*pdfcomponents(x,mu,sigma,p,N,fun),fpp*pdfcomponents(x,mu,mu,p,N,fun=gammai))
  if(take.sample) fold=sample(1:(2*N),size=1,prob=components)
  else fold=which.max(components)
  if(fold>N) return(NA) # in this case it's an fpp component
  else return(mu+(x-fold*mu)/sqrt(fold))
}
foldInterval=function(x,mu,sigma,p,fpp,N=5L,fun=normi,take.sample=F){
  if(length(x)>1) return(sapply(x,FUN=foldHelper,mu=mu,sigma=sigma,p=p,fpp=fpp,N=N,fun=fun,take.sample=take.sample))
  else return(foldHelper(x,mu,sigma,p,fpp,N,fun,take.sample))
}

#' Folds observed arrival intervals with missed observations back to their most likely fundamental interval
#' @title Folds observed arrival intervals to a fundamental interval
#' @param object an object of class \code{intRvals}, usually a result of a call to \link[intRvals]{estinterval}
#' @param take.sample when \code{TRUE} the number of folds of the fundamental interval is sampled randomly, taking into account the probability weight of each possibility. When \code{FALSE} the fold with the highest probability weight is taken.
#' @param sigma.within (optional) numeric value with an assumed within-group/subject standard deviation, or '\code{auto}' to estimate it automatically using \link[intRvals]{partition}.
#' @param silent logical, if \code{TRUE} print no text to console
#' @export
#' @details
#' Arrival intervals containing missed observations are folded to their most likely
#' fundamental interval according to a fit of the distribution of intervals by \link[intRvals]{estinterval}.
#'
#' There is inherent uncertainty on how many missed arrival events an observed interval contains, and therefore to
#' which fundamental interval it should be folded. Intervals folded to the fundamental
#' can therefore introduce extra unexplained variance.
#'
#' The default is to fold intervals to the
#' fundamental with the highest probability weight (\code{take.sample = F}). Alternatively, randomly sampled intervals
#' can be generated, that take into account the probability weights of each possible fold (\code{take.sample = T}).
#'
#' Intervals \code{x} are transformed to their fundamental interval according to
#' \deqn{\mu+(x-i*\mu)/\sqrt i}{\mu+(x-i*\mu)/\sqrt i}
#' with \code{i-1} the estimated number of missed observations within the interval. This transformation scales appropriately
#' with the expected broadening of the standard distributions \eqn{\phi(x | i \mu,\sqrt i \sigma)} with \code{i} in \link[intRvals]{intervalpdf}.
#'
#' When no \code{sigma.within} is provided, \eqn{\mu} equals the mean arrival rate, estimated by \link[intRvals]{estinterval}.
#'
#' When \code{sigma.within} is '\code{auto}', \code{sigma.within} is estimated using \link[intRvals]{partition}.
#'
#' When \code{sigma.within} is a user-specified numeric value or '\code{auto}', \eqn{\mu} is estimated for each group (
#' as specified in the group argument of \link[intRvals]{estinterval}),
#' by maximizing the log-likelihood of \link[intRvals]{intervalpdf}, with its \code{data} argument equals to the intervals of the group,
#' its \code{sigma} argument equal to \code{sigma.within}, and its remaining arguments taken from \code{object}.
#'
#' Intervals assigned to the \code{fpp} component (see \link[intRvals]{estinterval}) are not
#' folded, and return as \code{NA} values.
#'
#' @return numeric vector with intervals folded into the fundamental interval
#' @examples
#' dr=estinterval(goosedrop$interval,group=goosedrop$bout_id)
#' # fold assuming no within-group variation:
#' interval.fundamental=fold(dr)
#' # test whether there is evidence for within-group variation:
#' partition(dr)$`p<alpha`   #> TRUE
#' # there is evidence, therefore better to fold
#' # while accounting for within-group variation:
#' interval.fundamental=fold(dr,sigma.within='auto')
fold=function(object, take.sample=F, sigma.within=NA,silent=F){
  stopifnot(inherits(object, "intRvals"))
  if(!is.na(sigma.within)){
    if(length(object$group)==1 && is.na(object$group)) stop("no groups found in object")
    if(!(length(object$group)==length(object$data))) stop("'group' and 'data' are of unequal length")
    if(sigma.within=='auto'){
      sigma.within=partition(object,silent=silent)$sigma.within
      if(!silent) cat("   estimate sigma.within:",sigma.within,"\n\n")
    }
  }
  # set up distribution family:
  if(object$distribution=="normal"){
    funpdf=normi
    funcdf=normpi
  }
  else{
    funpdf=gammai
    funcdf=gammapi
  }
  # fold intervals:
  if(is.na(sigma.within)){
    # case sigma.within >> sigma.between
    output=foldInterval(object$data,object$mu,object$sigma,object$p,object$fpp,object$N,funpdf,take.sample=take.sample)
  }
  else{
    # case sigma.within same order sigma.between
    output=c()
    for(group in unique(object$group)) {
      intervals.group=object$data[object$group==group]
      foptim=function(mu) loglik(intervals.group,mu,sigma.within,object$p,object$N,funpdf,funcdf,0,object$trunc)
      # find the most likely within-group mean
      mu=optim(object$mu,foptim,method="Brent",lower=max(c(0,object$mu-3*abs(object$sigma))),upper=object$mu+3*abs(object$sigma),control=list(fnscale=-1))
      # fold the interval, using the optimized within-group mean
      output=c(output,foldInterval(intervals.group,mu$par,sigma.within,object$p,object$fpp,object$N,funpdf,take.sample=take.sample))
    }
  }
  output
}


#' summary method for class \code{intRvals}
#'
#' @param object An object of class \code{intRvals}, usually a result of a call to \link[intRvals]{estinterval}
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return The function \code{summary.intRvals} computes and returns a list of summary statistics
#' \describe{
#'   \item{\code{data}}{the interval data}
#'   \item{\code{mu}}{the modelled mean interval}
#'   \item{\code{mu.se}}{the modelled mean interval standard error}
#'   \item{\code{sigma}}{the modelled interval standard deviation}
#'   \item{\code{p}}{the modelled probability to not observe an arrival}
#'   \item{\code{fpp}}{the modelled fraction of arrivals following a random poisson process, see \link[intRvals]{intervalpdf}}
#'   \item{\code{N}}{the highest number of consecutive missed arrivals taken into account, see \link[intRvals]{intervalpdf}}
#'   \item{\code{convergence}}{convergence field of \link[stats]{optim}}
#'   \item{\code{counts}}{counts field of \link[stats]{optim}}
#'   \item{\code{loglik}}{vector of length 2, with first element the log-likelihood of the fitted model, and second element the log-likelihood of the model without a missed event probability (i.e. \code{p}=0)}
#'   \item{\code{df.residual}}{degrees of freedom, a 2-vector (1, number of intervals - \code{n.param})}
#'   \item{\code{n.param}}{number of optimized model parameters}
#'   \item{\code{distribution}}{assumed interval distribution, one of 'gamma' or 'normal'}
#'   \item{\code{trunc}}{interval range over which the interval pdf was truncated and normalized}
#'   \item{\code{fpp.method}}{A string equal to 'fixed' or 'auto'. When 'auto' fpp has been optimized as a free model parameter. When 'fixed' the model is fitted with a fixed value set by parameter \code{fpp}}
#'   \item{\code{deviance}}{deviance between the fitted model and a model without a missed event probability (i.e. \code{p}=0)}
#'   \item{\code{p.value}}{numeric vector with two elements. First element contains the p.value for a
#'   likelihood ratio (deviance) test between the fitted model and a model without a missed event probability (i.e. \code{p}=0).
#'   Second element contains the p.value for a likelihood ratio (deviance) test between the fitted model and a saturated null model.}
#' }
#' @examples
#' data(goosedrop)
#' dr=estinterval(goosedrop$interval)
#' summary(dr)
summary.intRvals=function(object, ...){
  stopifnot(inherits(object, "intRvals"))
  xx=object
  xx$deviance=c(2*(xx$loglik[1]-xx$loglik[2]),abs(2*xx$loglik[1]))
  xx$p.value=c(pchisq(xx$deviance[1], df=xx$df.residual[1], lower.tail=FALSE),pchisq(xx$deviance[2], df=xx$df.residual[2], lower.tail=FALSE))
  class(xx)="summary.intRvals"
  xx
}

#' print method for class \code{intRvals}
#'
#' @param x An object of class \code{intRvals}, usually a result of a call to \link[intRvals]{estinterval}
#' @keywords internal
#' @export
print.intRvals=function(x,digits = max(3L, getOption("digits") - 3L), ...){
  stopifnot(inherits(x, "intRvals"))
  cat("Analysis of arrival interval data with missed arrival observations\n\n")
  cat("          number of intervals: ",format(signif(length(x$data),digits)),"\n\n")
  cat("        mean arrival interval: ",format(signif(x$mu,digits))," ( s.e.",format(signif(x$mu.se,max(1,digits-1))),")\n")
  cat("           standard deviation: ",format(signif(x$sigma,digits)),"\n")
  cat("              fraction missed: ",format(signif(x$p)),"\n")
  if(x$fpp.method=="auto"){
  cat("fraction random poisson (fpp): ",format(signif(x$fpp)),"\n")
  }
}

#' print method for class \code{summary.intRvals}
#'
#' @param x An object of class \code{symmary.intRvals}, usually a result of a call to \link[intRvals]{summary.intRvals}
#' @keywords internal
#' @export
print.summary.intRvals=function(x,digits = max(3L, getOption("digits") - 3L), ...){
  stopifnot(inherits(x, "summary.intRvals"))
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("             mean arrival interval: ",format(signif(x$mu,digits)),"\n")
  cat("                standard deviation: ",format(signif(x$sigma,digits)),"\n")
  cat("   arrival observation probability: ",format(signif(1-x$p))," (1-miss probability)\n")
  cat(" baseline fraction Poisson process: ",format(signif(x$fpp))," (fpp)\n\n")

  cat("  fitted interval distribution: ",x$distribution,"\n")
  cat("           number of intervals: ",format(signif(length(x$data),digits)),"\n")
  cat("                Log-likelihood: ",format(signif(x$loglik[1],digits)),"\n")
  cat("             Residual deviance: ",format(signif(x$deviance[1]))," on ", x$df.residual[1], "degrees of freedom\n")
  cat("Likelihood ratio test, p-value: ",format(signif(x$p.value[1],digits)), "(against null model without miss probability)\n")
  cat("             Residual deviance: ",format(signif(x$deviance[2]))," on ", x$df.residual[2], "degrees of freedom\n")
  cat("Likelihood ratio test, p-value: ",format(signif(x$p.value[2],digits)), "(against saturated null model)\n")
}

#' Performs one and two sample t-tests on objects of class \code{intRvals}
#' @title Student's t-test to compare two means of objects of class \code{intRvals}
#' @param x an object of class \code{intRvals}, usually a result of a call to \link[intRvals]{estinterval}
#' @param y an (optional) object of class \code{intRvals}, usually a result of a call to \link[intRvals]{estinterval}
#' @param alternative a character string specifying the alternative hypothesis, must be one of "\code{two.sided}" (default), "\code{greater}" or "\code{less}". You can specify just the initial letter.
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test).
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param conf.level confidence level of the interval
#' @details \code{alternative = "greater"} is the alternative that \code{x} has a larger mean than \code{y}.
#' @details If the input data are effectively constant (compared to the larger of the two means) an error is generated.
#' @export
#' @return A list with class "\code{htest}" containing the same components as in \link[stats]{t.test}
#' @examples
#' data(goosedrop)
#' dr=estinterval(goosedrop$interval)
#' # perform a one-sample t-test
#' ttest(dr,mu=200)  #> FALSE, true mean not equal to 200
#' # two sample t-test
#' data.beforeMay=goosedrop[goosedrop$date<as.POSIXct('2013-05-01'),]
#' data.afterMay=goosedrop[goosedrop$date>as.POSIXct('2013-05-01'),]
#' dr.beforeMay=estinterval(data.beforeMay$interval)
#' dr.afterMay=estinterval(data.afterMay$interval)
#' ttest(dr.beforeMay,dr.afterMay)
#'
ttest = function (x, y = NULL, alternative = c("two.sided", "less", "greater"),
          mu = 0, var.equal = FALSE, conf.level = 0.95)
{
  stopifnot(inherits(x, "intRvals"))
  if (!is.null(y)){
    stopifnot(inherits(y, "intRvals"))
    if(length(x$data)!=length(y$data)){
      samedata=F
    }else{
      samedata=length(which(!(sort(x$data) == sort(y$data))))==0
    }
    if(samedata) stop("intRvals objects estimated on the same dataset")
  }
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    yok <- !is.na(y)
    xok <- !is.na(x)
    y <- y[yok]
  }
  else {
    dname <- deparse(substitute(x))
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  nx <- (1-x$fpp)*x$df.residual[2]
  mx <- x$mu
  vx <- (x$sigma)^2
  if (is.null(y)) {
    if (nx < 2)
      stop("not enough 'x' observations")
    df <- nx - 1
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx))
      stop("data are essentially constant")
    tstat <- (mx - mu)/stderr
    method <- "One Sample t-test"
    estimate <- setNames(mx, "mean of x")
  }
  else {
    ny <- (1-y$fpp)*y$df.residual[2]
    if (nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if (ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if (var.equal && nx + ny < 3)
      stop("not enough observations")
    my <- y$mu
    vy <- (y$sigma)^2
    method <- paste(if (!var.equal)
      "Welch", "Two Sample t-test")
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")
    if (var.equal) {
      df <- nx + ny - 2
      v <- 0
      if (nx > 1)
        v <- v + (nx - 1) * vx
      if (ny > 1)
        v <- v + (ny - 1) * vy
      v <- v/df
      stderr <- sqrt(v * (1/nx + 1/ny))
    }
    else {
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
      df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny -
                                                        1))
    }
    if (stderr < 10 * .Machine$double.eps * max(abs(mx),
                                                abs(my)))
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/stderr
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    cint <- c(-Inf, tstat + qt(conf.level, df))
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
  }
  else {
    pval <- 2 * pt(-abs(tstat), df)
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- if (!is.null(y))
    "difference in means"
  else "mean"
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               conf.int = cint, estimate = estimate, null.value = mu,
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}

#' Performs an F test to compare the variances of objects of class \code{intRvals}
#' @title F Test to compare two variances of objects of class \code{intRvals}
#' @param x an object of class \code{intRvals}, usually a result of a call to \link[intRvals]{estinterval}
#' @param y an (optional) object of class \code{intRvals}, usually a result of a call to \link[intRvals]{estinterval}
#' @param alternative a character string specifying the alternative hypothesis, must be one of "\code{two.sided}" (default), "\code{greater}" or "\code{less}". You can specify just the initial letter.
#' @param ratio the hypothesized ratio of the population variances of \code{x} and \code{y}.
#' @param conf.level confidence level for the returned confidence interval
#' @details The null hypothesis is that the ratio of the variances of the
#' data to which the models \code{x} and \code{y} were fitted, is equal to ratio.
#' @export
#' @return A list with class "\code{htest}" containing the same components as in \link[stats]{var.test}
#' @examples
#' data(goosedrop)
#' dr=estinterval(goosedrop$interval)
#' # split the interval data into two periods
#' data.beforeMay=goosedrop[goosedrop$date<as.POSIXct('2013-05-01'),]
#' data.afterMay=goosedrop[goosedrop$date>as.POSIXct('2013-05-01'),]
#' dr.beforeMay=estinterval(data.beforeMay$interval)
#' dr.afterMay=estinterval(data.afterMay$interval)
#' # perform an F test
#' vartest(dr.beforeMay,dr.afterMay)
#'
vartest=function (x, y, ratio = 1, alternative = c("two.sided", "less",
                                           "greater"), conf.level = 0.95)
{
  stopifnot(inherits(x, "intRvals") && inherits(y, "intRvals"))
  if(length(x$data)!=length(y$data)){
    samedata=F
  }else{
    samedata=length(which(!(sort(x$data) == sort(y$data))))==0
  }
  if(samedata) stop("intRvals objects estimated on the same dataset")
  if (!((length(ratio) == 1L) && is.finite(ratio) && (ratio >
                                                      0)))
    stop("'ratio' must be a single positive number")
  alternative <- match.arg(alternative)
  if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
        (conf.level > 0) && (conf.level < 1)))
    stop("'conf.level' must be a single number between 0 and 1")
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  DF.x <- (1-x$fpp)*x$df.residual[2]
  if (DF.x < 1L)
    stop("not enough 'x' observations")
  DF.y <- (1-y$fpp)*y$df.residual[2]
  if (DF.y < 1L)
    stop("not enough 'y' observations")
  V.x <- (x$sigma)^2
  V.y <- (y$sigma)^2

  ESTIMATE <- V.x/V.y
  STATISTIC <- ESTIMATE/ratio
  PARAMETER <- c(`num edf` = DF.x, `denom edf` = DF.y)
  PVAL <- pf(STATISTIC, DF.x, DF.y)
  if (alternative == "two.sided") {
    PVAL <- 2 * min(PVAL, 1 - PVAL)
    BETA <- (1 - conf.level)/2
    CINT <- c(ESTIMATE/qf(1 - BETA, DF.x, DF.y), ESTIMATE/qf(BETA,
                                                             DF.x, DF.y))
  }
  else if (alternative == "greater") {
    PVAL <- 1 - PVAL
    CINT <- c(ESTIMATE/qf(conf.level, DF.x, DF.y), Inf)
  }
  else CINT <- c(0, ESTIMATE/qf(1 - conf.level, DF.x, DF.y))
  names(STATISTIC) <- "F"
  names(ESTIMATE) <- names(ratio) <- "ratio of variances"
  attr(CINT, "conf.level") <- conf.level
  RVAL <- list(statistic = STATISTIC, parameter = PARAMETER,
               p.value = PVAL, conf.int = CINT, estimate = ESTIMATE,
               null.value = ratio, alternative = alternative, method = "F test to compare two variances",
               data.name = DNAME)
  attr(RVAL, "class") <- "htest"
  return(RVAL)
}

#' Compare model fits of \code{intRvals} objects estimated on the same data.
#' If one object is provided, the results of a deviance test against a model without a missed event probability 'p'
#' is reported. If two objects are provided, the results of a deviance test between the model fits of the two objects is given.
#' @title Compares model fits of \code{intRvals} objects
#' @param object an object of class \code{intRvals}, usually a result of a call to \link[intRvals]{estinterval}
#' @param y an (optional) object of class \code{intRvals}, usually a result of a call to \link[intRvals]{estinterval}
#' @param conf.level confidence level for the deviance test
#' @param digits the number of digits for printing to screen
#' @param ... other arguments to be passed to low level functions
#' @export
#' @return A list of class "\code{anova.intRvals}" with the best model (1 or 2), deviance statistic and test results
#' \describe{
#'   \item{\code{best.model}}{the index of the best model (1 is first argument, 2 is second)}
#'   \item{\code{deviance}}{the deviance between the two tested models}
#'   \item{\code{p.value}}{p-value for the deviance (likelihood-ratio) test}
#'   \item{\code{conf.level}}{assumed confidence level for the test}
#'   \item{\code{model1.call}}{call that generated model 1}
#'   \item{\code{model2.call}}{call that generated model 2}
#'   \item{\code{AIC}}{numeric 2-vector containg the AIC value for model 1 (first element) and model 2 (second element)}
#'   \item{\code{loglik}}{numeric 2-vector containg the log-likelihood value for model 1 (first element) and model 2 (second element)}
#' }
#' @examples
#' data(goosedrop)
#' model1=estinterval(goosedrop$interval,fun="gamma")
#' # visually inspect model1 fit:
#' plot(model1)
#' # The observed distribution has intervals near zero.
#' # We allow a small random baseline to reduce the effect
#' # of intervals near zero on the fit result.
#' model2=estinterval(goosedrop$interval,fun="gamma",fpp.method='auto')
#' # model2 performs better than model1:
#' anova(model1,model2)
anova.intRvals=function(object, y=NULL, conf.level = 0.95,digits = max(3L, getOption("digits") - 3L), ...)
{
  x=object
  stopifnot(inherits(x, "intRvals"))
  if(inherits(y,"intRvals") && inherits(conf.level,"intRvals")) stop("can compare up to two models only")
  if(!missing(y)){
    stopifnot(inherits(y, "intRvals"))
    samedata=length(which(!(sort(x$data) == sort(y$data))))==0
    if(!samedata) stop("intRvals objects not estimated on the same dataset")
    if(!(length(which(x$trunc!=y$trunc))==0)) stop("intRvals objects do not have the same truncation range 'trunc'")
    deviance=2*(x$loglik[1]-y$loglik[1])
    # CHECK BELOW!
    p.value=pchisq(abs(deviance), df=max(c(1,abs(x$n.param-y$n.param))), lower.tail=FALSE)
    AIC=c(2*x$n.param-2*x$loglik[1],2*y$n.param-2*y$loglik[1])
    if(AIC[1]<AIC[2]) bestmodel=1L else bestmodel=2L
    if(AIC[1]==AIC[2]) bestmodel=NA
    output=list(best.model=bestmodel,deviance=abs(deviance),p.value=p.value,conf.level=conf.level,model1.call=x$call,model2.call=y$call,AIC=AIC,loglik=c(x$loglik[1],y$loglik[1]))
  }
  else{
    deviance=2*(x$loglik[1]-x$loglik[2])
    p.value=pchisq(abs(deviance), df=x$df.residual[1], lower.tail=FALSE)
    AIC=c(2*x$n.param-2*x$loglik[1],2*(x$n.param-1)-2*x$loglik[2])
    if(AIC[1]<AIC[2]) bestmodel=1L else bestmodel=2L
    if(AIC[1]==AIC[2]) bestmodel=NA
    output=list(best.model=bestmodel,deviance=abs(deviance),p.value=p.value,conf.level=conf.level,model1.call=x$call,model2.call="model 1 with missed event probability p fixed to zero.",
                AIC=AIC,loglik=x$loglik)
  }
  attr(output, "class") <- "anova.intRvals"
  return(output)
}

#' print method for class \code{anova.intRvals}
#'
#' @param x An object of class \code{anova.intRvals}, usually a result of a call to \link[intRvals]{anova.intRvals}
#' @keywords internal
#' @export

print.anova.intRvals=function(x,digits = max(3L, getOption("digits") - 3L), ...){
  stopifnot(inherits(x, "anova.intRvals"))
  cat("Model 1 call: ",deparse(x$model1.call),"\n")
  cat("Model 2 call: ",deparse(x$model2.call),"\n\n")

  if(x$p.value<1-x$conf.level){
    if(x$best.model==1) cat(paste("Model 1 outperforms Model 2, Pr(>Chi) =",format(signif(x$p.value,digits)),"\n\n"))
    if(x$best.model==2) cat(paste("Model 2 outperforms Model 1, Pr(>Chi) =",format(signif(x$p.value,digits)),"\n\n"))
  } else{
    cat(paste("models are equivalent, Pr(>Chi) =",format(signif(x$p.value,digits)),"\n\n"))
  }
  cat(paste("Model 1:   AIC =",format(signif(x$AIC[1],digits)), "   Loglik =",format(signif(x$loglik[1],digits)),"\n"))
  cat(paste("Model 2:   AIC =",format(signif(x$AIC[2],digits)), "   Loglik =",format(signif(x$loglik[2],digits))))
}

#' Estimate within-group variation
#'
#' Estimate within-group variation in interval length
#' @export
#' @param x object inheriting from class \code{intRvals}
#' @param conf.level confidence level passed to function \link[intRvals]{fundamental}, used in selecting fundamental intervals
#' @param alpha significance level for differences within and between groups or subjects
#' @param silent logical, if \code{TRUE} print no text to console
#' @return A logical atomic vector indicating which intervals are fundamental.
#' \describe{
#'   \item{\code{sigma.within}}{within-group standard deviation in interval length, estimated on fundamental intervals with repeated measures only}
#'   \item{\code{sigma}}{the total standard deviation in interval length, copied from \code{x$sigma}}
#'   \item{\code{p.within}}{p-value form a likelihood-ratio test indicating whether there is evidence for a random effect of group or subject}
#'   \item{\code{n.within}}{average number of intervals per group}
#'   \item{\code{n.total}}{total number of intervals}
#'   \item{\code{n.repeat}}{number of fundamental intervals with repeated measures, the size of the dataset on which \code{sigma.within} was estimated}
#'   \item{\code{p<alpha}}{logical. Whether there was significant evidence for a difference in within- and between-group/subject variance}
#' }
#' @details
#' Within- and between-group variation is estimated on the subset of fundamental intervals only.
#'
#' The subset of fundamental intervals is selected using \link[intRvals]{fundamental}.
#'
#' We calculate \eqn{sigma.within = s_w n_{ind}/(n_{ind}+1)} with \eqn{s_w} the uncorrected sample standard deviation
#' of within-group centered values (obtained from subtracting the group's mean value from each observation value),
#' and \eqn{n_{ind}/(n_{ind}+1)} Bessel's correction with \eqn{n_{ind}} the average number of repeated measures
#' per group. Significance of within-group variation is determined by testing for a random effect
#' of group against a constant null model (van de Pol & Wright 2009),
#' using the R-package lme4 (Bates et al. 2015).
#' @examples
#' # select the group of intervals observed on Terschelling island
#' dropset=goosedrop[goosedrop$site=="terschelling",]
#' # estimate an interval model, with separate within- and between-group variation:
#' dr=estinterval(data=dropset$interval,group = dropset$bout_id)
#' # plot the model fit:
#' plot(dr)
#' # estimate within-group variation, and its significance:
#' output=partition(dr)
#' # print within-group standard deviation:
#' output$sigma.within
#' # is the model including within-group standard deviation signicant,
#' # relative to a null model without separate within-group sd,
#' # at the specified confidence level alpha?
#' output$`p<alpha`  #> TRUE
#' @references
#' van de Pol, M. & Wright, J. (2009). A simple method for distinguishing within- versus between-subject effects using mixed models. Animal Behaviour, 77, 753-758.
#'
#' Bates, D., M\"{a}chler, M., Bolker, B.M. & Walker, S.C. (2015). Fitting linear mixed-effects models using lme4. Journal of Statistical Software, 67, 1-48.
partition=function(x,conf.level=0.9,alpha=0.05,silent=F){
  stopifnot(inherits(x,"intRvals"))
  if(length(x$group)==1 && is.na(x$group)) stop("no groups found in object")
  if(!(length(x$group)==length(x$data))) stop("'group' and 'data' are of unequal length")
  #extract intervals in the fundamental
  index=which(fundamental(x,conf.level))
  #check that we have fundamental intervals
  if(length(index)<2){
    warning(paste("insufficient fundamental intervals at confidence level alpha =",conf.level))
    return(list('p<alpha'=F,sigma.within=NA,sigma=x$sigma,p.within=NA,n.within=NA,n.repeat=0,n.total=length(x$data)))
  }
  dat.fund=data.frame(interval=x$data[index],group_id=x$group[index])
  # suppress CRAN checks on invisible binding
  interval=NULL
  # calculate group statistics
  dat.fund=ddply(dat.fund, c("group_id"), transform, offset = scale(interval,scale=F), mean = mean(interval), n = length(interval))
  if(!silent) cat(paste("intervals in fundamental:",nrow(dat.fund),"of",length(x$data)))
  dat.fund=dat.fund[dat.fund$n>1,]
  if(!silent) cat(paste(", with repeated measures:",nrow(dat.fund),"\n"))
  # test with mixed model term whether within-individual variation is significant
  lmer0=lm(interval ~ 1,data=dat.fund)
  lmer1=lmer(interval ~ (1 | group_id),data=dat.fund,REML=F)
  p.within=anova(lmer1,lmer0)$`Pr(>Chisq)`[2]
  # mean number of individuals per group
  n.ind=length(dat.fund$group_id)/length(unique(dat.fund$group_id))
  # Bessel's correction of within-subject standard deviation
  # this is required to account for the effect of within-subject centering
  sd.within=sd(dat.fund$offset)*sqrt(n.ind/(n.ind-1))
  accept=p.within<alpha
  if(!accept) warning("no evidence for a random effect of group or subject")
  list('p<alpha'=accept,sigma.within=sd.within,sigma=x$sigma,p.within=p.within,n.within=n.ind,n.repeat=nrow(dat.fund),n.total=length(x$data))
}
