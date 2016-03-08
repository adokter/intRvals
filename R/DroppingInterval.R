## open problem convergence:
# this works:
# dr.schiermonnikoog=droprate(goosedrop[goosedrop$site=="schiermonnikoog",]$interval,fun="normal",mu=247,sigma=83)
# this not:
# dr.schiermonnikoog=droprate(goosedrop[goosedrop$site=="schiermonnikoog",]$interval,fun="normal")



terschelling=read.csv("~/Dropbox/metawad/veldwerk terschelling/dropping_intervals_terschelling.csv",colClasses=c("POSIXct","numeric"))
schiermonnikoog=read.csv("~/Dropbox/metawad/veldwerk terschelling/dropping_intervals_schier.csv",colClasses=c("POSIXct","numeric"))
terschelling$site="terschelling"
schiermonnikoog$site="schiermonnikoog"
goosedrop=rbind(terschelling,schiermonnikoog)
#save(goosedrop,file="~/git/R/droprate/data/goosedrop.RData")
load("~/git/R/droprate/data/goosedrop.RData")
#' This is data to be included in my package
#'
#' @author Adriaan Dokter \email{a.m.dokter@uva.nl}
#' @references \url{data_blah.com}
"goosedrop"

# normal distribution for (i-1) missed events component of the
# probability density function (PDF)
normi=  function(x,mu,sigma,p,i) (p^(i-1)-p^i)*dnorm(x,i*mu,sqrt(i)*sigma)
# gamma distribution for (i-1) missed events component of the PDF
gammai=  function(x,mu,sigma,p,i) (p^(i-1)-p^i)*dgamma(x,shape=i*mu^2/(sigma^2),scale=(sigma^2)/mu)

# normalized sum of all components of the PDF
normsum=function(x,mu,sigma,p,N,fun=normi) sapply(1:N,FUN=fun,x=x,p=p,mu=mu,sigma=sigma)
# the PDF
probdens=function(x,mu,sigma,p,N,fun=normi) {
  if(length(x)>1) return(rowSums(normsum(x,mu,sigma,p,N,fun=fun)))
  else return(sum(normsum(x,mu,sigma,p,N,fun=fun)))
}
# log-likelihood
loglik=function(data,mu,sigma,p,N,fun=normi) sum(log(probdens(data,mu,sigma,p,N,fun=fun)))
# log-likelihood, with parameters to be optimised in list params
loglikfn=function(params,data,N,fun=normi) sum(log(probdens(data,params[1],params[2],plogis(params[3]),N,fun=fun)))
# log-likelihood of a null model without miss chance (i.e. N=1, p=0)
logliknull=function(params,data,N,fun=normi) sum(log(probdens(data,params[1],params[2],0,N,fun=fun)))


#' Plot an interval histogram and fit of droprate object
#'
#' @param object An droprate class object
#' @param binsize Width of the histogram bins
#' @param line.col Color of the plotted curve for the model fit
#' @param main an overall title for the plot
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param ... Additional arguments to be passed to the low level plotting functions
#' @export
#' @return This function returns a list with data, corresponding to the model fit
#' @examples
#' data(goosedrop)
#' dr=droprate(goosedrop$interval)
#' plot(dr)
#' plot(dr,binsize=10,line.col='blue')
plot.droprate=function(object,binsize=20,xlab="Interval",ylab="Density",main="Interval histogram and fit", line.col='red', ...){
  stopifnot(inherits(object, "droprate"))
  if(object$distribution=="normal") fun2use=normi
  else fun2use=gammai
  hist(object$data,freq=F,breaks=seq(0,max(object$data)+binsize,binsize),xlab=xlab,ylab=ylab,main=main,...)
  curve(probdens(x,object$mean,object$stdev,object$fractionMissed,object$N,fun=fun2use),0,1500,col=line.col,add=T)
}

#' Analyse event rate with missed event observations
#'
#' @param data A numeric list of intervals.
#' @param mu Start value for the numeric optimization for the mean interval rate.
#' @param sigma Start value for the numeric optimization for the standard deviation of the interval rate.
#' @param p Start value for the numeric optimization for the chance to miss an observation of an event.
#' @param N Maximum number of missed observations to be taken into account.
#' @param fun Assumed distribution for the intervals, one of "\code{normal}" or "\code{gamma}", corresponding
#' to the \link[stats]{Normal} and \link[stats]{GammaDist} distributions
#' @export
#' @return This function returns an object of class \code{droprate}
#' @examples
#' data(goosedrop)
#' dr=droprate(goosedrop$interval)
#' plot(dr)
#' summary(dr)
#'
#' # let's analyze dropping intervals by site
#' dr1=droprate(goosedrop[goosedrop$site=="schiermonnikoog",]$interval)
#' dr2=droprate(goosedrop[goosedrop$site=="terschelling",]$interval)
#'
#' # The T test that accounts for missed observations shows the
#' # mean dropping interval is not significantly different at the two sites:
#' ttest(dr1,dr2)
#' # not accounting for missed observations leads to the false
#' # conclusion that the mean intervals are significantly different:
#' t.test(dr1$data,dr2$data)
#'
droprate=function(data,mu=mean(data),sigma=sd(data)/2,p=0.2,N=5,fun="normal"){
  call=match.call()
  if(min(data)<0) stop("data contains one or more negative intervals, only positive intervals allowed.")
  if (!(fun=="normal" || fun=="gamma")) stop("fun needs to be either 'normal' or 'gamma'")
  if(fun=="normal") fun2use=normi
  else fun2use=gammai
  p.logit=log(p/(1-p))
  opt=optim(par=c(mu,sigma,p.logit),loglikfn,data=data,N=N,fun=fun2use,control=list(fnscale=-1))
  optnull=optim(par=c(mu,sigma),logliknull,data=data,N=1,fun=fun2use,control=list(fnscale=-1))
  out=list(call=call,data=data,mean=opt$par[1],stdev=opt$par[2],fractionMissed=plogis(opt$par[3]),N=N,convergence=opt$convergence,counts=opt$counts,loglik=c(opt$value,optnull$value),df.residual=c(1,length(data)-3),distribution=fun)
  goodness.fit=pchisq(2*(out$loglik[1]-out$loglik[2]), df=out$df.residual[1], lower.tail=FALSE)
  if(goodness.fit>0.05){
    warning("model including a miss chance is not significant.\n\nCheck convergence using different starting values (arguments: mu,sigma and p), or change the interval distribution (argument fun)")
  }
  class(out)="droprate"
  out
}

#' summary method for class \code{droprate}
#'
#' @param x An object of class \code{droprate}, usually a result of a call to \link[droprate]{droprate}
#' @export
#' @return The function \code{summary.droprate} computes and returns a list of summary statistics
#' @examples
#' data(goosedrop)
#' dr=droprate(goosedrop$interval)
#' summary(dr)
summary.droprate=function(x){
  stopifnot(inherits(x, "droprate"))
  xx=x
  xx$deviance=c(2*(xx$loglik[1]-xx$loglik[2]),abs(2*xx$loglik[1]))
  xx$p.value=c(pchisq(xx$deviance[1], df=xx$df.residual[1], lower.tail=FALSE),pchisq(xx$deviance[2], df=xx$df.residual[2], lower.tail=FALSE))
  class(xx)="summary.droprate"
  xx
}

#' print method for class \code{droprate}
#'
#' @param x An object of class \code{droprate}, usually a result of a call to \link[droprate]{droprate}
#' @keywords internal
#' @export
print.droprate=function(x,digits = max(3L, getOption("digits") - 3L)){
  stopifnot(inherits(x, "droprate"))
  cat("Analysis of event interval data with missed event observations\n\n")
  cat("  number of intervals: ",format(signif(length(x$data),digits)),"\n\n")
  cat("      mean event rate: ",format(signif(x$mean,digits)),"\n")
  cat("   standard deviation: ",format(signif(x$stdev,digits)),"\n")
  cat("      fraction missed: ",format(signif(x$fractionMissed)),"\n")
}

#' print method for class \code{summary.droprate}
#'
#' @param x An object of class \code{symmary.droprate}, usually a result of a call to \link[droprate]{summary.droprate}
#' @keywords internal
#' @export
print.summary.droprate=function(x,digits = max(3L, getOption("digits") - 3L)){
  stopifnot(inherits(x, "summary.droprate"))
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("             distribution used: ",x$distribution,"\n")
  cat("           number of intervals: ",format(signif(length(x$data),digits)),"\n")
  cat("                Log-likelihood: ",format(signif(x$loglik[1],digits)),"\n")
  cat("             Residual deviance: ",format(signif(x$deviance[1]))," on ", x$df.residual[1], "degrees of freedom\n")
  cat("Likelihood ratio test, p-value: ",format(signif(x$p.value[1],digits)), "(against null model without miss chance)\n")
  cat("             Residual deviance: ",format(signif(x$deviance[2]))," on ", x$df.residual[2], "degrees of freedom\n")
  cat("Likelihood ratio test, p-value: ",format(signif(x$p.value[2],digits)), "(against saturated null model)\n\n")
  cat("               mean event rate: ",format(signif(x$mean,digits)),"\n")
  cat("            standard deviation: ",format(signif(x$stdev,digits)),"\n")
  cat("               fraction missed: ",format(signif(x$fractionMissed)),"\n")
}

#' Performs one and two sample t-tests on objects of class \code{droprate}
#' @title Student's t-test for objects of class \code{droprate}
#' @param x an object of class \code{droprate}, usually a result of a call to \link[droprate]{droprate}
#' @param y an (optional) object of class \code{droprate}, usually a result of a call to \link[droprate]{droprate}
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
#' dr=droprate(goosedrop$interval)
#' # perform a one-sample t-test
#' ttest(dr)
#' # two sample t-test
#' data.beforeMay=goosedrop[goosedrop$date<as.POSIXct('2013-05-01'),]
#' data.afterMay=goosedrop[goosedrop$date>as.POSIXct('2013-05-01'),]
#' dr.beforeMay=droprate(data.beforeMay$interval)
#' dr.afterMay=droprate(data.afterMay$interval)
#' ttest(dr.beforeMay,dr.afterMay)
#'
ttest = function (x, y = NULL, alternative = c("two.sided", "less", "greater"),
          mu = 0, var.equal = FALSE, conf.level = 0.95)
{
  stopifnot(inherits(x, "droprate"))
  if (!is.null(y)){
    stopifnot(inherits(y, "droprate"))
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
  nx <- length(x$data)
  mx <- x$mean
  vx <- (x$stdev)^2
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
    ny <- length(y)
    if (nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if (ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if (var.equal && nx + ny < 3)
      stop("not enough observations")
    my <- y$mean
    vy <- (y$stdev)^2
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

#' Performs an F test to compare the variances of objects of class \code{droprate}
#' @title F Test to compare two variances of objects of class \code{droprate}
#' @param x an object of class \code{droprate}, usually a result of a call to \link[droprate]{droprate}
#' @param y an (optional) object of class \code{droprate}, usually a result of a call to \link[droprate]{droprate}
#' @param alternative a character string specifying the alternative hypothesis, must be one of "\code{two.sided}" (default), "\code{greater}" or "\code{less}". You can specify just the initial letter.
#' @param ratio the hypothesized ratio of the population variances of \code{x} and \code{y}.
#' @param conf.level confidence level for the returned confidence interval
#' @details The null hypothesis is that the ratio of the variances of the
#' data to which the models \code{x} and \code{y} were fitted, is equal to ratio.
#' @export
#' @return A list with class "\code{htest}" containing the same components as in \link[stats]{var.test}
#' @examples
#' data(goosedrop)
#' dr=droprate(goosedrop$interval)
#' # split the interval data into two periods
#' data.beforeMay=goosedrop[goosedrop$date<as.POSIXct('2013-05-01'),]
#' data.afterMay=goosedrop[goosedrop$date>as.POSIXct('2013-05-01'),]
#' dr.beforeMay=droprate(data.beforeMay$interval)
#' dr.afterMay=droprate(data.afterMay$interval)
#' # perform an F test
#' vartest(dr.beforeMay,dr.afterMay)
#'
vartest=function (x, y, ratio = 1, alternative = c("two.sided", "less",
                                           "greater"), conf.level = 0.95)
{
  stopifnot(inherits(x, "droprate") && inherits(y, "droprate"))
  if (!((length(ratio) == 1L) && is.finite(ratio) && (ratio >
                                                      0)))
    stop("'ratio' must be a single positive number")
  alternative <- match.arg(alternative)
  if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
        (conf.level > 0) && (conf.level < 1)))
    stop("'conf.level' must be a single number between 0 and 1")
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  DF.x <- x$df.residual[2]
  if (DF.x < 1L)
    stop("not enough 'x' observations")
  DF.y <- y$df.residual[2]
  if (DF.y < 1L)
    stop("not enough 'y' observations")
  V.x <- (x$stdev)^2
  V.y <- (y$stdev)^2

  ESTIMATE <- V.x/V.y
  STATISTIC <- ESTIMATE/ratio
  PARAMETER <- c(`num df` = DF.x, `denom df` = DF.y)
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
