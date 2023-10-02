#' @title The probability density function (pdf) of a composite distribution with Pareto tail
#' @description \code{dcomppareto} returns the density of a composite distribution with a Pareto upper tail at a point x, with a specified distribution at the lower tail.
#' @param x A scalar or vector of positive values at which the density needs to be evaluated
#' @param spec The selection of the lower tail (head) distribution
#' @param alpha The shape parameter of the Pareto distribution
#' @param theta The scale parameter of Pareto, also serve as the location parameter of the composite model
#' @param log logical; if TRUE, probability p are given as log(p)
#' @param ... The parameter of the lower tail (head) distribution
#' @return an object of the same length of x as the density evaluated at x
#' @import actuar
#' @export
#' @examples
#' x<-1:100
#' dcomppareto(x, "lnorm", 0.4, 1, meanlog = 1, sdlog = 0.8)
#' dcomppareto(x, "weibull", alpha = 1.5, theta = 1.5, shape = 2, scale = 2)

dcomppareto<-function(x, spec, alpha = 1,theta=1, log=FALSE, ...)
{
  g1<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
  G1<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

  g2<-function (z) {actuar::dpareto(z,shape = alpha,scale = theta)}
  G2<-function (z) {actuar::ppareto(z,shape = alpha,scale = theta)}
  phi<-( g1(theta,...)*(1-G2(theta) ) / ( g2(theta)*G1(theta,...) ) )#continuous condition

  pdf<-x
  pdf[x<=theta&log==FALSE]<-g1(x[x<=theta],...) / ((1+phi)*G1(theta,...))
  pdf[x<=theta&log==TRUE]<-log(g1(x[x<=theta],...) / ((1+phi)*G1(theta,...)))
  pdf[x>theta&log==FALSE]<-phi*g2(x[x>theta])/((1+phi)*(1-G2(theta)))
  pdf[x>theta&log==TRUE]<-log(phi*g2(x[x>theta])/((1+phi)*(1-G2(theta))))

  return(pdf)
}

#' @title The cumulative distribution function (CDF) of a composite distribution with Pareto tail
#' @description \code{pcomppareto} returns the CDF of a composite distribution with a Pareto upper tail at x, with a specified distribution at the lower tail.
#' @param x A scalar or vector of positive values at which the CDF needs to be evaluated
#' @param spec The selection of the lower tail (head) distribution
#' @param alpha The shape parameter of the Pareto distribution
#' @param theta The scale parameter of Pareto, also serve as the location parameter of the composite model
#' @param log.p logical; if TRUE, probability p are given as log(p)
#' @param lower.tail logical; if FALSE, the upper tail probability is provided
#' @param ... The parameter of the lower tail (head) distribution
#' @return an object of the same length of x as the CDF evaluated at x
#' @import actuar
#' @export
#' @examples
#' x<-1:100
#' pcomppareto(x, "lnorm", 0.4, 1, meanlog = 1, sdlog = 0.8)
#' pcomppareto(x, "weibull", alpha = 1.5, theta = 1.5, shape = 2, scale = 2)

pcomppareto<-function(x, spec, alpha=1, theta=1, lower.tail=TRUE, log.p = FALSE , ...)
{
  g1<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
  G1<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

  g2<-function (z) {dpareto(z,shape = alpha,scale = theta)}
  G2<-function (z) {ppareto(z,shape = alpha,scale = theta)}
  phi<-( g1(theta,...)*(1-G2(theta) ) / ( g2(theta)*G1(theta,...) ) )#continuous condition

  cdf<-x
  cdf[x<=theta&lower.tail==TRUE&log.p==FALSE]<-G1(x[x<=theta],...) / ((1+phi)*G1(theta,...))
  cdf[x<=theta&lower.tail==FALSE&log.p ==FALSE]<-1-G1(x[x<=theta],...) / ((1+phi)*G1(theta,...))
  cdf[x>theta&lower.tail==TRUE&log.p==FALSE]<-1/(1+phi) + (phi/(1+phi))*(G2(x[x>theta])-G2(theta))/(1-G2(theta))
  cdf[x>theta&lower.tail==FALSE&log.p==FALSE]<-phi/(1+phi) - (phi/(1+phi))*(G2(x[x>theta])-G2(theta))/(1-G2(theta))

  cdf[x<=theta&lower.tail==TRUE&log.p==TRUE]<-log(G1(x[x<=theta],...) / ((1+phi)*G1(theta,...)))
  cdf[x<=theta&lower.tail==FALSE&log.p ==TRUE]<-log(1-G1(x[x<=theta],...) / ((1+phi)*G1(theta,...)))
  cdf[x>theta&lower.tail==TRUE&log.p==TRUE]<-log(1/(1+phi) + (phi/(1+phi))*(G2(x[x>theta])-G2(theta))/(1-G2(theta)))
  cdf[x>theta&lower.tail==FALSE&log.p==TRUE]<-log(phi/(1+phi) - (phi/(1+phi))*(G2(x[x>theta])-G2(theta))/(1-G2(theta)))
  return(cdf)
}

#' @title The quantile function of a composite distribution with Pareto tail
#' @description \code{qcomppareto} returns the quantile of a composite distribution with a Pareto upper tail given p, with a specified distribution at the lower tail.
#' @param p vector of probabilities
#' @param spec The selection of the lower tail (head) distribution
#' @param alpha The shape parameter of the Pareto distribution
#' @param theta The scale parameter of Pareto, also serve as the location parameter of the composite model
#' @param log.p logical; if TRUE, probability p are given as log(p)
#' @param ... The parameter of the lower tail (head) distribution
#' @return an object of the same length of x as the CDF evaluated at x
#' @import actuar
#' @export
#' @examples
#' p <-seq(0.01,0.99,b=0.01)
#' qcomppareto(p, "weibull", alpha = 1.5, theta = 1.5, shape = 2, scale = 2)

qcomppareto<-function(p, spec, alpha=1, theta=1, log.p=FALSE , ...)
{
  if (log.p==TRUE) {
    p = exp(p)
  }

  g1<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
  G1<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
  G1q<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

  g2<-function (z) {dpareto(z,shape = alpha,scale = theta)}
  G2<-function (z) {ppareto(z,shape = alpha,scale = theta)}
  G2q<-function (z) {qpareto(z,shape = alpha,scale = theta)}

  phi<-( g1(theta,...)*(1-G2(theta) ) / ( g2(theta)*G1(theta,...) ) )#continuous condition

  qf<-p
  qf[p<=1/(1+phi)]<-G1q( p[p<=1/(1+phi)]*(1+phi) * G1(theta,...),... )
  qf[p>1/(1+phi)]<-G2q( (p[p>(1/(1+phi))]*(1+phi)-1)/phi * (1-G2(theta))+ G2(theta))
  return(qf)
}

#' @title Generating random number from a discrete composite distribution with Pareto tail
#' @description \code{rcomppareto} returns a random sample of a composite distribution with a Pareto upper tail, with a specified distribution at the lower tail.
#' @param n number of observations
#' @param spec The selection of the lower tail (head) distribution
#' @param alpha The shape parameter of the Pareto distribution
#' @param theta The scale parameter of Pareto, also serve as the location parameter of the composite model
#' @param ... The parameter of the lower tail (head) distribution
#' @return an object of the same length of n
#' @import actuar
#' @importFrom stats runif
#' @export
#' @examples
#' n<-100
#' rcomppareto(n,"weibull", alpha = 1.5, theta = 1.5, shape = 2, scale = 2)
rcomppareto = function(n,spec,alpha = 1,theta=1, ...){
  rand = runif(n,0,1)
  rand_comppareto = qcomppareto(rand,spec,alpha = alpha,theta = theta, ...)
  return(rand_comppareto)
}


#' @title The probability mass function (pmf) of a discrete composite distribution with Pareto tail
#' @description \code{dwdcomppareto} returns the pmf of a discrete composite distribution with a Pareto upper tail at a point x, with a specified distribution at the lower tail.
#' @param x A scalar or vector of nonnegative integer values at which the probability mass needs to be evaluated
#' @param spec The selection of the lower tail (head) distribution
#' @param alpha The shape parameter of the Pareto distribution
#' @param theta The scale parameter of Pareto, also serve as the location parameter of the composite model
#' @param log logical; if TRUE, probability p are given as log(p)
#' @param ... The parameter of the lower tail (head) distribution
#' @return an object of the same length of x as the probability mass evaluated at x
#' @import actuar
#' @export
#' @examples
#' x<-1:100
#' dwdcomppareto(x, "lnorm", 0.4, 1, meanlog = 1, sdlog = 0.8)
#' dwdcomppareto(x, "weibull", alpha = 1.5, theta = 1.5, shape = 2, scale = 2)
dwdcomppareto <- function(x,spec,alpha,theta,log = FALSE,...){
  test <- c()
  for(i in 1:length(x)){
    test[i] <- x[i]%%1==0 & x[i]>=0
  }
  crit <- sum(test)
  if(crit==length(x)) {
    pmf<-pcomppareto(x+1,spec, alpha,theta,log.p = FALSE,...)-pcomppareto(x,spec, alpha,theta,log.p = FALSE,...)}
  else{
    warning("dwdcomppareto is only defined for nonnegative integers.")
  }
  if(log ==TRUE){
    pmf <- log(pmf)
  }
  return(pmf)
}


#' @title The cumulative distribution function (CDF) of a discrete composite distribution with Pareto tail
#' @description \code{pwdcomppareto} returns the CDF of a discrete composite distribution with a Pareto upper tail at x, with a specified distribution at the lower tail.
#' @param x A scalar or vector of positive values at which the CDF needs to be evaluated
#' @param spec The selection of the lower tail (head) distribution
#' @param alpha The shape parameter of the Pareto distribution
#' @param theta The scale parameter of Pareto, also serve as the location parameter of the composite model
#' @param log.p logical; if TRUE, probability p are given as log(p)
#' @param ... The parameter of the lower tail (head) distribution
#' @return an object of the same length of x as the CDF evaluated at x
#' @import actuar
#' @export
#' @examples
#' x<-1:100
#' pwdcomppareto(x, "lnorm", 0.4, 1, meanlog = 1, sdlog = 0.8)
#' pwdcomppareto(x, "weibull", alpha = 1.5, theta = 1.5, shape = 2, scale = 2)
pwdcomppareto <- function(x,spec,alpha,theta,log.p = FALSE,...){
  cdf<-pcomppareto(floor(x)+1,spec, alpha,theta,log.p = FALSE,...)
  if(log.p==TRUE){
    cdf <-log(cdf)
  }
  return(cdf)
}

#' @title The quantile function of a discrete composite distribution with Pareto tail
#' @description \code{qwdcomppareto} returns the quantile of a composite distribution with a Pareto upper tail given p, with a specified distribution at the lower tail.
#' @param p vector of probabilities
#' @param spec The selection of the lower tail (head) distribution
#' @param alpha The shape parameter of the Pareto distribution
#' @param theta The scale parameter of Pareto, also serve as the location parameter of the composite model
#' @param log logical; if TRUE, probability p are given as log(p)
#' @param ... The parameter of the lower tail (head) distribution
#' @return an object of the same length of x as the CDF evaluated at x
#' @import actuar
#' @export
#' @examples
#' p <-seq(0.1,0.9,b=0.1)
#' qcomppareto(p, "weibull", alpha = 1.5, theta = 1.5, shape = 2, scale = 2)
qwdcomppareto<-function(p,spec,alpha,theta,log = FALSE,...){
  if (log==TRUE) {
    p = exp(p)
  }
  x = c()
  for(i in 1:length(p)){
    if(p[i]>=0 & p[i]<1){
      x[i] <-0
      p0 =  pwdcomppareto(x[i],spec,alpha,theta,log.p = FALSE,...)
      while(p0<= p[i]){
        x[i] <- x[i]+1
        p0 <- pwdcomppareto(x[i],spec,alpha,theta,log.p = FALSE,...)
      }
    }
    else{
      warning("True probability should be between 0 and 1.")
    }
  }
  return(x)
}


#' @title Generating random number from a discrete composite distribution with Pareto tail
#' @description \code{rwdcomppareto} returns a random sample of a discrete composite distribution with a Pareto upper tail, with a specified distribution at the lower tail.
#' @param n number of observations
#' @param spec The selection of the lower tail (head) distribution
#' @param alpha The shape parameter of the Pareto distribution
#' @param theta The scale parameter of Pareto, also serve as the location parameter of the composite model
#' @param ... The parameter of the lower tail (head) distribution
#' @return an object of the same length of n
#' @import actuar
#' @importFrom stats runif
#' @export
#' @examples
#' n<-10
#' rcomppareto(n,"weibull", alpha = 1.5, theta = 1.5, shape = 2, scale = 2)
rwdcomppareto = function(n,spec,alpha = 1,theta=1, ...){
  p = runif(n)
  sample = qwdcomppareto(p,spec,alpha,theta,log.p = FALSE,...)
  return(sample)
}






