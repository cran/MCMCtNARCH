\name{MCMCtNARCH}
\alias{mctnarch}
\alias{plotpara}
\alias{aus2005}
\title{MCMC estimation of Box-Cox transformed nonlinear GARCH (1, 1) model with Student-t Distribution }
\description{
\code{MCMCtNARCH} will perform a MCMC estimation of Box-Cox transformed nonlinear GARCH (1, 1) model with Student-t Distribution.Simulation inefficiency factor (SIF) is calculated and shown in result which is estimated as the variance of the sample mean divided by the variance of the sample mean from a hypothetical sampler that draws independent random observations from the posterior distribution.
}

\usage{
mctnarch(xdata, m) 
plotpara(para.matrix, m)
}

\arguments{
  \item{xdata}{financial return data }
  \item{m}{MCMC sampling iteration times }
  \item{para.matrix}{Sampled matrix of parameters  for M iterations }
}

\value{
  
  \item{para.matrix }{ a matrix of returned values for the estimated parameters }
  \item{aus2005}{Percentage Return of AOI stock price data from 2005 to 2008 i.e.[LN(Pt)-LN(Pt-1)]*100}
}

\author{ Ken Ren \email{ken.ren@abs.gov.au}}

\examples{
##
##-- Test data: The following Return of AOI stock price data 2005-2008 is used to example purpose 
##-- \\MCMCTNARCH\DATA\aus2005.txt
##xt<-matrix(scan(file="\\MCMCTNARCH\DATA\ aus2005.txt"),ncol=1,byrow=T)

data(aus2005)
xt<-aus2005
m<-1000
mc=mctnarch(xt, m)
plotpara(mc, m)

}
\keyword{ distribution }
\keyword{ models }


