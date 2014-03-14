\name{fit.CLMM.NA}

\alias{fit.CLMM.NA}

\title{
  Clustering of Linear Mixed Models Method
}

\description{
  Fit a CLMM model for time course data allowing for missing value.
}

\usage{
  fit.CLMM.NA(data.y, data.x, data.z, na.flag, n.clst, n.start = 1)
}

\arguments{
  \item{data.y}{
    a three dimensional array of gene expression data, data.y[j, i, t] for gene j, sample i,  and time point t.
  }
  \item{data.x}{
    a three dimensional array of fixed effects (common for all genes), data.x[i, t, p] for sample i, time point t, and covariate p.
  }
  \item{data.z}{
    a three dimensional array of random effects (common for all genes), data.z[i, t, q] for sample i, time point t, and covariate p.
  }
  \item{na.flag}{
    a three dimensional array of indicator for missing data, 1 for missing, 0 for otherwise.
  }
  \item{n.clst}{
    an integer, number of clusters.
  }
  \item{n.start}{
    an integer used to get the starting value for the EM algorithm.
  }
}

\value{
  \item{u.hat}{
    a matrix containing the cluster membership for the genes.
  }
  \item{b.hat}{
    an array containing the estimated random effects.
  }
  \item{theta.hat}{
    a list comprised of five components: \emph{zeta.hat}, \emph{pi.hat}, \emph{D.hat}, \emph{sigma2.hat} and \emph{llh}. They are described as below:
  }
  \item{zeta.hat}{
    a matrix of the estimated fixed effects with one row for each cluster.
  }
  \item{pi.hat}{
    a vector of the relative frequency for each cluster.
  }
  \item{D.hat}{
    the estimated random effects variances for the clusters.
  }
  \item{sigma2.hat}{
    the estimated measurement error variances for the clusters.
  }
  \item{llh}{
    log likelihood for the model.
  }
}

\details{
  This function implements the Clustering of Linear Mixed Models Method
  of Qin and Self (2006).
}

\author{ Li-Xuan Qin 
  \email{qinl@mskcc.org}
}

\references{
 \itemize{
  \item Li-Xuan Qin and Steven G. Self (2006). The clustering of regression models
  method with applications in gene expression data.  \emph{Biometrics}
  62, 526-533.
  \item Li-Xuan Qin, Linda Breeden and Steven G. Self (2014). Finding gene clusters 
  for a replicated time course study. \emph{BMC Res Notes} 7:60.  
 }
}

\seealso{
  \code{\link{fit.CLM}},
  \code{\link{fit.CLMM}},
  \code{\link{fit.CLMM.2}},
  \code{\link{fit.CLMM.NA}},
  \code{\link{fit.CLMM.NA.2}}
}

\examples{
\dontrun{
#test data
  data(YeastCellCycle2)
  data.y <- YeastCellCycle2$normalizedData.WT
  data.x <- YeastCellCycle2$designMatrix.WT
  na.flag<- YeastCellCycle2$NA.flag.WT
#fit the model
  n.clst <- 8
  fit1   <- fit.CLMM.NA(data.y,data.x[,,1:9],data.x[,,1:9],na.flag,n.clst)
  fit1.u <- apply(fit1$u.hat, MARGIN=1, FUN=order, decreasing=TRUE)[1,]
  zeta.fitted <- fit1$theta.hat$zeta.hat
  profile.WT <- YeastCellCycle2$designMatrix.WT[1,,1:9] \%*\% t(zeta.fitted)
#display the profile of each cluster
 # remove bad time points for WTs
  n.time  <- 25
  time.WT <- (1:n.time)[-22]
  n.rpl.WT<- length(time.WT)
  n.gene.short<- dim(data.y)[1]
 # gene-specific profile: observed profiles averaged over replicates
  data.WT.mean  <- matrix(0, nrow=n.gene.short, ncol=n.rpl.WT)
  for(j in 1:n.gene.short){
   data.WT.mean[j,]	<- (data.y[j,1,]+data.y[j,2,])/2
  }
 # plot observed profiles by cluster
  col.panel=rainbow(8)
  par(mfrow=c(3, 3),mai=c(0.3,0.25,0.2,0.05))
  for(k in 1:n.clst){
  plot(5*(time.WT-1), profile.WT[,k], type="l", col=col.panel[k], ylim=c(-2,2),
       xlab="Time", ylab="Expression Value", main=paste("WT: cluster",k))
  i.plot 	<- (1:n.gene.short)[fit1.u==k]
  for(j in i.plot) lines(5*(time.WT-1), data.WT.mean[j,], lty=1)
  lines(5*(time.WT-1), profile.WT[,k], col=col.panel[k], lwd=2)
  text(125, -1.9, pos=2, paste(length(i.plot)," genes"))
  }
}
}

\keyword{Functions}

