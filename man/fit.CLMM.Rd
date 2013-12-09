\name{fit.CLMM}

\alias{fit.CLMM}

\title{
  Clustering of Linear Mixed Models Method
}

\description{
  Fit a CLMM model for single or replicated time course data.
}

\usage{
  fit.CLMM(data.y, data.x, data.z, n.clst, n.start = 1)
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
  \item{n.clst}{
    an integer, number of clusters.
  }
  \item{n.start}{
    an integer used to get the starting values for the EM algorithm.
  }
}

\value{
  \item{u.hat}{
    a matrix containing the cluster membership information for the genes.
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
  of Qin and Selft (2006). This method clusters genes based on
  the estimated parameters of the regression functions that are used to
  model the gene expression.
}

\author{ Li-Xuan Qin 
  \email{qinl@mskcc.org}
}

\references{
  Li-Xuan Qin and Steven G. Self (2006).  The clustering of regression models
  method with application in gene expression data.  \emph{Biometrics}
  62, 526-533.
  \url{ http://www.mskcc.org/mskcc/html/60448.cfm}  
}

\seealso{
  \code{\link{fit.CLM}},
  \code{\link{fit.CLMM}}
}

\examples{
#test data
  data(YeastCellCycle)
  data.y <- YeastCellCycle$normalizedData.sample
  data.x <- YeastCellCycle$designMatrix
#fit the model
  n.clst <- 6
  fit1 <- fit.CLMM(data.y, data.x, data.x, n.clst)

  fit1.u <- apply(fit1$u.hat, MARGIN=1, FUN=order, decreasing=TRUE)[1,]
  zeta.fitted <- fit1$theta.hat$zeta.hat
  profile <- data.x[1,,] \%*\% t(zeta.fitted)
#display the profile of each cluster
  n.knots <- 7
  plot.x <- n.knots*(1:dim(data.y)[3]-1)

  par(mfrow=c(2, ceiling((n.clst)/2)),mai=c(0.5,0.5,0.5,0.1),mgp=c(1,0.3,0))
  for(k in 1:n.clst){
  # plot the fitted cluster-specific profiles
    plot(plot.x,profile[,k],type="l", 
         ylim=c(-2,2), main=paste("Cluster",k), 
         xlab="time (min)", ylab=NA,xaxt="n",lwd=2)
    axis(side=1, at=plot.x[(1:8)*2-1], labels=paste(plot.x[(1:8)*2-1]), cex.axis=0.8)
  # plot the observed profiles for genes in this cluster
    i.plot <- (1:dim(data.y)[1])[fit1.u==k]
    for(j in i.plot) { lines(plot.x, data.y[j,1,], lty=3, lwd=1)}
    text(84,-1.9, paste(length(which(fit1.u==k)),"genes"))
  }

}

\keyword{Functions}

