\name{fit.CLM}

\alias{fit.CLM}

\title{
  Clustering of Linear Models Method
}

\description{
  Fit a CLM model for cross-sectional data.
}

\usage{
  fit.CLM(data.y, data.x, n.clst, n.start = 1)
}

\arguments{
  \item{data.y}{
    matrix of gene expression data, data.y[j, i] for sample i and gene j.
  }
  \item{data.x}{
    matrix of sample covariates, data.x[i, p] for sample i and covariate p.
  }
  \item{n.clst}{
    an integer, number of clusters .
  }
  \item{n.start}{
    an integer used to get the starting values for the EM algorithm.
  }
}

\value{
  \item{u.hat}{
    a matrix containing the cluster membership probability for each gene, whose row names are genes and column names are clusters.
  }
  \item{theta.hat}{
    a list comprised of four components: \emph{zeta.hat}, \emph{pi.hat}, \emph{sigma2.hat}, \emph{llh}. They are described as below:
  }
  \item{zeta.hat}{
      a matrix with the estimated regression parameters with one row for each cluster.
    }
  \item{pi.hat}{
      a vector with the relative frequency for each cluster.
    }
  \item{sigma2.hat}{
      a vector of variance parameters.
    }
  \item{llh}{
      log likelihood for the model.
    }
}

\details{
This function implements the Clustering of Linear Models Method of Qin (2006). This method clusters genes based on the estimated regression parameters that model the relation between gene expression and sample covariates. 
}

\author{ Li-Xuan Qin 
  \email{qinl@mskcc.org}
}

\references{
  Li-Xuan Qin and Steven G. Self (2006).  The clustering of regression models
  method with applications in gene expression data.  \emph{Biometrics}
  62, 526-533.
  \url{ http://www.mskcc.org/mskcc/html/60448.cfm}  
}

\seealso{
  \code{\link{fit.CLM}},
  \code{\link{fit.CLMM}}
}

\examples{
#test data
  data(BreastCancer)
  data.y <- BreastCancer$normalizedData
  data.x <- BreastCancer$designMatrix
#fit the model
  n.clst <- 9
  fit1   <- fit.CLM(data.y, data.x, n.clst)
  fit1.u <- apply(fit1$u.hat, MARGIN=1, FUN=order, decreasing=TRUE)[1,]
#display the results
  index.IDC <- which(data.x[,2]==0)
  index.ILC <- which(data.x[,2]==1)
  mean.IDC  <- apply(data.y[,index.IDC], MARGIN=1, FUN=mean, na.rm=TRUE)
  mean.ILC  <- apply(data.y[,index.ILC], MARGIN=1, FUN=mean, na.rm=TRUE)

  color  <- rainbow(n.clst)
  par(mai=c(1,1,0.5,0.1),cex.axis=0.8, cex.lab=1,mgp=c(1.5,0.5,0))
  plot((mean.IDC+mean.ILC)/2, 
       (mean.IDC-mean.ILC), 
       xlab="(IDC mean + ILC mean)/2",
       ylab="IDC mean - ILC mean",
       pch=paste(fit1.u),
       col=color[fit1.u],
       main=paste("K=",n.clst))
}

\keyword{Functions}
