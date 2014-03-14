\name{fit.CLMM.2}

\alias{fit.CLMM.2}

\title{
  Clustering of Linear Mixed Models Method
}

\description{
  Fit a CLMM model for replicated time course data allowing for two 
  sets of time points among biological or technical replicates. No 
  missing value is allowed.
}

\usage{
  fit.CLMM.2(data.y1,data.x1,data.z1,data.y2,data.x2,data.z2,n.clst, n.run = 1)
}

\arguments{
  \item{data.y1}{
    a three dimensional array of gene expression data, data.y1[j, i, t] for gene j, sample i,  and time point t.
  }
  \item{data.x1}{
    a three dimensional array of fixed effects (common for all genes), data.x1[i, t, p] for sample i, time point t, and covariate p.
  }
  \item{data.z1}{
    a three dimensional array of random effects (common for all genes), data.z1[i, t, q] for sample i, time point t, and covariate p.
  }
  \item{data.y2}{
    a three dimensional array of gene expression data, data.y2[j, i, t] for gene j, sample i,  and time point t.
  }
  \item{data.x2}{
    a three dimensional array of fixed effects (common for all genes), data.x2[i, t, p] for sample i, time point t, and covariate p.
  }
  \item{data.z2}{
    a three dimensional array of random effects (common for all genes), data.z2[i, t, q] for sample i, time point t, and covariate p.
  }
  \item{n.clst}{
    an integer, number of clusters.
  }
  \item{n.run}{
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
  Li-Xuan Qin, Linda Breeden and Steven G. Self (2014). Finding gene clusters 
  for a replicated time course study. \emph{BMC Res Notes} 7:60.  
}

\seealso{
  \code{\link{fit.CLM}},
  \code{\link{fit.CLMM}},
  \code{\link{fit.CLMM.2}},
  \code{\link{fit.CLMM.NA}},
  \code{\link{fit.CLMM.NA.2}}
}

\keyword{Functions}

