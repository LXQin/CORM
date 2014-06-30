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
    an integer used to get the starting value for the EM algorithm.
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
This function implements the Clustering of Linear Models Method of Qin and Self (2006). This method clusters genes based on the estimated regression parameters that model the relation between gene expression and sample covariates. 
}

\author{ Li-Xuan Qin 
  \email{qinl@mskcc.org}
}

\references{
 \itemize{
  \item Li-Xuan Qin and Steven G. Self (2006).  The clustering of regression models
  method with applications in gene expression data.  \emph{Biometrics}
  62, 526-533. 
  \item Li-Xuan Qin (2008). An integrative analysis of microRNA and mRNA expression - 
  a case study. \emph{Cancer Informatics} 6, 369-379.
 }
}

\seealso{
  \code{\link{fit.CLM}},
  \code{\link{fit.CLMM}},
  \code{\link{fit.CLMM.2}}
}

\examples{
#Example 1
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
 
\dontrun{
#Example 2
 #test data
  data(miRTargetGenes)
  data.y <- miRTargetGenes$normalizedData
  data.x <- miRTargetGenes$designMatrix
 #fit the model
  n.clst <- 9
  n.start<- 20
  fit2  	 <- fit.CLM(data.y, data.x, n.clst, n.start)
  fit2.u   <- apply(fit2$u.hat, MARGIN=1, FUN=order, decreasing=TRUE)[1,]
  fit2.u.o <- factor(fit2.u, levels=c(1,5,6,7,4,8,2,9,3), labels=1:9)
  library(limma)
  plot.y   <- lmFit(data.y, data.x)$coef \%*\% cbind(c(1,0,0,0),c(1,0,1,0),c(1,1,0,0),c(1,1,1,1))
  plot.x   <- 1:4
 #display the results
  color		 <- rainbow(n.clst)
  par(mfrow=c(3,4),mai=c(0.35, 0.4, 0.4, 0.2), mgp=c(1.6,0.4,0), tck=-0.01, las=2)
  for(k in 1:n.clst){
   plot(plot.x, plot.y[1,], type="n", xaxt="n", ylim=range(plot.y), 
        xlab="", ylab="gene expression")
   axis(1, plot.x, c("Normal \n","Normal \n +miRNA","Tumor \n","Tumor \n +miRNA"), 
        las=1, cex.axis=1, mgp=c(1.5,1.2,0))
   title(paste("cluster", k))
   abline(h=0, lty=2)
   for(j in which(fit2.u.o==k)) points(plot.x, plot.y[j,], type="b", col=color[k])
  }
}
}

\keyword{Functions}
