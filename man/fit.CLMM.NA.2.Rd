\name{fit.CLMM.NA.2}

\alias{fit.CLMM.NA.2}

\title{
  Clustering of Linear Mixed Models Method
}

\description{
  Fit a CLMM model for replicated time course data allowing for 
  missing value and allowing for two sets of time points among 
  the biological or technical replicates.
}

\usage{
  fit.CLMM.NA.2(data.y1,data.x1,data.z1,data.y2,data.x2,data.z2,
  na.flag1,na.flag2,n.clst, n.run = 1)
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
  \item{na.flag1}{
    NA status indicator for data.y1.
  }
  \item{na.flag2}{
    NA status indicator for data.y2.
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
 data.y1 <- YeastCellCycle2$normalizedData.WT
 data.x1 <- YeastCellCycle2$designMatrix.WT
 na.flag1<- YeastCellCycle2$NA.flag.WT
 data.y2 <- YeastCellCycle2$normalizedData.SM
 data.x2 <- YeastCellCycle2$designMatrix.SM
 na.flag2<- YeastCellCycle2$NA.flag.SM
 n.clst  <- 8
 fit1	<- fit.CLMM.NA.2(data.y1,data.x1,data.x1,data.y2,data.x2,data.x2,
                       na.flag1,na.flag2,n.clst)
 fit1.u  <- apply(fit1$u.hat, MARGIN=1, FUN=order, decreasing=TRUE)[1,]
 zeta.fitted <- fit1$theta.hat$zeta.hat
 profile.WT  <- YeastCellCycle2$designMatrix.WT[1,,] \%*\% t(zeta.fitted)
 profile.SM  <- YeastCellCycle2$designMatrix.SM[1,,] \%*\% t(zeta.fitted)
# remove bad time points for WTs and SMs
 n.time  <- 25
 time.WT <- (1:n.time)[-22]
 time.SM <- (1:n.time)[-c(6,9,12)]
 n.rpl.WT<- length(time.WT)
 n.rpl.SM<- length(time.SM)
 n.gene.short<-dim(YeastCellCycle2$normalizedData.WT)[1]
# gene-specific profile: observed profiles averaged over replicates
 data.WT.mean  <- matrix(0, nrow=n.gene.short, ncol=n.rpl.WT)
 data.SM.mean	<- matrix(0, nrow=n.gene.short, ncol=n.rpl.SM)
 for(j in 1:n.gene.short){
  data.WT.mean[j,] <- (YeastCellCycle2$normalizedData.WT[j,1,]+
                       YeastCellCycle2$normalizedData.WT[j,2,])/2
  data.SM.mean[j,] <- (YeastCellCycle2$normalizedData.SM[j,1,]+
                       YeastCellCycle2$normalizedData.SM[j,2,])/2
 }
# plot observed profiles by cluster -- wild type
 col.panel=rainbow(8)
 par(mai=c(0.3,0.25,0.2,0.05),mfrow=c(3,3))
 for(k in 1:n.clst){
  plot(5*(time.WT-1), profile.WT[,k], type="l", col=col.panel[k], ylim=c(-2,2),
       xlab="Time", ylab="Expression Value", main=paste("WT: cluster",k))
  i.plot <- (1:n.gene.short)[fit1.u==k]
  for(j in i.plot) lines(5*(time.WT-1), data.WT.mean[j,], lty=1)
  lines(5*(time.WT-1), profile.WT[,k], col=col.panel[k], lwd=2)
  text(125, -1.9, pos=2, paste(length(i.plot)," genes"))
}
# plot observed profiles by cluster -- single mutant
 par(mai=c(0.3,0.25,0.2,0.05),mfrow=c(3,3))
 for(k in 1:n.clst){
  plot(5*(time.SM-1), profile.SM[,k], type="l", col=col.panel[k], ylim=c(-2,2),
       xlab="Time", ylab="Expression Value", main=paste("SM: cluster",k))
  i.plot <- (1:n.gene.short)[fit1.u==k]
  for(j in i.plot) lines(5*(time.SM-1), data.SM.mean[j,], lty=1)
  lines(5*(time.SM-1), profile.SM[,k], col=col.panel[k], lwd=2)
  text(125, -1.9, pos=2, paste(length(i.plot)," genes"))
}
# plot fitted profiles by cluster
 par(mai=c(0.3,0.25,0.2,0.05),mfrow=c(3,3))
 for(k in 1:n.clst){ 
  plot(5*(time.WT-1), profile.WT[,k], type="l", ylim=c(-2,2), 
       xlab="Time", ylab="Expression Value", lwd=2)
  title(paste("Cluster", k))
  lines(5*(time.SM-1), profile.SM[,k], lty=3, lwd=2)
  if(k==1) legend(60, 2, c("WT", "SM"), lty=1:2, cex=0.8)
 }
}
}

\keyword{Functions}

