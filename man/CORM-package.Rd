\name{CORM}

\alias{CORM}

\docType{package}

\title{
  Clustering of Regression Models Method
}

\description{
  The clustering of regression models method 
  clusters genes based on the regression coefficients.
}

\details{
  \tabular{ll}{
    Package: \tab CORM\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2013-11-18\cr
    License: \tab GPL(>=2)\cr
  }
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

\keyword{Introduction}

\seealso{
  \code{\link{fit.CLM}},
  \code{\link{fit.CLMM}}
}

\examples{
  data(BreastCancer)
  data.y <- BreastCancer$normalizedData
  data.x <- BreastCancer$designMatrix
  n.clst <- 9
  fit1   <- fit.CLM(data.y, data.x, n.clst)
}
