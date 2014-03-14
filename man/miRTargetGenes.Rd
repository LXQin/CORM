\name{miRTargetGenes}

\docType{data}

\alias{miRTargetGenes}

\title{miR-let-7f Targets Data Set}

\description{Lu et al. (2005) profiled both miRNA expression and mRNA expression 
in multiple human cancer types. The miRNA profiles reflected the developmental lineage 
and differentiation state of the tumours. One of the profiled miRNA is let-7f. 
Its expression is highly correlated (Pearson correlation < -0.458) with the expression 
of 178 genes in tumors from 5 cancer types.
}

\usage{data(miRTargetGenes)}

\format{
  A list comprised of two components: \bold{normalizedData} and \bold{designMatrix}.
  \enumerate{
  \item \bold{normalizedData} is a matrix containing the normalized data of the miR-let-7f targets, whose row 
  names are gene IDs and column names indicate normal and tumor samples of 5 types of cancer.
  \item \bold{designMatrix} is the covariates matrix used to fit the clustering of
  linear models (CLM), whose row names are samples and column names are covariates.
 }
}

\references{
  Lu et al. (2005). MicroRNA expression profiles classify human cancers. 
  \emph{Nature}, 435, 834-838
}

\keyword{Datasets}
