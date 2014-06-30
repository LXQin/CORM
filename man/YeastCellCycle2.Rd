\name{YeastCellCycle2}
\docType{data}
\alias{YeastCellCycle2}
\title{Yeast Cell Cycle Data Sets of wild type and single mutant}
\description{
This data set is taken from a yeast cell cycle study in Qin et al. (2014) 
(See the reference below). It contains the normalized data of 256 cell cycle 
dependent genes in wild type(WT) and single mutant(SM) yeast, each 
mesured over time with two technical replicates.
}
\usage{data(YeastCellCycle2)}
\format{
  A list comprised of four components: \bold{normalizedData.WT}, \bold{normalizedData.SM}, \bold{designMatrix.WT}, and \bold{designMatrix.SM}.
 \enumerate{
  \item \bold{normalizedData.WT} and \bold{normalizedData.SM} are three-dimensional 
  arrays containing the normalized expression data for 256 genes in WT and SM yeast, respectively. Missing values are indicated by "NA". 
  \item \bold{designMatrix.WT} and \bold{designMatrix.SM} are three-dimensional arrays
   used to fit the clustering of linear mixed models (CLMM).
 }
}
\references{
Li-Xuan Qin, Linda Breeden and Steven G. Self (2014). 
Finding gene clusters for a replicated time course study. 
\emph{BMC Res Notes} 7:60.
}
\keyword{Datasets}
