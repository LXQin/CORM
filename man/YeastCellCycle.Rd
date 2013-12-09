\name{YeastCellCycle}

\docType{data}

\alias{YeastCellCycle}

\title{Yeast Cell Cycle Data Set}

\description{
This data set contains a subset of yeast cell cycle data taken from Spellman et al. (1998) (See the reference below). Spellman et al. (1998) monitored the genome-wide mRNA levels for 6108 yeast genes at 7-minute intervals for 119 minutes. A total of 256 genes were identified to oscillate significantly in at least two data sets. This example data set contains the log ratios of these 256 genes at the first 16 time points(from 0 min to 105 min).
}

\usage{data(YeastCellCycle)}

\format{
  A list comprised of three components: normalizedData, normalizedData.sample and designMatrix.
  normalizedData is a three dimensional array containing the normalized expression data of the 256 genes during yeast cell-cycle. 
  normalizedData.sample is a randomly selected sample from normalizedData. It only contains 64 genes.
  designMatrix is a three dimensional array used to fit the clustering of linear mixed models (CLMM).
}

\references{
  Spellman et al. (1998). Comprehensive identification of cell cycle-regulated
   genes of the yeast Saccharomyces cerevisiae by microarray hybridization. 
  \emph{Molecular Biology of the Cell},9, 3273-3297
}

\source{
  http://genome-www.stanford.edu/cellcycle/data/rawdata/
}

\keyword{Datasets}
