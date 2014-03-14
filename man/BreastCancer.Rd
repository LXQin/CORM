\name{BreastCancer}

\docType{data}

\alias{BreastCancer}

\title{Breast Cancer Data Set}

\description{Zhao et al. (2004) studied gene expression profiles of two types of breast cancer, invasive ductal carcinoma (IDC) and invasive lobular carcinoma (ILC). They analyzed the expression profiles of 38 IDC samples and 21 ILC samples, using cDNA arrays spotted for 42,000 clones. With the significance analysis of microarrays (SAM) method (Efron et al., 2001), they identified a total of 474 clones that were differentially expressed between IDCs and ILCs, representing 354 unique genes. This example data set contains the normalized expression value of these 354 genes for the 59 samples.
}

\usage{data(BreastCancer)}

\format{
 \enumerate{
  \item A list comprised of two components: \bold{normalizedData} and \bold{designMatrix}.
  \item \bold{normalizedData} is a matrix containing the normalized breast cancer data, whose row 
  names are gene IDs and column names indicate two subtypes of breast cancer (IDC vs. ILC).
  \item \bold{designMatrix} is the covariates matrix used to fit the clustering of
  linear models (CLM), whose row names are samples and column names are covariates.
 }
}

\references{
  Zhao et al. (2004).  Different gene expression patterns in invasive
  lobular and ductal carcinomas of the breast. \emph{Molecular Biology 
  of the Cell},15, 2523-2536
}

\keyword{Datasets}
