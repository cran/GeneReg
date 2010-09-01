\name{timedelay.lm.batch}
\alias{timedelay.lm.batch}
\title{time delay linear model batch}
\description{
Batch usage of time delay linear model.
}
\usage{
timedelay.lm.batch(bspline.data, expr.data, regulator.list, 
target.list=rownames(bspline.data), \dots)
}
\arguments{
  \item{bspline.data}{output of ts.bspline, an interpolated expression matrix with genes in the rows, predicted time points in the columns}
  \item{expr.data}{an observed expression matrix with genes in the rows, predicted time points in the columns}
  \item{regulator.list}{a candidate pool of potential regulators}
  \item{target.list}{the list of target genes that will be investigated}
  \item{ \dots}{Additional arguments, passed to \link{timedelay.lm}}
}
\references{
\cite{Tao Huang, Lei Liu, Ziliang Qian, Kang Tu, Yixue Li, Lu Xie: Using GeneReg to construct time delay gene regulatory networks. BMC Res Notes 2010, 3(1):142.}
\url{http://www.biomedcentral.com/1756-0500/3/142}
} 
