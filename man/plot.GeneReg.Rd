\name{plot.GeneReg}
\alias{plot.GeneReg}
\title{Plot gene regulatory network}
\description{
Plot gene regulatory network based on a series of time delay linear models.
}
\usage{
\method{plot}{GeneReg}(edge, \dots)
}
\arguments{
  \item{edge}{output of timedelay.lm.batch, which has 4 columns.
  The first column is the regulator, the second column is the target gene, 
  the third column is the regression coefficient, the fourth column is time delay.}
  \item{\dots}{Additional arguments, passed to plot.igraph.}
}
\references{
\cite{Tao Huang, Lei Liu, Ziliang Qian, Kang Tu, Yixue Li, Lu Xie: Using GeneReg to construct time delay gene regulatory networks. BMC Res Notes 2010, 3(1):142.}
\url{http://www.biomedcentral.com/1756-0500/3/142}
} 
