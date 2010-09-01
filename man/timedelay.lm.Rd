\name{timedelay.lm}
\alias{timedelay.lm}
\title{time delay linear model}
\description{
timedelay.lm is for time delay linear model.
}
\usage{
timedelay.lm(bspline.data, expr.data, target, regulator,   
maxdelay=ncol(bspline.data)*0.25, single.adj.r.squared=0.8,
multiple.adj.r.squared=0.9, min.coef=0.25, max.coef=4,
output=FALSE, topdf=FALSE, xlab='time point', ylab='log ratio')
}
\arguments{
  \item{bspline.data}{an interpolated expression matrix with genes in the rows, predicted time points in the columns}
  \item{expr.data}{an observed expression matrix with genes in the rows, predicted time points in the columns}
  \item{target}{target gene}
  \item{regulator}{potential regulators}
  \item{maxdelay}{max time points delay}
  \item{single.adj.r.squared}{the cutoff of adjusted R squared of univariate regression}
  \item{multiple.adj.r.squared}{the cutoff of adjusted R squared of multivariate regression}
  \item{min.coef}{the regression coefficient should be larger than min.coef}
  \item{max.coef}{the regression coefficient should be smaller than max.coef}  
  \item{output}{whether the model should be plotted}
  \item{topdf}{whether the model should be plotted into a pdf file }
  \item{xlab}{xlab of the plot}
  \item{ylab}{ylab of the plot}
}
\references{
\cite{Tao Huang, Lei Liu, Ziliang Qian, Kang Tu, Yixue Li, Lu Xie: Using GeneReg to construct time delay gene regulatory networks. BMC Res Notes 2010, 3(1):142.}
\url{http://www.biomedcentral.com/1756-0500/3/142}
} 
