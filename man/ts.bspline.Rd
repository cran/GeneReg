\name{ts.bspline}
\alias{ts.bspline}
\title{B spline interpolation}
\description{
ts.bspline is for B spline interpolation.
}
\usage{
ts.bspline(expr,ts.point=NULL, data.predict=100, df=round(ncol(expr)*0.8))
}
\arguments{
  \item{expr}{an expression matrix with genes in the rows, time points in the columns}
  \item{ts.point}{time points}
  \item{data.predict}{the number of time points to be predicted}
  \item{df}{degrees of freedom, used in B-Spline Basis for Polynomial Splines}
}  
