\name{ts.bspline}
\alias{ts.bspline}
\title{B spline interpolation}
\description{
ts.bspline is for B spline interpolation.
}
\usage{
ts.bspline(expr,ts.point=NULL,df=round(ncol(expr)*0.75),data.predict=100)
}
\arguments{
  \item{expr}{an expression matrix with genes in the rows, time points in the columns}
  \item{ts.point}{time points}
  \item{df}{degrees of freedom, used in B-Spline Basis for Polynomial Splines}
  \item{data.predict}{the number of time points to be predicted}
}  
