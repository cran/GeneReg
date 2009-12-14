\name{timedelay.univariate.lm}
\alias{timedelay.univariate.lm}
\title{Internal function for time delay linear model}
\description{
timedelay.univariate.lm is an internal function for time delay linear model.
}
\usage{
timedelay.univariate.lm(bspline.data, target, regulator,
maxdelay=ncol(bspline.data)*0.25,
single.adj.r.squared=0.8, min.coef=0.25, max.coef=4)
}
\arguments{
  \item{bspline.data}{an interpolated expression matrix with genes in the rows, predicted time points in the columns}
  \item{target}{target gene}
  \item{regulator}{potential regulators}
  \item{maxdelay}{max time points delay}
  \item{single.adj.r.squared}{the cutoff of adjusted R squared of univariate regression}
  \item{min.coef}{the regression coefficient should be larger than min.coef}
  \item{max.coef}{the regression coefficient should be smaller than max.coef}  
}
