\name{timedelay.univariate.lm}
\alias{timedelay.univariate.lm}
\title{Internal function for time delay linear model}
\description{
timedelay.univariate.lm is an internal function for time delay linear model.
}
\usage{
timedelay.univariate.lm(expr,target,regulator,maxdelay=ncol(expr)*0.25,
univariate.adj.r.squared=0.8,min.coef=0.25)
}
\arguments{
  \item{expr}{an expression matrix with genes in the rows, time points in the columns}
  \item{target}{target gene}
  \item{regulator}{potential regulators}
  \item{maxdelay}{max time points delay}
  \item{univariate.adj.r.squared}{the cutoff of adjusted R squared of univariate regression}
  \item{min.coef}{the regression coefficient should be larger than min.coef and smaller than 1/min.coef}
}
