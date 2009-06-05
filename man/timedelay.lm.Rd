\name{timedelay.lm}
\alias{timedelay.lm}
\title{time delay linear model}
\description{
timedelay.lm is for time delay linear model.
}
\usage{
timedelay.lm(expr,target,regulator,maxdelay=ncol(expr)*0.25,
univariate.adj.r.squared=0.8,adj.r.squared=0.9,
min.coef=0.25,output=FALSE,topdf=FALSE,
xlab='time point',ylab='log ratio')
}
\arguments{
  \item{expr}{an expression matrix with genes in the rows, time points in the columns}
  \item{target}{target gene}
  \item{regulator}{potential regulators}
  \item{maxdelay}{max time points delay}
  \item{univariate.adj.r.squared}{the cutoff of adjusted R squared of univariate regression}
  \item{adj.r.squared}{the cutoff of adjusted R squared of multivariate regression}
  \item{min.coef}{the regression coefficient should be larger than min.coef and smaller than 1/min.coef}
  \item{output}{whether the model should be plotted}
  \item{topdf}{whether the model should be plotted into a pdf file }
  \item{xlab}{xlab of the plot}
  \item{ylab}{ylab of the plot}
}

