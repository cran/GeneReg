\name{plot}
\alias{plot}
\title{Plot gene regulatory network}
\description{
Plot gene regulatory network based on a series of time delay linear models.
}
\usage{
plot(edge, \dots)
}
\arguments{
  \item{edge}{output of timedelay.lm.batch, which has 4 columns.
  The first column is the regulator, the second column is the target gene, 
  the third column is the regression coefficient, the fourth column is time delay.}
  \item{\dots}{Additional arguments, passed to plot.igraph.}
}

