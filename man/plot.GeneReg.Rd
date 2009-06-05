\name{plot.GeneReg}
\alias{plot.GeneReg}
\title{Plot gene regulatory network}
\description{
Plot gene regulatory network based on time delay linear model.
}
\usage{
\method{plot}{GeneReg}(edge, \dots)
}
\arguments{
  \item{edge}{output of timedelay.lm.batch or timedelay.lm, which has 4 columns.
  The first column is the regulator, the second column is the target gene, 
  the third column is the regression coefficient, the fourth column is time delay.}
  \item{\dots}{Additional arguments, passed to plot.igraph.}
}

