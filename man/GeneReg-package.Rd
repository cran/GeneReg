\name{GeneReg-package}
\alias{GeneReg}
\title{
GeneReg
}
\description{
\pkg{GeneReg} is an R package for inferring gene regulatory network with time delay 
using time course gene expression data.
}
\usage{
}
\arguments{
}
\details{
\tabular{ll}{
Package: \tab \pkg{GeneReg}\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-06-05\cr
License: \tab LGPL\cr
}
The goal of time delay linear model is to fit a linear regression model 
using a set of putative regulators to estimate the transcription pattern 
of a specific target gene. 
\cr \cr
\eqn{Eg = a1*Etf1(-t1) + a2*Etf2(-t2) + ... + ai*Etfi(-ti) + ... + an*Etfn(-tn)}
\cr \cr
Where Eg is the relative expression level of target gene, 
Etfi is the relative expression level of TFi, 
ti is time delay of TFi, 
ai is the regression coefficient of TFi.
\cr \cr
In the following example, we will demonstrate how to use this package to 
analysis a cell cycle datasets which is publicly available at GEO 
\url{http://www.ncbi.nlm.nih.gov/geo} under accession number GSE8799. 
The example time course gene expression profiles included 15 time points 
at 16 min resolution in wild type Saccharomyces cerevisiae. 
YEASTRACT (\url{http://www.yeastract.com/}) were used as potential gene regulation relationship.
\cr \cr
First, B spline interpolation was applied to estimate 100 time points 
according to the original 15 time points.
\cr \cr
\code{> data(expr)}
\cr \cr
\code{> expr<- log2(expr)}
\cr \cr
\code{> expr<- expr - expr[,1]}
\cr \cr
\code{> expr.bspline<- ts.bspline(expr, ts.point= as.numeric(colnames(expr)), }
\cr \cr
\code{+ data.predict=100)}
\cr \cr
Then, the time delay linear model was carried out based on the interpolated 
expression data and potential gene regulation relationship derived from YEASTRACT. 
univariate. adj.r.squared and adj.r.squared is the cutoff of Adjusted R2 of 
univariate regression and multivariate regression, separately. 
min.coef means the regression coefficient should be larger than 1/4 and smaller than 4. 
maxdelay means the time points delay should not exceed 1/4 of all time points.
\cr \cr
\code{> data(yeastract)}
\cr \cr
\code{> dir.create('model')}
\cr \cr
\code{> setwd('model')}
\cr \cr
\code{> models<-timedelay.lm.batch(expr.bspline, yeastract, }
\cr \cr
\code{+ univariate.adj.r.squared=0.8, adj.r.squared=0.9, }
\cr \cr
\code{+ min.coef=1/4, maxdelay=ncol(expr.bspline)/4, output=T,  }
\cr \cr
\code{+ topdf=T, xlab='time (minute)',ylab='log ratio')}
\cr \cr
\code{setwd('..')}
\cr \cr
Finally, the gene regulatory network can be plotted according to 
the time delay linear model of every target gene.
\cr \cr
\code{> pdf('network.pdf',width=70, height=70)}
\cr \cr
\code{> plot.GeneReg(model,vertex.size=2,layout=layout.fruchterman.reingold)}
\cr \cr
\code{> dev.off()}

}
\author{
Tao Huang

Maintainer: Tao Huang <\email{huangtao@sibs.ac.cn}>
}
\references{
\cite{Tao Huang, Ziliang Qian, Lei Liu, Yixue Li. GeneReg: R package for inferring gene regulatory network using time course gene expression data. 2009}
}  
