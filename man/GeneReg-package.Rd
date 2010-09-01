\name{GeneReg-package}
\Rdversion{1.1.2}
\alias{GeneReg}
\docType{package}
\title{
Construct time delay gene regulatory network
}
\description{
\pkg{GeneReg} is an R package for inferring time delay gene regulatory network 
using time course gene expression data.
}
\usage{
}
\arguments{
}
\details{

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
The regression models were optimized by the AIC forward selection procedure.
\cr \cr
In the following example, we will demonstrate how to use this package to 
analysis a cell cycle datasets which is publicly available at GEO 
\url{http://www.ncbi.nlm.nih.gov/geo} under accession number GSE8799. 
The data set included the gene expression profiles of wild-type cells and 
cyclin-mutant cells with 15 time points during two cell cycles. 
1271 periodic genes acquired from wild-type cells were analyzed to investigate 
how yeast regulators regulate target gene expression. These periodic genes formed 
the list of target genes. Then, a candidate pool of potential regulators 
of target genes was constructed by intersecting the periodic genes and 
the transcription factors in YEASTRACT \url{http://www.yeastract.com/}.   
\cr \cr
First, B spline interpolation was applied to estimate 100 time points 
according to the original 15 time points. 
\cr \cr
\code{> data(wt.expr.data)}
\cr \cr
\code{> wt.bspline.data<- ts.bspline(wt.expr.data,  }
\cr \cr
\code{+ ts.point=as.numeric(colnames(wt.expr.data)), data.predict=100)}
\cr \cr
Then, the time delay linear model was carried out based on the interpolated expression data. 
single.adj.r.squared and multiple .adj.r.squared were the cutoffs of adjusted R seqaured of 
single regulator regression and multiple regulator regression, respectively. 
tfs.list were the candidate pool of regulators. 
min.coef means the regression coefficient should be larger than 1/4 and smaller than 4. 
maxdelay means the time points delay should not exceed 1/4 of all time points.
\cr \cr
\code{> data(tf.list)}
\cr \cr
\code{> dir.create('wt.model')}
\cr \cr
\code{> setwd('wt.model')}
\cr \cr
\code{> wt.models<-timedelay.lm.batch(bspline.data=wt.bspline.data,  }
\cr \cr
\code{+ expr.data=wt.expr.data, regulator.list=tf.list,  }
\cr \cr
\code{+ target.list=rownames(wt.bspline.data),  }
\cr \cr
\code{+ single.adj.r.squared=0.8, multiple.adj.r.squared=0.9, }
\cr \cr
\code{+ maxdelay=ncol(wt.bspline.data)*0.1, min.coef=0.25, max.coef=4,  }
\cr \cr
\code{+ output=T, topdf=T, xlab='Time point (lifeline)',  }
\cr \cr
\code{+ ylab='Relative expression level (in log ratio)')}
\cr \cr
\code{> setwd('..')}
\cr \cr
Finally, the gene regulatory network can be plotted according to 
the series of time delay linear models.
\cr \cr
\code{> pdf('wt.network.pdf',width=70, height=70)}
\cr \cr
\code{> plot.GeneReg(wt.models,vertex.size=2,layout=layout.fruchterman.reingold)}
\cr \cr
\code{> dev.off()}
\cr \cr
The analysis of cyclin-mutant cells is similar.
}
\author{
Tao Huang <\email{huangtao@sibs.ac.cn}>
}
\references{
\cite{Tao Huang, Lei Liu, Ziliang Qian, Kang Tu, Yixue Li, Lu Xie: Using GeneReg to construct time delay gene regulatory networks. BMC Res Notes 2010, 3(1):142.}
\url{http://www.biomedcentral.com/1756-0500/3/142}
}  
