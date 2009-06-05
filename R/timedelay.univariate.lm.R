timedelay.univariate.lm <-
function(expr,target,regulator,maxdelay=ncol(expr)*0.25,univariate.adj.r.squared=0.8,min.coef=0.25) {
	y<-target.expr<-expr[target,]
	ts.point<-as.numeric(colnames(expr))
	finished=FALSE
	while (!finished) {
	max.adj.rsquare<-(-1)
	updated<-TRUE
	regulator.delay<-rep(0,length(regulator)); names(regulator.delay)<-regulator
	x<-regulator.expr<-expr[regulator,]
	if (length(regulator)==1) x<-regulator.expr<-t(as.matrix(regulator.expr))
	while (updated) {
		updated<-FALSE
		for (i in 1:nrow(x)) {
			regulator.delay.2<-regulator.delay
			for (delay in 0:maxdelay) {
				regulator.delay.2[i]<-delay
				
				x.to.fit<-matrix(0,length(regulator),length(y)-max(regulator.delay.2))
				rownames(x.to.fit)<-regulator
				for (j in 1:length(regulator)) {
					x.to.fit[j,]<-x[j,(max(regulator.delay.2)-regulator.delay.2[j]+1):(length(y)-regulator.delay.2[j])]
				}
				y.to.fit<-y[(max(regulator.delay.2)+1):length(y)]
				
				fit<-lm(y~.-1,data=data.frame(y=y.to.fit,t(x.to.fit)))
				adj.rsquare<-summary(fit)$adj.r.squared

				if (adj.rsquare>max.adj.rsquare) { 
					max.adj.rsquare<-adj.rsquare
					updated<-TRUE
					regulator.delay<-regulator.delay.2
					best.fit<-fit
					fit.coef<-summary(fit)$coef
				}   		
			}
		}
	}

	regulator2.ix<-which(abs(fit.coef[,1])>min.coef & abs(fit.coef[,1])<1/min.coef)
	finished<-(length(regulator)==length(regulator2.ix))
	if (length(regulator2.ix)==0) return(NULL)
	regulator<-regulator[regulator2.ix]
	}
	
	if (summary(best.fit)$adj.r.squared<=univariate.adj.r.squared) return(NULL)  
	return(list(delay=regulator.delay*(ts.point[2]-ts.point[1]),fit=best.fit))
}

