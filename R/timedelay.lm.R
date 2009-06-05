timedelay.lm <-
function(expr,target,regulator,maxdelay=ncol(expr)*0.25,univariate.adj.r.squared=0.8,adj.r.squared=0.9,min.coef=0.25,
output=FALSE,topdf=FALSE,xlab='time point',ylab='log ratio') {

	if (length(regulator)>1) {
		regulator.2<-NULL
		for (i in 1:length(regulator)) {
			regulator.i<-timedelay.univariate.lm(expr,target,regulator[i],maxdelay=maxdelay,univariate.adj.r.squared=univariate.adj.r.squared,min.coef=min.coef)  	
			if (!is.null(regulator.i)) {
				regulator.2<-c(regulator.2,regulator[i])
			}
		}
		regulator<-regulator.2
	}
	if (length(regulator)==0) return(NULL)

	y<-target.expr<-expr[target,]
	
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
	
	if (summary(best.fit)$adj.r.squared<=adj.r.squared) return(NULL)
	
	ts.point<-as.numeric(colnames(expr))
	if (output) {
		if (topdf) {
			pdf(paste(target,'.pdf',sep=''),width=8,height=6)
		} else {
			x11()
		}
		x.to.fit<-matrix(0,length(regulator),length(y)-max(regulator.delay))
		for (j in 1:length(regulator)) {
			x.to.fit[j,]<-x[j,(max(regulator.delay)-regulator.delay[j]+1):(length(y)-regulator.delay[j])]*fit.coef[j,1]
		}
		pal<-rainbow(round(length(regulator)*1.33))
		plot(ts.point[1:length(y)],y,type='l',col='black',ylim=c(min(y,x.to.fit,x),max(y,x.to.fit,x)),lwd=4,xlab=xlab,ylab=ylab)
		title(paste('adj.r.squared =',summary(best.fit)$adj.r.squared))
		y.to.fit<-y[1:(length(y)-max(regulator.delay))]
		lines(ts.point[(max(regulator.delay)+1):length(y)],predict(best.fit),lty='dashed',col='black',lwd=4)
		for (j in 1:length(regulator)) {
			lines(ts.point[1:length(y)],x[j,],col=pal[j],lwd=2)
			lines(ts.point[(max(regulator.delay)+1):length(y)],x.to.fit[j,],col=pal[j],lty='dashed',lwd=2)
		}
		legend.text<-c(target,paste('predicted',target))
		legend.col<-rep('black',2)
		legend.lty<-rep(c('solid','dashed'),length(regulator)+1)
		legend.lwd<-c(2,2,rep(2,length(regulator)*2))
		for (i in 1:length(regulator)) {
			legend.text<-c(legend.text,paste(regulator[i],': ',format(summary(best.fit)$coef[i,1],digits=3),sep=''),paste('delay: ',format(regulator.delay[i]*(ts.point[2]-ts.point[1]),digits=3)))
			legend.col<-c(legend.col,rep(pal[i],2))
		}
		legend('topleft',legend=legend.text,col=legend.col,lty=legend.lty,lwd=legend.lwd)
		if (topdf) dev.off()
	}
	
	return(list(delay=regulator.delay*(ts.point[2]-ts.point[1]),fit=best.fit))
}

