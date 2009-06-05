timedelay.lm.batch <-
function(expr,regulator2target,...) {
	regulator2target<-regulator2target[(regulator2target[,1] %in% rownames(expr))&(regulator2target[,2] %in% rownames(expr)),]
	
	target<-as.character(unique(regulator2target[,2]))
	result<-NULL  	
	c<-0   	
	for (i in 1:length(target)) {
		if ((i*100)%/%length(target)>=c) {
			cat(c,'% done.\n',sep='')
			c<-c+1
		}
		regulator<- as.character(regulator2target[regulator2target[,2]==target[i],1]) 
		result.i<-timedelay.lm(expr,target[i],regulator,...)
		if (!is.null(result.i)) {
			result<-rbind(result,data.frame(regulator=names(result.i$delay),target=rep(target[i],length(result.i$delay)),coef=summary(result.i$fit)$coef[,1],delay=result.i$delay))
		}
	}
	
	rownames(result)<-NULL
	return(result)
}

