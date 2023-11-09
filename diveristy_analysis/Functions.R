# you will need this for the alpha and beta diversity code.
#Code From Bolker 2009 and Fox 2002:

# Function to check for overdispersion (from http://glmm.wikidot.com/faq)
	overdisp_fun <- function(model) {
	## number of variance parameters in 
	##   an n-by-n variance-covariance matrix
	  vpars <- function(m) {
	    nrow(m)*(nrow(m)+1)/2
	  }
	  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
	  rdf <- nrow(model.frame(model))-model.df
	  rp <- residuals(model,type="pearson")
	  Pearson.chisq <- sum(rp^2)
	  prat <- Pearson.chisq/rdf
	  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
	  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
	}
	

#Function To Calculate Ratio of Model SS to Mean Parametric Bootstrap SS ('bias') from Harrison 2014 Peer J: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4194460/#supp-1
od<-function(bootobject){
	biasvals<-bootobject $t0/bootobject[2]$t
	bias<-mean(biasvals,na.rm=T)
	intervals<-quantile(biasvals,c(0.025,0.975),na.rm=T)
	dat<-c(bias,intervals)
	return(dat)
}


tryCatch.W.E <- function(expr){
     W <- NULL
     .number_of_warnings <- 0L
     frame_number <- sys.nframe()
     
     w.handler <- function(w){ # warning handler
       assign(".number_of_warnings",.number_of_warnings+1L,
              envir = sys.frame(frame_number))
       W <<- w
       invokeRestart("muffleWarning")
     }
     list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
 				     warning = w.handler),
          warning.no = .number_of_warnings,
 	        warning = W)
}

	
	
# Function to diagnose Poisson fits with overdispersion (from Dr. Paul Johnson, https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q3/020820.html)
# St dev used in Residuals discussed here: http://abdn.ac.uk/lambin-group/Papers/Paper%2041%202001%20Elston%20Tick%20aggregation%20Parasitology.pdf	
  plot.glmer<-function(mer.fit,type="pearson",overdispersion.term=NULL,title=NULL,...){
		require(AICcmodavg)
		if(is.null(overdispersion.term)){ 
			Fitted<-fitted(mer.fit)
			Residuals=resid(mer.fit,type)
		} else {
			response<-model.frame(mer.fit)[[1]]
			od.ranef<-lme4::ranef(mer.fit)[[overdispersion.term]][[1]]
			if(length(response)!=length(od.ranef) || fam.link.mer(mer.fit)$family!="poisson" || fam.link.mer(mer.fit)$link!="log")
			stop("Model is not lognormal-Poisson. Cannot use overdispersion term.")
			Fitted<-exp(log(fitted(mer.fit))-od.ranef)
			Residuals<-(response - Fitted)/sqrt(Fitted+(Fitted^2)*c(exp(lme4::VarCorr(mer.fit)[[overdispersion.term]])-1))
		}
		plot.data<-data.frame(Fitted=Fitted,Residuals=Residuals)
		plot.data$loess.line<-predict(loess(Residuals~Fitted,data=plot.data))
		plot.data<-plot.data[order(plot.data$Fitted),]
		sdabove = mean(Residuals)+sd(Residuals)*2.5
    sdbelow = mean(Residuals)-sd(Residuals)*2.5
    plot(plot.data[,c("Fitted","Residuals")],main=title,ylim=range(sdbelow,sdabove))
		abline(h=0)
    abline(h=sdabove,col="red",lty=2)
		abline(h=sdbelow,col="red",lty=2)
		points(plot.data[,c("Fitted","loess.line")],type="l",col="red")
		hist(plot.data$Residuals,xlab="Residuals",main="")
		}

# Function to detect outlier for glmm based on above code
outlier <- function(mer.fit,overdispersion.term=NULL,type="pearson"){
	require(AICcmodavg)
      		if(is.null(overdispersion.term)){ 
			      Fitted<-fitted(mer.fit)
			      Residuals=resid(mer.fit,type)
			      above=which((mean(Residuals)+sd(Residuals)*2.5)<Residuals)
      		  below=which((mean(Residuals)-sd(Residuals)*2.5)>Residuals)
      		  return(outliers=rbind(above=above,below=below))
		} else {      	
      	response<-model.frame(mer.fit)[[1]]
  		  od.ranef<-lme4::ranef(mer.fit)[[overdispersion.term]][[1]]
			  if(length(response)!=length(od.ranef) || fam.link.mer(mer.fit)$family!="poisson" || fam.link.mer(mer.fit)$link!="log")
			  stop("Model is not lognormal-Poisson. Cannot use overdispersion term.")
			  Fitted<-exp(log(fitted(mer.fit))-od.ranef)
			  Residuals<-(response - Fitted)/sqrt(Fitted+(Fitted^2)*c(exp(lme4::VarCorr(mer.fit)[[overdispersion.term]])-1))
      	above=which((mean(Residuals)+sd(Residuals)*2.5)<Residuals)
      	below=which((mean(Residuals)-sd(Residuals)*2.5)>Residuals)
      	return(outliers=rbind(above=above,below=below))
    }
    }

	
# Coefficients from glm (from Bolker 2009)
	plot.lmList <- function(object, ord.var = 1, ...) {
		require(reshape)
		cL <- coef(object)
		cL <- data.frame(grp = rownames(cL), cL)
		if (is.numeric(ord.var) & length(ord.var) == 1) {
			ord.var <- names(cL)[ord.var + 1]
			cL$grp <- reorder(cL$grp, cL[[ord.var]])
		}
		else cL$grp <- reorder(cL$grp, ord.var)
		dotplot(grp ~ value | variable, data = melt(cL), ...)
	}

# Assessing normality assumption of random effects (from Bolker 2009)
	qqmath.lmList <- function(object, ...) {
		require(reshape)
		qqmath(~value | variable, data = melt(coef(object)),
		prepanel = prepanel.qqmathline,
			panel = function(x, ...) {
	    	panel.qqmathline(x, ...)
    		panel.qqmath(x, ...)
		}, scale = list(y = list(relation = "free")))
}

plot_diversitytime<-function(index,data=tmp,whichones="all"){ 
# Plots diversity indices over time for each diet and for each individual within a diet
	for (div in index){

	if (sum(grepl(paste0(whichones,"+"),c("all","summary")))){
	par(mfrow=c(3,1),mai=c(0.75,0.75,0.5,0.5))
	
	for (dt in c("DF","HN","LN")){	
		plot(data[[div]]~week,type="n",data=data,subset=diet==dt[1]&ind==1,ylim=c(0.9*min(data[[div]]),1.1*max(data[[div]])),xlim=c(0.5,8.5),main=paste(div,"on Diet",dt),cex.lab=2,cex.main=2,cex=3,ylab=div,cex.axis=1.5)
			for (i in 1:max(unique(data[data$diet==dt,ind]))){
				lines(data[[div]]~week,type="o",pch=20,bg=i,data=data,col=i,subset=diet==dt&ind==i,cex=3)
			}
	}
	}

	if (sum(grepl(paste0(whichones,"+"),c("all","DF")))) print(xyplot(data[diet=="DF"][[div]]~week|individual,data=data[diet=="DF"],type="o",xlim=c(0,9),main=paste(div,"on Diet DF"),ylab=div))

	if (sum(grepl(paste0(whichones,"+"),c("all","HN")))) print(xyplot(data[diet=="HN"][[div]]~week|individual,data=data[diet=="HN"],type="o",xlim=c(0,9),main=paste(div,"on Diet HN"),ylab=div))
	
	if (sum(grepl(paste0(whichones,"+"),c("all","LN")))) print(xyplot(data[diet=="LN"][[div]]~week|individual,data=data[diet=="LN"],type="o",xlim=c(0,9), main=paste(div,"on Diet LN"),ylab=div))
	
	# Another way to make the same xyplots	
		# DF_grouped=groupedData(S~week|individual,data=diversity[diet=="DF"])
		# plot(DF_grouped)
		
	}	
	}
	
	
plot_coef<-function(data=diversity){	
# This is a preliminary scan of coefficient estimates and check of assumptions, and assumptions will need to be assessed with final model.

	##### Box plots following the Fox example ##########	
	# Couldn't get to work with div with lmList...
	# lmList produces individual lm fits for each level of grouping factor, in this case individual

# S
	df.list <- lmList(S ~ week | individual, subset = diet=="DF", data=data)
	hn.list <- lmList(S ~ week | individual, subset = diet=="HN", data=data)
	ln.list <- lmList(S ~ week | individual, subset = diet=="LN", data=data)

	#lmList works for nlme if week is numeric (doesn't work if it's a factor/character)
	df.coef <- coef(df.list)
	hn.coef <- coef(hn.list)
	ln.coef <- coef(ln.list)
	
	old <- par(mfrow=c(1,2))
	boxplot(df.coef[,1], hn.coef[,1], ln.coef[,1], main=paste("Intercepts of S, Family = Normal"),names=c("DF", "HN","LN"))
	boxplot(df.coef[,2], hn.coef[,2], ln.coef[,2], main=paste("Slopes of S, Family = Normal"),names=c("DF", "HN","LN"))
	par(old)
	rm(df.list,hn.list,ln.list,df.coef,hn.coef,ln.coef)
	
	# From Bolker 2009 GLMM for Ecologists
	data$wd<-as.factor(data$week):data$diet
	data$wd <- with(data, reorder(wd, S, mean))
	print(bwplot(S ~ wd, data = data, scales = list(x = list(rot = 90)),main="Variances Within Weeks and Diets"))
		
	data$log.S<-log(data$S)
	data$wd<-as.factor(data$week):data$diet
	data$wd <- with(data, reorder(wd, log.S, mean))
	print(bwplot(log.S ~ wd, data = data, scales = list(x = list(rot = 90)),main="Variances Within Weeks and Diets"))
	

	data$sqrt.S<-sqrt(data$S)
	data$wd<-as.factor(data$week):data$diet
	data$wd <- with(data, reorder(wd, sqrt.S, mean))
	print(bwplot(sqrt.S ~ wd, data = data, scales = list(x = list(rot = 90)),main="Variances Within Weeks and Diets"))
	
	# Plot coefficients with Poisson (from Bolker 2009)
	glm.lis <-lmList(S ~ week | individual, data = data, family = "poisson")
	try(print(plot(glm.lis, scale = list(x = list(relation = "free")),main="S, Family = Poisson")))
	try(print(qqmath(glm.lis)))
	
	grpMeans <- with(data, tapply(S, list(wd),mean))
	summary(grpMeans)

	grpVars <- with(data, tapply(S, list(wd),var))
	
	plot(grpVars ~ grpMeans, main="S, Var > Mean, Try quasiPoisson or negative binomial")
	abline(c(0, 1), lty = 2)

# Srar	
	df.list <- lmList(Srar ~ week | individual, subset = diet=="DF", data= data)
	hn.list <- lmList(Srar ~ week | individual, subset = diet=="HN", data= data)
	ln.list <- lmList(Srar ~ week | individual, subset = diet=="LN", data= data)

	df.coef <- coef(df.list)
	hn.coef <- coef(hn.list)
	ln.coef <- coef(ln.list)
	
	old <- par(mfrow=c(1,2))
	boxplot(df.coef[,1], hn.coef[,1], ln.coef[,1], main=paste("Intercepts of Srar"),names=c("DF", "HN","LN"))
	boxplot(df.coef[,2], hn.coef[,2], ln.coef[,2], main=paste("Slopes of Srar"),names=c("DF", "HN","LN"))
	par(old)
	rm(df.list,hn.list,ln.list,df.coef,hn.coef,ln.coef)
# H	
	df.list <- lmList(H ~ week | individual, subset = diet=="DF", data= data)
	hn.list <- lmList(H ~ week | individual, subset = diet=="HN", data= data)
	ln.list <- lmList(H ~ week | individual, subset = diet=="LN", data= data)

	df.coef <- coef(df.list)
	hn.coef <- coef(hn.list)
	ln.coef <- coef(ln.list)
	
	old <- par(mfrow=c(1,2))
	boxplot(df.coef[,1], hn.coef[,1], ln.coef[,1], main=paste("Intercepts of H"),names=c("DF", "HN","LN"))
	boxplot(df.coef[,2], hn.coef[,2], ln.coef[,2], main=paste("Slopes of H"),names=c("DF", "HN","LN"))
	par(old)
	rm(df.list,hn.list,ln.list,df.coef,hn.coef,ln.coef)
}
  
bootFun <-
	  function(x) {
	    th <- getME(x, "theta")
	    nvec <- sapply(getME(x, "cnms"), length)
	    scaleTh <- (isLMM(x) || isNLMM(x))
	    useSc <- as.logical(x@devcomp$dims["useSc"])
	    tn <- lme4:::tnames(x, old = FALSE, prefix = c("sd", "cor"))
	    if (scaleTh) {
	      ss <- setNames(lme4:::Cv_to_Sv(th, n = nvec, s = sigma(x)), 
	                     c(tn, "sigma"))
	    } else if (useSc) {
	      ss <- setNames(c(lme4:::Cv_to_Sv(th, n = nvec), sigma(x)), 
	                     c(tn, "sigma"))
	    } else {
	      ss <- setNames(lme4:::Cv_to_Sv(th, n = nvec), tn)
	    }
	    c(ss, fixef(x))
	  }
