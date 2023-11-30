# This code was written by Dr. Noelle Beckman

# MV Without Outliers: LN_3_6 & DF_5_4
# MV without LN weeks 5 - 8

library(data.table)
library(vegan)

location = " " # paste your location file path inside the quotes
file.name = c("KC_Filtered_QC_L5_20210325.txt")
data=read.table(file.path(location,file.name),header=TRUE, sep="\t")
redata=t(data[,-1])
colnames(redata)=data[,1]


# Without one row with outliers
redata2=redata[,-c(108, 109, 114:117, 36)]
redata2=t(redata2)
redata2=as.data.frame(redata2)
		# Abundant species can dominate. Common	transformations is the sqrt or the 4th root, that gives increasingly more weight to rare species
		# to give equal waiting for common and abundant species, can analyze the presence/absence data.
		# We do not want to standardize
			range(redata2)
			
			redata2_sqrt<-sqrt(redata2)
			
		# This makes the communities more similar to each other.	
			redata2_fourth<-sqrt(redata2_sqrt)
		
		# MDS (Non-metric multidimensional Scaling)
		# the distance metric defaults to Bray and common ecological data transformations are turned on by default. We want to turn these off.		

		# From Clark and Waricke
		# Stress<0.05 gives an excellent representation
		# Stress <0.1 corresponds to good ordination with no prospect of misleading interpretation; higher dimensions will not add additional infor to overall structure
		# Stress <0.2 gives potentially useful 2-dimensional pcitures, though for values at the upper end of this range too much relance shound not be placed on the detail of the plot
		# Stress >0.3 indicates that points are close to being arbitrarily places in the 2 D ordination space.

		NMDS.scree<-function(x,dis="bray",bin=FALSE) { #where x is the name of the data frame variable
			plot(rep(1,10),replicate(10,metaMDS(x,distance=dis,binary=bin,autotransform=F,k=1,engine="monoMDS")$stress),xlim=c(1,11),ylim=c(0,0.3),xlab="# of Dimensions",ylab="Stress",main="NMDS stress plot")
			for (i in 1:10) { #(nrow(x)-2)
				points(rep(i+1,10),replicate(10,metaMDS(x,distance=dis,binary=bin, autotransform=F,k=i+1)$stress))
			}
			abline(h=0.2, lty=4)
			abline(h=0.1, lty=2)
			abline(h=0.05,lty=3)	
		}	
			
		# For Bray-Curtis
			NMDS.scree(redata2,dis="bray",bin=T) # k=3 (stress= ca. 0.12) or k= 4 (stress= ca. 0.09), k=2 (stress = 0.16), k=7 (Stress=0.05)
			NMDS.scree(redata2_sqrt) # k=2 (stress = ca. 0.16), k=3 (stress= ca. 0.12), k=4 (stress = ca. 0.09), k=7 (stress=0.05) 
			NMDS.scree(redata2_fourth) # k=2 (stress = ca. 0.17), k=3 (stress= ca. 0.13), k=4 (stress = ca. 0.09), k=8 (stress = 0.05)
	
			mydata.mds0 <- metaMDS(redata2, distance = "bray", trymax=100,autotransform=FALSE,k=2)
			
			mydata.mds1 <- metaMDS(redata2_sqrt,distance="bray",trymax=100,autotransform=FALSE,k=2)
			names(mydata.mds1)
			mydata.mds1 # stress is 0.169
			mydata.mds = mydata.mds1

    # Sorenson's (presence/absence Bray-Curtis)
      mydata.mds2 <- metaMDS(redata2,distance="bray",trymax=100,autotransform=FALSE,k=2,binary=TRUE)
    # stress is nearly zero: you may have insufficient data

     meta.md<-metaMDS
     
     variableScores0 <- mydata.mds0$species
     sampleScores0 <- mydata.mds0$points
     
		variableScores1 <- mydata.mds1$species
   	sampleScores1 <- mydata.mds1$points
   	
   	variableScores2 <- mydata.mds2$species
   	sampleScores2 <- mydata.mds2$points
   	
   		# Open black circles correspond to samples and red crosses indicate taxa.  
   	par(mfrow = c(1, 3))		
   	plot(mydata.mds0,display="species", main = "Raw redata")
   	plot(mydata.mds1,display="species", main = "Redata_sqrt")
   	plot(mydata.mds2,display="species", main = "Redata Sorenson")
   			
   			#text(mydata.mds,display=c("site"))
  	
  			tmp=data.table(ID=rownames(redata2))
			tmp[,diet:=substring(ID,1,2)]
			tmp[,ind:=substring(ID,4,4)]
			tmp[,week:=substring(ID,6,6)]
			
			
			tmp[,ind:=as.numeric(ind)]
			tmp[diet=="DF",individual:=ind]
			tmp[diet=="HN",individual:=ind+7]
			tmp[diet=="LN",individual:=ind+13]
			
			tmp[,dietweek:=paste0(diet,week)]	

# We can draw convex hulls connecting the vertices of the points made by
# these communities on the plot			
			
			ordiplot(mydata.mds)
			ordihull(mydata.mds,choices=c(1,2),groups=tmp$diet,draw="polygon",col="grey90",label=F)
			points(mydata.mds,display="species",col="red")
			orditorp(mydata.mds,display="sites",col=c(rep("cyan4",5),rep("darkorange",5)),
			   air=0.01,cex=1)
			
# KC --> I this to look at differences between all 3
par(mfrow = c(1, 3))
ordiplot(mydata.mds0, main = "Raw Redata")
ordiplot(mydata.mds1, main = "Redata_Sqrt")
ordiplot(mydata.mds2, main = "Sorenson")


# One can also plot ellipses and "spider graphs" using the functions 
# `ordiellipse` and `orderspider` which emphasize the centroid of the 
# communities in each treatment

			 ordiplot(mydata.mds0)  
			 ordiellipse(mydata.mds,groups=tmp$diet)
			
# ordispider
			 ordiplot(mydata.mds0, type = "none", main = "Redata Raw")
			 ordispider(mydata.mds0,groups=tmp$diet)
			 
			 ordiplot(mydata.mds1)
			 ordispider(mydata.mds1,groups=tmp$diet, main = "Redata_Sqrt")
			 
			 ordiplot(mydata.mds2)
			 ordispider(mydata.mds2,groups=tmp$diet, main = "Sorenson")
			plot.new() 
# Another alternative is to plot a minimum spanning tree (from the 
# function `hclust`), which clusters communities based on their original 
# dissimilarities and projects the dendrogram onto the 2-D plot
ordiplot(mydata.mds)
orditorp(mydata.mds,display="species",col="red")
orditorp(mydata.mds,display="sites",col=c(rep("cyan4",5),rep("darkorange",5)),
         air=0.01,cex=1)
ordicluster(mydata.mds,hclust(vegdist(redata,"bray"))) 
# Note that clustering is based on Bray-Curtis distances
# This is one method suggested to check the 2-D plot for accuracy

# You could also plot the convex hulls, ellipses, spider plots, etc. colored based on the treatments
# First, create a vector of color values corresponding of the same length as the vector of treatment values
colors=c(rep("red",55),rep("cyan4",41),rep("darkorange",25))
ordiplot(mydata.mds,type = "none", xlim=c(-1,2),choices=c(1,2))
#Plot convex hulls with colors based on treatment
	for(i in unique(tmp$diet)) {
  		ordihull(mydata.mds$point[grep(i,tmp$diet),],choices=c(1,2),draw="polygon",groups=tmp$diet[tmp$diet==i],col=colors[grep(i,tmp$diet)],label=F) 
  	} 
	#points(mydata.mds,display="species",col="black",pch=19,cex=0.5)
	#points(data.frame(mydata.mds$points)$MDS1,data.frame(mydata.mds$points)$MDS2,display="sites",col=colors,pch=tmp$week,cex=1.5)
	#orditorp(mydata.mds,display="sites",col=colors,air=0.01,cex=1,choices=c(2,3))
    legend(x=1.3,y=1.3,legend=c("DF","HN","LN"),col=c("black","cyan4","darkorange"),pch=19,cex=1.5)
         #stressplot(mydata.mds)
 
 
 ################ ORDISPIDER Plots ##############
 #  ordispider draws a ‘spider’ diagram where each point is connected to the group centroid
#     DF (black), HN(cyan4), LN(dotted)
    
colors = tmp$diet
colors[colors=="DF"]="black"
colors[colors=="HN"]="cyan4"
colors[colors=="LN"]="darkorange"

fill = tmp$diet
fill[fill=="DF"]="black"
fill[fill=="HN"]="cyan4"
fill[fill =="LN"]="darkorange"

shape = tmp$diet
shape[shape=="DF"]=22
shape[shape=="HN"]=21
shape[shape=="LN"]=23

lt = tmp$diet
lt[lt=="DF"]="solid"
lt[lt=="HN"]="solid"
lt[lt=="LN"]="solid"

par(mfrow = c(1, 2))


# Raw Redata
plot(mydata.mds0, type = "none", main = "Bray Curtis Raw Abundance", xlim=c(-1.0, 1.5),ylim=c(-1.0, 1.0))
  
#Plot convex hulls with colors based on treatment
	for(i in unique(tmp$diet)) {
		ordispider(mydata.mds0,groups=tmp$diet,show.group=tmp$diet[tmp$diet==i],
		           col=colors[grep(i,tmp$diet)],
		           lty=lt[grep(i, tmp$diet)], label=F)  
  	} 

    #legend(x=1.5, y=0.5, legend=unique(tmp$diet),col= c("black","cyan4","darkorange"),
           #lty = c("solid", "solid", "solid"), pch=10,cex=1)

    mds0_points <- mydata.mds0[["points"]]
    mean(abs(mds0_points[1:54, ])) # 0.243468
    mean(abs(mds0_points[55:96, ])) # 0.3302821
    mean(abs(mds0_points[97:116, ])) # 0.4736498

   # distance between DF and HN or LN
    (mean(mds0_points[55:96, ])) - (mean(mds0_points[1:54, ])) # -0.1484951
    (mean(abs(mds0_points[97:116, ]))) - (mean(abs(mds0_points[1:54, ]))) # 0.2301818
    
# Square Root Redata
    plot(mydata.mds1, type = "none", main = "Bray Curtis Square Root", xlim=c(-0.5, 1.0),ylim=c(-0.5, 0.5))
    
    #Plot convex hulls with colors based on treatment
    for(i in unique(tmp$diet)) {
      ordispider(mydata.mds1,groups=tmp$diet,show.group=tmp$diet[tmp$diet==i],
                 col=colors[grep(i,tmp$diet)],
                 lty=lt[grep(i, tmp$diet)], label=F)  
    } 
    
    #legend(x=1.5, y=0, legend=unique(tmp$diet),col= c("black","cyan4","darkorange"),
          # lty = c("solid", "solid", "solid"), pch=10,cex=1)
    
    mds1_points <- mydata.mds1[["points"]]
    mean(abs(mds1_points[1:54, ])) # 0.1307379
    mean(abs(mds1_points[55:96, ])) # 0.1798059
    mean(abs(mds1_points[97:116, ])) # 0.2675978
    # distance between DF and HN or LN
    (mean(abs(mds1_points[55:96, ]))) - (mean(abs(mds1_points[1:54, ]))) # 0.04906805
    (mean(abs(mds1_points[97:116, ]))) - (mean(abs(mds1_points[1:54, ]))) # = 0.1368599
    
# Sorenson
    plot(mydata.mds2, type= "none", main = "Sorenson", xlim=c(-0.1, .1),ylim=c(-0.15, .125))
    
    #Plot convex hulls with colors based on treatment
    for(i in unique(tmp$diet)) {
      ordispider(mydata.mds2,groups=tmp$diet,show.group=tmp$diet[tmp$diet==i],
                 col=colors[grep(i,tmp$diet)],
                 lty=lt[grep(i, tmp$diet)], label=F)  
    } 
    
    plot(mydata.mds2, type= "none", main = "Legend", xlim=c(0, 1),ylim=c(0, 1),
         col.axis = "white", lty = 0, tcl = NA)
  
    legend(x= -0.05, y=-0.2, legend=unique(tmp$diet),col= c("black","cyan4","darkorange"),
           lty = c("solid", "solid", "solid"), pch=10,cex=0.5, border = "white")    

    mds2_points <- mydata.mds2[["points"]]
    mean(abs(mds2_points[1:54, ])) # 0.01596883
    mean(abs(mds2_points[55:96, ])) # 0.0273584
    mean(abs(mds2_points[97:116, ])) # 0.02251755
    # distance between DF and HN or LN
    (mean(mds2_points[55:96, ])) - (mean(mds2_points[1:54, ])) # 0.00130715
    (mean(mds2_points[97:116, ])) - (mean(mds2_points[1:54, ])) # -0.009976607
    
    
#################
## Permanova ##
#################
# From Clark and Waricke (pg 65-66 in the Primer manual):
# "It must be admitted that the correlation structure among samples through time, if any, 
# will be effectively ignored under permutation. Thus, differences in correlation structure
# through time among treatments (ie.e lack of spherecity) may, therefore, produce a statistically
# significant result. However, we consider the differences in correlation structure through time are
# indicative of (at least one type of) a treatment effect, so should warrant closer inspection...
# If the statistically significant results are obtained in a repeated measures analysis, the user
# may wish to accompnay the PERMANOVA with a separate test for spherecity (in the case of univariate data)
# or its analogue (using dissimilarities rather than differences) for multivariate responses, in order 
# to shed further light on the meaning of the results and appropriate inferences"
library(vegan)

# Adonis analysis suggests treatment effects are somewhat consistent throughout time for the square and fourth root but not for untransformed data.
# I should also do pairwise differences to see which treatments differ. Treatments differences could be due because differences in the means/centroids, 
# dispersion, or correlation structure through time among treatments (spherecity assumption).

# Diet and interaction between diet and week significant
    # I checked these results again and they match in the paper. There are small differences (i.e. p value in paper 0.0142 v 0.0156 when reran)
bc.adonis <- adonis2(redata2~tmp$diet*tmp$week,permutations=9999) # this is the results that we used in the paper
bc.adonis2 <- adonis2(redata2~tmp$week*tmp$diet,permutations=9999)

# Diet significiant and interaction between diet and week marginally significant
bc.adonis.sqrt <- adonis2(redata2_sqrt~tmp$diet*tmp$week,permutations=9999)
bc.adonis.sqrt2 <- adonis2(redata2_sqrt~tmp$week*tmp$diet,permutations=9999)

# Diet significiant and interaction between diet and week marginally significant
bc.adonis.fourthroot <- adonis2(redata2_fourth ~tmp$diet*tmp$week,permutations=9999)
bc.adonis.fourthroot2 <- adonis2(redata2_fourth ~tmp$week*tmp$diet,permutations=9999)

# Diet and interaction between diet and week significant
sor.adonis <- adonis2(redata2~tmp$diet*tmp$week,permutations=9999,binary=TRUE) # this is the result used in the paper
sor.adonis2 <- adonis2(redata2~tmp$week*tmp$diet,permutations=9999,binary=TRUE)

# To determine if spherecity assumption is met
bc.dist <- vegdist(redata2, method="bray")
bc.dist.sqrt<-vegdist(redata2_sqrt,method="bray")
bc.dist.fourth<-vegdist(redata2_fourth,method="bray")
sorenson<-vegdist(redata2,method="bray",binary=TRUE)


# Levene's test for spherecity suggests that homogeneity of variances among time points and that treatment effects are not caused by differences in the 
	# correlation structure through time among individuals.
k=1
g=1
group_length=653
	data.tmp<-bc.dist # Bray Curtis
	bc.matrix<-as.matrix(data.tmp)
	data.bc<-data.frame(bc=vector("numeric", group_length),time1=vector("integer", group_length),time2=vector("integer", group_length),group=vector("integer", group_length))
	for(i_time in 1:7){
		for(j_time in 2:8){
			test<-bc.matrix[substring(rownames(bc.matrix),6,6)==i_time,substring(colnames(bc.matrix),6,6)==j_time]
			for(j in 1:dim(test)[2]){
				for(i in 1:dim(test)[1]){
					if(substring(rownames(test),1,4)[i]==substring(colnames(test),1,4)[j]){
						data.bc$bc[k]<-test[i,j]
						data.bc$time1[k]<-i_time
						data.bc$time2[k]<-j_time
						data.bc$group[k]<-g
						k=k+1	
			}
		}
	}
			g=g+1
		}
	}
	
	rm(data.tmp)	

	library(car)
	
	# Spherecity assumptions holds up for all transformations
	leveneTest(bc ~ as.factor(group), data=data.bc[-which(data.bc$bc==0),]) 
	  # bc == 0 removes all values that are 0
	
k=1
g=1
group_length=653
	data.tmp<-sorenson #Sorenson
	bc.matrix<-as.matrix(data.tmp)
	data.bc<-data.frame(bc=vector("numeric", group_length),time1=vector("integer", group_length),time2=vector("integer", group_length),group=vector("integer", group_length))
	for(i_time in 1:7){
	  for(j_time in 2:8){
	    test<-bc.matrix[substring(rownames(bc.matrix),6,6)==i_time,substring(colnames(bc.matrix),6,6)==j_time]
	    for(j in 1:dim(test)[2]){
	      for(i in 1:dim(test)[1]){
	        if(substring(rownames(test),1,4)[i]==substring(colnames(test),1,4)[j]){
	          data.bc$bc[k]<-test[i,j]
	          data.bc$time1[k]<-i_time
	          data.bc$time2[k]<-j_time
	          data.bc$group[k]<-g
	          k=k+1	
	        }
	      }
	    }
	    g=g+1
	  }
	}
	
	rm(data.tmp)	
	
	library(car)
	
	# Spherecity assumptions holds up for all transformations
	leveneTest(bc ~ as.factor(group), data=data.bc[-which(data.bc$bc==0),])
	
	
	# Should do a dispersion test for each level of diet and week, but the PERMANOVA+Primer book says to view results with caution
	#if n<=10 in each level and avoided for n<=5, which we have here (pg. 96).
	#If there is an interaction, then dispersion for main effects should be interpreted cautiously 
	# as this could be do to the interaction and not dispersion.
	# For each transformation, dispersion varies with diet but not with week.
	anova(betadisper(bc.dist,tmp$week))
	permutest(betadisper(bc.dist,tmp$week),permutations=how(nperm=9999))
	TukeyHSD(betadisper(bc.dist,tmp$week))
	plot(TukeyHSD(betadisper(bc.dist,tmp$week)))
	plot(betadisper(bc.dist,tmp$week))
	
	untest<-betadisper(bc.dist,tmp$diet)
	permutest(untest,permutations=how(nperm=9999))
	TukeyHSD(untest)
	plot(TukeyHSD(untest))
	plot(untest)

	permutest(betadisper(bc.dist.sqrt,tmp$week),permutations=how(nperm=9999))
	untest2<-betadisper(bc.dist.sqrt,tmp$diet)
	permutest(untest2,permutations=how(nperm=9999))
	TukeyHSD(untest2)
	plot(TukeyHSD(untest2))
	plot(untest2)
	
	anova(betadisper(bc.dist.sqrt,tmp$week))
	anova(betadisper(bc.dist.sqrt,tmp$diet))


permutest(betadisper(sorenson,tmp$week),permutations=how(nperm=9999))
untest2<-betadisper(sorenson,tmp$diet)
permutest(untest2,permutations=how(nperm=9999))
TukeyHSD(untest2)
plot(TukeyHSD(untest2))
plot(untest2)
	
	anova(betadisper(bc.dist.fourth,tmp$week))
	anova(betadisper(bc.dist.fourth,tmp$diet))


############# Did not use in the paper, but can be used
# Similarity percentage, simper (Clarke 1993) is based on the decomposition of Bray-Curtis dissimilarity index (see vegdist, designdist).
# http://www.uwyo.edu/mdillon/hor/2013/richards.pd f


# Average dissimilarity between all pairs of inter-group samples 

contrib <- simper(redata2_sqrt,tmp$diet)
lapply(contrib, FUN=function(x){x$overall})
contrib

# Average broken down by contributions from each species
data_contrib <- summary(contrib)

write.table(x = data_contrib$DF_HN, file = "~/Frass_Paper/Sabree_L5/BC/Dissim_DF_HN.txt",sep = "\t")
write.table(x = data_contrib$HN_LN, file = "~/Frass_Paper/Sabree_L5/BC/Dissim_HN_LN.txt",sep = "\t")
write.table(x = data_contrib$DF_LN, file = "~/Frass_Paper/Sabree_L5/BC/Dissim_DF_LN.txt",sep = "\t")


  
location2 = "~/Frass_Paper/Sabree_L5/BC"
file.name2 = c("Dissim_DF_HN.txt")
Dissim_DF_HN=read.table(file.path(location2,file.name2),header=TRUE, sep="\t")


contrib2 <- simper(redata2,tmp$diet)
lapply(contrib2, FUN=function(x){x$overall})
contrib2

# Average broken down by contributions from each species
data_contrib2 <- summary(contrib2)


# Also look into this distance decay: 
	# http://onlinelibrary.wiley.com/store/10.1046/j.1365-2699.1999.00305.x/asset/j.1365-2699.1999.00305.x.pdf;jsessionid=26355632F369FE137257C20D1E20AF7F.f04t04?v=1&t=idhtdv91&s=8243f7f1c66de7f9e3268f7500fa7a193ced07c5
	#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3413688/#pone.0041938-Nekola1
