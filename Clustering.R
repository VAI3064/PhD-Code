function(data.input.file,fold.change.col,ID.col)
{
	all.data<-read.delim(data.input.file,header=T)
	all.data<-as.matrix(all.data)
	library(HSAUR)
	library(ggplot2)
	library(reshape2)
	library(fpc)
	library(foreign)
	library(psych)
	library(nFactors)
	library(gplots)
	library(lattice)
	
	#Create hierarchical clustering plot for each transcript#
	FCs<-all.data[,fold.change.col]
	FCs<-as.matrix(FCs)
	row.names(FCs)<-all.data[,ID.col]
	mode(FCs)<-'numeric'
	hclust.matrix<-hclust(dist(t((FCs)),method='euclidean'),method='complete')
	plot(hclust.matrix,axes=F)
	
	#Create hierarchical clustering plot for each array#
	d2 <- all.data[hclust.matrix$order,]
	hclust.matrix.2<-hclust(dist(t(d2)))
	plot(hclust.matrix.2,labels=headers)
	
	#Factor Analysis and determination of number of optimum groups#
	cor.all<-cor(FCs)
	ev <- eigen(cor.all)
	ap <- parallel(subject=nrow(all[,-c(1,2)]),var=ncol(all[,-c(1,2)]),rep=100,cent=.05)
	nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
	plotnScree(nS)
	groupsnumber<-nS$Components[[3]]
	factor.analysis<-fa(cor.all,nfactors=groupsnumber,rotate='varimax')
	factor.analysis
	plot(factor.analysis)
	fa.diagram(factor.analysis,simple=F,cex=0.7)

	#Display factor loadings#
	variable.group<-colnames(all[,-c(1,29,30,31,32)])
	melted <- cbind(variable.group, melt(factor.analysis$loadings[,1:groupsnumber]))
	ggplot(data=melted) +
  		geom_bar(aes(x=Var1, y=value, fill=variable.group), stat="identity") +
  		facet_wrap(~Var2)
  	
  	#Correlation heatmap of array#
  	subset.all.num<-data.input.file[,-c(1)]
	subset.all.num<-as.matrix(subset.all.num)
	mode(subset.all.num)<-'numeric'
	cor.matrix<-cor(subset.all.num)
	levelplot(cor.matrix)
}
