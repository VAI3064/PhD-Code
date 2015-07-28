function(all.fold.changes,hub.list)
{
	#Input all fold changes for all arrays, starting with ID and one column for each fold change, read in list of hub IDs that are the same as column one#
	alldata<-read.delim(all.fold.changes,header=T)
	alldata<-as.matrix(alldata)

	colours.vector<-colors(distinct = TRUE)
	colours.vector<-colours.vector[-c(137:234)]
	colours.vector<-colours.vector[-21]

	hubs<-read.delim(hub.list,header=F)
	hubs<-as.matrix(hubs)
	values<-length(hubs)
	colours.for.graph<-colours.vector[sample(1:403,values-1, replace=F)]
	colours.for.graph<-c(colours.for.graph)
	
	#Plot a graph of individual fold changes for each array#
	plot(1:ncol(all.fold.changes)-1,1:ncol(all.fold.changes)-1,col=colours.for.graph[1],ylim=c(-3,6),xlim=c(0,ncol(all.fold.changes)),type='n',ylab='Gene Centered FC',xlab='Index')
	for(i in 1:length(hubs))
	{
  		positions<-grep(hubs[i],alldata[,1])
  		plot.points<-alldata[positions[1],2:ncol(alldata)]
  		plot.points<-unname(plot.points)
  		mode(plot.points)<-'numeric'
  		points(1:ncol(all.fold.changes)-1,log2(plot.points),col=colours.for.graph[i])
  		lines(1:ncol(all.fold.changes)-1,log2(plot.points),col=colours.for.graph[i])
  		par(new=T)
	}
	legend('topright',hubs,cex=.45,col=colours.for.graph,lty=1, bty='n')
}	
  