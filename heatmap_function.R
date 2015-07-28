function(data.in,ID.search)
{
	data<-read.delim(data.in,header=T)
	data<-as.matrix(data)

	apriori<-read.delim(ID.search,header=F)
	apriori<-as.matrix(apriori)

	vector.out<-c()
	for(i in 1:nrow(apriori))
	{
  		positions<-grep(apriori[i,1],data[,1])
  		if(length(positions)>0)
  		{
   	 		data[positions[1],2]<-apriori[i,2]
  		}
	}

	names<-c()
	vector.out<-c()
	for(i in 1:nrow(apriori))
	{
  		positions<-grep(apriori[i,1],data[,1])
		if(length(grep(apriori[i,1],names))>=1)
		{
 			 next
		}
		else
		{
 		 	if(length(positions)==1)
  			{
    			vector.out<-rbind(vector.out,data[positions,2:ncol(data)])
    			names<-c(names,apriori[i,1])
  			}
  			if(length(positions)>1)
  			{
    			vector.out2<-c()
    			for(j in 2:ncol(data))
    				{
      					vector.out2<-c(vector.out2,mean(as.numeric(data[positions,j])))
   					 }
    			vector.out<-rbind(vector.out,vector.out2)
    			names<-c(names,apriori[i,1])
  			}
  			if(length(positions)==0)
  			{
    			next
    			names<-c(names,apriori[i,1])
  			}
		}
	}

	FCs<-vector.out
	mode(FCs)<-'numeric'

	library(gplots)
	colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
	heatmap.2(log2(FCs),col=colors,scale='none',key=T,keysize=1,density.info='none',trace='none',labCol=colnames(FCs),labRow=names,Rowv=F,Colv=F,symm=F,cexRow=0.5, cexCol=0.7)
}
