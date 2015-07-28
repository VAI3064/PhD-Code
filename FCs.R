FCs<- function (x,gaga_model,groups,names1,names2,workingdirectory,filename)
{
	log2FC<- rep(1,length(x[,1]))
	InverseLog<- rep(1,length(x[,1]))
	mpos<-posmeansGG(gaga_model,x=x,groups=groups,underpattern=1)
	Fold2<-""
	fc <- mpos[,1]-mpos[,2]
	for ( i in 1: length(x[,1]))
	{
		log2FC[i]<-fc[i]
		InverseLog[i]<-2^fc[i]
	}
	
	Fold<-cbind(names1,log2FC,InverseLog)
	for(i in 1:nrow(names2))
	{
		position<-grep(names2[i,1],Fold[,1])
		for(j in 1:length(position))
		{
			Fold2<-rbind(Fold2,Fold[position[j],])
		}
	}
	Fold2<-Fold2[-1,]
	write.table(Fold2,file = paste(workingdirectory,"/",filename,"/pattern1_fold_change.txt",sep=""), row.names = FALSE)
}
