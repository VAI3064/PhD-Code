match_drosophila2<- function (orthologs, descriptions)
{
	orthologs<- as.matrix(orthologs)
	descriptions<- as.matrix(descriptions)
	orthologs_and_descriptions<- ""
	orthologs_and_descriptions2<-""
        for(i in 1:length(orthologs[,1]))
        {
            f<-as.numeric(which(orthologs[i,2]==descriptions[,1]))
	    descriptions3<-c()
	    descriptions2<-c()
	    for(k in 1:length(f))
	    {
		descriptions3<-c(descriptions3,as.character(descriptions[f[k],1]))
		descriptions2<-c(descriptions2,as.character(descriptions[f[k],2]))
	    }
	    descriptions3<-paste(descriptions3,collapse=",",sep="")
	    descriptions2<-paste(descriptions2,colapse=",",sep="")
	    orthologs_and_descriptions<- rbind(orthologs_and_descriptions,c(as.character(orthologs[i,1]), descriptions3, descriptions2))
        }
	
	duplicated_elements<-duplicated(orthologs_and_descriptions[,1])
	orthologs_and_descriptions<-orthologs_and_descriptions[-1,]
	for(i in 1:nrow(orthologs_and_descriptions))
	{
		
		if(duplicated_elements[i]==TRUE)
		{
			next
		}
		else
		{
			f<-as.numeric(which(orthologs_and_descriptions[,1]==orthologs_and_descriptions[i,1]))
			descriptions5<-c()
			descriptions4<-c()
			for(k in 1:length(f))
			{
				descriptions4<-c(descriptions4,as.character(orthologs_and_descriptions[f[k],2]))
				descriptions5<-c(descriptions5,as.character(orthologs_and_descriptions[f[k],3]))
			}
		}
		descriptions4<-paste(descriptions4,collapse=",",sep="")
		descriptions5<-paste(descriptions5,collapse=",",sep="")
		orthologs_and_descriptions2<-rbind(orthologs_and_descriptions2,c(as.character(orthologs_and_descriptions[i,1]),descriptions4,descriptions5,"ortholog"))
	}
    return(orthologs_and_descriptions2)
}
