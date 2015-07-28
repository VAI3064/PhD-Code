orthologs2<- function (orthologs_list, annotated_genome,j)
{
	
descriptions<- ""

	for (i in 1:length(orthologs_list[,2]))
	{
	    descriptions2<-""
            starts<-1
	    ends<-10
                if(nchar(orthologs_list[i,2])>10)
		{
			numberinstring<-ceiling(nchar(orthologs_list[i,2])/11)
			for(m in 1:numberinstring)
			{
				probe<-substr(orthologs_list[i,2],starts,ends)
				starts<-starts+11
				ends<-ends+11
				f<-grep(probe,annotated_genome[,1],value=FALSE)
				descriptions2<-c(descriptions2,as.character(annotated_genome[f,2]))
			}
                    descriptions2<-descriptions2[-1]
                    descriptions3<-paste(descriptions2,sep=",",collapse="/")
                    descriptions<- rbind(descriptions, c(as.character(orthologs_list[i,1]),as.character(orthologs_list[i,2]),as.character(orthologs_list[i,3]), descriptions3, orthologs_list[i,4:(max(j)+4)]))
		}
                
                if(nchar(orthologs_list[i,2])>2 && nchar(orthologs_list[i,2])<10 )
                {
                    descriptions<- rbind(descriptions, c(as.character(orthologs_list[i,1]),as.character(orthologs_list[i,2]),as.character(orthologs_list[i,3]), "NONE", orthologs_list[i,4:(max(j)+4)]))
                }
                
                if(nchar(orthologs_list[i,2])==10)
                {
                    probe1<-orthologs_list[i,2]
                    f<-grep(probe1,annotated_genome[,1],value=FALSE)
                    descriptions<- rbind(descriptions, c(as.character(orthologs_list[i,1]),as.character(orthologs_list[i,2]),as.character(orthologs_list[i,3]), as.character(annotated_genome[f,2]), orthologs_list[i,4:(max(j)+4)]))
                }
                
                    if(is.na(orthologs_list[i,2]))
                {
                    descriptions<- rbind(descriptions, c(as.character(orthologs_list[i,1]),as.character(orthologs_list[i,2]),as.character(orthologs_list[i,3]), "NONE", orthologs_list[i,4:(max(j)+4)]))
                }
	}

	return (descriptions)
}
	