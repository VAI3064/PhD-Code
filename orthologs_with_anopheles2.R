orthologs_with_anopheles2<- function (orthologs_list, genes_of_interest,j)
{	
    organism_orthologs_list<- ""
    vecnames<-c()
    for(i in 1:length(genes_of_interest[,1]))
    {
        f<-which(genes_of_interest[i]==orthologs_list[,2])
        vecnames1<-rep(NA,length(f))
        if(length(f)>1)
        {
            for(i in 1:length(f))
            {
                vecnames1[i]<-as.character(orthologs_list[f[i],3])
            }
            entry0<-as.character(orthologs_list[f[1],2])
            entry1<-paste(vecnames1,collapse=",")
            entry2<-as.character(orthologs_list[f[1],4])
            organism_orthologs_list<- rbind(organism_orthologs_list, c(entry0,entry1,entry2, genes_of_interest[i,2:(max(j)+2)]))
        }
        
        if(identical(f,integer(0)))
        {
            entry0<-as.character(genes_of_interest[i])
            entry1<-c("None")
            entry3<-c("None")
            organism_orthologs_list<- rbind(organism_orthologs_list, c(entry0,entry1,entry3, genes_of_interest[i,2:(max(j)+2)]))
        }
        
        if(length(f)==1)
        {
        organism_orthologs_list<- rbind(organism_orthologs_list, c(as.character(orthologs_list[f,2]),as.character(orthologs_list[f,3]),as.character(orthologs_list[f,4]), genes_of_interest[i,2:(max(j)+2)]))
        }
    }

    return(organism_orthologs_list)
}
