all_orthologs_search<-function(workingdirectory,genes_vec,filename,num)
{
    orthologs<-read.delim(paste(workingdirectory,"/Ortholog_Searches/All_Homologs.txt",sep=""))
    orthologs<-as.matrix(orthologs)
    orthologs<-unname(orthologs)
    new_orthologs<-c()
    
    for(i in 1:length(genes_vec))
    {
        new_orthologs1<-c()
        if(nchar(genes_vec[i]>10))
        {
            ss<-substr(genes_vec[i],1,10)
            f<-which(orthologs==ss,arr.ind=T)
            if(length(f)==0)
            {
                next
            }
            if(length(f)>2)
            {
                length_f<-seq.int(1,length(f),2)
                for(j in 1:length(length_f))
                {
                    new_orthologs1<-orthologs[f[length_f[j]],]
                    new_orthologs<-rbind(new_orthologs,new_orthologs1)
                }
            }
            else
            {
                new_orthologs1<-orthologs[f[1],]
                new_orthologs<-rbind(new_orthologs,new_orthologs1)
            }
        }
        else
        {
            next
        }
    }
    
   write.table(new_orthologs, file = paste(workingdirectory,"/",filename,"/","pattern_",num-1,"_Large_Ortholog_Table.txt",sep=""), row.names = FALSE)
}