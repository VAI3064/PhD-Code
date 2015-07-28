pathways2 <- function(genes_vec,workingdirectory,all_orthologs)
{
    pathway<-read.delim(paste(workingdirectory,"/Functional_Searches/pathway_orthologs.txt",sep=""))
    pathway<-pathway[-1,]
    pathway<-as.matrix(pathway)
    pathway2<-c()
    out<-c()
    
    for(i in 1:length(genes_vec))
    {
        string<-genes_vec[i]
        if(nchar(string)>=10)
        {
            ss<-substr(string,1,10)
        }
        else
        {
            next
        }
        
        position<-grep(ss,pathway[,1])
        position2<-grep(ss,all_orthologs[,1])
        
        if(length(position)==0)
        {
            next
        }
        else
        {
            for(j in 1:length(position))
            {
                pathway0<-rbind(all_orthologs[position2,12:ncol(all_orthologs)])
                pathway1<-cbind(pathway[position[j],1],pathway[position[j],2],pathway0)
                pathway2<-rbind(pathway2,pathway1)
            }
        }
    }
    if(length(pathway2)==0)
    {
        print("No Enriched Pathways")
        return()
    }
    unique_pathways<-unique(pathway2[,2])
    unique_pathways<-as.vector(unique_pathways)
    unique_pathways<-unname(unique_pathways)
    mode(unique_pathways)<-'character'
    genes_in_pathway<-matrix(NA,ncol=10000,nrow=length(unique_pathways))
    
    for(i in 1:length(unique_pathways))
    {
        pp<-grep(unique_pathways[i],pathway2[,2])
        numbers<-c()
        genes_in_pathway1<-c()
        for(k in 1:length(pp))
        {
            numbers1<-pathway2[pp[k],3:ncol(pathway2)]
            numbers1<-unname(numbers1)
            numbers1<-as.matrix(numbers1)
            mode(numbers1)<-'numeric'
            numbers<-rbind(numbers,numbers1)
            pathwayxyz<-as.character(pathway2[pp[k],1])
            genes_in_pathway1<-c(genes_in_pathway1,pathwayxyz)
        }
        
        genes_in_pathway[i,1:length(genes_in_pathway1)]<-genes_in_pathway1
        
        mode(numbers)<-'numeric'
        aver<-c()
        aver1<-c()
        
        if(ncol(numbers)>1)
        {
            
            for(k in 1:ncol(numbers))
            {
                aver<-ave(numbers[,k])
                aver1<-cbind(aver1,aver[1])
            }   
        }
        else
        {
            aver1<-numbers
        }
        
        out<-rbind(out,aver1)
    }
    
    position_is_na<-c()
    for(b in 1:nrow(genes_in_pathway))
    {
        is_na<-is.na(genes_in_pathway[b,])
        which_na<-which(is_na==TRUE)
        position_is_na<-c(position_is_na,which_na[1])
    }
    
    cut_off<-max(position_is_na)
    genes_in_pathway<-genes_in_pathway[,-c((cut_off:ncol(genes_in_pathway)))]
    
    out<-as.matrix(out)
    mode(out)<-'numeric'
    genes_in_pathway<-cbind(out,unique_pathways,genes_in_pathway)
    
    
    return(genes_in_pathway)
}