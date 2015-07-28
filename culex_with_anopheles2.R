culex_with_anopheles2<- function(culex_anopheles, genes_of_interest_2,j)
{
    organism_orthologs_list<- ""
    for(i in 1:length(culex_anopheles[,1]))
    {
        f<-which(genes_of_interest_2[i]==culex_anopheles[,1])
        if(length(f)==0)
        {
            next
        }
        else
        {
                organism_orthologs_list<- rbind(organism_orthologs_list, c(culex_anopheles[f[1],1],culex_anopheles[f[1],2],culex_anopheles[f[1],3],"Culex2Anopheles", genes_of_interest_2[i,2:(max(j)+2)]))
        }
    }
    return (organism_orthologs_list)
}

