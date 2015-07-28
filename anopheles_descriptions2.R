anopheles_descriptions2<- function (anopheles_description_file, genes_of_interest_2,j)
{
    anopheles_description<- ""
    for(i in 1:length(genes_of_interest_2[,1]))
    {
        f<-grep(genes_of_interest_2[i],anopheles_description_file[,2],value=FALSE)
       if(length(f)>1)
       {
        anopheles_description<- rbind(anopheles_description, c(as.character(genes_of_interest_2[i]), as.character(anopheles_description_file[f[1],5]), genes_of_interest_2[i,2:(ncol(genes_of_interest_2)-1)]))
       }
       
       if(length(f)==1)
       {
        anopheles_description<- rbind(anopheles_description, c(as.character(genes_of_interest_2[i]), as.character(anopheles_description_file[f,5]), genes_of_interest_2[i,(ncol(genes_of_interest_2)-1)]))
       }
       
       if(identical(f,integer(0)))
       {
        anopheles_description<- rbind(anopheles_description, c(as.character(genes_of_interest_2[i]), "NA", genes_of_interest_2[i,2:(ncol(genes_of_interest_2)-1)]))
       }
    }
    return(anopheles_description)
}
