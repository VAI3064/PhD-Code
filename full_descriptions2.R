full_descriptions2<- function (genes_of_interest, new_GOI, anopheles_descriptions, aedes_description, culex_description, drosophila_description,j)
{
    all_descriptions<- matrix(" ", nrow = length(genes_of_interest[,1]), ncol = (max(j)+12))
    all_descriptions[,1]<- genes_of_interest[,1]
    all_descriptions[,12:(max(j)+12)]<- genes_of_interest[,2:(max(j)+2)]
    vector1<-rep("Empty",nrow(drosophila_description))
    drosophila_description<-cbind(drosophila_description,vector1)
        for(i in 1:length(all_descriptions[,1]))
        {
            f<-which(all_descriptions[i, 1]==anopheles_descriptions[, 1])
            if(length(f)==0)
            {
                all_descriptions[i,2]<- "NONE"  
            }
            else
            {
                all_descriptions[i,2]<- anopheles_descriptions[f,2]
            }
            
            g<-which(all_descriptions[i, 1]==aedes_description[, 1])
            
            if(length(g)==0)
            {
                nonevec<-c("NONE","NONE","NONE")
                all_descriptions[i,3:5]<- nonevec
            }
            else
            {
                all_descriptions[i,3:5]<- aedes_description[g,2:4]
            }
            
            h<-which(all_descriptions[i, 1]==culex_description[, 1])
            
            if(length(h)==0)
            {
                nonevec<-c("NONE","NONE","NONE")
                all_descriptions[i,6:8]<- nonevec
            }
            else
            {
              all_descriptions[i,6:8]<- culex_description[h,2:4]  
            }
            
            k<-which(new_GOI[i,1]==drosophila_description[,1])
            
            if(length(k)==0)
            {
                nonevec<-c("NONE","NONE","NONE")
                all_descriptions[i,9:11]<-nonevec
            }
            else
            {
               if(length(drosophila_description[k,2:4])==3)
               {
                all_descriptions[i,9:11]<- drosophila_description[k,2:4]
               }
               else
               {
                    nonevec<-c("NONE","NONE","NONE")
                    all_descriptions[i,9:11]<-nonevec
               }
            }
            
        }
    return(all_descriptions)
}
        