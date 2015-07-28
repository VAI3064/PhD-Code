heatmaps<- function (genes_vec, x2)
{
    heatmapinfo<-c()
    rows<-x2[,1]
    x3<-as.character(rows)
    for(i in 1:length(genes_vec))
    {
        f<-which(genes_vec[i]==rows)
        if(length(f)>1)
        {
            for(k in 1:length(f))
            {
               if(nchar(genes_vec[i])<10)
               {
                    name<-genes_vec[i]
                    heatmapinfo<-rbind(heatmapinfo,cbind(name,x2[f[k],12:ncol(x2)]))
               }
                if(nchar(genes_vec[i])>10)
                {
                    name<-substr(genes_vec[i],1,10)
                    heatmapinfo<-rbind(heatmapinfo,cbind(name,x2[f[k],12:ncol(x2)]))
                }
            }
        }
        
        else
        {
           if(nchar(genes_vec[i])<10)
           {
                name<-genes_vec[i]
                heatmapinfo<-rbind(heatmapinfo,cbind(name,x2[f,12:ncol(x2)]))
           }
           if(nchar(genes_vec[i])>10)
            {
                name<-substr(genes_vec[i],1,10)
               heatmapinfo<-rbind(heatmapinfo,cbind(name,x2[f,12:ncol(x2)]))
            }
        }
    }
    return(heatmapinfo)
}