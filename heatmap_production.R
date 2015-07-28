heatmap_production<-function(heatmapinfo3,num,groupsvector,patternsvector,col.label,col.map,filename,workingdirectory)
{
 
    num_groups<-unique(groupsvector)
    ends<-(length(num_groups)*num)
    starts<-(length(num_groups)*(num-1))+1
    new_patterns<-patternsvector[starts:ends]
    namesvec<-c(0,0,0)
    patternsvector1<-c()
    heatmapinfo<-heatmapinfo3[,2:ncol(heatmapinfo3)]
    heatmapinfo<-as.matrix(heatmapinfo)
    heatmapinfo<-unname(heatmapinfo)
    mode(heatmapinfo)<-'numeric'
    
    if(length(unique(new_patterns))>2)
    {
        pdf(paste(workingdirectory,"/",filename,"/","Pattern",num-1,".pdf",sep=""))
        heatmap.2(heatmapinfo,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,labRow=NA,ColSideColors=col.map,main=c(paste("pattern",num-1,sep="")," ",col.label))
        dev.off()
        return('Sorry, mean split is only available for a direct mean comparison between two groups')
    }
    
    for(i in 1:(length(unique(new_patterns))+1))
    {
        namesvec[i]<-paste("group",i,sep="")
    }
    
    if(length(new_patterns==2))
    {
        lengths<-length(which(groupsvector==1))
        patternsvector1<-c(patternsvector1,rep(0,lengths))
        lengths<-length(which(groupsvector==2))
        patternsvector1<-c(patternsvector1,rep(1,lengths))
    }
    else
    {
        k<-0
        for(i in 1:length(new_patterns))
        {   
            if(i==1)
            {
                lengths<-length(which(groupsvector==i))
                patternsvector1<-c(patternsvector1,rep(k,lengths))
                next
            }
            if(i==length(new_patterns) && patternsvector1)
            {
                lengths<-length(which(groupsvector==i))
                patternsvector1<-c(patternsvector1,rep(k,lengths))
                next
            }
            if(new_patterns[i]==new_patterns[i-1])
            {
                lengths<-length(which(groupsvector==i))
                patternsvector1<-c(patternsvector1,rep(k,lengths))
            }
            if(new_patterns[i]!=new_patterns[i-1])
            {
                lengths<-length(which(groupsvector==i))
                k<-k+1
                patternsvector1<-c(patternsvector1,rep(k,lengths))
            }
        }
    }
    
    k<-1
    for(i in 1:length(patternsvector1))
    {
        if(i == 1)
        {
            group_membership<-1
            next
        }
        if(patternsvector1[i]==patternsvector1[i-1])
        {
            group_membership<-cbind(group_membership,i)
        }
        if(patternsvector1[i]!=patternsvector1[i-1])
        {
            group_membership<-as.vector(group_membership)
            assign(namesvec[k],group_membership)
            k<-k+1
            group_membership<-c(i:length(patternsvector1))
            assign(namesvec[k],group_membership)
            break
        }
    }
    
        heatmapinfo1<-get(namesvec[1])
        heatmapinfo2<-get(namesvec[2])
        heatmapinfo_down<-c()
        heatmapinfo_up<-c()
        heatmapinfo_other<-c()
        names_down<-c()
        names_up<-c()
        names_other<-c()
        
        for(i in 1:nrow(heatmapinfo))
        {
            ave1<-heatmapinfo[i,heatmapinfo1]
            ave1<-as.vector(ave1)
            mode(ave1)<-'numeric'
            ave1<-ave(ave1)
            ave2<-heatmapinfo[i,heatmapinfo2]
            ave2<-as.vector(ave2)
            mode(ave2)<-'numeric'
            ave2<-ave(ave2)
            
            if(ave1[1]<0 && ave2[1]>0)
            {
                column_dow<-c(heatmapinfo[i,heatmapinfo1],heatmapinfo[i,heatmapinfo2])
                heatmapinfo_down<-rbind(heatmapinfo_down,column_dow)
                names_down<-c(names_down,as.character(heatmapinfo3[i,1]))
            }
            
            if(ave1[1]>0 && ave2[1]<0)
            {
                column_up<-c(heatmapinfo[i,heatmapinfo1],heatmapinfo[i,heatmapinfo2])
                heatmapinfo_up<-rbind(heatmapinfo_up,column_up)
                names_up<-c(names_up,as.character(heatmapinfo3[i,1]))
            }
            
            else
            {
                column_other<-c(heatmapinfo[i,heatmapinfo1],heatmapinfo[i,heatmapinfo2])
                heatmapinfo_other<-rbind(heatmapinfo_other,column_other)
                names_other<-c(names_other,as.character(heatmapinfo3[i,1]))
            }
        }
    heatmapinfo_down<-as.matrix(heatmapinfo_down)
    mode(heatmapinfo_down)<-'numeric'
    heatmapinfo_up<-as.matrix(heatmapinfo_up)
    mode(heatmapinfo_up)<-'numeric'
    
    pdf(paste(workingdirectory,"/",filename,"/","Pattern",num-1,"(1).pdf",sep=""))
    heatmap.2(heatmapinfo_down,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.5,labRow=names_down,ColSideColors=col.map,main=c(paste("pattern",num-1,"(1)",sep="")," ",col.label),labCol=1:ncol(heatmapinfo))
    dev.off()
    pdf(paste(workingdirectory,"/",filename,"/","Pattern",num-1,"(2).pdf",sep=""))
    heatmap.2(heatmapinfo_up,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.5,labRow=names_up,ColSideColors=col.map,main=c(paste("pattern",num-1,"(2)",sep="")," ",col.label),labCol=1:ncol(heatmapinfo))
    dev.off()
    pdf(paste(workingdirectory,"/",filename,"/","Pattern",num-1,"(3).pdf",sep=""))
    heatmap.2(heatmapinfo_other,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.5,labRow=names_other,ColSideColors=col.map,main=c(paste("pattern",num-1,"(3)",sep="")," ",col.label),labCol=1:ncol(heatmapinfo))
    dev.off()
    
    out_up<-cbind(names_up,heatmapinfo_up)
    out_down<-cbind(names_down,heatmapinfo_down)
    out_other<-cbind(names_other,heatmapinfo_other)
    namesvec2<-list()
    out_up<-unname(out_up)
    out_down<-unname(out_down)
    out_other<-unname(out_other)
    
    namesvec2[["up"]]<-out_up
    namesvec2[["down"]]<-out_down
    namesvec2[["other"]]<-out_other
    write.table(names_up, file = paste(workingdirectory,"/",filename,"/pattern",num-1,"(2).txt",sep=""), row.names = FALSE)
    write.table(names_down, file = paste(workingdirectory,"/",filename,"/pattern",num-1,"(1).txt",sep=""), row.names = FALSE)
    write.table(names_other, file = paste(workingdirectory,"/",filename,"/pattern",num-1,"(3).txt",sep=""), row.names = FALSE)
    
    return(namesvec2)
}