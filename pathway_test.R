pathway_test<-function(GO_terms,groups_length,col.label,col.map,num,filename,workingdirectory,up_reg,down_reg,other_reg)
{
    if(length(GO_terms)==0)
    {
        return()
    }
    
    pathway<-read.delim(paste(workingdirectory,"/Functional_Searches/in_pathways_WHOLEGENOME_numbers.txt",sep=""),header=F)
    pathway<-as.matrix(pathway)
    GO_terms<-unname(GO_terms)
    
    out<-c()
    
    for(i in 1:nrow(GO_terms))
    {
        is_na<-is.na(GO_terms[i,])
        is_na2<-which(is_na==TRUE)
        ends<-is_na2[1]-1
        starts<-(groups_length+3)
        if(is.na(starts) || is.na(ends))
        {
            vec<-GO_terms[i,(groups_length+2):ncol(GO_terms)]
        }
        else
        {
            vec<-GO_terms[i,starts:ends]
        }
        number<-length(vec)
        out1<-cbind(GO_terms[i,(groups_length+2)],number)
        out<-rbind(out,out1)
    }
    
    out1<-out[,2]
    mode(out1)<-'numeric'
    white_balls<-sum(out1)
    pathway1<-pathway[,2]
    mode(pathway1)<-'numeric'
    black_balls<-sum(pathway1)
    results<-c()
    
    for(i in 1:nrow(out))
    {
        start_numbers<-grep(out[i,1],pathway[,1])
        
        if(length(start_numbers)>1)
        {
            for(k in 1:length(start_numbers))
            {
                if(nchar(out[i,1])==nchar(pathway[start_numbers[k],1]))
                {
                    start_numbers<-start_numbers[k]
                    break
                }
            }
        }
        
        if(length(start_numbers)==0)
        {
            start_numbers<-which(pathway[, 1]==out[i, 1])
        }
        
        number_genome<-pathway[start_numbers,2]
        mode(number_genome)<-'numeric'
        number_sample<-out[i,2]
        mode(number_sample)<-'numeric'    
        result<-phyper(number_sample,number_genome,(black_balls-number_genome),white_balls,lower.tail=FALSE)
        results<-c(results,result)
    }
    
    results<-unname(results)
    out<-cbind(GO_terms[,1:(groups_length+2)],results)
    
    gtr2<-out[,ncol(out)]
    mode(gtr2)<-'numeric'
    sig_enrich<-which(gtr2<=0.05)
    heatmap_matrix<-c()
    row_names_matrix<-c()
    
    for(i in 1:length(sig_enrich))
    {
        heatmap_matrix1<-out[sig_enrich[i],1:groups_length]
        row_names_matrix1<-out[sig_enrich[i],(groups_length+2)]
        heatmap_matrix<-rbind(heatmap_matrix,heatmap_matrix1)
        row_names_matrix<-rbind(row_names_matrix,row_names_matrix1)
    }
    mode(heatmap_matrix)<-'numeric'
    heatmap_matrix = heatmap_matrix[rowSums(!is.na(heatmap_matrix))!=0, colSums(!is.na(heatmap_matrix))!=0]
    mode(heatmap_matrix)<-'numeric'
    if(length(heatmap_matrix)>2*groups_length)
    {
        pdf(paste(workingdirectory,"/",filename,"/","Pattern",num-1,"Significantly_enriched_pathways.pdf",sep=""))
        heatmap.2(heatmap_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.5,labRow=row_names_matrix,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched Pathways",sep="")," ",col.label))
        dev.off()
    }
    
    if(length(down_reg) == 1 && length(up_reg) == 1)
    {
        return("No Mean Split")
    }
    
    if(length(up_reg) == 1)
    {
        up_reg_matrix<-c()
        down_reg_matrix<-c()
        names_up<-c()
        names_down<-c()
        for(i in 1:nrow(down_reg))
        {
            if(length(grep(down_reg[i,1],GO_terms))>0)
            {
                down_reg_matrix<-rbind(down_reg_matrix,down_reg[i,2:ncol(up_reg)])
                names_down<-c(names_down,down_reg[i,1])
            }
            else
            {
                next
            }
        }
        mode(down_reg_matrix)<-'numeric'
        down_reg_matrix = down_reg_matrix[rowSums(!is.na(down_reg_matrix))!=0, colSums(!is.na(down_reg_matrix))!=0]
        mode(down_reg_matrix)<-'numeric'
        if(length(down_reg_matrix)>2*groups_length)
        {
            pdf(paste(workingdirectory,"/",filename,"/","Pattern",num-1,"Significantly_enriched_pathways_(1).pdf",sep=""))
            heatmap.2(down_reg_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.6,labRow=names_down,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched Pathways(2)",sep="")," ",col.label))
            dev.off()
        }
    }
    
    if(length(down_reg)==1)
    {
        up_reg_matrix<-c()
        down_reg_matrix<-c()
        names_up<-c()
        names_down<-c()
        for(i in 1:nrow(up_reg))
        {
            if(length(grep(up_reg[i,1],GO_terms))>0)
            {
                up_reg_matrix<-rbind(up_reg_matrix,up_reg[i,2:ncol(up_reg)])
                names_up<-c(names_up,up_reg[i,1])
            }
            else
            {
                next
            }
        }
        mode(up_reg_matrix)<-'numeric'
        up_reg_matrix = up_reg_matrix[rowSums(!is.na(up_reg_matrix))!=0, colSums(!is.na(up_reg_matrix))!=0]
        mode(up_reg_matrix)<-'numeric'
        if(length(up_reg_matrix)>2*groups_length)
        {
            pdf(paste(workingdirectory,"/",filename,"/","Pattern",num-1,"Significantly_enriched_pathways(1).pdf",sep=""))
            heatmap.2(up_reg_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.4,labRow=names_up,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched Pathways(1)",sep="")," ",col.label))
            dev.off()
        }
    }
    
    else
    {
        up_reg_matrix<-c()
        down_reg_matrix<-c()
        other_reg_matrix<-c()
        names_up<-c()
        names_down<-c()
        names_other<-c()
        
        for(i in 1:nrow(up_reg))
        {
            if(length(grep(up_reg[i,1],GO_terms))>0)
            {
                up_reg_matrix<-rbind(up_reg_matrix,up_reg[i,2:ncol(up_reg)])
                names_up<-c(names_up,as.character(up_reg[i,1]))
            }
            else
            {
                next
            }
        }
        
        for(i in 1:nrow(down_reg))
        {
            if(length(grep(down_reg[i,1],GO_terms))>0)
            {
                down_reg_matrix<-rbind(down_reg_matrix,down_reg[i,2:ncol(up_reg)])
                names_down<-c(names_down,as.character(down_reg[i,1]))
            }
            else
            {
                next
            }
        }
        
        for(i in 1:nrow(other_reg))
        {
            if(length(grep(other_reg[i,1],GO_terms))>0)
            {
                other_reg_matrix<-rbind(other_reg_matrix,other_reg[i,2:ncol(other_reg)])
                names_other<-c(names_other,as.character(other_reg[i,1]))
            }
            else
            {
                next
            }
        }
        
        mode(up_reg_matrix)<-'numeric'
        mode(down_reg_matrix)<-'numeric'
        mode(other_reg_matrix)<-'numeric'
        
        if(length(up_reg_matrix) == 0)
        {
            up_matrix<-matrix(0,1,1)
        }
        
        up_reg_matrix = up_reg_matrix[rowSums(!is.na(up_reg_matrix))!=0, colSums(!is.na(up_reg_matrix))!=0]
        mode(up_reg_matrix)<-'numeric'
        
        if(length(up_reg_matrix)>2*groups_length)
        {
            pdf(paste(workingdirectory,"/",filename,"/","Pattern",num-1,"Significantly_enriched_pathway_(1_).pdf",sep=""))
            heatmap.2(up_reg_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.4,labRow=names_up,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched Pathway (1)",sep="")," ",col.label))
            dev.off()
        }
        
        if(length(down_reg_matrix) == 0)
        {
            down_reg_matrix<-matrix(0,1,1)
        }
            down_reg_matrix = down_reg_matrix[rowSums(!is.na(down_reg_matrix))!=0, colSums(!is.na(down_reg_matrix))!=0]
            mode(down_reg_matrix)<-'numeric'
        
        if(length(down_reg_matrix)>2*groups_length)
        {
            pdf(paste(workingdirectory,"/",filename,"/","Pattern",num-1,"Significantly_enriched_pathway_(2).pdf",sep=""))
            heatmap.2(down_reg_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.4,labRow=names_down,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched Pathway (2)",sep="")," ",col.label))
            dev.off()
        }
        
        if(length(other_reg_matrix) == 0)
        {
            other_reg_matrix<-matrix(0,1,1)
        }
            other_reg_matrix = other_reg_matrix[rowSums(!is.na(other_reg_matrix))!=0, colSums(!is.na(other_reg_matrix))!=0]
            mode(other_reg_matrix)<-'numeric'
        
        if(length(other_reg_matrix)>2*groups_length)
        {
            pdf(paste(workingdirectory,"/",filename,"/","Pattern",num-1,"Significantly_enriched_pathway_(3).pdf",sep=""))
            heatmap.2(other_reg_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.4,labRow=names_other,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched Pathway (3)",sep="")," ",col.label))
            dev.off()
        }
    }
    
    out<-c(names_up,names_down,names_other)
    identity_gene<-read.delim(paste(workingdirectory,"/Functional_Searches/in_pathways_WHOLEGENOME.txt",sep=""))
    identity_gene<-as.matrix(identity_gene)
    identity_gene<-unname(identity_gene)
    out1<-c()
    for(i in 1:length(row_names_matrix))
    {
        position1<-grep(row_names_matrix[i],identity_gene[,1])
        for(j in 1:length(position1))
            {
                for(k in 1:length(out))
                {
                    is_it<-which(identity_gene[position1,]==out[k])
                    if(length(is_it)==0)
                    {
                        next
                    }
                    else
                    out1<-rbind(out1,cbind(identity_gene[position1,1],out[k]))
                }
            }
    }
    
    write.table(out1,file = paste(workingdirectory,"/",filename,"/Pathway_Test_Results.txt",sep=""), row.names = FALSE)
}