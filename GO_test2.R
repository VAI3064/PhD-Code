GO_test2<-function(pattern,alls,workingdirectory,num,all_orthologs,groups_length,col.map,col.label,up_reg,down_reg,other_reg,filename)
{
    pattern_p<-rep(0,length(pattern))
    names(pattern_p)<-pattern
    all1<-alls[-c(1:3)]
    
    for(i in 1:length(pattern))
    {
        index<-which(pattern[i]==all1)
        if(identical(index,integer(0)))
        {
            next
        }
        else
        {
            all1<-all1[-index]
        }
    }
    
    all_p<-rep(1,length(all1))
    names(all_p)<-all1
    
    overall<-c(pattern_p,all_p)
    names1<-names(overall)
    
    for(i in 1:length(names1))
    {
        ss<-substr(names1[i],1,10)
        names1[i]<-ss
    }
    
    names(overall)<-names1
    
    ont1<-c("CC","MF","BP")
    ont<-c("Cellular Component","Molecular Function","Biological Process")
    ont2<-c("/GO_CC.txt","/GO_MF.txt","/GO_BP.txt")
    
       
    logvals<-regexpr("AGAP",names(overall))!=-1
    logvals2<-c()
    overall2<-c()
    
    for(i in 1:length(logvals))
    {
        if(logvals[i]==TRUE)
        {
            overall2<-c(overall2,overall[i])
        }
    }
    
    overall<-overall2
    
    topDiffGenes <- function(allScore){
    return(allScore < 0.01)}
    
    
    genes<-names(overall)
    
    for(k in 1:3)
    {
        sampleGOdata<-new("topGOdata",description=paste("Pattern",num-1,sep=""),ontology=ont1[k],allGenes=overall,geneSel=topDiffGenes,annot=annFUN.file,file=paste(workingdirectory,"/Functional_Searches",ont2[k],sep=""),nodeSize=5)
        
        if(length(usedGO(sampleGOdata))>5)
        {
            classicFisher<-runTest(sampleGOdata,algorithm="classic",statistic="fisher")
            elimFisher<-runTest(sampleGOdata,algorithm="elim",statistic="fisher")
            weightFisher<-runTest(sampleGOdata,algorithm="weight",statistic="fisher")
            allRes <- GenTable(sampleGOdata, classicFisher = classicFisher,elimFisher=elimFisher,weightFisher = weightFisher, orderBy = "weightFisher", ranksOf = "elimFisher")
            write.table(allRes,file=paste(workingdirectory,"/",filename,"/pattern",num-1,"_",ont[k],".txt",sep=""))
            
            KSsig<-as.matrix(allRes[ncol(allRes)])
        
            if(length(grep("<",KSsig)) >= 1)
            {
                position<-grep("<",KSsig)
                for(m in 1:length(position))
                {
                    KSsig2<-substr(KSsig[position[m]],2,nchar(KSsig[position[m]]))
                    KSsig[position[m]]<-KSsig2
                }
            }
            mode(KSsig)<-'numeric'
            KSsig2<-which(KSsig<=0.05)
            names_heat<-c()
            names_heat2<-c()
            aver3<-c()
            aver4<-c()
            min_length<-(ncol(alls)-12)*2
            mode(min_length)<-'numeric'
            
            datas<-read.delim(paste(workingdirectory,"/Functional_Searches",ont2[k],sep=""))
            datas<-as.matrix(datas)
            datas<-unname(datas)
            
            for(i in 1:length(KSsig2))
            {
                term<-allRes[i,1]
                positions<-which(datas[,2]==term)
                search_terms2<-c()
                aver<-c()
                aver2<-c()
                
                if(identical(positions,integer(0)))
                {
                    next
                }
                else
                {
                    for(j in 1:length(positions))
                    {
                        search_terms<-datas[positions[j],1]
                        search_terms2<-c(search_terms2,search_terms)
                        names_heat<-c(names_heat,allRes[i,2])
                        names_heat2<-c(names_heat2,search_terms)
                    }
                }
                
                search_terms2<-unique(search_terms2)
                
                for(j in 1:length(search_terms2))
                {
                    in_all<-grep(search_terms2[j],all_orthologs[,1])
                    if(identical(in_all,integer(0)))
                    {
                        next
                    }
                    else
                    {
                        aver<-rbind(aver,all_orthologs[in_all,12:(ncol(all_orthologs)-1)])
                    }
                }
                
                aver<-as.matrix(aver)
                
                for(j in 1:ncol(aver))
                    {
                    aver_col<-aver[,j]
                    aver_col<-as.numeric(aver_col)
                    aver1<-ave(aver_col)
                    aver2<-cbind(aver2,aver1[1])
                }
                
                aver3<-rbind(aver3,aver2)
                aver4<-rbind(aver4,aver)
                
            }
            mode(aver3)<-"numeric"
            mode(aver4)<-'numeric'
            names_heat<-unique(names_heat)
            
            if(nrow(aver3)<2 || length(aver3)==0)
            {
                next
            }
                pdf(paste(workingdirectory,"/",filename,"/","Significantly_enriched_GO_terms_",ont[k],".pdf",sep=""))
                heatmap.2(aver3,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.4,labRow=names_heat,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched GO terms",ont[k],sep=" ")," ",col.label))
                dev.off()
                
            write.table(names_heat,file=paste(workingdirectory,"/",filename,"/pattern",num-1,"_Terms_",ont[k],".txt",sep=""))
            if(nrow(aver4)<2 || length(aver4) == 0)
            {
                next
            }
                pdf(paste(workingdirectory,"/",filename,"/","Significantly_enriched_GO_terms_",ont[k],"_Genes.pdf",sep=""))
                heatmap.2(aver4,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.4,labRow=names_heat2,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched GO terms",ont[k],"Genes",sep=" ")," ",col.label))
                dev.off()
            
            if(identical(up_reg,1) && identical(down_reg,1))
            {
                return("No Mean Split Value")
            }
            
            if(identical(up_reg,1))
            {
                up_reg_matrix<-c()
                names_up<-c()
                down_reg_matrix<-c()
                names_down<-c()
                
                names_heat3<-unique(names_heat2)
                
                for(m in 1:length(names_heat3))
                {
                    is_there2<-grep(names_heat3[m],down_reg[,1])
                    if(length(is_there2)>0)
                    {
                        down_reg_matrix<-rbind(down_reg_matrix,down_reg[is_there2,2:ncol(down_reg)])
                        names_down<-c(names_down,down_reg[is_there2,1])
                    }
                    else
                    {
                        next
                    }
                    
                }
                mode(down_reg_matrix)<-'numeric'
                output2<-c()
                
                if(identical(down_reg_matrix,numeric(0)) || length(down_reg_matrix)==0)
                {
                    print(output2,paste(ont[k],"significant genes in (2) has one member",names_down,sep=" "))
                    down_reg_matrix<-matrix(0,1,1)
                }
                if(length(down_reg_matrix) > groups_length*2)
                {
                    down_reg_matrix<-as.matrix(down_reg_matrix)
                    down_reg_matrix<-unname(down_reg_matrix)
                    mode(down_reg_matrix)<-'numeric'
                    pdf(paste(workingdirectory,"/",filename,"/","Significantly_enriched_GO_terms_",ont[k],"_(1).pdf",sep=""))
                    heatmap.2(down_reg_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.4,labRow=names_down,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched GO terms",ont[k],"(1)",sep=" ")," ",col.label))
                }
            }
            
            if(identical(down_reg,1))
            {
                up_reg_matrix<-c()
                names_up<-c()
                down_reg_matrix<-c()
                names_down<-c()
                
                names_heat3<-unique(names_heat2)
                
                for(m in 1:length(names_heat3))
                {
                    is_there2<-grep(names_heat3[m],up_reg[,1])
                    if(length(is_there2)>0)
                    {
                        up_reg_matrix<-rbind(up_reg_matrix,up_reg[is_there2,2:ncol(up_reg)])
                        names_up<-c(names_up,up_reg[is_there2,1])
                    }
                    else
                    {
                        next
                    }
                    
                }
                up_reg_matrix<-as.matrix(up_reg_matrix)
                up_reg_matrix<-unname(up_reg_matrix)
                mode(up_reg_matrix)<-'numeric'
                output2<-c()
                
                if(identical(up_reg_matrix,numeric(0)) || length(up_reg_matrix)==0)
                {
                    print(output2,paste(ont[k],"significant genes in (2) has one member",names_up,sep=" "))
                    up_reg_matrix<-matrix(0,1,1)
                }
                if(length(up_reg_matrix) > groups_length*2)
                {
                    pdf(paste(workingdirectory,"/",filename,"/","Significantly_enriched_GO_terms_",ont[k],"_(2).pdf",sep=""))
                    heatmap.2(up_reg_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.6,labRow=names_up,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched GO terms",ont[k],"(2)",sep=" ")," ",col.label))
                    dev.off()
                }
            }
            
            else
            {
                up_reg_matrix<-c()
                names_up<-c()
                down_reg_matrix<-c()
                names_down<-c()
                other_reg_matrix<-c()
                names_other<-c()
                
                names_heat3<-unique(names_heat2)
                
                for(m in 1:length(names_heat3))
                {
                    is_there<-grep(names_heat3[m],up_reg[,1])
                    is_there2<-grep(names_heat3[m],down_reg[,1])
                    is_there3<-grep(names_heat3[m],other_reg[,1])
                    if(length(is_there)>0)
                    {
                        up_reg_matrix<-rbind(up_reg_matrix,up_reg[is_there,2:ncol(up_reg)])
                        names_up<-c(names_up,up_reg[is_there,1])
                    }
                    if(length(is_there2)>0)
                    {
                        down_reg_matrix<-rbind(down_reg_matrix,down_reg[is_there2,2:ncol(down_reg)])
                        names_down<-c(names_down,down_reg[is_there2,1])
                    }
                    if(length(is_there3)>0)
                    {
                        other_reg_matrix<-rbind(other_reg_matrix,other_reg[is_there3,2:ncol(other_reg)])
                        names_other<-c(names_other,other_reg[is_there3,1])
                    }
                    else
                    {
                        next
                    }
                    
                }
            
                mode(up_reg_matrix)<-'numeric'
                mode(down_reg_matrix)<-'numeric'
                mode(other_reg_matrix)<-'numeric'
                output1<-c()
                output2<-c()
                output3<-c()
            
                if(identical(up_reg_matrix,numeric(0)) || length(up_reg_matrix)==0)
                {
                    print(paste(ont[k],"significant genes in (1) has one member",names_up,sep=" "))
                    up_reg_matrix<-matrix(0,1,1)
                }
                if(length(up_reg_matrix) > groups_length*2 && length(up_reg_matrix)>groups_length*2)
                {
                   up_reg_matrix<-as.matrix(up_reg_matrix)
                    up_reg_matrix<-unname(up_reg_matrix)
                    mode(up_reg_matrix)<-'numeric'
                    pdf(paste(workingdirectory,"/",filename,"/","Significantly_enriched_GO_terms_",ont[k],"_(_1).pdf",sep=""))
                   heatmap.2(up_reg_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.4,labRow=names_up,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched GO terms",ont[k],"(1)",sep=" ")," ",col.label))
                    dev.off()
                }
                
                if(identical(down_reg_matrix,numeric(0)) || length(down_reg_matrix)==0)
                {
                     print(paste(ont[k],"significant genes in (2) has one member",names_down,sep=" "))
                     down_reg_matrix<-matrix(0,1,1)
                }
                if(length(down_reg_matrix) > groups_length*2 && length(down_reg_matrix) > groups_length*2)
                {
                    down_reg_matrix<-as.matrix(down_reg_matrix)
                    down_reg_matrix<-unname(down_reg_matrix)
                    mode(down_reg_matrix)<-'numeric'
                    pdf(paste(workingdirectory,"/",filename,"/","Significantly_enriched_GO_terms_",ont[k],"_(2).pdf",sep=""))
                    heatmap.2(down_reg_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.4,labRow=names_down,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched GO terms",ont[k],"(2)",sep=" ")," ",col.label))
                    dev.off()
                }
                
                if(identical(other_reg_matrix,numeric(0)) || length(other_reg_matrix)==0)
                {
                     print(paste(ont[k],"significant genes in (3) has one member",names_down,sep=" "))
                     other_reg_matrix<-matrix(0,1,1)
                }
                if(length(other_reg_matrix) > groups_length*2 && length(other_reg_matrix)> groups_length*2)
                {
                    other_reg_matrix<-as.matrix(other_reg_matrix)
                    other_reg_matrix<-unname(other_reg_matrix)
                    mode(other_reg_matrix)<-'numeric'
                    pdf(paste(workingdirectory,"/",filename,"/","Significantly_enriched_GO_terms_",ont[k],"_(3).pdf",sep=""))
                    heatmap.2(other_reg_matrix,col=redgreen(75),scale="row",key=T,keysize=1.5,density.info="none",trace="none",cexCol=0.9,cexRow=0.4,labRow=names_other,labCol=1:groups_length,ColSideColors=col.map,main=c(paste("pattern ",num-1," Significantly Enriched GO terms",ont[k],"(2)",sep=" ")," ",col.label))
                    dev.off()
                }
                
                write.table(names_up,file=paste(workingdirectory,"/",filename,"/pattern",num-1,"(1)_up_",ont[k],".txt",sep=""))
                write.table(names_down,file=paste(workingdirectory,"/",filename,"/pattern",num-1,"(2)_down_",ont[k],".txt",sep=""))
                write.table(names_other,file=paste(workingdirectory,"/",filename,"/pattern",num-1,"(3)_other_",ont[k],".txt",sep=""))
                
            }
        }
    }    
}