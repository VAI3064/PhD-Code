fold_changes2<- function (x2, genes_list,j)
{
	fold_changes_matrix<- ""
	rows<- row.names(x2)
        
        for(i in 1:length(genes_list))
        {
            f<-grep(genes_list[i],rows,value=FALSE)
		fold_changes_matrix<- rbind(fold_changes_matrix, c(x2[f[1],1:(max(j)+1)]))	
        }
    
    return(fold_changes_matrix)
}
