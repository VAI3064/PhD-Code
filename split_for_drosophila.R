split_for_drosophila<- function(genes_of_interest)
{
	AGAPs<- genes_of_interest[,1]
	
	list1<- ""
	
	for ( i in 1: length(AGAPs))
	{
			if(nchar(AGAPs[i])<10)
		{
			list1<-c(list1,AGAPs[i])
			next
		}
		else
		{
			splits<- strsplit(AGAPs[i], "-")[[1]][1]
			list1<- c(list1, splits)
		}
	}
	list1 <- list1[2:length(list1)]
	genes_of_interest[,1]<- list1
	
	return(genes_of_interest)
}