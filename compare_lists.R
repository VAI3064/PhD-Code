compareLists<-function(list1, list2)
{
	list1<-as.vector(list1)
	list2<-as.vector(list2)
	
	inCommon<-""
	
	for (i in 1:length(list1)){
		for (j in 1:length(list2)){
			if (list1[i] == list2[j]){
				inCommon<-c(inCommon,list1[i])
				break
			}
		}
	}
	return (inCommon)
}
