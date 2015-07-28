function(all.fold.changes)
{
	#Outputs a file that can be read into Cytoscape and visualised as a network of co-correlated nodes#
	expr<-read.delim(all.fold.changes,header=T)
	expression.1<-expr[,-1]
	expression.1<-as.matrix(expression.1)
	mode(expression.1)<-'numeric'
	cor.matrix<-cor(t(expression.1))
	cor.matrix.2<-cor.matrix
	diag(cor.matrix.2)<-0
	cor.matrix.2[lower.tri(cor.matrix.2)]<-0
	cor.matrix.2<-abs(cor.matrix.2)
	positions<-which(cor.matrix.2>0.8,arr.ind=T)
	
	vector.out<-c()
	for(i in 1:nrow(positions))
	{
 		 vector.out<-rbind(vector.out,c(as.character(expr[positions[i,1],1]),as.character(expr[positions[i,2],1]),cor.matrix[positions[i,1],positions[i,2]]))
	}
	return(vector.out)
}
