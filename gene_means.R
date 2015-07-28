gene_means<- function (x)
{
	means<- rep(1,length(x[,1]))
for ( i in 1: length(x[,1]))
{
	means[i]<-(mean(x[i,]))
}
t_means<- t(means)
}
