colourmatrix <- function(colour,number)
{
    colmatrix<-c()
    for(i in 1:length(colour))
    {
        colmatrix<-c(colmatrix,paste(colour[i]," = ",number[i]))
    }
return(colmatrix)
}