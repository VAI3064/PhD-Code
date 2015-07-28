colouring <- function(expr)
{
    colour<-c()
    j <- 1
    colvec1<-colors()
    colvec2<-colvec1[17:140]
    colvec2<-c(colvec2,colvec1[365:652])
    colvec<-sample(colvec2,length(colvec2))
    for(i in 1:(length(expr)-1))
    {
        if(expr[i]==expr[i+1])
        {
            colour<-c(colour,colvec[j])
        }
        else
        {
            colour<-c(colour,colvec[j])
            j<-j+1
        }
        if(i == (length(expr)-1))
        {
            colour<-c(colour,colvec[j])
        }
    }
    return(colour)
}