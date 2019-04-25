#This is the function Iterat that find the m-tile of a sample
#L is the name of a sample and m is between 2 and m+1
#The function first sort the sample in ascending order
#If m =2 the iterat function just estimate the median
#if m=m+1 the iterat function estimate the mean
#if m is greater than 2 and less or equal to  m it estimate the m-tile
#The m-tile method uses the percentile method in picking the elements.
#L is the sample in a vector form eg. L=c(3,5,4,1)
#y is the sorted sample 
#n is the length of the sample
#u is the percentile 
#o gives the position of the elements to pick 
#q is the mean of the elements picked from the sample


iterat<-function(L, m){ 
  y<-sort(L)
  n<-length(L)
  u=100/(m-1)      
  o=round((u/100)*n)
  if (m==2){
    (a=median(L))
    return(a)
  }
  else if (m>2& m<=(n)){
    q<-y[seq((o), n, (o))] 
    a=(sum(q)/(length(q))) 
    return(a)
  }
  else if (m==n+1){
    a<-(mean(L))
    return(a)}
  else print("m is not a valid number for the length of your sample")
  
}  

#Examples 1
iterat(L=c(5,2,8,4,1,0), 2)

# Example 2 
L=c(5,8,4,7,2,1,3)
iterat(L, m=5)





