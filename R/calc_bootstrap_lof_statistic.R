Un.b.c<-function(x,x.b,phi,phi.b)
{
  # calculate the bootstrap lack-of-fit test statistic via C code
  # x: original sample
  # x.b: bootstrap sample of x
  # phi: results of funciton Phi with respect to x, in vector form
  # phi.star: results of function Phi with respect to x.star, in vector form
  # in this function, we use ||x_j-x_i||=sqrt(sum((x_j-x_i)^2))
  # in this function, we use c code to compute the result
  x<-as.matrix(x)
  x.b<-as.matrix(x.b)
  n<-nrow(x)
  p<-ncol(x)
  res<-.Call("unbc",x,x.b,phi,phi.b,as.integer(n),as.integer(p))
  return(res)
}
