Un<-function(x,phi)
{
  # calculate our lack-of-fit test statistic
  n = nrow(x)
  ind = combn(n,2)
  res<-apply(ind,2,function(ind,x,phi) {
    u<-prod(phi[ind])
    v<-(x[ind[1],]-x[ind[2],])^2
    v<-sqrt(sum(v))
    return(u*v)
  },x,phi)
  res = -sum(res)/ncol(ind)
  return(res)
}
