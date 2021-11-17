gen_bspline<-function(tau,breaks,degree)
{
  # generate b-splines for tau
  # tau : quantile level seq
  # breaks:
  # degree:
  splines_tau = bsplineS(tau,breaks,degree)
  return(splines_tau)
}

phi<-function(y,z_tau,tau,b,n)
{
  # calculate phi for test statistic
  kn = length(tau)-1
  Tau = matrix(rep(tau,each=n),nrow = n)
  psi = y - z_tau%*%b
  psi = matrix(psi,nrow = n)
  psi = Tau - (psi<0)
  return(rowSums(psi))
}

