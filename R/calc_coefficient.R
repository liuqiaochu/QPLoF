comp.B.MM.c<-function(ztau,Y,tau,breaks,basis.order,n,p,maxiter=1000,tol=10^-8,
                      epsilon=0.001)
{
  # compute the coefficient matrix B of B-spline using MM algorithm via C code
  # x,y: original data matrix
  # tau: quantile level vector
  # breaks: knots of B-spline
  # basis.order: order of B-spline
  # maxiter: max number of iteration of MM algorithm
  # tol: tolerance of convergence
  # epsilon: default value set to be 0.001
  kn<-length(tau)-1   # kn: k_n in the paper
  m<-length(breaks)+basis.order-2   # m: number of basis functions

  # prepare the data of standard form for MM algorithm
  theta<-solve(t(ztau)%*%ztau)%*%t(ztau)%*%Y   # theta: initial point of theta

  # call the MM algorithm
  theta<-.Call("compbmmc",ztau,Y,rep(tau,each=n),theta,as.integer(n),
               as.integer(p*m),as.integer(kn+1),as.integer(maxiter),tol,epsilon)
  R<-Y-ztau%*%theta

  B<-matrix(theta,nrow=p)
  R<-matrix(R,ncol=kn+1)
  res<-list(B=B,Residual=R)
  return(res)
}
