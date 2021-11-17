QPLoF<-function(X,y,n,p,beta0,gamma0,epsilon,tau,breaks,degree,B)
{
  # calculate the p-value of our lack-of-fit test through paired bootstrap
  Y = rep(y,length(tau))
  Phi_tau = gen_bspline(tau,breaks,degree)
  z_tau = kronecker(Phi_tau,X)
  result = comp.B.MM.c(z_tau,Y,tau,breaks,degree,n,p,maxiter=200,tol=10^-8,epsilon=0.01)
  B_hat = result$B
  b_hat = as.vector(B_hat)
  e_hat = result$Residual
  Phi = phi(Y,z_tau,tau,b_hat,n)
  U_star = Un(X,Phi)
  U_boot = vector(length = B)
  for(i in 1:B)
  {
    # resample by paired bootstrap
    Ind<-sample(1:n,n,replace = T)
    x_boot<-as.matrix(X)[Ind,]
    y_boot<-y[Ind]
    Y_boot = rep(y_boot,length(tau))
    Phi_tau = gen_bspline(tau,breaks,degree)
    z_tau_boot = kronecker(Phi_tau,x_boot)
    # calculate the process coefficient of the bootstrap samples
    B_boot<-comp.B.MM.c(z_tau_boot,Y_boot,tau,breaks,degree,n,p,maxiter=200,tol=10^-8,epsilon=0.01)$B
    Phi_boot = phi(Y_boot,z_tau_boot,tau,as.vector(B_boot),n)
    U_boot[i]<-Un.b.c(as.matrix(X),as.matrix(x_boot),Phi,Phi_boot)
  }
  pvalue = sum(U_star<U_boot)/B
  return(pvalue)
}
