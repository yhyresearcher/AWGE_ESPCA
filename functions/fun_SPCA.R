SPCA = function(X, k=5, kv=1, niter=1000, err=0.0001, Num.init=5){
  # --------------------------------------------
  # [n,p] = dim(X), n is the number of samples and p is the number of features
  # --------------------------------------------
  n = nrow(X); # n is the number of samples
  p = ncol(X); # p is the number of features
  print(k)
  if(length(kv)==1) kv = rep(kv,k)
  
  U=matrix(0,n,k);D=matrix(0,k,k);V=matrix(0,p,k)
  tX = X
  out = rank1.SPCA(tX, kv[1], niter, err, Num.init)
  U[,1] = out$u; V[,1] = out$v; D[1,1] = out$d
  if(k<2) return (list(U=U, D=D, V=V))
  # --------------------------------------------
  for(i in 2:k){
    tX = tX-c(out$d)*out$u%*%t(out$v); 
    UU = U%*%t(U)
    out = cycleFun1(tX, UU, kv[i], niter, err, Num.init)
    U[,i] = out$u; V[,i] = out$v; D[i,i] = out$d
  }
  return (list(U=U, D=D, V=V))
}
# ----------------------------------------------
rank1.SPCA = function(X, kv, niter=1000, err=0.0001, Num.init = 5){
  n = nrow(X); # n is the number of samples
  p = ncol(X); # p is the number of features
  d.opt = -100
  # set initial point
  for(ii in 1:Num.init){
    set.seed(ii*100)
    v0 = matrix(rnorm(p,0,1),ncol=1);v0 = v0/norm(v0,'E')
    u0 = matrix(rnorm(n,0,1),ncol=1);u0 = u0/norm(u0,'E')
    # Iterative algorithm to solve u and v values
    for (i in 1:niter){
      
      u = u.project1(X%*%v0) 
      v = l0project(t(X)%*%u, kv)
      
      # Algorithm termination condition norm(matrix(v),"E")
      if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)){break}
      else {
        u0 = u;v0 = v}
    }
    d =t(u)%*%X%*%v
    if(d>d.opt){
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}
# ----------------------------------------------
cycleFun1 = function(X, UU, kv, niter=1000, err=0.0001, Num.init = 5){
  
  n = nrow(X); # n is the number of samples
  p = ncol(X); # p is the number of features
  d.opt = -100
  # set five initial point
  for(ii in 1:Num.init){
    set.seed(ii*100)
    v0 = matrix(rnorm(p,0,1),ncol=1);v0 = v0/norm(v0,'E')
    u0 = matrix(rnorm(n,0,1),ncol=1);u0 = u0/norm(u0,'E')
    
    # Iterative algorithm to solve u and v values
    for (i in 1:niter){
      u = (diag(n) - UU)%*%(X%*%v0); u = u.project1(u)
      v = l0project(t(X)%*%u, kv)
      
      # Algorithm termination condition norm(matrix(v),"E")
      if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)){break}
      else {
        u0 = u;v0 = v}
    }
    d =t(u)%*%X%*%v
    if(d>d.opt){
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}
# ----------------------------------------------
u.project1 = function(z){  
  u = z
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    return(u)} 
}

select = function(x, k){
  if(k>=length(x)) return(x)
  x[-order(x,decreasing=T)[1:k]] = 0
  return(x)
}

l0project = function(z, k){  
  absz = abs(z);
  u = select(absz,k)
  
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    u = sign(z)*u
    return(u)} 
}
# ----------------------------------------------