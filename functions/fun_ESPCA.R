ESPCA = function(X, k=2, overlap.group, k.group=2, we=0.5, t = 0.1, niter=1000){
  # --------------------------------------------------------------------------
  # [n,p] = dim(X), n is the number of samples and p is the number of features
  # --------------------------------------------------------------------------
  n = nrow(X); # n is the number of samples
  p = ncol(X); # p is the number of features
  U = matrix(0,n,k); D = matrix(0,k,k); V = matrix(0,p,k)
  tX = X
  out = rank1.ESPCA(tX, overlap.group, k.group, we, t, niter)
  U[,1] = out$u; V[,1] = out$v; D[1,1] = out$d 
  if(k<2) return (list(U=U, D=D, V=V))
  # --------------------------------------------------------------------------
  
  for(i in 2:k){
    tX = tX-c(out$d)*out$u%*%t(out$v); 
    UU = U%*%t(U)
    #(we)
    out = cycleFun2(tX, UU, overlap.group, k.group,we, t, niter)
    U[,i] = out$u; V[,i] = out$v; D[i,i] = out$d
  }
  return (list(U=U, D=D, V=V))
}
# ----------------------------------------------------------------------------
rank1.ESPCA = function(X, overlap.group, k.group, we, t, niter=1000){
  n = nrow(X); # n is the number of samples
  p = ncol(X); # p is the number of features
  d.opt = -100
  # set five initial point
  for(ii in 1:5){
    we1 = we
    print("rank")
    set.seed(ii*100)
    v0 = matrix(rnorm(p,0,1),ncol=1);v0 = v0/norm(v0,'E')
    u0 = matrix(rnorm(n,0,1),ncol=1);u0 = u0/norm(u0,'E')
    # Iterative algorithm to solve u and v values
    for (i in 1:niter){
      u = u.project2(X%*%v0)
      v = overlap.group.penalty(t(X)%*%u, overlap.group, k.group, we1)
      if(we1 > 0){
        we1 = we1 -t
      }else{
        we1 = 0
      }
      # Algorithm termination condition norm(matrix(v),"E")
      if ((norm(u - u0,'E')<= 0.0001)&(norm(v - v0,'E')<= 0.0001)){break}
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
cycleFun2 = function(X, UU, overlap.group, k.group,we, t, niter=1000){
  n = nrow(X); # n is the number of samples
  p = ncol(X); # p is the number of features
  d.opt = -100
  for(ii in 1:5){
    print("cyc")
    we1 = we
    set.seed(ii*100)
    v0 = matrix(rnorm(p,0,1),ncol=1);v0 = v0/norm(v0,'E')
    u0 = matrix(rnorm(n,0,1),ncol=1);u0 = u0/norm(u0,'E')
    # Iterative algorithm to solve u and v values
    for(i in 1:niter){
      u = (diag(n) - UU)%*%(X%*%v0); u = u.project2(u)
      v = overlap.group.penalty(t(X)%*%u, overlap.group, k.group, we1)
      if(we1 > 0){
        we1 = we1 -t
      }else{
        we1 = 0
      }
      # Algorithm termination condition norm(matrix(v),"E")
      if ((norm(u - u0,'E')<= 0.0001)&(norm(v - v0,'E')<= 0.0001)){break}
      else {
        u0 = u;v0 = v}
    }
    d = t(u)%*%X%*%v
    if(d>d.opt){
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}
# ----------------------------------------------
u.project2 = function(z){  
  u = z
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    return(u)} 
}
# ----------------------------------------------
overlap.group.penalty = function(u, overlap.group, k0 = 1.0, we=0.5){
  # Greedy k-edges-sparse projection
  # k0 : the number of varibales
  # overlap.group edges的集合
  #-----------------------------
  group.num = length(overlap.group)
  #求边的个数
  group.norm = rep(0,group.num)
  #建立二维矩阵，q权重，边
  for(i in 1:group.num){
    g.set = overlap.group[[i]]
    w_i = 1/sqrt(length(g.set))
    #print(g.set)
    group.norm[i] = norm(w_i*u[g.set],"2")
  }
  if(we > 0){
    k1 = k0*(1+we)
    k1 = ceiling(k1)
  }else{
    k1 = k0
  }
  if(k1 < k0){
    k1 = k0
  }
  #print(k1)
  if(we != 0){
    ID1 = order(group.norm, decreasing = TRUE)[1:k1]
    ID = sample(ID1, k0)
  }else {
    ID = order(group.norm, decreasing = TRUE)[1:k0]
  }
  
  #print("k number")
  #print(ID)
  #从大到小排序，保留K0个边
  select.features = c()
  #c() 建立一个空的向量集合
  for(i in 1:length(ID)){
    temp = overlap.group[[ID[i]]]
    select.features = c(select.features, temp)
  }
  #使用向量记录保留边的ID
  index= sort(unique(select.features))
  #去除可能重复的元素
  x.opt = u; x.opt[-c(index)] = 0
  if(sum(abs(x.opt))==0) return(x.opt)
  else return(x.opt/norm(x.opt,"E")) 
}