AWGE = function(X, k=2, overlap.group, k.group=2, we=0.5, t = 0.1, niter=1000){
  # --------------------------------------------------------------------------
  # [n,p] = dim(X), n is the number of samples and p is the number of features
  # --------------------------------------------------------------------------
  n = nrow(X); # n is the number of samples
  p = ncol(X); # p is the number of features
  U = matrix(0,n,k); D = matrix(0,k,k); V = matrix(0,p,k)
  tX = X
  out = rank1.AWGE(tX, overlap.group, k.group, we, t, niter)
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
rank1.AWGE = function(X, overlap.group, k.group, we, t, niter=1000){
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
boxplot_outlier_correction = function(X, LQ = 0.35){
  Q1 <- quantile(X, probs=LQ)
  Q3 <- quantile(X, probs=1-LQ)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1 * IQR
  upper_bound <- Q3 + 0.5 * IQR
  X[X < lower_bound] <- 0
  X[X > upper_bound] <- 0
  return(X)
}
# ?޸ĺ??? overlap.group.penalty ????
# 加载权重
load_weights = function(filepath) {
  weights = as.matrix(read.table(filepath, header = FALSE))
  return(weights)
}

# 应用权重到数据
apply_weights <- function(data, weights) {
    # 确保数据是矩阵
    if (!is.matrix(data)) {
        stop("Data should be a matrix.")
    }
    
    # 确保权重是单列矩阵
    if (!(is.matrix(weights) && ncol(weights) == 1)) {
        stop("Weights should be a single-column matrix.")
    }
    
    # 获取数据的列数
    num_columns <- ncol(data)
    
    # 创建权重矩阵，复制每行以匹配数据的列数
    weights_matrix <- matrix(weights, nrow = nrow(data), ncol = num_columns, byrow = TRUE)
    
    # 应用权重到数据，逐元素乘法
    weighted_data <- data * weights_matrix
    
    return(weighted_data)
}


# 更新权重以匹配剩余数据
update_weights = function(weights, removed_indices) {
  weights[-removed_indices, , drop = FALSE]
}
# 修改后的overlap.group.penalty函数
overlap.group.penalty = function(u, overlap.group, k0 = 1.0, we=0.5){
  # 读取权重
  weights = load_weights("power.txt")

  # 使用箱形图算法对u数据进行修正
  u_corrected = boxplot_outlier_correction(u)

  # 应用权重
  u_weighted = apply_weights(u_corrected, weights)
  # ??????ԭ?еĺ???????
  group.num = length(overlap.group)
  group.norm = rep(0,group.num)
  for(i in 1:group.num){
    g.set = overlap.group[[i]]
    w_i = 1/sqrt(length(g.set))
    group.norm[i] = norm(w_i*u_weighted[g.set],"2")
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
  #?Ӵ???С???򣬱???K0????
  select.features = c()
  #c() ??  һ???յ???  ????
  for(i in 1:length(ID)){
    temp = overlap.group[[ID[i]]]
    select.features = c(select.features, temp)
  }
  #ʹ????  ??¼?????ߵ?ID
  index= sort(unique(select.features))
  #ȥ???????ظ???Ԫ??
  x.opt = u; x.opt[-c(index)] = 0
  if(sum(abs(x.opt))==0) return(x.opt)
  else return(x.opt/norm(x.opt,"E")) 
}
