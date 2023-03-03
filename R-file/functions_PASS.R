DataSplit <- function(data, n, n_test, n_cal, n_rest){
  if(n_test>0){
    index_test <- sample(1:n, n_test, replace = FALSE)
    data_test <- data[index_test,]
    data_rest2 <- data[-index_test,]
    data_rest <- data_rest2[sample(1:dim(data_rest2)[1], n_rest),]
    index_cal <- sample(1:dim(data_rest)[1], n_cal)
    data_train <- data_rest[-index_cal,]
    data_cal <- data_rest[index_cal,]
    return(list(data_train = data_train, data_cal = data_cal, data_test = data_test, data_rest = data_rest))
  }else{
    data_rest <- data[sample(1:dim(data)[1], n_rest, replace = FALSE),]
    index_cal <- sample(1:dim(data_rest)[1], n_cal)
    data_train <- data_rest[-index_cal,]
    data_cal <- data_rest[index_cal,]
    return(list(data_train = data_train, data_cal = data_cal, data_test = 0, data_rest = data_rest))
  }
}


confomalPvalue <- function(W_cal, W_test, Null_cal, Value){
  Phi_cal <- -ScoreCompute(W_cal, Value)
  Phi_test <- -ScoreCompute(W_test, Value)
  Phi_Null <- Phi_cal[Null_cal]
  n2 <- length(Phi_Null)
  pvalue <- sapply(Phi_test, function(t){
    (sum(Phi_Null<=t)+1)/(n2+1)
  })
  return(pvalue)
}


NullIndex <- function(y, Value){
  if(Value$type=="==A"){
    index <- which(y==Value$v)
  }else if(Value$type=="<=A"){
    index <- which(y<=Value$v)
  }else if(Value$type==">=B"){
    index <- which(y>=Value$v)
  }else if(Value$type=="<=A|>=B"){
    index <- which(y<=Value$v[1]|y>=Value$v[2])
  }else if(Value$type==">=A&<=B"){
    index <- which(y>=Value$v[1]&y<=Value$v[2])
  }
  return(index)
}


ClassPred <- function(pred, algo){
  if(algo@name=='RFc'|algo@name=='NN'){
    return(2*as.numeric(pred>0.5)-1)
  }else if(algo@name=='SVM'){
    return(2*as.numeric(pred>0)-1)
  }
}


ScoreCompute <- function(pred, Value){
  if(Value$type=="==A"|Value$type=="<=A"){
    Phi <- pred
  }else if(Value$type==">=B"){
    Phi <- -pred
  }else if(Value$type=="<=A|>=B"){
    Phi <- pmin(pred-Value$v[1], Value$v[2]-pred)
  }else if(Value$type==">=A&<=B"){
    Phi <- pmax(Value$v[1]-pred, pred-Value$v[2])
  }
  return(Phi)
}


SolvePi <- function(TN, m, S, Screened, alpha, constraint = TRUE){
  c <- rep(0, length(Screened))
  if(constraint){
    dT <- diag((1-TN)[Screened])
    H <- 2*dT%*%S[Screened, Screened]%*%dT
    A <- rbind(rep(1, length(Screened)), TN[Screened])
    b <- c(1, 0)
    r <- c(0, alpha)
  }else{
    H <- 2*S[Screened, Screened]
    A <- rep(1, length(Screened))
    b <- c(1)
    r <- c(0)
  }
  l <- rep(0, length(Screened))
  u <- rep(1/m, length(Screened))
  
  Solu <- ipop(c, H, A, b, l, u, r)
  
  Pi <- Solu@primal
  Pi[Pi<1e-8] <- 0
  Pi <- Pi/sum(Pi)
  
  return(Pi)
}


Tcompute <- function(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = FALSE, IsCali = FALSE){
  if(h1==0){
    h1 <- density(W_cal[Null_cal])$bw
  }
  if(h2==0){
    h2 <- density(W_cal)$bw
  }
  if(IsSame){
    h1 <- h2
  }
  f0 <- kde(W_cal[Null_cal], h = h1, eval.points = W_test)
  f <- kde(W_cal, h = h2, eval.points = W_test)
  phat <- 1-length(Null_cal)/length(W_cal)
  TN <- (1-phat)*f0$estimate/f$estimate
  TN[which(TN>1)] <- 1
  TN[which(TN<0)] <- 0
  if(IsCali){
    Ind <- order(W_test)
    TNW <- TN[Ind]
    for (i in (length(TNW)-20):1) {
      if(TNW[i]<TNW[i+1]){TNW[i] <- TNW[i+1]}
    }
    TN[Ind] <- TNW
  }
  return(TN)
}


KfoldTcompute <- function(X_rest, Y_rest, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = FALSE, IsCali = FALSE){
  Folds <- createFolds(1:length(Y_rest), K)
  m <- length(W_test_whole)
  W_cross <- rep(0, length(Y_rest))
  for (i in 1:K) {
    model <- fitting(algo, X_rest[-Folds[[i]],], Y_rest[-Folds[[i]]], lambda = lambda)
    lens <- length(Folds[[i]])
    if(lens==1){
      W_cross[Folds[[i]]] <- Pred(algo, model, t(X_rest[Folds[[i]],]))
    }else{
      W_cross[Folds[[i]]] <- Pred(algo, model, X_rest[Folds[[i]],])
    } 
  }
  TN_cross <- Tcompute(W_cross, W_test_whole, Null_rest, algo, h1, h2, IsSame, IsCali)
  return(TN_cross)
}


L2Normsq <- function(x){
  return(sum(x^2))
}


Scompute <- function(g, Z_test, d){
  if(g=='RBF'){
    S <- Z_test%*%t(Z_test)
    for (i in 1:dim(Z_test)[1]){
      if(d>1){
        S[i,] <- apply((Z_test[i,]-t(Z_test)), 2, L2Normsq)
      }
      if(d==1){
        S[i,] <- apply((Z_test[i]-t(Z_test)), 2, L2Normsq)
      }
    }
    S <- exp(-S)
  }
  if(g=='COS'){
    S <- (Z_test/sqrt(apply(Z_test, 1, L2Normsq)))%*%t(Z_test/sqrt(apply(Z_test, 1, L2Normsq)))
  }
  return(S)
}


CriterionCompute <- function(Pi, Index, Null_test, m, S, Retime = 1){
  record <- sapply(1:Retime, function(i){
    Sampled <- sample(Index, m, replace = TRUE, Pi)
    FDP <- sum(is.element(Sampled, Null_test))/m
    H1 <- Sampled[!is.element(Sampled, Null_test)]
    SIM <- sum(S[H1, H1])/sum(!is.element(Sampled, Null_test))^2
    return(c(FDP, SIM))
    })
  result <- apply(record, 1, mean)
  return(data.frame(FDP = result[1], SIM = result[2]))
}