#' Solving the quadratic optimization problem for the sampling probability vector
#'
#' @param TN Vector. Estimated local FDR of the unlabeled individuals.
#' @param m Integral. The sampling budget, that is, how many unlabeled individuals should be sampled.
#' @param S Matrix. The similarity matrix of the unlabeled covariates.
#' @param Screened Vector. Indices of the screened individuals.
#' @param alpha Numeric. The pre-specified error rate level.
#'
#' @return The subsampling probability vector.
#' @export
SolvePi <- function(TN, m, S, Screened, alpha){
  c <- rep(0, length(Screened))
  dT <- diag((1-TN)[Screened])
  H <- 2*dT%*%S%*%dT
  A <- rbind(rep(1, length(Screened)), TN[Screened])
  b <- c(1, 0)
  r <- c(0, alpha)
  l <- rep(0, length(Screened))
  u <- rep(1/m, length(Screened))

  Solu <- kernlab::ipop(c, H, A, b, l, u, r)

  Pi <- Solu@primal
  Pi[Pi<1e-8] <- 0
  Pi <- Pi/sum(Pi)

  return(Pi)
}



#' Estimating the local FDR of unlabeled individuals
#'
#' @param W_cal Vector. Prediction values on the calibration data given by the model trained on the training set.
#' @param W_test Vector. Prediction values on the unlabeled data given by the model trained on the training set.
#' @param Null_cal Vector. Indices of the individuals not of interest among the calibration data.
#' @param classification Logical. Should be TRUE if under classification settings. FALSE as default.
#'
#' @return Vector containing the estimated local FDR of the unlabeled individuals.
#' @export
Tcompute <- function(W_cal, W_test, Null_cal, classification = FALSE){
  h1 <- stats::density(W_cal[Null_cal])$bw
  h2 <- stats::density(W_cal)$bw
  if(classification){
    h1 <- h2
  }
  f0 <- ks::kde(W_cal[Null_cal], h = h1, eval.points = W_test)
  f <- ks::kde(W_cal, h = h2, eval.points = W_test)
  phat <- 1-length(Null_cal)/length(W_cal)
  TN <- (1-phat)*f0$estimate/f$estimate
  TN[which(TN>1)] <- 1
  TN[which(TN<0)] <- 0
  if(classification){
    Ind <- order(W_test)
    TNW <- TN[Ind]
    for (i in (length(TNW)-20):1) {
      if(TNW[i]<TNW[i+1]){TNW[i] <- TNW[i+1]}
    }
    TN[Ind] <- TNW
  }
  return(TN)
}



#' Computing the similarity matrix
#'
#' @param g Character. The criterion function used for measuring similarity. Take "RBF" or "COS" and "RBF" as default.
#' @param Z_test Matrix. Normalized covariate matrix of the unlabeled data.
#' @param d Integral. Dimensionality of the covariates.
#'
#' @return The similarity matrix of the unlabeled data set.
#' @export
Scompute <- function(g, Z_test, d){
  if(g=='RBF'){
    S <- Z_test%*%t(Z_test)
    for (i in 1:dim(Z_test)[1]){
      if(d>1){
        S[i,] <- apply((Z_test[i,]-t(Z_test)), 2, function(x){sum(x^2)})
      }
      if(d==1){
        S[i,] <- apply((Z_test[i]-t(Z_test)), 2, function(x){sum(x^2)})
      }
    }
    S <- exp(-S)
  }
  if(g=='COS'){
    S <- (Z_test/sqrt(apply(Z_test, 1, function(x){sum(x^2)})))%*%t(Z_test/sqrt(apply(Z_test, 1, function(x){sum(x^2)})))
  }
  return(S)
}



#' Estimating the local FDR of unlabeled individuals under PASSC
#'
#' @param X_label Matrix. The covariate matrix of the labeled data.
#' @param Y_label Vector. Responses of the labeled data. Should take +-1 for binary Y.
#' @param X_test Matrix. The covariate matrix of the unlabeled data.
#' @param Null_label Vector. Indices of the individuals not of interest among the labeled data.
#' @param W_test_whole Vector. Prediction values on the unlabeled data given by the model trained on the whole labeled data set without splitting.
#' @param K Integral. Number of folds used for PASSC. Five as default.
#' @param lambda Numeric. The tuning parameter used for fitting the predictive model.
#' @param algo Class. A class specifying the algorithms used for fitting and prediction. The class should be well-defined with two functions. One named "fitting", which takes algorithm class, covariate matrix, responses and possibly a tuning parameter as inputs and outputs a model object. The other named "Pred", which takes algorithm class, a model object and a new covariate matrix as inputs and outputs a vector of predciton values.
#' @param classification classification Logical. If the response is discrete then set this to TRUE. FALSE as default.
#' @param shift Logical. If there is possible distribution shift between labelled and unlabelled data then set this to TRUE. FALSE as default.
#' @param region List. Should contain two attributes. The first attribute named 'type' takes values among "==A", "<=A", ">=B", "<=A|>=B" and ">=A&<=B". This characterizes the shape of the null region, where "==A" corresponds to classification settings. The second attribute named "v" characterizes the boundary of the null region. It should be a single number or a vector of length 2. If Y is binary then set v=-1.
#' @param tau Numeric. Threshold of the Storey estimation of the signal proportion under distribution shift.
#'
#' @return Vector containing the estimated local FDR of the unlabeled individuals under PASSC.
#' @export
KfoldTcompute <- function(X_label, Y_label, X_test, Null_label, W_test_whole, K, lambda, algo, classification = FALSE, shift = FALSE, region = NULL, tau = 0.5){
  Folds <- caret::createFolds(1:length(Y_label), K)
  W_cross <- rep(0, length(Y_label))
  for (i in 1:K) {
    model <- fitting(algo, X_label[-Folds[[i]],], Y_label[-Folds[[i]]], lambda = lambda)
    lens <- length(Folds[[i]])
    if(lens==1){
      W_cross[Folds[[i]]] <- Pred(algo, model, t(X_label[Folds[[i]],]))
    }else{
      W_cross[Folds[[i]]] <- Pred(algo, model, X_label[Folds[[i]],])
    }
  }
  if(shift){
    phat <- pEstimate(W_cross, W_test_whole, Null_label, region, tau)
    TN_cross <- Tcompute_shift(W_cross, W_test_whole, Null_label, classification, phat)
  }else{
    TN_cross <- Tcompute(W_cross, W_test_whole, Null_label, classification)
  }
  return(TN_cross)
}



#' Estimating the local FDR of unlabeled individuals under distribution shift
#'
#' @param W_cal Vector. Prediction values on the calibration data given by the model trained on the training set.
#' @param W_test Vector. Prediction values on the unlabeled data given by the model trained on the training set.
#' @param Null_cal Vector. Indices of the individuals not of interest among the calibration data.
#' @param classification Logical. Should be TRUE if under classification settings. FALSE as default.
#' @param phat Numeric. Estimated value of the signal proportion.
#'
#' @return Vector containing the estimated local FDR of the unlabeled individuals.
#' @export
Tcompute_shift <- function(W_cal, W_test, Null_cal, classification = FALSE, phat){
  h1 <- stats::density(W_cal[Null_cal])$bw
  h2 <- stats::density(W_test)$bw
  if(classification){
    h2 <- h1
  }
  f0 <- ks::kde(W_cal[Null_cal], h = h1, eval.points = W_test)
  f <- ks::kde(W_test, h = h2, eval.points = W_test)
  TN <- (1-phat)*f0$estimate/f$estimate
  TN[which(TN>1)] <- 1
  TN[which(TN<0)] <- 0
  if(classification){
    Ind <- order(W_test)
    TNW <- TN[Ind]
    for (i in (length(TNW)-20):1) {
      if(TNW[i]<TNW[i+1]){TNW[i] <- TNW[i+1]}
    }
    TN[Ind] <- TNW
  }
  return(TN)
}



#' Estimating the signal proportion by Storey method.
#'
#' @param W_cal Vector. Prediction values on the calibration data given by the model trained on the training set.
#' @param W_test Vector. Prediction values on the unlabeled data given by the model trained on the training set.
#' @param Null_cal Vector. Indices of the individuals not of interest among the calibration data.
#' @param region List. Should contain two attributes. The first attribute named 'type' takes values among "==A", "<=A", ">=B", "<=A|>=B" and ">=A&<=B". This characterizes the shape of the null region, where "==A" corresponds to classification settings. The second attribute named "v" characterizes the boundary of the null region. It should be a single number or a vector of length 2. If Y is binary then set v=-1.
#' @param tau Numeric. Threshold of the Storey estimation of the signal proportion under distribution shift.
pEstimate <- function(W_cal, W_test, Null_cal, region, tau = 0.5){
  pval <- confomalPvalue(W_cal, W_test, Null_cal, region)
  Fn <- stats::ecdf(pval)
  phat <- max(1-(1-Fn(tau))/(1-tau), 0)
  return(phat)
}



#' Computing the non-conformity score
#'
#' @param pred Vector. Prediction values.
#' @param region List. Should contain two attributes. The first attribute named 'type' takes values among "==A", "<=A", ">=B", "<=A|>=B" and ">=A&<=B". This characterizes the shape of the null region, where "==A" corresponds to classification settings. The second attribute named "v" characterizes the boundary of the null region. It should be a single number or a vector of length 2. If Y is binary then set v=-1.
ScoreCompute <- function(pred, region){
  if(region$type=="==A"|region$type=="<=A"){
    Phi <- pred
  }else if(region$type==">=B"){
    Phi <- -pred
  }else if(region$type=="<=A|>=B"){
    Phi <- pmin(pred-region$v[1], region$v[2]-pred)
  }else if(region$type==">=A&<=B"){
    Phi <- pmax(region$v[1]-pred, pred-region$v[2])
  }
  return(Phi)
}



#' Computing conformal p-values
#'
#' @param W_cal Vector. Prediction values on the calibration data given by the model trained on the training set.
#' @param W_test Vector. Prediction values on the unlabeled data given by the model trained on the training set.
#' @param Null_cal Vector. Indices of the individuals not of interest among the calibration data.
#' @param region List. Should contain two attributes. The first attribute named 'type' takes values among "==A", "<=A", ">=B", "<=A|>=B" and ">=A&<=B". This characterizes the shape of the null region, where "==A" corresponds to classification settings. The second attribute named "v" characterizes the boundary of the null region. It should be a single number or a vector of length 2. If Y is binary then set v=-1.
confomalPvalue <- function(W_cal, W_test, Null_cal, region){
  Phi_cal <- -ScoreCompute(W_cal, region)
  Phi_test <- -ScoreCompute(W_test, region)
  Phi_Null <- Phi_cal[Null_cal]
  n2 <- length(Phi_Null)
  pvalue <- sapply(Phi_test, function(t){
    (sum(Phi_Null<=t)+1)/(n2+1)
  })
  return(pvalue)
}



#' Implementation of the PASS subsampling method
#'
#' @param X_label Matrix. The covariate matrix of the labeled data.
#' @param Y_label Vector. Responses of the labeled data. Should take +-1 for binary Y.
#' @param X_test Matrix. The covariate matrix of the unlabeled data.
#' @param Null_label Vector. Indices of the individuals not of interest among the labeled data.
#' @param g Character. The criterion function used for measuring similarity. Take "RBF" or "COS" and "RBF" as default.
#' @param m Integral. The sampling budget, that is, how many unlabeled individuals should be sampled.
#' @param alpha Numeric. The pre-specified error rate level.
#' @param algo Class. A class specifying the algorithms used for fitting and prediction. The class should be well-defined with two functions. One named "fitting", which takes algorithm class, covariate matrix, responses and possibly a tuning parameter as inputs and outputs a model object. The other named "Pred", which takes algorithm class, a model object and a new covariate matrix as inputs and outputs a vector of predciton values.
#' @param cv Logical. Set to TRUE if PASSC is to be implemented. FALSE as default.
#' @param K Integral. Number of folds used for PASSC. Five as default.
#' @param lambda Numeric. The tuning parameter used for fitting the predictive model.
#' @param classification Logical. If the response is discrete then set this to TRUE. FALSE as default.
#' @param shift Logical. If there is possible distribution shift between labelled and unlabelled data then set this to TRUE. FALSE as default.
#' @param region List. Should contain two attributes. The first attribute named 'type' takes values among "==A", "<=A", ">=B", "<=A|>=B" and ">=A&<=B". This characterizes the shape of the null region, where "==A" corresponds to classification settings. The second attribute named "v" characterizes the boundary of the null region. It should be a single number or a vector of length 2. If Y is binary then set v=-1.
#' @param tau Numeric. Threshold of the Storey estimation of the signal proportion under distribution shift.
#'
#' @return A list containing a vector of sampling probabilities supported on the unlabeled data set and a vector indicating which individuals are sampled.
#' @export
PASS_sampling <- function(X_label, Y_label, X_test, Null_label, g = "RBF", m, alpha, algo, cv = FALSE, shift = FALSE, K = 5, lambda = 500, classification = FALSE, region, tau = 0.5){
  d <- dim(data.frame(X_label))[2]
  n <- dim(data.frame(X_label))[1]
  N <- dim(data.frame(X_test))[1]
  X_test_scale <- scale(X_test, center = TRUE, scale = TRUE)
  if(cv){
    model_whole <- fitting(algo, X_label, Y_label, lambda = lambda)
    W_test_whole <- Pred(algo, model_whole, X_test)
    TN_cv <- KfoldTcompute(X_label, Y_label, X_test, Null_label, W_test_whole, K, lambda, algo, classification, shift, region, tau)
    Screened_cv <- which(TN_cv<=0.5)
    S <- Scompute(g, X_test_scale[Screened_cv, ], d)
    Pi_cv <- SolvePi(TN_cv, m, S, Screened_cv, alpha)
    Pi_recover_cv <- rep(0, N)
    Pi_recover_cv[Screened_cv] <- Pi_cv
    Sampled_cv <- sample(1:N, m, replace = TRUE, Pi_recover_cv)
    return(list(prob = Pi_recover_cv, Indices = Sampled_cv))
  }else{
    TrainIndex <- sample(1:n, floor(n/2))
    data_label <- data.frame(x = X_label, y = Y_label)
    X_train <- as.matrix(data_label[TrainIndex,])[, 1:d]
    X_cal <- as.matrix(data_label[-TrainIndex,])[, 1:d]
    Y_train <- as.matrix(data_label[TrainIndex,])[, d+1]
    Y_cal <- as.matrix(data_label[-TrainIndex,])[, d+1]
    X_test <- as.matrix(data.frame(x = X_test))
    model <- fitting(algo, X_train, Y_train, lambda = lambda)
    W_cal <- Pred(algo, model, X_cal)
    W_test <- Pred(algo, model, X_test)
    Null_cal <- which(setdiff(1:n, TrainIndex)%in%Null_label)
    if(shift){
      phat <- pEstimate(W_cal, W_test, Null_cal, region)
      TN <- Tcompute_shift(W_cal, W_test, Null_cal, classification = T, phat)
    }else{
      TN <- Tcompute(W_cal, W_test, Null_cal, classification)
    }
    Screened <- which(TN<=0.5)
    S <- Scompute(g, X_test_scale[Screened, ], d)
    Pi <- SolvePi(TN, m, S, Screened, alpha)
    Pi_recover <- rep(0, N)
    Pi_recover[Screened] <- Pi
    Sampled <- sample(1:N, m, replace = TRUE, Pi_recover)
    return(list(prob = Pi_recover, Indices = Sampled))
  }
}


