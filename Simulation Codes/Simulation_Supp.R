library(kernlab)
library(MASS)
library(ks)
library(foreach)
library(randomForest)
library(doParallel)
require(tidyverse)
require(ggplot2)
library(ggpubr)
library(nnet)
library(caret)

source("functions_PASS.R")
source("algoclass_PASS.R")


########## Select bandwidth by cross-validation
### Revised function for CV-selected bandwidths
Tcompute_cv <- function(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = FALSE, IsCali = FALSE){
  if(h1==0){
    h1 <- density(W_cal[Null_cal], bw = 'ucv')$bw
  }
  if(h2==0){
    h2 <- density(W_cal, bw = 'ucv')$bw
  }
  if(IsSame){
    h1 <- h2
  }
  f0 <- kde(W_cal[Null_cal], h = h1, eval.points = W_test)
  f <- kde(W_cal, h = h2, eval.points = W_test)
  qhat <- 1-length(Null_cal)/length(W_cal)
  TN <- (1-qhat)*f0$estimate/f$estimate
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

N <- 5000 ### unlabeled data size
n <- 400 ### labeled data size
alpha <- 0.1 ### pre-specified FSR level
K <- 5 ### number of folds for PASSC
n_train <- round(n/2) ### training set size for single splitting
n_cal <- n-n_train ### calibration set size for single splitting
d <- 4 ### dimension of covariates
sgm <- 1 ### sigma value of the RBF criterion
algoarray <- c(new('RF'), new('SVM-R'), new('NN-R')) ### algorithms used for fitting and prediction
lambda <- 500 ### number opf trees
m <- 200 ### sampling budget
g <- 'RBF' ### similarity criterion, RBF or COS
ns <- 120 ### time of replications

### determining the region A
X <- matrix(runif((n+N)*d, -1, 1), nrow = n+N, ncol = d)
Y <- 3*X[,1]+6*X[,2]^2+7*exp(X[,3]+2*X[,4]^2)+rnorm(n+N, sd = 1)
l1 <- quantile(Y, 0.1)
l2 <- quantile(Y, 0.9)

cl <- makeCluster(8) ### number of parallel threads
registerDoParallel(cl)
ResultA <- foreach(iter = 1:ns, .combine = "rbind", .packages = c('MASS', 'glmnet', "caret", "kernlab", "randomForest", "ks", "nnet"), .errorhandling = "remove")%dopar% {

  info <- data.frame()

  ### data generation and splitting
  X <- matrix(runif((n+N)*d, -1, 1), nrow = n+N, ncol = d)
  Y <- 3*X[,1]+6*X[,2]^2+7*exp(X[,3]+2*X[,4]^2)+rnorm(n+N, sd = 1)
  data <- data.frame(x = X, y = Y)
  Value <- list(type = ">=A&<=B", v = c(l1, l2))

  datawork <- DataSplit(data, n+N, N, n_cal)
  data_train <- datawork$data_train
  data_cal <- datawork$data_cal
  data_rest <- datawork$data_rest
  data_test <- datawork$data_test

  X_train <- as.matrix(data_train[colnames(data_train)[-d-1]])
  Y_train <- as.matrix(data_train$y)
  X_cal <- as.matrix(data_cal[colnames(data_cal)[-d-1]])
  Y_cal <- as.matrix(data_cal$y)
  X_rest <- as.matrix(data_rest[colnames(data_rest)[-d-1]])
  Y_rest <- as.matrix(data_rest$y)
  X_test <- as.matrix(data_test[colnames(data_test)[-d-1]])
  Y_test <- as.matrix(data_test$y)

  Null_cal <- NullIndex(data_cal$y, Value)
  Null_rest <- NullIndex(data_rest$y, Value)
  Null_test <- NullIndex(data_test$y, Value)
  Alter_test <- setdiff(1:length(data_test$y), Null_test)

  ### rescaling the unlabeled covariates
  X_test_scale <- scale(X_test, center = TRUE, scale = TRUE)

  ### implementation via different algorithms
  for (algo in algoarray) {

    ### model fitting and prediction
    model <- fitting(algo, X_train, Y_train, lambda)
    model_whole <- fitting(algo, X_rest, Y_rest, lambda)
    W_cal <- Pred(algo, model, X_cal)
    W_test <- Pred(algo, model, X_test)
    W_test_whole <- Pred(algo, model_whole, X_test)

    ### implementation of different methods
    #---PASS---#
    TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0)
    Screened <- which(TN<=0.5)
    S <- Scompute(g, X_test_scale[Screened,], d, sgm = as.numeric(sgm))
    Pi_PASS <- SolvePi(TN, m, S, Screened = Screened, alpha)
    Record_PASS <- CriterionCompute(Pi_PASS, Screened, Null_test, m, S, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_PASS$FDP, SIM = Record_PASS$SIM, Method = "PASS", Algorithm = algo@name, datasetting = "Setting A"
    ))

    #---PASSC---#
    TN_cv <- KfoldTcompute(X_rest, Y_rest, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0)
    Screened_cv <- which(TN_cv<=0.5)
    S_cv <- Scompute(g, X_test_scale[Screened_cv,], d, sgm = as.numeric(sgm))
    Pi_PASSC <- SolvePi(TN_cv, m, S_cv, Screened_cv, alpha)
    Record_PASSC <- CriterionCompute(Pi_PASSC, Screened_cv, Null_test, m, S_cv, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_PASSC$FDP, SIM = Record_PASSC$SIM, Method = "PASSC", Algorithm = algo@name, datasetting = "Setting A"
    ))

    #---CP---#
    pval <- confomalPvalue(W_cal, W_test, Null_cal, Value)
    rej <- sort(pval)<((1:N)/N)*alpha
    rejnum <- max(which(rej==T))
    reject <- which(pval<=sort(pval)[rejnum])
    S_cp <- Scompute(g, X_test_scale[reject,], d, sgm = as.numeric(sgm))
    if(rejnum==-Inf){
      Record_CP <- data.frame(FDP = 0, SIM = NA)
    }else{
      Record_CP <- CriterionCompute(rep(1/rejnum, rejnum), reject, Null_test, m, S_cp, Retime = 50)
    }

    info <- rbind(info, list(
      FDP = Record_CP$FDP, SIM = Record_CP$SIM, Method = "CP", Algorithm = algo@name, datasetting = "Setting A"
    ))

    #---SSD---#
    Screened_pred <- setdiff(1:N, NullIndex(W_test_whole, Value))
    S_pred <- Scompute(g, X_test_scale[Screened_pred,], d, sgm = as.numeric(sgm))
    Pi_SSD <- SolvePi(TN_cv, m, S_pred, Screened_pred, alpha, constraint = FALSE)
    Record_SSD <- CriterionCompute(Pi_SSD, Screened_pred, Null_test, m, S_pred, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_SSD$FDP, SIM = Record_SSD$SIM, Method = "SSD", Algorithm = algo@name, datasetting = "Setting A"
    ))

    #---SS---#
    Screened_pred <- setdiff(1:N, NullIndex(W_test_whole, Value))
    Pi_SS <- rep(1/length(Screened_pred), length(Screened_pred))
    Record_SS <- CriterionCompute(Pi_SS, Screened_pred, Null_test, m, S_pred, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_SS$FDP, SIM = Record_SS$SIM, Method = "SS", Algorithm = algo@name, datasetting = "Setting A"
    ))

    #---LocF---#
    TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0)
    Screened <- which(TN<=0.5)
    Pi_loc <- rep(1/length(Screened), length(Screened))
    Record_loc <- CriterionCompute(Pi_loc, Screened, Null_test, m, S, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_loc$FDP, SIM = Record_loc$SIM, Method = "LocF", Algorithm = algo@name, datasetting = "Setting A"
    ))
  }

  return(info)
}
stopCluster(cl)

algoarray <- c(new('RFc'), new('SVM'), new('NN')) ### algorithms used for fitting and prediction
sigma <- 0.8 ### proportion of null individuals

cl <- makeCluster(8) ### number of parallel threads
registerDoParallel(cl)
ResultB <- foreach(iter = 1:ns, .combine = "rbind", .packages = c('MASS', 'glmnet', "caret", "kernlab", "randomForest", "ks", "nnet"), .errorhandling = "remove")%dopar% {

  info <- data.frame()

  ### data generation and splitting
  n1 <- round((n+N)*sigma)
  X1 <- mvrnorm(n1, c(rep(3, round(d/4)), rep(0, d-round(d/4))), diag(rep(1, d)))
  X2 <- mvrnorm((n+N)-n1, c(rep(0, d/2), rep(2, d/2)), diag(rep(1, d)))
  X <- rbind(X1, X2)
  data <- data.frame(x = X, y = c(rep(-1, n1), rep(1, (n+N)-n1)))
  data <- data[sample(1:dim(data)[1], dim(data)[1]),]
  Value <- list(type = "==A", v = -1)

  datawork <- DataSplit(data, n+N, N, n_cal, n)
  data_train <- datawork$data_train
  data_cal <- datawork$data_cal
  data_rest <- datawork$data_rest
  data_test <- datawork$data_test

  X_train <- as.matrix(data_train[colnames(data_train)[-d-1]])
  Y_train <- as.matrix(data_train$y)
  X_cal <- as.matrix(data_cal[colnames(data_cal)[-d-1]])
  Y_cal <- as.matrix(data_cal$y)
  X_rest <- as.matrix(data_rest[colnames(data_rest)[-d-1]])
  Y_rest <- as.matrix(data_rest$y)
  X_test <- as.matrix(data_test[colnames(data_test)[-d-1]])
  Y_test <- as.matrix(data_test$y)

  Null_cal <- NullIndex(data_cal$y, Value)
  Null_rest <- NullIndex(data_rest$y, Value)
  Null_test <- NullIndex(data_test$y, Value)
  Alter_test <- setdiff(1:length(data_test$y), Null_test)

  ### rescaling the unlabeled covariates
  X_test_scale <- scale(X_test, center = TRUE, scale = TRUE)

  ### implementation via different algorithms
  for (algo in algoarray) {

    ### model fitting and prediction
    if(algo@name=='NN'){
      model <- fitting(algo, X_train, (Y_train+1)/2, lambda)
      model_whole <- fitting(algo, X_rest, (Y_rest+1)/2, lambda)
    }else{
      model <- fitting(algo, X_train, Y_train, lambda)
      model_whole <- fitting(algo, X_rest, Y_rest, lambda)
    }
    W_cal <- Pred(algo, model, X_cal)
    W_test <- Pred(algo, model, X_test)
    W_test_whole <- Pred(algo, model_whole, X_test)
    W_test_class <- ClassPred(W_test_whole, algo)

    ### implementation of different methods
    #---PASS---#
    TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
    Screened <- which(TN<=0.5)
    S <- Scompute(g, X_test_scale[Screened,], d, sgm = sgm)
    Pi_PASS <- SolvePi(TN, m, S, Screened = Screened, alpha)
    Record_PASS <- CriterionCompute(Pi_PASS, Screened, Null_test, m, S, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_PASS$FDP, SIM = Record_PASS$SIM, Method = "PASS", Algorithm = algo@name, datasetting = "Setting B"
    ))

    #---PASSC---#
    if(algo@name=='NN'){
      TN_cv <- KfoldTcompute(X_rest, (Y_rest+1)/2, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
    }else{
      TN_cv <- KfoldTcompute(X_rest, Y_rest, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
    }
    Screened_cv <- which(TN_cv<=0.5)
    S_cv <- Scompute(g, X_test_scale[Screened_cv,], d, sgm = sgm)
    Pi_PASSC <- SolvePi(TN_cv, m, S_cv, Screened_cv, alpha)
    Record_PASSC <- CriterionCompute(Pi_PASSC, Screened_cv, Null_test, m, S_cv, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_PASSC$FDP, SIM = Record_PASSC$SIM, Method = "PASSC", Algorithm = algo@name, datasetting = "Setting B"
    ))

    #---CP---#
    pval <- confomalPvalue(W_cal, W_test, Null_cal, Value)
    rej <- sort(pval)<((1:N)/N)*alpha
    rejnum <- max(which(rej==T))
    reject <- which(pval<=sort(pval)[rejnum])
    S_cp <- Scompute(g, X_test_scale[reject,], d, sgm = sgm)
    if(rejnum==-Inf){
      Record_CP <- data.frame(FDP = 0, SIM = NA)
    }else{
      Record_CP <- CriterionCompute(rep(1/rejnum, rejnum), reject, Null_test, m, S_cp, Retime = 50)
    }

    info <- rbind(info, list(
      FDP = Record_CP$FDP, SIM = Record_CP$SIM, Method = "CP", Algorithm = algo@name, datasetting = "Setting B"
    ))

    #---SSD---#
    Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value))
    S_pred <- Scompute(g, X_test_scale[Screened_pred,], d, sgm = sgm)
    Pi_SSD <- SolvePi(TN_cv, m, S_pred, Screened_pred, alpha, constraint = FALSE)
    Record_SSD <- CriterionCompute(Pi_SSD, Screened_pred, Null_test, m, S_pred, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_SSD$FDP, SIM = Record_SSD$SIM, Method = "SSD", Algorithm = algo@name, datasetting = "Setting B"
    ))

    #---SS---#
    Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value))
    Pi_SS <- rep(1/length(Screened_pred), length(Screened_pred))
    Record_SS <- CriterionCompute(Pi_SS, Screened_pred, Null_test, m, S_pred, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_SS$FDP, SIM = Record_SS$SIM, Method = "SS", Algorithm = algo@name, datasetting = "Setting B"
    ))

    #---LocF---#
    TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
    Screened <- which(TN<=0.5)
    Pi_loc <- rep(1/length(Screened), length(Screened))
    Record_loc <- CriterionCompute(Pi_loc, Screened, Null_test, m, S, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_loc$FDP, SIM = Record_loc$SIM, Method = "LocF", Algorithm = algo@name, datasetting = "Setting B"
    ))
  }

  return(info)
}
stopCluster(cl)

ResultA <- ResultA[1:1800,]
ResultB <- ResultB[1:1800,]

ResultA$Algorithm[ResultA$Algorithm=="NN-R"] <- "NN"
ResultA$Algorithm[ResultA$Algorithm=="SVM-R"] <- "SVM"
ResultA$Method <- factor(ResultA$Method, levels = c("PASS", "PASSC", "SS", "SSD", "LocF", "CP"))
ResultA$Algorithm <- factor(ResultA$Algorithm, levels = c("RF", "SVM", "NN"))

P1 <- ggplot(data = ResultA, aes(x = Method, y = FDP, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_continuous(name = "FSR") +
  scale_x_discrete(name = "Method") +
  theme_bw() +
  geom_hline(aes(yintercept = alpha), colour = "#AA0000") +
  stat_summary(mapping = aes(group = Method),
               fun = "mean",
               geom = "point", shape = 23, size = 1.1, fill = "red",
               position = position_dodge(0.8)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Algorithm, scales = "free")

P2 <- ggplot(data = ResultA, aes(x = Method, y = SIM, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_continuous(name = "ES") +
  scale_x_discrete(name = "Method") +
  theme_bw() +
  stat_summary(mapping = aes(group = Method),
               fun = "mean",
               geom = "point", shape = 23, size = 1.1, fill = "red",
               position = position_dodge(0.8)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Algorithm, scales = "free")

PP <- ggarrange(P1, P2)
PPR <- annotate_figure(PP, top = text_grob("Scenario A", color = "black", face = "bold", size = 14))

ResultB$Algorithm[ResultB$Algorithm=="RFc"] <- "RF"
ResultB$Method <- factor(ResultB$Method, levels = c("PASS", "PASSC", "SS", "SSD", "LocF", "CP"))
ResultB$Algorithm <- factor(ResultB$Algorithm, levels = c("RF", "SVM", "NN"))

P1 <- ggplot(data = ResultB, aes(x = Method, y = FDP, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_continuous(name = "FSR") +
  scale_x_discrete(name = "Method") +
  theme_bw() +
  geom_hline(aes(yintercept = alpha), colour = "#AA0000") +
  stat_summary(mapping = aes(group = Method),
               fun = "mean",
               geom = "point", shape = 23, size = 1.1, fill = "red",
               position = position_dodge(0.8)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Algorithm, scales = "free")

P2 <- ggplot(data = ResultB, aes(x = Method, y = SIM, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_continuous(name = "ES") +
  scale_x_discrete(name = "Method") +
  theme_bw() +
  stat_summary(mapping = aes(group = Method),
               fun = "mean",
               geom = "point", shape = 23, size = 1.1, fill = "red",
               position = position_dodge(0.8)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Algorithm, scales = "free")

PP <- ggarrange(P1, P2)
PPC <- annotate_figure(PP, top = text_grob("Scenario B", color = "black", face = "bold", size = 14))

PPRC <- ggarrange(PPC, PPR, nrow = 2)




########## Different values of sigma in the RBF kernel
N <- 5000 ### unlabeled data size
n <- 400 ### labeled data size
alpha <- 0.1 ### pre-specified FSR level
K <- 5 ### number of folds for PASSC
n_train <- round(n/2) ### training set size for single splitting
n_cal <- n-n_train ### calibration set size for single splitting
d <- 4 ### dimension of covariates
sgmarray <- c('1', '2', '5', '10') ### sigma value of the RBF criterion
algoarray <- c(new('RFc')) ### algorithms used for fitting and prediction
lambda <- 500 ### hyper-parameter used for different algorithms
m <- 200 ### sampling budget
g <- 'RBF' ### similarity criterion, RBF or COS
sigma <- 0.8 ### proportion of null individuals
ns <- 110 ### time of replications

cl <- makeCluster(8) ### number of parallel threads
registerDoParallel(cl)
Result <- foreach(iter = 1:ns, .combine = "rbind", .packages = c('MASS', 'glmnet', "caret", "kernlab", "randomForest", "ks", "nnet"), .errorhandling = "remove")%dopar% {

  info <- data.frame()

  ### data generation and splitting
  n1 <- round((n+N)*sigma)
  X1 <- mvrnorm(n1, c(rep(3, round(d/4)), rep(0, d-round(d/4))), diag(rep(1, d)))
  X2 <- mvrnorm((n+N)-n1, c(rep(0, d/2), rep(2, d/2)), diag(rep(1, d)))
  X <- rbind(X1, X2)
  data <- data.frame(x = X, y = c(rep(-1, n1), rep(1, (n+N)-n1)))
  data <- data[sample(1:dim(data)[1], dim(data)[1]),]
  Value <- list(type = "==A", v = -1)

  datawork <- DataSplit(data, n+N, N, n_cal, n)
  data_train <- datawork$data_train
  data_cal <- datawork$data_cal
  data_rest <- datawork$data_rest
  data_test <- datawork$data_test

  X_train <- as.matrix(data_train[colnames(data_train)[-d-1]])
  Y_train <- as.matrix(data_train$y)
  X_cal <- as.matrix(data_cal[colnames(data_cal)[-d-1]])
  Y_cal <- as.matrix(data_cal$y)
  X_rest <- as.matrix(data_rest[colnames(data_rest)[-d-1]])
  Y_rest <- as.matrix(data_rest$y)
  X_test <- as.matrix(data_test[colnames(data_test)[-d-1]])
  Y_test <- as.matrix(data_test$y)

  Null_cal <- NullIndex(data_cal$y, Value)
  Null_rest <- NullIndex(data_rest$y, Value)
  Null_test <- NullIndex(data_test$y, Value)
  Alter_test <- setdiff(1:length(data_test$y), Null_test)

  ### rescaling the unlabeled covariates
  X_test_scale <- scale(X_test, center = TRUE, scale = TRUE)

  for (sgm in sgmarray) {

    ### implementation via different algorithms
    for (algo in algoarray) {

      ### model fitting and prediction
      if(algo@name=='NN'){
        model <- fitting(algo, X_train, (Y_train+1)/2, lambda)
        model_whole <- fitting(algo, X_rest, (Y_rest+1)/2, lambda)
      }else{
        model <- fitting(algo, X_train, Y_train, lambda)
        model_whole <- fitting(algo, X_rest, Y_rest, lambda)
      }
      W_cal <- Pred(algo, model, X_cal)
      W_test <- Pred(algo, model, X_test)
      W_test_whole <- Pred(algo, model_whole, X_test)
      W_test_class <- ClassPred(W_test_whole, algo)

      ### implementation of different methods
      #---PASS---#
      TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
      Screened <- which(TN<=0.5)
      S <- Scompute(g, X_test_scale[Screened,], d, sgm = as.numeric(sgm))
      Pi_PASS <- SolvePi(TN, m, S, Screened = Screened, alpha)
      Record_PASS <- CriterionCompute(Pi_PASS, Screened, Null_test, m, S, Retime = 50)

      info <- rbind(info, list(
        FDP = Record_PASS$FDP, SIM = Record_PASS$SIM, Method = "PASS", Algorithm = algo@name, datasetting = "CLA", sgm = sgm
      ))

      #---PASSC---#
      if(algo@name=='NN'){
        TN_cv <- KfoldTcompute(X_rest, (Y_rest+1)/2, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
      }else{
        TN_cv <- KfoldTcompute(X_rest, Y_rest, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
      }
      Screened_cv <- which(TN_cv<=0.5)
      S_cv <- Scompute(g, X_test_scale[Screened_cv,], d, sgm = as.numeric(sgm))
      Pi_PASSC <- SolvePi(TN_cv, m, S_cv, Screened_cv, alpha)
      Record_PASSC <- CriterionCompute(Pi_PASSC, Screened_cv, Null_test, m, S_cv, Retime = 50)

      info <- rbind(info, list(
        FDP = Record_PASSC$FDP, SIM = Record_PASSC$SIM, Method = "PASSC", Algorithm = algo@name, datasetting = "CLA", sgm = sgm
      ))

      #---CP---#
      pval <- confomalPvalue(W_cal, W_test, Null_cal, Value)
      rej <- sort(pval)<((1:N)/N)*alpha
      rejnum <- max(which(rej==T))
      reject <- which(pval<=sort(pval)[rejnum])
      S_cp <- Scompute(g, X_test_scale[reject,], d, sgm = as.numeric(sgm))
      if(rejnum==-Inf){
        Record_CP <- data.frame(FDP = 0, SIM = NA)
      }else{
        Record_CP <- CriterionCompute(rep(1/rejnum, rejnum), reject, Null_test, m, S_cp, Retime = 50)
      }

      info <- rbind(info, list(
        FDP = Record_CP$FDP, SIM = Record_CP$SIM, Method = "CP", Algorithm = algo@name, datasetting = "CLA", sgm = sgm
      ))

      #---SSD---#
      Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value))
      S_pred <- Scompute(g, X_test_scale[Screened_pred,], d, sgm = as.numeric(sgm))
      Pi_SSD <- SolvePi(TN_cv, m, S_pred, Screened_pred, alpha, constraint = FALSE)
      Record_SSD <- CriterionCompute(Pi_SSD, Screened_pred, Null_test, m, S_pred, Retime = 50)

      info <- rbind(info, list(
        FDP = Record_SSD$FDP, SIM = Record_SSD$SIM, Method = "SSD", Algorithm = algo@name, datasetting = "CLA", sgm = sgm
      ))

      #---SS---#
      Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value))
      Pi_SS <- rep(1/length(Screened_pred), length(Screened_pred))
      Record_SS <- CriterionCompute(Pi_SS, Screened_pred, Null_test, m, S_pred, Retime = 50)

      info <- rbind(info, list(
        FDP = Record_SS$FDP, SIM = Record_SS$SIM, Method = "SS", Algorithm = algo@name, datasetting = "CLA", sgm = sgm
      ))

      #---LocF---#
      TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
      Screened <- which(TN<=0.5)
      Pi_loc <- rep(1/length(Screened), length(Screened))
      Record_loc <- CriterionCompute(Pi_loc, Screened, Null_test, m, S, Retime = 50)

      info <- rbind(info, list(
        FDP = Record_loc$FDP, SIM = Record_loc$SIM, Method = "LocF", Algorithm = algo@name, datasetting = "CLA", sgm = sgm
      ))
    }
  }

  return(info)
}
stopCluster(cl)

### summary
Result <- Result[1:3000,]
pp <- Result%>%
  group_by(Method, Algorithm, datasetting, sgm)%>%
  dplyr::summarize(FDR = mean(FDP), SIMmean = mean(SIM), sdFDR = sd(FDP)/sqrt(100), SIMsd = sd(SIM)/sqrt(100))
pp

### plotting
Result$datasetting <- "Setting B"
Result$Algorithm[Result$Algorithm=="RFc"] <- "RF"
Result$Method <- factor(Result$Method, levels = c("PASS", "PASSC", "SS", "SSD", "LocF", "CP"))
Result$Algorithm <- factor(Result$Algorithm, levels = c("RF", "SVM", "NN"))
Result$sgm[Result$sgm=='1'] <- 'aaa'
Result$sgm[Result$sgm=='2'] <- 'bbb'
Result$sgm[Result$sgm=='5'] <- 'ccc'
Result$sgm[Result$sgm=='10'] <- 'ddd'
Result$sgm <- factor(Result$sgm, levels = c('aaa', 'bbb', 'ccc', 'ddd'))
sgmname <- c(aaa = "sigma = 1", bbb = "sigma = 2", ccc = "sigma = 5", ddd ="sigma = 10")

P1 <- ggplot(data = Result, aes(x = Method, y = FDP, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_continuous(name = "FSR") +
  scale_x_discrete(name = "Method") +
  theme_bw() +
  geom_hline(aes(yintercept = alpha), colour = "#AA0000") +
  stat_summary(mapping = aes(group = Method),
               fun = "mean",
               geom = "point", shape = 23, size = 1.1, fill = "red",
               position = position_dodge(0.8)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~sgm, scales = "free", labeller = labeller(sgm = sgmname))

P2 <- ggplot(data = Result, aes(x = Method, y = SIM, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_continuous(name = "ES") +
  scale_x_discrete(name = "Method") +
  theme_bw() +
  stat_summary(mapping = aes(group = Method),
               fun = "mean",
               geom = "point", shape = 23, size = 1.1, fill = "red",
               position = position_dodge(0.8)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~sgm, scales = "free", labeller = labeller(sgm = sgmname))

PP <- ggarrange(P1, P2, nrow = 2)
PPC <- annotate_figure(PP, top = text_grob("Scenario B", color = "black", face = "bold", size = 14))
PPC




########## Implementation under distribution shift
### Use default fitting for NN
setClass("NN", slots = list(name = "character"),
         prototype = list(name = "NN"))
setMethod("fitting", "NN", function(obj, X, Y, lambda){
  datawork <- data.frame(X, y = Y)
  nnet(y~., data = datawork, size = 10, decay = 5e-4, maxit = 2000)
})
setMethod("Pred", "NN", function(obj, model, X_test, class = F){
  if(class==T){
    predict(model, as.data.frame(x = X_test), type = 'class')
  }else{
    predict(model, as.data.frame(x = X_test))
  }
})

N <- 5000 ### unlabeled data size
n <- 800 ### labeled data size
alpha <- 0.1 ### pre-specified FSR level
K <- 5 ### number of folds for PASSC
n_train <- round(n/2) ### training set size for single splitting
n_cal <- n-n_train ### calibration set size for single splitting
d <- 4 ### dimension of covariates
sgm <- 1 ### sigma value of the RBF criterion
algoarray <- c(new('RFc'), new('SVM'), new('NN')) ### algorithms used for fitting and prediction
lambda <- 500 ### number of trees
m <- 150 ### sampling budget
g <- 'RBF' ### similarity criterion, RBF or COS
sigma <- 0.8 ### proportion of labelled null individuals
alpUn <- 0.05 ### proportion of unlabelled non-null individuals
tau <- 0.2 ### threshold of the Storey estimation

ns <- 250 ### time of replications
cl <- makeCluster(8) ### number of parallel threads
registerDoParallel(cl)
Result <- foreach(iter = 1:ns, .combine = "rbind", .packages = c('MASS', 'glmnet', "caret", "kernlab", "randomForest", "ks", "nnet", "evd"), .errorhandling = "remove")%dopar% {

  info <- data.frame()

  ### data generation and splitting
  X1 <- rbind(mvrnorm(round(n*(1-sigma)), c(rep(0, d/2), 2, 2), diag(rep(1, d))), mvrnorm(round(N*alpUn), c(rep(0, d/2), 2.5, 1.5), diag(rep(1, d))))
  X2 <- mvrnorm(n*sigma+N*(1-alpUn), c(rep(3, round(d/4)), rep(0, d-round(d/4))), diag(rep(1, d)))
  X <- rbind(X1, X2)
  data <- data.frame(x = X, y = c(rep(1, n*(1-sigma)+N*alpUn), rep(-1, n*(sigma)+N*(1-alpUn))))
  Value <- list(type = "==A", v = -1)

  data_label <- data[c(1:((1-sigma)*n),(n*(1-sigma)+N*alpUn+1):(n+N*alpUn)),]
  data_test <- data[-c(1:((1-sigma)*n),(n*(1-sigma)+N*alpUn+1):(n+N*alpUn)),]
  data_label <- data_label[sample(1:dim(data_label)[1], dim(data_label)[1]),]
  data_train <- data_label[1:n_train,]
  data_cal <- data_label[(n_train+1):n,]

  X_train <- as.matrix(data_train[colnames(data_train)[-d-1]])
  Y_train <- as.matrix(data_train$y)
  X_cal <- as.matrix(data_cal[colnames(data_cal)[-d-1]])
  Y_cal <- as.matrix(data_cal$y)
  X_label <- as.matrix(data_label[colnames(data_label)[-d-1]])
  Y_label <- as.matrix(data_label$y)
  X_test <- as.matrix(data_test[colnames(data_test)[-d-1]])
  Y_test <- as.matrix(data_test$y)

  Null_cal <- NullIndex(data_cal$y, Value)
  Null_label <- NullIndex(data_label$y, Value)
  Null_test <- NullIndex(data_test$y, Value)
  Alter_test <- setdiff(1:length(data_test$y), Null_test)

  ### rescaling the unlabeled covariates
  X_test_scale <- scale(X_test, center = TRUE, scale = TRUE)

  ### implementation via different algorithms
  for (algo in algoarray) {

    ### model fitting and prediction
    if(algo@name=='NN'){
      model <- fitting(algo, X_train, (Y_train+1)/2, lambda)
      model_whole <- fitting(algo, X_label, (Y_label+1)/2, lambda)
    }else{
      model <- fitting(algo, X_train, Y_train, lambda)
      model_whole <- fitting(algo, X_label, Y_label, lambda)
    }
    W_cal <- Pred(algo, model, X_cal)
    W_test <- Pred(algo, model, X_test)
    W_test_whole <- Pred(algo, model_whole, X_test)
    W_test_class <- ClassPred(W_test_whole, algo)

    pval <- Pvalue(W_cal, W_test, Null_cal, classification = T)
    qhat <- max(1-sum(pval>tau)/length(pval)/(1-tau), 0)

    ### implementation of different methods
    #---PASS---#
    TN <- rTcompute(W_cal, W_test, Null_cal, classification = T)
    Screened <- which(TN<=0.5)
    S <- Scompute(g, X_test_scale[Screened,], d, sgm = as.numeric(sgm))
    Pi_PASS <- SolvePi(TN, m, S, Screened = Screened, alpha)
    Record_PASS <- CriterionCompute(Pi_PASS, Screened, Null_test, m, S, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_PASS$FDP, SIM = Record_PASS$SIM, Method = "rPASS", Algorithm = algo@name, datasetting = "CLA"
    ))

    #---CP---#
    pval <- confomalPvalue(W_cal, W_test, Null_cal, Value)
    rej <- sort(pval)<((1:N)/N)*alpha
    rejnum <- max(which(rej==T))
    reject <- which(pval<=sort(pval)[rejnum])
    if(rejnum==-Inf){
      Record_CP <- data.frame(FDP = 0, SIM = NA)
    }else{
      S_cp <- Scompute(g, X_test_scale[reject,], d, sgm = as.numeric(sgm))
      Record_CP <- CriterionCompute(rep(1/rejnum, rejnum), reject, Null_test, m, S_cp, Retime = 50)
    }

    info <- rbind(info, list(
      FDP = Record_CP$FDP, SIM = Record_CP$SIM, Method = "CP", Algorithm = algo@name, datasetting = "CLA"
    ))

    #---SSD---#
    Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value))
    S_pred <- Scompute(g, X_test_scale[Screened_pred,], d, sgm = as.numeric(sgm))
    Pi_SSD <- SolvePi(TN_cv, m, S_pred, Screened_pred, alpha, constraint = FALSE)
    Record_SSD <- CriterionCompute(Pi_SSD, Screened_pred, Null_test, m, S_pred, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_SSD$FDP, SIM = Record_SSD$SIM, Method = "SSD", Algorithm = algo@name, datasetting = "CLA"
    ))

    #---SS---#
    Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value))
    Pi_SS <- rep(1/length(Screened_pred), length(Screened_pred))
    Record_SS <- CriterionCompute(Pi_SS, Screened_pred, Null_test, m, S_pred, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_SS$FDP, SIM = Record_SS$SIM, Method = "SS", Algorithm = algo@name, datasetting = "CLA"
    ))

    #---LocF---#
    Screened <- which(TN<=0.5)
    Pi_loc <- rep(1/length(Screened), length(Screened))
    Record_loc <- CriterionCompute(Pi_loc, Screened, Null_test, m, S, Retime = 50)

    info <- rbind(info, list(
      FDP = Record_loc$FDP, SIM = Record_loc$SIM, Method = "LocF", Algorithm = algo@name, datasetting = "CLA"
    ))
  }

  return(info)
}
stopCluster(cl)

### summary
Result <- Result[1:1800,]
pp <- Result%>%
  group_by(Method, Algorithm, datasetting)%>%
  dplyr::summarize(FDR = mean(FDP), SIMmean = mean(SIM), sdFDR = sd(FDP)/sqrt(100), SIMsd = sd(SIM)/sqrt(100))
pp

### plotting
Result$datasetting <- "Setting B"
Result$Algorithm[Result$Algorithm=="RFc"] <- "RF"
Result$Method <- factor(Result$Method, levels = c("PASS", "SS", "SSD", "LocF", "CP"))
Result$Algorithm <- factor(Result$Algorithm, levels = c("RF", "SVM", "NN", "GLM"))

P1 <- ggplot(data = Result, aes(x = Method, y = FDP, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_continuous(name = "FSR", limits = c(0,0.5)) +
  scale_x_discrete(name = "Method") +
  theme_bw() +
  geom_hline(aes(yintercept = alpha), colour = "#AA0000") +
  stat_summary(mapping = aes(group = Method),
               fun = "mean",
               geom = "point", shape = 23, size = 1.1, fill = "red",
               position = position_dodge(0.8)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Algorithm, scales = "free")

P2 <- ggplot(data = Result, aes(x = Method, y = SIM, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_continuous(name = "ES") +
  scale_x_discrete(name = "Method") +
  theme_bw() +
  stat_summary(mapping = aes(group = Method),
               fun = "mean",
               geom = "point", shape = 23, size = 1.1, fill = "red",
               position = position_dodge(0.8)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Algorithm, scales = "free")

PP <- ggarrange(P1, P2)
PPC <- annotate_figure(PP, top = text_grob("Scenario B", color = "black", face = "bold", size = 14))
PPC



########## Computation time
N <- 5000 ### unlabeled data size
n <- 400 ### labeled data size
alpha <- 0.1 ### pre-specified FSR level
K <- 5 ### number of folds for PASSC
n_train <- round(n/2) ### training set size for single splitting
n_cal <- n-n_train ### calibration set size for single splitting
d <- 4 ### dimension of covariates
sgm <- 1 ### sigma value of the RBF criterion
algoarray <- c(new('RFc'), new('SVM'), new('NN')) ### algorithms used for fitting and prediction
lambda <- 500 ### number of trees
m <- 200 ### sampling budget
g <- 'RBF' ### similarity criterion, RBF or COS
sigma <- 0.8 ### proportion of null individuals

ns <- 100 ### time of replications
cl <- makeCluster(8) ### number of parallel threads
registerDoParallel(cl)
Result <- foreach(iter = 1:ns, .combine = "rbind", .packages = c('MASS', 'glmnet', "caret", "kernlab", "randomForest", "ks", "nnet"), .errorhandling = "remove")%dopar% {

  info <- data.frame()


  ### data generation and splitting
  n1 <- round((n+N)*sigma)
  X1 <- mvrnorm(n1, c(rep(3, round(d/4)), rep(0, d-round(d/4))), diag(rep(1, d)))
  X2 <- mvrnorm((n+N)-n1, c(rep(0, d/2), rep(2, d/2)), diag(rep(1, d)))
  X <- rbind(X1, X2)
  data <- data.frame(x = X, y = c(rep(-1, n1), rep(1, (n+N)-n1)))
  data <- data[sample(1:dim(data)[1], dim(data)[1]),]
  Value <- list(type = "==A", v = -1)

  datawork <- DataSplit(data, n+N, N, n_cal, n)
  data_train <- datawork$data_train
  data_cal <- datawork$data_cal
  data_rest <- datawork$data_rest
  data_test <- datawork$data_test

  X_train <- as.matrix(data_train[colnames(data_train)[-d-1]])
  Y_train <- as.matrix(data_train$y)
  X_cal <- as.matrix(data_cal[colnames(data_cal)[-d-1]])
  Y_cal <- as.matrix(data_cal$y)
  X_rest <- as.matrix(data_rest[colnames(data_rest)[-d-1]])
  Y_rest <- as.matrix(data_rest$y)
  X_test <- as.matrix(data_test[colnames(data_test)[-d-1]])
  Y_test <- as.matrix(data_test$y)

  Null_cal <- NullIndex(data_cal$y, Value)
  Null_rest <- NullIndex(data_rest$y, Value)
  Null_test <- NullIndex(data_test$y, Value)
  Alter_test <- setdiff(1:length(data_test$y), Null_test)


  ### rescaling the unlabeled covariates
  X_test_scale <- scale(X_test, center = TRUE, scale = TRUE)

  ### implementation via different algorithms
  for (algo in algoarray) {

    ### model fitting and prediction

    if(algo@name=='NN'){
      t_train_1 <- Sys.time()
      model <- fitting(algo, X_train, (Y_train+1)/2, lambda)
      t_train_2 <- Sys.time()
      model_whole <- fitting(algo, X_rest, (Y_rest+1)/2, lambda)
    }else{
      t_train_1 <- Sys.time()
      model <- fitting(algo, X_train, Y_train, lambda)
      t_train_2 <- Sys.time()
      model_whole <- fitting(algo, X_rest, Y_rest, lambda)
    }
    t_train_3 <- Sys.time()
    t_pred_1 <- t_train_3
    W_cal <- Pred(algo, model, X_cal)
    t_pred_2<-Sys.time()
    W_test <- Pred(algo, model, X_test)
    t_pred_3 <- Sys.time()
    W_test_whole <- Pred(algo, model_whole, X_test)
    t_pred_4 <- Sys.time()
    W_test_class <- ClassPred(W_test_whole, algo)
    t_pred_5 <- Sys.time()

    ### implementation of different methods
    #---PASS---#
    t_pass_1 <- Sys.time()
    TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
    t_pass_2 <- Sys.time()

    Screened <- which(TN<=0.5)

    t_pass_3 <- Sys.time()
    S <- Scompute(g, X_test_scale[Screened,], d, sgm = as.numeric(sgm))
    t_pass_4 <- Sys.time()

    Pi_PASS <- SolvePi(TN, m, S, Screened = Screened, alpha)
    t_pass_5 <- Sys.time()

    Record_PASS <- CriterionCompute(Pi_PASS, Screened, Null_test, m, S, Retime = 50)

    time_pass_estloc <- (t_train_2-t_train_1)+(t_pred_2-t_pred_1)+(t_pred_3-t_pred_2)+(t_pass_2-t_pass_1)+(t_pass_3-t_pass_2)
    time_pass_S <- t_pass_4-t_pass_3
    time_pass_opt <- t_pass_5-t_pass_4

    info <- rbind(info, list(
      FDP = Record_PASS$FDP, SIM = Record_PASS$SIM, Method = "PASS", Algorithm = algo@name, datasetting = "CLA",
      time_estloc=time_pass_estloc,time_S=time_pass_S,time_opt=time_pass_opt
    ))

    #---PASSC---#
    t_passc_1<-Sys.time()
    if(algo@name=='NN'){
      TN_cv <- KfoldTcompute(X_rest, (Y_rest+1)/2, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
    }else{
      TN_cv <- KfoldTcompute(X_rest, Y_rest, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
    }
    t_passc_2 <- Sys.time()

    Screened_cv <- which(TN_cv<=0.5)
    t_passc_3 <- Sys.time()
    S_cv <- Scompute(g, X_test_scale[Screened_cv,], d, sgm = as.numeric(sgm))
    t_passc_4 <- Sys.time()
    Pi_PASSC <- SolvePi(TN_cv, m, S_cv, Screened_cv, alpha)
    t_passc_5 <- Sys.time()
    Record_PASSC <- CriterionCompute(Pi_PASSC, Screened_cv, Null_test, m, S_cv, Retime = 50)

    time_passc_estloc <- (t_train_3-t_train_2)+(t_pred_4-t_pred_3)+(t_passc_2-t_passc_1)+(t_passc_3-t_passc_2)
    time_passc_S <- t_passc_4-t_passc_3
    time_passc_opt <- t_passc_5-t_passc_4

    info <- rbind(info, list(
      FDP = Record_PASSC$FDP, SIM = Record_PASSC$SIM, Method = "PASSC", Algorithm = algo@name, datasetting = "CLA",
      time_estloc=time_passc_estloc,time_S=time_passc_S,time_opt=time_passc_opt

    ))

    #---CP---#
    t_cp_1 <- Sys.time()
    pval <- confomalPvalue(W_cal, W_test, Null_cal, Value)
    rej <- sort(pval)<((1:N)/N)*alpha
    rejnum <- max(which(rej==T))
    reject <- which(pval<=sort(pval)[rejnum])
    t_cp_2 <- Sys.time()

    S_cp <- Scompute(g, X_test_scale[reject,], d, sgm = as.numeric(sgm))
    if(rejnum==-Inf){
      Record_CP <- data.frame(FDP = 0, SIM = NA)
    }else{
      Record_CP <- CriterionCompute(rep(1/rejnum, rejnum), reject, Null_test, m, S_cp, Retime = 50)
    }

    ### be careful!!! the optimization time of CP is the whole time of CP
    time_cp_whole <- (t_train_2-t_train_1)+(t_pred_2-t_pred_1)+(t_pred_3-t_pred_2)+(t_cp_2-t_cp_1)

    info <- rbind(info, list(
      FDP = Record_CP$FDP, SIM = Record_CP$SIM, Method = "CP", Algorithm = algo@name, datasetting = "CLA",
      time_estloc=0,time_S=0,time_opt=time_cp_whole

    ))

    #---SSD---#
    t_ssd_1 <- Sys.time()
    Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value))
    t_ssd_2 <- Sys.time()

    S_pred <- Scompute(g, X_test_scale[Screened_pred,], d, sgm = as.numeric(sgm))
    t_ssd_3 <- Sys.time()

    Pi_SSD <- SolvePi(TN_cv, m, S_pred, Screened_pred, alpha, constraint = FALSE)
    t_ssd_4 <- Sys.time()

    Record_SSD <- CriterionCompute(Pi_SSD, Screened_pred, Null_test, m, S_pred, Retime = 50)

    time_ssd_estloc <- (t_train_3-t_train_2)+(t_pred_4-t_pred_3)+(t_pred_5-t_pred_4)+(t_passc_2-t_passc_1)+(t_ssd_2-t_ssd_1)
    time_ssd_S <- t_ssd_3-t_ssd_2
    time_ssd_opt <- t_ssd_4-t_ssd_3

    info <- rbind(info, list(
      FDP = Record_SSD$FDP, SIM = Record_SSD$SIM, Method = "SSD", Algorithm = algo@name, datasetting = "CLA",
      time_estloc=time_ssd_estloc,time_S=time_ssd_S,time_opt=time_ssd_opt

    ))

    #---SS---#
    t_ss_1 <- Sys.time()

    Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value))
    t_ss_2 <- Sys.time()

    Pi_SS <- rep(1/length(Screened_pred), length(Screened_pred))
    t_ss_3 <- Sys.time()

    Record_SS <- CriterionCompute(Pi_SS, Screened_pred, Null_test, m, S_pred, Retime = 50)

    time_ss_whole <- (t_train_3-t_train_2)+(t_pred_4-t_pred_3)+(t_pred_5-t_pred_4)+(t_ss_3-t_ss_1)

    info <- rbind(info, list(
      FDP = Record_SS$FDP, SIM = Record_SS$SIM, Method = "SS", Algorithm = algo@name, datasetting = "CLA",
      time_estloc=0,time_S=0,time_opt=time_ss_whole
    ))

    #---LocF---#
    t_locf_1=Sys.time()

    TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
    t_locf_2 <- Sys.time()

    Screened <- which(TN<=0.5)
    t_locf_3 <- Sys.time()

    Pi_loc <- rep(1/length(Screened), length(Screened))
    t_locf_4 <- Sys.time()

    Record_loc <- CriterionCompute(Pi_loc, Screened, Null_test, m, S, Retime = 50)

    time_locf_whole <- (t_train_2-t_train_1)+(t_pred_2-t_pred_1)+(t_pred_3-t_pred_2)+t_locf_4-t_locf_1

    info <- rbind(info, list(
      FDP = Record_loc$FDP, SIM = Record_loc$SIM, Method = "LocF", Algorithm = algo@name, datasetting = "CLA",
      time_estloc=0,time_S=0,time_opt=time_locf_whole

    ))
  }

  return(info)
}
stopCluster(cl)

### summary
pp <- Result%>%
  group_by(Method, Algorithm, datasetting)%>%
  dplyr::summarize(Time_estlfdr=mean(time_estloc),Time_S=mean(time_S),Time_opt=mean(time_opt),Time_whole=mean(time_estloc+time_S+time_opt),FDR = mean(FDP), SIMmean = mean(SIM))
pp
