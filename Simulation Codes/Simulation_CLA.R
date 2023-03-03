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


N <- 5000 ### unlabeled data size
n <- 400 ### labeled data size
alpha <- 0.1 ### pre-specified FSR level
K <- 5 ### number of folds for PASSC
n_train <- round(n/2) ### training set size for single splitting
n_cal <- n-n_train ### calibration set size for single splitting
d <- 4 ### dimension of covariates
algoarray <- c(new('NN'), new('RFc'), new('SVM')) ### algorithms used for fitting and prediction
lambda <- 500 ### hyper-parameter used for different algorithms
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
  
  
  ### rescaling the unlabeled covariates and compute the similarity matrix
  X_test_scale <- scale(X_test, center = TRUE, scale = TRUE)
  S <- Scompute(g, X_test_scale, d)
  
  
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
    Pi_PASS <- SolvePi(TN, m, S, Screened = Screened, alpha)
    Record_PASS <- CriterionCompute(Pi_PASS, Screened, Null_test, m, S, Retime = 50)
    
    info <- rbind(info, list(
      FDP = Record_PASS$FDP, SIM = Record_PASS$SIM, Method = "PASS", Algorithm = algo@name, datasetting = "CLA"
    ))
    
    
    #---PASSC---#
    if(algo@name=='NN'){
      TN_cv <- KfoldTcompute(X_rest, (Y_rest+1)/2, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
    }else{
      TN_cv <- KfoldTcompute(X_rest, Y_rest, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
    }
    Screened_cv <- which(TN_cv<=0.5)
    Pi_PASSC <- SolvePi(TN_cv, m, S, Screened_cv, alpha)
    Record_PASSC <- CriterionCompute(Pi_PASSC, Screened_cv, Null_test, m, S, Retime = 50)
    
    info <- rbind(info, list(
      FDP = Record_PASSC$FDP, SIM = Record_PASSC$SIM, Method = "PASSC", Algorithm = algo@name, datasetting = "CLA"
    ))
    
    
    #---CP---#
    pval <- confomalPvalue(W_cal, W_test, Null_cal, Value)
    rej <- sort(pval)<((1:N)/N)*alpha
    rejnum <- max(which(rej==T))
    reject <- which(pval<=sort(pval)[rejnum])
    if(rejnum==-Inf){
      Record_CP <- data.frame(FDP = 0, SIM = 0)
    }else{
      Record_CP <- CriterionCompute(rep(1/rejnum, rejnum), reject, Null_test, m, S, Retime = 50)
    }
    
    info <- rbind(info, list(
      FDP = Record_CP$FDP, SIM = Record_CP$SIM, Method = "CP", Algorithm = algo@name, datasetting = "CLA"
    ))
    
    
    #---SSD---#
    Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value)) 
    Pi_SSD <- SolvePi(TN_cv, m, S, Screened_pred, alpha, constraint = FALSE)
    Record_SSD <- CriterionCompute(Pi_SSD, Screened_pred, Null_test, m, S, Retime = 50)
    
    info <- rbind(info, list(
      FDP = Record_SSD$FDP, SIM = Record_SSD$SIM, Method = "SSD", Algorithm = algo@name, datasetting = "CLA"
    ))
    
    
    #---SS---#
    Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value))
    Pi_SS <- rep(1/length(Screened_pred), length(Screened_pred))
    Record_SS <- CriterionCompute(Pi_SS, Screened_pred, Null_test, m, S, Retime = 50)
    
    info <- rbind(info, list(
      FDP = Record_SS$FDP, SIM = Record_SS$SIM, Method = "SS", Algorithm = algo@name, datasetting = "CLA"
    ))
    
    
    #---LocF---#
    TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
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
pp <- Result%>%
  group_by(Method, Algorithm, datasetting)%>%
  dplyr::summarize(FDR = mean(FDP), SIMmean = mean(SIM), sdFDR = sd(FDP)/sqrt(ns), SIMsd = sd(SIM)/sqrt(ns))
pp


### plotting
Result$datasetting <- "Setting B"
Result$Algorithm[Result$Algorithm=="RFc"] <- "RF"
Result$Method <- factor(Result$Method, levels = c("PASS", "PASSC", "SS", "SSD", "LocF", "CP"))
Result$Algorithm <- factor(Result$Algorithm, levels = c("RF", "SVM", "NN"))

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
  facet_grid(.~Algorithm, scales = "free")

P2 <- ggplot(data = Result, aes(x = Method, y = SIM, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_continuous(name = "ES", limits = c(0.03, 0.125)) +
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
