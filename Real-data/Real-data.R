library(kernlab)
library(MASS)
library(ks)
library(randomForest)
library(doParallel)
require(tidyverse)
require(ggplot2)
library(ggpubr)
library(caret)
library(nnet)
library(ROSE)
library(ranger)
library(RColorBrewer)


source("functions_PASS.R")
source("algoclass_PASS.R")

### loading and pre-processing the data
data0 <- read.table("adult.csv", header = T, sep = ",")
data0[data0=="?"] <- NA
data1 <- na.omit(data0)

USindex <- which(data1$native.country=="United-States") ### only consider american individuals
data1 <- select(data1[USindex,], select = -c(native.country))

varcontinue <- c("age", "fnlwgt", "education.num", "capital.gain", "capital.loss", "hours.per.week")  
colname <- colnames(data1)
y <- as.numeric(data1$income==factor(">50K", levels = c("<=50K", ">50K")))
y <- 2*y-1
data1 <- cbind(lapply(data1[, varcontinue], function(x){as.numeric(as.character(x))}),
               as.data.frame(lapply(data1[, setdiff(colname, varcontinue)], function(x){factor(x)})))

dummy <- dummyVars(" ~ .", data = data1[, -length(colname)]) ### one-hot encoding
data <- data.frame(predict(dummy, newdata = data1)) 
data$y <- y

Number <- dim(data)[1] ### sample size
N <- 5000 ### unlabeled sample size
n <- 1000 ### labeled sample size
alpha <- 0.2 ### pre-specified FSR level
K <- 5 ### number of folds for PASSC
n_train <- round(n/2) ### training set size for single splitting
n_cal <- n-n_train ### calibration set size for single splitting
d <- dim(data)[2]-1 ### dimension of covariates
algo <- new("RFc") ### algorithm used for fitting and prediction
lambda <- 500 ### number of trees
m <- 200 ### sampling budget
g <- "RBF" ### similarity criterion, RBF or COS
sigma <- 1-sum(data$y==1)/Number ### proportion of null individuals

### Specify indices of individuals of interest
Value <- list(type = "==A", v = -1)
Null <- NullIndex(data$y, Value)
Alter <- setdiff(1:Number, Null)


ns <- 250
cl <- makeCluster(8) ### number of parallel threads
registerDoParallel(cl)
Result <- foreach(iter = 1:ns, .combine = "rbind", .packages = c('MASS', 'glmnet', "caret", "kernlab", "randomForest", "ks", "nnet", "dplyr"), .errorhandling = "remove")%dopar% {
  info <- data.frame()
  
  
  ### sampling from the whole big data set (data generation)
  IndexSample <- c(sample(Null, round((n+N)*sigma), replace = FALSE), sample(Alter, n+N-round((n+N)*sigma), replace = FALSE))
  newdata <- data[sample(IndexSample, n+N, replace = FALSE),]
  
  datawork <- DataSplit(newdata, n+N, N, n_cal, n)
  data_train <- datawork$data_train
  data_cal <- datawork$data_cal
  data_label <- datawork$data_rest
  data_test <- datawork$data_test
  ClassData <- data1[rownames(data_test),]
  
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
  X_test_scale[, 7:13] <- X_test_scale[, 7:13]/7
  X_test_scale[, 14:29] <- X_test_scale[, 14:29]/16
  X_test_scale[, 30:36] <- X_test_scale[, 30:36]/7
  X_test_scale[, 37:50] <- X_test_scale[, 37:50]/14
  X_test_scale[, 51:56] <- X_test_scale[, 51:56]/6
  X_test_scale[, 57:61] <- X_test_scale[, 57:61]/5
  X_test_scale[, 62:63] <- X_test_scale[, 62:63]/2
  
  
  ### model fitting and prediction
  model <- fitting(algo, X_train, Y_train, lambda)
  model_whole <- fitting(algo, X_label, Y_label, lambda)
  W_cal <- Pred(algo, model, X_cal)
  W_test <- Pred(algo, model, X_test)
  W_test_whole <- Pred(algo, model_whole, X_test)
  W_test_class <- ClassPred(W_test_whole, algo)
  
  
  #---PASS---# 
  TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
  Screened <- which(TN<=0.5)
  S <- Scompute(g, X_test_scale[Screened,], d, sgm = 1)
  Pi_PASS <- SolvePi(TN, m, S, Screened = Screened, alpha)
  
  ### Counting and calculating proportions
  WorkTable <- data.frame()
  for (i in 1:100) {
    S_Pi <- sample(Screened, m, replace = TRUE, Pi_PASS)
    S_Pi <- S_Pi[is.element(S_Pi, Alter_test)]
    m0 <- length(S_Pi)
    Male <- sum(ClassData$sex[S_Pi]=='Male')/m0
    Female <- sum(ClassData$sex[S_Pi]=='Female')/m0
    White <- sum(ClassData$race[S_Pi]=='White')/m0
    Black <- sum(ClassData$race[S_Pi]=='Black')/m0
    Other <- sum(ClassData$race[S_Pi]=='Other')/m0
    AIE <- sum(ClassData$race[S_Pi]=='Amer-Indian-Eskimo')/m0
    API <- sum(ClassData$race[S_Pi]=='Asian-Pac-Islander')/m0
    WorkTable <- rbind(WorkTable, list(Male = Male, Female = Female, W = White, B = Black, O = Other, AIE = AIE, API = API))
  }
  pp <- WorkTable%>%
    group_by()%>%
    dplyr::summarize(Ma = mean(Male), Fe = mean(Female), Wh = mean(W), Bl = mean(B), Ot = mean(O), AI = mean(AIE), AP = mean(API))
  info <- rbind(info, list(Male = pp$Ma, Female = pp$Fe, W = pp$Wh, B = pp$Bl, O = pp$Ot, AIE = pp$AI, API = pp$AP, method = 'PASS'))
  
  #---PASSC---#
  TN_cv <- KfoldTcompute(X_label, Y_label, X_test, Null_label, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
  Screened_cv <- which(TN_cv<=0.5)
  S_cv <- Scompute(g, X_test_scale[Screened_cv,], d, sgm = 1)
  Pi_PASSC <- SolvePi(TN_cv, m, S_cv, Screened_cv, alpha)
  
  ### Counting and calculating proportions
  WorkTable <- data.frame()
  for (i in 1:100) {
    S_Pi <- sample(Screened_cv, m, replace = TRUE, Pi_PASSC)
    S_Pi <- S_Pi[is.element(S_Pi, Alter_test)]
    m0 <- length(S_Pi)
    Male <- sum(ClassData$sex[S_Pi]=='Male')/m0
    Female <- sum(ClassData$sex[S_Pi]=='Female')/m0
    White <- sum(ClassData$race[S_Pi]=='White')/m0
    Black <- sum(ClassData$race[S_Pi]=='Black')/m0
    Other <- sum(ClassData$race[S_Pi]=='Other')/m0
    AIE <- sum(ClassData$race[S_Pi]=='Amer-Indian-Eskimo')/m0
    API <- sum(ClassData$race[S_Pi]=='Asian-Pac-Islander')/m0
    WorkTable <- rbind(WorkTable, list(Male = Male, Female = Female, W = White, B = Black, O = Other, AIE = AIE, API = API))
  }
  pp <- WorkTable%>%
    group_by()%>%
    dplyr::summarize(Ma = mean(Male), Fe = mean(Female), Wh = mean(W), Bl = mean(B), Ot = mean(O), AI = mean(AIE), AP = mean(API))
  info <- rbind(info, list(Male = pp$Ma, Female = pp$Fe, W = pp$Wh, B = pp$Bl, O = pp$Ot, AIE = pp$AI, API = pp$AP, method = 'PASSC'))
  
  #---CP---#
  pval <- confomalPvalue(W_cal, W_test, Null_cal, Value)
  rej <- sort(pval)<((1:N)/N)*alpha
  rejnum <- max(which(rej==T))
  reject <- which(pval<=sort(pval)[rejnum])
  Pi_CP <- rep(1/rejnum, rejnum)
  
  ### Counting and calculating proportions
  WorkTable <- data.frame()
  for (i in 1:100) {
    S_Pi <- sample(reject, m, replace = TRUE, Pi_CP)
    S_Pi <- S_Pi[is.element(S_Pi, Alter_test)]
    m0 <- length(S_Pi)
    Male <- sum(ClassData$sex[S_Pi]=='Male')/m0
    Female <- sum(ClassData$sex[S_Pi]=='Female')/m0
    White <- sum(ClassData$race[S_Pi]=='White')/m0
    Black <- sum(ClassData$race[S_Pi]=='Black')/m0
    Other <- sum(ClassData$race[S_Pi]=='Other')/m0
    AIE <- sum(ClassData$race[S_Pi]=='Amer-Indian-Eskimo')/m0
    API <- sum(ClassData$race[S_Pi]=='Asian-Pac-Islander')/m0
    WorkTable <- rbind(WorkTable, list(Male = Male, Female = Female, W = White, B = Black, O = Other, AIE = AIE, API = API))
  }
  pp <- WorkTable%>%
    group_by()%>%
    dplyr::summarize(Ma = mean(Male), Fe = mean(Female), Wh = mean(W), Bl = mean(B), Ot = mean(O), AI = mean(AIE), AP = mean(API))
  info <- rbind(info, list(Male = pp$Ma, Female = pp$Fe, W = pp$Wh, B = pp$Bl, O = pp$Ot, AIE = pp$AI, API = pp$AP, method = 'CP'))
  
  return(info)
}
stopCluster(cl)

### Summarizing and plotting
Result <- Result[1:300,]
pp <- Result%>%
  group_by(method)%>%
  dplyr::summarize(Ma = mean(Male), Fe = mean(Female), Wh = mean(W), Bl = mean(B), Ot = mean(O), AI = mean(AIE), AP = mean(API))
pp

ppPlotG <- data.frame(class = rep(c('Female', 'Male'), each = 3), prop = c(pp$Fe, pp$Ma), method = rep(c('CP', 'PASS', 'PASSC'), 2))
ppPlotG$method <- factor(ppPlotG$method, levels = c("PASS", "PASSC", "CP"))
ppPlotR <- data.frame(class = rep(c("Amer-Indian-Eskimo", "Asian-Pac-Islander", "Black", "Other", "White"), each = 3), prop = c(pp$AI, pp$AP, pp$Bl, pp$Ot, pp$Wh), method = rep(c('CP', 'PASS', 'PASSC'), 5))
ppPlotR$method <- factor(ppPlotR$method, levels = c("PASS", "PASSC", "CP"))
ppPlotR$class <- factor(ppPlotR$class, levels = c("Amer-Indian-Eskimo", "Asian-Pac-Islander", "Black", "Other", "White"))

PG <- ggplot(data = ppPlotG, aes(x = method, y = prop, fill = class)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, width = 0.8) +
  scale_fill_manual(values = c("#E31A1C", "#1F78B4"), name = "") +
  scale_y_continuous(name = "Proportion") +
  scale_x_discrete(name = "Gender", labels = c("PASS", "PASSC", "CP")) +
  theme_bw() +
  theme(legend.position = "top", axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 1))

PR <- ggplot(data = ppPlotR, aes(x = method, y = prop, fill = class)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.8) +
  scale_fill_manual(values = c("#B15928", "#FF7F00", "#E31A1C", "#33A02C", "#1F78B4"), name = "", labels = c("AIE", "API", "B", "O", "W")) +
  scale_y_continuous(name = "Race") +
  scale_x_discrete(name = "Race", labels = c("PASS", "PASSC", "CP")) +
  theme_bw() +
  theme(legend.position = "top", axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 1))

PPGR <- ggarrange(PG, PR, common.legend = FALSE)
PPGR


### Education length and age
ns <- 200
cl <- makeCluster(8) ### number of parallel threads
registerDoParallel(cl)
Result <- foreach(iter = 1:ns, .combine = "rbind", .packages = c('MASS', 'glmnet', "caret", "kernlab", "randomForest", "ks", "nnet", "dplyr"), .errorhandling = "remove")%dopar% {
  info <- data.frame()
  
  
  ### sampling from the whole big data set (data generation)
  IndexSample <- c(sample(Null, round((n+N)*sigma), replace = FALSE), sample(Alter, n+N-round((n+N)*sigma), replace = FALSE))
  newdata <- data[sample(IndexSample, n+N, replace = FALSE),]
  
  datawork <- DataSplit(newdata, n+N, N, n_cal, n)
  data_train <- datawork$data_train
  data_cal <- datawork$data_cal
  data_label <- datawork$data_rest
  data_test <- datawork$data_test
  ClassData <- data1[rownames(data_test),]
  
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
  
  
  ### rescaling the unlabeled covariates and compute the similarity matrix
  X_test_scale <- scale(X_test, center = TRUE, scale = TRUE)
  X_test_scale[, 7:13] <- X_test_scale[, 7:13]/7
  X_test_scale[, 14:29] <- X_test_scale[, 14:29]/16
  X_test_scale[, 30:36] <- X_test_scale[, 30:36]/7
  X_test_scale[, 37:50] <- X_test_scale[, 37:50]/14
  X_test_scale[, 51:56] <- X_test_scale[, 51:56]/6
  X_test_scale[, 57:61] <- X_test_scale[, 57:61]/5
  X_test_scale[, 62:63] <- X_test_scale[, 62:63]/2
  
  
  ### model fitting and prediction
  model <- fitting(algo, X_train, Y_train, lambda)
  model_whole <- fitting(algo, X_label, Y_label, lambda)
  W_cal <- Pred(algo, model, X_cal)
  W_test <- Pred(algo, model, X_test)
  W_test_whole <- Pred(algo, model_whole, X_test)
  W_test_class <- ClassPred(W_test_whole, algo)
  
  #---PASS---# 
  TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
  Screened <- which(TN<=0.5)
  S <- Scompute(g, X_test_scale[Screened,], d, sgm = 1)
  Pi_PASS <- SolvePi(TN, m, S, Screened = Screened, alpha)
  
  #---CP---#
  pval <- confomalPvalue(W_cal, W_test, Null_cal, Value)
  rej <- sort(pval)<((1:N)/N)*alpha
  rejnum <- max(which(rej==T))
  reject <- which(pval<=sort(pval)[rejnum])
  Pi_CP <- rep(1/rejnum, rejnum)
  
  
  ### Counting and calculating proportions
  WorkTableAge <- data.frame()
  for (i in 1:100) {
    S_Pi <- sample(reject, m, replace = TRUE, Pi_CP)
    S_Pi <- S_Pi[is.element(S_Pi, Alter_test)]
    m0 <- length(S_Pi)
    HistReAge <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], age)
    WorkTableAge <- rbind(WorkTableAge, list(class = as.character(HistReAge$age), numb = HistReAge$n/m0,
                                             type = rep('CP', dim(HistReAge)[1]), Times = rep(as.character(i), dim(HistReAge)[1])))
    S_Pi <- sample(Screened, m, replace = TRUE, Pi_PASS)
    S_Pi <- S_Pi[is.element(S_Pi, Alter_test)]
    m0 <- length(S_Pi)
    HistReAge <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], age)
    WorkTableAge <- rbind(WorkTableAge, list(class = as.character(HistReAge$age), numb = HistReAge$n/m0,
                                             type = rep('PASS', dim(HistReAge)[1]), Times = rep(as.character(i), dim(HistReAge)[1])))
  }
  ppAge <- WorkTableAge%>%
    group_by(class, type)%>%
    dplyr::summarize(va = sum(numb)/100)
  ppAge$rt <- iter
  ppAge$att <- 'age'
  info <- rbind(info, ppAge)
  
  WorkTableEdu <- data.frame()
  for (i in 1:100) {
    S_Pi <- sample(reject, m, replace = TRUE, Pi_CP)
    S_Pi <- S_Pi[is.element(S_Pi, Alter_test)]
    m0 <- length(S_Pi)
    HistReEdu <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], education.num)
    WorkTableEdu <- rbind(WorkTableEdu, list(class = as.character(HistReEdu$education.num), numb = HistReEdu$n/m0,
                                             type = rep('CP', dim(HistReEdu)[1]), Times = rep(as.character(i), dim(HistReEdu)[1])))
    S_Pi <- sample(Screened, m, replace = TRUE, Pi_PASS)
    S_Pi <- S_Pi[is.element(S_Pi, Alter_test)]
    m0 <- length(S_Pi)
    HistReEdu <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], education.num)
    WorkTableEdu <- rbind(WorkTableEdu, list(class = as.character(HistReEdu$education.num), numb = HistReEdu$n/m0,
                                             type = rep('PASS', dim(HistReEdu)[1]), Times = rep(as.character(i), dim(HistReEdu)[1])))
  }
  ppEdu <- WorkTableEdu%>%
    group_by(class, type)%>%
    dplyr::summarize(va = sum(numb)/100)
  ppEdu$rt <- iter
  ppEdu$att <- 'edu'
  info <- rbind(info, ppEdu)
  
  return(info)
}
stopCluster(cl)


### Summarizing and plotting
rt <- intersect(Result$rt, Result$rt)
Result <- Result[Result$rt<=rt[100],]
ResultEdu <- Result[Result$att=='edu',]
ResultAge <- Result[Result$att=='age',]
pp <- ResultEdu%>%
  group_by(type, class)%>%
  dplyr::summarize(va = sum(va)/100)

pp$va[pp$type=="CP"] <- -pp$va[pp$type=="CP"]
pp$class <- as.numeric(pp$class)
pp$type <- factor(pp$type, levels = c("PASS", "CP"))

Pedu <- pp%>%ggplot(aes(x = class, y = va, fill = type)) +
  geom_col(alpha = 0.8) +
  scale_y_continuous(name = "Proportion") +
  theme_bw() +
  scale_x_continuous(name = "Education length") +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C"), labels = c('PASS', 'CP'), name = 'Method') +
  theme(axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"))

pp1 <- pp[pp$class<=8,]
Pedu1 <- pp1%>%ggplot(aes(x = class, y = va, fill = type)) +
  geom_col(alpha = 0.8) +
  scale_y_continuous(name = "Proportion") +
  theme_bw() +
  scale_x_continuous(name = "Education length") +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C"), labels = c('PASS', 'CP'), name = 'Method') +
  theme(axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"))

PPedu <- ggarrange(Pedu, Pedu1, common.legend = T, widths = c(1,1.04))
PPedu


pp <- ResultAge%>%
  group_by(type, class)%>%
  dplyr::summarize(va = sum(va)/ns)
pp$va[pp$type=="CP"] <- -pp$va[pp$type=="CP"]
pp$class <- as.numeric(pp$class)
pp$type <- factor(pp$type, levels = c("PASS", "CP"))

Page <- pp%>%ggplot(aes(x = class, y = va, fill = type)) +
  geom_col(width = 0.7, alpha = 0.8) +
  scale_y_continuous(name = "Proportion", limits = c(-0.053, 0.053)) +
  theme_bw() +
  scale_x_continuous(name = "Age") +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C"), labels = c('PASS', 'CP'), name = 'Method') +
  theme(axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"))
Page

pp1 <- pp[pp$class<30,]
Page1 <- pp1%>%ggplot(aes(x = class, y = va, fill = type)) +
  geom_col(width = 0.7, alpha = 0.8) +
  scale_y_continuous(name = "Proportion") +
  theme_bw() +
  scale_x_continuous(name = "Age", breaks = c(20, 23, 26, 29)) +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C"), labels = c('PASS', 'CP'), name = 'Method') +
  theme(axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"))

pp2 <- pp[pp$class>60,]
Page2 <- pp2%>%ggplot(aes(x = class, y = va, fill = type)) +
  geom_col(width = 0.7, alpha = 0.8) +
  scale_y_continuous(name = "Proportion") +
  theme_bw() +
  scale_x_continuous(name = "Age") +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C"), labels = c('PASS', 'CP'), name = 'Method') +
  theme(axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"))

pp3 <- pp[pp$class>=30&pp$class<=60,]
Page3 <- pp3%>%ggplot(aes(x = class, y = va, fill = type)) +
  geom_col(width = 0.7, alpha = 0.8) +
  scale_y_continuous(name = "Proportion") +
  theme_bw() +
  scale_x_continuous(name = "Age") +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C"), labels = c('PASS', 'CP'), name = 'Method') +
  theme(axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"))

PPage <- ggarrange(Page1, Page3, Page2, common.legend = T, ncol = 3)
PPage


