library(kernlab)
library(MASS)
library(ks)
library(randomForest)
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


### sampling from the whole big data set (data generation)
Value <- list(type = "==A", v = -1)
Null <- NullIndex(data$y, Value)
Alter <- setdiff(1:Number, Null)
IndexSample <- c(sample(Null, round((n+N)*sigma), replace = FALSE), sample(Alter, n+N-round((n+N)*sigma), replace = FALSE))
newdata <- data[sample(IndexSample, n+N, replace = FALSE),]

datawork <- DataSplit(newdata, n+N, N, n_cal, n)
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
X_test_scale[, 7:13] <- X_test_scale[, 7:13]/7
X_test_scale[, 14:29] <- X_test_scale[, 14:29]/16
X_test_scale[, 30:36] <- X_test_scale[, 30:36]/7
X_test_scale[, 37:50] <- X_test_scale[, 37:50]/14
X_test_scale[, 51:56] <- X_test_scale[, 51:56]/6
X_test_scale[, 57:61] <- X_test_scale[, 57:61]/5
X_test_scale[, 62:63] <- X_test_scale[, 62:63]/2
S <- Scompute(g, X_test_scale, d)


### model fitting and prediction
model <- fitting(algo, X_train, Y_train, lambda)
model_whole <- fitting(algo, X_rest, Y_rest, lambda)
W_cal <- Pred(algo, model, X_cal)
W_test <- Pred(algo, model, X_test)
W_test_whole <- Pred(algo, model_whole, X_test)
W_test_class <- ClassPred(W_test_whole, algo)


### implementation of different methods
#---PASS---# 
TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
Screened <- which(TN<=0.5)
Pi_PASS <- SolvePi(TN, m, S, Screened = Screened, alpha)


#---PASSC---#
TN_cv <- KfoldTcompute(X_rest, Y_rest, X_test, Null_rest, W_test_whole, K, lambda, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
Screened_cv <- which(TN_cv<=0.5)
Pi_PASSC <- SolvePi(TN_cv, m, S, Screened_cv, alpha)


#---CP---#
pval <- confomalPvalue(W_cal, W_test, Null_cal, Value)
rej <- sort(pval)<((1:N)/N)*alpha
rejnum <- max(which(rej==T))
reject <- which(pval<=sort(pval)[rejnum])
Pi_CP <- rep(1/rejnum, rejnum)


#---SSD---#
Screened_pred <- setdiff(1:N, NullIndex(W_test_class, Value)) 
Pi_SSD <- SolvePi(TN, m, S, Screened_pred, alpha, constraint = FALSE)


#---SS---#
Pi_SS <- rep(1/length(Screened_pred), length(Screened_pred))


#---LocF---#
Pi_loc <- rep(1/length(Screened), length(Screened))


### plotting
ClassData <- data1[rownames(data_test),]


### gender
Pi1 <- Pi_PASSC
Pi2 <- Pi_CP
Index1 <- Screened_cv
Index2 <- reject
Medname1 <- "PASSC"
Medname2 <- "CP"
WorkTable <- data.frame()
for (i in 1:100) {
  S_Pi <- sample(Index1, m, replace = TRUE, Pi1)
  HistRe <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], sex)
  WorkTable <- rbind(WorkTable, list(class = as.character(HistRe$sex), numb = HistRe$n,
                                 type = rep(Medname1, dim(HistRe)[1]), Times = rep(i, dim(HistRe)[1])))
  S_Pi <- sample(Index2, m, replace = TRUE, Pi2)
  HistRe <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], sex)
  WorkTable <- rbind(WorkTable, list(class = as.character(HistRe$sex), numb = HistRe$n,
                                 type = rep(Medname2, dim(HistRe)[1]), Times = rep(as.character(i), dim(HistRe)[1])))
  S_Pi <- sample(Screened, m, replace = TRUE, Pi_PASS)
  HistRe <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], sex)
  WorkTable <- rbind(WorkTable, list(class = as.character(HistRe$sex), numb = HistRe$n,
                                 type = rep("PASS", dim(HistRe)[1]), Times = rep(as.character(i), dim(HistRe)[1])))
}
pp2 <- WorkTable%>%
  group_by(class, type)%>%
  dplyr::summarize(va = mean(numb))
pp2$va[c(1, 4)] <- pp2$va[c(1, 4)]/sum(pp2$va[c(1, 4)])
pp2$va[c(2, 5)] <- pp2$va[c(2, 5)]/sum(pp2$va[c(2, 5)])
pp2$va[c(3, 6)] <- pp2$va[c(3, 6)]/sum(pp2$va[c(3, 6)])
pp2$type <- factor(pp2$type, levels = c("PASS", "PASSC", "CP"))

P2 <- ggplot(data = pp2, aes(x = type, y = va, fill = class)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, width = 0.8) +
  scale_fill_manual(values = c("#E31A1C", "#1F78B4"), name = "") +
  scale_y_continuous(name = "Proportion") +
  scale_x_discrete(name = "Gender", labels = c("PASS", "PASSC", "CP")) +
  theme_bw() +
  theme(legend.position = "top", axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 1))


### race
WorkTable <- data.frame()
for (i in 1:100) {
  S_Pi <- sample(Index1, m, replace = TRUE, Pi1)
  HistRe <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], race)
  WorkTable <- rbind(WorkTable, list(class = as.character(HistRe$race), numb = HistRe$n,
                                 type = rep(Medname1, dim(HistRe)[1]), Times = rep(i, dim(HistRe)[1])))
  S_Pi <- sample(Index2, m, replace = TRUE, Pi2)
  HistRe <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], race)
  WorkTable <- rbind(WorkTable, list(class = as.character(HistRe$race), numb = HistRe$n,
                                 type = rep(Medname2, dim(HistRe)[1]), Times = rep(as.character(i), dim(HistRe)[1])))
  S_Pi <- sample(Screened, m, replace = TRUE, Pi_PASS)
  HistRe <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], race)
  WorkTable <- rbind(WorkTable, list(class = as.character(HistRe$race), numb = HistRe$n,
                                 type = rep("PASS", dim(HistRe)[1]), Times = rep(as.character(i), dim(HistRe)[1])))
 }
pp3 <- WorkTable%>%
  group_by(class, type)%>%
  dplyr::summarize(va = mean(numb))
pp3$type <- factor(pp3$type, levels = c("PASS", "PASSC", "CP"))
pp3$class <- factor(pp3$class, levels = c("Amer-Indian-Eskimo", "Asian-Pac-Islander", "Black", "Other", "White"))
pp3$va[c(1, 4, 7, 10, 13)] <- pp3$va[c(1, 4, 7, 10, 13)]/sum(pp3$va[c(1, 4, 7, 10, 13)])
pp3$va[c(2, 5, 8, 11, 14)] <- pp3$va[c(2, 5, 8, 11, 14)]/sum(pp3$va[c(2, 5, 8, 11, 14)])
pp3$va[c(3, 6, 9, 12, 15)] <- pp3$va[c(3, 6, 9, 12, 15)]/sum(pp3$va[c(3, 6, 9, 12, 15)])

P3 <- ggplot(data = pp3, aes(x = type, y = va, fill = class)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c("#B15928", "#FF7F00", "#E31A1C", "#33A02C", "#1F78B4"), name = "", labels = c("AIE", "API", "B", "O", "W")) +
  scale_y_continuous(name = "Porpotion") +
  theme_bw() +
  scale_x_discrete(name = "Method", labels = c("PASS",  "PASSC", "CP")) +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 1))
P3


PP1 <- ggarrange(P2, P3, common.legend = FALSE)
PP1


##education length
Pi1 <- Pi_PASS
Pi2 <- Pi_CP
Index1 <- Screened
Index2 <- reject
Medname1 <- "PASS"
Medname2 <- "CP"
WorkTable <- data.frame()
for (i in 1:100) {
  S_Pi <- sample(Index2, m, replace = TRUE, Pi2)
  HistRe <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], education.num)
  WorkTable <- rbind(WorkTable, list(class = as.character(HistRe$education.num), numb = HistRe$n,
                                 type = rep(Medname2, dim(HistRe)[1]), Times = rep(as.character(i), dim(HistRe)[1])))
  S_Pi <- sample(Index1, m, replace = TRUE, Pi1)
  HistRe <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], education.num)
  WorkTable <- rbind(WorkTable, list(class = as.character(HistRe$education.num), numb = HistRe$n,
                                 type = rep(Medname1, dim(HistRe)[1]), Times = rep(i, dim(HistRe)[1])))
}
pp4 <- WorkTable%>%
  group_by(class, type)%>%
  dplyr::summarize(va = mean(numb))
pp4$va[pp4$type=="CP"] <- -pp4$va[pp4$type=="CP"]
pp4$class <- as.numeric(pp4$class)
pp4$type <- factor(pp4$type, levels = c("PASS", "CP"))

P4 <- pp4%>%ggplot(aes(x = class, y = va, fill = type)) +
  geom_col(alpha = 0.8) +
  scale_y_continuous(name = "Number") +
  theme_bw() +
  scale_x_continuous(name = "Education length") +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C"), labels = c('PASS', 'CP'), name = 'Method') +
  theme(axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"))
P4


##age
WorkTable <- data.frame()
for (i in 1:100) {
  S_Pi <- sample(Index2, m, replace = TRUE, Pi2)
  HistRe <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], age)
  WorkTable <- rbind(WorkTable, list(class = as.character(HistRe$age), numb = HistRe$n,
                                 type = rep(Medname2, dim(HistRe)[1]), Times = rep(as.character(i), dim(HistRe)[1])))
  S_Pi <- sample(Index1, m, replace = TRUE, Pi1)
  HistRe <- count(ClassData[S_Pi,][is.element(S_Pi, Alter_test),], age)
  WorkTable <- rbind(WorkTable, list(class = as.character(HistRe$age), numb = HistRe$n,
                                 type = rep(Medname1, dim(HistRe)[1]), Times = rep(i, dim(HistRe)[1])))
}
pp5 <- WorkTable%>%
  group_by(class, type)%>%
  dplyr::summarize(va = mean(numb))
pp5$va[pp5$type=="CP"] <- -pp5$va[pp5$type=="CP"]
pp5$class <- as.numeric(pp5$class)
pp5$type <- factor(pp5$type, levels = c("PASS", "CP"))

P5 <- pp5%>%ggplot(aes(x = class, y = va, fill = type)) +
  geom_col(width = 0.7, alpha = 0.8) +
  scale_y_continuous(name = "Number") +
  theme_bw() +
  scale_x_continuous(name = "Age") +
  scale_fill_manual(values = c("#1F78B4", "#E31A1C"), labels = c('PASS', 'CP'), name = 'Method') +
  theme(axis.title.x = element_text(colour = "grey20", size = 14, face = "bold"))
P5

PP2 <- ggarrange(P4, P5, ncol = 2, common.legend = TRUE)
PP2
