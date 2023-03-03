library(MASS)
library(ks)
require(ggplot2)
library(randomForest)
library(kernlab)
library(tidyverse)


source("functions_PASS.R")
source("algoclass_PASS.R")


### adjust the number of variables randomly sampled as candidates at each split of random forest
setMethod("fitting", "RFc", function(obj, X_rest, Y_rest, lambda){
  randomForest(X_rest, as.factor(Y_rest), mtry = 2, ntree = lambda)
})


N <- 1500 ### unlabeled data size
n <- 400 ### labeled data size
alpha <- 0.1 ### pre-specified FSR level
n_train <- round(n/2) ### training set size for single splitting
n_cal <- n-n_train ### calibration set size for single splitting
d <- 4 ### dimension of covariates
algo <- new("RFc") ### algorithm used for fitting and prediction
lambda <- 500 ### number of trees
m <- 100 ### sampling budget
g <- 'RBF' ### similarity criterion, RBF or COS
sigma <- 0.8 ### proportion of null individuals


### data generation and splitting
n1 <- round((n+N)*sigma)
X1 <- mvrnorm(n1, c(rep(3, round(d/4)), rep(0, d-round(d/4))), diag(rep(1, d)))
X2 <- mvrnorm((n+N)-n1, rep(0, d), diag(rep(1, d)))
X <- rbind(X1, X2)
data <- data.frame(x = X, y = c(rep(-1, n1), rep(1, (n+N)-n1)))
data <- data[sample(1:dim(data)[1], dim(data)[1]),]
Value <- list(type = "==A", v = -1)

datawork <- DataSplit(data[1:n,], n, 0, n_cal)
data_train <- datawork$data_train
data_cal <- datawork$data_cal
data_rest <- datawork$data_rest
data_test <- data[(n+1):(n+N),]

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


### model fitting and prediction
model <- fitting(algo, X_train, Y_train, lambda)
model_whole <- fitting(algo, X_rest, Y_rest, lambda)
W_cal <- Pred(algo, model, X_cal)
W_test <- Pred(algo, model, X_test)
W_test_whole <- Pred(algo, model_whole, X_test)


### implementation of PASS
TN <- Tcompute(W_cal, W_test, Null_cal, algo, h1 = 0, h2 = 0, IsSame = T, IsCali = T)
Screened <- which(TN<=0.5)
Pi_PASS <- SolvePi(TN, m, S, Screened = Screened, alpha)
Sampled <- sample(Screened, m, replace = TRUE, Pi_PASS)


### plotting
Type <- rep("Not-of-interest", N)
Type[-Null_test] <- "Of-interest"
Type[intersect(Sampled, Null_test)] <- "Falsely-sampled"
Type[intersect(Sampled, Alter_test)] <- "Correctly-sampled"
Type <- factor(Type, levels = c("Not-of-interest", "Of-interest", "Falsely-sampled", "Correctly-sampled"))
dataP <- data.frame(X_test1 = X_test[, 1], X_test2 = X_test[, 2], Types = Type)
dataP <- dataP[order(dataP$Types),]

P <- dataP %>% 
  ggplot(aes(x = X_test1, y = X_test2, color = Types, size = Types, shape = Types)) +
  geom_point() +
  scale_y_continuous(name = "X2", position = "left") +
  scale_x_continuous(name = "X1") +
  scale_color_manual(values = c('sandybrown', "lightgreen", "red", "darkgreen"), name = "", 
                     labels = c("Not-of-interest", "Of-interest", "Falsely-sampled", "Correctly-sampled")) +
  scale_size_manual(values = c(0.5, 0.5, 1.5, 1.5), labels = c("Not-of-interest", "Of-interest", "Falsely-sampled", "Correctly-sampled")) +
  scale_shape_manual(values = c(16, 16, 17, 16), name = "", labels = c("Not-of-interest", "Of-interest", "Falsely-sampled", "Correctly-sampled")) +
  theme(plot.title = element_text(color = 'black',hjust = 0.5),
        legend.position = "top", legend.spacing.x = unit(0.4, 'cm')) +
  guides(size = "none", color = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "top")
P
