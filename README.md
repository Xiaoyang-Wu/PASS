# PASS
Codes for reproducing the numerical results in *Optimal Subsampling via Predictive Inference (2023)*. And contents of the R package implementing the proposed method.

# Description of PASS
In the big data era, subsampling or sub-data selection techniques are often adopted to extract a fraction of informative individuals from the massive data. Existing subsampling algorithms focus mainly on obtaining a representative subset to achieve the best estimation accuracy under a given class of models. In this paper, we consider a semi-supervised setting wherein a small or moderate sized “labeled” data is available in addition to a much larger sized “unlabeled” data. The goal is to sample from the unlabeled data with a given budget to obtain informative individuals that are characterized by their unobserved responses. We propose an optimal subsampling procedure
that is able to maximize the diversity of the selected subsample and control the false selection rate (FSR) simultaneously, allowing us to explore reliable information as
much as possible. The key ingredients of our method are the use of predictive inference for quantifying the uncertainty of response predictions and a reformulation of the objective into a constrained optimization problem.

# Contents

- **R**: Package function definitions.
- **man**: Package .Rd fils.
- **Simulation Codes**: Codes for reproducing the illustration example plot and the simulated results.
- **Real-data**: Codes for reproducing the real-data experiment and the data set used.

# How to install the PASS package
Use the R package `devtools` to install:
```
devtools::install_github("Xiaoyang-Wu/PASS")  
library(PASS)
```

# Notes on usage of the PASS package
For direct implementation of the method, just use the `PASS_sampling` function in the package and all arguments needed can be found in the help documentation `?PASS_sampling`.  
Note that the `algo` argument is of type `class` and should be defined by the user depending on the specific fitting and predicting algorithm chosen. Here is an example for the definition of the algorithm class:
```
library(PASS)
library(MASS)
library(randomForest)

setGeneric("fitting", function(obj, ...) standardGeneric("fitting"))
setGeneric("Pred", function(obj, ...) standardGeneric("pred"))
###random forest classification
setClass("RFc", slots = list(name = "character"), prototype = list(name = "RFc"))
setMethod("fitting", "RFc", function(obj, X_rest, Y_rest, lambda){
  randomForest(X_rest, as.factor(Y_rest), ntree = lambda)
})
setMethod("Pred", "RFc", function(obj, model, X_test, type = 'decision'){
  if(type=="class"){
    predraw <- as.numeric(as.character(predict(model, X_test)))
  }else{
    predraw <- predict(model, X_test, type = "prob")[, 2]
  }
  return(predraw)
})


N <- 5000 ### unlabeled data size
n <- 400 ### labeled data size
alpha <- 0.1 ### pre-specified FSR level
K <- 5 ### number of folds for PASSC
n_train <- round(n/2) ### training set size for single splitting
n_cal <- n-n_train ### calibration set size for single splitting
d <- 4 ### dimension of covariates
sgm <- 1 ### sigma value of the RBF criterion
algo <- new('RFc') ### algorithm used for fitting and prediction
lambda <- 500 ### number of trees
m <- 200 ### sampling budget
g <- 'RBF' ### similarity criterion, RBF or COS
sigma <- 0.8 ### proportion of null individuals

n1 <- round((n+N)*sigma)
X1 <- mvrnorm(n1, c(rep(3, round(d/4)), rep(0, d-round(d/4))), diag(rep(1, d)))
X2 <- mvrnorm((n+N)-n1, c(rep(0, d/2), rep(2, d/2)), diag(rep(1, d)))
X <- rbind(X1, X2)
data <- data.frame(x = X, y = c(rep(-1, n1), rep(1, (n+N)-n1)))
data <- data[sample(1:dim(data)[1], dim(data)[1]),]
region <- list(type = "==A", v = -1)
data_label <- data[1:n,]
data_test <- data[(n+1):(n+N),]
X_label <- as.matrix(data_label[colnames(data_label)[-d-1]])
Y_label <- as.matrix(data_label$y)
X_test <- as.matrix(data_test[colnames(data_test)[-d-1]])
Null_label <- which(Y_label==-1)

Sampling_results <- PASS_sampling(X_label, Y_label, X_test, Null_label, g, m, alpha, algo = algo, lambda = lambda, classification = TRUE, region = region)
```
The function `fitting` should take the algorithm class object, a covariate matrix, the corresponding response vector and optionally a tuning parameter as inputs and outputs a model object. The function `Pred` should take the algorithm class object, the fitted model object and a new covariate matrix as inputs and outputs a vector or prediction values.
