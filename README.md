# PASS
Codes for reproducing the numerical results in *Optimal Subsampling via Predictive Inference (2023)* by Xiaoyang Wu, Yuyang Huo, Haojie Ren and Changliang Zou. And contents of the R package implementing the proposed method.

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
#SVM regression  
setClass("SVM-R", slots = list(name = "character"),  
         prototype = list(name = "SVM-R"))  
setMethod("fitting", "SVM-R", function(obj, X, Y, lambda){  
  datawork <- data.frame(X, y = Y)  
  ksvm(y~., data = datawork, C = lambda)  
})  
setMethod("Pred", "SVM-R", function(obj, model, X_test){  
  predict(model, X_test)  
})  
algo <- new('SVM-R')  
Sampling_results <- PASS_sampling(..., algo = algo, ...)
```
The function `fitting` should take the algorithm class object, a covariate matrix, the corresponding response vector and optionally a tuning parameter as inputs and outputs a model object. The function `Pred` should take the algorithm class object, the fitted model object and a new covariate matrix as inputs and outputs a vector or prediction values.
