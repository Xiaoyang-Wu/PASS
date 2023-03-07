library(methods)


setGeneric("fitting", function(obj, ...) standardGeneric("fitting"))
setGeneric("Pred", function(obj, ...) standardGeneric("Pred"))
setGeneric("Tuning", function(obj, ...) standardGeneric("Tuning"))


###ridge
setClass("ridge", slots = list(name = "character", alpha = "numeric", family = "character"), prototype = list(name = "RR", alpha = 0, family = "gaussian") )
setMethod("fitting", "ridge", function(obj, X, Y, lambda){
  glmnet(x = X, y = Y, family = obj@family, alpha = obj@alpha, lambda = lambda)
})
setMethod("Pred", "ridge", function(obj, model, X_test){
  predict(model, X_test)
})


###lasso
setClass("lasso", slots = list(name = "character", alpha = "numeric", family = "character"), prototype = list(name = "Lasso", alpha = 1, family = "gaussian") )
setMethod("fitting", "lasso", function(obj, X, Y, lambda){
  glmnet(x = X, y = Y, family ='binomial', alpha = obj@alpha, lambda = lambda)
})
setMethod("Pred", "lasso", function(obj, model, X_test){
  predict(model, X_test)
})


###random forest regression
setClass("RF", slots = list(name = "character"), prototype = list(name = "RF"))
setMethod("fitting", "RF", function(obj, X, Y, lambda){
  datawork <- data.frame(X, y = Y)
  randomForest(y~., data = datawork, mtry = 3, ntree = lambda)
})
setMethod("Pred", "RF", function(obj, model, X_test){
  predict(model, as.data.frame(x = X_test))
})


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


###NN regression
setClass("NN-R", slots = list(name = "character"),
         prototype = list(name = "NN-R"))
setMethod("fitting", "NN-R", function(obj, X, Y, lambda){
  datawork <- data.frame(X, y = Y)
  nnet(y~., data = datawork, size = 10, linout = T, maxit = 2000)
})
setMethod("Pred", "NN-R", function(obj, model, X_test){
  predict(model, as.data.frame(x = X_test))})


###NN classification
setClass("NN", slots = list(name = "character"),
         prototype = list(name = "NN"))
setMethod("fitting", "NN", function(obj, X, Y, lambda){
  datawork <- data.frame(X, y = Y)
  nnet(y~., data = datawork, size = 10, decay = 5e-4, entropy = TRUE, maxit = 2000)
})
setMethod("Pred", "NN", function(obj, model, X_test, class = F){
  if(class==T){
    predict(model, as.data.frame(x = X_test), type = 'class')
  }else{
    predict(model, as.data.frame(x = X_test))
  }
})


#SVM regression
setClass("SVM-R", slots = list(name = "character"),
         prototype = list(name = "SVM-R"))
setMethod("fitting", "SVM-R", function(obj, X, Y, lambda){
  datawork <- data.frame(X, y = Y)
  ksvm(y~., data = datawork, C = lambda)
})
setMethod("Pred", "SVM-R", function(obj, model, X_test, type = "decision"){
  predict(model, X_test)
})


###SVM classification
setClass("SVM", slots = list(name = "character"),
         prototype = list(name = "SVM"))
setMethod("fitting", "SVM", function(obj, X, Y, lambda){
  datawork <- data.frame(X, y = as.factor(Y))
  ksvm(y~., data = datawork, type = "C-svc", C = lambda)
})
setMethod("Pred", "SVM", function(obj, model, X_test, type = "decision"){
  if(type=="class"){
    predraw <- as.numeric(as.character(predict(model, X_test)))
  }else{
    predraw <- predict(model, X_test, type = "decision")
  }
  return(predraw)
})


###linear regression
setClass("LRs", slots = list(name = "character"), prototype = list(name = "LR-standard") )
setMethod("fitting", "LRs", function(obj, X, Y, lambda){
  if(dim(X)[1]>=dim(X)[2]){
    datawork <- data.frame(X, y = Y)
    return(lm(y~., datawork))
  }else{
    X_a1 <- cbind(rep(1, dim(X)[1]), X)
    return(ginv(t(X_a1)%*%X_a1)%*%t(X_a1)%*%Y)
  }
})
setMethod("Pred", "LRs", function(obj, model, X_test){
  if(class(model)=="lm"){
    predict(model, as.data.frame(x = X_test))
  }else{
    X_a1 <- cbind(rep(1, dim(X_test)[1]), X_test)
    X_a1%*%model
  }
})


###GLM(logistic regression)
setClass("glml", slots = list(name = "character", family = "character"), prototype = list(name = "GLM", family = "binomial"))
setMethod("fitting", "glml", function(obj, X, Y, lambda){
  datawork <- data.frame(X, y = Y)
  glm(y~., datawork, family = binomial())
})
setMethod("Pred", "glml", function(obj, model, X_test){
  predict(model, as.data.frame(x = X_test), type = "response")
})
