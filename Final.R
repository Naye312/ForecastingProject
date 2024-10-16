## Fix dataset by removing unneeded variables and removing outliers

data <- read.table("C:/Users/nyanl/Desktop/UCSD/Econ 178/Final/data_tr.txt", head = T)[,-1]

data_fixed <- subset(data, select= -c(hval, nohs)) ## removed hmort and no highschool

## Remove outliers 
Q <- quantile(data_fixed$tw, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(data_fixed$tw)
data_fixed <- subset(data_fixed, data_fixed$tw > (Q[1] - 1.5*iqr) & data_fixed$tw < (Q[2]+1.5*iqr))

## Histograms before and after outlier removal 
hist(data$tw)
hist(data_fixed$tw)

## Outcome 
y <- data_fixed$tw

set.seed(123)
library(MASS)
library(glmnet)

n <- length(y)
k <- 5
ii <- sample(rep(1:k, length= n))
pr.stepwise_backward <- pr.stepwise_forward <- pr.lasso <- pr.ridge <- pr.ols1 <- pr.ols2 <- pr.ols3 <- rep(NA, length(y))

## No transform
for (j in 1:k){
  hold <- (ii == j)
  train <- (ii != j)
  
  ## OLS
  regcross1 <- lm(tw ~ age + marr + inc + educ, data= data_fixed[train,])
  regcross2 <- lm(tw ~ 1, data= data_fixed[train,])
  regcross3 <- lm(tw ~ ., data= data_fixed[train,])
  
  
  pr.ols1[hold] <- predict(regcross1, newdata=data_fixed[hold,])
  pr.ols2[hold] <- predict(regcross2, newdata=data_fixed[hold,])
  pr.ols3[hold] <- predict(regcross3, newdata=data_fixed[hold,])
  
  ## Stepwise 
  full <- lm(tw ~ ., data=data_fixed[train,])
  null <- lm(tw ~ 1, data=data_fixed[train,])
  a <- stepAIC(null, scope=list(lower=null, upper=full), trace = FALSE, direction='forward')
  b <- stepAIC(full, scope=list(lower=null, upper=full), trace = FALSE, direction='backward')
  pr.stepwise_backward[hold] <- predict(b, newdata=data_fixed[hold,])
  pr.stepwise_forward[hold] <- predict(a, newdata=data_fixed[hold,])
  
  ## lasso and ridge
  xx.tr <- data_fixed[train,-1]
  y.tr <-  y[train]
  xx.te <- data_fixed[hold,-1]
  ridge.cv <- cv.glmnet(x=as.matrix(xx.tr), y=y.tr, nfolds=k, alpha=0)
  lasso.cv <- cv.glmnet(x=as.matrix(xx.tr), y=y.tr, nfolds=k, alpha=1)
  ridge <- glmnet(x = xx.tr, y = y.tr, lambda = ridge.cv$lambda.min, alpha = 0)
  lasso <- glmnet(x = xx.tr, y = y.tr, lambda = lasso.cv$lambda.min, alpha = 1)
  pr.lasso[hold] <- predict(lasso, newx=as.matrix(xx.te))
  pr.ridge[hold] <- predict(ridge, newx=as.matrix(xx.te))
}

mspe_step_backward <- mean((pr.stepwise_backward-y)^2)
mspe_step_forward <- mean((pr.stepwise_forward-y)^2)
mspe.Lasso <- mean((pr.lasso-y)^2)
mspe.ridge <- mean((pr.ridge-y)^2)
mspe.ols1 <- mean((pr.ols1-y)^2)
mspe.ols2 <- mean((pr.ols2-y)^2)
mspe.ols3 <- mean((pr.ols3-y)^2)

## Check the best
c(mspe.ols1, mspe.ols2, mspe.ols3, mspe_step_backward, mspe_step_forward, mspe.Lasso, mspe.ridge)
# 1739601805 2455317801  317342442  317357485  317387315  317383348  326766943
# OLS with all the predictors gives the best result when comparing just the linear models. 

## Polynomials
plot(data_fixed$nifa, data_fixed$tw)
plot(data_fixed$inc, data_fixed$tw)
X <- poly(data_fixed$nifa,3) 
X_poly <- data.frame(X)
y <- data_fixed$tw

n <- length(y) 
k <- 5
set.seed(123) 
id <- sample(rep(1:k, length= n)) 

p <- 10
mspe_degree <- mspe_degree2 <- mspe_degree_simple <- rep(NA, p)

## We'll start off with GAMs instead of single predictors since the degree results can be different.
for (deg in 1:p) {
  mspe <- vector(length = k)
  mspe2 <- vector(length = k)
  mspe3 <- vector(length = k)
  
  for (f in 1:k) {
    test <- (id == f)
    train <- (id != f)
    
    # Polynomial fit with deg degrees
    poly_reg <- lm(tw ~ . - nifa + poly(nifa,deg), data= data_fixed[train,])
    poly_reg2 <- lm(tw ~ . - inc + poly(inc,deg), data= data_fixed[train,])
    poly_reg_simple <- lm(tw ~ poly(nifa,deg), data= data_fixed[train,])
    
    # Make predictions
    pred <- predict(poly_reg, newdata = data_fixed[test,])
    pred2 <- predict(poly_reg2, newdata = data_fixed[test,])
    pred_simple <- predict(poly_reg_simple, newdata = data_fixed[test,])
    # Calculate MSPE
    mspe[f] <- mean((y[test] - pred)^2)
    mspe2[f] <- mean((y[test] - pred2)^2)
    mspe3[f] <- mean((y[test] - pred_simple)^2)
  }
  mspe_degree[deg] <- mean(mspe)
  mspe_degree2[deg] <- mean(mspe2)
  mspe_degree_simple[deg] <- mean(mspe3)
}

which.min(mspe_degree) # 5 is the best degree for nifa
which.min(mspe_degree2) # 4 is the best degree for inc 
which.min(mspe_degree_simple) # 10 is the best degree for nifa if the nifa polynomial is the only predictor


MSPE.ols <- MSPE.ols_poly <- MSPE.ols_single <- rep(NA, k)

for (f in 1:k){
  test <- (id == f)
  train <- (id != f)
  
  regcrosscontrol <- lm(tw ~ ., data= data_fixed[train,])
  prcontrol <- predict(regcrosscontrol, newdata=data_fixed[test,])
  MSPE.ols[f] <- mean((y[test] - prcontrol)^2)
  
  regcrosspoly <- lm(tw ~ . - nifa + poly(nifa,5) - inc + poly(inc,4), data= data_fixed[train,])
  prpoly <- predict(regcrosspoly, newdata=data_fixed[test,])
  MSPE.ols_poly[f] <- mean((y[test] - prpoly)^2)
  
  regcross <- lm(tw ~ poly(nifa,5), data= data_fixed[train,])
  pr <- predict(regcross, newdata=data_fixed[test,])
  MSPE.ols_single[f] <- mean((y[test] - pr)^2)
}

mean(MSPE.ols_poly) # 313854010 ## With transformed predictors
mean(MSPE.ols) # 317339936 ## This is the control when regressing with just plain predictors
mean(MSPE.ols_single) # 1650943961 ## This is another control with only a single transformed predictor, results in low accuracy. 

## Splines
plot(data_fixed$ira, data_fixed$tw)
plot(data_fixed$hmort, data_fixed$tw)

library(splines)
forward <- stepAIC(null, scope=list(lower=null, upper=full), trace = FALSE, direction='forward')
summary(forward)

y <- data_fixed$tw
n <- length(y)
set.seed(613)
k <- 5
id <- sample(rep(1:k, length= n))

# max number of df
p <- 20
mspe_knots <- mspe_knots2 <- rep(NA, p)

# finding the best degree
for (deg in 1:p) {
  mspe <- vector(length = k)
  mspe2 <-  vector(length = k)
  
  for (f in 1:k) {
    test <- (id == f)
    train <- (id != f)
    spline_reg <- lm(tw ~ . - hmort + bs(hmort, df = deg), data = data_fixed[train,])
    pred <- predict(spline_reg, newdata = data_fixed[test,])
    mspe[f] <- mean((y[test] - pred)^2)
    
    spline_reg2 <- lm(tw ~ . - ira + bs(ira, df = deg), data = data_fixed[train,])
    pred2 <- predict(spline_reg2, newdata = data_fixed[test,])
    mspe2[f] <- mean((y[test] - pred2)^2)
  }
  mspe_knots[deg] <- mean(mspe)
  mspe_knots2[deg] <- mean(mspe2)
} 
which.min(mspe_knots) # 3 degrees of freedom are the best for hmort
which.min(mspe_knots2) # 11 degrees of freedom are the best for ira

mspe_OLS <- mspe_OLS2 <- mspe_OLS_control <- vector(length = k)
for (f in 1:k) {
  test <- (id == f)
  train <- (id != f)
  
  ## OLS
  spline_reg <- lm(tw ~ . - hmort + bs(hmort, 3), data = data_fixed[train,])
  spline_reg2 <- lm(tw ~ . - ira + bs(ira, 11), data = data_fixed[train,])
  control_reg <- lm(tw ~ . , data = data_fixed[train,])
  pred <- predict(spline_reg, newdata = data_fixed[test,])
  pred2 <- predict(spline_reg2, newdata = data_fixed[test,])
  control_pred <- predict(control_reg, newdata = data_fixed[test,])
  mspe_OLS[f] <- mean((y[test] - pred)^2)
  mspe_OLS2[f] <- mean((y[test] - pred2)^2)
  mspe_OLS_control[f] <- mean((y[test] - control_pred)^2)
}

mean(mspe_OLS) # 317729111 ## Using B-spline on hmort worsens the accuracy
mean(mspe_OLS2) # 316370133 ## Using B-spline on ira increases accuracy 
mean(mspe_OLS_control) # 317547401 ## Control model without any splines

## Natural Cubic Spline 
p <- 20
mspe_knots <- rep(NA, p)

# finding the best degree
for (deg in 1:p) {
  mspe <- vector(length = k)
  
  for (f in 1:k) {
    test <- (id == f)
    train <- (id != f)
    spline_reg <- lm(tw ~ . - hmort + ns(hmort, df = deg), data = data_fixed[train,])
    pred <- predict(spline_reg, newdata = data_fixed[test,])
    mspe[f] <- mean((y[test] - pred)^2)
  }
  mspe_knots[deg] <- mean(mspe)
} 
which.min(mspe_knots) # 2 degrees of freedom are the best for hmort

mspe_OLS <- mspe_OLS_control <- vector(length = k)
for (f in 1:k) {
  test <- (id == f)
  train <- (id != f)
  
  ## OLS
  spline_reg <- lm(tw ~ . - hmort + ns(hmort, 2), data = data_fixed[train,])
  control_reg <- lm(tw ~ . , data = data_fixed[train,])
  pred <- predict(spline_reg, newdata = data_fixed[test,])
  control_pred <- predict(control_reg, newdata = data_fixed[test,])
  mspe_OLS[f] <- mean((y[test] - pred)^2)
  mspe_OLS_control[f] <- mean((y[test] - control_pred)^2)
}

mean(mspe_OLS) # 317514110 ## Using natural spline on hmort improves the accuracy
mean(mspe_OLS_control) # 317547401 ## Control

## Final Model
y <- data_fixed$tw
n <- length(y)
set.seed(613)
k <- 5
id <- sample(rep(1:k, length= n))
mspe_control <- mspe_OLS <- mspe_forward <- mspe_ridge <- mspe_lasso <- vector(length = k)
pr.stepwise_backward <- pr.stepwise_forward <- pr.lasso <- pr.ridge <- rep(NA, length(y))

hmort.spline <- as.data.frame(ns(data_fixed$hmort, 2))
ira.spline <- as.data.frame(bs(data_fixed$ira, 11))
inc.poly <- as.data.frame(poly(data_fixed$inc, 4))
nifa.poly <-as.data.frame(poly(data_fixed$nifa, 5))
names(hmort.spline) <- paste0("hmort_spline", 1:2)
names(ira.spline) <- paste0("hmort_spline", 1:11)
names(inc.poly) <- paste0("inc.poly", 1:4)
names(nifa.poly) <- paste0("nifa.poly", 1:5)


# Combine dataset with hmort.spline
data_transform <- cbind(data_fixed, hmort.spline)
data_transform <- cbind(data_transform, ira.spline)
data_transform <- cbind(data_transform, inc.poly)
data_transform <- cbind(data_transform, nifa.poly)

# Remove hmort from transformed data
data_transform <- subset(data_transform, select = -hmort)
data_transform <- subset(data_transform, select = -ira)
data_transform <- subset(data_transform, select = -inc)
data_transform <- subset(data_transform, select = -nifa)


for (f in 1:k) {
  test <- (id == f)
  train <- (id != f)
  
  ## OLS
  reg <- lm(tw ~ . - hmort + ns(hmort, 2) - ira + bs(ira, 11) - inc + poly(inc, 4) - nifa + poly(nifa, 5), data = data_fixed[train,])
  pred <- predict(reg, newdata = data_fixed[test,])
  mspe_OLS[f] <- mean((y[test] - pred)^2)
  
  reg2 <- lm(tw ~ ., data = data_fixed[train,])
  pred2 <- predict(reg2, newdata = data_fixed[test,])
  mspe_control[f] <- mean((y[test] - pred2)^2)
  
  ## Stepwise 
  full <- lm(tw ~ . - hmort + ns(hmort, 2) - ira + bs(ira, 11) - inc + poly(inc, 4) - nifa + poly(nifa, 5), data=data_fixed[train,])
  null <- lm(tw ~ 1, data=data_fixed[train,])
  a <- stepAIC(null, scope=list(lower=null, upper=full), trace = FALSE, direction='forward')
  b <- stepAIC(full, scope=list(lower=null, upper=full), trace = FALSE, direction='backward')
  pr.stepwise_backward[test] <- predict(b, newdata=data_fixed[test,])
  pr.stepwise_forward[test] <- predict(a, newdata=data_fixed[test,])
  
  ## lasso and ridge
  xx.tr <- data_transform[train,-1]
  y.tr <-  y[train]
  xx.te <- data_transform[test,-1]
  ridge.cv <- cv.glmnet(x=as.matrix(xx.tr), y=y.tr, nfolds=k, alpha=0)
  lasso.cv <- cv.glmnet(x=as.matrix(xx.tr), y=y.tr, nfolds=k, alpha=1)
  ridge <- glmnet(x = xx.tr, y = y.tr, lambda = ridge.cv$lambda.min, alpha = 0)
  lasso <- glmnet(x = xx.tr, y = y.tr, lambda = lasso.cv$lambda.min, alpha = 1)
  pr.lasso[test] <- predict(lasso, newx=as.matrix(xx.te))
  pr.ridge[test] <- predict(ridge, newx=as.matrix(xx.te))
}

mean(mspe_OLS) ## 312393646
mean((pr.stepwise_backward-y)^2) ## 312320914
mean((pr.stepwise_forward-y)^2) ## 312320914
mean((pr.lasso-y)^2) ## 312580608
mean((pr.ridge-y)^2) ## 322128696
mean(mspe_control) ## 317547401

## Stepwise regressions with polynomial transformations of income and nifa, along with splines of hmort and ira, is the best model.
full <- lm(tw ~ . - hmort + ns(hmort, 2) - ira + bs(ira, 11) - inc + poly(inc, 4) - nifa + poly(nifa, 5), data=data_fixed)
null <- lm(tw ~ 1, data=data_fixed)
my_model <- stepAIC(null, scope=list(lower=null, upper=full), trace = FALSE, direction='forward')

data_te <- read.table("C:/Users/nyanl/Desktop/UCSD/Econ 178/Final/data_for_prediction.txt", header = TRUE, sep = "\t", dec = ".")[,-1]
data_te <- subset(data_te, select= -c(hval, nohs)) ## removed hmort and no highschool

my_predictions <- predict(my_model, newdata = data.frame(data_te))
length(my_predictions)
write.table(my_predictions, file = 'C:/Users/nyanl/Desktop/UCSD/Econ 178/Final/my_predictions.txt')
