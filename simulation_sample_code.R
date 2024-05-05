
library(ggplot2)
library(dplyr)
library(MASS)
library(splines)
library(lme4)
# library(lmerTest)
library(lsa)
library(simml)
library(lme4)
library(MASS)
library(lsa)
library(simml)
library(splines)
library(dplyr)
library(e1071)
library(randomForest)
library(datasets)
library(DTRlearn2)

# Load functions
source('data_generation.R')
source('functions.R')

################################################
# Parameter Setting
################################################

# Dimension of covariates
p <- 2
# Theta angle
theta_angle <- 5
# Number of subjects 
n <- 200
ntest <- 100
# Random error
error <- 1
# Mean value of covariates
mux <- p:1 / sqrt(sum((1:p)^2))
mux[1:floor(p/2)] <- -mux[1:floor(p/2)]
# Covariance matrix of covariates
sigmax <- matrix(0, p, p)
sigmax <- matrix((0.9 / p)^(abs(row(sigmax) - col(sigmax))), p, p)
diag(sigmax) <- 0.5 
# Covariate names
covar_names <- paste('X.', 1:p, sep='')
# Observation time 
t <- seq(0, 7, 1)
# Design matrix 
Xij <- cbind(1, t, t^2, t^3, t^4)
Zij <- cbind(1, t, t^2)
# True alpha
alpha <- 1:p
alpha <- alpha / sqrt(sum(alpha^2))
# Covariance matrix D
D.drg <- matrix(
  c(0.050, -0.010, -0.001,   
    -0.010, 0.050, -0.001,  
    -0.001, -0.001, 0.001),  
  nrow = 3,                  
  ncol = 3,                
  byrow = FALSE             
)
D.pbo <- matrix(
  c(0.040, -0.012, -0.001, 
    -0.012, 0.050, -0.001,  
    -0.001, -0.001, 0.001), 
  nrow = 3,                 
  ncol = 3,                
  byrow = FALSE         
)
beta.drg <- matrix(c(20, 3, -0.5), 3, 1)
beta.pbo <- matrix(c(20, 2.3, -0.4), 3, 1)
gamma.pbo <- matrix(c(0, cos(theta_angle/ 180 * pi), 
                      sin(theta_angle/ 180 * pi)), 3, 1)
gamma.drg <- matrix(c(0, cos(theta_angle/ 180 * pi), 
                      -sin(theta_angle/ 180 * pi)), 3, 1)

# # Non-quadratic trajectory
# pbo = function(t, w){
#   return(30*(cos(0.2*pi*t) - sin(0.2*pi*t*w))+30)
# }
# drg = function(t, w){
#   return(30*(cos(0.2*pi*t) + cos(0.1*pi*t*w) - 0.003*t))
# }

# Quadratic trajectory
pbo <- function(t, w){
  return(cbind(1, t, t^2) %*% (beta.pbo + gamma.pbo * c(w)))
}
drg <- function(t, w){
  return(cbind(1, t, t^2) %*% (beta.drg + gamma.drg * c(w)))
}

################################################
# Simulate
# one training and one test data set
################################################

# Initialize result vectors
result_lingem = alpha_lingem = c()
result_simml = alpha_simml = c()
result_owl = result_svmr = result_rf = c()  
result_mle = alpha_mle = c()  
result_pats = alpha_pats = c()
result_kld = alpha_kld = c()
result_npats = alpha_npats = c()

niters = 1 # 100

for(seed in niters){
  
  set.seed(seed)
  
  # Simulate training data
  data <- simulateData(n, mux, sigmax, alpha, D.pbo, D.drg, error, t, p)
  cv_train = data
  
  # Compute change scores
  uni.data <- computeChangeScores(data)
  
  # Get the estimated mean and covariance matrix for X
  covarlist = c()
  for(i in 1:p){
    covarlist = cbind(covarlist, data[,paste('X.',i, sep = '')])
  }
  estmux = matrix(apply(unique(covarlist), 2, mean), p, 1)
  estsigmax = cov(unique(covarlist))

  # Simulate test data
  data.test <- simulateData(ntest, mux, sigmax, alpha, D.pbo, D.drg, error, t, p)
  # Compute change scores
  data.test <- computeChangeScores(data.test)
  xtest = data.test[,paste('X.', 1:p, sep='')]
  xtest = as.matrix(xtest)
  
  ## Display results
  # head(data.test)
  # range(data.test$cs)
  
  ################################################
  # Utilize different methods to calculate 
  # PCD and Value
  ################################################
  
  ########################
  # Linear GEM
  ########################
  
  # Fit linear model
  lm_formula <- as.formula('cs ~ -1 + X.1')
  for (i in paste('X.', 2:p, sep = '')) {
    lm_formula <- update(lm_formula, paste('~', i))
  }
  
  # Fit linear models
  l_est_drg <- lm(lm_formula, data = uni.data[uni.data$group == 'drg', ])
  l_est_pbo <- lm(lm_formula, data = uni.data[uni.data$group == 'pbo', ])
  
  # Test results
  test_drg <- predict(l_est_drg, newdata = data.test[, paste('X.', 1:p, sep = '')])
  test_pbo <- predict(l_est_pbo, newdata = data.test[, paste('X.', 1:p, sep = '')])
  lin_group <- ifelse(test_drg > test_pbo, 'pbo', 'drg')
  
  # Compute alpha
  beta1 <- coef(l_est_drg)
  beta2 <- coef(l_est_pbo)
  alphatemp <- (beta1 - beta2) / sqrt(t(beta1 - beta2) %*% estsigmax %*% (beta1 - beta2))[1][1]
  alphatemp <- alphatemp / sqrt(sum(alphatemp^2))
  alphatemp <- matrix(alphatemp, length(alphatemp), 1)
  
  # Store alpha
  alpha_lingem <- rbind(alpha_lingem, c(alphatemp))
  colnames(alpha_lingem) <- paste('a', 1:p, sep = '')
  
  # Store results
  result_lingem <- rbind(result_lingem, c(
    sum(lin_group == data.test$truegroup.cs) / nrow(data.test),
    abs(cosine(c(alpha), c(alphatemp))),
    mean(data.test[data.test$group == lin_group, ]$cs, na.rm = TRUE)
  ))
  colnames(result_lingem) <- c('PCD', 'Cosine', 'Value')
  
  
  ########################
  # SIMML
  ########################
  
  # Extracting data
  yy <- uni.data$cs
  trt <- ifelse(uni.data$group == 'drg', 1, 2)
  xx <- matrix(unlist(uni.data[, paste('X.', 1:p, sep ='')]), dim(uni.data)[1], p)
  xx_test <- matrix(unlist(data.test[, paste('X.', 1:p, sep ='')]), dim(data.test)[1], p)
  
  # Fitting simml
  glm_fit <- glm(yy ~ xx, family = "gaussian")
  mu_hat <- as.vector(predict(glm_fit, newX = xx, type = "link"))
  simml_obj <- simml(yy, trt, xx, Xm = mu_hat, family = "gaussian")
  
  # Predicting treatment rules
  simml_trt_rule <- pred.simml(simml_obj, newX = xx_test, maximize = FALSE)$trt.rule
  simml_trt_rule <- ifelse(simml_trt_rule == 1, 'drg', 'pbo')
  
  # Storing results
  temp <- simml_obj$beta.coef / sqrt(sum(simml_obj$beta.coef^2))
  result_simml <- rbind(result_simml, c(
    sum(simml_trt_rule == data.test$truegroup.cs) / nrow(data.test),
    cosine(c(alpha), c(temp)),
    mean(data.test[data.test$group == simml_trt_rule, ]$cs, na.rm = TRUE)
  ))
  alpha_simml <- rbind(alpha_simml, c(temp))
  colnames(alpha_simml) <- paste('a', 1:p, sep = '')
  colnames(result_simml) = c('PCD', 'Cosine', 'Value')
  
  
  ########################
  # OWL-Gaussian
  ########################
  
  # Prepare treatment indicator
  treatment_indicator <- ifelse(uni.data$group == 'drg', 1, -1)
  
  # Response variable
  response_variable <- uni.data$cs
  
  # Predictor matrix
  predictor_matrix <- as.matrix(uni.data[, paste('X.', 1:p, sep = '')])
  
  # Fit OWL-Gaussian model
  owl_model <- owl(H = predictor_matrix, 
                   AA = treatment_indicator,
                   RR = response_variable, 
                   n = length(response_variable),
                   K = 1,
                   pi = rep(0.5, n),
                   res.lasso = TRUE, loss = 'hinge', kernel = 'linear',
                   augment = FALSE, c = 2^(-2:2), 
                   sigma = c(0.01, 0.1, 1, 10, 100),
                   s = 2^(-2:2), m = 4)
  
  # Prepare test predictor matrix
  test_predictor_matrix <- as.matrix(data.test[, paste('X.', 1:p, sep = '')])
  
  # Make predictions
  predictions <- predict(owl_model, H = test_predictor_matrix, K = 1)
  treatment_rule <- ifelse(predictions$treatment[[1]] == -1, 'drg', 'pbo')
  
  # Compute performance metrics
  result_owl <- rbind(result_owl, c(
    mean(data.test$truegroup.cs == treatment_rule, na.rm = TRUE), 
    mean(data.test[data.test$group == treatment_rule, ]$cs, na.rm = TRUE)
  ))
  colnames(result_owl) <- c('PCD', 'Value')

  
  ########################
  # SVM-Radial
  ########################
  
  # Subset data for 'pbo' group
  pbo_data <- uni.data[uni.data$group == 'pbo', c(paste('X.', 1:p, sep = ''), 'cs')]
  
  # Subset data for 'drg' group
  drg_data <- uni.data[uni.data$group == 'drg', c(paste('X.', 1:p, sep = ''), 'cs')]
  
  # Fit SVM model for 'pbo' group
  svm_fit_pbo <- svm(cs ~ ., data = pbo_data, kernel = "radial",
                     cost = 10, scale = FALSE)
  psvm_pbo <- predict(svm_fit_pbo, xtest)
  
  # Fit SVM model for 'drg' group
  svm_fit_drg <- svm(cs ~ ., data = drg_data, kernel = "radial",
                     cost = 10, scale = FALSE)
  psvm_drg <- predict(svm_fit_drg, xtest)
  
  # Determine treatment rule based on SVM predictions
  psvm <- ifelse(psvm_pbo > psvm_drg, 'drg', 'pbo')
  
  result_svmr = rbind(result_svmr, c(
    mean(data.test$truegroup.cs == psvm, na.rm = TRUE),
    mean(data.test[data.test$group == psvm, ]$cs, na.rm = TRUE)
  ))
  
  colnames(result_svmr) = c('PCD', 'Value')

  
  ########################
  # Random Forest
  ########################
  
  # Fit Random Forest model for 'pbo' group
  rf_fit_pbo <- randomForest(cs ~ ., data = pbo_data, proximity = TRUE)
  prf_pbo <- predict(rf_fit_pbo, xtest)
  
  # Fit Random Forest model for 'drg' group
  rf_fit_drg <- randomForest(cs ~ ., data = drg_data, proximity = TRUE)
  prf_drg <- predict(rf_fit_drg, xtest)
  
  # Determine treatment rule based on Random Forest predictions
  prf <- ifelse(prf_pbo > prf_drg, 'drg', 'pbo')
  
  result_rf = rbind(result_rf, c(
    mean(data.test$truegroup.cs == prf, na.rm = TRUE), 
    mean(data.test[data.test$group == prf, ]$cs, na.rm = TRUE)
  ))
  
  colnames(result_rf) = c('PCD', 'Value')

  ########################
  # MLE
  ########################

  # Generate random initial alpha and normalize
  temp <- rnorm(2)
  temp <- temp / sqrt(sum(temp^2))
  
  # Calculate alpha using Maximum Likelihood Estimation (MLE)
  temp <- alphaLL(temp)
  temp <- temp / sqrt(sum(temp^2))
  alphatemp = temp
  
  # Calculate performance metrics
  result_tmp <- pcdcalculation(temp, data, data.test)
  
  # Store alpha values
  alpha_mle <- rbind(alpha_mle, matrix(temp, 1, p))
  colnames(alpha_mle) <- paste('a', 1:p, sep = '')
  
  # Store results
  result_mle <- rbind(result_mle, c(result_tmp$cs, abs(cosine(alpha, temp)), result_tmp$value))
  colnames(result_mle) <- c('PCD', 'Cosine', 'Value')
  
  
  ########################
  # PATS
  ########################

  # Optimize alpha using PATS
  resats1 <- optim(alphatemp, alphaATS, control = list(reltol = 1e-4, maxit = 100))
  
  # Extract optimized alpha and normalize
  temp <- resats1$par
  temp <- temp / sqrt(sum(temp^2))
  
  # Calculate performance metrics
  result_tmp <- pcdcalculation(temp, data, data.test)
  
  # Store alpha values
  alpha_pats <- rbind(alpha_pats, matrix(temp,1,p))
  colnames(alpha_pats) <- paste('a', 1:p, sep = '')
  
  # Store results
  result_pats <- rbind(result_pats, c(result_tmp$cs, abs(cosine(alpha, temp)), result_tmp$value))
  colnames(result_pats) <- c('PCD', 'Cosine', 'Value')
  
  ########################
  # LS-KLD
  ########################
  
  # Optimize alpha using LS-KLD
  reskld <- optim(alphatemp, alphaKL, control = list(reltol = 1e-2, maxit = 100))
  
  # Extract optimized alpha
  temp <- reskld$par
  temp <- temp / sqrt(sum(temp^2))
  
  # Calculate performance metrics using KLD-optimized alpha
  result_tmp <- pcdcalculation(temp, data, data.test)
  
  # Store alpha values
  alpha_kld <- rbind(alpha_kld, matrix(temp,1,p))
  colnames(alpha_kld) <- paste('a', 1:p, sep = '')
  
  # Store results
  result_kld <- rbind(result_kld, c(result_tmp$cs, abs(cosine(alpha, temp)), result_tmp$value))
  colnames(result_kld) <- c('PCD', 'Cosine', 'Value')
  
  
  ########################
  # NPATS
  ########################
  
  # Prepare data
  datatemp <- data
  
  # Optimize alpha using NPATS
  npats <- optim(alphatemp, alphaNPATS, control = list(reltol = 1e-2, maxit = 100))
  temp <- npats$par
  temp <- temp / sqrt(sum(temp^2))
  
  # Calculate performance metrics using NPATS
  result_tmp <- pcdcalculation_nonquad(temp, data, data.test)
  
  # Store alpha values
  alpha_npats <- rbind(alpha_npats, matrix(temp,1,p))
  colnames(alpha_npats) <- paste('a', 1:p, sep = '')
  
  # Store results
  result_npats <- rbind(result_npats, c(result_tmp$cs, abs(cosine(alpha, temp)), result_tmp$value))
  colnames(result_npats) <- c('PCD', 'Cosine', 'Value')

  
}

result_lingem
result_simml
result_owl
result_mle
result_pats
result_kld
result_npats 

