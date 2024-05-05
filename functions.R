library(dplyr)
library(MASS)
library(splines)
library(lme4)
library(lsa)
library(simml)
library(ggplot2)

# Function to control and handle error
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

# Function for optimizing alpha with PATS
alphaATS = function(alpha){
  
  # Normalize alpha vector
  alpha = alpha / sqrt(sum(alpha^2))
  alpha = matrix(alpha, p, 1)
  
  # Prepare data and compute model fits for each treatment group
  datatemp = data
  datatemp$W = covarlist%*%alpha
  dat_pbo_est = datatemp[datatemp$group == 'pbo', ]
  dat_drg_est = datatemp[datatemp$group == 'drg', ]
  
  fit_pbo_est = myTryCatch(lmer(Y ~ week + I(week^2) + W + W * week +
                                  W * I(week^2) + (week+I(week^2)|subj),
                                data = dat_pbo_est, REML = FALSE))
  fit_drg_est = myTryCatch(lmer(Y ~ week + I(week^2) + W + W * week +
                                  W * I(week^2) + (week+I(week^2)|subj),
                                data = dat_drg_est, REML = FALSE))
  
  # Extract fixed effects coefficients
  beta1 = as.matrix(fixef(fit_drg_est$value))[1:3]
  gamma1 = as.matrix(fixef(fit_drg_est$value))[4:6] 
  # D1 = as.matrix(VarCorr(fit_drg_est$value)$subj)[2:3, 2:3] 
  beta2 = as.matrix(fixef(fit_pbo_est$value))[1:3] 
  gamma2 = as.matrix(fixef(fit_pbo_est$value))[4:6] 
  # D2 = as.matrix(VarCorr(fit_pbo_est$value)$subj)[2:3, 2:3] 
  
  # Compute the parameters
  t = sort(unique((data$week))) 
  gdiff = (cbind(1, t, t^2)[length(t),] - cbind(1, t, t^2)[1,])/ (max(t) - min(t))
  gdiff = matrix(gdiff, 3, 1)
  c1 = t(beta1 - beta2) %*% (gdiff) %*% t(gdiff) %*% (beta1 - beta2)
  c2 = 2*t(beta1 - beta2) %*% (gdiff) %*% t(gdiff) %*% (gamma1 - gamma2)
  c3 = t(gamma1 - gamma2) %*% (gdiff) %*% t(gdiff) %*% (gamma1 - gamma2)
  
  # Compute the criterion to maximize
  res = c1 + c2 * c(t(estmux) %*% alpha) + c3 * c(
    t(alpha) %*% (estmux %*% t(estmux)+ estsigmax) %*% alpha)
  
  # print(res)
  # print(c(alpha))
  return(-res) # Return the negative of the criterion value for optimization
}

# Function for optimizing alpha with KLD
alphaKL <- function(alpha) {
  
  # Normalize alpha vector
  alpha <- alpha / sqrt(sum(alpha^2))
  alpha <- matrix(alpha, nrow = p, ncol = 1)
  
  # Prepare the training data
  datatemp <- cv_train
  datatemp$W <- as.matrix(covarlist) %*% alpha
  
  # Subset data for pbo and drg groups
  dat_pbo_est <- datatemp[datatemp$group == 'pbo',]
  dat_drg_est <- datatemp[datatemp$group == 'drg',]
  
  # Model fitting for pbo and drg groups
  fit_pbo_est <- myTryCatch(lmer(Y ~ week + I(week^2) + W + W * week +
                                   W * I(week^2) + (week + I(week^2) | subj),
                                 data = dat_pbo_est, REML = FALSE))
  fit_drg_est <- myTryCatch(lmer(Y ~ week + I(week^2) + W + W * week +
                                   W * I(week^2) + (week + I(week^2) | subj),
                                 data = dat_drg_est, REML = FALSE))
  
  # Extract fixed effects and covariance matrices
  beta1 <- fixef(fit_drg_est$value)[2:3]
  gamma1 <- fixef(fit_drg_est$value)[5:6]
  D1 <- VarCorr(fit_drg_est$value)$subj[2:3, 2:3]
  
  beta2 <- fixef(fit_pbo_est$value)[2:3]
  gamma2 <- fixef(fit_pbo_est$value)[5:6]
  D2 <- VarCorr(fit_pbo_est$value)$subj[2:3, 2:3]
  
  # Calculate KL divergence components
  a1 <- -dim(D1)[1] + 0.5 * (sum(diag(solve(D1) %*% D2)) + sum(diag(solve(D2) %*% D1))) +
    0.5 * t(beta1 - beta2) %*% (solve(D1) + solve(D2)) %*% (beta1 - beta2)
  a2 <- t(beta1 - beta2) %*% (solve(D1) + solve(D2)) %*% (gamma1 - gamma2)
  a3 <- 0.5 * t(gamma1 - gamma2) %*% (solve(D1) + solve(D2)) %*% (gamma1 - gamma2)
  
  # Calculate and print result
  res <- a1 + a2 * (t(estmux) %*% alpha) + a3 * (t(alpha) %*% (estmux %*% t(estmux) + estsigmax) %*% alpha)
  print(res)
  
  return(-res)
}

# Function for optimizing alpha with MLE
alphaLL <- function(alpha) {
  
  # Normalize alpha vector
  alpha <- alpha / sqrt(sum(alpha^2))
  alpha <- matrix(alpha, nrow = p, ncol = 1)
  
  # Data preparation
  datatemp <- cv_train
  idpbo <- unique(datatemp[datatemp$group == 'pbo', ]$subj)
  iddrg <- unique(datatemp[datatemp$group == 'drg', ]$subj)
  
  for (j in 1:100) {
    
    # Store the old alpha for convergence checking
    alpha.old <- alpha
    
    # Update weights in the data
    datatemp$W <- as.matrix(covarlist) %*% alpha
    dat_pbo_est <- datatemp[datatemp$group == 'pbo', ]
    dat_drg_est <- datatemp[datatemp$group == 'drg', ]
    
    # Fit models for pbo and drg groups
    fit_pbo_est <- myTryCatch(lmer(Y ~ week + I(week^2) + W + W * week + W * I(week^2) + (week + I(week^2) | subj),
                                   data = dat_pbo_est, REML = FALSE))
    fit_drg_est <- myTryCatch(lmer(Y ~ week + I(week^2) + W + W * week + W * I(week^2) + (week + I(week^2) | subj),
                                   data = dat_drg_est, REML = FALSE))
    
    # Extract model estimates
    beta1 <- fixef(fit_pbo_est$value)[1:3]
    gamma1 <- fixef(fit_pbo_est$value)[4:6]
    D1 <- VarCorr(fit_pbo_est$value)$subj[1:3, 1:3]
    sigma1 <- summary(fit_pbo_est$value)$sigma
    
    beta2 <- fixef(fit_drg_est$value)[1:3]
    gamma2 <- fixef(fit_drg_est$value)[4:6]
    D2 <- VarCorr(fit_drg_est$value)$subj[1:3, 1:3]
    sigma2 <- summary(fit_drg_est$value)$sigma
    
    part1 <- 0
    part2 <- 0
    
    # Calculate parts for the pbo group
    for (i in idpbo) {
      t <- datatemp[datatemp$subj == i, 'week']
      Xij <- cbind(1, t, t^2)
      Sigma.yj <- Xij %*% D1 %*% t(Xij) + diag(length(t)) * sigma1^2
      
      yij <- datatemp[datatemp$subj == i, 'Y']
      xij <- matrix(unlist(datatemp[datatemp$subj == i, covar_names][1,]), p, 1)
      part1 <- part1 + c(t(gamma1) %*% t(Xij) %*% solve(Sigma.yj) %*% (yij - Xij %*% beta1)) * t(xij)
      part2 <- part2 + c(t(gamma1) %*% t(Xij) %*% solve(Sigma.yj) %*% Xij %*% gamma1) * xij %*% t(xij)
    }
    
    # Calculate parts for the drg group
    for (i in iddrg) {
      t <- datatemp[datatemp$subj == i, 'week']
      Xij <- cbind(1, t, t^2)
      Sigma.yj <- Xij %*% D2 %*% t(Xij) + diag(length(t)) * sigma2^2
      
      yij <- datatemp[datatemp$subj == i, 'Y']
      xij <- matrix(unlist(datatemp[datatemp$subj == i, covar_names][1,]), p, 1)
      part1 <- part1 + c(t(gamma2) %*% t(Xij) %*% solve(Sigma.yj) %*% (yij - Xij %*% beta2)) * t(xij)
      part2 <- part2 + c(t(gamma2) %*% t(Xij) %*% solve(Sigma.yj) %*% Xij %*% gamma2) * xij %*% t(xij)
    }
    
    # Update alpha and check for convergence
    alpha <- solve(part2) %*% t(part1)
    alpha <- alpha / sqrt(sum(alpha^2))
    
    if (abs(cosine(as.vector(alpha), as.vector(alpha.old))) > 0.9999) {
      break
    }
  }
  
  return(alpha)
}

# Function for optimizing alpha with NPATS
alphaNPATS = function(alphatmp) {
  
  # Compute W matrix
  datatemp$W = covarlist %*% alphatmp
  
  # Subset data for each treatment group
  dat_pbo_est = datatemp[datatemp$group == 'pbo', ]
  dat_drg_est = datatemp[datatemp$group == 'drg', ]
  
  # Fit linear mixed-effects models for each group
  fit.pbo.formula = lmer(Y ~ bs(week, knots = c(3, 5)) * bs(W, knots = c(-1, 1)) + (week + I(week^2) | subj),
                         data = dat_pbo_est, REML = FALSE)
  fit.drg.formula = lmer(Y ~ bs(week, knots = c(3, 5)) * bs(W, knots = c(-1, 1)) + (week + I(week^2) | subj),
                         data = dat_drg_est, REML = FALSE)
  
  # Initialize variables for calculations
  part1 = part2 = 0
  eta.pbo = eta.drg = numeric()
  
  # Perform calculations for each subject in pbo group
  for (i in unique(dat_pbo_est$subj)) {
    temp = dat_pbo_est[dat_pbo_est$subj == i, ]
    wi_pbo = c(t(alphatmp) %*% matrix(unlist(temp[1, paste('X.', 1:p, sep = '')]), p, 1))
    Amatrix = c()
    tt = sort(temp$week)
    for (ii in tt) {
      Amatrix = rbind(Amatrix, t(A(ii)))
    }
    X_alphai_pbo = Amatrix %x% t(Bu(wi_pbo))
    Yi_pbo = matrix(temp$Y, length(temp$Y), 1)
    
    Z.pbo = cbind(1, tt, tt^2)
    Dhat.pbo = as.matrix(VarCorr(fit.pbo.formula)$subj)
    Z.drg = cbind(1, tt, tt^2)
    Dhat.drg = as.matrix(VarCorr(fit.drg.formula)$subj)
    Vhat.pbo = Z.pbo %*% Dhat.pbo %*% t(Z.pbo) + summary(fit.pbo.formula)$sigma^2 * diag(dim(Z.pbo)[1])
    Vhat.drg = Z.drg %*% Dhat.drg %*% t(Z.drg) + summary(fit.drg.formula)$sigma^2 * diag(dim(Z.drg)[1])
    
    part1 = part1 + t(X_alphai_pbo) %*% solve(Vhat.pbo) %*% X_alphai_pbo
    part2 = part2 + t(X_alphai_pbo) %*% solve(Vhat.pbo) %*% Yi_pbo
  }
  
  # Calculate eta for pbo group
  eta.pbo = ginv(part1) %*% part2
  
  # Reset variables for calculations
  part1 = part2 = 0
  
  # Perform calculations for each subject in drg group
  for (i in unique(dat_drg_est$subj)) {
    temp = dat_drg_est[dat_drg_est$subj == i, ]
    wi_drg = c(t(alphatmp) %*% matrix(unlist(temp[1, paste('X.', 1:p, sep = '')]), p, 1))
    Amatrix = c()
    tt = sort(temp$week)
    for (ii in tt) {
      Amatrix = rbind(Amatrix, t(A(ii)))
    }
    X_alphai_drg = Amatrix %x% t(Bu(wi_drg))
    Yi_drg = matrix(temp$Y, length(temp$Y), 1)
    
    Z.pbo = cbind(1, tt, tt^2)
    Dhat.pbo = as.matrix(VarCorr(fit.pbo.formula)$subj)
    Z.drg = cbind(1, tt, tt^2)
    Dhat.drg = as.matrix(VarCorr(fit.drg.formula)$subj)
    Vhat.pbo = Z.pbo %*% Dhat.pbo %*% t(Z.pbo) + summary(fit.pbo.formula)$sigma^2 * diag(dim(Z.pbo)[1])
    Vhat.drg = Z.drg %*% Dhat.drg %*% t(Z.drg) + summary(fit.drg.formula)$sigma^2 * diag(dim(Z.drg)[1])
    
    part1 = part1 + t(X_alphai_drg) %*% solve(Vhat.drg) %*% X_alphai_drg
    part2 = part2 + t(X_alphai_drg) %*% solve(Vhat.drg) %*% Yi_drg
  }
  
  # Calculate eta for drg group
  eta.drg = ginv(part1) %*% part2
  
  # Calculate Astar
  tm = max(tt); t1 = min(tt)
  Astar = (A(tm) - A(t1)) %x% t(A(tm) - A(t1))
  
  # Define function G
  G = function(u, alphatmp = alphatmp, eta.drg = eta.drg, eta.pbo = eta.pbo) {
    alphatmp = matrix(alphatmp, p, 1)
    tmp = myTryCatch(t(eta.drg - eta.pbo) %*% (Astar %x% (Bu(u) %*% t(Bu(u)))) %*% (eta.drg - eta.pbo) * f_alpha_fun(alphatmp, u))
    if (is.null(tmp$error) == FALSE) {
      tmp = 0
    } else {
      tmp = tmp$value
    }
    if (length(tmp) > 1) {
      print(tmp)
    }
    return(c(tmp))
  }
  
  # Print eta.drg and alphatmp
  print(eta.drg[1:3])
  alphatmp = alphatmp / sqrt(sum(alphatmp^2))
  alphatmp = matrix(alphatmp, p, 1)
  print(alphatmp)
  
  # Integrate function G
  res = myTryCatch(integrate(Vectorize(G, vectorize.args = 'u'), upper = Inf, lower = -Inf, alphatmp = alphatmp, eta.drg = eta.drg, eta.pbo = eta.pbo, subdivisions = 10000))
  tmp = res$value$value
  print(tmp)
  
  # Return result
  return(-tmp)
}

# Function for calculating performance metrics (pcd)
pcdcalculation = function(alpha, data, data.test) {
  # Prepare data for estimation
  datatemp = data
  datatemp$W = covarlist %*% alpha
  dat_pbo_est = datatemp[datatemp$group == 'pbo', ]
  dat_drg_est = datatemp[datatemp$group == 'drg', ]
  
  # Fit linear mixed-effects models for each treatment group
  fit_pbo_est = myTryCatch(
    lmer(Y ~ week + I(week^2) + W + W * week + W * I(week^2) + (week + I(week^2) | subj),
         data = dat_pbo_est, REML = FALSE)
  )
  fit_drg_est = myTryCatch(
    lmer(Y ~ week + I(week^2) + W + W * week + W * I(week^2) + (week + I(week^2) | subj),
         data = dat_drg_est, REML = FALSE)
  )
  
  # Extract fixed effects coefficients for each group
  beta1 = as.matrix(fixef(fit_drg_est$value))[1:3]
  gamma1 = as.matrix(fixef(fit_drg_est$value))[4:6] 
  beta2 = as.matrix(fixef(fit_pbo_est$value))[1:3] 
  gamma2 = as.matrix(fixef(fit_pbo_est$value))[4:6] 
  
  # Calculate difference in coefficients between the groups
  t = sort(unique(data$week)) # Sort and get unique week values
  gdiff = (cbind(1, t, t^2)[length(t),] - cbind(1, t, t^2)[1,]) / (max(t) - min(t))
  gdiff = matrix(gdiff, 3, 1)
  
  # Estimate treatment group for each subject in the test data
  estgroup = c()
  for (i in unique(data.test$subj)) {
    temp = c(matrix(unlist(data.test[i, paste('X.', 1:p, sep = '')]), 1, p) %*% alpha)
    if (t(gdiff) %*% (beta1 + temp * gamma1) > t(gdiff) %*% (beta2 + temp * gamma2)) {
      estgroup = c(estgroup, 'pbo')
    } else {
      estgroup = c(estgroup, 'drg')
    }
  }
  
  # Calculate performance metrics
  result = list(
    cs = sum(estgroup == data.test$truegroup.cs) / length(estgroup),
    value = mean(data.test[data.test$group == estgroup, ]$cs, na.rm = TRUE)
  )
  
  return(result)
}

pcdcalculation_nonquad = function(alpha, data, data.test) {
  # Prepare data for estimation
  xtest = data.test[,paste('X.', 1:p, sep='')]
  xtest = as.matrix(xtest)
  estW = xtest %*% alpha
  datatemp <- data
  datatemp$W <- covarlist %*% alpha
  dat_pbo_est <- datatemp[datatemp$group == 'pbo', ]
  dat_drg_est <- datatemp[datatemp$group == 'drg', ]
  
  # Fit linear mixed-effects models for each treatment group
  fit.pbo_form <- lmer(Y ~ bs(week, knots = c(3,5)) *
                         bs(W, knots = c(-1,1)) +
                         (week + I(week^2) | subj),
                       data = dat_pbo_est, REML = FALSE)
  fit.drg_form <- lmer(Y ~ bs(week, knots = c(3,5)) *
                         bs(W, knots = c(-1,1)) +
                         (week + I(week^2) | subj),
                       data = dat_drg_est, REML = FALSE)
  
  # Create design matrices for week 0 and week 7
  dat_for_design <- data.frame(week = rep(t, ntest),
                               Y = rnorm(length(t) * ntest),
                               W = rep(unique(estW), each = length(t)))
  design <- model.matrix(lm(Y ~ bs(week, knots = c(3,5)) *
                              bs(W, knots = c(-1,1)),
                            data = dat_for_design))
  design <- cbind(rep(t, ntest), design)
  design0 <- design[design[,1] == 0, 2:dim(design)[2]]
  design7 <- design[design[,1] == 7, 2:dim(design)[2]]
  rownames(design0) <- rownames(design7) <- NULL
  
  design0 <- data.frame(t(design0))
  design0$name <- rownames(design0)
  design7 <- data.frame(t(design7))
  design7$name <- rownames(design7)
  
  # Extract coefficients
  tmpcoef_drg <- data.frame(name = rownames(summary(fit.drg_form)$coefficient),
                            coef = summary(fit.drg_form)$coefficient[,1])
  tmpcoef_pbo <- data.frame(name = rownames(summary(fit.pbo_form)$coefficient),
                            coef = summary(fit.pbo_form)$coefficient[,1])
  tmpdf_drg0 <- merge(design0, tmpcoef_drg, by = 'name')
  tmpdf_drg7 <- merge(design7, tmpcoef_drg, by = 'name')
  
  # Calculate ats for drg
  ats_drg0 <- ats_drg7 <- c()
  for(i in 2:(dim(tmpdf_drg0)[2] - 1)){
    ats_drg0 <- c(ats_drg0, c(t(tmpdf_drg0[,i]) %*% tmpdf_drg0[,'coef']))
    ats_drg7 <- c(ats_drg7, c(t(tmpdf_drg7[,i]) %*% tmpdf_drg7[,'coef']))
  }
  ats_drg <- ats_drg7 - ats_drg0
  
  # Calculate ats for pbo
  tmpdf_pbo0 <- merge(design0, tmpcoef_pbo, by = 'name')
  tmpdf_pbo7 <- merge(design7, tmpcoef_pbo, by = 'name')
  
  ats_pbo0 <- ats_pbo7 <- c()
  for(i in 2:(dim(tmpdf_pbo0)[2] - 1)){
    ats_pbo0 <- c(ats_pbo0, c(t(tmpdf_pbo0[,i]) %*% tmpdf_pbo0[,'coef']))
    ats_pbo7 <- c(ats_pbo7, c(t(tmpdf_pbo7[,i]) %*% tmpdf_pbo7[,'coef']))
  }
  ats_pbo <- ats_pbo7 - ats_pbo0
  
  # Determine estimated treatment group
  estgroup <- ifelse(ats_drg > ats_pbo, 'pbo', 'drg')
  
  # Compute performance metrics
  res_npats_tmp <- data.frame(cs = sum(estgroup == data.test$truegroup.cs) / length(estgroup), 
                              value = mean(data.test[data.test$group == estgroup,]$cs, na.rm = TRUE))
  return(res_npats_tmp)
}

# Function for computing log-likelihood values for alpha
alphaLLvalue = function(alpha) {
  # Normalize alpha
  alpha = alpha / sqrt(sum(alpha^2))
  alpha = matrix(alpha, p, 1)
  
  # Compute W matrix
  datatemp = data
  datatemp$W = covarlist %*% alpha
  
  # Extract unique subject IDs for each treatment group
  idpbo = unique(datatemp[datatemp$group == 'pbo', ]$subj)
  iddrg = unique(datatemp[datatemp$group == 'drg', ]$subj)
  
  # Subset data for each treatment group
  dat_pbo_est = datatemp[datatemp$group == 'pbo', ]
  dat_drg_est = datatemp[datatemp$group == 'drg', ]
  
  # Fit linear mixed-effects models for each group
  fit_pbo_est = myTryCatch(lmer(Y ~ week + I(week^2) + W + W * week + W * I(week^2) + (week + I(week^2) | subj),
                                data = dat_pbo_est, REML = FALSE))
  fit_drg_est = myTryCatch(lmer(Y ~ week + I(week^2) + W + W * week + W * I(week^2) + (week + I(week^2) | subj),
                                data = dat_drg_est, REML = FALSE))
  
  # Return log-likelihood values
  return(c(summary(fit_pbo_est$value)$logLik + summary(fit_drg_est$value)$logLik))
}

# Function for computing f_alpha function
f_alpha_fun = function(alpha, u){
  tmp.mux = c(t(alpha) %*% estmux)
  tmp.sigmax = c(t(alpha) %*% estsigmax %*% alpha)
  return(1/sqrt(2*pi*tmp.sigmax^2) * exp(-1/(2*tmp.sigmax^2)*(u-tmp.mux)^2))
}

# Calculate basis function for A()
x = seq(0,7,0.01)
basis = bs(x, knots = c(3,5))
part1 = 1:length(seq(0,3,0.01)); 
part2 = length(seq(0,3,0.01)):length(seq(0,5,0.01)); 
part3 = length(seq(0,5,0.01)):length(seq(0,7,0.01))
x1 = x[part1]; x2 = x[part2]; x3 = x[part3]

fit_basis1_part1 = lm(basis[part1, 1] ~ x1 + I(x1^2) + I(x1^3))
fit_basis2_part1 = lm(basis[part1, 2] ~ x1 + I(x1^2) + I(x1^3))
fit_basis3_part1 = lm(basis[part1, 3] ~ x1 + I(x1^2) + I(x1^3))
fit_basis4_part1 = lm(basis[part1, 4] ~ x1 + I(x1^2) + I(x1^3))
fit_basis5_part1 = lm(basis[part1, 5] ~ x1 + I(x1^2) + I(x1^3))

fit_basis1_part2 = lm(basis[part2, 1] ~ x2 + I(x2^2) + I(x2^3))
fit_basis2_part2 = lm(basis[part2, 2] ~ x2 + I(x2^2) + I(x2^3))
fit_basis3_part2 = lm(basis[part2, 3] ~ x2 + I(x2^2) + I(x2^3))
fit_basis4_part2 = lm(basis[part2, 4] ~ x2 + I(x2^2) + I(x2^3))
fit_basis5_part2 = lm(basis[part2, 5] ~ x2 + I(x2^2) + I(x2^3))

fit_basis1_part3 = lm(basis[part3, 1] ~ x3 + I(x3^2) + I(x3^3))
fit_basis2_part3 = lm(basis[part3, 2] ~ x3 + I(x3^2) + I(x3^3))
fit_basis3_part3 = lm(basis[part3, 3] ~ x3 + I(x3^2) + I(x3^3))
fit_basis4_part3 = lm(basis[part3, 4] ~ x3 + I(x3^2) + I(x3^3))
fit_basis5_part3 = lm(basis[part3, 5] ~ x3 + I(x3^2) + I(x3^3))

coefmatrix = rbind(fractions(coef(fit_basis1_part1)),
                   fractions(coef(fit_basis1_part2)),
                   fractions(coef(fit_basis1_part3)),
                   fractions(coef(fit_basis2_part1)),
                   fractions(coef(fit_basis2_part2)),
                   fractions(coef(fit_basis2_part3)),
                   fractions(coef(fit_basis3_part1)),
                   fractions(coef(fit_basis3_part2)),
                   fractions(coef(fit_basis3_part3)),
                   fractions(coef(fit_basis4_part1)),
                   fractions(coef(fit_basis4_part2)),
                   fractions(coef(fit_basis4_part3)),
                   fractions(coef(fit_basis5_part1)),
                   fractions(coef(fit_basis5_part2)),
                   fractions(coef(fit_basis5_part3)))

f_b1_p1 = function(x){
  return(coefmatrix[1,1] + coefmatrix[1,2]*x + coefmatrix[1,3]*x^2 + coefmatrix[1,4]*x^3)
}
f_b1_p2 = function(x){
  return(coefmatrix[2,1] + coefmatrix[2,2]*x + coefmatrix[2,3]*x^2 + coefmatrix[2,4]*x^3)
}
f_b1_p3 = function(x){
  return(coefmatrix[3,1] + coefmatrix[3,2]*x + coefmatrix[3,3]*x^2 + coefmatrix[3,4]*x^3)
}

f_b2_p1 = function(x){
  tmpdif = 3
  return(coefmatrix[1+tmpdif,1] + coefmatrix[1+tmpdif,2]*x + 
           coefmatrix[1+tmpdif,3]*x^2 + coefmatrix[1+tmpdif,4]*x^3)
}
f_b2_p2 = function(x){
  tmpdif = 3
  return(coefmatrix[2+tmpdif,1] + coefmatrix[2+tmpdif,2]*x + 
           coefmatrix[2+tmpdif,3]*x^2 + coefmatrix[2+tmpdif,4]*x^3)
}
f_b2_p3 = function(x){
  tmpdif = 3
  return(coefmatrix[3+tmpdif,1] + coefmatrix[3+tmpdif,2]*x + 
           coefmatrix[3+tmpdif,3]*x^2 + coefmatrix[3+tmpdif,4]*x^3)
}

f_b3_p1 = function(x){
  tmpdif = 6
  return(coefmatrix[1+tmpdif,1] + coefmatrix[1+tmpdif,2]*x + 
           coefmatrix[1+tmpdif,3]*x^2 + coefmatrix[1+tmpdif,4]*x^3)
}
f_b3_p2 = function(x){
  tmpdif = 6
  return(coefmatrix[2+tmpdif,1] + coefmatrix[2+tmpdif,2]*x + 
           coefmatrix[2+tmpdif,3]*x^2 + coefmatrix[2+tmpdif,4]*x^3)
}
f_b3_p3 = function(x){
  tmpdif = 6
  return(coefmatrix[3+tmpdif,1] + coefmatrix[3+tmpdif,2]*x + 
           coefmatrix[3+tmpdif,3]*x^2 + coefmatrix[3+tmpdif,4]*x^3)
}

f_b4_p1 = function(x){
  tmpdif = 9
  return(coefmatrix[1+tmpdif,1] + coefmatrix[1+tmpdif,2]*x + 
           coefmatrix[1+tmpdif,3]*x^2 + coefmatrix[1+tmpdif,4]*x^3)
}
f_b4_p2 = function(x){
  tmpdif = 9
  return(coefmatrix[2+tmpdif,1] + coefmatrix[2+tmpdif,2]*x + 
           coefmatrix[2+tmpdif,3]*x^2 + coefmatrix[2+tmpdif,4]*x^3)
}
f_b4_p3 = function(x){
  tmpdif = 9
  return(coefmatrix[3+tmpdif,1] + coefmatrix[3+tmpdif,2]*x + 
           coefmatrix[3+tmpdif,3]*x^2 + coefmatrix[3+tmpdif,4]*x^3)
}

f_b5_p1 = function(x){
  tmpdif = 12
  return(coefmatrix[1+tmpdif,1] + coefmatrix[1+tmpdif,2]*x + 
           coefmatrix[1+tmpdif,3]*x^2 + coefmatrix[1+tmpdif,4]*x^3)
}
f_b5_p2 = function(x){
  tmpdif = 12
  return(coefmatrix[2+tmpdif,1] + coefmatrix[2+tmpdif,2]*x + 
           coefmatrix[2+tmpdif,3]*x^2 + coefmatrix[2+tmpdif,4]*x^3)
}
f_b5_p3 = function(x){
  tmpdif = 12
  return(coefmatrix[3+tmpdif,1] + coefmatrix[3+tmpdif,2]*x + 
           coefmatrix[3+tmpdif,3]*x^2 + coefmatrix[3+tmpdif,4]*x^3)
}

A = function(x){
  if(x < 3){
    tmp = c(f_b1_p1(x), f_b2_p1(x), f_b3_p1((x)), f_b4_p1(x), f_b5_p1(x))
  }else if(x >= 3 & x < 5){
    tmp = c(f_b1_p2(x), f_b2_p2(x), f_b3_p2((x)), f_b4_p2(x), f_b5_p2(x))
  }else{
    tmp = c(f_b1_p3(x), f_b2_p3(x), f_b3_p3((x)), f_b4_p3(x), f_b5_p3(x))
  }
  tmp = c(1, tmp)
  tmp = matrix(tmp, length(tmp), 1)
  return(tmp)
}


# Basis function for weight
x = seq(-3,3,0.01)
basis = bs(x, knots = c(-1,1))
part1 = 1:length(seq(-3,-1,0.01)); 
part2 = length(seq(-3,-1,0.01)):length(seq(-3,1,0.01)); 
part3 = length(seq(-3,1,0.01)):length(seq(-3,3,0.01))

x1 = x[part1]; x2 = x[part2]; x3 = x[part3]

fit_basis1_part1 = lm(basis[part1, 1] ~ x1 + I(x1^2) + I(x1^3))
fit_basis2_part1 = lm(basis[part1, 2] ~ x1 + I(x1^2) + I(x1^3))
fit_basis3_part1 = lm(basis[part1, 3] ~ x1 + I(x1^2) + I(x1^3))
fit_basis4_part1 = lm(basis[part1, 4] ~ x1 + I(x1^2) + I(x1^3))
fit_basis5_part1 = lm(basis[part1, 5] ~ x1 + I(x1^2) + I(x1^3))

fit_basis1_part2 = lm(basis[part2, 1] ~ x2 + I(x2^2) + I(x2^3))
fit_basis2_part2 = lm(basis[part2, 2] ~ x2 + I(x2^2) + I(x2^3))
fit_basis3_part2 = lm(basis[part2, 3] ~ x2 + I(x2^2) + I(x2^3))
fit_basis4_part2 = lm(basis[part2, 4] ~ x2 + I(x2^2) + I(x2^3))
fit_basis5_part2 = lm(basis[part2, 5] ~ x2 + I(x2^2) + I(x2^3))

fit_basis1_part3 = lm(basis[part3, 1] ~ x3 + I(x3^2) + I(x3^3))
fit_basis2_part3 = lm(basis[part3, 2] ~ x3 + I(x3^2) + I(x3^3))
fit_basis3_part3 = lm(basis[part3, 3] ~ x3 + I(x3^2) + I(x3^3))
fit_basis4_part3 = lm(basis[part3, 4] ~ x3 + I(x3^2) + I(x3^3))
fit_basis5_part3 = lm(basis[part3, 5] ~ x3 + I(x3^2) + I(x3^3))

coefmatrix = rbind(fractions(coef(fit_basis1_part1)),
                   fractions(coef(fit_basis1_part2)),
                   fractions(coef(fit_basis1_part3)),
                   fractions(coef(fit_basis2_part1)),
                   fractions(coef(fit_basis2_part2)),
                   fractions(coef(fit_basis2_part3)),
                   fractions(coef(fit_basis3_part1)),
                   fractions(coef(fit_basis3_part2)),
                   fractions(coef(fit_basis3_part3)),
                   fractions(coef(fit_basis4_part1)),
                   fractions(coef(fit_basis4_part2)),
                   fractions(coef(fit_basis4_part3)),
                   fractions(coef(fit_basis5_part1)),
                   fractions(coef(fit_basis5_part2)),
                   fractions(coef(fit_basis5_part3)))

f_b1_p1 = function(x){
  return(coefmatrix[1,1] + coefmatrix[1,2]*x + coefmatrix[1,3]*x^2 + coefmatrix[1,4]*x^3)
}
f_b1_p2 = function(x){
  return(coefmatrix[2,1] + coefmatrix[2,2]*x + coefmatrix[2,3]*x^2 + coefmatrix[2,4]*x^3)
}
f_b1_p3 = function(x){
  return(coefmatrix[3,1] + coefmatrix[3,2]*x + coefmatrix[3,3]*x^2 + coefmatrix[3,4]*x^3)
}

f_b2_p1 = function(x){
  tmpdif = 3
  return(coefmatrix[1+tmpdif,1] + coefmatrix[1+tmpdif,2]*x + 
           coefmatrix[1+tmpdif,3]*x^2 + coefmatrix[1+tmpdif,4]*x^3)
}
f_b2_p2 = function(x){
  tmpdif = 3
  return(coefmatrix[2+tmpdif,1] + coefmatrix[2+tmpdif,2]*x + 
           coefmatrix[2+tmpdif,3]*x^2 + coefmatrix[2+tmpdif,4]*x^3)
}
f_b2_p3 = function(x){
  tmpdif = 3
  return(coefmatrix[3+tmpdif,1] + coefmatrix[3+tmpdif,2]*x + 
           coefmatrix[3+tmpdif,3]*x^2 + coefmatrix[3+tmpdif,4]*x^3)
}

f_b3_p1 = function(x){
  tmpdif = 6
  return(coefmatrix[1+tmpdif,1] + coefmatrix[1+tmpdif,2]*x + 
           coefmatrix[1+tmpdif,3]*x^2 + coefmatrix[1+tmpdif,4]*x^3)
}
f_b3_p2 = function(x){
  tmpdif = 6
  return(coefmatrix[2+tmpdif,1] + coefmatrix[2+tmpdif,2]*x + 
           coefmatrix[2+tmpdif,3]*x^2 + coefmatrix[2+tmpdif,4]*x^3)
}
f_b3_p3 = function(x){
  tmpdif = 6
  return(coefmatrix[3+tmpdif,1] + coefmatrix[3+tmpdif,2]*x + 
           coefmatrix[3+tmpdif,3]*x^2 + coefmatrix[3+tmpdif,4]*x^3)
}

f_b4_p1 = function(x){
  tmpdif = 9
  return(coefmatrix[1+tmpdif,1] + coefmatrix[1+tmpdif,2]*x + 
           coefmatrix[1+tmpdif,3]*x^2 + coefmatrix[1+tmpdif,4]*x^3)
}
f_b4_p2 = function(x){
  tmpdif = 9
  return(coefmatrix[2+tmpdif,1] + coefmatrix[2+tmpdif,2]*x + 
           coefmatrix[2+tmpdif,3]*x^2 + coefmatrix[2+tmpdif,4]*x^3)
}
f_b4_p3 = function(x){
  tmpdif = 9
  return(coefmatrix[3+tmpdif,1] + coefmatrix[3+tmpdif,2]*x + 
           coefmatrix[3+tmpdif,3]*x^2 + coefmatrix[3+tmpdif,4]*x^3)
}

f_b5_p1 = function(x){
  tmpdif = 12
  return(coefmatrix[1+tmpdif,1] + coefmatrix[1+tmpdif,2]*x + 
           coefmatrix[1+tmpdif,3]*x^2 + coefmatrix[1+tmpdif,4]*x^3)
}
f_b5_p2 = function(x){
  tmpdif = 12
  return(coefmatrix[2+tmpdif,1] + coefmatrix[2+tmpdif,2]*x + 
           coefmatrix[2+tmpdif,3]*x^2 + coefmatrix[2+tmpdif,4]*x^3)
}
f_b5_p3 = function(x){
  tmpdif = 12
  return(coefmatrix[3+tmpdif,1] + coefmatrix[3+tmpdif,2]*x + 
           coefmatrix[3+tmpdif,3]*x^2 + coefmatrix[3+tmpdif,4]*x^3)
}

Bu = function(x){
  if(x < -1){
    tmp = c(f_b1_p1(x), f_b2_p1(x), f_b3_p1((x)), f_b4_p1(x), f_b5_p1(x))
  }else if(x >= -1 & x < 1){
    tmp = c(f_b1_p2(x), f_b2_p2(x), f_b3_p2((x)), f_b4_p2(x), f_b5_p2(x))
  }else{
    tmp = c(f_b1_p3(x), f_b2_p3(x), f_b3_p3((x)), f_b4_p3(x), f_b5_p3(x))
  }
  tmp = c(1, tmp)
  tmp = matrix(tmp, length(tmp), 1)
  return(tmp)
}
