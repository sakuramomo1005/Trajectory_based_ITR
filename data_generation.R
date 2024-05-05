
library(dplyr)
library(MASS)
library(splines)
library(lme4)
library(lsa)
library(simml)
library(ggplot2)

#########################################################################################
# Simulate data for a clinical trial comparing two treatment groups.
#########################################################################################

# Function to simulate training and testing data
simulateData <- function(n, mux, sigmax, alpha, D.pbo, D.drg, error, t, p, missing = 'No') {
  
  # Initialize variables
  subj <- X <- W0 <- week <- truegroup.cs <- truegroup.fix <- truegroup.ran <- y.drg <- y.pbo <- c()
  y.drg.fix <- y.pbo.fix <- y.drg.ran <- y.pbo.ran <- group <- c()
  
  # Simulate test data
  for(i in 1:n) {
    xij <- mvrnorm(1, mux, sigmax)
    w <- c(matrix(xij, 1, p) %*% matrix(alpha, p, 1))
    W0 <- c(W0, rep(w, length(t)))
    X <- rbind(X, matrix(xij, length(t), p, byrow = TRUE))

    b.i.pbo <- mvrnorm(1, rep(0,3), D.pbo)
    b.i.drg <- mvrnorm(1, rep(0,3), D.drg)

    epsilon.ij.drg <- rnorm(length(t), 0, error)
    epsilon.ij.pbo <- rnorm(length(t), 0, error)

    temp.pbo <- pbo(t, w)
    temp.drg <- drg(t, w)

    yij.drg.fix <- temp.drg
    yij.pbo.fix <- temp.pbo 
    yij.drg.ran <- temp.drg + cbind(1, t, t^2) %*% b.i.drg 
    yij.pbo.ran <- temp.pbo + cbind(1, t, t^2) %*% b.i.pbo

    yij.pbo <- yij.pbo.ran + epsilon.ij.pbo
    yij.drg <- yij.drg.ran + epsilon.ij.drg

    y.pbo <- c(y.pbo, yij.pbo)
    y.drg <- c(y.drg, yij.drg)
    y.pbo.fix <- c(y.pbo.fix, yij.pbo.fix)
    y.pbo.ran <- c(y.pbo.ran, yij.pbo.ran)
    y.drg.fix <- c(y.drg.fix, yij.drg.fix)
    y.drg.ran <- c(y.drg.ran, yij.drg.ran)

    # Determine true groups
    truegroup.fix = c(truegroup.fix, rep(if ((yij.pbo.fix[length(t)] - yij.pbo.fix[1]) < (yij.drg.fix[length(t)] - yij.drg.fix[1])) 'pbo' else 'drg', length(t)))
    truegroup.ran = c(truegroup.ran, rep(if ((yij.pbo.ran[length(t)] - yij.pbo.ran[1]) < (yij.drg.ran[length(t)] - yij.drg.ran[1])) 'pbo' else 'drg', length(t)))
    truegroup.cs = c(truegroup.cs, rep(if ((yij.pbo[length(t)] - yij.pbo[1]) < (yij.drg[length(t)] - yij.drg[1])) 'pbo' else 'drg', length(t)))
    
    group <- c(group, rep(sample(c('drg','pbo'))[1], length(t)))
    subj <- c(subj, rep(i, length(t)))
    week <- c(week, t)
  }

  # Create data frame
  data <- data.frame(subj = subj, week = week, W0 = W0, X = X, truegroup.cs = truegroup.cs,
                          truegroup.fix = truegroup.fix, truegroup.ran = truegroup.ran, group = group,
                          y.drg = y.drg, y.pbo = y.pbo, y.drg.ran = y.drg.ran, y.pbo.ran = y.pbo.ran,
                          y.drg.fix = y.drg.fix, y.pbo.fix = y.pbo.fix)
  data$Y <- ifelse(data$group == 'pbo', data$y.pbo, data$y.drg)

  if(missing == 'MCAR'){
    
    mis = c()
    for(ii in 1:(n)){
      mis = c(mis, 1, 1, sample(c(0,1),6, replace = TRUE))
    }
    mis = ifelse(mis == 0, NA, mis)
    data$mis = mis
    data = na.omit(data)
    
  }else if(missing == 'Dropout'){
    
    mistemp =sample(8:4, size = n, replace = TRUE,
                    prob = c(0.5,0.3,0.1,0.05,0.05))
    mis = c()
    for(i in 1:(n)){
      if(mistemp[i] == 8){
        temp = rep(1, 8)
      }else if(mistemp[i] == 7){
        temp = c(rep(1,7),0)
      }else if(mistemp[i] == 6){
        temp = c(rep(1,6), rep(0,2))
      }else if(mistemp[i] == 5){
        temp = c(rep(1,5), rep(0,3))
      }else if(mistemp[i] == 4){
        temp = c(rep(1,4), rep(0,4))
      }
      mis = c(mis, temp)
    }
    mis = ifelse(mis == 0, NA, mis)
    data$mis = mis
    data = na.omit(data)
    
  }
  
  return(data)
}

# Function to compute change scores
computeChangeScores <- function(data) {
  tmp1 <- data.frame(data %>% group_by(subj) %>% slice(1) %>% mutate(yfirst = Y))
  tmp2 <- data.frame(data %>% group_by(subj) %>% slice(n()) %>% mutate(ylast = Y, totalweek = week))
  tmp <- data.frame(subj = tmp1$subj, cs = (tmp2$ylast - tmp1$yfirst) / tmp2$totalweek)
  unique_data <- unique(data[, c('subj', 'truegroup.cs', 'group', paste('X.', 1:p, sep=''), 'W0')])
  unique_data <- merge(unique_data, tmp, by = 'subj')
  return(unique_data)
}

# # Simulate training data
# data <- simulateData(n, mux, sigmax, alpha, D.pbo, D.drg, error, t, p)
# cv_train = data
# # Compute change scores
# unidata <- computeChangeScores(data)
# # Visualization of results using ggplot2
# ggplot(data, aes(x = week, y = Y, group = subj)) + geom_line(col = as.numeric(factor(data$group)))
# 
# # Simulate test data
# data.test <- simulateData(ntest, mux, sigmax, alpha, D.pbo, D.drg, error, t, p)
# # Compute change scores
# data.test <- computeChangeScores(data.test)
# # Display results
# head(data.test)
# range(data.test$cs)
