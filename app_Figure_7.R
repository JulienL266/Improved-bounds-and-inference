library(gplots)
library(nnet)
library(boot)
## Loading the data
data1 <- read.csv("~/Documents/Github/Bounds-and-Simulation/card.csv")
set.seed(1)


# Mean values to replace NAs(following https://arxiv.org/pdf/1611.09925.pdf)
data <- data1
for(i in 1:nrow(data1)){
  for(j in 1:ncol(data1)){
    if(is.na(data1[i,j])){
      data[i,j] <- mean(data1[,j], na.rm = TRUE)
    }
  }
}
data1 <- data

#Binning age, education of parents and iq
for( i in 1:nrow(data)){
  if(data1$age[i] <= 29){
    data$age[i] <- "14-19"
  }else if(data1$age[i] <= 34){
    data$age[i] <- "19-24"
  }
  if(data1$fatheduc[i] <= 12){
    data$fatheduc[i] <- "High-school"
  }else if(data1$fatheduc[i] <= 16){
    data$fatheduc[i] <- "Undergraduate"
  }else{
    data$fatheduc[i] <- "Graduate"
  }
  if(data1$motheduc[i] <= 12){
    data$motheduc[i] <- "High-school"
  }else if(data1$motheduc[i] <= 16){
    data$motheduc[i] <- "Undergraduate"
  }else{
    data$motheduc[i] <- "Graduate"
  }
  if(data1$iq[i] < 70){
    data$iq[i] <- "50-69"
  } 
  if(data1$iq[i] < 90  & data1$iq[i] >= 70){
    data$iq[i] <- "70-89"
  } 
  if(data1$iq[i] < 110 & data1$iq[i] >= 90){
    data$iq[i] <- "90-109"
  } 
  if(data1$iq[i] < 130 & data1$iq[i] >= 110){
    data$iq[i] <- "110-129"
  }
  if(data1$iq[i] < 150 & data1$iq[i] >= 130){
    data$iq[i] <- "130-149"
  }
}
data$age <- as.factor(data$age)
data$fatheduc <- as.factor(data$fatheduc)
data$motheduc <- as.factor(data$motheduc)
data$iq <- as.factor(data$iq)

## Creating all the needed variables
data$Y <- rep(NA, 3010)
cutoff <- median(data$lwage)
for(i in 1:3010){
  if(data$lwage[i] > cutoff){
    data$Y[i] <- 1
  }else{
    data$Y[i] <- 0
  }
}
Y <- data$Y
Z <- as.factor(data$nearc4)
L <- data[,c("age", "black","fatheduc","motheduc","south","smsa","iq")]
# Transforming A into an indicator of education after high-school
A <- data$educ
for(i in 1:length(A)){
  if(A[i] > 12){
    A[i] <- 1
  }else{
    A[i] <- 0
  }
}
n <- length(A)
dat <- data.frame(Z = Z, A = A, Y = Y)
dat <- cbind(dat,L)

l <- L[1,] 
l[1,] <- c("14-19", 0,"High-school","High-school",0,1,"70-89")
l$black <- as.integer(l$black)
l$south <- as.integer(l$south)
l$smsa <- as.integer(l$smsa)


##Uncomment upper and lower bound fcts depending on the estimand of interest
## lower bound functions
#ATE
#p.l1 <- function(pi) { pi[8] + pi[1] - 1 }
#p.l2 <- function(pi) { pi[4] + pi[5] - 1 }
#p.l3 <- function(pi) { -pi[6] - pi[7] }
#p.l4 <- function(pi) { -pi[2] - pi[3] }
#p.l5 <- function(pi) { pi[4] - pi[8] - pi[7] - pi[2] - pi[3] }
#p.l6 <- function(pi) { pi[8] - pi[4] - pi[3] - pi[6] - pi[7] }
#p.l7 <- function(pi) { pi[5] - pi[6] - pi[7] - pi[2] - pi[1] }
#p.l8 <- function(pi) { pi[1] - pi[2] - pi[3] - pi[6] - pi[5] }
#a = 0 intervention
#p.l1 <- function(pi) { pi[7] }
#p.l2 <- function(pi) { pi[3] }
#p.l3 <- function(pi) { pi[3] + pi[4] - pi[5] - pi[8] }
#p.l4 <- function(pi) { pi[2] + pi[3] - pi[5] - pi[6] }
#p.l6 <- function(pi) { -5  } 
#p.l5 <- function(pi) { -5}
#p.l7 <- function(pi) {-5  } 
#p.l8 <- function(pi) { -5  }

#a = 1 intervention
p.l1 <- function(pi) { pi[4] }
p.l2 <- function(pi) { pi[8] }
p.l3 <- function(pi) { pi[5] + pi[8] - pi[1] - pi[2] }
p.l4 <- function(pi) { pi[7] + pi[8] - pi[2] - pi[3] }
p.l5 <- function(pi) { -5 }
p.l6 <- function(pi) { -5  } 
p.l7 <- function(pi) {-5  } 
p.l8 <- function(pi) { -5  }

gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                               p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
arg_gamma.l <- function(pi) {which.max(c(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                         p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)))}

## upper bound functions
#ATE
#p.u1 <- function(pi) { 1 - pi[6] - pi[3] }
#p.u2 <- function(pi) { 1 - pi[2] - pi[7] }
#p.u3 <- function(pi) { pi[8] + pi[5] }
#p.u4 <- function(pi) { pi[4] + pi[1] }
#p.u5 <- function(pi) { -pi[2] + pi[6] + pi[5] + pi[4] + pi[1] }
#p.u6 <- function(pi) { -pi[6] + pi[2] + pi[1] + pi[8] + pi[5] }
#p.u7 <- function(pi) { -pi[7] + pi[8] + pi[5] + pi[4] + pi[3] }
#p.u8 <- function(pi) { -pi[3] + pi[4] + pi[1] + pi[8] + pi[7] }
#a = 0 intervention
#p.u1 <- function(pi) { 1 - pi[5] }
#p.u2 <- function(pi) { 1 - pi[1] }
#p.u3 <- function(pi) { pi[2] + pi[3] + pi[7] + pi[8] }
#p.u4 <- function(pi) { pi[3] + pi[4] + pi[6] + pi[7] }
#p.u5 <- function(pi) { 5 }
#p.u6 <- function(pi) { 5  } 
#p.u7 <- function(pi) { 5 } 
#p.u8 <- function(pi) { 5 }

#a = 1 intervention
p.u1 <- function(pi) { 1 - pi[6] }
p.u2 <- function(pi) { 1 - pi[2] }
p.u3 <- function(pi) { pi[1] + pi[4] + pi[7] + pi[8] }
p.u4 <- function(pi) { pi[3] + pi[4] + pi[5] + pi[8] }
p.u5 <- function(pi) { 5  }
p.u6 <- function(pi) { 5  } 
p.u7 <- function(pi) { 5 } 
p.u8 <- function(pi) { 5 }

gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                               p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
arg_gamma.u <- function(pi) {which.min(c(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                         p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)))}


dat$Y.A <- factor(dat$Y):factor(dat$A)
pi.ya.0 <- multinom(Y.A ~ ., data = dat[dat$Z==0, -(1:3)], family = "binomial")
pi.ya.1 <- multinom(Y.A ~ ., data = dat[dat$Z==1, -(1:3)], family = "binomial")

pi.hat <- data.frame(predict(pi.ya.0, type = 'probs', newdata = l))
pi.hat <- data.frame(t(pi.hat))
colnames(pi.hat) <-
  c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
## compile pi.hat_{ya.1} estimates
colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
  c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
pi.hat <- as.numeric(pi.hat)
round(c(gamma.l(pi.hat), gamma.u(pi.hat)), 3)

##Jiang & Ding procedure(optimal)
###bootstrap to get variances
sd.l1.hat <- rep(NA, 8)
sd.u1.hat <- rep(NA, 8)
data = data.frame(Y = Y, A = A, Z = Z)
data <- cbind(data,L)
fct.l1 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.l1(pi.hat))
}
fct.l2 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.l2(pi.hat))
}
fct.l3 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.l3(pi.hat))
}



fct.l4 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.l4(pi.hat))
}

fct.l5 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.l5(pi.hat))
}

fct.l6 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.l6(pi.hat))
}

fct.l7 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.l7(pi.hat))
}

fct.l8 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.l8(pi.hat))
}
fct.u1 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.u1(pi.hat))
}
fct.u2 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.u2(pi.hat))
}
fct.u3 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.u3(pi.hat))
}
fct.u4 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.u4(pi.hat))
}
fct.u5 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.u5(pi.hat))
}
fct.u6 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.u6(pi.hat))
}
fct.u7 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.u7(pi.hat))
}
fct.u8 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3,12)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,11)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,11)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  return(p.u8(pi.hat))
}
set.seed(2023)
R <- 1000
sd.l1.hat[1] <- sd(boot(dat, fct.l1, R)$t)
sd.l1.hat[2] <- sd(boot(dat, fct.l2, R)$t)
sd.l1.hat[3] <- sd(boot(dat, fct.l3, R)$t)
sd.l1.hat[4] <- sd(boot(dat, fct.l4, R)$t)
sd.l1.hat[5] <- sd(boot(dat, fct.l5, R)$t)
sd.l1.hat[6] <- sd(boot(dat, fct.l6, R)$t)
sd.l1.hat[7] <- sd(boot(dat, fct.l7, R)$t)
sd.l1.hat[8] <- sd(boot(dat, fct.l8, R)$t)

sd.u1.hat[1] <- sd(boot(dat, fct.u1, R)$t)
sd.u1.hat[2] <- sd(boot(dat, fct.u2, R)$t)
sd.u1.hat[3] <- sd(boot(dat, fct.u3, R)$t)
sd.u1.hat[4] <- sd(boot(dat, fct.u4, R)$t)
sd.u1.hat[5] <- sd(boot(dat, fct.u5, R)$t)
sd.u1.hat[6] <- sd(boot(dat, fct.u6, R)$t)
sd.u1.hat[7] <- sd(boot(dat, fct.u7, R)$t)
sd.u1.hat[8] <- sd(boot(dat, fct.u8, R)$t)

alpha = 0.05
f <- function(x){
  return(pnorm(x + (gamma.u(pi.hat) - gamma.l(pi.hat))/(max(sd.l1.hat[arg_gamma.l(pi.hat)], sd.u1.hat[arg_gamma.u(pi.hat)]))) - pnorm(-x) - (1-alpha))
}
C <- uniroot(f, interval = c(-3, 3), extendInt = "yes")
C <- C$root
BP.bounds <- c(gamma.l(pi.hat) - C*sd.l1.hat[arg_gamma.l(pi.hat)], gamma.u(pi.hat) + C*sd.u1.hat[arg_gamma.u(pi.hat)])

#Uncomment quantities below depending on the estimate of interest

#BP.bounds.ATE <- BP.bounds
#BP.bounds0 <- BP.bounds
BP.bounds1 <- BP.bounds

#gam.lATE <- gamma.l(pi.hat)
#gam.l0 <- gamma.l(pi.hat)
gam.l1 <- gamma.l(pi.hat)

#gam.uATE <- gamma.u(pi.hat)
#gam.u0 <- gamma.u(pi.hat)
gam.u1 <- gamma.u(pi.hat)

##Plot(uncomment different arrow commands depending on the estimand)
par(mar = c(3,4,1,1))
plot(NULL, ylim = c(0,1), xlim = c(0.2,0.8), 
     xlab = "", ylab = "", xaxt = 'n')
axis(1, at = c(0.35, 0.65), labels = c("a = 0", "a = 1"))
abline(h = 0)
#arrows(x0=0.5, y0=BP.bounds[1], x1=0.5, 
#     y1= BP.bounds[2], 
#  code=3, angle=90, length=0.05, lwd=2, col = 'red')
#arrows(x0=0.5, y0=gamma.l(pi.hat), x1=0.5, y1=gamma.u(pi.hat), 
#     code=3, angle=90, length=0.025, lwd=2, col = 'blue')
arrows(x0=0.35, y0=BP.bounds0[1], x1=0.35, 
       y1=BP.bounds0[2], 
       code=3, angle=90, length=0.05, lwd=2, col = 'red')
arrows(x0=0.35, y0=gam.l0, x1=0.35, y1=gam.u0, 
       code=3, angle=90, length=0.025, lwd=2, col = 'blue')
arrows(x0=0.65, y0=BP.bounds1[1], x1=0.65, 
       y1=BP.bounds1[2], 
       code=3, angle=90, length=0.05, lwd=2, col = 'red')
arrows(x0=0.65, y0=gam.l1, x1=0.65, y1=gam.u1, 
       code=3, angle=90, length=0.025, lwd=2, col = 'blue')

#Decision criteria(optimal)
##Maximax
as.numeric(BP.bounds1[2] > BP.bounds0[2])
##Maximin
as.numeric(BP.bounds1[1] > BP.bounds0[1])
##Minimax
as.numeric(BP.bounds.ATE[1] > 0 || (BP.bounds.ATE[2] > 0 && abs(BP.bounds.ATE[2]) > abs(BP.bounds.ATE[1])))
#Healthcare
as.numeric(BP.bounds.ATE[1] > 0)



#Balke-Pearl bounds(superoptimal)
#Uncomment upper and lower bound fcts depending on the estimand of interest
#CATE(A = 0)
#p.l1 <- function(pi) { pi[4] }
#p.l2 <- function(pi) { pi[8] }
#p.l3 <- function(pi) { pi[5] + pi[8] - pi[1] - pi[2] }
#p.l4 <- function(pi) { pi[7] + pi[8] - pi[2] - pi[3] }
#p.l5 <- function(pi) { -5 }
#p.l6 <- function(pi) { -5  } 
#p.l7 <- function(pi) { -5  }
#p.l8 <- function(pi) { -5  }

#CATE(A = 1)
#p.l1 <- function(pi) { pi[7] }
#p.l2 <- function(pi) { pi[3] }
#p.l3 <- function(pi) { pi[3] + pi[4] - pi[5] - pi[8] }
#p.l4 <- function(pi) { pi[2] + pi[3] - pi[5] - pi[6] }
#p.l5 <- function(pi) { -5}
#p.l6 <- function(pi) { -5  } 
#p.l7 <- function(pi) {-5  } 
#p.l8 <- function(pi) { -5  }

#a = 0, A = 1 intervention
#p.l1 <- function(pi) { pi[7] }
#p.l2 <- function(pi) { pi[3] }
#p.l3 <- function(pi) { pi[3] + pi[4] - pi[5] - pi[8] }
#p.l4 <- function(pi) { pi[2] + pi[3] - pi[5] - pi[6] }
#p.l5 <- function(pi) { -5}
#p.l6 <- function(pi) { -5  } 
#p.l7 <- function(pi) {-5  } 
#p.l8 <- function(pi) { -5  }

#a = 1, A = 0 intervention
p.l1 <- function(pi) { pi[4] }
p.l2 <- function(pi) { pi[8] }
p.l3 <- function(pi) { pi[5] + pi[8] - pi[1] - pi[2] }
p.l4 <- function(pi) { pi[7] + pi[8] - pi[2] - pi[3] }
p.l5 <- function(pi) { -5 }
p.l6 <- function(pi) { -5  } 
p.l7 <- function(pi) { -5  } 
p.l8 <- function(pi) { -5  }


gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                               p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
arg_gamma.l <- function(pi) {which.max(c(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                         p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)))}

## upper bound functions
#CATE(A = 0)
#p.u1 <- function(pi) { 1 - pi[6] }
#p.u2 <- function(pi) { 1 - pi[2] }
#p.u3 <- function(pi) { pi[1] + pi[4] + pi[7] + pi[8] }
#p.u4 <- function(pi) { pi[3] + pi[4] + pi[5] + pi[8] }
#p.u5 <- function(pi) { 5  }
#p.u6 <- function(pi) { 5 } 
#p.u7 <- function(pi) { 5 } 
#p.u8 <- function(pi) { 5}

#CATE(A = 1)
#p.u1 <- function(pi) { 1 - pi[5] }
#p.u2 <- function(pi) { 1 - pi[1] }
#p.u3 <- function(pi) { pi[2] + pi[3] + pi[7] + pi[8] }
#p.u4 <- function(pi) { pi[3] + pi[4] + pi[6] + pi[7] }
#p.u5 <- function(pi) { 5 }
#p.u6 <- function(pi) { 5  } 
#p.u7 <- function(pi) { 5 } 
#p.u8 <- function(pi) { 5 }

#a = 0, A = 1 intervention
#p.u1 <- function(pi) { 1 - pi[5] }
#p.u2 <- function(pi) { 1 - pi[1] }
#p.u3 <- function(pi) { pi[2] + pi[3] + pi[7] + pi[8] }
#p.u4 <- function(pi) { pi[3] + pi[4] + pi[6] + pi[7] }
#p.u5 <- function(pi) { 5 }
#p.u6 <- function(pi) { 5  } 
#p.u7 <- function(pi) { 5 } 
#p.u8 <- function(pi) { 5 }

#a = 1, A = 0 intervention
p.u1 <- function(pi) { 1 - pi[6] }
p.u2 <- function(pi) { 1 - pi[2] }
p.u3 <- function(pi) { pi[1] + pi[4] + pi[7] + pi[8] }
p.u4 <- function(pi) { pi[3] + pi[4] + pi[5] + pi[8] }
p.u5 <- function(pi) { 5  }
p.u6 <- function(pi) { 5 } 
p.u7 <- function(pi) { 5 } 
p.u8 <- function(pi) { 5 }


gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                               p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
arg_gamma.u <- function(pi) {which.min(c(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                         p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)))}

#Balke Pearl bounds(superoptimal)
dat$Y.A <- factor(dat$Y):factor(dat$A)
pi.ya.0 <- multinom(Y.A ~ ., data = dat[dat$Z==0, -(1:3)], family = "binomial")
pi.ya.1 <- multinom(Y.A ~ ., data = dat[dat$Z==1, -(1:3)], family = "binomial")
pi.hat <- data.frame(predict(pi.ya.0, type = 'probs', newdata = l))
pi.hat <- data.frame(t(pi.hat))
colnames(pi.hat) <-
  c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
## compile pi.hat_{ya.1} estimates
colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
  c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
pi.hat <- as.numeric(pi.hat)

#models for Y and A based on L
fmY <- glm(Y~., data = dat[, -c(1,2,ncol(dat))], family = "binomial")
EY <- predict(fmY, newdata = l, type = "response")
fmY.A0 <- glm(Y~., data = dat[which(A == 0),-c(1,2,ncol(dat))], family = "binomial")
EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
fmY.A1 <- glm(Y~., data = dat[which(A == 1),-c(1,2,ncol(dat))], family = "binomial")
EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
fmA <- glm(A~., data = dat[, -c(1,3,ncol(dat))], family = "binomial")
pA <- predict(fmA, newdata = l, type = "response")
##Uncomment following depending on the estimand of interest
#CATE(A = 0)
#round(c((gamma.l(pi.hat) - EY)/(1-pA), (gamma.u(pi.hat) - EY)/(1-pA)), 3)
#CATE(A = 1)
#round(c( (-gamma.u(pi.hat) + EY)/pA, (- gamma.l(pi.hat) + EY)/pA), 3)
#a = 0, A = 1 intervention
#round(c((gamma.l(pi.hat) - EY.A0*(1-pA))/pA, (gamma.u(pi.hat) - EY.A0*(1-pA))/pA),3)
#a = 1, A = 0 intervention
#round(c((gamma.l(pi.hat) - EY.A1*pA)/(1-pA), (gamma.u(pi.hat) - EY.A1*pA)/(1-pA)),3)

##Jiang & Ding procedure(superoptimal)
###bootstrap to get variances
##Uncomment upper and lower bound fcts depending on the estimand of interest
sd.l1.hat <- rep(NA, 8)
sd.u1.hat <- rep(NA, 8)
fct.l1 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.l1(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.l1(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.l1(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.l1(pi.hat) - EY.A0*pA)/(1-pA))
}
fct.l2 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1 )], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.l2(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.l2(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.l2(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.l2(pi.hat) - EY.A0*pA)/(1-pA))
}
fct.l3 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.l3(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.l3(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.l3(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.l3(pi.hat) - EY.A0*pA)/(1-pA))
}

fct.l4 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) -1,ncol(dat_boot) )], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.l4(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.l4(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.l4(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.l4(pi.hat) - EY.A0*pA)/(1-pA))
}


fct.l5 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.l5(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.l5(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.l5(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.l5(pi.hat) - EY.A0*pA)/(1-pA))
}

fct.l6 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.l6(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.l6(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.l6(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.l6(pi.hat) - EY.A0*pA)/(1-pA))
}

fct.l7 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) -1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) -1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.l7(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.l7(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.l7(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.l7(pi.hat) - EY.A0*pA)/(1-pA))
}

fct.l8 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.l8(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.l8(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.l8(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.l8(pi.hat) - EY.A0*pA)/(1-pA))
}
fct.u1 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.u1(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.u1(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.u1(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.u1(pi.hat) - EY.A0*pA)/(1-pA))
}
fct.u2 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.u2(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.u2(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.u2(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.u2(pi.hat) - EY.A0*pA)/(1-pA))
}
fct.u3 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.u3(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.u3(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.u3(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.u3(pi.hat) - EY.A0*pA)/(1-pA))
}
fct.u4 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.u4(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.u4(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.u4(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.u4(pi.hat) - EY.A0*pA)/(1-pA))
}
fct.u5 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.u5(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.u5(pi.hat) - EY)/pA)
  #a = 0, A = 1 intervention
  #return((p.u5(pi.hat) - EY.A0*(1-pA))/pA)
  #a = 1, A = 0 intervention
  return((p.u5(pi.hat) - EY.A0*pA)/(1-pA))
}
fct.u6 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot) )], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.u6(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.u6(pi.hat) - EY)/pA)
}
fct.u7 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.u7(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.u7(pi.hat) - EY)/pA)
}
fct.u8 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  pi.ya.0 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==0, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  pi.ya.1 <- multinom(Y.A_boot ~ ., data = dat_boot[dat_boot$Z_boot==1, -c(1,2,3,ncol(dat_boot) - 1)], family = "binomial")
  
  pi.hat <- data.frame(t(predict(pi.ya.0, type = 'probs', newdata = l)))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  ## compile pi.hat_{ya.1} estimates
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  #models for Y and A based on L
  fmY <- glm(Y_boot~., data = dat_boot[, -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y_boot~., data = dat_boot[which(A_boot == 0),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y_boot~., data = dat_boot[which(A_boot == 1),-c(1,2,ncol(dat_boot) -1,ncol(dat_boot) )], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A_boot~., data = dat_boot[, -c(1,3,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  #CATE(A=0)
  #return((p.u8(pi.hat) - EY)/(1-pA))
  #CATE(A = 1)
  #return((p.u8(pi.hat) - EY)/pA)
  return((p.u8(pi.hat) - EY.A0*pA)/(1-pA))
}
set.seed(2023)
R <- 1000
sd.l1.hat[1] <- sd(boot(dat, fct.l1, R)$t)
sd.l1.hat[2] <- sd(boot(dat, fct.l2, R)$t)
sd.l1.hat[3] <- sd(boot(dat, fct.l3, R)$t)
sd.l1.hat[4] <- sd(boot(dat, fct.l4, R)$t)
sd.l1.hat[5] <- sd(boot(dat, fct.l5, R)$t)
sd.l1.hat[6] <- sd(boot(dat, fct.l6, R)$t)
sd.l1.hat[7] <- sd(boot(dat, fct.l7, R)$t)
sd.l1.hat[8] <- sd(boot(dat, fct.l8, R)$t)

sd.u1.hat[1] <- sd(boot(dat, fct.u1, R)$t)
sd.u1.hat[2] <- sd(boot(dat, fct.u2, R)$t)
sd.u1.hat[3] <- sd(boot(dat, fct.u3, R)$t)
sd.u1.hat[4] <- sd(boot(dat, fct.u4, R)$t)
sd.u1.hat[5] <- sd(boot(dat, fct.u5, R)$t)
sd.u1.hat[6] <- sd(boot(dat, fct.u6, R)$t)
sd.u1.hat[7] <- sd(boot(dat, fct.u7, R)$t)
sd.u1.hat[8] <- sd(boot(dat, fct.u8, R)$t)

alpha = 0.05
f <- function(x){
  return(pnorm(x + (gamma.u(pi.hat) - gamma.l(pi.hat))/(max(sd.l1.hat[arg_gamma.l(pi.hat)], sd.u1.hat[arg_gamma.u(pi.hat)]))) - pnorm(-x) - (1-alpha))
}
C <- uniroot(f, interval = c(-3, 3), extendInt = "yes")
C <- C$root


#Uncomment following depending on estimand of interest
#CATE(A = 0)
#BP.bounds <- c((gamma.l(pi.hat)- EY)/(1-pA) - C*sd.l1.hat[arg_gamma.l(pi.hat)], (gamma.u(pi.hat) - EY)/(1-pA) + C*sd.u1.hat[arg_gamma.u(pi.hat)])
#CATE(A = 1)
#BP.bounds.CATE1 <- c( (-gamma.u(pi.hat) + EY)/pA - C*sd.u1.hat[arg_gamma.u(pi.hat)], (-gamma.l(pi.hat)+ EY)/pA + C*sd.l1.hat[arg_gamma.l(pi.hat)])


#BP.bounds.CATE0 <- BP.bounds



#gam.lCATE1 <- gamma.l(pi.hat)
#gam.lCATE0 <- gamma.l(pi.hat)



#gam.uCATE1 <- gamma.u(pi.hat)
#gam.uCATE0 <- gamma.u(pi.hat)









#Combined CATE plots
##Setting line width
width = 3
##Combined plot
par(mar = c(4,5,1,1), cex.axis = 1.5, cex.lab = 1.5)
plot(NULL, ylim = c(-1,1), xlim = c(0.2,0.8), 
     xlab = "", ylab = "", xaxt = 'n')
axis(1, at = c(0.35, 0.5, 0.65), labels = c("l-CATE", "(l,0)-CATE", "(l,1)-CATE"))
abline(h = 0)
###L-CATE
arrows(x0=0.35, y0=BP.bounds.ATE[1], x1=0.35, 
       y1= BP.bounds.ATE[2], 
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
arrows(x0=0.35, y0=gam.lATE, x1=0.35, y1=gam.uATE, 
       code=3, angle=90, length=0.025, lwd=width, col = 'blue')
###A = 0
arrows(x0=0.5, y0=BP.bounds.CATE0[1], x1=0.5, 
       y1=BP.bounds.CATE0[2], 
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
arrows(x0=0.5, y0=(gam.lCATE0 - EY)/(1-pA), x1=0.5, y1=(gam.uCATE0 - EY)/(1-pA), 
       code=3, angle=90, length=0.025, lwd=width, col = 'blue')
###A = 1
arrows(x0=0.65, y0=BP.bounds.CATE1[1], x1=0.65, 
       y1=BP.bounds.CATE1[2], 
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
arrows(x0=0.65, y0=(-gam.uCATE1 + EY)/pA, x1=0.65, y1=(-gam.lCATE1 + EY)/pA, 
       code=3, angle=90, length=0.025, lwd=width, col = 'blue')
