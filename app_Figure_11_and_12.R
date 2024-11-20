library(gplots)
library(nnet)
library(boot)
## Loading the data
data1 <- read.csv("card.csv")
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
Z <- data$nearc4
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


#L-optimal analysis is the same as Y^a is independent of Z | L


#Balke-Pearl bounds((L,A)-optimal)
#a = 0, A = 1 intervention
p.l1 <- function(pi) { pi[7] }
p.l2 <- function(pi) { pi[3] }
p.l3 <- function(pi) { pi[3] + pi[4] - pi[5] - pi[8] }
p.l4 <- function(pi) { pi[2] + pi[3] - pi[5] - pi[6] }
p.l5 <- function(pi) { -5}
p.l6 <- function(pi) { -5  } 
p.l7 <- function(pi) {-5  } 
p.l8 <- function(pi) { -5  }


gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                               p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
arg_gamma.l <- function(pi) {which.max(c(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                         p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)))}

## upper bound functions
#a = 0, A = 1 intervention
p.u1 <- function(pi) { 1 - pi[5] }
p.u2 <- function(pi) { 1 - pi[1] }
p.u3 <- function(pi) { pi[2] + pi[3] + pi[7] + pi[8] }
p.u4 <- function(pi) { pi[3] + pi[4] + pi[6] + pi[7] }
p.u5 <- function(pi) { 5 }
p.u6 <- function(pi) { 5  } 
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

#models for Y and A based on L, Z = 0
fmY.Z0 <- glm(Y~., data = dat[which(Z == 0), -c(1,2,ncol(dat) - 1,ncol(dat))], family = "binomial")
EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
fmY.A0Z0 <- glm(Y~., data = dat[which((1-A)*(1-Z) == 1),-c(1,2,ncol(dat) - 1,ncol(dat))], family = "binomial")
EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
fmY.A1Z0 <- glm(Y~., data = dat[which(A*(1-Z) == 1),-c(1,2,ncol(dat)-1,ncol(dat))], family = "binomial")
EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
fmA.Z0 <- glm(A~., data = dat[which(Z == 0), -c(1,3,ncol(dat) -1,ncol(dat))], family = "binomial")
pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
#models for Y and A based on L, Z = 1
fmY.Z1 <- glm(Y~., data = dat[which(Z == 1), -c(1,2,ncol(dat) - 1,ncol(dat))], family = "binomial")
EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
fmY.A0Z1 <- glm(Y~., data = dat[which((1-A)*Z == 1),-c(1,2,ncol(dat) - 1,ncol(dat))], family = "binomial")
EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
fmY.A1Z1 <- glm(Y~., data = dat[which(A*Z == 1),-c(1,2,ncol(dat)-1,ncol(dat))], family = "binomial")
EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
fmA.Z1 <- glm(A~., data = dat[which(Z == 1), -c(1,3,ncol(dat) -1,ncol(dat))], family = "binomial")
pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")


##Jiang & Ding procedure(superoptimal)
###bootstrap to get variances
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.l1(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)

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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.l2(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)

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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.l3(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.l4(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.l5(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.l6(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.l7(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.l8(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.u1(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.u2(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.u3(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.u4(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.u5(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.u6(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
  
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.u7(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 0, A = 1 intervention | Z = 0
  return((p.u8(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0)

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
#a = 0, A = 1 intervention | Z = 0
BP.bounds.a0A1.Z0 <- c((gamma.l(pi.hat)- EY.A0Z0*(1-pA.Z0))/pA.Z0 - C*sd.l1.hat[arg_gamma.l(pi.hat)], (gamma.u(pi.hat) - EY.A0Z0*(1-pA.Z0))/pA.Z0 + C*sd.u1.hat[arg_gamma.u(pi.hat)])

gam.la0A1 <- gamma.l(pi.hat)

gam.ua0A1 <- gamma.u(pi.hat)

#Balke-Pearl bounds((L,A)-optimal)
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

#models for Y and A based on L, Z = 0
fmY.Z0 <- glm(Y~., data = dat[which(Z == 0), -c(1,2,ncol(dat) - 1,ncol(dat))], family = "binomial")
EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
fmY.A0Z0 <- glm(Y~., data = dat[which((1-A)*(1-Z) == 1),-c(1,2,ncol(dat) - 1,ncol(dat))], family = "binomial")
EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
fmY.A1Z0 <- glm(Y~., data = dat[which(A*(1-Z) == 1),-c(1,2,ncol(dat)-1,ncol(dat))], family = "binomial")
EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
fmA.Z0 <- glm(A~., data = dat[which(Z == 0), -c(1,3,ncol(dat) -1,ncol(dat))], family = "binomial")
pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
#models for Y and A based on L, Z = 1
fmY.Z1 <- glm(Y~., data = dat[which(Z == 1), -c(1,2,ncol(dat) - 1,ncol(dat))], family = "binomial")
EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
fmY.A0Z1 <- glm(Y~., data = dat[which((1-A)*Z == 1),-c(1,2,ncol(dat) - 1,ncol(dat))], family = "binomial")
EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
fmY.A1Z1 <- glm(Y~., data = dat[which(A*Z == 1),-c(1,2,ncol(dat)-1,ncol(dat))], family = "binomial")
EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
fmA.Z1 <- glm(A~., data = dat[which(Z == 1), -c(1,3,ncol(dat) -1,ncol(dat))], family = "binomial")
pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")


##Jiang & Ding procedure(superoptimal)
###bootstrap to get variances
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.l1(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
  
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.l2(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
  
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.l3(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.l4(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.l5(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.l6(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.l7(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.l8(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.u1(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.u2(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.u3(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
  #models for Y and A based on L, Z = 0
  fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
  #models for Y and A based on L, Z = 1
  fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
  pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
  #a = 1, A = 0 intervention / Z = 0
  return((p.u4(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
    #models for Y and A based on L, Z = 0
    fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
    fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
    fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
    EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
    fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
    pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
    #models for Y and A based on L, Z = 1
    fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
    fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
    fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
    EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
    fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
    pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
    #a = 1, A = 0 intervention / Z = 0
    return((p.u5(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
    #models for Y and A based on L, Z = 0
    fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
    fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
    fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
    EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
    fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
    pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
    #models for Y and A based on L, Z = 1
    fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
    fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
    fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
    EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
    fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
    pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
    #a = 1, A = 0 intervention / Z = 0
    return((p.u6(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
    #models for Y and A based on L, Z = 0
    fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
    fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
    fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
    EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
    fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
    pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
    #models for Y and A based on L, Z = 1
    fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
    fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
    fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
    EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
    fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
    pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
    #a = 1, A = 0 intervention / Z = 0
    return((p.u7(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
    #models for Y and A based on L, Z = 0
    fmY.Z0 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 0), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.Z0 <- predict(fmY.Z0, newdata = l, type = "response")
    fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
    fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
    EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
    fmA.Z0 <- glm(A_boot~., data = dat_boot[which(Z_boot == 0), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
    pA.Z0 <- predict(fmA.Z0, newdata = l, type = "response")
    #models for Y and A based on L, Z = 1
    fmY.Z1 <- glm(Y_boot~., data = dat_boot[which(Z_boot == 1), -c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.Z1 <- predict(fmY.Z1, newdata = l, type = "response")
    fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
    EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
    fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot)-1,ncol(dat_boot))], family = "binomial")
    EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
    fmA.Z1 <- glm(A_boot~., data = dat_boot[which(Z_boot == 1), -c(1,3,ncol(dat_boot) -1,ncol(dat_boot))], family = "binomial")
    pA.Z1 <- predict(fmA.Z1, newdata = l, type = "response")
    #a = 1, A = 0 intervention / Z = 0
    return((p.u8(pi.hat) - EY.A0Z0*pA.Z0)/(1-pA.Z0))
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
  #a = 1, A = 0 intervention | Z = 0
  BP.bounds.a1A0.Z0 <- c((gamma.l(pi.hat)- EY.A1Z0*pA.Z0)/(1-pA.Z0) - C*sd.l1.hat[arg_gamma.l(pi.hat)], (gamma.u(pi.hat) - EY.A1Z0*pA.Z0)/(1-pA.Z0) + C*sd.u1.hat[arg_gamma.u(pi.hat)])
  
  gam.la1A0 <- gamma.l(pi.hat)
  
  gam.ua1A0 <- gamma.u(pi.hat)
  

#conditional means estimation
fctA0Z0 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  fmY.A0Z0 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z0 <- predict(fmY.A0Z0, newdata = l, type = "response")
  return(EY.A0Z0)
}
fctA0Z1 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  fmY.A0Z1 <- glm(Y_boot~., data = dat_boot[which((1-A_boot)*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A0Z1 <- predict(fmY.A0Z1, newdata = l, type = "response")
  return(EY.A0Z1)
}
fctA1Z0 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  fmY.A1Z0 <- glm(Y_boot~., data = dat_boot[which(A_boot*(1-Z_boot) == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1Z0 <- predict(fmY.A1Z0, newdata = l, type = "response")
  return(EY.A1Z0)
}
fctA1Z1 <- function(data, ind){
  Y_boot <- data$Y[ind]
  A_boot <- data$A[ind]
  Z_boot <- data$Z[ind]
  L_boot <- data[ind, -c(1,2,3)]
  dat_boot <- data.frame(Z_boot = Z_boot, A_boot = A_boot, Y_boot = Y_boot)
  dat_boot <- cbind(dat_boot, L_boot)
  Y.A_boot <- data$Y.A[ind]
  dat_boot <- cbind(dat_boot, Y.A_boot)
  colnames(dat_boot)[ncol(dat_boot)] <- c("Y.A_boot")
  fmY.A1Z1 <- glm(Y_boot~., data = dat_boot[which(A_boot*Z_boot == 1),-c(1,2,ncol(dat_boot) - 1,ncol(dat_boot))], family = "binomial")
  EY.A1Z1 <- predict(fmY.A1Z1, newdata = l, type = "response")
  return(EY.A1Z1)
}
meanA0Z0 <- EY.A0Z0
meanA1Z0 <- EY.A1Z0
meanA0Z1 <- EY.A0Z1
meanA1Z1 <- EY.A1Z1
set.seed(2023)
sd.A0Z0 <- sd(boot(dat, fctA0Z0, R)$t)
sd.A1Z0 <- sd(boot(dat, fctA1Z0, R)$t)
sd.A0Z1 <- sd(boot(dat, fctA0Z1, R)$t)
sd.A1Z1 <- sd(boot(dat, fctA1Z1, R)$t)
Bounds.A0Z0 <- c(meanA0Z0 - qnorm(0.975)*sd.A0Z0, meanA0Z0 + qnorm(0.975)*sd.A0Z0)
Bounds.A1Z0 <- c(meanA1Z0 - qnorm(0.975)*sd.A1Z0, meanA1Z0 + qnorm(0.975)*sd.A1Z0)
Bounds.A0Z1 <- c(meanA0Z1 - qnorm(0.975)*sd.A0Z1, meanA0Z1 + qnorm(0.975)*sd.A0Z1)
Bounds.A1Z1 <- c(meanA1Z1 - qnorm(0.975)*sd.A1Z1, meanA1Z1 + qnorm(0.975)*sd.A1Z1)


#Decision criteria for Figure 12
##Maximax(A = 0, Z = 0 )
as.numeric(BP.bounds.a1A0.Z0[2] > Bounds.A0Z0[2])
##Maximin(A = 0, Z = 0)
as.numeric(BP.bounds.a1A0.Z0[1] > Bounds.A0Z0[1])
##Maximax(A = 1, Z = 0)
as.numeric(Bounds.A1Z0[2] > BP.bounds.a0A1.Z0[2])
##Maximin(A = 1, Z = 0)
as.numeric(Bounds.A1Z0[1] > BP.bounds.a0A1.Z0[1])





#Combined plots for Figure 11, Z = 0
##Setting line width
width = 3
##Combined plot
par(mar = c(3,4,1,1))
plot(NULL, ylim = c(0,1), xlim = c(0.1,0.9), 
     xlab = "", ylab = "", xaxt = 'n')
axis(1, at = c(0.2, 0.4, 0.6, 0.8), labels = c( "E(Y | A = 0, Z = 0)", "E(Y^1|A = 0, Z = 0)", "E(Y^0 | A = 1, Z = 0)", "E(Y|A = 1, Z = 0)"))
abline(h = 0)
###a = 0, A = 0, Z = 0
arrows(x0=0.2, y0=max(0,Bounds.A0Z0[1]), x1=0.2, 
       y1=min(1,Bounds.A0Z0[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
arrows(x0=0.2, y0=meanA0Z0, x1=0.2, y1=meanA0Z0, 
       code=3, angle=90, length=0.025, lwd=width, col = 'blue')
###a = 1, A = 0, Z = 0
arrows(x0=0.4, y0=max(0,BP.bounds.a1A0.Z0[1]), x1=0.4, 
       y1=min(1,BP.bounds.a1A0.Z0[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
arrows(x0=0.4, y0=max(0,(gam.la1A0 - EY.A1Z0*pA.Z0)/(1-pA.Z0)), x1=0.4, y1=min(1,(gam.ua1A0- EY.A1Z0*pA.Z0)/(1-pA.Z0)), 
       code=3, angle=90, length=0.025, lwd=width, col = 'blue')
###a = 0, A = 1, Z = 0
arrows(x0=0.6, y0=max(0,BP.bounds.a0A1.Z0[1]), x1=0.6, 
       y1=min(1,BP.bounds.a0A1.Z0[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
arrows(x0=0.6, y0=(gam.la0A1- EY.A0Z0*(1-pA.Z0))/pA.Z0, x1=0.6, y1=(gam.ua0A1- EY.A0Z0*(1-pA.Z0))/pA.Z0, 
       code=3, angle=90, length=0.025, lwd=width, col = 'blue')
#a = 1, A = 1, Z = 0
arrows(x0=0.8, y0=max(0,Bounds.A1Z0[1]), x1=0.8, 
       y1=min(1,Bounds.A1Z0[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
arrows(x0=0.8, y0=meanA1Z0, x1=0.8, y1=meanA1Z0, 
       code=3, angle=90, length=0.025, lwd=width, col = 'blue')


