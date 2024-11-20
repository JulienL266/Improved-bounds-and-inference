library(rpart)
library(randomForest)
library(nnet)


#Card example (remove '#'s if you want to do the analysis with the NLSYM dataset)
data1 <- read.csv("card.csv")


# Mean values to replace NAs (following https://arxiv.org/pdf/1611.09925.pdf)
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
}else {
 data$fatheduc[i] <- "Higher education"
}
if(data1$motheduc[i] <= 12){
 data$motheduc[i] <- "High-school"
}else {
 data$motheduc[i] <- "Higher education"
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
# Transforming A into a variable about education after high-school
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

# Data analysis
set.seed(532)
dat$Y.A <- factor(dat$Y):factor(dat$A)
sample.split <- sample(1:n, floor(0.2*n), replace = FALSE)
dat.t <- dat[sample.split, ]
dat.te <- dat[-sample.split, ]
M <- 10 ## number of data splits for cross-fitting
test.indices <- list()
remaining <- 1:nrow(dat.te)
for (m in 1:(M-1)) {
  test.indices[[m]] <- sample(remaining, floor(nrow(dat.te)/M), replace = F)
  remaining <- remaining[! (remaining %in% test.indices[[m]])]
}
test.indices[[M]] <- remaining

bound_regime <- function(g, dat.train, dat.test){
  #Bounds are estimated using Levis procedure, but the bounds used to find the 
  #decision criteria are found in the same way as in the Card example
  hist.incr.l <- c()
  hist.incr.u <- c()
  hist.l <- c(0)
  hist.u <- c(0)
  lower_bound <- 0
  upper_bound <- 0
  #define A,Z and L for below from dat.test
  n <- nrow(dat.test)
  A <- dat.test$A
  Z <- dat.test$Z
  Y <- dat.test$Y
  L <- dat.test[,-c(1,2,3,ncol(dat.test))]
  #models needed
  ## Y models
  EY.gA0Z0 <- rpart(as.factor(Y)~., data = dat.train[which((1-dat.train$Z)*(1-dat.train$A) == 1), -c(1,2,ncol(dat.train))])
  EY.gA0Z1 <- rpart(as.factor(Y)~., data = dat.train[which(dat.train$Z*(1-dat.train$A) == 1), -c(1,2,ncol(dat.train))])
  EY.gA1Z0 <- rpart(as.factor(Y)~., data = dat.train[which((1-dat.train$Z)*dat.train$A == 1), -c(1,2,ncol(dat.train))])
  EY.gA1Z1 <- rpart(as.factor(Y)~., data = dat.train[which(dat.train$Z*dat.train$A == 1), -c(1,2,ncol(dat.train))])
  ## A:Z models
  A.Z <- factor(dat.train$A):factor(dat.train$Z)
  pAZ <- randomForest(A.Z ~., data = dat.train[, -c(1,2,3,ncol(dat.train))])
  ## Y:A models
  pYA.gZ0 <- rpart(Y.A~., data = dat.train[dat.train$Z == 0, -(1:3)])
  pYA.gZ1 <- rpart(Y.A~., data = dat.train[dat.train$Z == 1, -(1:3)])
  ## Z model
  lambda <- glm(Z~., data = dat.train[, -c(2,3,ncol(dat.train))], family = "binomial") #correctly specified if we only have categorical data
  ## A models
  pAZ0 <- rpart(as.factor(A)~., data = dat.train[dat.train$Z == 0, -c(1,3,ncol(dat.train))])
  pAZ1 <- rpart(as.factor(A)~., data = dat.train[dat.train$Z == 1, -c(1,3,ncol(dat.train))])
  pZ <- rpart(as.factor(Z)~., data = dat.train[, -c(2,3,ncol(dat.train))])
  for(i in 1:nrow(dat.test)){
    #model fits
    pAZ.hat <- predict(pAZ, type = "prob", newdata = L[i,]) #00,01,10,11
    pYA.gZ0.hat <- predict(pYA.gZ0, type = "prob", newdata = L[i,])
    pYA.gZ1.hat <- predict(pYA.gZ1, type = "prob", newdata = L[i,])
    pi.hats <- cbind(pYA.gZ0.hat, pYA.gZ1.hat)
    lambda.1 <- predict(lambda, type = "response", newdata = L[i,])
    pAZ0.hat <- predict(pAZ0, type = "prob", newdata = L[i,])
    pAZ1.hat <- predict(pAZ1, type = "prob", newdata = L[i,])
    EY.gA0Z0.hat <- predict(EY.gA0Z0, type = "prob", newdata = L[i,])
    EY.gA0Z1.hat <- predict(EY.gA0Z1, type = "prob", newdata = L[i,])
    EY.gA1Z0.hat <- predict(EY.gA1Z0, type = "prob", newdata = L[i,])
    EY.gA1Z1.hat <- predict(EY.gA0Z0, type = "prob", newdata = L[i,])
    pZ.hat <- predict(pZ, type = "prob", newdata = L[i,])
    if(g(A[i], L[i,], Z[i]) == 0){
      if(A[i] == 0){
        #E(Y | A = 0, L = l, Z = z) bounds
        if(Z[i] == 0){
          EY.gA0Z0.hat <- predict(EY.gA0Z0, type = "prob", newdata = L[i,])[,2]
          #lower_bound = lower_bound + max(0,min(1,(EY.gA0Z0.hat + Y[i] - EY.gA0Z0.hat)))
          #upper_bound = upper_bound + max(0,min(1,(EY.gA0Z0.hat + Y[i] - EY.gA0Z0.hat)))
          lower_bound = lower_bound + (EY.gA0Z0.hat + Y[i] - EY.gA0Z0.hat)
          upper_bound = upper_bound + (EY.gA0Z0.hat + Y[i] - EY.gA0Z0.hat)
        }else{
          EY.gA0Z1.hat <- predict(EY.gA0Z1, type = "prob", newdata = L[i,])[,2]
          #lower_bound = lower_bound + max(0,min(1,(EY.gA0Z1.hat + Y[i] - EY.gA0Z1.hat)))
          #upper_bound = upper_bound + max(0,min(1,(EY.gA0Z1.hat + Y[i] - EY.gA0Z1.hat)))
          lower_bound = lower_bound + (EY.gA0Z1.hat + Y[i] - EY.gA0Z1.hat)
          upper_bound = upper_bound + (EY.gA0Z1.hat + Y[i] - EY.gA0Z1.hat)
          
        }
      }else{
        #E(Y^0 | A = 1, L = l, Z = z) bounds
        if(Z[i] == 0){
          ## bound fcts
          p.l1 <- function(pi) { pi[7] }
          p.l2 <- function(pi) { pi[3] }
          p.l3 <- function(pi) { pi[3] + pi[4] - pi[5] - pi[8] }
          p.l4 <- function(pi) { pi[2] + pi[3] - pi[5] - pi[6] }
          p.l5 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) }
          p.l6 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
          p.l7 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi))  } 
          p.l8 <- function(pi) {pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi))  }
          gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                         p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
          p.u1 <- function(pi) { 1 - pi[5] }
          p.u2 <- function(pi) { 1 - pi[1] }
          p.u3 <- function(pi) { pi[2] + pi[3] + pi[7] + pi[8] }
          p.u4 <- function(pi) { pi[3] + pi[4] + pi[6] + pi[7] }
          p.u5 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 }
          p.u6 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 } 
          p.u7 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 } 
          p.u8 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 }
          gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                         p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
          ## constructing psis and pis
          psi_00.0 <- pYA.gZ0.hat[1] + 
            (1 - Z[i]) * ((1 - Y[i]) * (1 - A[i]) - pYA.gZ0.hat[1]) / (1 - lambda.1)
          psi_01.0 <- pYA.gZ0.hat[2] + 
            (1 - Z[i]) * ((1 - Y[i]) * A[i] - pYA.gZ0.hat[2]) / (1 - lambda.1)
          psi_10.0 <- pYA.gZ0.hat[3] + 
            (1 - Z[i]) * (Y[i] * (1 - A[i]) - pYA.gZ0.hat[3]) / (1 - lambda.1)
          psi_11.0 <- pYA.gZ0.hat[4] + 
            (1 - Z[i]) * (Y[i] * A[i] - pYA.gZ0.hat[4]) / (1 - lambda.1)
          psi_00.1 <- pYA.gZ1.hat[1] + 
            Z[i] * ((1 - Y[i]) * (1 - A[i]) - pYA.gZ1.hat[1]) / lambda.1
          psi_01.1 <- pYA.gZ1.hat[2] + 
            Z[i] * ((1 - Y[i]) * A - pYA.gZ1.hat[2]) / lambda.1
          psi_10.1 <- pYA.gZ1.hat[3] + 
            Z[i] * (Y[i] * (1 - A[i]) - pYA.gZ1.hat[3]) / lambda.1
          psi_11.1 <- pYA.gZ1.hat[4] + 
            Z[i] * (Y[i] * A[i] - pYA.gZ1.hat[4]) / lambda.1
          psi.hats <- c(psi_00.0, psi_01.0, psi_10.0, psi_11.0, 
                        psi_00.1, psi_01.1, psi_10.1, psi_11.1)
          
          ## computing L's and p's
          p.l.hats <- cbind(p.l1(pi.hats)[1], p.l2(pi.hats)[1], p.l3(pi.hats)[1],
                            p.l4(pi.hats)[1], p.l5(pi.hats)[1], p.l6(pi.hats)[1],
                            p.l7(pi.hats)[1], p.l8(pi.hats)[1])
          p.u.hats <- cbind(p.u1(pi.hats)[1], p.u2(pi.hats)[1], p.u3(pi.hats)[1],
                            p.u4(pi.hats)[1], p.u5(pi.hats)[1], p.u6(pi.hats)[1],
                            p.u7(pi.hats)[1], p.u8(pi.hats)[1])
          L.hats <- cbind(p.l1(psi.hats)[1], p.l2(psi.hats)[1], p.l3(psi.hats)[1],
                          p.l4(psi.hats)[1], p.l5(psi.hats)[1], p.l6(psi.hats)[1],
                          p.l7(psi.hats)[1], p.l8(psi.hats)[1])
          U.hats <- cbind(p.u1(psi.hats)[1], p.u2(psi.hats)[1], p.u3(psi.hats)[1],
                          p.u4(psi.hats)[1], p.u5(psi.hats)[1], p.u6(psi.hats)[1],
                          p.u7(psi.hats)[1], p.u8(psi.hats)[1])
          
          ## influence fcts
          argmax.p.l <- apply(p.l.hats, MARGIN = 1, FUN = which.max)
          
          argmin.p.u <- apply(p.u.hats, MARGIN = 1, FUN = which.min)
          
          IF.l <- L.hats[argmax.p.l]
          IF.u <- U.hats[argmin.p.u]
          
          #centered influence functions
          pi_1.l.hat <- p.l.hats[argmax.p.l]
          pi_1.u.hat <- p.u.hats[argmin.p.u]
          UIF.l <- IF.l - pi_1.l.hat
          UIF.u <- IF.u - pi_1.u.hat
          
          #superopt influence functions
          supUIF.l0 <- (UIF.l*pAZ0.hat[2] - pi_1.l.hat*(A[i] - pAZ0.hat[2]) 
                        -((1-A[i])*(Y[i] - EY.gA0Z0.hat[2]) 
                          + EY.gA0Z0.hat[,2]*((1-dat.test$A[i]) - pAZ0.hat[,1]))*pAZ0.hat[,2]
                        + EY.gA0Z0.hat[,2]*pAZ0.hat[,1]*(dat.test$A[i] - pAZ0.hat[,2]))/(pAZ0.hat[,2]^2)
          supUIF.u0 <- (UIF.u*pAZ0.hat[,2] - pi_1.u.hat*(dat.test$A[i] - pAZ0.hat[,2]) 
                        -((1-dat.test$A[i])*(dat.test$Y[i] - EY.gA0Z0.hat[,2]) 
                          + EY.gA0Z0.hat[,2]*((1-dat.test$A[i]) - pAZ0.hat[,1]))*pAZ0.hat[,2]
                        + EY.gA0Z0.hat[,2]*pAZ0.hat[,1]*(dat.test$A[i] - pAZ0.hat[,2]))/(pAZ0.hat[,2]^2)
          supIF.l0 <- (pi_1.l.hat - EY.gA0Z0.hat[,2]*pAZ0.hat[,1])/pAZ0.hat[,2] + supUIF.l0
          supIF.u0 <- (pi_1.u.hat - EY.gA0Z0.hat[,2]*pAZ0.hat[,1])/pAZ0.hat[,2] + supUIF.u0
          #lower_bound <- lower_bound + max(0, min(1,(supIF.l0)))
          #upper_bound <- upper_bound + max(0,min(1,(supIF.u0)))
          lower_bound <- lower_bound + (supIF.l0)
          upper_bound <- upper_bound + (supIF.u0)
        }else{
          ## bound fcts
          p.l1 <- function(pi) { pi[7] }
          p.l2 <- function(pi) { pi[3] }
          p.l3 <- function(pi) { pi[3] + pi[4] - pi[5] - pi[8] }
          p.l4 <- function(pi) { pi[2] + pi[3] - pi[5] - pi[6] }
          p.l5 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) }
          p.l6 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
          p.l7 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi))  } 
          p.l8 <- function(pi) {pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi))  }
          gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                         p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
          p.u1 <- function(pi) { 1 - pi[5] }
          p.u2 <- function(pi) { 1 - pi[1] }
          p.u3 <- function(pi) { pi[2] + pi[3] + pi[7] + pi[8] }
          p.u4 <- function(pi) { pi[3] + pi[4] + pi[6] + pi[7] }
          p.u5 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))  + 1 }
          p.u6 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))  + 1 } 
          p.u7 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))  + 1 } 
          p.u8 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))  + 1 }
          gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                         p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
          ## constructing psis and pis
          psi_00.0 <- pYA.gZ0.hat[1] + 
            (1 - Z[i]) * ((1 - Y[i]) * (1 - A[i]) - pYA.gZ0.hat[1]) / (1 - lambda.1)
          psi_01.0 <- pYA.gZ0.hat[2] + 
            (1 - Z[i]) * ((1 - Y[i]) * A[i] - pYA.gZ0.hat[2]) / (1 - lambda.1)
          psi_10.0 <- pYA.gZ0.hat[3] + 
            (1 - Z[i]) * (Y[i] * (1 - A[i]) - pYA.gZ0.hat[3]) / (1 - lambda.1)
          psi_11.0 <- pYA.gZ0.hat[4] + 
            (1 - Z[i]) * (Y[i] * A[i] - pYA.gZ0.hat[4]) / (1 - lambda.1)
          psi_00.1 <- pYA.gZ1.hat[1] + 
            Z[i] * ((1 - Y[i]) * (1 - A[i]) - pYA.gZ1.hat[1]) / lambda.1
          psi_01.1 <- pYA.gZ1.hat[2] + 
            Z[i] * ((1 - Y[i]) * A[i] - pYA.gZ1.hat[2]) / lambda.1
          psi_10.1 <- pYA.gZ1.hat[3] + 
            Z[i] * (Y[i] * (1 - A[i]) - pYA.gZ1.hat[3]) / lambda.1
          psi_11.1 <- pYA.gZ1.hat[4] + 
            Z[i] * (Y[i] * A[i] - pYA.gZ1.hat[4]) / lambda.1
          psi.hats <- c(psi_00.0, psi_01.0, psi_10.0, psi_11.0, 
                        psi_00.1, psi_01.1, psi_10.1, psi_11.1)
          
          ## computing L's and p's
          p.l.hats <- cbind(p.l1(pi.hats)[1], p.l2(pi.hats)[1], p.l3(pi.hats)[1],
                            p.l4(pi.hats)[1], p.l5(pi.hats)[1], p.l6(pi.hats)[1],
                            p.l7(pi.hats)[1], p.l8(pi.hats)[1])
          p.u.hats <- cbind(p.u1(pi.hats)[1], p.u2(pi.hats)[1], p.u3(pi.hats)[1],
                            p.u4(pi.hats)[1], p.u5(pi.hats)[1], p.u6(pi.hats)[1],
                            p.u7(pi.hats)[1], p.u8(pi.hats)[1])
          L.hats <- cbind(p.l1(psi.hats)[1], p.l2(psi.hats)[1], p.l3(psi.hats)[1],
                          p.l4(psi.hats)[1], p.l5(psi.hats)[1], p.l6(psi.hats)[1],
                          p.l7(psi.hats)[1], p.l8(psi.hats)[1])
          U.hats <- cbind(p.u1(psi.hats)[1], p.u2(psi.hats)[1], p.u3(psi.hats)[1],
                          p.u4(psi.hats)[1], p.u5(psi.hats)[1], p.u6(psi.hats)[1],
                          p.u7(psi.hats)[1], p.u8(psi.hats)[1])
          
          ## influence fcts
          #argmax.p.l <- apply(p.l.hats, MARGIN = 1, FUN = which.max)
          argmax.p.l <- which.max(p.l.hats)
          
          #argmin.p.u <- apply(p.u.hats, MARGIN = 1, FUN = which.min)
          argmin.p.u <- which.min(p.u.hats)
          
          IF.l <- L.hats[argmax.p.l]
          IF.u <- U.hats[argmin.p.u]
          
          #centered influence functions
          pi_1.l.hat <- p.l.hats[argmax.p.l]
          pi_1.u.hat <- p.u.hats[argmin.p.u]
          UIF.l <- IF.l - pi_1.l.hat
          UIF.u <- IF.u - pi_1.u.hat
          
          #superopt influence functions
          supUIF.l0 <- (UIF.l*pAZ1.hat[2] - pi_1.l.hat*(A[i] - pAZ1.hat[2]) 
                        -((1-A[i])*(Y[i] - EY.gA0Z1.hat[2]) 
                          + EY.gA0Z1.hat[,2]*((1-dat.test$A[i]) - pAZ1.hat[,1]))*pAZ1.hat[,2]
                        + EY.gA0Z1.hat[,2]*pAZ1.hat[,1]*(dat.test$A[i] - pAZ1.hat[,2]))/(pAZ1.hat[,2]^2)
          supUIF.u0 <- (UIF.u*pAZ1.hat[,2] - pi_1.u.hat*(dat.test$A[i] - pAZ1.hat[,2]) 
                        -((1-dat.test$A[i])*(dat.test$Y[i] - EY.gA0Z1.hat[,2]) 
                          + EY.gA0Z1.hat[,2]*((1-dat.test$A[i]) - pAZ1.hat[,1]))*pAZ1.hat[,2]
                        + EY.gA0Z1.hat[,2]*pAZ1.hat[,1]*(dat.test$A[i] - pAZ1.hat[,2]))/(pAZ1.hat[,2]^2)
          supIF.l0 <- (pi_1.l.hat - EY.gA0Z1.hat[,2]*pAZ1.hat[,1])/pAZ1.hat[,2] + supUIF.l0
          supIF.u0 <- (pi_1.u.hat - EY.gA0Z1.hat[,2]*pAZ1.hat[,1])/pAZ1.hat[,2] + supUIF.u0
          #lower_bound <- lower_bound + max(0,min(1,supIF.l0))
          #upper_bound <- upper_bound + max(0,min(1,supIF.u0))
          lower_bound <- lower_bound + (supIF.l0)
          upper_bound <- upper_bound + (supIF.u0)
        }
      }
    }else{
      if(A[i] == 0){
        #E(Y^1 | A = 0, L = l, Z = z) bounds
        if(Z[i] == 0){
          p.l1 <- function(pi) { pi[4] }
          p.l2 <- function(pi) { pi[8] }
          p.l3 <- function(pi) { pi[5] + pi[8] - pi[1] - pi[2] }
          p.l4 <- function(pi) { pi[7] + pi[8] - pi[2] - pi[3] }
          p.l5 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) }
          p.l6 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
          p.l7 <- function(pi) {pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
          p.l8 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi))  }
          
          gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                         p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
          
          
          #a = 1 intervention
          p.u1 <- function(pi) { 1 - pi[6] }
          p.u2 <- function(pi) { 1 - pi[2] }
          p.u3 <- function(pi) { pi[1] + pi[4] + pi[7] + pi[8] }
          p.u4 <- function(pi) { pi[3] + pi[4] + pi[5] + pi[8] }
          p.u5 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1  }
          p.u6 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1  }
          p.u7 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) +1} 
          p.u8 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) +1}
          
          gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                         p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
          
          ## constructing psis and pis
          psi_00.0 <- pYA.gZ0.hat[1] + 
            (1 - Z[i]) * ((1 - Y[i]) * (1 - A[i]) - pYA.gZ0.hat[1]) / (1 - lambda.1)
          psi_01.0 <- pYA.gZ0.hat[2] + 
            (1 - Z[i]) * ((1 - Y[i]) * A[i] - pYA.gZ0.hat[2]) / (1 - lambda.1)
          psi_10.0 <- pYA.gZ0.hat[3] + 
            (1 - Z[i]) * (Y[i] * (1 - A[i]) - pYA.gZ0.hat[3]) / (1 - lambda.1)
          psi_11.0 <- pYA.gZ0.hat[4] + 
            (1 - Z[i]) * (Y[i] * A[i] - pYA.gZ0.hat[4]) / (1 - lambda.1)
          psi_00.1 <- pYA.gZ1.hat[1] + 
            Z[i] * ((1 - Y[i]) * (1 - A[i]) - pYA.gZ1.hat[1]) / lambda.1
          psi_01.1 <- pYA.gZ1.hat[2] + 
            Z[i] * ((1 - Y[i]) * A[i] - pYA.gZ1.hat[2]) / lambda.1
          psi_10.1 <- pYA.gZ1.hat[3] + 
            Z[i] * (Y[i] * (1 - A[i]) - pYA.gZ1.hat[3]) / lambda.1
          psi_11.1 <- pYA.gZ1.hat[4] + 
            Z[i] * (Y[i] * A[i] - pYA.gZ1.hat[4]) / lambda.1
          psi.hats <- c(psi_00.0, psi_01.0, psi_10.0, psi_11.0, 
                        psi_00.1, psi_01.1, psi_10.1, psi_11.1)
          
          ## computing L's and p's
          p.l.hats <- cbind(p.l1(pi.hats)[1], p.l2(pi.hats)[1], p.l3(pi.hats)[1],
                            p.l4(pi.hats)[1], p.l5(pi.hats)[1], p.l6(pi.hats)[1],
                            p.l7(pi.hats)[1], p.l8(pi.hats)[1])
          p.u.hats <- cbind(p.u1(pi.hats)[1], p.u2(pi.hats)[1], p.u3(pi.hats)[1],
                            p.u4(pi.hats)[1], p.u5(pi.hats)[1], p.u6(pi.hats)[1],
                            p.u7(pi.hats)[1], p.u8(pi.hats)[1])
          L.hats <- cbind(p.l1(psi.hats)[1], p.l2(psi.hats)[1], p.l3(psi.hats)[1],
                          p.l4(psi.hats)[1], p.l5(psi.hats)[1], p.l6(psi.hats)[1],
                          p.l7(psi.hats)[1], p.l8(psi.hats)[1])
          U.hats <- cbind(p.u1(psi.hats)[1], p.u2(psi.hats)[1], p.u3(psi.hats)[1],
                          p.u4(psi.hats)[1], p.u5(psi.hats)[1], p.u6(psi.hats)[1],
                          p.u7(psi.hats)[1], p.u8(psi.hats)[1])
          
          ## influence fcts
          argmax.p.l <- apply(p.l.hats, MARGIN = 1, FUN = which.max)
          
          argmin.p.u <- apply(p.u.hats, MARGIN = 1, FUN = which.min)
          
          IF.l <- L.hats[argmax.p.l]
          IF.u <- U.hats[argmin.p.u]
          
          #centered influence functions
          pi_1.l.hat <- p.l.hats[argmax.p.l]
          pi_1.u.hat <- p.u.hats[argmin.p.u]
          UIF.l <- IF.l - pi_1.l.hat
          UIF.u <- IF.u - pi_1.u.hat
          supUIF.l1 <- (UIF.l*pAZ0.hat[,1] - pi_1.l.hat*((1-dat.test$A[i]) - pAZ0.hat[,1]) 
                        -((dat.test$A[i])*(dat.test$Y[i] - EY.gA1Z0.hat[,2]) 
                          + EY.gA1Z0.hat[,2]*((dat.test$A[i]) - pAZ0.hat[,2]))*pAZ0.hat[,1]
                        + EY.gA1Z0.hat[,2]*pAZ0.hat[,2]*((1-dat.test$A[i]) - pAZ0.hat[,1]))/(pAZ0.hat[,1]^2)
          supUIF.u1 <- (UIF.u*pAZ0.hat[,1] - pi_1.u.hat*((1-dat.test$A[i]) - pAZ0.hat[,1]) 
                        -((dat.test$A[i])*(dat.test$Y[i] - EY.gA1Z0.hat[,2]) 
                          + EY.gA1Z0.hat[,2]*((dat.test$A[i]) - pAZ0.hat[,2]))*pAZ0.hat[,1]
                        + EY.gA1Z0.hat[,2]*pAZ0.hat[,2]*((1-dat.test$A[i]) - pAZ0.hat[,1]))/(pAZ0.hat[,1]^2)
          supIF.l1 <- (pi_1.l.hat - EY.gA1Z0.hat[,2]*pAZ0.hat[,2])/pAZ0.hat[,1] + supUIF.l1
          supIF.u1 <- (pi_1.u.hat - EY.gA1Z0.hat[,2]*pAZ0.hat[,2])/pAZ0.hat[,1] + supUIF.u1
          #lower_bound <- lower_bound + max(0,min(1,(supIF.l1)))
          #upper_bound <- upper_bound + max(0,min(1,(supIF.u1)))
          lower_bound <- lower_bound + (supIF.l1)
          upper_bound <- upper_bound + (supIF.u1)
        }else{
          p.l1 <- function(pi) { pi[4] }
          p.l2 <- function(pi) { pi[8] }
          p.l3 <- function(pi) { pi[5] + pi[8] - pi[1] - pi[2] }
          p.l4 <- function(pi) { pi[7] + pi[8] - pi[2] - pi[3] }
          p.l5 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) }
          p.l6 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
          p.l7 <- function(pi) {pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
          p.l8 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi))  }
          
          gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                         p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
          
          
          #a = 1 intervention
          p.u1 <- function(pi) { 1 - pi[6] }
          p.u2 <- function(pi) { 1 - pi[2] }
          p.u3 <- function(pi) { pi[1] + pi[4] + pi[7] + pi[8] }
          p.u4 <- function(pi) { pi[3] + pi[4] + pi[5] + pi[8] }
          p.u5 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- max(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1  }
          p.u6 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- max(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1  }
          p.u7 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- max(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 } 
          p.u8 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- max(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 }
          
          gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                         p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
          
          ## constructing psis and pis
          psi_00.0 <- pYA.gZ0.hat[1] + 
            (1 - Z[i]) * ((1 - Y[i]) * (1 - A[i]) - pYA.gZ0.hat[1]) / (1 - lambda.1)
          psi_01.0 <- pYA.gZ0.hat[2] + 
            (1 - Z[i]) * ((1 - Y[i]) * A[i] - pYA.gZ0.hat[2]) / (1 - lambda.1)
          psi_10.0 <- pYA.gZ0.hat[3] + 
            (1 - Z[i]) * (Y[i] * (1 - A[i]) - pYA.gZ0.hat[3]) / (1 - lambda.1)
          psi_11.0 <- pYA.gZ0.hat[4] + 
            (1 - Z[i]) * (Y[i] * A[i] - pYA.gZ0.hat[4]) / (1 - lambda.1)
          psi_00.1 <- pYA.gZ1.hat[1] + 
            Z[i] * ((1 - Y[i]) * (1 - A[i]) - pYA.gZ1.hat[1]) / lambda.1
          psi_01.1 <- pYA.gZ1.hat[2] + 
            Z[i] * ((1 - Y[i]) * A[i] - pYA.gZ1.hat[2]) / lambda.1
          psi_10.1 <- pYA.gZ1.hat[3] + 
            Z[i] * (Y[i] * (1 - A[i]) - pYA.gZ1.hat[3]) / lambda.1
          psi_11.1 <- pYA.gZ1.hat[4] + 
            Z[i] * (Y[i] * A[i] - pYA.gZ1.hat[4]) / lambda.1
          psi.hats <- c(psi_00.0, psi_01.0, psi_10.0, psi_11.0, 
                        psi_00.1, psi_01.1, psi_10.1, psi_11.1)
          
          ## computing L's and p's
          p.l.hats <- cbind(p.l1(pi.hats)[1], p.l2(pi.hats)[1], p.l3(pi.hats)[1],
                            p.l4(pi.hats)[1], p.l5(pi.hats)[1], p.l6(pi.hats)[1],
                            p.l7(pi.hats)[1], p.l8(pi.hats)[1])
          p.u.hats <- cbind(p.u1(pi.hats)[1], p.u2(pi.hats)[1], p.u3(pi.hats)[1],
                            p.u4(pi.hats)[1], p.u5(pi.hats)[1], p.u6(pi.hats)[1],
                            p.u7(pi.hats)[1], p.u8(pi.hats)[1])
          L.hats <- cbind(p.l1(psi.hats)[1], p.l2(psi.hats)[1], p.l3(psi.hats)[1],
                          p.l4(psi.hats)[1], p.l5(psi.hats)[1], p.l6(psi.hats)[1],
                          p.l7(psi.hats)[1], p.l8(psi.hats)[1])
          U.hats <- cbind(p.u1(psi.hats)[1], p.u2(psi.hats)[1], p.u3(psi.hats)[1],
                          p.u4(psi.hats)[1], p.u5(psi.hats)[1], p.u6(psi.hats)[1],
                          p.u7(psi.hats)[1], p.u8(psi.hats)[1])
          
          ## influence fcts
          argmax.p.l <- apply(p.l.hats, MARGIN = 1, FUN = which.max)
          
          argmin.p.u <- apply(p.u.hats, MARGIN = 1, FUN = which.min)
          
          IF.l <- L.hats[argmax.p.l]
          IF.u <- U.hats[argmin.p.u]
          
          #centered influence functions
          pi_1.l.hat <- p.l.hats[argmax.p.l]
          pi_1.u.hat <- p.u.hats[argmin.p.u]
          UIF.l <- IF.l - pi_1.l.hat
          UIF.u <- IF.u - pi_1.u.hat
          supUIF.l1 <- (UIF.l*pAZ1.hat[,1] - pi_1.l.hat*((1-dat.test$A[i]) - pAZ1.hat[,1]) 
                        -((dat.test$A[i])*(dat.test$Y[i] - EY.gA1Z1.hat[,2]) 
                          + EY.gA1Z1.hat[,2]*((dat.test$A[i]) - pAZ1.hat[,2]))*pAZ1.hat[,1]
                        + EY.gA1Z1.hat[,2]*pAZ1.hat[,2]*((1-dat.test$A[i]) - pAZ1.hat[,1]))/(pAZ1.hat[,1]^2)
          supUIF.u1 <- (UIF.u*pAZ1.hat[,1] - pi_1.u.hat*((1-dat.test$A[i]) - pAZ1.hat[,1]) 
                        -((dat.test$A[i])*(dat.test$Y[i] - EY.gA1Z1.hat[,2]) 
                          + EY.gA1Z1.hat[,2]*((dat.test$A[i]) - pAZ0.hat[,2]))*pAZ1.hat[,1]
                        + EY.gA1Z1.hat[,2]*pAZ1.hat[,2]*((1-dat.test$A[i]) - pAZ1.hat[,1]))/(pAZ1.hat[,1]^2)
          supIF.l1 <- (pi_1.l.hat - EY.gA1Z1.hat[,2]*pAZ1.hat[,2])/pAZ1.hat[,1] + supUIF.l1
          supIF.u1 <- (pi_1.u.hat - EY.gA1Z1.hat[,2]*pAZ1.hat[,2])/pAZ1.hat[,1] + supUIF.u1
          #lower_bound <- lower_bound + max(0,min(1,(supIF.l1)))
          #upper_bound <- upper_bound + max(0,min(1,(supIF.u1)))
          lower_bound <- lower_bound + (supIF.l1)
          upper_bound <- upper_bound + (supIF.u1)
        }
      }else{
        #E(Y | A = 1, L = l, Z = z) bounds
        if(Z[i] == 0){
          EY.gA1Z0.hat <- predict(EY.gA1Z0, type = "prob", newdata = L[i,])[,2]
          #lower_bound = lower_bound + max(0,min(1,(EY.gA1Z0.hat + Y[i] - EY.gA1Z0.hat)))
          #upper_bound = upper_bound + max(0,min(1,(EY.gA1Z0.hat + Y[i] - EY.gA1Z0.hat)))
          lower_bound = lower_bound + (EY.gA1Z0.hat + Y[i] - EY.gA1Z0.hat)
          upper_bound = upper_bound + (EY.gA1Z0.hat + Y[i] - EY.gA1Z0.hat)
        }else{
          EY.gA1Z1.hat <- predict(EY.gA1Z1, type = "prob", newdata = L[i,])[,2]
          #lower_bound = lower_bound + max(0,min(1,(EY.gA1Z1.hat + Y[i] - EY.gA1Z1.hat)))
          #upper_bound = upper_bound + max(0,min(1,(EY.gA1Z1.hat + Y[i] - EY.gA1Z1.hat)))
          lower_bound = lower_bound + (EY.gA1Z1.hat + Y[i] - EY.gA1Z1.hat)
          upper_bound = upper_bound + (EY.gA1Z1.hat + Y[i] - EY.gA1Z1.hat)
        }
      }
    }
    hist.l <- c(hist.l, lower_bound)
    hist.incr.l <- c(hist.incr.l, hist.l[length(hist.l)] - hist.l[length(hist.l) - 1])
    hist.u <- c(hist.u, upper_bound)
    hist.incr.u <- c(hist.incr.u, hist.u[length(hist.u)] - hist.u[length(hist.u) - 1])
  }
  lower_bound <- lower_bound/nrow(dat.test)
  upper_bound <- upper_bound/nrow(dat.test)
  return(c(lower = lower_bound, upper = upper_bound, lower.var = var(hist.incr.l), upper.var = var(hist.incr.u)))
}


#######################################################################################################
#decision criteria
compute_bounds <- function(dat, l){
  #ATE
  p.l1 <- function(pi) { pi[8] + pi[1] - 1 }
  p.l2 <- function(pi) { pi[4] + pi[5] - 1 }
  p.l3 <- function(pi) { -pi[6] - pi[7] }
  p.l4 <- function(pi) { -pi[2] - pi[3] }
  p.l5 <- function(pi) { pi[4] - pi[8] - pi[7] - pi[2] - pi[3] }
  p.l6 <- function(pi) { pi[8] - pi[4] - pi[3] - pi[6] - pi[7] }
  p.l7 <- function(pi) { pi[5] - pi[6] - pi[7] - pi[2] - pi[1] }
  p.l8 <- function(pi) { pi[1] - pi[2] - pi[3] - pi[6] - pi[5] }
  gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                 p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
  arg_gamma.l <- function(pi) {which.max(c(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                           p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)))}
  p.u1 <- function(pi) { 1 - pi[6] - pi[3] }
  p.u2 <- function(pi) { 1 - pi[2] - pi[7] }
  p.u3 <- function(pi) { pi[8] + pi[5] }
  p.u4 <- function(pi) { pi[4] + pi[1] }
  p.u5 <- function(pi) { -pi[2] + pi[6] + pi[5] + pi[4] + pi[1] }
  p.u6 <- function(pi) { -pi[6] + pi[2] + pi[1] + pi[8] + pi[5] }
  p.u7 <- function(pi) { -pi[7] + pi[8] + pi[5] + pi[4] + pi[3] }
  p.u8 <- function(pi) { -pi[3] + pi[4] + pi[1] + pi[8] + pi[7] }
  gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                 p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
  arg_gamma.u <- function(pi) {which.min(c(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                           p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)))}
  pi.ya.0 <- multinom(Y.A ~ ., data = dat[dat$Z==0, -(1:3)], family = "binomial")
  pi.ya.1 <- multinom(Y.A ~ ., data = dat[dat$Z==1, -(1:3)], family = "binomial")
  pi.hat <- data.frame(predict(pi.ya.0, type = 'probs', newdata = l))
  pi.hat <- data.frame(t(pi.hat))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  ATE_bounds <- c(gamma.l(pi.hat), gamma.u(pi.hat))
  #CATE with A = 0
  p.l1 <- function(pi) { pi[4] }
  p.l2 <- function(pi) { pi[8] }
  p.l3 <- function(pi) { pi[5] + pi[8] - pi[1] - pi[2] }
  p.l4 <- function(pi) { pi[7] + pi[8] - pi[2] - pi[3] }
  p.l5 <- function(pi) { -1 }
  p.l6 <- function(pi) { -1  } 
  p.l7 <- function(pi) { -1  }
  p.l8 <- function(pi) { -1  }
  gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                 p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
  arg_gamma.l <- function(pi) {which.max(c(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                           p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)))}
  p.u1 <- function(pi) { 1 - pi[6] }
  p.u2 <- function(pi) { 1 - pi[2] }
  p.u3 <- function(pi) { pi[1] + pi[4] + pi[7] + pi[8] }
  p.u4 <- function(pi) { pi[3] + pi[4] + pi[5] + pi[8] }
  p.u5 <- function(pi) { 1  }
  p.u6 <- function(pi) { 1 } 
  p.u7 <- function(pi) { 1 } 
  p.u8 <- function(pi) { 1}
  gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                 p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
  arg_gamma.u <- function(pi) {which.min(c(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                           p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)))}
  pi.ya.0 <- multinom(Y.A ~ ., data = dat[dat$Z==0, -(1:3)], family = "binomial")
  pi.ya.1 <- multinom(Y.A ~ ., data = dat[dat$Z==1, -(1:3)], family = "binomial")
  pi.hat <- data.frame(predict(pi.ya.0, type = 'probs', newdata = l))
  pi.hat <- data.frame(t(pi.hat))
  colnames(pi.hat) <-
    c("pi_00.0", "pi_01.0", "pi_10.0", "pi_11.0")
  pi.hat <- cbind.data.frame(pi.hat, data.frame(t(predict(pi.ya.1, type = 'probs', newdata = l))))
  colnames(pi.hat)[(ncol(pi.hat) - 3):ncol(pi.hat)] <-
    c("pi_00.1", "pi_01.1", "pi_10.1", "pi_11.1")
  pi.hat <- as.numeric(pi.hat)
  
  fmY <- glm(Y~., data = dat[, -c(1,2,ncol(dat))], family = "binomial")
  EY <- predict(fmY, newdata = l, type = "response")
  fmY.A0 <- glm(Y~., data = dat[which(A == 0),-c(1,2,ncol(dat))], family = "binomial")
  EY.A0 <- predict(fmY.A0, newdata = l, type = "response")
  fmY.A1 <- glm(Y~., data = dat[which(A == 1),-c(1,2,ncol(dat))], family = "binomial")
  EY.A1 <- predict(fmY.A1, newdata = l, type = "response")
  fmA <- glm(A~., data = dat[, -c(1,3,ncol(dat))], family = "binomial")
  pA <- predict(fmA, newdata = l, type = "response")
  CATE0_bounds <- c((gamma.l(pi.hat) - EY)/(1-pA), (gamma.u(pi.hat) - EY)/(1-pA))
  #CATE with A = 1
  p.l1 <- function(pi) { pi[7] }
  p.l2 <- function(pi) { pi[3] }
  p.l3 <- function(pi) { pi[3] + pi[4] - pi[5] - pi[8] }
  p.l4 <- function(pi) { pi[2] + pi[3] - pi[5] - pi[6] }
  p.l5 <- function(pi) { -1}
  p.l6 <- function(pi) { -1  } 
  p.l7 <- function(pi) { -1  } 
  p.l8 <- function(pi) { -1  }
  gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                 p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
  arg_gamma.l <- function(pi) {which.max(c(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                           p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)))}
  p.u1 <- function(pi) { 1 - pi[5] }
  p.u2 <- function(pi) { 1 - pi[1] }
  p.u3 <- function(pi) { pi[2] + pi[3] + pi[7] + pi[8] }
  p.u4 <- function(pi) { pi[3] + pi[4] + pi[6] + pi[7] }
  p.u5 <- function(pi) { 1  }
  p.u6 <- function(pi) { 1  } 
  p.u7 <- function(pi) { 1 } 
  p.u8 <- function(pi) { 1 }
  gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                 p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
  arg_gamma.u <- function(pi) {which.min(c(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                           p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)))}
  #Balke Pearl bounds(superoptimal)
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
  
  CATE1_bounds <- c((-gamma.u(pi.hat) + EY)/pA, (- gamma.l(pi.hat) + EY)/pA)
  #E(Y^0 | L)
  ## lower bound functions
  #a = 0 intervention
  p.l1 <- function(pi) { pi[7] }
  p.l2 <- function(pi) { pi[3] }
  p.l3 <- function(pi) { pi[3] + pi[4] - pi[5] - pi[8] }
  p.l4 <- function(pi) { pi[2] + pi[3] - pi[5] - pi[6] }
  p.l5 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -1 }
  p.l6 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -1 } 
  p.l7 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -1  } 
  p.l8 <- function(pi) {pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -1  }
  
  gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                 p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
  
  ## upper bound functions
  #a = 0 intervention
  p.u1 <- function(pi) { 1 - pi[5] }
  p.u2 <- function(pi) { 1 - pi[1] }
  p.u3 <- function(pi) { pi[2] + pi[3] + pi[7] + pi[8] }
  p.u4 <- function(pi) { pi[3] + pi[4] + pi[6] + pi[7] }
  p.u5 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 }
  p.u6 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 } 
  p.u7 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 } 
  p.u8 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 }
  
  gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                 p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
  
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
  Y0_bounds <- c(gamma.l(pi.hat), gamma.u(pi.hat))
  #E(Y^1 | L)
  ## lower bound functions
  #a = 1 intervention
  p.l1 <- function(pi) { pi[4] }
  p.l2 <- function(pi) { pi[8] }
  p.l3 <- function(pi) { pi[5] + pi[8] - pi[1] - pi[2] }
  p.l4 <- function(pi) { pi[7] + pi[8] - pi[2] - pi[3] }
  p.l5 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) }
  p.l6 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
  p.l7 <- function(pi) {pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
  p.l8 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi))  }
  
  gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                 p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
  
  ## upper bound functions
  #a = 1 intervention
  p.u1 <- function(pi) { 1 - pi[6] }
  p.u2 <- function(pi) { 1 - pi[2] }
  p.u3 <- function(pi) { pi[1] + pi[4] + pi[7] + pi[8] }
  p.u4 <- function(pi) { pi[3] + pi[4] + pi[5] + pi[8] }
  p.u5 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))-pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1  }
  p.u6 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))-pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1  }
  p.u7 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))-pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 } 
  p.u8 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))-pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 }
  
  gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                 p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
  
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
  Y1_bounds <- c(gamma.l(pi.hat), gamma.u(pi.hat))
  #E(Y^1 | A = 0, L)
  ## lower bound functions
  #a = 1 intervention
  p.l1 <- function(pi) { pi[4] }
  p.l2 <- function(pi) { pi[8] }
  p.l3 <- function(pi) { pi[5] + pi[8] - pi[1] - pi[2] }
  p.l4 <- function(pi) { pi[7] + pi[8] - pi[2] - pi[3] }
  p.l5 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) }
  p.l6 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
  p.l7 <- function(pi) {pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
  p.l8 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi))  }
  
  gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                 p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
  
  ## upper bound functions
  #a = 1 intervention
  p.u1 <- function(pi) { 1 - pi[6] }
  p.u2 <- function(pi) { 1 - pi[2] }
  p.u3 <- function(pi) { pi[1] + pi[4] + pi[7] + pi[8] }
  p.u4 <- function(pi) { pi[3] + pi[4] + pi[5] + pi[8] }
  p.u5 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1  }
  p.u6 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1  }
  p.u7 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 } 
  p.u8 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi))- pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 }
  
  gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                 p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
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
  Y1.A0_bounds <- c((gamma.l(pi.hat) - EY.A1*pA)/(1-pA), (gamma.u(pi.hat) - EY.A1*pA)/(1-pA))
  #E(Y^0 | A = 1, L)
  ## lower bound functions
  #a = 0 intervention
  p.l1 <- function(pi) { pi[7] }
  p.l2 <- function(pi) { pi[3] }
  p.l3 <- function(pi) { pi[3] + pi[4] - pi[5] - pi[8] }
  p.l4 <- function(pi) { pi[2] + pi[3] - pi[5] - pi[6] }
  p.l5 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) }
  p.l6 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) } 
  p.l7 <- function(pi) { pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi))  } 
  p.l8 <- function(pi) {pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi)) -pmin(p.l1(pi), p.l2(pi),p.l3(pi), p.l4(pi))  }
  
  gamma.l <- function(pi) { pmax(p.l1(pi), p.l2(pi), p.l3(pi), p.l4(pi),
                                 p.l5(pi), p.l6(pi), p.l7(pi), p.l8(pi)) }
  
  ## upper bound functions
  #a = 0 intervention
  p.u1 <- function(pi) { 1 - pi[5] }
  p.u2 <- function(pi) { 1 - pi[1] }
  p.u3 <- function(pi) { pi[2] + pi[3] + pi[7] + pi[8] }
  p.u4 <- function(pi) { pi[3] + pi[4] + pi[6] + pi[7] }
  p.u5 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 }
  p.u6 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 } 
  p.u7 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 } 
  p.u8 <- function(pi) { pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) - pmax(p.u1(pi), p.u2(pi),p.u3(pi), p.u4(pi)) + 1 }
  
  gamma.u <- function(pi) { pmin(p.u1(pi), p.u2(pi), p.u3(pi), p.u4(pi),
                                 p.u5(pi), p.u6(pi), p.u7(pi), p.u8(pi)) }
  
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
  Y0.A1_bounds <- c((gamma.l(pi.hat) - EY.A0*(1-pA))/pA, (gamma.u(pi.hat) - EY.A0*(1-pA))/pA)
  
  return(data.frame(Y0_bounds = Y0_bounds, Y1_bounds = Y1_bounds,
                    Y1.A0_bounds = Y1.A0_bounds, Y0.A1_bounds = Y0.A1_bounds,
                    ATE_bounds = ATE_bounds, CATE0_bounds = CATE0_bounds,
                    CATE1_bounds = CATE1_bounds))
}

#Computing decision criteria regimes

maximax_opt <- function(a,l,z){
  bounds <- compute_bounds(dat.t, l)
  if(bounds$Y1_bounds[2] > bounds$Y0_bounds[2]){
    return(1)
  }else{
    return(0)
  }
}

maximin_opt <- function(a,l,z){
  bounds <- compute_bounds(dat.t, l)
  if(bounds$Y1_bounds[1] > bounds$Y0_bounds[1]){
    return(1)
  }else{
    return(0)
  }
}

minimax_opt <- function(a,l,z){
  bounds <- compute_bounds(dat.t, l)
  if(bounds$ATE_bounds[2] < 0){
    return(0)
  }else{
    if(bounds$ATE_bounds[1] >= 0){
      return(1)
    }else{
      if(abs(bounds$ATE_bounds[2]) > abs(bounds$ATE_bounds[1])){
        return(1)
      }else{
        return(0)
      }
    }
  }
}

healthcare_opt <- function(a,l,z){
  bounds <- compute_bounds(dat.t, l)
  if(bounds$ATE_bounds[1] >= 0){
    return(1)
  }else{
    return(0)
  }
}

maximax_superopt <- function(a,l,z){
  bounds <- compute_bounds(dat.t, l)
  fmY.A0 <- glm(Y~., data = dat.t[dat.t$A == 0, -c(1,2,ncol(dat.t))]) 
  EY.A0 <- predict(fmY.A0,l, type = "response")
  fmY.A1 <- glm(Y~., data = dat.t[dat.t$A == 1, -c(1,2,ncol(dat.t))]) 
  EY.A1 <- predict(fmY.A1,l, type = "response")
  if(a == 0){
    if(bounds$Y1.A0_bounds[2] > EY.A0){
      return(1)
    }else{
      return(0)
    }
  }else{
    if(EY.A1 > bounds$Y0.A1_bounds[2]){
      return(1)
    }else{
      return(0)
    }
  }
}

maximin_superopt <- function(a,l,z){
  bounds <- compute_bounds(dat.t, l)
  fmY.A0 <- glm(Y~., data = dat[dat$A == 0, -c(1,2,ncol(dat.t))]) 
  EY.A0 <- predict(fmY.A0,l, type = "response")
  fmY.A1 <- glm(Y~., data = dat[dat$A == 1, -c(1,2,ncol(dat.t))]) 
  EY.A1 <- predict(fmY.A1,l, type = "response")
  if(a == 0){
    if(bounds$Y1.A0_bounds[1] > EY.A0){
      return(1)
    }else{
      return(0)
    }
  }else{
    if(EY.A1 > bounds$Y0.A1_bounds[1]){
      return(1)
    }else{
      return(0)
    }
  }
}

minimax_superopt <- function(a,l,z){
  bounds <- compute_bounds(dat.t, l)
  if(a == 0){
    if(bounds$CATE0_bounds[2] < 0){
      return(0)
    }else{
      if(bounds$CATE0_bounds[1] >= 0){
        return(1)
      }else{
        if(abs(bounds$CATE0_bounds[2]) > abs(bounds$CATE0_bounds[1])){
          return(1)
        }else{
          return(0)
        }
      }
    }
  }else{
    if(bounds$CATE1_bounds[2] < 0){
      return(0)
    }else{
      if(bounds$CATE1_bounds[1] >= 0){
        return(1)
      }else{
        if(abs(bounds$CATE1_bounds[2]) > abs(bounds$CATE1_bounds[1])){
          return(1)
        }else{
          return(0)
        }
      }
    }
  }
}

healthcare_superopt <- function(a,l,z){
  bounds <- compute_bounds(dat.t, l)
  if(a == 0){
    if(bounds$CATE0_bounds[1] >= 0){
      return(1)
    }else{
      return(0)
    }
  }else{
    if(bounds$CATE1_bounds[1] >= 0){
      return(1)
    }else{
      return(0)
    }
  }
}

healthcare_obs_baseline <- function(a,l,z){
  bounds <- compute_bounds(dat.t, l)
  if(a == 0){
    if(bounds$CATE0_bounds[1] > 0){
      return(1)
    }else{
      return(0)
    }
  }else{
    if(bounds$CATE1_bounds[2] < 0){
      return(0)
    }else{
      return(1)
    }
  }
}

healthcare_obs_baseline.opt <- function(a,l,z){
  bounds <- compute_bounds(dat.t, l)
  fmY <- glm(Y~., data = dat[-c(1,2,ncol(dat.t))]) 
  EY <- predict(fmY,l, type = "response")
  I_1 <- 0
  I_2 <- 0
  I_3 <- 0
  if(bounds$Y0_bounds[1] > EY){
    I_1 <- 1
  }
  if(bounds$Y1_bounds[1] <= bounds$Y0_bounds[2]){
    I_2 <- 1
  }
  if(bounds$Y1_bounds[1] > EY){
    I_3 <- 1
  }
  return(a*(1-I_1*I_2) + (1-a)*(I_3*(1-I_2)*I_1 + I_3*(1-I_1)))
}

######################################################################################################
#maximax (optimal)
results.maximax_opt <- lapply(test.indices, function(inds) {
  bound_regime(g = maximax_opt, dat.train = dat.te[-inds,], dat.test = dat.te[inds,])
})
lower_bound.maximax_opt <-  mean(sapply(results.maximax_opt, function(x) {x[1]}, simplify = 0))
upper_bound.maximax_opt <-  mean(sapply(results.maximax_opt, function(x) {x[2]}, simplify = 0))
maximax_opt.low.var <- mean(sapply(results.maximax_opt, function(x) {x[3]}, simplify = 0))
maximax_opt.upp.var <- mean(sapply(results.maximax_opt, function(x) {x[4]}, simplify = 0))
maximax_CI <- c(lower_bound.maximax_opt - qnorm(0.975)*sqrt(maximax_opt.low.var/n), upper_bound.maximax_opt + qnorm(0.975)*sqrt(maximax_opt.upp.var/n))

#maximin (optimal)
results.maximin_opt <- lapply(test.indices, function(inds) {
  bound_regime(g = maximin_opt, dat.train = dat.te[-inds,], dat.test = dat.te[inds,])
})
lower_bound.maximin_opt <-  mean(sapply(results.maximin_opt, function(x) {x[1]}, simplify = 0))
upper_bound.maximin_opt <-  mean(sapply(results.maximin_opt, function(x) {x[2]}, simplify = 0))
maximin_opt.low.var <- mean(sapply(results.maximin_opt, function(x) {x[3]}, simplify = 0))
maximin_opt.upp.var <- mean(sapply(results.maximin_opt, function(x) {x[4]}, simplify = 0))
maximin_CI <- c(lower_bound.maximin_opt - qnorm(0.975)*sqrt(maximin_opt.low.var/n), upper_bound.maximin_opt + qnorm(0.975)*sqrt(maximin_opt.upp.var/n))


#minimax (optimal)
results.minimax_opt <- lapply(test.indices, function(inds) {
  bound_regime(g = minimax_opt, dat.train = dat.te[-inds,], dat.test = dat.te[inds,])
})
lower_bound.minimax_opt <-  mean(sapply(results.minimax_opt, function(x) {x[1]}, simplify = 0))
upper_bound.minimax_opt <-  mean(sapply(results.minimax_opt, function(x) {x[2]}, simplify = 0))
minimax_opt.low.var <- mean(sapply(results.minimax_opt, function(x) {x[3]}, simplify = 0))
minimax_opt.upp.var <- mean(sapply(results.minimax_opt, function(x) {x[4]}, simplify = 0))
minimax_CI <- c(lower_bound.minimax_opt - qnorm(0.975)*sqrt(minimax_opt.low.var/n), upper_bound.minimax_opt + qnorm(0.975)*sqrt(minimax_opt.upp.var/n))


#healthcare (optimal)
results.healthcare_opt <- lapply(test.indices, function(inds) {
  bound_regime(g = healthcare_opt, dat.train = dat.te[-inds,], dat.test = dat.te[inds,])
})
lower_bound.healthcare_opt <-  mean(sapply(results.healthcare_opt, function(x) {x[1]}, simplify = 0))
upper_bound.healthcare_opt <-  mean(sapply(results.healthcare_opt, function(x) {x[2]}, simplify = 0))
#healthcare_opt.low.var <- mean(sapply(results.maximin_opt, function(x) {x[3]}, simplify = 0))
#healthcare_opt.upp.var <- mean(sapply(results.maximin_opt, function(x) {x[4]}, simplify = 0))
healthcare_opt.low.var <- mean(sapply(results.healthcare_opt, function(x) {x[3]}, simplify = 0))
healthcare_opt.upp.var <- mean(sapply(results.healthcare_opt, function(x) {x[4]}, simplify = 0))
healthcare_CI <- c(lower_bound.healthcare_opt - qnorm(0.975)*sqrt(healthcare_opt.low.var/n), upper_bound.healthcare_opt + qnorm(0.975)*sqrt(healthcare_opt.upp.var/n))


#maximax (superoptimal)
results.maximax_superopt <- lapply(test.indices, function(inds) {
  bound_regime(g = maximax_superopt, dat.train = dat.te[-inds,], dat.test = dat.te[inds,])
})
lower_bound.maximax_superopt <-  mean(sapply(results.maximax_superopt, function(x) {x[1]}, simplify = 0))
upper_bound.maximax_superopt <-  mean(sapply(results.maximax_superopt, function(x) {x[2]}, simplify = 0))
maximax_superopt.low.var <- mean(sapply(results.maximax_superopt, function(x) {x[3]}, simplify = 0))
maximax_superopt.upp.var <- mean(sapply(results.maximax_superopt, function(x) {x[4]}, simplify = 0))
maximax_superopt_CI <- c(lower_bound.maximax_superopt - qnorm(0.975)*sqrt(maximax_superopt.low.var/n), upper_bound.maximax_superopt + qnorm(0.975)*sqrt(maximax_superopt.upp.var/n))


#maximin (superoptimal)
results.maximin_superopt <- lapply(test.indices, function(inds) {
  bound_regime(g = maximin_superopt, dat.train = dat.te[-inds,], dat.test = dat.te[inds,])
})
lower_bound.maximin_superopt <-  mean(sapply(results.maximin_superopt, function(x) {x[1]}, simplify = 0))
upper_bound.maximin_superopt <-  mean(sapply(results.maximin_superopt, function(x) {x[2]}, simplify = 0))
maximin_superopt.low.var <- mean(sapply(results.maximin_superopt, function(x) {x[3]}, simplify = 0))
maximin_superopt.upp.var <- mean(sapply(results.maximin_superopt, function(x) {x[4]}, simplify = 0))
maximin_superopt_CI <- c(lower_bound.maximin_superopt - qnorm(0.975)*sqrt(maximin_superopt.low.var/n), upper_bound.maximin_superopt + qnorm(0.975)*sqrt(maximin_superopt.upp.var/n))


#minimax (superoptimal)
results.minimax_superopt <- lapply(test.indices, function(inds) {
  bound_regime(g = minimax_superopt, dat.train = dat.te[-inds,], dat.test = dat.te[inds,])
})
lower_bound.minimax_superopt <-  mean(sapply(results.minimax_superopt, function(x) {x[1]}, simplify = 0))
upper_bound.minimax_superopt <-  mean(sapply(results.minimax_superopt, function(x) {x[2]}, simplify = 0))
minimax_superopt.low.var <- mean(sapply(results.minimax_superopt, function(x) {x[3]}, simplify = 0))
minimax_superopt.upp.var <- mean(sapply(results.minimax_superopt, function(x) {x[4]}, simplify = 0))
minimax_superopt_CI <- c(lower_bound.minimax_superopt - qnorm(0.975)*sqrt(minimax_superopt.low.var/n), upper_bound.minimax_superopt + qnorm(0.975)*sqrt(minimax_superopt.upp.var/n))


#healthcare (superoptimal)
results.healthcare_superopt <- lapply(test.indices, function(inds) {
  bound_regime(g = healthcare_superopt, dat.train = dat.te[-inds,], dat.test = dat.te[inds,])
})
lower_bound.healthcare_superopt <-  mean(sapply(results.healthcare_superopt, function(x) {x[1]}, simplify = 0))
upper_bound.healthcare_superopt <-  mean(sapply(results.healthcare_superopt, function(x) {x[2]}, simplify = 0))
healthcare_superopt.low.var <- mean(sapply(results.healthcare_superopt, function(x) {x[3]}, simplify = 0))
healthcare_superopt.upp.var <- mean(sapply(results.healthcare_superopt, function(x) {x[4]}, simplify = 0))
healthcare_superopt_CI <- c(lower_bound.healthcare_superopt - qnorm(0.975)*sqrt(healthcare_superopt.low.var/n), upper_bound.healthcare_superopt + qnorm(0.975)*sqrt(healthcare_superopt.upp.var/n))

#healthcare with observed baseline
results.healthcare_obs_baseline <- lapply(test.indices, function(inds) {
  bound_regime(g = healthcare_obs_baseline, dat.train = dat.te[-inds,], dat.test = dat.te[inds,])
})
lower_bound.healthcare_obs_baseline <-  mean(sapply(results.healthcare_obs_baseline, function(x) {x[1]}, simplify = 0))
upper_bound.healthcare_obs_baseline <-  mean(sapply(results.healthcare_obs_baseline, function(x) {x[2]}, simplify = 0))
healthcare_obs_baseline.low.var <- mean(sapply(results.healthcare_obs_baseline, function(x) {x[3]}, simplify = 0))
healthcare_obs_baseline.upp.var <- mean(sapply(results.healthcare_obs_baseline, function(x) {x[4]}, simplify = 0))
healthcare_obs_baseline_CI <- c(lower_bound.healthcare_obs_baseline - qnorm(0.975)*sqrt(healthcare_obs_baseline.low.var/n), upper_bound.healthcare_obs_baseline + qnorm(0.975)*sqrt(healthcare_obs_baseline.upp.var/n))

#healthcare with observed baseline (opt version)
results.healthcare_obs_baseline.opt <- lapply(test.indices, function(inds) {
  bound_regime(g = healthcare_obs_baseline.opt, dat.train = dat.te[-inds,], dat.test = dat.te[inds,])
})
lower_bound.healthcare_obs_baseline.opt <-  mean(sapply(results.healthcare_obs_baseline.opt, function(x) {x[1]}, simplify = 0))
upper_bound.healthcare_obs_baseline.opt <-  mean(sapply(results.healthcare_obs_baseline.opt, function(x) {x[2]}, simplify = 0))
healthcare_obs_baseline.opt.low.var <- mean(sapply(results.healthcare_obs_baseline.opt, function(x) {x[3]}, simplify = 0))
healthcare_obs_baseline.opt.upp.var <- mean(sapply(results.healthcare_obs_baseline.opt, function(x) {x[4]}, simplify = 0))
healthcare_obs_baseline_CI.opt <- c(lower_bound.healthcare_obs_baseline.opt - qnorm(0.975)*sqrt(healthcare_obs_baseline.opt.low.var/n), upper_bound.healthcare_obs_baseline.opt + qnorm(0.975)*sqrt(healthcare_obs_baseline.opt.upp.var/n))


#Plots
width = 3
par(mar = c(5,5,1,1), cex.axis = 1.5, cex.lab = 1.5)
plot(NULL, ylim = c(0,1), xlim = c(0,0.95), 
     xlab = "", ylab = expression("E(Y"^"g"*")"), xaxt = 'n')
axis(1, at = c(0.05, 0.15, 0.3, 0.4, 0.55, 0.65, 0.8, 0.9),
     labels = c(expression("g"["opt"]), expression("g"["sup"]), expression("g"["opt"]), expression("g"["sup"]), expression("g"["opt"]), expression("g"["sup"]), expression("g"["opt"]), expression("g"["sup"])))
axis(1, at = c(0.1, 0.35, 0.6, 0.85),
     labels = c("Optimist", "Pessimist", "Opportunist", "Healthcare"), tick = FALSE, line = 1)

abline(h = 0)
abline(h = mean(dat$Y), col = "red")
#Optimist
arrows(x0=0.05, y0=max(0,maximax_CI[1]), x1=0.05, 
       y1=min(1,maximax_CI[2]),
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
arrows(x0=0.05, y0=max(0,lower_bound.maximax_opt), x1=0.05, y1=min(1,upper_bound.maximax_opt), 
       code=3, angle=90, length=0.025, lwd=width, col = 'blue')
arrows(x0=0.15, y0=max(0,maximax_superopt_CI[1]), x1=0.15, 
       y1=min(1,maximax_superopt_CI[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
arrows(x0=0.15, y0=max(0,lower_bound.maximax_superopt), x1=0.15, y1=min(1,upper_bound.maximax_superopt), 
       code=3, angle=90, length=0.025, lwd=width, col = 'blue')
#Pessimist
arrows(x0=0.3, y0=max(0,maximin_CI[1]), x1=0.3, 
       y1=min(1,maximin_CI[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'green')
arrows(x0=0.3, y0=max(0,lower_bound.maximin_opt), x1=0.3, y1=min(1,upper_bound.maximin_opt), 
       code=3, angle=90, length=0.025, lwd=width, col = 'black')
arrows(x0=0.4, y0=max(0,maximin_superopt_CI[1]), x1=0.4, 
       y1=min(1,maximin_superopt_CI[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'green')
arrows(x0=0.4, y0=max(0,lower_bound.maximin_superopt), x1=0.4, y1=min(1,upper_bound.maximin_superopt), 
       code=3, angle=90, length=0.025, lwd=width, col = 'black')
#Opportunist
arrows(x0=0.55, y0=max(0,minimax_CI[1]), x1=0.55, 
       y1=min(1,minimax_CI[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'red4')
arrows(x0=0.55, y0=max(0,lower_bound.minimax_opt), x1=0.55, y1=min(1,upper_bound.minimax_opt), 
       code=3, angle=90, length=0.025, lwd=width, col = 'cyan')
arrows(x0=0.65, y0=max(0,minimax_superopt_CI[1]), x1=0.65, 
       y1=min(1,minimax_superopt_CI[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'red4')
arrows(x0=0.65, y0=max(0,lower_bound.minimax_superopt), x1=0.65, 
       y1=min(1,upper_bound.minimax_superopt), 
       code=3, angle=90, length=0.025, lwd=width, col = 'cyan')
#Healthcare
arrows(x0=0.8, y0=max(0,healthcare_CI[1]), x1=0.8, 
       y1=min(1,healthcare_CI[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'darkgreen')
arrows(x0=0.8, y0=max(0,lower_bound.healthcare_opt), x1=0.8, 
       y1=min(1,upper_bound.healthcare_opt), 
       code=3, angle=90, length=0.025, lwd=width, col = 'azure3')
arrows(x0=0.9, y0=max(0,healthcare_superopt_CI[1]), x1=0.9, 
       y1=min(1,healthcare_superopt_CI[2]), 
       code=3, angle=90, length=0.05, lwd=width, col = 'darkgreen')
arrows(x0=0.9, y0=max(0,lower_bound.healthcare_superopt), x1=0.9, 
       y1=min(1,upper_bound.healthcare_superopt), 
       code=3, angle=90, length=0.025, lwd=width, col = 'azure3')
