#Construct data
set.seed(2023)
n <- 2000
#n <- 10000
U <- rbinom(n,1, 0.5)
X <- runif(n, min = -2, max = 2)
sigma <- function(x){
  return(1/(1 + exp(-x)))
}
e_t <- function(x){
  return(sigma(0.75*x + 0.5))
}
alpha_t <- function(x, Gamma, t = 1){
  e_t.val <- e_t(x)
  if(t == 0){
    e_t.val <- 1 - e_t.val
  }
  return(1/(Gamma*e_t.val) + 1 - 1/Gamma)
}
beta_t <- function(x, Gamma, t = 1){
  e_t.val <- e_t(x)
  if(t == 0){
    e_t.val <- 1 - e_t.val
  }
  return(Gamma/e_t.val+ 1 - Gamma)
}
Gamma_st <- exp(1)
e <- function(x,u){
  return(u/alpha_t(x,Gamma_st) + (1-u)/beta_t(x, Gamma_st))
}
A <- rep(NA,n)
for(i in 1:n){
  A[i] <- rbinom(1,1, e(X[i], U[i]))
}
epsilon <- rnorm(n)
Y_0 = -X - 1 - 2*sin(2*(-1)*X) - 2*(2*U - 1)*(1 + 0.5*X) + epsilon
Y_1 = X + 1 - 2*sin(2*(+1)*X) - 2*(2*U - 1)*(1 + 0.5*X) + epsilon
Y <- rep(NA, n)
for(i in 1:n){
  if(A[i] == 0){
    Y[i] <- Y_0[i]
  }else{
    Y[i] <- Y_1[i]
  }
}

#Estimation

#re-order data
indices <- order(Y)
Y <- Y[indices]
A <- A[indices]
X <- X[indices]

e_t.model <- glm(A~X, family = "binomial")
alpha_t.estim <- function(t,x, Gamma){
  e_t.estim <- predict(e_t.model, newdata = data.frame(X = x), type = "response")
  if(t == 0){
    e_t.estim <- 1 - e_t.estim
  }
  return(1/(Gamma*e_t.estim) + 1 - 1/Gamma)
}
beta_t.estim <- function(t,x, Gamma){
  e_t.estim <- predict(e_t.model, newdata = data.frame(X = x), type = "response")
  if(t == 0){
    e_t.estim <- 1 - e_t.estim
  }
  return(Gamma/e_t.estim + 1 - Gamma)
}
#choose h
h0 <- bw.bcv(X[which(A == 0)], nb = length(X[which(A == 0)]))

h1 <- bw.bcv(X[which(A == 1)], nb = length(X[which(A == 1)]))

lower.boundary <- (min(X) - max(X))/min(c(h0,h1))
upper.boundary <- (max(X) - min(X))/min(c(h0,h1))
norm.factor <- (pnorm(upper.boundary) - pnorm(lower.boundary))*sqrt(2*pi)
K <- function(u){
  if(u < lower.boundary | u > upper.boundary){
    return(0)
  }else{
    return(exp(-(u^2)/2)/norm.factor) 
  }
}
alpha_tilde <- function(i,t,x,Gamma){
  alpha <- 0
  if(A[i] == t){
    if(t == 0){
      h <- h0
    }else{
      h <- h1
    }
    alpha <- alpha_t.estim(t,X[i], Gamma)*K((X[i] - x)/h)
  }
  return(alpha)
}
beta_tilde <- function(i,t,x,Gamma){
  beta <- 0
  if(A[i] == t){
    if(t == 0){
      h <- h0
    }else{
      h <- h1
    }
    beta <- beta_t.estim(t,X[i], Gamma)*K((X[i] - x)/h)
  }
  return(beta)
}
lamba_up <- function(t,k,x,Gamma){
  s1 <- 0
  s2 <- 0
  s3 <- 0
  s4 <- 0
  if(k >= 1){
    for(i in 1:k){
      s1 <- s1 + alpha_tilde(i,t,x,Gamma)*Y[i]
      s3 <- s3 + alpha_tilde(i,t,x,Gamma)
    }
  }
  if(k + 1 <= n){
    for(i in (k + 1):n){
      s2 <- s2 + beta_tilde(i,t,x,Gamma)*Y[i]
      s4 <- s4 + beta_tilde(i,t,x,Gamma) 
    }
  }
  return((s1 + s2)/(s3 + s4))
}
lamba_down <- function(t,k,x,Gamma){
  s1 <- 0
  s2 <- 0
  s3 <- 0
  s4 <- 0
  if(k >= 1){
    for(i in 1:k){
      s1 <- s1 + beta_tilde(i,t,x,Gamma)*Y[i]
      s3 <- s3 + beta_tilde(i,t,x,Gamma)
    }
  }
  if(k + 1 <= n){
    for(i in (k + 1):n){
      s2 <- s2 + alpha_tilde(i,t,x,Gamma)*Y[i]
      s4 <- s4 + alpha_tilde(i,t,x,Gamma) 
    }
  }
  return((s1 + s2)/(s3 + s4))
}
k_H <- function(t,x, Gamma){
  for(k in 0:(n-1)){
    if(lamba_up(t,k,x,Gamma) >lamba_up(t,k + 1,x,Gamma)){ #strict inequality in new Kallus paper
      return(k)
    }
  }
  return(n)
}
k_L <- function(t,x, Gamma){
  for(k in 0:(n-1)){
    if(lamba_down(t,k,x,Gamma) < lamba_down(t,k + 1,x,Gamma)){ #strict inequality in new Kallus paper
      return(k)
    }
  }
  return(n)
}
mu_up <- function(t,x, Gamma){
  return(lamba_up(t,k_H(t,x, Gamma), x, Gamma))
}
mu_down <- function(t,x, Gamma){
  return(lamba_down(t,k_L(t,x, Gamma), x, Gamma))
}

#Inference

#Gamma <- exp(0.5)
Gamma <- exp(1)
#Gamma <- exp(1.5)

#choose x
EY.A.X <- function(a,x){
  term_U0 <- x + 1 - 2*sin(2*x)
  if(a == 0){
    term_U0 <- - term_U0
  }
  term_U0 <- term_U0 - 2*(-1)*(1 + 0.5*x)
  if(a == 0){
    term_U0 <- term_U0*(1-e(x,0))*0.5/(0.5*(1-e(x,0) + 1-e(x,1)))
  }else{
    term_U0 <- term_U0*e(x,0)*0.5/(0.5*(e(x,0) + e(x,1)))
  }
  
  term_U1 <- x + 1 - 2*sin(2*x)
  if(a == 0){
    term_U1 <- - term_U1
  }
  term_U1 <- term_U1 - 2*(1)*(1 + 0.5*x)
  if(a == 0){
    term_U0 <- term_U0*(1-e(x,1))*0.5/(0.5*(1-e(x,0) + 1-e(x,1)))
  }else{
    term_U0 <- term_U0*e(x,1)*0.5/(0.5*(e(x,0) + e(x,1)))
  }
  return(0.5*(term_U0 + term_U1))
}
EY.X <- function(x){
  return(EY.A.X(0,x)*0.5*(1-e(x,0) + 1-e(x,1)) + EY.A.X(1,x)*0.5*(e(x,0) + e(x,1)))
}
EY1.X <- function(x){
  return(x + 1 - 2*sin(2*x))
}
EY0.X <- function(x){
  return(-x - 1 + 2*sin(2*x))
}
#pick x with reversal of effects
for(x_test in seq(-2,2,0.1)){
  if(EY1.X(x_test) > EY.X(x_test)){
    if(EY0.X(x_test) > EY.X(x_test)){
      x <- x_test
    }
  }
}

Bounds.EY1 <- c(mu_down(1,x,Gamma), mu_up(1,x,Gamma))
Bounds.EY0 <- c(mu_down(0,x,Gamma), mu_up(0,x,Gamma))
Bounds.L.CATE <- c(Bounds.EY1[1] - Bounds.EY0[2], Bounds.EY1[2] - Bounds.EY0[1])

#(L,A)-optimal regime bounds(well-specified model version)
e_estim <- function(x,u, Gamma){
  return(u/alpha_t.estim(1,x,Gamma) + (1-u)/beta_t.estim(1,x, Gamma))
}
EY.A.X_estim <- function(a,x, Gamma){
  term_U0 <- x + 1 - 2*sin(2*x)
  if(a == 0){
    term_U0 <- - term_U0
  }
  term_U0 <- term_U0 - 2*(-1)*(1 + 0.5*x)
  if(a == 0){
    term_U0 <- term_U0*(1-e_estim(x,0, Gamma))*0.5/(0.5*(1-e_estim(x,0, Gamma) + 1-e_estim(x,1, Gamma)))
  }else{
    term_U0 <- term_U0*e_estim(x,0, Gamma)*0.5/(0.5*(e_estim(x,0, Gamma) + e_estim(x,1, Gamma)))
  }
  
  term_U1 <- x + 1 - 2*sin(2*x)
  if(a == 0){
    term_U1 <- - term_U1
  }
  term_U1 <- term_U1 - 2*(1)*(1 + 0.5*x)
  if(a == 0){
    term_U0 <- term_U0*(1-e_estim(x,1, Gamma))*0.5/(0.5*(1-e_estim(x,0, Gamma) + 1-e_estim(x,1, Gamma)))
  }else{
    term_U0 <- term_U0*e_estim(x,1, Gamma)*0.5/(0.5*(e_estim(x,0, Gamma) + e_estim(x,1, Gamma)))
  }
  return(0.5*(term_U0 + term_U1))
}
EY.X_estim <- function(x, Gamma){
  return(EY.A.X_estim(0,x, Gamma)*0.5*(1-e_estim(x,0, Gamma) + 1-e_estim(x,1, Gamma)) + EY.A.X_estim(1,x, Gamma)*0.5*(e_estim(x,0, Gamma) + e_estim(x,1, Gamma)))
}
Bounds.L1.CATE <- c(EY.X_estim(x, Gamma) - Bounds.EY0[2], EY.X_estim(x, Gamma) - Bounds.EY0[1])/(1-predict(e_t.model, newdata = data.frame(X = x), type = "response"))
Bounds.L0.CATE <- c(-EY.X_estim(x, Gamma) + Bounds.EY1[1], -EY.X_estim(x, Gamma) + Bounds.EY1[2])/(predict(e_t.model, newdata = data.frame(X = x), type = "response"))

#(L,A)-optimal regime bounds
library(rpart)
np.model <- rpart(Y~X)
np.EY.X_estim <- predict(np.model, newdata = data.frame(X = x), type = "vector")
sp.Bounds.L1.CATE <- c(np.EY.X_estim - Bounds.EY0[2], np.EY.X_estim - Bounds.EY0[1])/(1-predict(e_t.model, newdata = data.frame(X = x), type = "response"))
sp.Bounds.L0.CATE <- c(-np.EY.X_estim + Bounds.EY1[1], -np.EY.X_estim + Bounds.EY1[2])/(predict(e_t.model, newdata = data.frame(X = x), type = "response"))

#True CATEs
True.L.CATE <- EY1.X(x)- EY0.X(x)
True.EY.X <- EY.X(x) 
True.L0.CATE <- (EY1.X(x) - EY.X(x))/e_t(x)
True.L1.CATE <- (EY.X(x) - EY0.X(x))/(1 - e_t(x)) 
#superopt bounds correctly identify superoptimal regime


#Combined plot
## L-CATE
width = 3
par(mfrow = c(1,3) , cex=1, cex.axis = 1.5, cex.lab = 1.5)
plot(NULL, ylim = c(-4,4), xlim = c(0.2,0.8),xlab  = "", ylab = "L-CATE", xaxt = 'n')
axis(1, at = c(0.35, 0.5, 0.65),
     labels = c("exp(0.5)", "exp(1)", "exp(1.5)"))
abline(h = 0)
abline( h=True.L.CATE, col = 'blue')
# Gamma = exp(1)
arrows(x0=0.5, y0=Bounds.L.CATE[1], x1=0.5, 
       y1=Bounds.L.CATE[2],
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
# Gamma = exp(0.5)
Gamma = exp(0.5)
Bounds.EY1 <- c(mu_down(1,x,Gamma), mu_up(1,x,Gamma))
Bounds.EY0 <- c(mu_down(0,x,Gamma), mu_up(0,x,Gamma))
Bounds.L.CATE.L.G5 <- c(Bounds.EY1[1] - Bounds.EY0[2], Bounds.EY1[2] - Bounds.EY0[1])
arrows(x0=0.35, y0=Bounds.L.CATE.L.G5[1], x1=0.35, 
       y1=Bounds.L.CATE.L.G5[2],
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
# Gamma = exp(1.5)
Gamma = exp(1.5)
Bounds.EY1 <- c(mu_down(1,x,Gamma), mu_up(1,x,Gamma))
Bounds.EY0 <- c(mu_down(0,x,Gamma), mu_up(0,x,Gamma))
Bounds.L.CATE.L.G15 <- c(Bounds.EY1[1] - Bounds.EY0[2], Bounds.EY1[2] - Bounds.EY0[1])
arrows(x0=0.65, y0=Bounds.L.CATE.L.G15[1], x1=0.65, 
       y1=Bounds.L.CATE.L.G15[2],
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
## (L,0)-CATE
plot(NULL, ylim = c(-4,4), xlim = c(0.2,0.8), xlab = "", ylab = "(L,0)-CATE", xaxt = 'n')
axis(1, at = c(0.35, 0.5, 0.65),
     labels = c("exp(0.5)", "exp(1)", "exp(1.5)"))
abline(h = 0)
abline( h=True.L0.CATE, col = 'blue')
# Gamma = exp(1)
arrows(x0=0.5, y0=Bounds.L0.CATE[1], x1=0.5, 
       y1=Bounds.L0.CATE[2],
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
# Gamma = exp(0.5)
Gamma = exp(0.5)
Bounds.EY1 <- c(mu_down(1,x,Gamma), mu_up(1,x,Gamma))
Bounds.L0.CATE.G5 <- c(-EY.X_estim(x, Gamma) + Bounds.EY1[1], -EY.X_estim(x, Gamma) + Bounds.EY1[2])/(predict(e_t.model, newdata = data.frame(X = x), type = "response"))
arrows(x0=0.35, y0=Bounds.L0.CATE.G5[1], x1=0.35, 
       y1=Bounds.L0.CATE.G5[2],
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
# Gamma = exp(1.5)
Gamma = exp(1.5)
Bounds.EY1 <- c(mu_down(1,x,Gamma), mu_up(1,x,Gamma))
Bounds.L0.CATE.G15 <- c(-EY.X_estim(x, Gamma) + Bounds.EY1[1], -EY.X_estim(x, Gamma) + Bounds.EY1[2])/(predict(e_t.model, newdata = data.frame(X = x), type = "response"))
arrows(x0=0.65, y0=Bounds.L0.CATE.G15[1], x1=0.65, 
       y1=Bounds.L0.CATE.G15[2],
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
## (L,1)-CATE
plot(NULL, ylim = c(-18,18), xlim = c(0.2,0.8), xlab = "", ylab = "(L,1)-CATE", xaxt = 'n')
axis(1, at = c(0.35, 0.5, 0.65),
     labels = c("exp(0.5)", "exp(1)", "exp(1.5)"))
abline(h = 0)
abline( h=True.L1.CATE, col = 'blue')
# Gamma = exp(1)
arrows(x0=0.5, y0=Bounds.L1.CATE[1], x1=0.5, 
       y1=Bounds.L1.CATE[2],
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
# Gamma = exp(0.5)
Gamma = exp(0.5)
Bounds.EY0 <- c(mu_down(0,x,Gamma), mu_up(0,x,Gamma))
Bounds.L1.CATE.G5 <- c(EY.X_estim(x, Gamma) - Bounds.EY0[2], EY.X_estim(x, Gamma) - Bounds.EY0[1])/(1-predict(e_t.model, newdata = data.frame(X = x), type = "response"))
arrows(x0=0.35, y0=Bounds.L1.CATE.G5[1], x1=0.35, 
       y1=Bounds.L1.CATE.G5[2],
       code=3, angle=90, length=0.05, lwd=width, col = 'red')
# Gamma = exp(1.5)
Gamma = exp(1.5)
Bounds.EY0 <- c(mu_down(0,x,Gamma), mu_up(0,x,Gamma))
Bounds.L1.CATE.G15 <- c(EY.X_estim(x, Gamma) - Bounds.EY0[2], EY.X_estim(x, Gamma) - Bounds.EY0[1])/(1-predict(e_t.model, newdata = data.frame(X = x), type = "response"))
arrows(x0=0.65, y0=Bounds.L1.CATE.G15[1], x1=0.65, 
       y1=Bounds.L1.CATE.G15[2],
       code=3, angle=90, length=0.05, lwd=width, col = 'red')