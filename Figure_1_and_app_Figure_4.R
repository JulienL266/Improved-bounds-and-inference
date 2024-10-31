library(bpbounds)

is.feasible <- function(pZ, p_rd.ry){ #p_rd.ry is the joint probability distribution
  p_rd <- c(sum(p_rd.ry[1:4]), sum(p_rd.ry[5:8]), sum(p_rd.ry[9:12]), sum(p_rd.ry[13:16]))
  names(p_rd) <- c("0", "1", "2", "3")
  p_a.z <- c(p_rd[1] + p_rd[2], p_rd[3] + p_rd[4], p_rd[1] + p_rd[3], p_rd[2] + p_rd[4] )
  names(p_a.z) <- c("0.0", "1.0", "0.1", "1.1")
  p_ya.z <- c(p_rd.ry["0.0"] + p_rd.ry["0.1"] + p_rd.ry["1.0"] + p_rd.ry["1.1"],
              p_rd.ry["2.0"] + p_rd.ry["2.2"] + p_rd.ry["3.0"] + p_rd.ry["3.2"],
              p_rd.ry["0.2"] + p_rd.ry["0.3"] + p_rd.ry["1.2"] + p_rd.ry["1.3"],
              p_rd.ry["2.1"] + p_rd.ry["2.3"] + p_rd.ry["3.1"] + p_rd.ry["3.3"],
              p_rd.ry["0.0"] + p_rd.ry["0.1"] + p_rd.ry["2.0"] + p_rd.ry["2.1"],
              p_rd.ry["1.0"] + p_rd.ry["1.2"] + p_rd.ry["3.0"] + p_rd.ry["3.2"],
              p_rd.ry["0.2"] + p_rd.ry["0.3"] + p_rd.ry["2.2"] + p_rd.ry["2.3"],
              p_rd.ry["1.1"] + p_rd.ry["1.3"] + p_rd.ry["3.1"] + p_rd.ry["3.3"])
  names(p_ya.z) <- c("00.0", "01.0", "10.0", "11.0",
                     "00.1", "01.1", "10.1", "11.1")
  tabp = as.table(array(
    p_ya.z,
    dim = c(2, 2, 2),
    dimnames = list(
      x = c(0, 1),
      y = c(0, 1),
      z = c(0, 1)
    )
  ))
  bounds <- bpbounds(tabp)
  EY <- (p_ya.z["10.0"] + p_ya.z["11.0"])*(1-pZ) + (p_ya.z["10.1"] + p_ya.z["11.1"])*pZ
  if(max(bounds$p11upp, bounds$p10upp) < EY && bounds$bplb < 0 && bounds$bpub > 0){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
set.seed(2023)
n <- 10000
done <- FALSE
for(pZ in seq(from = 0, to = 1, length.out = n)){
  for(rd1 in seq(from = 0, to = 1, length.out = n)){
    p_rd.ry <- rep(NA,16)
    p_rd.ry[1] <- rd1
    for(j in 2:15){
      if(sum(p_rd.ry, na.rm = TRUE) == 1){
        p_rd.ry[j] <- 0
      }else{
        p_rd.ry[j] <- sample(seq(from = 0, to = 1-sum(p_rd.ry, na.rm = TRUE), length.out = n), 1)
      }
    }
    p_rd.ry[16] <- 1 - sum(p_rd.ry, na.rm = TRUE)
    names(p_rd.ry) <- c("0.0", "0.1", "0.2", "0.3",
                        "1.0", "1.1", "1.2", "1.3",
                        "2.0", "2.1", "2.2", "2.3",
                        "3.0", "3.1", "3.2", "3.3")
    if(is.feasible(pZ, p_rd.ry)){
      print(pZ)
      print(p_rd.ry)
      pZ_ex <- pZ
      p_rd.ry_ex <- p_rd.ry
      done <- TRUE
      break
    }
  }
  if(done){break}
}
##Plot(superoptimal, CATE)
width = 4
#computing bounds
p_rd <- c(sum(p_rd.ry_ex[1:4]), sum(p_rd.ry_ex[5:8]), sum(p_rd.ry_ex[9:12]), sum(p_rd.ry_ex[13:16]))
names(p_rd) <- c("0", "1", "2", "3")
p_a.z <- c(p_rd[1] + p_rd[2], p_rd[3] + p_rd[4], p_rd[1] + p_rd[3], p_rd[2] + p_rd[4] )
names(p_a.z) <- c("0.0", "1.0", "0.1", "1.1")
p_ya.z <- c(p_rd.ry_ex["0.0"] + p_rd.ry_ex["0.1"] + p_rd.ry_ex["1.0"] + p_rd.ry_ex["1.1"],
            p_rd.ry_ex["2.0"] + p_rd.ry_ex["2.2"] + p_rd.ry_ex["3.0"] + p_rd.ry_ex["3.2"],
            p_rd.ry_ex["0.2"] + p_rd.ry_ex["0.3"] + p_rd.ry_ex["1.2"] + p_rd.ry_ex["1.3"],
            p_rd.ry_ex["2.1"] + p_rd.ry_ex["2.3"] + p_rd.ry_ex["3.1"] + p_rd.ry_ex["3.3"],
            p_rd.ry_ex["0.0"] + p_rd.ry_ex["0.1"] + p_rd.ry_ex["2.0"] + p_rd.ry_ex["2.1"],
            p_rd.ry_ex["1.0"] + p_rd.ry_ex["1.2"] + p_rd.ry_ex["3.0"] + p_rd.ry_ex["3.2"],
            p_rd.ry_ex["0.2"] + p_rd.ry_ex["0.3"] + p_rd.ry_ex["2.2"] + p_rd.ry_ex["2.3"],
            p_rd.ry_ex["1.1"] + p_rd.ry_ex["1.3"] + p_rd.ry_ex["3.1"] + p_rd.ry_ex["3.3"])
names(p_ya.z) <- c("00.0", "01.0", "10.0", "11.0",
                   "00.1", "01.1", "10.1", "11.1")
tabp = as.table(array(
  p_ya.z,
  dim = c(2, 2, 2),
  dimnames = list(
    x = c(0, 1),
    y = c(0, 1),
    z = c(0, 1)
  )
))
bounds <- bpbounds(tabp)

#Computing EY and pA
EY <- pZ_ex*sum(p_ya.z[c("10.1", "11.1")]) + (1-pZ_ex)*sum(p_ya.z[c("10.0", "11.0")])
pA <- pZ*p_a.z["1.1"] + (1-pZ)*p_a.z["1.0"]

#Numbers in Figure 4 of the Appendix
print(p_rd.ry_ex)

##Plots for Figure 1 in the main text
par(mar = c(3,4,1,1), cex.axis = 2, cex.lab = 2)
plot(NULL, ylim = c(-1,1), xlim = c(0.2,0.8), 
     xlab = "", ylab = "", xaxt = 'n')
axis(1, at = c(0.25, 0.5, 0.75), labels = c("ATE", "A = 0", "A = 1"))
abline(h = 0)
arrows(x0=0.25, y0=bounds$bplb, x1=0.25, y1=bounds$bpub, 
       code=3, angle=90, length=0.025, lwd=width, col = 'blue')
arrows(x0=0.5, y0=(bounds$p11low - EY)/(1-pA), x1=0.5, y1=(bounds$p11upp - EY)/(1-pA), 
       code=3, angle=90, length=0.025, lwd=width, col = 'red')
arrows(x0=0.75, y0=(-bounds$p10upp + EY)/pA, x1=0.75, y1=(-bounds$p10low + EY)/pA, 
       code=3, angle=90, length=0.025, lwd=width, col = 'red')

