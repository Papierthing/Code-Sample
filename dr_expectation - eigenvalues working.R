#This was the format of the code used in the Financial Mathematics Team Challenge. It involved doing an eignedecomposition of 
#the variance-covariance matrix of stock returns on the Dow-Jones index. Small eigenvalues and their associated vectors were
#dropped, and the rest used to implement a technique known as DR-expectation


library(MASS)
#library(mvtnorm)
#library(GLDEX)
library(stats)

MyData <- read.csv("F:/FMTC/djia constituent log returns.csv", header = TRUE, sep = ";")

daysused <- 300
daystested <- 252
startindex <- daystested + 1
endindex <- daystested + daysused 
n <- 30
kprime <- 8
k <- 22
lambda <- 30
eigensused <- 15

thirtystockdata <- MyData[startindex:endindex,2:(n+1)]

#muest <- colMeans(thirtystockdata)
Vest <- var(thirtystockdata)
ones <- integer(n)+1

eigenest <- eigen(Vest, symmetric = TRUE)
eigensused <- 0
while (eigensused < 20)
{
  if (eigenest$values[eigensused+1] > 1e-6)
  {
    eigensused = eigensused + 1
  }
  else
  {
    break
  }
}
eigenvectors <- eigenest$vectors[,1:eigensused]
eigenest <- eigenest$values[1:eigensused]
muest<-rowMeans(ginv(eigenvectors)%*%as.matrix(t(thirtystockdata)))

estimate <- c(muest*1e4,eigenest*1e4+1)

index1 <- eigensused+1
index2 <- length(estimate)

LL <- function(theta)
{
  mu <- as.vector(theta[1:eigensused])*1e-4
  eigenvalues <- as.vector(theta[index1:index2])*1e-4
  #sigma <- eigenvectors%*%diag(eigenvalues)%*%t(eigenvectors)
  sigmainv <- eigenvectors%*%diag(1/eigenvalues)%*%t(eigenvectors)
  mu<-eigenvectors%*%mu
  #datanorm=as.matrix(thirtystockdata-mu)%*%sigmainv%*%t(as.matrix(thirtystockdata-mu))
  datanorm = mahalanobis(thirtystockdata,mu,sigmainv, inverted = T)
  (sum(datanorm)/2+sum(log(eigenvalues))*daysused/2)
}

mleest = optim(estimate, LL, method = "BFGS", hessian = TRUE)
#thetahat <- mleest$par
thetahat <- mleest$par*1e-4
muhat <- thetahat[1:eigensused]
muhat<-eigenvectors%*%muhat

#muhat <- thetahat[1:n]*1e-4
eigenvalueshat <- thetahat[index1:index2]
sigmahat <- eigenvectors%*%diag(eigenvalueshat)%*%t(eigenvectors)

#eigenvalueshat <- thetahat[index1:index2]*1e-4


Vinv <- ginv(sigmahat)
gamma <- (2*lambda-t(muhat)%*%Vinv%*%ones)/(t(ones)%*%Vinv%*%(ones))
pistar = 1/(2*lambda)*Vinv%*%(muhat+c(gamma)*(ones))

hes <- mleest$hessian
#suplik <- mleest$value
#library(ggplot2)
#library(reshape2)
#ggplot(data=melt(log(abs(hes))),  aes(x=Var1, y=Var2, fill=value)) + geom_tile()

Rpenalty <- function(theta)
{
  lldif <- 1/2*t(theta-thetahat)%*%hes%*%(theta-thetahat)*1e8
  (1/k*lldif)^kprime
}

#Rpenalty <- function(theta)
#{
#  lldif <- LL(theta)-suplik
#  (1/k*lldif)^kprime
#}


#epsilon = function(pi1, tole=1e-4)
#{
#  pi <- c(pi1, 1-sum(pi1))
#  innerep = function(theta)
#  {
#    mu <- theta[1:n]*1e-4
#    eigenvalues <- theta[index1:index2]*1e-4
#    sigma <- eigenvectors%*%diag(eigenvalues)%*%t(eigenvectors)
#    mu<-eigenvectors%*%t(eigenvectors)%*%mu
#    
#  }
#  
#  est = optim(estimate, innerep,method = "BFGS", control=list(maxit=30,reltol=tole,trace=2))
#  est$value
#}


eps2=function(theta){
  mu <- theta[1:eigensused]*1e-4
  eigenvalues <- theta[index1:index2]*1e-4
  sigma <- eigenvectors%*%diag(eigenvalues)%*%t(eigenvectors)
  mu<-eigenvectors%*%mu
  Vinv <- eigenvectors%*%diag(1/eigenvalues)%*%t(eigenvectors)
  gamma <- (2*lambda-t(mu)%*%Vinv%*%ones)/(t(ones)%*%Vinv%*%(ones))
  pistar = 1/(2*lambda)*Vinv%*%(mu+c(gamma)*(ones))
  t(pistar)%*%mu - lambda*t(pistar)%*%sigma%*%pistar+Rpenalty(theta*1e-4)
}

thetaminimax=optim(estimate, eps2,method = "BFGS", control=list(maxit = 20000, trace=2))

thetahat <- thetaminimax$par*1e-4
muhat <- thetahat[1:eigensused]
muhat<-eigenvectors%*%muhat

#muhat <- thetahat[1:n]*1e-4
eigenvalueshat <- thetahat[index1:index2]
sigmahat <- eigenvectors%*%diag(eigenvalueshat)%*%t(eigenvectors)

#eigenvalueshat <- thetahat[index1:index2]*1e-4


Vinv <- ginv(sigmahat)
gamma <- (2*lambda-t(muhat)%*%Vinv%*%ones)/(t(ones)%*%Vinv%*%(ones))
piminimax = 1/(2*lambda)*Vinv%*%(muhat+c(gamma)*(ones))

#guess <- pistar[1:(n-1)]
#guess <- integer(29)+1/n
#piest = optim(guess,epsilon, method = "BFGS", control = list(fnscale = -1,maxit=20, reltol=1e-4, trace=6))
#piweights <- c(piest$par, 1-sum(piest$par))
evenpi <- integer(n)+1/n


ClosingPrices <- read.csv("F:/FMTC/DJIA constituent prices.csv", header = TRUE, sep = ";")
testdata <- ClosingPrices[2:daystested,2:(n+1)]
marketpi <- c(t(ClosingPrices[startindex,2:(n+1)]/sum(ClosingPrices[startindex,2:31])))

  

testdata2 <- MyData[daystested:2,2:(n+1)]
#testdata2 <- t(testdata2)

#value1 <- piweights%*%t(testdata)
#value1 <- value1/value1[199]
#value2 <- t(pistar)%*%t(testdata)
#value2 <- value2/value2[199]

#hist(value1-value2)

returns1 <- t(piminimax)%*%t(testdata2)
returns2 <- t(pistar)%*%t(testdata2)
returns3 <- evenpi%*%t(testdata2)
returns4 <- marketpi%*%t(testdata2)

plot(cumsum((returns1[1:daystested-1])),,'l')
lines(cumsum(returns2[1:daystested-1]), col='red')
lines(cumsum(returns3[1:daystested-1]), col='blue')
lines(cumsum(returns4[1:daystested-1]), col='green')

returns <- t(rbind(returns1,returns2,returns3,returns4))
filename <- c(k, kprime, lambda)
filename = paste(as.character(filename), collapse = ',')
write.table(returns, file = filename)
returnstest <- read.csv(filename, header = T, sep = " ")
