## E-STEP FUNCTIONS

## prob of ND given Z
cpNDz <- function(z, pyfit){
  predict(pyfit, newdata=data.frame(gavg=z), type="response")
}

## joint prob of ND,Z
pNDZ <- function(z, mu, s2, pyfit){
  cpNDz(z, pyfit)*dnorm(z, mu, sqrt(s2))
}

## marginal prob of ND
pND <- function(mu, s2, pyfit){
  f <- function(z) pNDZ(z, mu, s2, pyfit)
  integrate(f, lower=mu-10*sqrt(s2), upper=mu+10*sqrt(s2))$value
}

## E[Z]
EZ <- function(mu, s2, pyfit){
  f <- function(z) z*pNDZ(z, mu, s2, pyfit)/pND(mu, s2, pyfit)
  integrate(f, lower=mu-10*sqrt(s2), upper=mu+10*sqrt(s2))$value
}

## E[Z^2]
EZ2 <- function(mu, s2, pyfit){
  f <- function(z) (z^2)*pNDZ(z, mu, s2, pyfit)/pND(mu, s2, pyfit)
  integrate(f, lower=mu-10*sqrt(s2), upper=mu+10*sqrt(s2))$value
}

#######################

## M-STEP FUNCTIONS

updateS2 <- function(Ct, thetaVec, dj, ez, ez2, i.nd, ngene, nperts){
  s2Vec <- vector(length=length(Ct))
  s2Mat <- matrix(nrow=max(ngene), ncol=max(nperts))
  for(i in 1:max(ngene)){
    for(j in 1:max(nperts)){
      ind <- which(ngene == i)
      indD <- intersect(ind, which(!i.nd))
      indND <- intersect(ind, which(i.nd))
      p1 <- sum((Ct[indD]-(thetaVec[indD]+dj[indD]))^2)
      p2 <- sum(ez2[indND]-2*ez[indND]*(thetaVec[indND]+dj[indND])
                +((thetaVec[indND]+dj[indND])^2))
      s2Vec[ind] <- s2Mat[i,j] <- (p1+p2)/length(ind)
    }
  }
  return(list(s2Vec, s2Mat))
}

updateTheta <- function(Ct, ez, dj, i.nd, ngene, nperts){
  thetaVec <- vector(length=length(Ct))
  thetaMat <- matrix(nrow=max(ngene), ncol=max(nperts))
  for(i in 1:max(ngene)){
    for(j in 1:max(nperts)){
      ind <- intersect(which(ngene==i), which(nperts==j))
      indD <- intersect(ind, which(!i.nd))
      indND <- intersect(ind, which(i.nd))
      num <- sum(Ct[indD]-dj[indD])+sum(ez[indND]-dj[indND])
      denom <- length(ind)
      thetaVec[ind] <- thetaMat[i,j] <- num/denom
    }
  }
  return(list(thetaVec, thetaMat))
}

logLik <- function(Ct, ez, ez2, s2Vec, thetaVec, dj, i.nd){
  p1 <- -sum(log(2*pi*s2Vec)/2)
  p2 <- -sum(((Ct[which(!i.nd)]-(thetaVec[which(!i.nd)]+dj[which(!i.nd)]))^2)
             / (2*s2Vec[which(!i.nd)]))
  p3 <- -sum((ez2[which(i.nd)]-2*ez[which(i.nd)]
              *(thetaVec[which(i.nd)]+dj[which(i.nd)])+
              (thetaVec[which(i.nd)]+dj[which(i.nd)])^2)
             / (2*s2Vec[which(i.nd)]))
  p1+p2+p3
}
