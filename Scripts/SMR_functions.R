library(spatstat)
library(tidyverse)
library(coda)
library(rgdal)
library(maptools)
library(parallel)

# Convert WGS to projected CRS ----
WGS84toCRS <- function(x,y){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84") 
  res <- spTransform(xy, crs)
  return(as.data.frame(res))
}

# SMR function ----
SMR_dens<- function(n, X, y, M, mmax, obsmod = c("pois", "bern"), niters, 
                    npics, region, sigma.prior, lam0.prior, inits, delta) 
{
  obsmod <- match.arg(obsmod)
  J <- nrow(n)
  K <- ncol(n)
  S <- inits$S
  D <- calc.dist(S, X) # NB these are SQUARED distances
  A<- area.owin(region)
  sigma <- inits$sigma
  lam0 <- inits$lam0
  ## NB: HHN would be specified with -d^2 in secr, but D = d^2 here
  lam <- lam0 * exp(-D/(2 * sigma * sigma)) 
  nObs <- nrow(y)
  Y <- array(0, c(M + mmax, J, K))
  Y[1:nObs, , ] <- y
  marked <- rep(FALSE, M + mmax)
  marked[1:mmax] <- TRUE
  psi <- inits$psi
  psim <- inits$psim
  z <- rbinom(M + mmax, 1, psi)
  z[1:nObs] <- 1
  for (j in 1:J) {
    for (k in 1:K) {
      if(is.na(n[j,k])) {
        Y[!marked, j, k] <- NA #if missing
        next
      }
      if (n[j, k] == 0) {
        Y[!marked, j, k] <- 0
        next
      }
      nUnknown <- n[j, k]
      probs <- lam[!marked, j] * z[!marked]
      probs <- probs/sum(probs)
      if (identical(obsmod, "pois")) 
        Y[!marked, j, k] <- rmultinom(1, nUnknown, probs)
      else if (identical(obsmod, "bern")) {
        Y[!marked, j, k] <- 0
        guys <- sample((mmax + 1):(M + mmax), nUnknown, prob = probs)
        Y[guys, j, k] <- 1
      }
    }
  }
  cr <- rep(1, M + mmax)
  if (missing(npics)) {
    crat <- 1
  }
  else {
    crat <- npics[1]/npics[2]
  }
  cr[marked] <- crat
  
  out <- matrix(NA, nrow = niters, ncol = 8)
  colnames(out) <- c("sigma", "lam0", "c", "psi", "psim", "m", "N", "D")
  cat("\nstarting values =", c(sigma, lam0, crat, psi, psim, 
                               sum(z[marked]), sum(z)), sum(z)/A, "\n\n")
  for (iter in 1:niters) {
    if (iter%%100 == 0) {
      cat("iter", iter, format(Sys.time(), "%H:%M:%S"), 
          "\n")
      cat("   current =", out[iter - 1, ], "\n")
    }
    if (identical(obsmod, "pois")) {
      ll <- sum(dpois(Y, lam * cr * z, log = TRUE), na.rm=TRUE)
    }
    else if (identical(obsmod, "bern")) {
      ll <- sum(dbinom(Y, 1, lam * cr * z, log = TRUE), na.rm=TRUE)
    }
    if (!missing(npics)) {
      crat <- rbeta(1, 1 + npics[1], 1 + npics[2] - npics[1])
      cr[marked] <- crat
    }
    # Update sigma
    sigma.cand <- rnorm(1, sigma, delta[1])
    sig.prior.curr<- prior.density(sigma, sigma.prior,logdens=T)
    sig.prior.cand<- prior.density(sigma.cand, sigma.prior,logdens=T)
    
    if (sigma.cand > 0) {
      lam.cand <- lam0 * exp(-D/(2 * sigma.cand * sigma.cand))
      if (identical(obsmod, "pois")) {
        ll <- sum(dpois(Y, lam * cr * z, log = TRUE), na.rm=TRUE)
        llcand <- sum(dpois(Y, lam.cand * cr * z, log = TRUE), na.rm=TRUE)
      }
      else if (identical(obsmod, "bern")) {
        ll <- sum(dbinom(Y, 1, lam * cr * z, log = TRUE), na.rm=TRUE)
        llcand <- sum(dbinom(Y, 1, lam.cand * cr * z, log = TRUE), na.rm=TRUE)
      }
      R1<- exp((llcand + sig.prior.cand) - (ll + sig.prior.curr))
      if(is.na(R1)) R1<- -1   
      if (runif(1) < R1) {
        ll <- llcand
        lam <- lam.cand
        sigma <- sigma.cand
      }
    }
    # Update lam0 
    lam0.cand <- rnorm(1, lam0, delta[2])
    lam0.prior.curr<- prior.density(lam0, lam0.prior,logdens=T)
    lam0.prior.cand<- prior.density(lam0.cand, lam0.prior,logdens=T)
    if(is.finite(lam0.prior.cand)) {  # make sure proposal has prior support
      lam.cand <- lam0.cand * exp(-D/(2 * sigma * sigma))
      if (identical(obsmod, "pois")) 
        llcand <- sum(dpois(Y, lam.cand * cr * z, log = TRUE), na.rm=TRUE)
      else if (identical(obsmod, "bern")) 
        llcand <- sum(dbinom(Y, 1, lam.cand * cr * z, log = TRUE), na.rm=TRUE)
      R2<- exp(llcand - ll)
      if(is.na(R2)) R2<- -1  
      if (runif(1) < R2) {
        ll <- llcand
        lam0 <- lam0.cand
        lam <- lam.cand
      }
    }
    # Update z for marked
    zUpsm <- zUps <- 0
    for (i in (nObs + 1):mmax) {
      zcand <- ifelse(z[i] == 0, 1, 0)
      if (identical(obsmod, "pois")) {
        llz <- sum(dpois(Y[i, , ], lam[i, ] * cr[i] * z[i], log = TRUE), na.rm=TRUE)
        llcandz <- sum(dpois(Y[i, , ], lam[i, ] * cr[i] * zcand, log = TRUE), na.rm=TRUE)
      }
      else if (identical(obsmod, "bern")) {
        llz <- sum(dbinom(Y[i, , ], 1, lam[i, ] * cr[i] * z[i], log = TRUE), na.rm=TRUE)
        llcandz <- sum(dbinom(Y[i, , ], 1, lam[i, ] * cr[i] * zcand, log = TRUE), na.rm=TRUE)
      }
      prior <- dbinom(z[i], 1, psim, log = TRUE)
      prior.cand <- dbinom(zcand, 1, psim, log = TRUE)
      if (runif(1) < exp((llcandz + prior.cand) - (llz + prior))) {
        z[i] <- zcand
        zUpsm <- zUpsm + 1
      }
    }
    # Update z for unmarked
    seen <- apply(Y > 0, 1, any, na.rm=TRUE)
    for (i in (mmax + 1):(M + mmax)) {
      if (seen[i]) 
        next
      zcand <- ifelse(z[i] == 0, 1, 0)
      if (identical(obsmod, "pois")) {
        ll <- sum(dpois(Y[i, , ], lam[i, ] * z[i], log = TRUE), na.rm=TRUE)
        llcand <- sum(dpois(Y[i, , ], lam[i, ] * zcand, log = TRUE), na.rm=TRUE)
      }
      else if (identical(obsmod, "bern")) {
        ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * z[i], log = TRUE), na.rm=TRUE)
        llcand <- sum(dbinom(Y[i, , ], 1, lam[i, ] * zcand, log = TRUE), na.rm=TRUE)
      }
      prior <- dbinom(z[i], 1, psi, log = TRUE)
      prior.cand <- dbinom(zcand, 1, psi, log = TRUE)
      rat <- (llcand + prior.cand) - (ll + prior)
      if (runif(1) < exp(rat)) {
        z[i] <- zcand
        zUps <- zUps + 1
      }
    }
    for (j in 1:J) {
      zip <- lam[!marked, j] * z[!marked]
      for (k in 1:K) {
        if (is.na(n[j, k])) {
          Y[!marked, j, k] <- NA
          next
        }
        if (n[j, k] == 0) {
          Y[!marked, j, k] <- 0
          next
        }
        nUnknown <- n[j, k]
        probs <- zip/sum(zip)
        if (identical(obsmod, "pois")) 
          Y[!marked, j, k] <- rmultinom(1, nUnknown, probs)
        else if (identical(obsmod, "bern")) {
          Y[!marked, j, k] <- 0
          guy <- sample((mmax + 1):(M + mmax), nUnknown, prob = probs)
          Y[guy, j, k] <- 1
        }
      }
    }
    # Gibbs sample psi for marked and unmarked 
    psim <- rbeta(1, 1 + sum(z[marked]), 1 + mmax - sum(z[marked]))
    psi <- rbeta(1, 1 + sum(z[!marked]), 1 + M - sum(z[!marked]))
    
    # Udate HR locations
    Sups <- 0
    for (i in 1:(M + mmax)) {
      Scand <- c(rnorm(1, S[i, 1], delta[3]), rnorm(1, S[i, 2], delta[3]))
      inregion<- inside.owin(Scand[1],Scand[2], region)
      if(!inregion)
        next
      dtmp <- (Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2 # SQUARED distance to each cam
      lam.cand <- lam0 * exp(-dtmp/(2 * sigma * sigma))
      if (identical(obsmod, "pois")) {
        ll <- sum(dpois(Y[i, , ], lam[i, ] * cr[i] * z[i], log = TRUE), na.rm=TRUE)
        llcand <- sum(dpois(Y[i, , ], lam.cand * cr[i] * z[i], log = TRUE), na.rm=TRUE)
      }
      else if (identical(obsmod, "bern")) {
        ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * cr[i] * z[i], log = TRUE), na.rm=TRUE)
        llcand <- sum(dbinom(Y[i, , ], 1, lam.cand * cr[i] * z[i], log = TRUE), na.rm=TRUE)
      }
      if (runif(1) < exp(llcand - ll)) {
        ll <- llcand
        S[i, ] <- Scand
        lam[i, ] <- lam.cand
        D[i, ] <- dtmp
        Sups <- Sups + 1
      }
    }
    if (iter%%100 == 0) {
      cat("   Acceptance rates\n")
      cat("     z =", zUps/M, "\n")
      cat("     zm =", zUpsm/mmax, "\n")
      cat("     S =", Sups/(M + mmax), "\n")
    }
    out[iter, ] <- c(sigma, lam0, crat, psi, psim, sum(z[marked]), sum(z), sum(z)/A)
  }
  return(out)
}

if(1==2){
# Build a capthist file from a matrix of detections ----
build.capt<- function(history, ID, sessid) {
  # cells of history with multiple ID must be separated by ","
  n<- length(ID)
  J<- dim(history)[1]
  K<- dim(history)[2]
  na.ind<- which(is.na(history), arr.ind=TRUE)
  chist<- list()
  ind<- 1
  for(i in 1:n) {
    id<- ID[i]
    for(j in 1:J) {
      for(k in 1:K) {
        ss<- strsplit(as.character(history[j,k]),",")[[1]]
        if(any(grepl(paste0("^",id,"$"),ss))) {
          chist[[ind]]<- data.frame(sessid,id,k,j)
          ind<- ind + 1
        }
      }
    }
  }
  chist.df<- do.call('rbind', chist)
  names(chist.df)<- c("Session","ID","Occasion","Detector")
  list(chist=chist.df,na.ind=na.ind)
}

# Build a unmarked sighting file from a matrix of detections ----
build.sight<- function(history, ID, sessid) {
  # cells of history with multiple ID must be separated by ","
  n<- length(ID)
  J<- dim(history)[1]
  K<- dim(history)[2]
  Yu<- matrix(0, J, K)
  na.mat<- is.na(history)
  Yu[na.mat]<- NA
  for(i in 1:n){
    id<- ID[i]
    for(j in 1:J) {
      for(k in 1:K) {
        ss<- strsplit(as.character(history[j,k]),",")[[1]]
        if(any(grepl(paste0("^",id,"$"),ss))) Yu[j,k]<- Yu[j,k] + 1
      }
    }
  }
  Yu
}


# Spread TrapID dataframe into detection array ----
capt.array<- function(ch, nocc, ntraps) {
  # ch is output from build.capt()
  trapid<- ch$chist
  id<- sort(unique(trapid$ID))
  n<- length(id)
  Yk<- array(0, c(n, ntraps, nocc))
  Yk[, ch$na.ind[,1], ch$na.ind[,2]]<- NA # Fill in missing cells
  for(i in 1:n) {
    tmpid<- trapid[trapid$ID == id[i],]
    K<- nrow(tmpid)
    for(k in 1:K) {
      Yk[i, tmpid$Detector[k], tmpid$Occasion[k]]<- 1
    }}
  Yk
}
}

# Calculate distance ----
calc.dist<- function(A, D) {
  #calculates distance matrix between devices and HR centres
  dx<- outer(A[,1], D[,1], FUN="-")
  dy<- outer(A[,2], D[,2], FUN="-")
  return(dx^2 + dy^2)
}

# Sample prior ----
sample.prior<- function(prior) {
  name<- prior[[1]]
  switch(EXPR = name, 
         uniform = {param<- runif(1, prior[[2]],prior[[3]])}, 
         normal = {param<- rnorm(1, prior[[2]],prior[[3]])}, 
         lognormal = {param<- rlnorm(1, prior[[2]],prior[[3]])}, 
         gamma = {param<- rgamma(1, prior[[2]],prior[[3]])},
         beta = {param<- rbeta(1, prior[[2]],prior[[3]])})
  if(!name %in% c("uniform","normal","lognormal","gamma","beta"))
    stop("Prior distribution not recognised")
  param
}

# Prior density ----
prior.density<- function(val, prior, logdens=FALSE) {
  name<- prior[[1]]
  switch(EXPR = name, 
         uniform = {param<- dunif(val, prior[[2]],prior[[3]],log=logdens)}, 
         normal = {param<- dnorm(val, prior[[2]],prior[[3]],log=logdens)}, 
         lognormal = {param<- dlnorm(val, prior[[2]],prior[[3]],log=logdens)}, 
         gamma = {param<- dgamma(val, prior[[2]],prior[[3]],log=logdens)},
         beta = {param<- dbeta(val, prior[[2]],prior[[3]],log=logdens)})
  if(!name %in% c("uniform","normal","lognormal","gamma","beta"))
    stop("Prior distribution not recognised")
  param
}



