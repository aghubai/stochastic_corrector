
sample1F <- function (n, prob=rep(1, n)) {
  # size = 1, replace = FALSE
  # if (length(n) > 0) {stop("This algorithm is for length(n) == 1")}
  # will not throw error if sum(prob > 0) < 1
  
  pos <- prob > 0
  if (any(pos))
    which(pos)[.Internal(sample(sum(pos), 1, FALSE, prob[pos]))]
  else .Internal(sample(n, 1, FALSE, NULL))
}

sampleSF <- function (n, size, prob=rep(1, n)) {
  # size >= 1, replace = FALSE
  # if (length(n) > 0) {stop("This algorithm is for length(n) == 1")}
  # will not throw error if sum(prob > 0) < size
  
  pos <- prob > 0
  s <- sum(pos)
  d <- s - size
  if (d > 0)
    which(pos)[.Internal(sample(sum(pos), size, FALSE, prob[pos]))]
  else
    c(which(pos), which(!pos)[.Internal(sample(n-s, -d, FALSE, NULL))])[
      .Internal(sample(size, size, FALSE, NULL))]
}

fitnessC <- function (vs, tau, vmax, eps) {
  # vs: vector of copy numbers
  # tau: number of genes
  # vmax: maximum number of individual molecules
  # eps: exponent of quality
  
  v <- sum(vs)
  if (v > 0) {
    qual <- prod(vs * tau/v)
    quan <- v/vmax
    return(qual^eps * quan)
  } else {return(0)}
}

fitnessQ <- function (vs, tau, vmax, eps) {
  v <- sum(vs)
  if (v > 0) {
    qual <- prod(vs * tau/v)
    return(qual)
  } else {return(0)}
}

fitnessR <- function (vs, tau, vmax, eps) {
  return(prod(vs * tau/vmax))
}

fitness <- fitnessR

fissionS <- function (vs, tau) { # stochastic
  # vs: vector of copy numbers
  # tau: number of genes
  
  vs1 <- rbinom(tau, vs, 0.5)
  return(rbind(vs1, vs - vs1))
}

fissionN <- function (vs, tau, n) { # stochastic, multiple offspring
  # vs: vector of copy numbers
  # tau: number of genes
  # n: number of offspring vesicles
  
  res <- sapply(vs, function (i) tabulate(sample(n, i, TRUE), n))
  return(res)
}

fission <- fissionS

growS <- function (vs, vmax, tau, a, mu, tp) { # stochastic
  # vs: vector of copy numbers
  # vmax: maximum number of individual molecules
  # tau: number of genes
  # a: vector of gene affinities
  # mu: error rate
  # tp: tau+1 (w/ parasites)
  
  if (sum(vs) > 0) {
    nc <- vmax - sum(vs)
    for (i in seq_len(nc)) {
      x <- sample1F(tp, a*vs)
      if (runif(1) < mu) x <- 1
      vs[x] <- vs[x] +1
    }
  }
  
  return(vs)
}

growsl <- function (vs, vmax, tau, a, mu, tp) { # stochastic, slow
  # see 'growS'
  
  if (sum(vs) > 0) {
    x <- sample1F(tp, a*vs)
    if (runif(1) < mu) x <- 1
    vs[x] <- vs[x] +1
  }
  
  return(vs)
}

grow <- growS

regr <- function (y, x=seq_len(length(y))) {
  xsum <- sum(as.numeric(x))
  ysum <- sum(as.numeric(y))
  prodsum <- sum(as.numeric(x*y))
  x2sum <- sum(as.numeric(x^2))
  y2sum <- sum(as.numeric(y^2))
  n <- length(x)
  
  d <- (n*x2sum - xsum^2)
  m <- (n*prodsum - xsum*ysum) / d
  b <- (x2sum*ysum - xsum*prodsum) / d
  return(c(b, m))
}

lineys <- function (x, y, slope, xs) {
  # give line ys corresponding to xs
  #   based on slope and a single (x, y) point
  
  b <- y - x * slope
  ys <- b + xs * slope
  return(ys)
}

# checking time series equilibrium
isequi <- function (series, debug=FALSE, pos="topright") {
  if (!is.na(tail(series, 1))) {
    
    lgser <- log(series)
    l <- length(series)
    
    ph <- l %/% 2
    h2 <- seq(ph, l)
    l2 <- l - ph + 1
    q3 <- mean(c(ph, l))
    
    st <- lgser[1]
    s2 <- lgser[h2]
    m <- mean(s2)
    sl <- (m - st) / (q3 - 1)
    
    rs <- regr(s2, h2)
    ls <- rs[1] + h2 * rs[2]
    le <- rep(m, l2)
    ln <- lineys(q3, m, sl, h2)
    de <- sum((ls - le)^2)
    dn <- sum((ls - ln)^2)
    
    if (debug) {
      plot(lgser, type="l")
      lines(h2, ls, col=3, lw=3)
      lines(h2, le, col=4, lw=3)
      lines(h2, ln, col=2, lw=3)
      legend(pos, , c("regression line",
                      paste("I(0) error:", round(de, 0)),
                      paste("I(1) error:", round(dn, 0))),
             inset=0.05, bty="n", col=c(3, 4, 2), lwd=3)}
    
    return(dn > de *20)
  } else return(NA)
}

sim <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- sampleSF(N*2, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

simsl <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop <- array(rep(vid, each=N), c(N, tp))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  fitn <- apply(vspop[, -1, drop=FALSE], 1,
                function (i) fitness(i, tau, vmax, eps))
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop[, 1] / rowSums(vspop))
  hst[1, ] <- c(N, fitn[1], pr)
  
  while (inf.loss == FALSE && j <= g*N) {
    
    vsums <- rowSums(vspop) # tmp: no need to recalculate every time
    vsm <- max(vsums)
    while (vsm < vmax) { # max(vsums) < vmax
      i <- sample.int(N, 1, , fitn)
      vspop[i, ] <- growsl(vspop[i, , drop=FALSE], vmax, tau, a, mu, tp)
      fitn[i] <- fitness(vspop[i, -1, drop=FALSE], tau, vmax, eps)
      vsums[i] <- vsm <- vsums[i] +1
    }
    
    h <- sample.int(N, 1) # outflux
    vspop[c(i, h), ] <- fission(vspop[i, , drop=FALSE], tp)
    fitn[i] <- fitness(vspop[i, -1, drop=FALSE], tau, vmax, eps)
    fitn[h] <- fitness(vspop[h, -1, drop=FALSE], tau, vmax, eps)
    
    sf <- sum(fitn > 0)
    if (sf > 0) {
      if (j %% N == 0) {
        pr <- mean(vspop[, 1] / rowSums(vspop), na.rm=TRUE)
        hst[j %/% N +1, ] <- c(sf, mean(fitn), pr)}
      j <- j +1
      
      if (j > g*N & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst))
}

# no selection
sim_nsel <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    surv <- sampleSF(N*2, N)
    sf <- sum(fitn[surv] > 0)
    if (sf > 0) {
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

# single offspring
sim_soff <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- seq(0, N-1)*2 + apply(matrix(fitn, 2), 2, which.max)
      #surv <- sampleSF(N*2, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- sum(fitn[surv] > 0)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

# no stochasticity in fission
sim_nsto <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is deterministic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (inf.loss == FALSE && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp) # also see growD
      vspop2[i*2 - 1:0, ] <- fissionD(vspop1[i, , drop=FALSE], tp)
    }
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- sampleSF(2*N, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

# adjustable number of offsprings (at fission)
sim_noff <- function (g, N, tau, v, mu, o=2, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*o
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*o, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*o - seq_len(o) +1, ] <- fissionN(vspop1[i, , drop=FALSE], tp, o)
    }
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- sampleSF(N*o, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

# deciding whether we need more repetitions
repn <- function (a, b, err=0.05){
  sm <- a + b
  #if (sm == 0) {return(TRUE)}
  mn <- min(a/sm, b/sm)
  return(mn >= err)
}

# do a (max) rep number of simulations
repsim <- function (reps, rerr, rmin, ...) {
  a <- b <- i <- 0
  cont <- TRUE
  ares <- array(, c(reps, 3))  # array of results
  lhst <- vector("list", 0) # list of histories
  
  while (i < reps & cont) {
    i <- i +1
    simout <- simq(...)
    ares[i, ] <- simout[[1]]
    lhst[[i]] <- simout[[2]]
    if (ares[i, 1] > 0)
      a <- a +1 else b <- b +1
    if (i %% rmin == 0)
      cont <- repn(a, b, rerr)}
  
  cont <- repn(a, b, rerr)
  aratio <- a/(a+b)
  mres <- colMeans(ares, na.rm=TRUE)
  
  out <- list(c(a > b, aratio, mres, i), lhst)
  if (cont) out[[1]][1] <- 2 #  # result is uncertain
  return(out)
}

# calculate master node ID and x for whole range
calcndcoos <- function (xrange, scalemax, wlog=10, int=FALSE) {
  m <- 2^scalemax +1
  ilog <- !is.na(wlog)
  if (ilog) xrange <- log(xrange, wlog)
  # 3× faster than: xs <- seq(xrange[1], xrange[2], , m)
  xs <- (xrange[2] - xrange[1]) * (0:(m-1))/(m-1) + xrange[1]
  if (ilog) xs <- wlog^xs
  
  if (int) xs <- round(xs)
  sclvl <- c(0, 0)
  for (i in seq_len(scalemax)) sclvl <- c(rbind(i, sclvl))[-1]
  
  xdupl <- array(, m)
  for (i in seq_len(m)) {
    xwhd <- which(xs == xs[i])
    xdupl[i] <- xwhd[which.min(sclvl[xwhd])]}
  
  return(cbind(xdupl, xs))
}

# find master node
whmaster <- function (ndIDs, m=length(xds[, 1])) {
  ndx <- (ndIDs-1) %%  m +1
  ndy <- (ndIDs-1) %/% m +1
  xm <- xds[ndx, 1] # master node's index
  ym <- yds[ndy, 1]
  mID <- (ym -1) * m + xm
  return(cbind(mID, xds[ndx, 2], yds[ndy, 2]))
}

# convert node ID into given scale
ndT <- function (ndID, scalefr, scaleto) {
  m <- 2^scalefr +1
  ndx <- (ndID-1) %% m
  ndy <- (ndID-1) %/% m
  mT <- 2^scaleto +1
  scd <- 2^(scaleto - scalefr)
  return(ndy*scd*mT + ndx*scd +1)
}

# check equality of an array of numbers
same <- function (i) {
  return(all(i == i[1]))}

# calculate sqL, sqL.nodes, sq.nodesL & whcheck
#   w/ converting node IDs into given scale
sqdivT <- function (sqID, scalefr, scaleto) {
  m <- 2^scalefr
  sqx <- (sqID-1) %% m
  sqy <- (sqID-1) %/% m
  mL <- 2^(scalefr +1)
  mT <- 2^scaleto +1
  scd <- 2^(scaleto - (scalefr +1))
  
  # old squares -> 4 squares (scalefr +1)
  bsq <- sqy*2*mL + sqx*2 +1
  p <- c(0, 1, mL, mL+1)
  sq <- rep(bsq, each=4) + p
  
  # old squares -> 9 nodes (scaleto)
  bnd <- sqy*2*scd*mT + sqx*2*scd
  p <- rep(scd*mT * 0:2, each=3) + scd * 0:2 +1
  nd <- rep(bnd, each=9) + p
  mnd <- matrix(nd, 9)
  
  # new squares' nodes (4 × 4 / old square)
  nsn <- mnd[c(1, 2, 4, 5, 2, 3, 5, 6,
               4, 5, 7, 8, 5, 6, 8, 9), ]
  mnsn <- matrix(nsn, 4)
  
  # unique internal nodes (scaleto)
  ndi <- unique(c(mnd[c(2, 4, 5, 6, 8), ]))
  
  return(list(sq, mnd, mnsn, ndi))
}

# scan a finite parameter space
dynscan <- function (rangex, rangey, scalemax=7, reps=100, rerr=0.2,
                     int=0, tlog=3, lazy=TRUE, debug=FALSE,
                     lbx="x", lby="y", lbm="", draw=FALSE, ...) {
  range1 <<- rangex
  range2 <<- rangey
  
  xds <<- calcndcoos(range1, scalemax, int=int %%  2)
  yds <<- calcndcoos(range2, scalemax, int=int %/% 2)
  
  rmin <- 1/rerr
  
  scale <- 0
  ndc <- 9
  sqc <- 7
  nodes <- matrix(, 0, ndc, dimnames=list(
    NULL, c("nodeID", "xpos", "ypos", "surv", "srel", "snum.m", "mfit.m", "equi.m", "reps")))
  squares <- matrix(, 0, sqc, dimnames=list(
    NULL, c("xmin", "ymin", "xmax", "ymax", "surv", "scale", "sqID")))
  ndhash <- hash()
  histories <- vector("list", 0)
  
  nd0 <- ndT(1:4, scale, scalemax)
  out <- list(1, NA, matrix(nd0, , 1), nd0)
  
  if (draw & debug)
    cmplot(range1, range2, nodes, squares, , tlog,
           main=lbm, xlab=lbx, ylab=lby)
  
  stime <- Sys.time()
  while (TRUE) {
    
    # lazy: ignoring midpoints of homogeneous edges
    if (lazy & scale > 0) {
      
      # finding master for every 9 nodes
      nd9 <- matrix(whmaster(out[[2]])[, 1], 9)
      
      # finding homogeneous edges
      nd4ov <- matrix(nodes[sapply(c(nd9[c(1, 3, 7, 9), ]), function (i)
        ndhash[[as.character(i)]]), 4], 4)
      ishmg <- nd4ov[c(1, 1, 2, 3), ] == nd4ov[c(2, 3, 4, 4), ]
      st <- nd4ov[c(1, 1, 2, 3), ]
      
      # checking whether midpoints are duplicates
      nd4e <- nd9[c(2, 4, 6, 8), ]
      isnew <- sapply(nd4e, function (i)
        is.null(ndhash[[as.character(i)]]))
      newhmg <- nd4e[isnew & ishmg]
      newhet <- nd4e[isnew & !ishmg]
      
      # storing midpoint values on homogeneous edges
      if (length(newhmg) > 0) {
        nhst <- unique(cbind(whmaster(newhmg), st[isnew & ishmg]))
        ndvalx <- cbind(nhst, NA, NA, NA, NA, 0)
        length(histories) <- length(histories) + nrow(ndvalx)
        ndhash[ndvalx[, 1]] <- nrow(nodes) + seq_len(nrow(ndvalx))
        nodes <- rbind(nodes, ndvalx)}
      
      # remaining nodes: 'newhet' and central nodes
      nd1 <- nd9[5, ]
      isnew1 <- sapply(nd1, function (i)
        is.null(ndhash[[as.character(i)]]))
      ndnew <- unique(whmaster(c(newhet, nd1[isnew1])))
      
    } else {
      
      # new nodes: master is missing from hash
      ndIDs <- out[[4]]
      whmst <- unique(whmaster(ndIDs))
      isnew <- sapply(whmst[, 1], function (i)
        is.null(ndhash[[as.character(i)]]))
      ndnew <- whmst[isnew, ] # nodeID, x, y
    }
    
    # calculate and store value of new nodes
    if (nrow(ndnew) > 0) {
      ndval <- array(, c(0, ndc))
      for (i in seq_len(nrow(ndnew))) {
        v <- repsim(reps, rerr, rmin, tau=ndnew[i, 2], v=ndnew[i, 3], ...)
        ndval <- rbind(ndval, c(ndnew[i, ], v[[1]]))
        histories[[length(histories) +1]] <- v[[2]]
        
        if (draw & debug)
          points(ndnew[i, 2], ndnew[i, 3], col=v[[1]][1] +2, pch=16)}
      ndhash[ndval[, 1]] <- nrow(nodes) + seq_len(nrow(ndval))
      nodes <- rbind(nodes, ndval)}
    
    # new squares: node masters do not overlap
    sqsnd <- out[[3]]
    sqmst <- apply(sqsnd, 2, function (i) whmaster(i)[, 1])
    whnew <- which(apply(sqmst, 2, function (i)
      length(unique(i)) == 4))
    
    # checking homogeneity of new squares
    sqind <- apply(sqmst[, whnew, drop=FALSE], 2, function (i)
      sapply(i, function (j) ndhash[[as.character(j)]]))
    if (length(sqind) == 0) break()
    homg <- apply(sqind, 2, function (i) same(nodes[i, 4]))
    sqIDs <- out[[1]][whnew[!homg]]
    
    # storing coordinates of homogeneous squares
    if (sum(homg)) {
      hsq <- sqind[, homg, drop=FALSE]
      newsq <- t(apply(hsq, 2, function (i)
        c(nodes[i, 2:4])[c(1, 5, 4, 8, 9)]))
      squares <- rbind(squares, cbind(newsq, scale, whnew[homg]))}
    
    if (debug) {
      outnd <<- nodes
      outsq <<- squares
      ouths <<- histories
      print(Sys.time() - stime)
      if (draw) cmplot(range1, range2, nodes, squares,
                       tlog=tlog, main=lbm, xlab=lbx, ylab=lby)}
    
    # if possible, divide squares
    if (scale < scalemax & length(sqIDs) > 0) {
      out <- sqdivT(sqIDs, scale, scalemax)
      scale <- scale +1
    } else {break()}
  }
  print(Sys.time() - stime)
  
  if (draw & !debug) cmplot(range1, range2, nodes, squares,
                            tlog=tlog, main=lbm, xlab=lbx, ylab=lby)
  
  return(list(nodes, squares, histories))
}

stascan <- function (rangex, rangey, scalemax=7, reps=100, rerr=0.2,
                     int=0, tlog=3, lazy=NA, debug=FALSE,
                     lbx="x", lby="y", lbm="", draw=FALSE, ...) {
  
  range1 <<- rangex
  range2 <<- rangey
  
  xds <<- calcndcoos(range1, scalemax, int=int %%  2)
  yds <<- calcndcoos(range2, scalemax, int=int %/% 2)
  
  rmin <- 1/rerr
  
  ndc <- 9
  nodes <- matrix(, 0, ndc, dimnames=list(
    NULL, c("nodeID", "xpos", "ypos", "surv", "srel", "snum.m", "mfit.m", "equi.m", "reps")))
  ndhash <- hash()
  histories <- vector("list", 0)
  
  if (draw & debug)
    ndplot(range1, range2, nodes, , tlog,
           main=lbm, xlab=lbx, ylab=lby)
  
  stime <- Sys.time()
  for (i in seq_len(nrow(xds))) {
    for (j in seq_len(nrow(yds))) {
      id <- (yds[j, 1] -1) * nrow(xds) + xds[i, 1]
      if (is.null(ndhash[[as.character(id)]])) {
        p1 <- xds[i, 2]
        p2 <- yds[j, 2]
        v <- repsim(reps, rerr, rmin, tau=p1, v=p2, ...)
        ndhash[id] <- 1
        nodes <- rbind(nodes, c(id, p1, p2, v[[1]]))
        histories[[length(histories) +1]] <- v[[2]]
        
        if (draw & debug)
          points(p1, p2, col=v[[1]][1] +2, pch=16)
      }}}
  
  if (debug) {
    outnd <<- nodes
    ouths <<- histories}
  
  squares <- nd2sq(range1, range2, scalemax, nodes, int)
  if (debug) outsq <<- squares
  print(Sys.time() - stime)
  
  if (draw & !debug) ndplot(range1, range2, nodes,
                            tlog=tlog, main=lbm, xlab=lbx, ylab=lby)
  
  return(list(nodes, squares, histories))
}

simwrapper <- function (inp) {
  
  type <- strsplit(inp, "")[[1]]
  if (! type[1] %in% c("B", "S", "X", "Y")) stop("invalid 1st character")
  if (! type[2] %in% c("C", "Q", "R"))      stop("invalid 2nd character")
  if (! type[3] %in% c("N", "P", "1", "2", "3", "4", "5"))  stop("invalid 3rd character")
  if (! type[4] %in% c("B", "C", "D", "E", "F", "G"))  stop("invalid 4th character")
  if (! type[5] %in% c(NA, "L"))            stop("invalid 5th character")
  if (! type[6] %in% c(NA, 0:9))            stop("invalid 6th character")
  
  if (type[1] == "B") simq <<- sim else {         # basic
    if (type[1] == "S") simq <<- simsl else       # slow: molecules replicate individually
      if (type[1] == "X") simq <<- sim_nsel else  # no selection: survival is fitness-independent
        if (type[1] == "Y") simq <<- sim_nsto}    # no stochasticity: fission is symmetric
  
  if (type[2] == "C") fitness <<- fitnessC else { # complex: both quality and quantity
    if (type[2] == "Q") fitness <<- fitnessQ else # qualitative: vesicle size does not matter
      if (type[2] == "R") fitness <<- fitnessR}   # rate: simplified formula (Q & Q)
  
  if (type[3] == "N") mu <- 0 else {              # no parasites
    if (type[3] == "P") mu <- 0.05 else
      if (type[3] == "1") mu <- 0.01 else
        if (type[3] == "2") mu <- 0.1 else
          if (type[3] == "3") mu <- 0.001 else
            if (type[3] == "4") mu <- 0.0001 else
              if (type[3] == "5") mu <- 0.5}
  
  apar <- 1.1
  if (type[4] == "B") apar <- 2/3 else
    if (type[4] == "C") apar <- 0.95 else
      if (type[4] == "D") apar <- 1 else
        if (type[4] == "E") apar <- 1.1 else
          if (type[4] == "F") apar <- 1.5 else
            if (type[4] == "G") apar <- 2
  
  l <- FALSE
  if (length(type) > 4)
    if (type[5] == "L") l <- TRUE                 # do not scan homogeneous edges
  
  o <- 2
  if (!is.na(o <- as.integer(type[6]) +2)) simq <<- sim_noff
  
  cat(paste0("Simulation ", inp, " started at ", strftime(Sys.time()), "\n"))
  
  g <- 100 # number of generations
  N <- 100 # number of cells
  range1 <- c(1, 100) # number of types (not counting the parasite)
  range2 <- c(1, 100) # number of copies (at start)
  
  out <- dynscan(range1, range2, lbx="genes", lby="copies", lbm=inp,
                 int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o, apar=apar)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0("ndsq_second_", inp, ".out")
  save(range1, range2, nodes, squares, file=fname)
  
  fnamel <- paste0("ndsq_second_", inp, "_large.out")
  save(range1, range2, nodes, squares, histories, file=fnamel)
}

# interpolate values along a line
ndc_inter <- function (vs) {
  vs2 <- c(-1, vs, -1)
  l <- length(vs2)
  m <- which(is.na(vs2))
  if (length(m) > 0) {
    vb <- sapply(m, function (i) vs2[1:i][tail(which(!is.na(vs2[1:i])), 1)])
    va <- sapply(m, function (i) vs2[i:l][head(which(!is.na(vs2[i:l])), 1)])
    e <- va == vb
    vs2[m[e]] <- va[e]
    vs2[m[!e]] <- 2
    u <- va == -1 | vb == -1
    vs2[m[u]] <- NA
    vs <- vs2[-c(1, l)]
  }
  return(vs)
}

# interpolate values of nodes table
ndclean <- function (nodes, scalefr=NA, scaleto=NA, int=NA) {
  
  if (!is.na(scaleto) && !is.na(int)) {
    xds <- calcndcoos(range1, scaleto, int=int %%  2)
    yds <- calcndcoos(range2, scaleto, int=int %/% 2)
    
  } else if (exists("xds") && exists("yds")
             && length(xds) == length(yds))
    scaleto <- log(length(xds[, 1]) -1, 2) else
      stop("Either (scaleto & int) or matching global (xds & yds) is needed.")
  
  # step 1: creating a grid
  xs <- sort(unique(nodes[, 2]))
  ys <- sort(unique(nodes[, 3]))
  nx <- length(xs)
  ny <- length(ys)
  
  xhash <- hash()
  xhash[xs] <- seq_len(nx)
  posx <- sapply(nodes[, 2], function (i) xhash[[as.character(i)]])
  
  yhash <- hash()
  yhash[ys] <- seq_len(ny)
  posy <- sapply(nodes[, 3], function (i) yhash[[as.character(i)]])
  
  ndtable <- matrix(, nx, ny)
  for (i in seq_len(nrow(nodes))) ndtable[posx[i], posy[i]] <- nodes[i, 4]
  
  # step 2: interpolating missing values
  #   points on the boundary
  ndtable[ 1, ] <- ndc_inter(ndtable[ 1, ])
  ndtable[nx, ] <- ndc_inter(ndtable[nx, ])
  ndtable[,  1] <- ndc_inter(ndtable[,  1])
  ndtable[, ny] <- ndc_inter(ndtable[, ny])
  
  #   remaining points
  ndv <- t(apply(ndtable, 1, ndc_inter))
  ndh <- apply(ndtable, 2, ndc_inter)
  e <- ndv == ndh
  ndtable[e] <- ndv[e]
  ndtable[!e] <- 2
  
  # step 3: expanding the grid
  xs2 <- sort(unique(c(nodes[, 2], xds[, 2])))
  ys2 <- sort(unique(c(nodes[, 3], yds[, 2])))
  nx2 <- length(xs2)
  ny2 <- length(ys2)
  
  ndtable2 <- matrix(, nx2, ny2)
  ndtable2[match(xs, xs2), match(ys, ys2)] <- ndtable
  
  # step 4: further interpolation
  #   points with known neighbours
  ndv2 <- t(apply(ndtable2, 1, ndc_inter))
  ndh2 <- apply(ndtable2, 2, ndc_inter)
  v <- !is.na(ndv2) & is.na(ndh2)
  h <- is.na(ndv2) & !is.na(ndh2)
  ndtable2[v] <- ndv2[v]
  ndtable2[h] <- ndh2[h]
  
  #   remaining points
  ndv2b <- t(apply(ndtable2, 1, ndc_inter))
  ndh2b <- apply(ndtable2, 2, ndc_inter)
  e <- ndv2b == ndh2b
  ndtable2[e] <- ndv2b[e]
  ndtable2[!e] <- 2
  
  # step 5: keeping only the new grid
  xs3 <- sort(unique(xds[, 2]))
  ys3 <- sort(unique(yds[, 2]))
  nx3 <- length(xs3)
  ny3 <- length(ys3)
  
  ndtable3 <- matrix(, nx3, ny3)
  ndtable3 <- ndtable2[match(xs3, xs2), match(ys3, ys2), drop=FALSE]
  
  # step 6: creating new nodes array
  nx4 <- length(xds[, 1])
  ny4 <- length(yds[, 1])
  ux <- unique(xds[, 1])
  uy <- unique(yds[, 1])
  
  nodescl <- matrix(, nx3*ny3, 9)
  for (i in seq_len(nx3))
    for (j in seq_len(ny3)) {
      r <- (j-1)*nx3 + i
      id <- (uy[j]-1)*ny4 + ux[i]
      nodescl[r, ] <- c(id, xs3[i], ys3[j], ndtable3[i, j], NA, NA, NA, NA, 0)
    }
  
  # step 7: copying data from original nodes array
  if (is.na(scalefr)) scalefr <- log(sqrt(max(nodes[, 1])) -1, 2)
  scmax <- max(scalefr, scaleto)
  m <- match(ndT(nodes[, 1], scalefr, scmax),
             ndT(nodescl[, 1], scaleto, scmax))
  nodescl[m[!is.na(m)], 5:9] <- nodes[!is.na(m), 5:9]
  
  return(nodescl)
}

# create squares from nodes
nd2sq <- function (range1, range2, scalemax, nodes, int=0,
                   greedy=FALSE, draw=FALSE, clean=TRUE, tlog=3) {
  
  xds <<- calcndcoos(range1, scalemax, int=int %%  2)
  yds <<- calcndcoos(range2, scalemax, int=int %/% 2)
  
  scale <- 0
  sqc <- 7
  squares <- matrix(, 0, sqc, dimnames=list(
    NULL, c("xmin", "ymin", "xmax", "ymax", "surv", "scale", "sqID")))
  
  if (clean) nodes <- ndclean(nodes, scalemax)
  
  ndhash <- hash()
  for (i in seq_len(nrow(nodes)))
    ndhash[nodes[i, 1]] <- i
  
  nd0 <- ndT(1:4, scale, scalemax)
  out <- list(1, NA, matrix(nd0, , 1))
  
  while (TRUE) {
    
    # new squares: node masters do not overlap
    sqsnd <- out[[3]]
    sqmst <- apply(sqsnd, 2, function (i) whmaster(i)[, 1])
    whnew <- which(apply(sqmst, 2, function (i)
      length(unique(i)) == 4))
    
    # checking homogeneity of new squares
    sqind <- apply(sqmst[, whnew, drop=FALSE], 2, function (i)
      sapply(i, function (j)
        if (is.null(ndhash[[as.character(j)]]))
          NA else ndhash[[as.character(j)]]))
    if (length(sqind) == 0) break()
    homg <- apply(sqind, 2, function (i) same(nodes[i, 4]))
    nhomg <- !homg
    
    hna <- is.na(homg)
    homg[hna] <- nhomg[hna] <- FALSE
    
    if (!greedy & any(homg)) {
      q <- which(homg)
      qnd <- sqmst[, whnew[homg], drop=FALSE]
      for (j in seq_len(ncol(qnd))) {
        qcoo <- t(sapply(qnd[c(1, 4), j], whmaster))
        qin <- qcoo[1, 2] < nodes[, 2] & nodes[, 2] < qcoo[2, 2] &
          qcoo[1, 3] < nodes[, 3] & nodes[, 3] < qcoo[2, 3]
        homg[q[j]] <- same(nodes[qin, 4])}
    }
    
    # storing coordinates of homogeneous squares
    if (sum(homg)) {
      hsq <- sqind[, homg, drop=FALSE]
      newsq <- t(apply(hsq, 2, function (i)
        c(nodes[i, 2:4])[c(1, 5, 4, 8, 9)]))
      squares <- rbind(squares, cbind(newsq, scale, whnew[homg]))
    }
    
    sqIDs <- out[[1]][whnew[nhomg]]
    
    # if possible, divide squares
    if (scale < scalemax & length(sqIDs) > 0) {
      out <- sqdivT(sqIDs, scale, scalemax)
      scale <- scale +1
    } else {break()}
  }
  
  if (draw) cmplot(range1, range2, nodes, squares, tlog=tlog)
  
  return(squares)
}

cols.cr <- function (vals, pnts=101) {
  rng <- range(vals, na.rm=TRUE)
  if (diff(rng) > 0)
    v <- round((vals - rng[1]) / diff(rng) * (pnts -1) + 1, 0) else
      v <- vals / rng[1] * (pnts / 2)
    plt <- topo.colors(pnts)
    return(plt[v])
}

# revisit the node's state based on the fitness histories of repeated runs
hsrevisit <- survrev <- function (ndhist, draw=FALSE, inds=NULL, rerr=0.2,
                                  longest=FALSE, limitfail=TRUE) {
  
  reps <- length(ndhist)
  rlen <- sapply(ndhist, nrow)
  
  hsf <- array(, c(reps, max(rlen)))
  rfinx <- rfiny <- requi <- array(, reps)
  for (j in seq_len(reps)) {
    len <- rlen[j]
    fit <- hsf[j, seq_len(len)] <- ndhist[[j]][, 2]
    fin <- rfinx[j] <- sum(!is.na(fit))                      # final surviving generation
    rfiny[j] <- fit[fin]
    requi[j] <- ifelse(fin < len, 0, (2:1)[isequi(fit) +1])  # equilibrium state
  }
  
  requiR <- requi
  uls <- sort(unique(rfinx[requi == 1]))
  for (j in seq_len(length(uls))) {
    
    # lowest fitness of each surviving run
    jsurv <- which(requi == 1 & rfinx == uls[j])
    smins <- apply(hsf[jsurv, , drop=FALSE], 1, min, na.rm=TRUE)
    
    # highest fitness of outlasting parts of failing runs
    jfail <- which(requi == 0 & rfinx > uls[j])
    if (length(jfail) == 0) fmax <- 0 else
      fmax <- max(apply(hsf[jfail, -seq_len(uls[j]), drop=FALSE], 1, max, na.rm=TRUE))
    
    # new uncertainty: more fit runs have failed, when inspected for longer
    uncR <- jsurv[smins < fmax]
    requiR[uncR] <- 3
  }
  
  # revision of survival state
  freq <- freqorig <- tabulate(requiR +1, 4)
  
  if (limitfail) {
    # a failed run can only cancel out a single equilibrium run
    srem <- max(freq[4] - freq[1], 0)
    freq[2] <- freq[2] + srem
    freq[4] <- freq[4] - srem
  }
  
  ncert <- sum(freq[1:2])
  srel <- freq[2] / ncert
  
  stR <- 2
  if (ncert > 0)
    if (srel < rerr) stR <- 0 else
      if (srel > 1 - rerr) stR <- 1
  
  if (longest) {
    # longest runs survive (w/ non-equilibrium cases excluded)
    cert <- which(requiR < 2)
    lcert <- rfinx[cert]
    stcertlong <- requiR[cert[lcert == max(lcert)]]
    
    leff <- TRUE
    if (all(stcertlong == 1) & stR == 0) stR <- 2 else
      leff <- FALSE
  }
  
  if (draw) {
    hsplot(list(ndhist), 1, log="y", showtitle=FALSE)
    if (length(inds) == 2)
      title(main=paste0("Repeated runs at tau=", inds[1], ", c=", inds[2]))
    
    cols <- c(2, 3, 8, 4)
    points(rfinx -1, rfiny, pch=16, col=cols[requiR +1])
    
    # legend of runs' state
    rtext <- c("failure", "success", "non-equilibrium", "uncertain")
    leg <- paste(rtext, freqorig, sep=": ")
    o <- c(2, 1, 4, 3)
    legend("topright", , leg[o], bty="n", inset=0.05, col=cols[o], pch=16)
    
    # legend of node's state
    ntext <- c("gets lost", "survives", "may survive")
    res <- paste("information", ntext[stR+1])
    if (longest && leff) res <- c(res, "(based on longest run)")
    if (limitfail && srem > 0) res <-
      c(res, paste("(less uncertainty: ", paste(freq[o], collapse=" "), ")", sep=""))
    legend("bottomright", , res, bty="n", inset=0.05)
  }
  
  return(stR)
}

# draw nodes (options: integer values, log scale)
ndplot <- function (range1, range2, nodes, squares=NULL, int=0, tlog=3, border=NULL,
                    lazy=TRUE, mirror=FALSE, repcol=9, ...) {
  tlg <- c("", "x", "y", "xy")[tlog %% 4 +1]
  coo <- list(integer(0), 2, 3, 2:3)[[int %% 4 + 1]]
  nodes[, coo] <- round(nodes[, coo])
  
  if (lazy) nodes <- nodes[which(nodes[, repcol] > 0), , drop=FALSE]
  
  if (mirror) {
    rngs <- rbind(range1, range2)
    range1 <- rngs[2, ]
    range2 <- rngs[1, ]
    nodes[, 2:3] <- nodes[, 3:2]
    squares[, 1:4] <- squares[, c(2, 1, 4, 3)]}
  
  plot(range1, range2, type="n", log=tlg, ...)
  cols <- nodes[, 4] + 2
  cols[which(cols %% 1 != 0)] <- 4
  points(nodes[, 2], nodes[, 3], col=cols, pch=16)
}

# draw squares (options: integer values, log scale, without borders)
sqplot <- function (range1, range2, squares, int=0, tlog=3, border=TRUE, ...) {
  tlg <- c("", "x", "y", "xy")[tlog %% 4 +1]
  coo <- list(integer(0), c(1, 3), c(2, 4), 1:4)[[int %% 4 + 1]]
  squares[, coo] <- round(squares[, coo])
  
  plot(range1, range2, type="n", log=tlg, ...)
  ifelse(border, brd <- par("fg"), brd <- NA)
  rect(squares[, 1], squares[, 2], squares[, 3], squares[, 4],
       col=squares[, 5] + 2, border=brd)
}

# draw nodes with gradient colouring
ndplot_grad <- function (range1, range2, nodes, int=0, tlog=3, pvar=4, cols=NULL, ...) {
  tlg <- c("", "x", "y", "xy")[tlog %% 4 +1]
  coo <- list(integer(0), 2, 3, 2:3)[[int %% 4 + 1]]
  nodes[, coo] <- round(nodes[, coo])
  
  if (is.null(cols)) {
    cols <- cols.cr(nodes[, pvar])
    ttl <- paste0("Colouring based on \'", colnames(nodes)[pvar], "\'")
  } else ttl <- "Custom colouring"
  
  plot(range1, range2, type="n", log=tlg, main=ttl, ...)
  points(nodes[, 2], nodes[, 3], col=cols, pch=16)
}

# draw both nodes and squares (options: integer values, log scale, without borders)
cmplot <- function (range1, range2, nodes, squares, int=0, tlog=3, border=FALSE,
                    lazy=TRUE, mirror=FALSE, repcol=9, ...) {
  tlg <- c("", "x", "y", "xy")[tlog %% 4 +1]
  coo <- list(integer(0), c(1, 3), c(2, 4), 1:4)[[int %% 4 + 1]]
  nodes[, coo] <- round(nodes[, coo])
  squares[, coo] <- round(squares[, coo])
  if (lazy) nodes <- nodes[which(nodes[, repcol] > 0), , drop=FALSE]
  
  if (mirror) {
    rngs <- rbind(range1, range2)
    range1 <- rngs[2, ]
    range2 <- rngs[1, ]
    nodes[, 2:3] <- nodes[, 3:2]
    squares[, 1:4] <- squares[, c(2, 1, 4, 3)]}
  
  sqnodes <- unique(rbind(squares[, c(1, 2)], squares[, c(1, 4)],
                          squares[, c(3, 2)], squares[, c(3, 4)]))
  ndpos_sq <- paste(sqnodes[, 1], sqnodes[, 2])
  ndpos_nd <- paste(nodes[, 2], nodes[, 3])
  redund <- ndpos_nd %in% ndpos_sq
  snodes <- nodes[!redund, , drop=FALSE]
  
  uncert <- nodes[, repcol] > min(nodes[, repcol], Inf, na.rm=TRUE)
  unodes <- nodes[uncert, , drop=FALSE]
  
  rnodes <- nodes[redund & !uncert, , drop=FALSE]
  
  plot(range1, range2, type="n", log=tlg, ...)
  ifelse(border, brd <- par("fg"), brd <- NA)
  rect(squares[, 1], squares[, 2], squares[, 3], squares[, 4],
       col=squares[, 5] + 2, border=brd)
  
  cols <- snodes[, 4] + 2
  points(snodes[, 2], snodes[, 3], col=cols, pch=16)  # separate nodes
  points(unodes[, 2], unodes[, 3])                    # uncertain nodes
  points(rnodes[, 2], rnodes[, 3], pch=".", cex=2)    # rest
}

# draw histories (of a single node)
hsplot <- function (histories, nodenum=NULL, repnum=NULL, series=1, add=FALSE,
                    nodes=NULL, x=NULL, y=NULL, showequi=FALSE, shownums=FALSE,
                    showtitle=TRUE, ...) {
  # if repnum is 'NULL', shows history of all runs
  # series 0, 1 & 2 are respectively surv, fitn & par
  
  n <- nodenum
  nd <- nodes
  if (all(c(is.null(n), !is.null(nd), !is.null(x), !is.null(y))))
    n <- which(nd[, 2] == x & nd[, 3] == y)
  
  if (length(n) == 1 && length(histories[[n]]) > 0) {
    h <- histories[[n]]
    nrep <- length(h)
    if (is.null(repnum)) rs <- seq_len(nrep) else rs <- repnum
    
    if (!is.null(nd)) {
      e <- nd[n, ]
      ttl <- paste0("Node: x=", e[2], ", y=", e[3], ", ID=", e[1],
                    ", row=", n, ", reps=", nrep)
    } else ttl <- paste0("Node: row=", n, ", reps=", nrep)
    
    ty <- series+1
    l <- max(unlist(lapply(h[rs], function (i) length(i[, ty]))))
    lns <- simplify2array(lapply(h[rs], function (i) {
      v <- i[, ty]; length(v) <- l; return(v)}))
    
    xs <- seq(0, l-1)
    yrng <- range(lns, na.rm=TRUE)
    nrs <- length(rs)
    if (nrs > 1) lcl <- rgb(0, 0, 0, 1/log(nrs, 2)) else lcl <- 1
    
    if (!showtitle) ttl <- ""
    if (!add) plot(c(0, l-1), yrng, type="n", main=ttl,
                   xlab="generations", ylab=c("#survivors", "fitness", "parasite ratio")[ty], ...)
    invisible(apply(lns, 2, function (i) lines(xs, i, col=lcl)))
    
    if (showequi)
      invisible(lapply(h[rs], function (i) {
        d <- i[, ty]
        last <- sum(!is.na(d))
        ie <- ifelse(last < length(d), 2, (4:3)[isequi(d) +1])
        points(last-1, d[last], pch=16, col=ie)}))
    
    if (shownums)
      invisible(sapply(seq(nrs), function (j) {
        i <- lns[, j]; last <- sum(!is.na(i)); text(last-1, i[last], j)}))
  }
}

# compare (longest) fitness histories of several nodes
hscompareplot <- function (histories, nodes, pos, horiz=TRUE,
                           dimname=NULL, qmin=NULL, qmax=NULL) {
  
  if (horiz) d <- 3 else d <- 2
  sel <- nodes[, 9] > 0
  if (!is.null(qmin)) sel <- sel & nodes[, 5-d] >= qmin
  if (!is.null(qmax)) sel <- sel & nodes[, 5-d] <= qmax
  
  rnodes <- nodes[sel, ]
  coos <- sort(rnodes[rnodes[, d] == pos, 5-d])
  lc <- length(coos)
  
  if (lc > 0) {
    lhist <- vector("list", lc)
    for (j in 1:lc) {
      h <- histories[[which(nodes[, 5-d] == coos[j] & nodes[, d] == pos)]]
      lh <- unlist(lapply(h, function (k) length(k[, 2])))
      lhist[[j]] <- h[[which.max(lh)]][, 2]
    }
    
    if (!is.null(dimname)) dname <- dimname else {
      if (horiz) dname <- "y" else dname <- "x"}
    
    xrng <- c(1, max(unlist(lapply(lhist, function (k) length(k)))))
    yrng <- range(unlist(lapply(lhist, function (k) range(k, na.rm=TRUE))))
    plot(xrng, yrng, type="n", log="y", ylab="mean fitness of population",
         xlab="generations", main=paste("Comparing nodes,", dname, "=", pos))
    legend("bottomright", , c("outer: node state", "inner: run state"), bty="n", inset=0.05)
    
    for(j in 1:lc) {
      hlc <- lhist[[j]]
      lines(1:length(hlc), hlc)
      e <- sum(!is.na(hlc))
      text(e + 0.025*xrng[2], hlc[e], coos[j])
      eq1 <- max(isequi(hlc[1:e]), 2*(e == 1001))
      eq2 <- nodes[nodes[, 5-d] == coos[j] & nodes[, d] == pos, 4]
      points(e, hlc[e], col=(2:4)[eq2 +1], pch=16, cex=1.5)
      points(e, hlc[e], col=0, pch=16)
      points(e, hlc[e], col=(2:4)[eq1 +1], pch=20)
    }
  }
}

# drawing phylogenetic tree
phylplot <- function (survs, o=2, draw=4, ...) {
  N <- ncol(survs)
  g <- nrow(survs)
  Ns <- seq_len(N)
  gs <- seq_len(g)
  plt <- rainbow(trunc(N*1.1))[Ns]
  
  so <- array(, c(g, N))
  for (i in gs)
    so[i, ] <- (survs[i, ] -1) %/% o +1
  
  sc <- array(, c(g, N))
  sc[1, ] <- so[1, ]
  for (i in seq(2, g))
    sc[i, ] <- sc[i-1, so[i, ]]
  
  for (i in gs)
    sc[i, ] <- sort(sc[i, ], na.last=TRUE)
  
  sp <- array(, c(g, N))
  sp[1, ] <- sort(lnext <- so[1, ])
  for (i in seq(2, g))
    sp[i, ] <- sort(lnext <- order(order(lnext))[so[i, ]], na.last=TRUE)
  
  if (1 %in% draw) {
    # parent-offspring links
    plot(c(1, N), c(0, g), type="n", xlab="population", ylab="time", ylim=c(g, 0), ...)
    for (i in gs)
      for (j in Ns)
        segments(sp[i, j], i-1, j, i, col=plt[sc[i, j]])}
  
  if (2 %in% draw) {
    # horizontal lines
    plot(c(1, N), c(0, g), type="n", xlab="population", ylab="time", ylim=c(g, 0), ...)
    for (i in gs)
      for (j in unique(sc[i, ]))
        lines(y <- which(sc[i, ] == j), rep(i, length(y)), col=plt[j])}
  
  if (3 %in% draw) {
    # vertical boundaries
    plot(c(1, N), c(0, g), type="n", xlab="population", ylab="time", ylim=c(g, 0), ...)
    for (i in Ns)
      lines(apply(sc, 1, function (j) max(which(j == i), -Inf) +0.5), gs)}
  
  if (4 %in% draw) {
    # polygons
    plot(c(1, N), c(0, g), type="n", xlab="population", ylab="time", ylim=c(g, 0), ...)
    for (i in Ns) {
      xmin <- apply(sc, 1, function (j) min(which(j == i),  Inf))
      xmax <- apply(sc, 1, function (j) max(which(j == i), -Inf))
      bnd <- rbind(c(i, xmin, rev(xmax)), c(0, gs, rev(gs)))
      bnd <- bnd[, abs(bnd[1, ]) != Inf, drop=FALSE]
      polygon(bnd[1, ], bnd[2, ], col=plt[i], border=8)}}
  
  if (5 %in% draw) {
    # polygons with links
    plot(c(1, N), c(0, g), type="n", xlab="population", ylab="time", ylim=c(g, 0), ...)
    for (i in Ns) {
      xmin <- apply(sc, 1, function (j) min(which(j == i),  Inf))
      xmax <- apply(sc, 1, function (j) max(which(j == i), -Inf))
      bnd <- rbind(c(i, xmin, rev(xmax)), c(0, gs, rev(gs)))
      bnd <- bnd[, abs(bnd[1, ]) != Inf, drop=FALSE]
      polygon(bnd[1, ], bnd[2, ], col=plt[i], border=NA)}
    for (i in gs)
      for (j in Ns)
        segments(sp[i, j], i-1, j, i, col=8)}
}

# drawing value on phylogenetic tree
phylvalplot <- function (vals, survs, o=2, valtransf=FALSE, logcolours=TRUE,
                         showtree=FALSE, showbest=FALSE, ...) {
  N <- ncol(survs)
  g <- nrow(survs)
  Ns <- seq_len(N)
  gs <- seq_len(g)
  plt <- rainbow(trunc(N*1.1))[Ns]
  
  if (valtransf) {
    valsAll <- vals
    vals <- array(, c(g, N))
    for (i in gs) vals[i, ] <- valsAll[i, survs[i, ]]
  }
  
  rowMaxs <- apply(vals, 1, max, na.rm=TRUE)
  vrel <- vals / rowMaxs
  
  if (logcolours) {
    vrel <- pmin(15, vrel * 14 +1)
    vrel <- 1 - trunc(log(16 - vrel, 2))/4
    
  } else
    vrel <- pmin(1, trunc(vrel * 4 +1) / 4)
  
  vrel <- vrel^2
  vrel[vals == 0] <- 0
  vrel <- array(vrel, c(g, N))
  
  # ID of survivors' parents in the same generation
  so <- array(, c(g, N))
  for (i in gs)
    so[i, ] <- (survs[i, ] -1) %/% o +1
  
  # ID of survivors' parents in the first generation
  sc <- array(, c(g, N))
  sc[1, ] <- so[1, ]
  for (i in seq(2, g))
    sc[i, ] <- sc[i-1, so[i, ]]
  
  for (i in gs)
    sc[i, ] <- sort(sc[i, ], na.last=TRUE)
  
  sp <- array(, c(g, N))
  sp[1, ] <- sort(lnext <- so[1, ])
  vrel[1, ] <- vrel[1, order(so[1, ])]
  for (i in seq(2, g)) {
    vrel[i, ] <- vrel[i, order(lnext)]
    sp[i, ] <- sort(lnext <- order(order(lnext))[so[i, ]], na.last=TRUE)}
  vrel[is.na(vrel)] <- 0
  
  plot(c(1, N), c(0, g), type="n", xlab="population", ylab="time", ylim=c(g, 0), ...)
  
  if (showtree)
    for (i in gs)
      for (j in Ns)
        segments(sp[i, j], i-1, j, i, col=8)
  
  xs <- sort(unique(c(sc)))
  for (i in xs) {
    clr <- col2rgb(plt[i]) / 255
    pos <- which(sc == i) -1
    points(pos %/% g +1, pos %% g +1, pch=16,
           col=rgb(clr[1], clr[2], clr[3], vrel[pos]))}
  
  if (showbest) {
    b <- which.max(table(sc[g, ]))
    best <- ifelse(length(b) > 0, as.integer(names(b)), NULL)
    title(sub=paste0("(", best, ")"))
    return(best)}
}

# drawing the compositions of all vesicles
popplot <- function (vspop, vmax, cols) {
  
  N <- nrow(vspop)
  tp <- ncol(vspop)
  
  if (missing(cols))
    cols <- c("red", "blue", "grey10", "grey50", "grey30", "grey70")
  
  csum <- t(rbind(0, apply(vspop, 1, cumsum)))
  mcol <- array(3:4, c(N, tp))
  mcol[, -seq(1, N, 2)] <- mcol[, -seq(1, N, 2)] + 2
  mcol[, 1] <- 2
  ext <- array(, c(0, 2))
  
  plot(c(0, N), c(0, vmax), type="n",
       main="Distribution of genes", xlab="population", ylab="molecules")
  for (i in 1:N)
    for (j in 1:tp)
      if (vspop[i, j] > 0)
        rect(i-1, csum[i, j], i, csum[i, j+1], col=cols[mcol[i, j]], border=NA) else
          ext <- rbind(ext, c(i-0.5, csum[i, j+1]))
  
  if (nrow(ext) > 0)
    points(ext[, 1], ext[, 2], pch=16, col=cols[1])
}

# drawing the rank-abundance curves of vesicular compositions
distplot <- function (taus, vs, vspoplist, verbose=TRUE,
                      draw=TRUE, log="xy", wholepop=TRUE, cols) {
  
  n <- length(vspoplist)
  N <- dim(vspoplist[[1]])[2]
  
  vpl_equi <- lapply(vspoplist, function (i) i[dim(i)[1], , ])
  vpl_sorted <- lapply(vpl_equi, function (i) apply(i[, -1], 1, sort))
  vpl_median <- lapply(vpl_sorted, function (i) apply(i, 1, median))
  
  if (draw) {
    if (missing(cols)) cols <- rainbow(n, alpha=0.05)
    cols_o <- substr(cols, 1, 7)
    if (wholepop) lcols <- rep(1, n) else lcols <- cols_o
    leg <- sprintf("%*.0f, %*.0f",
                   max(nchar(taus)), taus, max(nchar(vs)), vs)
    
    plot(range(1, taus), c(1, range(vpl_sorted)[2]), type="n", log=log,
         main="Vesicular distribution of molecules", xlab="gene rank", ylab="abundance")
    if (wholepop)
      for (i in seq_len(n))
        for (j in seq_len(N))
          lines(rev(seq_len(taus[i])), vpl_sorted[[i]][, j], col=cols[i])
    for (i in seq_len(n))
      lines(rev(vpl_median[[i]]), lwd=2, col=lcols[i])
    legend("bottomright", , rev(leg), bty="n", inset=0.05,
           title="genes, copies", col=rev(cols_o), pch=16, text.font=10)
  }
  
  if (verbose)
    return(list(vpl_sorted=vpl_sorted, vpl_median=vpl_median))
}

fitnplot <- function (fitnsAll, survs, o=2, coloured=FALSE,
                      parents=NULL, first=NULL, ...) {
  
  N <- ncol(survs)
  Ns <- seq_len(N)
  N2s <- seq_len(N*2)
  prn <- (N2s-1) %/% o +1
  alpha <- ifelse(N > 1, 1/log(N, 2), 1)
  
  g <- nrow(survs)
  e <- 0
  if (!is.null(first[1])) {
    fitnsAll <- rbind(first[1], fitnsAll)
    survs <- rbind(Ns*2, survs)
    g <- g +1
    e <- 1}
  gs <- seq_len(g)
  
  so <- array(, c(g, N))
  for (i in gs)
    so[i, ] <- (survs[i, ] -1) %/% o +1
  
  vals <- fitnsAll
  yrng <- range(vals[vals!=0], na.rm=TRUE)
  vals[vals==0] <- NA
  
  plot(range(gs)-e, yrng, type="n", xlab="time", ylab="fitness", log="y", ...)
  
  sc <- array(, c(g, N))
  sc[1, ] <- so[1, ]
  for (i in seq(2, g))
    sc[i, ] <- sc[i-1, so[i, ]]
  
  if (is.null(parents[1])) n <- N else n <- length(parents)
  
  if (!coloured) plt <- rep(rgb(0, 0, 0, alpha), n) else
    plt <- rainbow(trunc(n*1.1), alpha=alpha)[seq_len(n)]
  
  if (is.null(parents[1]))
    for (i in seq(2, g))
      for (j in N2s) {
        yp <- vals[i-1, survs[i-1, prn[j]]]
        cl <- plt[sc[i-1, prn[j]]]
        if (!is.na(v <- vals[i, j]))
          segments(i-e, v, i-e-1, yp, col=cl) else points(i-e-1, yp, col=cl)}
  
  if (!is.null(parents[1]))
    for (i in seq(2, g))
      for (j in N2s)
        if (!is.na(p <- match(sc[i-1, prn[j]], parents))) {
          yp <- vals[i-1, survs[i-1, prn[j]]]
          cl <- plt[p]
          if (!is.na(v <- vals[i, j]))
            segments(i-e, v, i-e-1, yp, col=cl) else points(i-e-1, yp, col=cl)}
}

noffplot <- function (vals, survs, o=2, valtransf=FALSE,
                      showgens=TRUE, coloured=FALSE, alpha=NULL, ...) {
  
  N <- ncol(survs)
  Ns <- seq_len(N)
  g <- nrow(survs)
  gs <- seq_len(g)
  
  if (valtransf) {
    valsAll <- vals
    vals <- array(, c(g, N))
    for (i in gs) vals[i, ] <- valsAll[i, survs[i, ]]
  }
  
  so <- noffs <- array(, c(g, N))
  for (i in seq(2, g)) {
    so[i, ] <- (survs[i, ] -1) %/% o +1
    noffs[i, ] <- tabulate(so[i, ], N)
  }
  
  ys <- seq(0, o)
  pts <- array(, c(g, o+1))
  for (i in seq(2, g))
    for (j in ys)
      pts[i, j+1] <- mean(vals[i-1, noffs[i, ] == j], na.rm=TRUE)
  
  plot(range(pts, na.rm=TRUE), c(-0.1, o + 0.1), type="n",
       xlab="Mean parent quality", ylab="Number of offsprings", ...)
  
  if (is.null(alpha)) alpha <- 1/log(g +1, 2)
  
  if (showgens) {
    
    if (!coloured)
      plt <- rep(rgb(1/2, 1/2, 1/2, alpha), g) else
        plt <- rainbow(trunc(g*1.1), alpha=alpha)[gs]
      
      for (i in seq(2, g)) {
        points(pts[i, ], ys, col=plt[i], pch=16)
        
        if (abs(cor(ys, pts[i, ])) > 0.5)
          abline(lm(ys ~ pts[i, ]), col=plt[i])}
  }
  
  mpts <- colMeans(pts[-1, ])
  
  points(mpts, ys, pch=16)
  abline(lm(ys ~ mpts), lwd=2)
}

survplot <- function (fitnsAll, hst, ...) {
  
  N <- ncol(fitnsAll) %/% 2
  nsurv <- apply(fitnsAll, 1, function (i) sum(i > 0))
  s <- rbind(nsurv, 2*N - nsurv)
  g <- ncol(s)
  
  par(mar = c(5, 4.5, 4, 4) + 0.1)
  barplot(s, col=c(0, 8), space=0, border=NA, main="Survivors among offsprings",
          xlab="generations", ylab="number of survivors", ...)
  axis(1, seq(0, g, , 5), pos=0)
  sapply(1:2, function (i) lines(c(0, g), c(i, i)*N, lty=2))
  
  maxf <- max(hst[, 2], na.rm=TRUE)
  mfitn <- hst[, 2] * 2*N / maxf
  lines(seq_len(length(mfitn)) -1.5, mfitn, col=4, lwd=2)
  axis(4, seq(0, 2*N, , 5), format(seq(0, maxf, , 5), digits=2, sci=TRUE))
  mtext("mean fitness", 4, 2.5)
  
  legend("topright", , c("# survivors", "mean fitness"), bty="n", inset=0.05,
         col=c(1, 4), lwd=c(NA, 2), pch=c(22, NA), pt.bg=c("white", NA))
  par(mar = c(5, 4, 4, 2) + 0.1)
}

# cmplot of the 100% certain nodes (to serve as a conservative estimate)
certplot <- function (range1, range2, nodes, squares, int=0, tlog=3, border=FALSE,
                      lazy=TRUE, repcol=9, scalemax=7, greedy=FALSE, clean=TRUE, ...) {
  
  unc <- which(nodes[, 5] %% 1 != 0)
  nodes[unc, 4] <- 2
  nodes[unc, 9] <- min(nodes[-which(nodes[, 9] == 0), 9])
  squares <- nd2sq(range1, range2, scalemax, nodes, int, greedy, FALSE, clean, tlog)
  cmplot(range1, range2, nodes, squares, int, tlog, border, lazy, repcol, ...)
}

# select point on the two sides of the threshold
thrpoints <- function (nodes, valabove=0, valbelow=1, repcol=9) {
  # valabove, valbelow: value below and above the threshold
  
  # searching for the rightmost 1s & the leftmost 0s (in horizontal slices)
  ys <- sort(unique(nodes[, 3]))
  ptsa <- array(, c(0, 4))
  ptsb <- array(, c(0, 4))
  for (ypos in ys) {
    nodes_sub <- nodes[nodes[, 3] == ypos, ]
    if (all(c(valabove, valbelow) %in% unique(nodes_sub[, 4]))) {
      o <- order(nodes_sub[, 2])
      b <- cbind(ypos, nodes_sub[o, c(2, 4, repcol)])
      sts <- nodes_sub[o, 4]
      
      va <- head(which(sts == valabove), 1)
      vb <- tail(which(sts == valbelow), 1)
      if (va > vb) {
        ptsa <- rbind(ptsa, b[va, ])
        ptsb <- rbind(ptsb, b[vb, ])
      }
    }
  }
  
  # adding topmost and bottommost missing slice
  yrng <- range(ptsa[, 1], ptsb[, 1])
  range2 <- range(nodes[, 3])
  
  if (range2[1] < yrng[1]) {
    ybot <- max(-Inf, ys[ys < yrng[1]])
    nodes_sub <- nodes[nodes[, 3] == ybot, ]
    o <- order(nodes_sub[, 2])
    ndbot <- c(ybot, nodes_sub[o[1], c(2, 4, repcol)])
    if (ndbot[3] == valabove) {
      ptsa <- rbind(c(ndbot), ptsa)
      ptsb <- rbind(c(ndbot[1], NA, NA, NA), ptsb)
    }
  }
  
  if (range2[2] > yrng[2]) {
    ytop <- min(Inf, ys[ys > yrng[2]])
    nodes_sub <- nodes[nodes[, 3] == ytop, ]
    o <- order(nodes_sub[, 2], decreasing=TRUE)
    ndtop <- c(ytop, nodes_sub[o[1], c(2, 4, repcol)])
    if (ndtop[3] == valbelow) {
      ptsa <- rbind(ptsa, c(ndtop[1], NA, NA, NA))
      ptsb <- rbind(ptsb, c(ndtop))
    }
  }
  
  # outputting result
  colnames(ptsa) <- colnames(ptsb) <-
    c("ypos", "xpos", "surv", "reps")
  
  out <- list(pbelow=ptsb, pabove=ptsa)
  return(out)
}

# finding spline with smallest error for df (via incremental weight optimization)
optispline <- function (xs, ya, yb, epsylon=0, fixedends=TRUE, draw=FALSE) {
  
  dfmin <- 2
  dfmax <- length(unique(xs))
  
  xall <- rep(xs, 2)
  yall <- c(ya, yb)
  wsta <- rep(1, length(xs))
  wsta[c(1, length(xs))] <- 1000 # against edge artifacts
  
  if (fixedends) {
    l <- length(xs)
    ya[1] <- yb[1] <- mean(c(ya[1], yb[1]))
    ya[l] <- yb[l] <- mean(c(ya[l], yb[l]))
    yall <- c(ya, yb)
  }
  
  i <- dfmin
  err_i <- err_best <- Inf
  wpair <- wsta
  
  while (TRUE) {
    ispline <- smooth.spline(xall, yall, rep(wpair, 2), i)
    yspl <- predict(ispline, xs)$y
    off <- yspl > ya | yspl < yb
    wpair[off] <- wpair[off] *2
    err <- sum(pmax(0, yspl - ya, na.rm=TRUE),
               pmax(0, yb - yspl, na.rm=TRUE))
    
    if (err < err_i) {
      splbest_i <- ispline
      err_i <- err
      
    } else {
      
      if (err_i < err_best) {
        splbest_best <- splbest_i
        err_best <- err_i
        dfout <- i
      }
      
      if (i >= dfmax || err <= epsylon) break()
      
      i <- i +1
      splbest_i <- NA
      err_i <- Inf
      wpair <- wsta
    }
  }
  
  if (draw) {
    yspl <- predict(splbest_best, xs)$y
    plot(range(xs), range(yall), type="n", log=c("xy", "")[logscale +1],
         xlab="x", ylab="y", main=paste0("df = ", dfout, ", epsylon = ", epsylon))
    legend("bottomright", , c("above", "below"), bty="n", inset=0.05, col=2:3, pch=16)
    points(xs, ya, pch=16, col=2)
    points(xs, yb, pch=16, col=3)
    lines(xs, yspl, lwd=2, col=4)
  }
  
  return(list(spl=splbest_best, error=err_best, df=dfout))
}

# calculate threshold in parameter space scan ('inner' spline)
calcthr <- function (filename, epsylon=0.05, x, thralongx=TRUE, logscale=TRUE, repcol=9,
                     spline.hi=NULL, spline.lo=NULL, draw=FALSE, add=FALSE, mirror=TRUE, pos) {
  # 'epyslon': maximum error of spline
  # 'thralongx': boundary point pairs are horizontal (vs. vertical)
  #   cf. y -> x is unambiguous (vs. x -> y)
  # 'logscale': small x & y values have a greater impact on spline fit
  # 'spline.hi' and 'spline.lo' should be a spline created with calcthr
  #   using the same 'thralongx' and 'logscale' settings
  # 'draw', 'mirror' and 'pos' make no difference in the output
  
  load(filename) # used: range1, range2, nodes
  if (!thralongx) nodes[, 2:3] <- nodes[, 3:2]
  valabove <- nodes[which(nodes[, 2]==range1[2] &
                            nodes[, 3]==range2[1]), 4]
  valbelow <- nodes[which(nodes[, 2]==range1[1] &
                            nodes[, 3]==range2[2]), 4]
  
  # finding threshold points
  thr <- thrpoints(nodes, valabove, valbelow, repcol)
  xs <- thr$pabove[, 1]
  ya <- thr$pabove[, 2]
  yb <- thr$pbelow[, 2]
  
  # filling two sets
  r <- min(ya / yb, na.rm=TRUE)
  w <- which(is.na(ya))
  ya[w] <- yb[w] * r
  w <- which(is.na(yb))
  yb[w] <- ya[w] / r
  
  if (logscale) {
    xs <- log(xs)
    ya <- log(ya)
    yb <- log(yb)
  }
  
  # setting bounds
  exc <- c(-1, -length(xs))
  if (!is.null(spline.hi))
    ya[exc] <- pmin(ya[exc], predict(spline.hi, xs[exc])$y)
  if (!is.null(spline.lo))
    yb[exc] <- pmax(yb[exc], predict(spline.lo, xs[exc])$y)
  
  ospl <- optispline(xs, ya, yb, epsylon)
  
  # x: points for spline interpolation
  if (missing(x)) x <- seq(min(xs), max(xs), , 100) else
    if (logscale) x <- log(x)
  y <- predict(ospl$spl, x)$y
  
  if (logscale) {
    x <- exp(x)
    y <- exp(y)
  }
  
  if (thralongx) {
    yn <- x
    x <- y
    y <- yn
  }
  
  # points
  if (!thralongx) {
    ptsa <- list(x=thr$pbelow[, 1], y=thr$pbelow[, 2])
    ptsb <- list(x=thr$pabove[, 1], y=thr$pabove[, 2])
    
  } else{
    ptsa <- list(x=thr$pabove[, 2], y=thr$pabove[, 1])
    ptsb <- list(x=thr$pbelow[, 2], y=thr$pbelow[, 1])
  }
  
  # draw
  if (draw) {
    mn <- paste0("Threshold spline, df=", ospl$df, ", epsylon=", epsylon)
    
    if (!mirror) {
      if (!add) {
        plot(range1, range2, log="xy", type="n",
             main=mn, xlab="genes", ylab="copies")
        points(ptsa$x, ptsa$y, pch=16, col=2)
        points(ptsb$x, ptsb$y, pch=16, col=3)
        
        if (missing(pos)) pos <- "bottomright"
        legend(pos, , c("right", "left"), col=2:3, pch=16, bty="n", inset=0.05)
      }
      lines(x, y, lwd=2, col=4)
      
    } else {
      if (!add) {
        plot(range2, range1, log="xy", type="n",
             main=mn, xlab="copies", ylab="genes")
        points(ptsa$y, ptsa$x, pch=16, col=2)
        points(ptsb$y, ptsb$x, pch=16, col=3)
        
        if (missing(pos)) pos <- "topleft"
        legend(pos, , c("above", "below"), col=2:3, pch=16, bty="n", inset=0.05)
      }
      lines(y, x, lwd=2, col=4)
    }
  }
  
  return(list(line=list(x=x, y=y), pts.above=ptsa, pts.below=ptsb,
              spline=ospl$spl, spline.df=ospl$df, spline.err=ospl$error))
}

# plot thresholds, possibly bounded by each other
thrplot <- function (filenames, varname, varvalues, epsylons=NULL, boundid,
                     boundhi=TRUE, areadown=TRUE, verbose=FALSE, plot=TRUE,
                     mirror=TRUE, linejoin=TRUE, joinright=TRUE, pos=NULL, ...) {
  # 'mirror', linejoin', 'joinright' and 'pos' are aesthetic features
  
  n <- length(filenames)
  
  if (is.null(epsylons)) epsylons <- rep(0, n)
  if (length(epsylons) == 1) epsylons <- rep(epsylons, n)
  if (missing(boundid)) boundid <- rep(NA, n)
  if (length(boundhi) == 1) boundhi <- rep(boundhi, n)
  if (is.null(pos)) pos <- "topleft"
  
  load(filenames[1]) # to know range1 and range2
  
  if (!mirror) {
    rangeX <- range1
    rangeY <- range2
  } else {
    rangeX <- range2
    rangeY <- range1
  }
  
  xs <- exp(seq(log(rangeX[1]), log(rangeX[2]), , 100))
  thrs <- lns <- nds <- vector("list", n)
  for (i in 1:n) {
    bid <- boundid[i]
    f <- filenames[i]
    e <- epsylons[i]
    
    if (is.na(bid)) thr <- calcthr(f, e, x=xs, ...) else
      if (boundhi[i]) thr <- calcthr(f, e, x=xs, spline.hi=thrs[[bid]], ...) else
        thr <- calcthr(f, e, x=xs, spline.lo=thrs[[bid]], ...)
    
    thrs[[i]] <- thr$spline
    lns[[i]] <- thr$line
    nds[[i]] <- list(xa=thr$pts.above$x, ya=thr$pts.above$y,
                     xb=thr$pts.below$x, yb=thr$pts.below$y)
  }
  rm(thrs)
  
  b <- 0.4
  areaback <- rgb(1, b, b)
  areacols <- rgb(b, seq(100, 255, , n+1)/255, b)[-1]
  
  ys <- vector("list", n)
  for (i in 1:n)
    ys[[i]] <- lns[[i]]$x
  
  if (linejoin) {
    for (i in which(!is.na(boundid))) {
      bid <- boundid[i]
      yto <- ys[[bid]]
      yi <- ys[[i]]
      
      if (joinright)
        sgn <- sign(tail(yto - yi, 1)) else
          sgn <- sign(head(yto - yi, 1))
      
      if (sgn > 0)
        closeto <- which(yto / yi < 20/19) else
          closeto <- which(yto / yi > 19/20)
      
      if (length(closeto) > 0) {
        if (joinright)
          whcl <- seq(max(closeto), 1) else
            whcl <- seq(1, max(closeto))
          d <- abs(yto - yi)[whcl]
          ys[[i]][whcl] <- yto[whcl] - sgn * cummin(d)
      }
    }
  }
  
  # outline areas
  areas <- vector("list", n)
  for (i in 1:n) {
    if (areadown) {
      x <- c(xs, rangeX[2], rangeX[2], rangeX[1])
      y <- c(ys[[i]], rangeY[2], rangeY[1], rangeY[1])
      
    } else {
      x <- c(xs, rangeX[2], rangeX[1], rangeX[1])
      y <- c(ys[[i]], rangeY[2], rangeY[2], rangeY[1])
    }
    
    areas[[i]] <- list(x=x, y=y)
  }
  
  # plotting
  if (plot) {
    plot(rangeX, rangeY, type="n", log="xy",
         xlab="copies", ylab="genes", main="Sustainability limit")
    drawlogticks(rangeX, 1)
    drawlogticks(rangeY, 2)
    
    rect(rangeX[1], rangeY[1], rangeX[2], rangeY[2], col=areaback, border=NA)
    
    for (i in 1:n)
      polygon(areas[[i]]$x, areas[[i]]$y, col=areacols[i], border=NA)
    
    for (i in 1:n)
      lines(xs, ys[[i]], lty=i, lwd=2)
    
    legend(pos, , varvalues, lty=1:n, title=varname, lwd=2, bty="n", inset=0.05)
    
    # margin
    xoff <- c(rangeX[1] /2, rangeX, rangeX[2] *2)
    yoff <- c(rangeY[1] /2, rangeY, rangeY[2] *2)
    rect(xoff[1], yoff[1], xoff[2], yoff[4], col=0, border=NA)
    rect(xoff[3], yoff[1], xoff[4], yoff[4], col=0, border=NA)
    rect(xoff[1], yoff[1], xoff[4], yoff[2], col=0, border=NA)
    rect(xoff[1], yoff[3], xoff[4], yoff[4], col=0, border=NA)
    box()
  }
  
  if (verbose) return(list(range1=range1, range2=range2,
                           areas=areas, nodes=nds, lines=lns))
}

maxgenesplot <- function (files, varname, varvalues, epsylons, below=TRUE,
                          verbose=FALSE, plot=TRUE, plotplace=FALSE, ...) {
  
  x <- as.numeric(varvalues)
  thr <- thrplot(files[x > 0], varname, varvalues[x > 0], epsylons[x > 0],
                 verbose=TRUE, plot=plotplace, ...)
  
  x <- x[x > 0]
  ymaxs <- array(, c(length(x), 3)) # ymax, x @ ymax, is y monotonous and growing
  for (i in seq_along(x)) {
    bnd <- thr$nodes[[i]] # points above and below boundary (below: surviving)
    if (below) {
      yi <- bnd$xb; xi <- bnd$yb} else {
        yi <- bnd$xa; xi <- bnd$ya}
    
    wm <- which.max(yi)
    mongrow <- all(yi == cummax(yi)) &
      tail(yi, 1) > head(yi, 1)
    ymaxs[i, ] <- c(yi[wm], xi[wm], mongrow)
  }
  
  y <- ymaxs[, 1]
  g <- ymaxs[, 3]
  if (plotplace)
    points(ymaxs[, 2], y, pch=c(16, 1)[g+1])
  
  # plot
  if (plot) {
    range1 <- thr$range1
    range2 <- range(x)
    mn <- paste("Effect of", varname, "on information integration")
    plot(x, y, xlim=range2, ylim=range1, log="xy", type="o", pch=c(16, 1)[g+1],
         main=mn, xlab=varname, ylab="max. genome size")
    drawlogticks(range2, 1)
    drawlogticks(range1, 2)
    
    if (below) ttl <- NULL else ttl <- "extinction side of boundary"
    leg <- c("approximation", "lower bound")
    legend("topright", , leg, title=ttl, pch=c(16, 1), bty="n", inset=0.05)
  }
  
  if (verbose) return(list(molerror=x, maxgenes=ymaxs[, 1], below=below))
}

drawlogticks <- function (values, side=1, base=10) {
  
  values <- values[values > 0]
  rng <- range(values)
  
  logrng <- log(rng, base)
  if (any(logrng %% 1 != 0))
    logrng <- c(floor(logrng[1]), ceiling(logrng[2]))
  
  lmin <- logrng[1]
  if (lmin < 0) logrng <- logrng - lmin
  
  ordmagn <- base^seq(logrng[1], logrng[2])
  tickpos <- c(1:9 %o% ordmagn)
  
  if (lmin < 0) tickpos <- tickpos * base^lmin
  tickpos <- tickpos[which(tickpos >= rng[1] & tickpos <= rng[2])]
  
  axis(side, tickpos, label=FALSE)
}

simvideo <- function (tau, v, g=100, N=100, mu=0.05, apar=1.1,
                      path="images/popvideos", folder=1) {
  
  fname <- paste0(path, "/", folder)
  while (file.exists(fname)) {
    folder <- folder +1
    fname <- paste0(path, "/", folder)
  }
  dir.create(fname)
  
  png(paste0(fname, "/pop_1_%03d.png"), 600, 600, , 25)
  out <- sim(g, N, tau, v, mu, apar=apar, draw=TRUE, drawpop=TRUE)
  dev.off()
}

calchistmeans <- function (histories, nseries=3, snames=c("msurv", "mfitn", "mpar")) {
  
  lh <- length(histories)
  m2N <- matrix(, lh, nseries, dimnames=list(NULL, snames))
  for (i in seq_len(lh)) {
    
    hi <- histories[[i]]
    lhi <- length(hi)
    mnode <- matrix(, lhi, nseries)
    for (j in seq_len(lhi)) {
      
      hij <- hi[[j]]
      nr <- nrow(hij)
      if (!is.na(hij[nr, 1])) {
        half <- seq(nr %/% 2 +1, nr)
        mnode[j, ] <- colMeans(hij[half, , drop=FALSE], na.rm=TRUE)}
    }
    m2N[i, ] <- colMeans(mnode, na.rm=TRUE)
  }
  
  return(m2N)
}

# history data (hst) describes the 2N population
sim2N <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N*2, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- sampleSF(N*2, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      pr <- mean(vspop2[, 1] / rowSums(vspop2), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N*2), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

# extra output: fitness of 2N offsprings, parasites of N survivors
sim_verbose <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  fitns2N <- array(, c(g, N*2))
  paras <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    fitns2N[j, ] <- fitn
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- sampleSF(N*2, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      paras[j, ] <- vspop1[, 1] / rowSums(vspop1)
      pr <- mean(paras[j, ], na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))
        fitns2N <- rbind(fitns2N, array(, c(gsta, N*2)))
        paras <- rbind(paras, array(, c(gsta, N)))
      }
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs, fitns2N, paras))
}

sim_test <- function (g, N, tau, v, mu) {
  out <- log(tau)^2 < v
  hst <- 1:100 + tau
  if (exists("errq"))
    if (runif(1) < errq) out <- !out
  return(list(c(out, NA, NA), hst))
}

sim_rand <- function (g, N, tau, v, mu) {
  s <- round(v*100, 0) * 100 + round(tau*100, 0)
  if (exists("const")) s <- s + const
  set.seed(s)
  out <- sample(c(T, F), 1)
  hst <- 1:100
  return(list(c(out, NA, NA), hst))
}

drawline <- function (x1, y1, x2, y2, ...) {
  b <- (y2 - y1) / (x2 - x1)
  a <- y1 - x1 * b
  abline(a, b, untf=TRUE, ...)
  return(c(a, b))
}

cmplot_vmax <- function (range1, range2, nodes, squares, int=0, tlog=3, border=FALSE,
                         lazy=TRUE, mirror=FALSE, repcol=9, ...) {
  nodes[, 3] <- nodes[, 2] * nodes[, 3] * 2
  range2 <- range1 * range2 * 2
  labs <- c("genes", "vmax")
  
  tlg <- c("", "x", "y", "xy")[tlog %% 4 +1]
  coo <- list(integer(0), c(1, 3), c(2, 4), 1:4)[[int %% 4 + 1]]
  nodes[, coo] <- round(nodes[, coo])
  squares[, coo] <- round(squares[, coo])
  if (lazy) nodes <- nodes[which(nodes[, repcol] > 0), , drop=FALSE]
  
  if (mirror) {
    rngs <- rbind(range1, range2)
    range1 <- rngs[2, ]
    range2 <- rngs[1, ]
    nodes[, 2:3] <- nodes[, 3:2]
    labs <- labs[2:1]}
  
  sqnodes <- unique(rbind(squares[, c(1, 2)], squares[, c(1, 4)],
                          squares[, c(3, 2)], squares[, c(3, 4)]))
  ndpos_sq <- paste(sqnodes[, 1], sqnodes[, 2])
  ndpos_nd <- paste(nodes[, 2], nodes[, 3])
  redund <- ndpos_nd %in% ndpos_sq
  snodes <- nodes[!redund, , drop=FALSE]
  
  uncert <- nodes[, repcol] > min(nodes[, repcol], Inf, na.rm=TRUE)
  unodes <- nodes[uncert, , drop=FALSE]
  
  rnodes <- nodes[redund & !uncert, , drop=FALSE]
  
  plot(range1, range2, type="n", log=tlg,
       xlab=labs[1], ylab=labs[2], ...)
  ifelse(border, brd <- par("fg"), brd <- NA)
  for (i in seq_len(nrow(squares))) {
    sx <- c(squares[i, 1], squares[i, 3])[c(1, 1, 2, 2)]
    sy <- c(squares[i, 2], squares[i, 4])[c(1, 2, 2, 1)]
    
    if (!mirror)
      polygon(sx, sy * sx * 2, col=squares[i, 5] + 2, border=brd) else
        polygon(sy * sx * 2, sx, col=squares[i, 5] + 2, border=brd)
  }
  
  cols <- snodes[, 4] + 2
  points(snodes[, 2], snodes[, 3], col=cols, pch=16)  # separate nodes
  points(unodes[, 2], unodes[, 3])                    # uncertain nodes
  points(rnodes[, 2], rnodes[, 3], pch=".", cex=3)    # rest
}

# ---
# additional functions for variant models

sim_degrade <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, d=0.05, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times, and degrade
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    vspop2 <- vspop2 - sapply(vspop2, function (i) rbinom(1, i, d))
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- sampleSF(N*2, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

sim_vesmort <- function (g, N, tau, v, mu, o, mort=0, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    if (mort > 0)
      vspop2[sample(2*N, rbinom(1, N*2, mort)), ] <- 0
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- sampleSF(N*2, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

sim_parasitoid <- function (g, N, tau, v, mu, o, eps=1, npar=0, ndiff=1, apar=1.1, adiff=1.1,
                            pop=0.1, tra=0.1, mort=0, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+2
  vmax <- v*tau*2
  a <- c(apar, adiff, rep(1, tau))
  
  if (npar == -1) npar <- v
  if (ndiff == -1) ndiff <- v
  vid <- c(npar, ndiff, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 4))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  pd <- mean(vspop1[, 2] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-(1:2)], tau, vmax, eps), pr, pd)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    if (pop > 0 || mort > 0) {
      # parasitoids bursting vesicles
      m <- sapply(vspop2[, 2], function (i) rbinom(1, i, pop)) > 0
      m[sample(N*2, rbinom(1, N*2, mort))] <- TRUE
      if (sum(m) > 0) {
        pool <- colSums(vspop2[m, , drop=FALSE])
        vspop2[m, ] <- 0
        
        # pool molecules getting back into vesicles
        l <- sum(!m)
        if (l > 0) {
          mpool <- unlist(lapply(seq_len(tp), function (i) rep(i, pool[i])))
          mpool <- mpool[runif(length(mpool)) < tra]
          nmp <- length(mpool)
          vsaim <- which(!m)[sample(l, nmp, TRUE)]
          for (k in seq_len(nmp))
            vspop2[vsaim[k], mpool[k]] <- vspop2[vsaim[k], mpool[k]] +1
          
          # vesicles above critical size fission
          os <- rowSums(vspop2) > vmax
          if ((no <- sum(os)) > 0) {
            wo <- which(os)
            vstofission <- cbind(1, vspop2[wo, , drop=FALSE])
            vsfissioned <- array(, c(0, tp+1))
            
            while (length(vstofission) > 0) {
              vsf <- fission(vstofission[1, , drop=FALSE], tp+1)
              wvsf <- rowSums(vsf) <= vmax
              vsfissioned <- rbind(vsfissioned, vsf[wvsf, ])
              vstofission <- rbind(vstofission[-1, ], vsf[!wvsf, ])
            }
            
            # one offspring vesicle per parent always survives
            wf <- which(vsfissioned[, 1] == 1)
            vspop2[wo, ] <- vsfissioned[wf, -1]
            vsfissioned <- vsfissioned[-wf, -1, drop=FALSE]
            
            # the remaining offsprings may take the place of burst vesicles
            nvf <- nrow(vsfissioned)
            vsnew <- min(sum(m), nvf)
            vpi <- which(m)[sample(sum(m), vsnew)]
            voi <- sample(nvf, vsnew)
            vspop2[vpi, ] <- vsfissioned[voi, ]
          }
        }
      }
    }
    
    fitn <- apply(vspop2[, -(1:2), drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- sampleSF(N*2, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      pd <- mean(vspop1[, 2] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr, pd)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 4)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio", "parasitoid ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
    plot(0:g, hst[, 4], xlab=ttl[1], ylab=ttl[5], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

fusion <- function (vs, tau, vmax) {
  vs2 <- fission(colSums(vs), tau)
  
  # oversized vesicle
  o <- which(rowSums(vs2) > vmax)
  while (length(o) > 0) {
    vs2[o, ] <- fission(vs2[o, , drop=FALSE], tau)[sample(2, 1), ]
    o <- which(rowSums(vs2) > vmax)}
  
  return(vs2)
}

sim_fusion <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE, pfus=0) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    # fusion & fission
    nfus <- min(N, rbinom(1, N, pfus))
    ifus <- matrix(sample(seq_len(N*2)), 2)[, seq_len(nfus), drop=FALSE]
    for (i in seq_len(nfus))
      vspop2[ifus[, i], ] <- fusion(vspop2[ifus[, i], , drop=FALSE], tp, vmax)
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- sampleSF(N*2, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

checkdiff <- function (vs, a, mu, test) {
  vp <- vs[1]
  ve <- vs[-1]
  ap <- a[1]
  ae <- a[2]
  ve_n <- ve * exp(ae * (1 - mu) * test)
  co <- sum(ve)*ae*mu / (ae*(1-mu) - ap)
  vp_n <- (vp - co) * exp(ap * test) + co * exp(ae * (1-mu) * test)
  return(c(vp_n, ve_n))
}

growD <- function (vs, vmax, tau, a, mu, tp, err=1e-5) {
  
  tout <<- tau
  vsout <<- vs
  vmaxout <<- vmax
  
  if (0 < sum(vs) && sum(vs) < vmax) {
    ts <- log(vmax / sum(vs)) / c(max(a), min(a) * (1 - mu))
    ds <- vmax - rowSums(rbind(checkdiff(vs, a, mu, ts[1]),
                               checkdiff(vs, a, mu, ts[2])))
    
    # finding t
    while(!any(abs(ds) < err)) {
      tnew <- mean(ts)
      dnew <- vmax - sum(checkdiff(vs, a, mu, tnew))
      
      i <- which(sign(ds) == sign(dnew))
      if (length(i) == 0) stop("@growD_diff")
      ts[i] <- tnew
      ds[i] <- dnew
    }
    
    i <- which.min(abs(ds))
    vs <- trunc(checkdiff(vs, a, mu, ts[i]))
  }
  
  return(vs)
}

fissionD <- function (vs, tau) { # deterministic
  # see 'fission'
  
  vs2 <- vs %/% 2
  return(rbind(vs2, vs2))
}

sim_selectbest <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- tail(order(fitn), N)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

growS_proc <- function (vs, proc, tau, a, mu, tp) { # stochastic, processive
  # vs: vector of copy numbers
  # proc: number of molecular replications
  # tau: number of genes
  # a: vector of gene affinities
  # mu: error rate
  # tp: tau+1 (w/ parasites)
  
  if (sum(vs) > 0) {
    for (i in seq_len(proc)) {
      x <- sample1F(tp, a*vs)
      if (runif(1) < mu) x <- 1
      vs[x] <- vs[x] +1
    }
  }
  
  return(vs)
}

sim_proc <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2, proc is v*tau)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  proc <- v*tau
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitness(vid[-1], tau, vmax, eps), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], proc, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitness(i, tau, vmax, eps))
    sf <- sum(fitn > 0)
    if (sf > 0) {
      surv <- sampleSF(N*2, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}

simsl_viable <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop <- array(rep(vid, each=N), c(N, tp))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  fitn <- apply(vspop[, -1, drop=FALSE], 1,
                function (i) fitness(i, tau, vmax, eps))
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop[, 1] / rowSums(vspop))
  hst[1, ] <- c(N, fitn[1], pr)
  
  while (inf.loss == FALSE && j <= g*N) {
    
    vsums <- rowSums(vspop) # tmp: no need to recalculate every time
    vsm <- max(vsums)
    while (vsm < vmax) { # max(vsums) < vmax
      i <- sample.int(N, 1, , fitn)
      vspop[i, ] <- growsl(vspop[i, , drop=FALSE], vmax, tau, a, mu, tp)
      fitn[i] <- fitness(vspop[i, -1, drop=FALSE], tau, vmax, eps)
      vsums[i] <- vsm <- vsums[i] +1
    }
    
    if (any(fitn == 0))
      h <- which(fitn == 0)[1] else
        h <- sample.int(N, 1) # outflux
    vspop[c(i, h), ] <- fission(vspop[i, , drop=FALSE], tp)
    fitn[i] <- fitness(vspop[i, -1, drop=FALSE], tau, vmax, eps)
    fitn[h] <- fitness(vspop[h, -1, drop=FALSE], tau, vmax, eps)
    
    sf <- sum(fitn > 0)
    if (sf > 0) {
      if (j %% N == 0) {
        pr <- mean(vspop[, 1] / rowSums(vspop), na.rm=TRUE)
        hst[j %/% N +1, ] <- c(sf, mean(fitn), pr)}
      j <- j +1
      
      if (j > g*N & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst))
}

fitnessNE_min <- function (vs, tau, vmax, eps, b0=0.1) {
  return(prod(pmax(b0, vs) * tau/vmax))
}

fitnessNE_sum <- function (vs, tau, vmax, eps, b0=0.1) {
  return(prod( (vs + b0)/(vmax + tau*b0) * tau))
}

simNE <- function (g, N, tau, v, mu, o, eps=1, npar=0, apar=1.1, adiff=NULL, b0=0.1, draw=FALSE) {
  # g: number of generations
  # N: number of vesicles
  # tau: number of genes
  # v: 'starting' copy number (vmax is v*tau*2)
  # eps: exponent of quality
  # adiff: if !NULL, the 'one-different' affinity (base affinity is 1)
  # both fission and growth is stochastic
  
  gsta <- g
  gmax <- g * 10
  tp <- tau+1
  vmax <- v*tau*2
  if (npar == -1) npar <- v
  vid <- c(npar, rep(v, tau))
  vspop1 <- array(rep(vid, each=N), c(N, tp))
  vspop2 <- array(, c(N*2, tp))
  survs <- array(, c(g, N))
  a <- c(apar, rep(1, tau))
  if (!is.null(adiff)) a[2] <- adiff
  
  inf.loss <- FALSE # information loss ends the simulation
  j <- 1            # otherwise, simulation goes on for 'g' generations
  
  hst <- array(, c(g+1, 3))
  pr <- mean(vspop1[, 1] / rowSums(vspop1))
  hst[1, ] <- c(N, fitnessNE(vid[-1], tau, vmax, eps, b0), pr)
  
  while (!inf.loss && j <= g) {
    
    # instant growth and fission, N times
    for (i in seq_len(N)) {
      vspop1[i, ] <- grow(vspop1[i, , drop=FALSE], vmax, tau, a, mu, tp)
      vspop2[i*2 - 1:0, ] <- fission(vspop1[i, , drop=FALSE], tp)
    }
    
    fitn <- apply(vspop2[, -1, drop=FALSE], 1,
                  function (i) fitnessNE(i, tau, vmax, eps, b0))
    cons <- apply(vspop2[, -1, drop=FALSE] > 0, 1, prod)
    sf <- sum(cons > 0)
    if (sf > 0) {
      surv <- sampleSF(N*2, N, fitn)
      survs[j, ] <- surv
      vspop1 <- vspop2[surv, , drop=FALSE]
      sf <- min(N, sf)
      pr <- mean(vspop1[, 1] / rowSums(vspop1), na.rm=TRUE)
      hst[j+1, ] <- c(sf, mean(fitn[surv]), pr)
      j <- j +1
      
      if (j > g & gmax > g & !isequi(hst[, 2])) {
        g <- g + gsta
        hst <- rbind(hst, array(, c(gsta, 3)))
        survs <- rbind(survs, array(, c(gsta, N)))}
    } else {inf.loss <- TRUE}
  }
  
  if (draw) {
    ttl <- c("generations", "survivors", "fitness", "parasite ratio")
    frng <- range(hst[, 2], na.rm=TRUE)
    plot(0:g, hst[, 1], xlab=ttl[1], ylab=ttl[2], ylim=c(0, N), type="l")
    plot(0:g, hst[, 2], xlab=ttl[1], ylab=ttl[3], ylim=frng, type="l", log="y")
    plot(0:g, hst[, 3], xlab=ttl[1], ylab=ttl[4], ylim=c(0, 1), type="l")
  }
  
  fts <- hst[, 2]
  efit <- mean(fts[seq(g %/% 2, g) +1]) # NA if inf.loss
  eq <- isequi(fts)
  return(list(c(sf, efit, eq), hst, survs))
}
