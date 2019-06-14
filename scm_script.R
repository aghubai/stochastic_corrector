
rm(list=ls())
library("hash")

setwd(choose.dir())
source("scm_functions.R")

simtype <- 0

# ---

# basic model, changing mutation rate
if (simtype == 1) {
  simwrapper("BRNEL")
  simwrapper("BRPEL")
  simwrapper("BR1EL")
  simwrapper("BR2EL")
  simwrapper("BR3EL")
  simwrapper("BR4EL")
}

# basic model, changing parasite replication affinity
if (simtype == 2) {
  simwrapper("BRPBL")
  simwrapper("BRPCL")
  simwrapper("BRPDL")
  simwrapper("BRPEL")
  simwrapper("BRPFL")
  simwrapper("BRPGL")
}

# basic model, rate difference between templates (one vs the rest)
if (simtype == 3) {
  simq <- sim
  fitness <- fitnessR
  mu <- 0.05
  apar <- 1
  
  n <- 1
  adiff <- c(2/3, 3/2, 0.5, 2)[n]
  l <- TRUE
  o <- 2
  
  inp <- paste0("BRPDL_adiff", n)
  cat(paste0("Simulation ", inp, " started at ", strftime(Sys.time()), "\n"))
  
  g <- 100 # number of generations
  N <- 100 # number of cells
  range1 <- c(1, 100) # number of types (not counting the parasite)
  range2 <- c(1, 100) # number of copies (at start)
  
  out <- dynscan(range1, range2, lbx="genes", lby="copies", lbm=inp,
                 int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o, apar=apar, adiff=adiff)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, file=fnamel)
}

# basic model, population size
if (simtype == 4) {
  simq <- sim
  
  mu <- 0.05
  apar <- 1.1
  
  l <- TRUE
  o <- 2
  
  N <- 10 # number of cells
  inp <- paste0("BRPEL_popsize_N", N)
  cat(paste0("Simulation ", inp, " started at ", strftime(Sys.time()), "\n"))
  
  g <- 100 # number of generations
  range1 <- c(1, 100) # number of types (not counting the parasite)
  range2 <- c(1, 100) # number of copies (at start)
  
  out <- dynscan(range1, range2, lbx="genes", lby="copies", lbm=inp,
                 int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o, apar=apar)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, file=fnamel)
}

# molecular degradation
if (simtype == 5) {
  simq <- sim_degrade
  
  g <- 100 # number of generations
  N <- 100 # number of cells
  range1 <- c(1, 100) # number of types (not counting the parasite)
  range2 <- c(1, 100) # number of copies (at start)
  
  l <- TRUE
  mu <- 0.05
  apar <- 1.1
  o <- 2
  
  d <- 0.01
  inp <- paste0("degrade_d", d)
  cat(paste0("Simulation ", inp, " started at ", strftime(Sys.time()), "\n"))
  
  out <- dynscan(range1, range2, 7, lbx="genes", lby="copies", lbm=inp,
                 int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o, apar=apar, d=d)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, mu, apar, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, mu, apar, file=fnamel)
}

# burst variant w/o reuptake
if (simtype == 6) {
  simq <- sim_vesmort
  
  g <- 100 # number of generations
  N <- 100 # number of cells
  range1 <- c(1, 100) # number of types (not counting the parasite)
  range2 <- c(1, 100) # number of copies (at start)
  
  l <- TRUE
  mu <- 0.05 # 0.01
  apar <- 1.1 # 1.5, 2
  o <- 2
  
  mort <- 0.01
  inp <- paste0("vesmort_d", mort)
  
  out <- dynscan(range1, range2, 7, lbx="genes", lby="copies", lbm=inp,
                 int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o, apar=apar, mort=mort)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, mu, apar, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, mu, apar, file=fnamel)
}

# burst variant w/ reuptake
if (simtype == 7) {
  simq <- sim_parasitoid
  
  g <- 100 # number of generations
  N <- 100 # number of cells
  range1 <- c(1, 100) # number of types (not counting the parasite)
  range2 <- c(1, 100) # number of copies (at start)
  
  l <- TRUE
  mu <- 0.05
  apar <- 1.1
  o <- 2
  
  npar <- 0
  ndiff <- 0    # no parasitoids
  
  vmort <- c(0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  vtra <- c(0, 0, 0.5, 1, 0, 0.5, 1, 1, 1, 1, 1, 1, 1, 1)
  
  i <- 1
  mort <- vmort[i]
  tra <- vtra[i]
  inp <- paste0("vesmortupt_nopd_d", mort, "_r", tra)
  
  if (FALSE) {
    g <- 1000
    inp <- paste0("vesmortupt_nopd_g1000_d", mort, "_r", tra)
  }
  
  out <- dynscan(range1, range2, 7, lbx="genes", lby="copies", lbm=inp,
                 int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o,
                 npar=npar, ndiff=ndiff, apar=apar, tra=tra, mort=mort)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, mu, apar, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, mu, apar, file=fnamel)
}

# fusion variant
if (simtype == 8) {
  g <- 100 # number of generations
  N <- 100 # number of cells
  mu <- 0.05 # error rate per molecule
  pfus <- 0.5
  simq <- sim_fusion
  
  range1 <- c(1, 100) # tau
  range2 <- c(1, 100) # v
  
  inp <- paste0("f", round(pfus*100), "_mu", round(mu*100))
  out <- dynscan(range1, range2, scalemax=7, int=3, draw=TRUE, debug=TRUE,
                 g=g, N=N, mu=mu, pfus=pfus, lbx="genes", lby="copies", lbm=inp)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, mu, pfus, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, mu, pfus, file=fnamel)
}

# offspring variant
if (simtype == 9) {
  simwrapper("BRPEL0")
  simwrapper("BRPEL3")
  simwrapper("BRPEL8")
}

# deterministic variant
if (simtype == 10) {
  
  # deterministic growth
  if (simsubtype == "a") {
    inp <- "gDfS"
  }
  
  # deterministic fission
  if (simsubtype == "b") {
    inp <- "gSfD"
  }
  
  # deterministic growth and fission
  if (simsubtype == "c") {
    inp <- "gDfD"
  }
  
  simq <- sim
  
  gftype <- strsplit(strsplit(inp, "g")[[1]][2], "f")[[1]]
  grow <- get(paste0("grow", gftype[1]))
  fission <- get(paste0("fission", gftype[2]))
  
  g <- 100 # number of generations
  N <- 100 # number of cells
  range1 <- c(1, 100) # number of types (not counting the parasite)
  range2 <- c(1, 100) # number of copies (at start)
  
  l <- TRUE
  mu <- 0.05 # 0.01
  apar <- 1.1 # 1.5, 2
  o <- 2
  
  out <- dynscan(range1, range2, 7, lbx="genes", lby="copies", lbm=inp,
                 int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o, apar=apar)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, mu, apar, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, mu, apar, file=fnamel)
}

# non-probabilistic selection (cf. "prospective value")
if (simtype == 11) {
  simq <- sim_selectbest
  fitness <- fitnessR
  mu <- 0.05
  apar <- 1.1
  l <- TRUE
  o <- 2
  
  inp <- "BRPEL_selectbest"
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
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, file=fnamel)
}

# fission based on processivity
if (simtype == 12) {
  simq <- sim_proc
  grow <- growS_proc
  mu <- 0.05
  apar <- 1.1
  l <- TRUE
  o <- 2
  
  inp <- "BRPEL_proc"
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
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, file=fnamel)
}

# selection type variant
if (simtype == 13) {
  # a) growth-based selection
  simwrapper("SRNEL")
  simwrapper("SRPEL")
  
  # b) growth-based selection w/ viable advantage
  if (simsubtype == "b") {
    simq <- simsl_viable
    
    g <- 100 # number of generations
    N <- 100 # number of cells
    range1 <- c(1, 100) # number of types (not counting the parasite)
    range2 <- c(1, 100) # number of copies (at start)
    
    l <- TRUE
    mu <- 0.05
    apar <- 1.1
    o <- 2
    
    inp <- "SRPEL_viable"
    cat(paste0("Simulation ", inp, " started at ", strftime(Sys.time()), "\n"))
    
    out <- dynscan(range1, range2, 7, lbx="genes", lby="copies", lbm=inp,
                   int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o, apar=apar)
    
    nodes <- out[[1]]
    squares <- out[[2]]
    histories <- out[[3]]
    
    fname <- paste0(inp, ".out")
    save(range1, range2, nodes, squares, mu, apar, file=fname)
    
    fnamel <- paste0(inp, "_large.out")
    save(range1, range2, nodes, squares, histories, mu, apar, file=fnamel)
  }
  
  # c) no selection
  simwrapper("XRNEL")
  simwrapper("XRPEL")
}

# nonessential templates
if (simtype == 14) {
  
  # substitutive
  if (simsubtype == "a") {
    fitnessNE <- fitnessNE_min; fn <- "_min"
  }
  
  # additive
  if (simsubtype == "b") {
    fitnessNE <- fitnessNE_sum; fn <- "_sum"
  }
  
  simq <- simNE
  b0 <- 0.1
  
  mu <- 0.05
  apar <- 1.1
  l <- TRUE
  o <- 2
  
  inp <- paste0("BRPEL_nonesn_", b0, fn)
  cat(paste0("Simulation ", inp, " started at ", strftime(Sys.time()), "\n"))
  
  g <- 100 # number of generations
  N <- 100 # number of cells
  range1 <- c(1, 100) # number of types (not counting the parasite)
  range2 <- c(1, 100) # number of copies (at start)
  
  out <- dynscan(range1, range2, lbx="genes", lby="copies", lbm=inp,
                 int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o, apar=apar, b0=b0)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, file=fnamel)
}

# mean-field variant
if (simtype == 15) {
  simq <- sim_nsel
  
  if (TRUE) {
    mu <- 0.05
    apar <- 1.1
    inp <- "mfBRPEL"
  }
  
  if (FALSE) {
    mu <- 0.05
    apar <- 1
    inp <- "mfBRPDL"
  }
  
  if (FALSE) {
    mu <- 0
    apar <- 1.1
    inp <- "mfBRNEL"
  }
  
  l <- TRUE
  o <- 2
  
  cat(paste0("Simulation ", inp, " started at ", strftime(Sys.time()), "\n"))
  
  g <- 10000 # number of generations
  N <- 1 # number of cells
  range1 <- c(1, 100) # number of types (not counting the parasite)
  range2 <- c(1, 100)*100 # number of copies (at start)
  
  out <- dynscan(range1, range2, lbx="genes", lby="copies", lbm=inp,
                 int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o, apar=apar)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, file=fnamel)
}

# extra redundancy
if (simtype == 16) {
  simq <- sim
  fitness <- fitnessR
  mu <- 0.05
  apar <- 1.2
  
  l <- TRUE
  o <- 2
  
  inp <- "BRPHL_redund"
  cat(paste0("Simulation ", inp, " started at ", strftime(Sys.time()), "\n"))
  
  g <- 100 # number of generations
  N <- 100 # number of cells
  range1 <- c(1, 25) # number of types (not counting the parasite)
  range2 <- c(1, 500) # number of copies (at start)
  
  out <- dynscan(range1, range2, lbx="genes", lby="copies", lbm=inp,
                 int=3, lazy=l, draw=TRUE, debug=TRUE, g=g, N=N, mu=mu, o=o, apar=apar)
  
  nodes <- out[[1]]
  squares <- out[[2]]
  histories <- out[[3]]
  
  fname <- paste0(inp, ".out")
  save(range1, range2, nodes, squares, file=fname)
  
  fnamel <- paste0(inp, "_large.out")
  save(range1, range2, nodes, squares, histories, file=fnamel)
}
