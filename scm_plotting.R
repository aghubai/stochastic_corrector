
rm(list=ls())

setwd(choose.dir())
source("scm_functions.R")

# ---

f <- paste0("demo data/", c("reference",
                            "a_p.67", "a_p.95", "a_p1", "a_p1.5", "a_p2",
                            "pop3", "pop10", "pop30", "pop300", "pop1000"), ".out")

vn <- c("reference", rep("parasite replication rate", 5), rep("population size", 5))

vv <- c(NA, "0.67", "0.95", "1", "1.5", "2", "3", "10", "30", "300", "1000")

# ---
# basic model, changing parasite replication affinity

n <- c(2:4, 1, 5:6)
files <- f[n]; varname <- vn[n[1]]; varvalues <- vv[n]
varvalues[which(n==1)] <- "0.05"

bnd <- c(NA, 1, 2, 3, 4, 5)
thrplot(files, varname, varvalues, boundid=bnd)

# ---
# population size

n <- rev(c(7:9, 1, 10:11))
files <- f[n]; varname <- vn[n[1]]; varvalues <- vv[n]
varvalues[which(n==1)] <- "100"

thrplot(files, varname, varvalues)

maxgenes()

draweach()
