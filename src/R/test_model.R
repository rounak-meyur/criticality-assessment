#install.packages("", repos="http://cran.r-project.org")
library(lpSolveAPI)
library(stringr)
library(pracma)
library(jsonlite)
library(rjson)
library(RJSONIO)
library(RCurl)
library(httr)
library(dplyr)
library(igraph)
library(bigmemory)
library(Matrix)
library(spatstat.utils)
library(base)

#---------------------------------------------------------------------------------------------------------------------------------------------------
get.data = function(bus, gen, branch) {
  gens = matrix(nrow = length(unique(gendat$bus_id)), ncol = 3)
  colnames(gens) = c("bus_id", "pmax", "pmin")
  index = 1
  gens[1, 1] = gendat$bus_id[1]
  gens[1, 2] = gendat$pmax[1]
  gens[1, 3] = gendat$pmin[1]
  for (l in 2:nrow(gendat)) {
    if (gendat$bus_id[l] == gendat$bus_id[l - 1]) {
      gens[index, 2] = gens[index, 2] + gendat$pmax[l]
      gens[index, 3] = gens[index, 3] + gendat$pmin[l]
    } else {
      index = index + 1
      gens[index, 1] = gendat$bus_id[l]
      gens[index, 2] = gendat$pmax[l]
      gens[index, 3] = gendat$pmin[l]
    }
  }
  
  g_max = c()
  g_min = c()
  index = 1
  for (l in 1:nrow(busdat)) {
    if (index <= 485) {
      if (busdat$bus_id[l] == gens[index, 1]) {
        g_max = append(g_max, gens[index, 2] / 100.0)
        g_min = append(g_min, gens[index, 3] / 100.0)
        index = index + 1
      } else {
        g_max = append(g_max, 0)
        g_min = append(g_min, 0)
      }
    } else {
      g_max = append(g_max, 0)
      g_min = append(g_min, 0)
    }
  }
  
  d_min = rep(0, nrow(bus))
  d_max = bus$pd / 100.0
  
  flim = branch$rateA / 100.0
  
  return(list("gmin" = g_min, "gmax" = g_max, "dmax" = d_max, "dmin" = d_min, "flim" = flim))
}

calculate.s = function(bus, branch, line_r) {
  g = graph.data.frame(branchdat[, 1:2], directed = FALSE)
  size = gsize(g)
  
  if (line_removal != 0) {
    g = delete_edges(g, c(line_removal))
  }
  
  gs = decompose(g, min.vertices = 2)
  comp = length(gs)
  components = append(components, comp)

  for (t in 1:comp) {
    nodes = as.numeric(as_ids(V(gs[[t]])))
    if (line_r != 0) {
      for (x in 1:nrow(gendat)) {
        if (which(bus$type == 2)[x] %in% nodes) {
          ind = which(bus$type == 2)[x]
          bus[ind, 2] = 3
          break
        }
      }
    }
  }
  
  noslack = (1:nrow(bus))[bus$type != 3]
  
  colind = c(match(branch$fbus, bus$bus_id), match(branch$tbus, bus$bus_id))
  rowind = c(seq(1, nrow(branch)), seq(1, nrow(branch)))
  entries = c(rep(1, nrow(branch)), rep(-1, nrow(branch)))
  A = sparseMatrix(i = rowind, j = colind, x = entries, dims = c(nrow(branch), nrow(bus)))
  
  b = (branch$status) / (branch$x)
  
  Bf = sparseMatrix(i = rowind, j = colind, x = c(b, -b), dims = c(nrow(branch), nrow(bus)))
  Br = t(A) %*% Bf
  
  S = matrix(0, nrow(branch), nrow(bus))
  S[,noslack] = as.matrix(Bf[, noslack]) %*% inv(as.matrix(Br[noslack, noslack]))
  
  return(S)
}

#---------------------------------------------------------------------------------------------------------------------------------------------------
ptm <- proc.time()

branchdat = read.csv("~/branchdat.csv")
gendat = read.csv("~/gendat.csv")
busdat = read.csv("~/busdat.csv")

iterations = 1000

criticality = zeros(1, nrow(branchdat) + 1)
components = c()
data = get.data(busdat, gendat, branchdat)

for (line_removal in 0:0) { #nrow(branchdat)
  svMisc::progress(line_removal, nrow(branchdat))
  
  branchdat = read.csv("~/branchdat.csv")
  gendat = read.csv("~/gendat.csv")
  busdat = read.csv("~/busdat.csv")
  
  if (line_removal != 0) {
    branchdat[line_removal, 11] = 0
  }
  
  S = calculate.s(busdat, branchdat, line_removal)
  
  for (a in 1:iterations) {
    svMisc::progress(a, iterations)
    
    d_i = c()
    for (l in 1:nrow(busdat)) {
      d_i = append(d_i, runif(1, data$dmin[l], data$dmax[l]))
    }
    
    g_i = ((data$gmax - data$gmin) / sum(data$gmax - data$gmin) * (sum(d_i) - sum(data$gmin))) + data$gmin
    p = matrix(g_i) - matrix(d_i)
    
    f = S %*% p
    
    fail = 0
    success = 0
    for (l in 1:nrow(branchdat)) {
      range = c(0, data$flim[l] * 1.5)
      if (inside.range(abs(f[l]), range)) {
        success = success + 1
      }
      else {
        fail = fail + 1
      }
    }
    
    if (fail == 0) {
      criticality[line_removal + 1] = criticality[line_removal + 1] + 1
    }
  }
}
rm(a, l)
criticality = criticality / iterations
proc.time() - ptm
# criticality
# sort(criticality)
# rank(criticality)
