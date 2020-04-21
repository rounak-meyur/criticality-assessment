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
library(matrixStats)

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

gen.bus.con = function(bus, gen) {
  rowind = match(gen$bus_id,bus$bus_id)
  colind = 1:nrow(gen)
  x = rep(1,nrow(gen))
  
  C = sparseMatrix(i = rowind, j = colind, x = x, dims = c(nrow(bus), nrow(gen)))
  return(C)
}

identify.slack = function(bus, gen, nlist) {
  C = gen.bus.con(bus, gen)
  
  busind = match(nlist, bus$bus_id)
  bus.pmax = as.matrix(C %*% gen$pmax)[busind]
  
  return(which.max(bus.pmax))
}

handle.islands = function(bus, branch, gen, verbose = FALSE) {
  g = graph.data.frame(branch[, 1:2], directed = FALSE)
  
  line_r = which(branch$status == 0)
  g = delete_edges(g, line_r)
  
  gs = decompose(g, min.vertices = 1)
  for(t in 1:length(gs))
  {
    nodes = as.numeric(as_ids(V(gs[[t]])))
    if(verbose) {
      cat("Component",t,":",nodes,"\n")
    }
    busind = match(nodes, bus$bus_id)
    btypes = bus$type[busind]
    
    if(!(3 %in% btypes)) {
      if(length(nodes) == 1) {
        if(verbose) {
          cat("Isolated node\n")
        }
        bus$type[busind[1]] = 3
        bus$pd[busind] = 0
      } else if(2 %in% btypes) {
        if(verbose) {
          cat("At least one PV bus\n")
        }
        slack = identify.slack(bus, gen, nodes)
        bus$type[busind[slack$max]] = 3
      } else {
        if(verbose) {
          cat("No PV bus\n")
        }
        bus$type[busind[1]] = 3
        bus$pd[busind] = rep(0,length(nodes))
      }
    } else {
      if(verbose) {
        cat("Slack bus already present\n")
      }
    }
    if(verbose) {
      cat("Slack bus:", nodes[which(bus$type[busind] == 3)], "\n")
    }

    C = gen.bus.con(bus,gen)
    if(sum(bus$pd[busind]) < sum(as.matrix(C%*%gen$pmin)[busind])) {
      if(verbose) {
        cat("Maximum demand:", sum(bus$pd[busind])," Minimum generation:", sum(as.matrix(C%*%gen$pmin)[busind]))
        cat("\nInfeasible island. Zeroing demands and generations\n")
      }
      bus$pd[busind] = rep(0, length(nodes))
      genind = which(colSums(as.matrix(C)[busind,]) == 1)
      gen$pmin[genind] = rep(0, length(genind))
      gen$pmax[genind] = rep(0, length(genind))
    }
  }
  return(list("bus" = bus, "gen" = gen))
}

calculate.s = function(bus, branch, gen, dat) {
  all = handle.islands(bus,branch,gen)
  bus = all$bus

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
  
  return(list("gmin" = dat$gmin, "gmax" = dat$gmax, "dmax" = dat$dmax, "dmin" = dat$dmin, "flim" = dat$flim, "S" = S, "bus" = all$bus, "gen" = all$gen))
}

#---------------------------------------------------------------------------------------------------------------------------------------------------
ptm = proc.time()

branchdat = read.csv("~/branchdat.csv")
gendat = read.csv("~/gendat.csv")
busdat = read.csv("~/busdat.csv")

iterations = 1000

criticality = zeros(1, nrow(branchdat) + 1)
data = get.data(busdat, gendat, branchdat)

for (line_removal in 1:2) { #nrow(branchdat)
  svMisc::progress(line_removal, nrow(branchdat))
  
  branchdat = read.csv("~/branchdat.csv")
  gendat = read.csv("~/gendat.csv")
  busdat = read.csv("~/busdat.csv")
  
  if (line_removal != 0) {
    branchdat[line_removal, 11] = 0
  }
  
  S_data = calculate.s(busdat, branchdat, gendat, data)
  
  ind = c()
  for (h in 1:nrow(S_data$gen)) {
    ind = append(ind, which(S_data$bus$bus_id == S_data$gen$bus_id[h]))
  }
  
  C = gen.bus.con(S_data$bus, S_data$gen)
  
  for (a in 1:iterations) {
    svMisc::progress(a, iterations)
    
    d_i = c()
    for (l in 1:nrow(S_data$bus)) {
      d_i = append(d_i, runif(1, S_data$dmin[l], S_data$dmax[l]))
    }
    
    g_i = ((S_data$gmax - S_data$gmin) / sum(S_data$gmax - S_data$gmin) * (sum(d_i) - sum(S_data$gmin))) + S_data$gmin
    
    if (sum(d_i) > sum(g_i)) {
      g_i = ((S_data$gmax - S_data$gmin) / sum(S_data$gmax - S_data$gmin) * (sum(d_i) - sum(S_data$gmin))) + S_data$gmin
    }
    
    g_i = g_i[c(ind)]
    
    node_gen = C %*% g_i
    p = node_gen - d_i

    f = S_data$S %*% p
    
    fail = 0
    success = 0
    for (l in 1:nrow(branchdat)) {
      range = c(0, S_data$flim[l] * 1.5)
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
  print(proc.time() - ptm)
}
rm(a, l)
criticality = criticality / iterations
# criticality
# sort(criticality)
# rank(criticality)
