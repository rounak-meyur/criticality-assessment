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
  g_max = gen$pmax / 100.0
  g_min = gen$pmin / 100.0
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

identify.slack.node = function(bus, gen, nlist) {
  C = gen.bus.con(bus, gen)
  bus_ind = match(nlist, bus$bus_id)
  bus.pmax = as.matrix(C %*% gen$pmax)[bus_ind]
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
    bus_ind = match(nodes, bus$bus_id)
    btypes = bus$type[bus_ind]
    
    if(!(3 %in% btypes)) {
      if(length(nodes) == 1) {
        if(verbose) {
          cat("Isolated node\n")
        }
        bus$type[bus_ind[1]] = 3
        bus$pd[bus_ind] = 0
      } else if(2 %in% btypes) {
        if(verbose) {
          cat("At least one PV bus\n")
        }
        slack = identify.slack.node(bus, gen, nodes)
        bus$type[bus_ind[slack$max]] = 3
      } else {
        if(verbose) {
          cat("No PV bus\n")
        }
        bus$type[bus_ind[1]] = 3
        bus$pd[bus_ind] = rep(0,length(nodes))
      }
    } else {
      if(verbose) {
        cat("Slack bus already present\n")
      }
    }
    if(verbose) {
      cat("Slack bus:", nodes[which(bus$type[bus_ind] == 3)], "\n")
    }

    C = gen.bus.con(bus,gen)
    if(sum(bus$pd[bus_ind]) < sum(as.matrix(C%*%gen$pmin)[bus_ind])) {
      if(verbose) {
        cat("Maximum demand:", sum(bus$pd[bus_ind])," Minimum generation:", sum(as.matrix(C%*%gen$pmin)[bus_ind]))
        cat("\nInfeasible island. Zeroing demands and generations\n")
      }
      bus$pd[bus_ind] = rep(0, length(nodes))
      gen_ind = which(colSums(as.matrix(C)[bus_ind,]) == 1)
      gen$pmin[gen_ind] = rep(0, length(gen_ind))
      gen$pmax[gen_ind] = rep(0, length(gen_ind))
    }
  }
  return(list("bus" = bus, "gen" = gen, "graphs" = gs))
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
  
  return(list("gmin" = dat$gmin, "gmax" = dat$gmax, "dmax" = dat$dmax, "dmin" = dat$dmin, "flim" = dat$flim, "S" = S, "bus" = all$bus, "gen" = all$gen, "graphs" = all$graphs))
}

#---------------------------------------------------------------------------------------------------------------------------------------------------
ptm = proc.time()

branchdat = read.csv("~/branchdat.csv")
gendat = read.csv("~/gendat.csv")
busdat = read.csv("~/busdat.csv")

iterations = 1000

criticality = zeros(1, nrow(branchdat) + 1)
data = get.data(busdat, gendat, branchdat)

for (line_removal in 0:1) { #nrow(branchdat)
  svMisc::progress(line_removal, nrow(branchdat))
  
  branchdat = read.csv("~/branchdat.csv")
  gendat = read.csv("~/gendat.csv")
  busdat = read.csv("~/busdat.csv")
  
  if (line_removal != 0) {
    branchdat[line_removal, 11] = 0
  }
  
  S_data = calculate.s(busdat, branchdat, gendat, data)
  
  gs = S_data$graphs
  
  for (graph in gs) {
    nodes = as.numeric(as_ids(V(graph)))
    bus_ind = which(S_data$bus$bus_id %in% nodes)
    C = as.matrix(gen.bus.con(S_data$bus, S_data$gen))[bus_ind,]
    gen_ind = which(colSums(C) == 1)
    
    d_max = S_data$dmax[bus_ind]
    d_min = S_data$dmin[bus_ind]
    g_max = S_data$gmax[gen_ind]
    g_min = S_data$gmin[gen_ind]
    
    mat = as.matrix(as_edgelist(graph))
    mat = cbind(mat, 0)
    for (i in 1:nrow(mat)) {
      index = which(branchdat[which(branchdat$fbus == mat[i, 1]), 2] == mat[i, 2])[1]
      if (is.na(index)) {
        index = which(branchdat[which(branchdat$fbus == mat[i, 2]), 2] == mat[i, 1])[1]
        mat[i, 3] = which(branchdat$fbus == mat[i, 2])[index]
      } else {
        mat[i, 3] = which(branchdat$fbus == mat[i, 1])[index]
      }
    }
    for (i in 1:(nrow(mat)-1)) {
      if (as.numeric(mat[i, 3]) == as.numeric(mat[i+1,3])) {
        mat[i+1,3] = as.numeric(mat[i, 3]) + 1
      }
    }
    
    edge_ind = as.numeric(as.vector(mat[,3]))
    S.comp = S_data$S[edge_ind, bus_ind]
    flim = S_data$flim[edge_ind]
    
    for (a in 1:iterations) {
      svMisc::progress(a, iterations)
      
      d_i = c()
      for (l in 1:length(d_max)) {
        d_i = append(d_i, runif(1, d_min[l], d_max[l]))
      }
      
      if (sum(d_i) > sum(g_min)) {
        g_i = ((g_max - g_min) / sum(g_max - g_min) * (sum(d_i) - sum(g_min))) + g_min
        p = S %*% g_i - d_i
        f = S.comp %*% p
        
        fail = 0
        success = 0
        for (l in 1:length(edge_ind)) {
          range = c(0, flim[l] * 1.5)
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
      } else {
        a = a - 1
      }
    }
  }
}
print(proc.time() - ptm)
rm(a, l, i, index)
criticality = criticality / iterations
# criticality
# sort(criticality)
# rank(criticality)
