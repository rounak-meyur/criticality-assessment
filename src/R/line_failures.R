#install.packages("", repos="http://cran.r-project.org")
library(pracma)
library(igraph)
library(Matrix)
library(spatstat.utils)
library(matrixStats)
library(graphics)
library(doParallel)

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
  rowind = match(gen$bus_id, bus$bus_id)
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

handle.islands = function(bus, branch, gen, dat) {
  g = graph.data.frame(branch[, 1:2], directed = FALSE)
  g = delete_edges(g, which(branch$status == 0))
  gs = decompose(g, min.vertices = 1) #generate connected components
  
  for(t in 1:length(gs)) {
    nodes = as.numeric(as_ids(V(gs[[t]])))
    bus_ind = match(nodes, bus$bus_id)
    b_types = bus$type[bus_ind]
    
    if(!(3 %in% b_types)) {
      if(length(nodes) == 1) { #isolated node
        bus$type[bus_ind] = 3
        dat$dmax[bus_ind] = 0
      } else if(2 %in% b_types) { #one or more generators
        slack = identify.slack.node(bus, gen, nodes)
        bus$type[bus_ind[slack]] = 3
      } else { #no generators
        bus$type[bus_ind[1]] = 3
        dat$dmax[bus_ind] = 0
      }
    }
    
    C = gen.bus.con(bus, gen)
    if(sum(dat$dmax[bus_ind]) < sum(as.matrix(C %*% dat$gmin)[bus_ind])) {
      dat$dmax[bus_ind] = 0
      if (length(bus_ind) == 1) {
        if (bus_ind %in% gen$bus_id) {
          gen_ind = bus_ind
          dat$gmax[gen_ind] = 0
          dat$gmin[gen_ind] = 0
        }
      } else {
        gen_ind = which(colSums(as.matrix(C)[bus_ind,]) == 1)
        dat$gmax[gen_ind] = 0
        dat$gmin[gen_ind] = 0
      }
    }
  }
  return(list("bus" = bus, "gen" = gen, "graphs" = gs, "gmin" = dat$gmin, "gmax" = dat$gmax, "dmax" = dat$dmax, "dmin" = dat$dmin, "flim" = dat$flim))
}

calculate.s = function(bus, branch, gen, dat) {
  all_data = handle.islands(bus, branch, gen, dat)
  bus = all_data$bus
  gen = all_data$gen
  
  noslack = (1:nrow(bus))[(bus$type != 3)]
  
  colind = c(match(branch$fbus, bus$bus_id), match(branch$tbus, bus$bus_id))
  rowind = c(seq(1, nrow(branch)), seq(1, nrow(branch)))
  entries = c(rep(1, nrow(branch)), rep(-1, nrow(branch)))
  A = sparseMatrix(i = rowind, j = colind, x = entries, dims = c(nrow(branch), nrow(bus)))
  
  b = (branch$status) / (branch$x)
  
  Bf = sparseMatrix(i = rowind, j = colind, x = c(b, -b), dims = c(nrow(branch), nrow(bus)))
  Br = t(A) %*% Bf
  
  S = matrix(0, nrow(branch), nrow(bus))
  S[,noslack] = as.matrix(Bf[, noslack]) %*% inv(as.matrix(Br[noslack, noslack]))
  
  return(list("gmin" = all_data$gmin, "gmax" = all_data$gmax, "dmax" = all_data$dmax, "dmin" = all_data$dmin, "flim" = all_data$flim,
              "S" = S, "bus" = bus, "gen" = gen, "graphs" = all_data$graphs))
}

calculate.di = function(d_max, Sg, g_max) {
  u_i = c()
  for (l in 1:length(d_max)) {
    u_i = append(u_i, runif(1, 0, 1))
  }
  
  S_r = sum(u_i * d_max)
  
  if (S_r >= Sg) {
    d_i = c()
    for (l in 1:length(g_max)) {
      d_i = append(d_i, u_i[l] * d_max * Sg / S_r)
    }
  } else {
    delta = 1
    K = 1
    while ((delta == 1 | (delta * K > 1)) & delta >= 0) {
      delta = delta - 0.01
      c = delta * S_r / Sg
      
      d_i = rep(0, length(d_max))
      ind = which(u_i >= c)
      d_i[ind] = d_max[ind]
      
      if (Sg < sum(d_i)) {
        return(c())
      }
      
      K = (Sg - sum(d_max[which(u_i >= c)])) / (Sg / S_r * sum(u_i[which(u_i < c)] * d_max[which(u_i < c)]))
    }
    ind = which(u_i < c)
    d_i[ind] = u_i[ind] * d_max[ind] * K * Sg / S_r
  }
  return(d_i)
}

base_case = function(bus, gen, branch) {
  data = get.data(bus, gen, branch)
  S_data = calculate.s(bus, branch, gen, data)
  gs = S_data$graphs
  graph = gs[[1]]
  
  nodes = as.numeric(as_ids(V(graph)))
  total_nodes = length(nodes)
  bus_ind = which(S_data$bus$bus_id %in% nodes)
  C = as.matrix(gen.bus.con(S_data$bus, S_data$gen))[bus_ind,]
  gen_ind = which(colSums(C) == 1)
  # C = C[,gen_ind]
  
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
  for (i in 1:(nrow(mat) - 1)) {
    if (nrow(mat) > 1) {
      if (as.numeric(mat[i, 3]) == as.numeric(mat[i + 1,3])) {
        mat[i + 1,3] = as.numeric(mat[i, 3]) + 1
      }
    }
  }
  
  edge_ind = as.numeric(as.vector(mat[,3]))
  S.comp = S_data$S[edge_ind, bus_ind]
  flim = S_data$flim[edge_ind]
  return(list("C" = C, "S.comp" = S.comp, "flim" = flim, "edge_ind" = edge_ind, "bus_ind" = bus_ind, "gen_ind" = gen_ind))
}

calculate.crit = function(edge_ind, flim, f) {
  fail = 0
  success = 0
  crit = 0
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
    crit = 1
  }
  return(crit)
}

out.of.bounds = function(edge_ind, flim, f) {
  res = matrix(nrow = 3206, ncol = 1)
  for (l in 1:length(edge_ind)) {
    range = c(0, flim[l] * 1.5)
    if (inside.range(abs(f[l]), range)) {
      res[edge_ind[l], 1] = 0
    }
    else {
      res[edge_ind[l], 1] = 1
    }
  }
  return(res)
}

#---------------------------------------------------------------------------------------------------------------------------------------------------
branchdat = read.csv("~/branchdat.csv")
gendat = read.csv("~/gendat.csv")
busdat = read.csv("~/busdat.csv")

results = matrix(nrow = 3206, ncol = 1000)

base_case_list = base_case(busdat, gendat, branchdat)

iterations = 1000

criticality = rep(0.0, nrow(branchdat) + 1)
data = get.data(busdat, gendat, branchdat)

cl = makeCluster(28)
registerDoParallel(cl)

c = foreach (line_removal = 3007, .packages = c("igraph", "Matrix", "pracma", "spatstat.utils"))  %dopar% {
  svMisc::progress(line_removal, nrow(branchdat))
  
  branchdat = read.csv("~/branchdat.csv")
  gendat = read.csv("~/gendat.csv")
  busdat = read.csv("~/busdat.csv")
  
  branchdat[line_removal, 11] = 0
  S_data = calculate.s(busdat, branchdat, gendat, data)
  gs = S_data$graphs
  total_nodes = 0
  
  for (graph in gs) {
    nodes = as.numeric(as_ids(V(graph)))
    total_nodes = total_nodes + length(nodes)
    if (length(nodes) > 1) {
      bus_ind = which(S_data$bus$bus_id %in% nodes)
      C = as.matrix(gen.bus.con(S_data$bus, S_data$gen))[bus_ind,]
      gen_ind = which(colSums(C) == 1)
      if (!isempty(gen_ind)) {
        C = C[,gen_ind]
        
        d_max = S_data$dmax[bus_ind]
        d_min = S_data$dmin[bus_ind]
        g_max = S_data$gmax[gen_ind]
        g_min = S_data$gmin[gen_ind]
        
        mat = as.matrix(as_edgelist(graph))
        mat = cbind(mat, 0)
        ch = rep(0, 3206)
        for (i in 1:nrow(mat)) {
          u = 1
          ind = intersect(which(branchdat$fbus == mat[i, 1]), which(branchdat$tbus == mat[i, 2]))
          index = intersect(which(branchdat$fbus == mat[i, 1]), which(branchdat$tbus == mat[i, 2]))[u]
          if (is.na(index)) {
            u = 1
            index = intersect(which(branchdat$fbus == mat[i, 2]), which(branchdat$tbus == mat[i, 1]))[u]
            while (ch[index] != 0) {
              u = u + 1
              index = intersect(which(branchdat$fbus == mat[i, 2]), which(branchdat$tbus == mat[i, 1]))[u]
            }
            ch[index] = 1
            mat[i, 3] = index
          } else {
            if (ch[index] == 1 & length(ind) == 1) {
              u = 1
              index = intersect(which(branchdat$fbus == mat[i, 2]), which(branchdat$tbus == mat[i, 1]))[u]
              while (ch[index] != 0) {
                u = u + 1
                index = intersect(which(branchdat$fbus == mat[i, 2]), which(branchdat$tbus == mat[i, 1]))[u]
              }
              ch[index] = 1
              mat[i, 3] = index
            } else {
              while (ch[index] != 0) {
                u = u + 1
                index = intersect(which(branchdat$fbus == mat[i, 1]), which(branchdat$tbus == mat[i, 2]))[u]
              }
              ch[index] = 1
              mat[i, 3] = index
            }
          }
        }
        
        edge_ind = as.numeric(as.vector(mat[,3]))
        S.comp = S_data$S[edge_ind, bus_ind]
        flim = S_data$flim[edge_ind]
        
        crit = 0
        crit_base = 0
        bad = 0
        iter = 1
        
        while (iter <= iterations) {
          svMisc::progress(iter, iterations)
          
          C_base = base_case_list$C
          S.comp_base = base_case_list$S.comp
          flim_base = base_case_list$flim
          edge_ind_base = base_case_list$edge_ind
          bus_ind_base = base_case_list$bus_ind
          gen_ind_base = base_case_list$gen_ind
          
          g_i = c()
          for (l in 1:length(g_max)) {
            g_i = append(g_i, runif(1, g_min[l], g_max[l]))
          }
          Sg = sum(g_i)
          
          g_i_base = rep(0, length(gen_ind_base))
          g_i_base[gen_ind] = g_i
          
          d_i = c()
          while (isempty(d_i)) {
            d_i = calculate.di(d_max, Sg, g_max)
          }
          
          d_i_base = rep(0, length(bus_ind_base))
          d_i_base[bus_ind] = d_i
          
          p = (C %*% g_i) - d_i
          p_base = (C_base %*% g_i_base) - d_i_base
          
          f = S.comp %*% p
          f_base = S.comp_base %*% p_base
          
          base = calculate.crit(edge_ind_base, flim_base, f_base)
          if (base == 1) {
            iter = iter + 1
            crit_base = crit_base + 1
            crit = crit + calculate.crit(edge_ind, flim, f)
            results[,crit_base] = out.of.bounds(edge_ind, flim, f)
          }
        }
        criticality[line_removal + 1] = criticality[line_removal + 1] + (length(nodes) * crit / crit_base)
      }
    }
  }
  criticality[line_removal + 1] = criticality[line_removal + 1] / total_nodes
  criticality[line_removal + 1] = 1 - criticality[line_removal + 1]
  
}

for (i in 1:(length(criticality) - 1)) {
  criticality[i + 1] = c[[i]]
}
#---------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(results_3007, "results_3007.csv")
