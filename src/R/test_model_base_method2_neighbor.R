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

  d_max = S_data$dmax[bus_ind]
  d_min = S_data$dmin[bus_ind]
  g_max = S_data$gmax[gen_ind]
  g_min = S_data$gmin[gen_ind]
  
  mat = as.matrix(as_edgelist(graph))
  mat = cbind(mat, 0)
  for (i in 1:nrow(mat)) {
    index = which(branch[which(branch$fbus == mat[i, 1]), 2] == mat[i, 2])[1]
    if (is.na(index)) {
      index = which(branch[which(branch$fbus == mat[i, 2]), 2] == mat[i, 1])[1]
      mat[i, 3] = which(branch$fbus == mat[i, 2])[index]
    } else {
      mat[i, 3] = which(branch$fbus == mat[i, 1])[index]
    }
  }
  for (i in 1:(nrow(mat) - 1)) {
    if (nrow(mat) > 1) {
      if (as.numeric(mat[i, 3]) == as.numeric(mat[i + 1,3])) {
        mat[i + 1, 3] = as.numeric(mat[i, 3]) + 1
      }
    }
  }
  
  edge_ind = as.numeric(as.vector(mat[,3]))
  S.comp = S_data$S[edge_ind, bus_ind]
  flim = S_data$flim[edge_ind]
  return(list("C" = C, "S.comp" = S.comp, "flim" = flim, "edge_ind" = edge_ind, "bus_ind" = bus_ind, "gen_ind" = gen_ind))
}

calculate.crit = function(edge_ind, flim, f, mul) {
  fail = 0
  success = 0
  crit = 0
  for (l in 1:length(edge_ind)) {
    range = c(0, flim[l] * mul)
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

neighbor_check = function(branch, mul, st) {
  check.criticality = read.delim2(st, header = FALSE, stringsAsFactors = FALSE)
  line_id = c()
  for (line in 1:nrow(check.criticality)) {
    rem = check.criticality[line, 2]
    rem = strsplit(rem, ",")
    rem = rem[[1]]
    line_id = append(line_id, rem)
  }
  line_id = unique(line_id)
  for (l in 1:length(line_id)) {
    branch[as.integer(line_id[l]), 6] = branch[as.integer(line_id[l]), 6] * mul
  }
  return(branch)
}

load_check = function(bus, mul, st) {
  check.criticality = read.delim2(st, header = FALSE, stringsAsFactors = FALSE)
  line_id = c()
  for (line in 1:nrow(check.criticality)) {
    rem = check.criticality[line, 2]
    rem = strsplit(rem, ",")
    rem = rem[[1]]
    bus_id = append(bus_id, rem)
  }
  bus_id = unique(bus_id)
  for (l in 1:length(bus_id)) {
    bus[as.integer(bus_id[l]), 3] = bus[as.integer(bus_id[l]), 3] * mul
  }
  return(bus)
}

full_model = function(rate, rateA_all, neighbor, neighbor_file, load, load_file) {
  criticality_final = zeros(4, 3207)
  
  for (v in rate) {

    branchdat = read.csv("~/branchdat.csv")
    gendat = read.csv("~/gendat.csv")
    busdat = read.csv("~/busdat.csv")
    
    if (neighbor) {
      branchdat = neighbor_check(branchdat, v, neighbor_file)
    }
    
    if (load) {
      busdat = load_check(busdat, v, load_file)
    }
    
    base_case_list = base_case(busdat, gendat, branchdat)
    
    iterations = 1000
    
    data = get.data(busdat, gendat, branchdat)
    
    require("doParallel")
    require("foreach")
    
    cl = makeCluster(28)
    registerDoParallel(cl)
    
    print(v)
    
    ptime = system.time({
      c = foreach (line_removal = 1:nrow(branchdat), .packages = c("igraph", "Matrix", "pracma", "spatstat.utils"), .combine = 'c')  %dopar% {

        branchdat = read.csv("~/branchdat.csv")
        gendat = read.csv("~/gendat.csv")
        busdat = read.csv("~/busdat.csv")
        
        if (neighbor) {
          branchdat = neighbor_check(branchdat, v, neighbor_file)
        }
        
        if (load) {
          busdat = load_check(busdat, v, load_file)
        }
        
        criticality = 0
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
                  if (as.numeric(mat[i, 3]) == as.numeric(mat[i + 1, 3])) {
                    mat[i + 1, 3] = as.numeric(mat[i, 3]) + 1
                  }
                }
              }
              
              edge_ind = as.numeric(as.vector(mat[,3]))
              S.comp = S_data$S[edge_ind, bus_ind]
              flim = S_data$flim[edge_ind]
              
              crit = 0
              bad = 0
              iter = 1
              
              while (iter <= iterations) {

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
                
                base = calculate.crit(edge_ind_base, flim_base, f_base, rateA_all)
                if (base == 1) {
                  iter = iter + 1
                  crit = crit + calculate.crit(edge_ind, flim, f, rateA_all)
                }
              }
              criticality = criticality + (length(nodes) * crit / iterations)
            }
          }
        }
        criticality = criticality / total_nodes
        criticality = 1 - criticality
        criticality
      }})
    
    stopCluster(cl)
    
    for (i in 1:(nrow(branchdat))) {
      criticality_final[which(rate == v), i + 1] = c[[i]]
    }
  }
  
  return(criticality_final)
}

#---------------------------------------------------------------------------------------------------------------------------------------------------
criticality1.5N = full_model(c(3, 4, 5, 6), 1.5, TRUE, "~/check-criticality-1.5.txt", FALSE, "")
criticality1.25N = full_model(c(3, 4, 5, 6), 1.25, TRUE, "~/check-criticality-1.25.txt", FALSE, "")
criticality1.5L = full_model(c(0.9, 0.8, 0.7, 0.6), 1.5, FALSE, "", TRUE, "~/check-load-criticality-1.5.txt")
criticality1.25L = full_model(c(0.9, 0.8, 0.7, 0.6), 1.25, FALSE, "", TRUE, "~/check-load-criticality-1.25.txt")

write.csv(criticality1.5N, "criticality_base_1.5_2_neighbors.csv")
write.csv(criticality1.25N, "criticality_base_1.25_2_neighbors.csv")
write.csv(criticality1.5L, "criticality_base_1.5_2_loads.csv")
write.csv(criticality1.25L, "criticality_base_1.25_2_loads.csv")

#---------------------------------------------------------------------------------------------------------------------------------------------------
h = hist(criticality1.5N[1,], breaks = 20, main = "Transmission Line Criticality Frequencies", xlab = "Criticality", ylab = "Count", 
         col = "darkmagenta", freq = TRUE)
text(h$mids, h$counts, labels = h$counts, adj = c(0.5, -0.5))

#---------------------------------------------------------------------------------------------------------------------------------------------------
# criticality_base_1.5_2 = read.csv("~/criticality_base_1.5_2.csv")
# criticality_base_1.5_2_neighbors = read.csv("~/criticality_base_1.5_2_neighbors.csv")
# e = read.delim2("~/check-criticality-1.5.txt", header = FALSE, stringsAsFactors = FALSE)
# x = e[,1]
# x = x + 1
# plot(criticality_base_1.5_2[x, 2], ylim = c(0, 1), col = "black")
# points(criticality_base_1.5_2_neighbors[2, x], col = "red")
# points(criticality_base_1.5_2_neighbors[3, x], col = "purple")
# points(criticality_base_1.5_2_neighbors[4, x], col = "blue")
# points(criticality_base_1.5_2_neighbors[5, x], col = "green")

criticality_base_1.25_2 = read.csv("~/criticality_base_1.25_2.csv")
criticality_base_1.25_2_neighbors = read.csv("~/criticality_base_1.25_2_neighbors.csv")
e = read.delim2("~/check-criticality-1.25.txt", header = FALSE, stringsAsFactors = FALSE)
x = e[,1]
# x = x + 1
plot(criticality_base_1.25_2[x, 2], ylim = c(0, 1), col = "black")
points(criticality_base_1.25_2_neighbors[2, x], col = "red")
points(criticality_base_1.25_2_neighbors[3, x], col = "purple")
points(criticality_base_1.25_2_neighbors[4, x], col = "blue")
points(criticality_base_1.25_2_neighbors[5, x], col = "green")

# criticality_base_1.5_2 = read.csv("~/criticality_base_1.5_2.csv")
# criticality_base_1.5_2_loads = read.csv("~/criticality_base_1.5_2_loads.csv")
# e = read.delim2("~/check-criticality-1.5.txt", header = FALSE, stringsAsFactors = FALSE)
# x = e[,1]
# # x = x + 1
# plot(criticality_base_1.5_2[x, 2], ylim = c(0, 1), col = "black")
# points(criticality_base_1.5_2_loads[2, x], col = "red")
# points(criticality_base_1.5_2_loads[3, x], col = "purple")
# points(criticality_base_1.5_2_loads[4, x], col = "blue")
# points(criticality_base_1.5_2_loads[5, x], col = "green")

# criticality_base_1.25_2 = read.csv("~/criticality_base_1.5_2.csv")
# criticality_base_1.25_2_loads = read.csv("~/criticality_base_1.5_2_loads.csv")
# e = read.delim2("~/check-criticality-1.25.txt", header = FALSE, stringsAsFactors = FALSE)
# x = e[,1]
# x = x + 1
# plot(criticality_base_1.25_2[x, 2], ylim = c(0, 1), col = "black")
# points(criticality_base_1.25_2_loads[2, x], col = "red")
# points(criticality_base_1.25_2_loads[3, x], col = "purple")
# points(criticality_base_1.25_2_loads[4, x], col = "blue")
# points(criticality_base_1.25_2_loads[5, x], col = "green")