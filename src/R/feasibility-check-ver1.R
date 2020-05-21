# Import necessary libraries
library(Matrix)
library(spatstat.utils)
library(pracma)
library(igraph)


gen.bus.con <- function(bus,gen)
{
  ##################################################
  # Creates a sparse matrix of dimension nb-by-ng  #
  # which maps each generator index to index of the#
  # bus where it is connected                      #
  #                                                #
  # Inputs:                                        #
  #   bus: csv table of bus data                   #
  #   gen: csv table of generator data             #
  #                                                #
  # Returns:                                       #
  #   C: a sparse matrix with connection info      #
  #      given by zeros and ones                   #
  #                                                #
  ##################################################
  ng <- nrow(gen)
  nb <- nrow(bus)
  
  rowind <- match(gen$bus_id,bus$bus_id)
  colind <- 1:ng
  x <- rep(1,ng)
  
  C <- sparseMatrix(i=rowind,j=colind,
                    x=x,dims=c(nb,ng))
  return(C)
}


#####################################################
# Sampling strategies
sample.dg <- function(limits)
{
  ng <- length(limits$gmax)
  dmax <- limits$dmax; dmin <- limits$dmin
  gmax <- limits$gmax; gmin <- limits$gmin
  
  gi = c()
  for(i in 1:ng)
  {
    gi <- append(gi,runif(1,gmin[i],gmax[i]))
  }
  if(sum(gi)<=sum(dmax))
  {
    di <- dmax/sum(dmax)*sum(gi)
    return(list("success"=T,"di"=di,"gi"=gi))
  }
  else 
  {
    return(list("success"=F,"di"=NA,"gi"=NA))
  }
}

sample.g <- function(limits)
{
  ng <- length(limits$gmax)
  gmax <- limits$gmax; gmin <- limits$gmin
  
  gi <- c()
  for(i in 1:ng)
  {
    gi <- c(gi,runif(1,gmin[i],gmax[i]))
  }
  return(gi)
}

sample.d <- function(limits,gi)
{
  nd <- length(limits$dmax)
  dmax <- limits$dmax; dmin <- rep(0,nd)
  di <- rep(0,nd)
  
  # Generate uniform random variables between 0 and 1
  ui <- runif(nd, 0, 1)
  
  # Compute sum of randomized demands
  Sd <- sum(ui*dmax)
  Sg <- sum(gi)
  if(Sd>=Sg)
  {
    di <- ui*dmax*(Sg/Sd)
    return(list("success"=T,"di"=di,"gi"=gi))
  }
  else
  {
    delta <- 1.0; K <- 2
    while (delta*K > 1)
    {
      delta <- delta - 0.01
      # Define a threshold c for the random fractions
      c <- delta*(Sd/Sg)
      # Assign maximum demands to those fractions greater than c
      ind.greater <- which(ui >= c)
      ind.lesser <- which(ui < c)
      di[ind.greater] = dmax[ind.greater]
      # If sum of generation is less then the random fractions do not work
      if (Sg<sum(di))
      {
        return(list("success"=F,"di"=NA,"gi"=NA))
      }
      else
      {
        K = (Sd/Sg)*(Sg - sum(dmax[ind.greater]))/(sum(ui[ind.lesser]*dmax[ind.lesser]))
      }
    }
    # Assign demands to the remaining loads
    di[ind.lesser] = ui[ind.lesser]*dmax[ind.lesser]*K*Sg/Sd
    return(list("success"=T,"di"=di,"gi"=gi))
  }
}

###########################################################
# Feasibility evaluation

norm.vec <- function(x) sqrt(sum(x^2))


buildIslandMatrix <- function(bus,comps)
{
  ##################################################
  # Creates a sparse matrix of dimension nb-by-nc  #
  # which maps each node index to index of the     #
  # network component where it is present          #
  #                                                #
  # Inputs:                                        #
  #   bus: csv table of bus data                   #
  #   gen: list of network components              #
  #                                                #
  # Returns:                                       #
  #   M: a sparse matrix with connection info      #
  #      given by zeros and ones                   #
  #                                                #
  ##################################################
  nb <- nrow(bus)
  nc <- length(comps)
  
  rowind <- c()
  colind <- c()
  x <- rep(1,nb)
  
  for(i in 1:nc)
  {
    nodes <- as.numeric(as_ids(V(comps[[i]])))
    colind <- c(colind,match(nodes,bus$bus_id))
    rowind <- c(rowind,rep(i,length(nodes)))
  }
  M <- sparseMatrix(i=rowind,j=colind,
                    x=x,dims=c(nc,nb))
}

sample.base <- function(bus,gen,base=100.0,
                        sampling.method=2,
                        max.iter=1000,verbose=FALSE)
{
  ###################################################
  # Samples generation and demand for the base case #
  # where there is no contingency.                  #
  #                                                 #
  # Inputs:                                         #
  #   bus: csv table of bus data                    #
  #   gen: list of network components               #
  #   base: base MVA, default is 100.0 MVA          #
  #   sampling.method:                              #
  #     1 for randomized generation only            #
  #     2 for randomized generation and demands     #
  #   max.iter: total number of scenarios           #
  #   verbose: print outputs                        #
  #                                                 #
  # Returns:                                        #
  #   P:  a matrix with entries denoting the power  #
  #       injection at a bus in a scenario          #
  #                                                 #
  ###################################################
  
  # Get the sparse generator connection matrix
  C <- as.matrix(gen.bus.con(bus,gen))
  
  # Get the limits
  limits <- list("dmax"=bus$pd/base,
                 "dmin"=rep(0,nrow(bus)),
                 "gmax"=gen$pmax/base,
                 "gmin"=gen$pmin/base)
  
  # Sample the generations and demands
  P <- array(0,dim=c(nrow(bus),max.iter))
  iter <- 0
  while(iter<max.iter)
  {
    if (sampling.method==1)
    {
      sample <- sample.dg(limits)
    }
    else
    {
      gi <- sample.g(limits)
      sample <- sample.d(limits,gi)
    }
    if(sample$success==T)
    {
      if(verbose) cat("Sampling is successful. Iteration count:",
                      iter+1,"\n")
      di <- sample$di
      gi <- sample$gi
      iter <- iter+1
      
      # Compute power injection vector
      P[,iter] <- (C%*%gi)-di
    }
    else
    {
      if(verbose) cat("Sampling is unsuccessful\n")
    }
  }
  return(P)
}

eval.feasibility <- function(bus,gen,branch,brind,
                             base=100.0,max.iter=1000,
                             verbose=FALSE)
{
  ##################################################
  # Evaluates the criticality of a scenario in the #
  # power grid network.                            #
  #                                                #
  # Inputs:                                        #
  #   bus: csv table of bus data                   #
  #   gen: csv table of generator data             #
  #   branch: csv table of branch data             #
  #   brind: index of branch to be taken out       #
  #   base: base MVA: default to 100.0 MVA         #
  #   max.iter: number of random iterations        # 
  #   verbose: print outputs                       #
  #                                                #
  ##################################################
  
  # Create a graph from the branch data
  g <- graph.data.frame(branch[, 1:2],directed=FALSE)
  g <- delete_edges(g, brind)
  # Decompose the graph into components
  comps <- decompose(g, min.vertices = 1)
  M <- buildIslandMatrix(bus,comps)
  P <- sample.base(bus,gen,base=base,max.iter=max.iter,
                   verbose=verbose)
  dept.matrix <- as.matrix(M)%*%P
  dept <- apply(dept.matrix,2,norm.vec)
  return(dept)
}



###############################################################
# Directory paths
cur.dir <- setwd("~/GitHub/criticality-assessment/src/R")
path.dir <- dirname(cur.dir)
data.dir <- file.path(path.dir,"case")

# Load the csv data for the IEEE 39 bus system
branchdat <- read.csv(file.path(data.dir,"branchdat.csv"))
gendat <- read.csv(file.path(data.dir,"gendat.csv"))
busdat <- read.csv(file.path(data.dir,"busdat.csv"))


# Plot histogram of feasibility
par(mfrow=c(2,2))
dept <- eval.feasibility(busdat,gendat,branchdat,20)
hist(dept,xlab="Departure",ylab="Density",freq = FALSE,
     breaks=100,col='pink',
     main="Histogram of departure from feasibility with branch 20 removed")

dept <- eval.feasibility(busdat,gendat,branchdat,27)
hist(dept,xlab="Departure",ylab="Density",freq = FALSE,
     breaks=100,col='violet',
     main="Histogram of departure from feasibility with branch 27 removed")
dept <- eval.feasibility(busdat,gendat,branchdat,1)
hist(dept,xlab="Departure",ylab="Density",freq = FALSE,
     breaks=100,col='cyan',
     main="Histogram of departure from feasibility with branch 1 removed")

dept <- eval.feasibility(busdat,gendat,branchdat,43)
hist(dept,xlab="Departure",ylab="Density",freq = FALSE,
     breaks=100,col='green',
     main="Histogram of departure from feasibility with branch 43 removed")



###################################################################
I <- matrix(0,nrow=nrow(branchdat),ncol=2)
for(i in 1:nrow(branchdat))
{
  g <- graph.data.frame(branchdat[, 1:2],directed=FALSE)
  g <- delete_edges(g, i)
  # Decompose the graph into components
  comps <- decompose(g, min.vertices = 1)
  for(t in 1:length(comps))
  {
    nodes <- as.numeric(as_ids(V(comps[[t]])))
    I[i,t] <- length(nodes)
  }
}
