Criticality Assessment
================

## Initialization

First we initailize the file paths. As an example we show the analysis
on IEEE 39 bus power system. The system has 39 buses or nodes and 46
edges or branches. There are 10 generator buses in the system.

``` r
# Directory paths
path.dir <- getwd()
data.dir <- file.path(path.dir,"src/case")

# Load the csv data for the IEEE 39 bus system
branchdat <- read.csv(file.path(data.dir,"branchdat39.csv"))
gendat <- read.csv(file.path(data.dir,"gendat39.csv"))
busdat <- read.csv(file.path(data.dir,"busdat39.csv"))
```

## Bus to generator connection matrix

``` r
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
```

## Handling islands in power system

When analyzing cascading failures due to multiple trips in the power
system, the following cases might occur:

  - Case 1: Isolated nodes in the network. In this case, we consider
    each of the isolated nodes as disconnected. For simplicity in
    further computation , we can assume them to be slack buses. The
    generations and demands at these nodes are also set to zero. This is
    done by setting the limits at the buses to zero.
  - Case 2: Disconnected components in the network. When dealing with
    such disconnected islands, we may face the following scenarios.
      - an island with no feasible demand/generation combination. This
        appears if there is no generation in the island or the total
        maximum demand is less than the total minimum generation.
      - an island with at least one feasible demand/generation
        combination.

<!-- end list -->

``` r
identify.slack <- function(bus,gen,nlist)
{
  # Get the bus-to-generator connection matrix
  C <- gen.bus.con(bus,gen)
  
  # Get the indices of buses in nodelist
  busind <- match(nlist,bus$bus_id)
  bus.pmax <- as.matrix(C%*%gen$pmax)[busind]
  
  # Returns the node index with maximum pmax
  return(which.max(bus.pmax))
}


handle.islands = function(bus,branch,gen,verbose=FALSE)
{
  ##################################################
  # Handles the islands in the network. Identifies #
  # new slack buses in each component.             #
  #                                                #
  # Inputs:                                        #
  #   bus: csv table of bus data                   #
  #   branch: csv table of branch data             #
  #   verbose: print outputs                       #
  #                                                #
  ##################################################
  
  # Create a graph from the branch data
  g <- graph.data.frame(branch[, 1:2], directed = FALSE)
  
  # Remove the edge(s) whose status is/are equal to 0
  line_r <- which(branch$status==0)
  g <- delete_edges(g, line_r)
  
  # Decompose the graph into components and iterate
  gs <- decompose(g, min.vertices = 1)
  for(t in 1:length(gs))
  {
    nodes <- as.numeric(as_ids(V(gs[[t]])))
    if(verbose) cat("Component",t,":",nodes,"\n")
    busind <- match(nodes,bus$bus_id)
    btypes <- bus$type[busind]
    
    if(!(3%in%btypes))
    {
      if(length(nodes)==1)
      {
        if(verbose) cat("Isolated node\n")
        bus$type[busind[1]] <- 4
        bus$pd[busind] <- 0
      }
      else if(2%in%btypes)
      {
        if(verbose) cat("At least one PV bus\n")
        slack <- identify.slack(bus,gen,nodes)
        bus$type[busind[slack]] <- 3
      }
      else
      {
        if(verbose) cat("No PV bus\n")
        bus$type[busind[1]] <- 4
        bus$pd[busind] <- rep(0,length(nodes))
      }
    }
    else
    {
      if(verbose) cat("Slack bus already present\n")
    }
    if(verbose) cat("Slack bus:",
        nodes[which(bus$type[busind]==3)],
        "Disconnected bus:",
        nodes[which(bus$type[busind]==4)],"\n")
    
    # Check for the other infeasible case: 
    # sum(dmax)<sum(gmin)
    max.demand <- sum(bus$pd[busind])
    C <- gen.bus.con(bus,gen)
    min.generation <- sum(as.matrix(C%*%gen$pmin)[busind])
    if(max.demand<min.generation)
    {
      if(verbose){
        cat("Maximum demand:",max.demand," Minimum generation:",min.generation)
        cat("\nInfeasible island. Zeroing demands and generations\n")
      }
      bus$type[busind] <- rep(4,length(nodes))
      bus$pd[busind] <- rep(0,length(nodes))
      genind <- which(colSums(as.matrix(C)[busind,])==1)
      gen$pmin[genind] <- rep(0,length(genind))
      gen$pmax[genind] <- rep(0,length(genind))
    }
  }
  return(list("bus"=bus,"gen"=gen))
}
```

We now show examples of handling these islands by creating line trips
and settinf up new slack buses. We first divide the network into
disconnected components and search for slack bus. If a slack bus is
already present, we leave it as it is. If not, we search for a generator
bus and select the first among them as the slack bus for the island. If
there is no generator bus, we convert each of the buses to a slack bus.

    ## [1] "#### Example 1: Branch 27 outage ####"

    ## Component 1 : 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 21 22 23 25 26 28 29 39 30 18 11 31 32 24 27 35 36 37 38 
    ## Slack bus already present
    ## Slack bus: 31 Disconnected bus:  
    ## Component 2 : 19 20 33 34 
    ## At least one PV bus
    ## Slack bus: 33 Disconnected bus:

    ## [1] "#### Example 2: Branch 20 outage ####"

    ## Component 1 : 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 19 20 21 22 23 25 26 28 29 39 30 18 11 31 24 27 33 34 35 36 37 38 
    ## Slack bus already present
    ## Slack bus: 31 Disconnected bus:  
    ## Component 2 : 32 
    ## Isolated node
    ## Slack bus:  Disconnected bus: 32

    ## [1] "#### Example 3: Branch 20,27 outage ####"

    ## Component 1 : 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 21 22 23 25 26 28 29 39 30 18 11 31 24 27 35 36 37 38 
    ## Slack bus already present
    ## Slack bus: 31 Disconnected bus:  
    ## Component 2 : 19 20 33 34 
    ## At least one PV bus
    ## Slack bus: 33 Disconnected bus:  
    ## Component 3 : 32 
    ## Isolated node
    ## Slack bus:  Disconnected bus: 32

    ## [1] "#### Example 4: Branch 44,45 outage, minimum generation of Gen 9: 300 MW ####"

    ## Component 1 : 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 19 20 21 22 23 25 26 28 39 30 18 11 31 32 24 27 33 34 35 36 37 
    ## Slack bus already present
    ## Slack bus: 31 Disconnected bus:  
    ## Component 2 : 29 38 
    ## At least one PV bus
    ## Slack bus: 38 Disconnected bus:  
    ## Maximum demand: 283.5  Minimum generation: 300
    ## Infeasible island. Zeroing demands and generations

## Power Transfer Distribution Factor (PTDF)

We compute the PTDF matrix for a given network. This matrix relates the
power injections at each node in the network to flows along edges of the
network.

``` r
make.PTDF = function(bus,branch,gen)
{
  ##################################################
  # Computes the PTDF matrix for a given network.  #
  #                                                #
  # Inputs:                                        #
  #   bus: csv table of bus data                   #
  #   branch: csv table of branch data             #
  #                                                #
  ##################################################
  nb <- nrow(bus)
  nl <- nrow(branch)
  busnum <- bus$bus_id
  
  # Handle islands without any type 3 bus
  bus <- handle.islands(bus,branch,gen)
  
  # The columns and rows to be chosen
  noslack <- (1:nb)[bus$type!=3]
  
  # Compute the PTDF matrix
  colind <- c(match(branch$fbus,busnum),
              match(branch$tbus,busnum))
  rowind <- c(seq(1,nl),seq(1,nl))
  entries <- c(rep(1,nl),rep(-1,nl))
  A <- sparseMatrix(i=rowind,j=colind,
                    x=entries,dims=c(nl,nb))
  
  
  b <- (branch$status)/(branch$x)
  
  # These are inv(X)*A and t(A)*inv(X)*A respectively
  Bf <- sparseMatrix(i=rowind,j=colind,
                     x=c(b,-b),dims=c(nl,nb))
  B <- t(A)%*%Bf
  
  # Create the PTDF matrix
  S <- matrix(0,nl,nb)
  S[,noslack] = as.matrix(Bf[,noslack]) %*% inv(as.matrix(B[noslack,noslack]))
  
  return(S)
}
```

## Sample space and scenario generation

For a N-bus power system with Ng generators, we define a
(N+Ng)-dimensional subspace of load demands and generations. We sample
the demands and generations such that the total demand satisfies the
total generation in the network and the marginal distribution of each
demand/generation follows a uniform distribution. For this purpose, we
require to compute the limits of demands and generation at each bus in
the system. Thereafter, the demands/generations are randomly sampled to
generate test scenarios.

``` r
sample.dg <- function(limits)
{
  ng <- length(limits$gmax)
  dmax <- limits$dmax; dmin <- limits$dmin
  gmax <- limits$dmax; gmin <- limits$dmin
  
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

eval.scenario <- function(nodes,edgeind,bus,gen,branch,
                          base=100.0,max.iter=100,
                          verbose=FALSE)
{
  ##################################################
  # Evaluates the criticality of a scenario for a  #
  # component in the power grid network.           #
  #                                                #
  # Inputs:                                        #
  #   nodes: list of nodes in component            #
  #   edgeind: list of branch indices in component #
  #   bus: csv table of bus data                   #
  #   gen: csv table of generator data             #
  #   branch: csv table of branch data             #
  #   base: base MVA: default to 100.0 MVA         #
  #   max.iter: number of random iterations        # 
  #   verbose: print outputs                       #
  #                                                #
  ##################################################
  
  # Get indices of buses and generators in the component
  busind <- match(nodes,bus$bus_id)
  C <- as.matrix(gen.bus.con(bus,gen))[busind,]
  genind <- which(colSums(C)==1)
  
  # Get the limits
  limits <- list("dmax"=bus$pd[busind]/base,
                 "dmin"=rep(0,length(busind)),
                 "gmax"=gen$pmax[genind]/base,
                 "gmin"=gen$pmin[genind]/base)
  
  # Get S matrix slice for the component
  S.comp <- S[edgeind,busind]
  flim <- branch$rateA[edgeind]
  
  # Perform the iterations
  criticality <- 0
  iter <- 0
  while(iter<max.iter)
  {
    sample <- sample.dg(limits)
    if(sample$success==T)
    {
      if(verbose) cat("Sampling is successful. Iteration count:",
                      iter+1,"\n")
      di <- sample$di
      gi <- sample$gi
      iter <- iter+1
      
      # Compute power injection and flows
      p.comp <- (C%*%gi)-di
      f.comp <- S.comp%*%p.comp
      
      # Success/failure of line limit
      fail <- 0
      success <- 0
      for(l in 1:length(edgeind))
      {
        range <- c(0,flim[l]*1)
        if(inside.range(abs(f.comp[l]),range)) success <- success + 1
        else fail <- fail + 1
      }
      
      # Compute criticality
      if(fail==0) criticality <- criticality + 1
      # Other methods
      # criticality <- criticality + success/(fail+success)
    }
    else
    {
      if(verbose) cat("Sampling is unsuccessful\n")
    }
  }
  return(criticality)
}
```

``` r
# Compute PTDF matrix
S <- make.PTDF(busdat,branchdat,gendat)

# nodes and edges
nodelist <- busdat$bus_id
edgeind <- 1:nrow(branchdat)

crit <- eval.scenario(nodelist,edgeind,
                      busdat,gendat,branchdat)

cat("Criticality: Successes:",crit,"Failures:",100-crit)
```

    ## Criticality: Successes: 100 Failures: 0
