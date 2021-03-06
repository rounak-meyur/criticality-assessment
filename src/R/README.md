Criticality Assessment
================

## Initialization

First we initailize the file paths. As an example we show the analysis
on IEEE 39 bus power system. The system has 39 buses or nodes and 46
edges or branches. There are 10 generator buses in the system.

``` r
# Directory paths
path.dir <- dirname(getwd())
data.dir <- file.path(path.dir,"case")

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
    
    if(length(nodes)==1)
    {
      if(verbose) cat("Isolated node\n")
      bus$type[busind] <- 4
      bus$pd[busind] <- 0
    }
    else
    {
      if(3%in%btypes)
      {
        if(verbose) cat("Slack bus already present\n")
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
        bus$type[busind] <- 4
        bus$pd[busind] <- rep(0,length(nodes))
      }
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
    if((length(nodes)>1) & (max.demand<min.generation))
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
  noslack <- (1:nb)[(bus$type!=3)|(bus$type!=4)]
  
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
```

## Estimating probability of reaching a feasible case

The first aspect of criticality assessment arises regarding the
load-generation balance in every disconnected components in the network
(also termed as islands). We need to make sure that the load-generation
balance is maintained in each island. It is easy to show that the
feasible space shrinks after a contingency, i.e., the load-generation
balance constraint gets stringent after the loss of a branch. We aim to
evaluate this departure in feasibility through the computation of the
L-2 norm of the load-generation mismatch. First, we evaluate multiple
feasible scenarios for the base case and then we apply these scenarios
for a given contingency. Thereafter, we compute the load-generation
mismatch in each generated island. The departure in feasibility for each
scenario for a given contingency is evaluated as the L-2 norm of the
mismatch over all islands.

``` r
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
```

## Test Cases

First we evaluate the departure in feasibility for four test cases. We
remove branches 20, 27, 1 and 43 respecively in each case and compute
the departure in feasibility over 1000 scenarios.

![](README_files/figure-gfm/test.feasibility-1.png)<!-- -->

## Criticality evaluation

We evaluate the criticality for each component of the power grid
separately. For this purpose, first the coponents are identified and
then random test scenarios are generated for each component. The
criticality of each component is evaluated as the fraction of the total
number of random scenarios where it is feasible. This component
criticality is weighed by the size of the component and the weighted
average criticality of the entire power grid is computed. Note that for
isolated nodes and coponents which are infeasible, the criticality os
considered to be zero.

``` r
eval.component <- function(S,nodes,edgeind,bus,gen,branch,
                          base=100.0,max.iter=100,
                          verbose=FALSE)
{
  ##################################################
  # Evaluates the criticality of a component in the#
  # power grid network.                            #
  #                                                #
  # Inputs:                                        #
  #   S: PTDF matrix of the system                 #
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
  C <- as.matrix(C[,genind])
  
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


eval.contingency <- function(bus,branch,gen,verbose=F,
                          max.iter=100,base=100.0)
{
  ##################################################
  # Evaluates the criticality of a contingency in  #
  # power grid network.                            #
  #                                                #
  # Inputs:                                        #
  #   bus: csv table of bus data                   #
  #   gen: csv table of generator data             #
  #   branch: csv table of branch data             #
  #   base: base MVA: default to 100.0 MVA         #
  #   max.iter: number of random iterations        # 
  #   verbose: print outputs                       #
  #                                                #
  ##################################################
  
  # Initialize return output
  net.crit <- 0
  net.node <- 0
  
  # Handle islands
  dat <- handle.islands(bus,branch,gen,verbose=verbose)
  bus <- dat$bus
  gen <- dat$gen
  # Compute PTDF matrix
  S <- make.PTDF(bus,branch,gen)
  
  # Create a graph from the branch data
  g <- graph.data.frame(branch[, 1:2],directed=FALSE)
  # Remove the edge(s) whose status is/are equal to 0
  line_r <- which(branch$status==0)
  g <- delete_edges(g, line_r)
  # Decompose the graph into components and iterate
  gs <- decompose(g, min.vertices = 1)
  
  # Evaluate criticality for each island
  for(t in 1:length(gs))
  {
    nodes <- as.numeric(as_ids(V(gs[[t]])))
    nodeind <- match(nodes,bus$bus_id)
    edgeind <- (1:nrow(branch))[branch$status==1]
    nodetype <- bus$type[nodeind]
    
    if((length(nodes)>1)&(3 %in% nodetype))
    {
      crit <- eval.component(S,nodes,edgeind,bus,gen,branch,
                            max.iter=max.iter,base=base)
    }
    else crit <- 0
    
    cat("Island:",t,"-- Nodes:",length(nodes),
        " -- Criticality: Successes:",crit,"Failures:",100-crit,"\n")
    
    net.crit <- net.crit + crit*length(nodes)
    net.node <- net.node + length(nodes)
  }
  if (verbose) cat("Net criticality:",net.crit/net.node)
  return(net.crit/net.node)
}
```

## Examples of test cases

We now show examples of handling these islands by creating line trips
and settinf up new slack buses. We first divide the network into
disconnected components and search for slack bus. If a slack bus is
already present, we leave it as it is. If not, we search for a generator
bus and select the first among them as the slack bus for the island. If
there is no generator bus, we convert each of the buses to a slack bus.

### Test case 1: base case with no outages

``` r
# Basecase example
c <- eval.contingency(busdat,branchdat,gendat,
                   verbose=T)
```

    ## Component 1 : 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 19 20 21 22 23 25 26 28 29 39 30 18 11 31 32 24 27 33 34 35 36 37 38 
    ## Slack bus already present
    ## Slack bus: 31 Disconnected bus:  
    ## Island: 1 -- Nodes: 39  -- Criticality: Successes: 100 Failures: 0 
    ## Net criticality: 100

### Test case 2: outage in branch number 27

``` r
# Example of the function with branch 27 removed
alt1.branchdat <- branchdat
alt1.branchdat$status[27]=0

c <- eval.contingency(busdat,alt1.branchdat,gendat,
                   verbose=T)
```

    ## Component 1 : 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 21 22 23 25 26 28 29 39 30 18 11 31 32 24 27 35 36 37 38 
    ## Slack bus already present
    ## Slack bus: 31 Disconnected bus:  
    ## Component 2 : 19 20 33 34 
    ## At least one PV bus
    ## Slack bus: 33 Disconnected bus:  
    ## Island: 1 -- Nodes: 35  -- Criticality: Successes: 100 Failures: 0 
    ## Island: 2 -- Nodes: 4  -- Criticality: Successes: 100 Failures: 0 
    ## Net criticality: 100

### Test case 3: outage in branch number 20

``` r
# Example of the function with branch 20 removed
alt2.branchdat <- branchdat
alt2.branchdat$status[20]=0

c <- eval.contingency(busdat,alt2.branchdat,gendat,
                   verbose=T)
```

    ## Component 1 : 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 19 20 21 22 23 25 26 28 29 39 30 18 11 31 24 27 33 34 35 36 37 38 
    ## Slack bus already present
    ## Slack bus: 31 Disconnected bus:  
    ## Component 2 : 32 
    ## Isolated node
    ## Slack bus:  Disconnected bus: 32 
    ## Island: 1 -- Nodes: 38  -- Criticality: Successes: 100 Failures: 0 
    ## Island: 2 -- Nodes: 1  -- Criticality: Successes: 0 Failures: 100 
    ## Net criticality: 97.4359

### Test case 4: outage in branch number 20 and 27

``` r
# Example of the function with branches 20,27 removed
alt3.branchdat <- branchdat
alt3.branchdat$status[c(20,27)]=0

c <- eval.contingency(busdat,alt3.branchdat,gendat,
                   verbose=T)
```

    ## Component 1 : 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 21 22 23 25 26 28 29 39 30 18 11 31 24 27 35 36 37 38 
    ## Slack bus already present
    ## Slack bus: 31 Disconnected bus:  
    ## Component 2 : 19 20 33 34 
    ## At least one PV bus
    ## Slack bus: 33 Disconnected bus:  
    ## Component 3 : 32 
    ## Isolated node
    ## Slack bus:  Disconnected bus: 32 
    ## Island: 1 -- Nodes: 34  -- Criticality: Successes: 100 Failures: 0 
    ## Island: 2 -- Nodes: 4  -- Criticality: Successes: 100 Failures: 0 
    ## Island: 3 -- Nodes: 1  -- Criticality: Successes: 0 Failures: 100 
    ## Net criticality: 97.4359

### Test case 5: outage in branch number 44 and 45

We further alter the generator data such that the minimum generation of
Generator 9 at bus 38 is 300 MW. This results in an infeasible island
where the minimum generation exeeds the maximum demand. We check the
feasibility of our algorithm in this scenario.

``` r
# Example of the function with branches 44,45 removed
alt4.branchdat <- branchdat
alt4.branchdat$status[c(44,45)]=0
alt4.gendat <- gendat
alt4.gendat$pmin[9]=300.0

c <- eval.contingency(busdat,alt4.branchdat,alt4.gendat,
                   verbose=T)
```

    ## Component 1 : 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 19 20 21 22 23 25 26 28 39 30 18 11 31 32 24 27 33 34 35 36 37 
    ## Slack bus already present
    ## Slack bus: 31 Disconnected bus:  
    ## Component 2 : 29 38 
    ## At least one PV bus
    ## Slack bus: 38 Disconnected bus:  
    ## Maximum demand: 283.5  Minimum generation: 300
    ## Infeasible island. Zeroing demands and generations
    ## Island: 1 -- Nodes: 37  -- Criticality: Successes: 100 Failures: 0 
    ## Island: 2 -- Nodes: 2  -- Criticality: Successes: 0 Failures: 100 
    ## Net criticality: 94.87179
