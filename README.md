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

## Scenario generation

For a N-bus power system, we define a 2N-dimensional subspace of load
demands and generations at each of the N buses respectively. We sample
the demands and generations such that the total demand satisfies the
total generation in the network and the marginal distribution of each
demand/generation follows a uniform distribution. For this purpose, we
require to compute the limits of demands and generation at each bus in
the system. Thereafter, the demands/generations are randomly sampled to
generate test scenarios.

``` r
get.data <- function(bus,gen,branch,base=100.0)
{
  ##################################################
  # Computes the generation/demand limits at each  #
  # bus in the power system and returns as a list  #
  # object. It also returns the flow limit at each #
  # branch in the network.                         #
  #                                                #
  # Inputs:                                        #
  #   bus: csv table of bus data                   #
  #   gen: csv table of generator data             #
  #   branch: csv table of branch data             #
  #   base: base MVA, default value is 100.0 MVA   #
  #                                                #
  # Returns:                                       #
  #   data: a list object with gmin, gmax, dmin,   #
  #         dmax, flim representing minimum and    #
  #         maximum generation, demand and maximum #
  #         flow limits.                           #
  #                                                #
  ##################################################
  # Construct the gmin and gmax vectors
  busnum <- bus$bus_id
  genbus <- gen$bus_id
  nb <- nrow(bus)
  ng <- nrow(gen)
  gmax <- rep(0,nb)
  gmin <- rep(0,nb)
  for(i in 1:ng)
  {
    busind <- match(genbus[i],busnum)
    gmin[busind] = gmin[busind] + gen$pmin[i]/base
    gmax[busind] = gmax[busind] + gen$pmax[i]/base
  }
  
  # Construct the dmin and dmax vectors
  dmin <- rep(0,nb)
  dmax <- bus$pd/base
  
  # Construct the line rating vector
  flim <- branch$rateA/base
  
  # Return as an object
  return(list("gmin"=gmin,"gmax"=gmax,"dmax"=dmax,
              "dmin"=dmin,"flim"=flim))
}

# Get the limits
data <- get.data(busdat,gendat,branchdat)
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
handle.islands = function(bus,branch,verbose=FALSE)
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
    if(verbose) cat("Component",t,"\n")
    nodes <- as.numeric(as_ids(V(gs[[t]])))
    busind <- match(nodes,bus$bus_id)
    btypes <- bus$type[busind]
    
    if(!(3%in%btypes))
    {
      # No slack bus in component
      if(2%in%btypes)
      {
        # At least one PV bus in component
        slack <- which(btypes==2)[1]
      }
      else
      {
        # No PV bus in component
        slack <- which(btypes==1)
      }
      # Change bus type to slack
      bus$type[busind[slack]] <- 3
      
      if(verbose) cat("Following buses are converted to slack:",
          nodes[slack],"\n")
    }
    if(verbose) cat("Slack bus(es)",
        nodes[which(bus$type[busind]==3)],"\n")
  }
  bus
}
```

We now show examples of handling these islands by creating line trips
and settinf up new slack buses. We first divide the network into
disconnected components and search for slack bus. If a slack bus is
already present, we leave it as it is. If not, we search for a generator
bus and select the first among them as the slack bus for the island. If
there is no generator bus, we convert each of the buses to a slack bus.

    ## [1] "############## Example 1 ##############"

    ## Component 1 
    ## Slack bus(es) 31 
    ## Component 2 
    ## Following buses are converted to slack: 33 
    ## Slack bus(es) 33

    ## [1] "############## Example 2 ##############"

    ## Component 1 
    ## Slack bus(es) 31 
    ## Component 2 
    ## Following buses are converted to slack: 32 
    ## Slack bus(es) 32

    ## [1] "############## Example 3 ##############"

    ## Component 1 
    ## Slack bus(es) 31 
    ## Component 2 
    ## Following buses are converted to slack: 33 
    ## Slack bus(es) 33 
    ## Component 3 
    ## Following buses are converted to slack: 32 
    ## Slack bus(es) 32

## Power Transfer Distribution Factor (PTDF)

We compute the PTDF matrix for a given network. This matrix relates the
power injections at each node in the network to flows along edges of the
network.

``` r
make.PTDF = function(bus,branch)
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
  bus <- handle.islands(bus,branch)
  
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

# Compute PTDF matrix
S <- make.PTDF(busdat,branchdat)
```
