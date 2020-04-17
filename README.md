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
  ## INSERT CODE ##
  
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
