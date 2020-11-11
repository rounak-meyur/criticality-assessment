library(Matrix)
library(pracma)

###############################################################
# Directory paths
cur.dir <- setwd("C:/Users/rouna/Documents/GitHub/criticality-assessment/src/R")
path.dir <- dirname(cur.dir)
data.dir <- file.path(path.dir,"case")
result.dir <- file.path(path.dir,"results")

# Load the csv data for the IEEE 39 bus system
branchdat <- read.csv(file.path(data.dir,"branchdat.csv"))
gendat <- read.csv(file.path(data.dir,"gendat.csv"))
busdat <- read.csv(file.path(data.dir,"busdat.csv"))




# Functions to compute PTDF and LODF
make.LODF = function(bus,branch)
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
  #bus <- handle.islands(bus,branch,gen)
  
  # The columns and rows to be chosen
  noslack <- (1:nb)[(bus$type==1)|(bus$type==2)]
  
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
  PTDF <- matrix(0,nl,nb)
  PTDF[,noslack] = as.matrix(Bf[,noslack]) %*% inv(as.matrix(B[noslack,noslack]))
  
  H <- PTDF %*% t(A)
  h = diag(H)
  C = matrix(1,nl,nl) - (matrix(1,nl,1) %*% t(h))
  LODF = H / C
  LODF = LODF - diag(diag(LODF)) - diag(nl)
  
  return(LODF)
}


eval.overloads <- function(i)
{
  file.name <- paste("flows_1.5_",toString(i),".csv", sep = "")
  flowdat <- read.csv(file.path(result.dir,file.name))
  flow.mean <- apply(flowdat,1,mean)
  flow.rate <- 1.5*branchdat$rateA/100.0
  # flow.buffer <- flow.rate-flow.mean
  L <- make.LODF(busdat,branchdat)
  LODF.interest <- L[,i]
  print(sort.int(LODF.interest,index.return=TRUE,
                 decreasing=TRUE)$ix[1:3])
  # flow.cont <- flow.mean[i]
  # 
  # rhs <- abs(flow.cont)*abs(LODF.interest)
  lhs <- flow.mean + LODF.interest*flow.mean[i]
  print(sort.int(abs(lhs)/flow.rate,index.return=TRUE,
                 decreasing=TRUE)$ix[1:3])
  # check <- flow.buffer-rhs
  check <- flow.rate - abs(lhs)
  
  # plot(check,type='h',xlab="Transmission lines in the network",
  #      ylab="Buffer - Redistribution")
  
  overloaded <- (1:nrow(branchdat))[check<0.0]
  return(overloaded)
}

over.85 <- eval.overloads(85)
print(over.85)
over.87 <- eval.overloads(87)
print(over.87)
over.464 <- eval.overloads(464)
print(over.464)
over.3007 <- eval.overloads(3007)
print(over.3007)



# # Check individually
# ind = 3007
# flowdat <- read.csv(file.path(result.dir,
#                               paste("flows_1.5_",toString(ind),
#                                     ".csv",sep="")))
# flow.check <- as.matrix(flowdat[ind,])
# L <- make.LODF(busdat,branchdat)
# LODF.check <- as.matrix(L[,ind])
# redistributed.flow <- LODF.check%*%flow.check
# post.flow.check <- flowdat+redistributed.flow
# flow.rate <- repcol(1.5*branchdat$rateA/100.0,1000)
# flag.check <- 