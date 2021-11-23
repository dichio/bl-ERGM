#!/usr/bin/env Rscript
#### 1 INPUT

## 1.1 INPUT-EXT
library("rjson")
config <- fromJSON(file = "config.json")

namefile = as.character(config[1])
my_formula = as.character(config[2])
nsim_gof <- as.integer(config[3])
unfiltered = as.logical(config[4])

## 1.2 INPUT-INT

# 1.2.1 Network
#ECO Neural Network Filter----
ECO <- function(data, k=3, directed = FALSE)
{
  #correlation matrix
  if(nrow(data)==ncol(data)){cormat<-data
  }else{cormat<-cor(data)}
  C<-cormat
  n<-ncol(C)
  S<-C
  if(directed)
  {
    numcon<-k*n
    ind<-which(C!=0)
  }else{C<-upper.tri(C,diag=TRUE)
  numcon<-k/2*n
  ind<-which(upper.tri(C,diag=TRUE)!=0)}
  S<-ifelse(C==1,S,0)
  if(numcon>length(ind))
  {
    stop("Input matrix is too sparse")
  }
  sorind<-matrix(0,nrow=length(ind),ncol=2)
  G<-S
  S<-abs(S)
  x<-S[ind]
  y<-ind
  h<-cbind(ind,S[ind])
  sorind<-h[order(-h[,2]),]
  C[sorind[(numcon+1):nrow(sorind),1]]<-0
  if(directed)
  {W<-C}else{W<-C+t(C)
  diag(W)<-1}
  J<-G+t(G)
  diag(J)<-1
  W<-ifelse(W!=0,J,0)
  W<-as.data.frame(W)
  colnames(W)<-colnames(data)
  row.names(W)<-colnames(data)
  W<-as.matrix(W)
  return(W)
}
#----

json_data <- fromJSON(file=namefile)

n_nodes = length(json_data$graph$nodes)
n_attr = length(json_data$graph$nodes[[1]]$metadata)
n_ecov = length(json_data$graph$edges[[1]]$metadata) - 1

attrs = matrix(nrow = n_nodes, ncol = n_attr)
for(r in 1:n_nodes){
  for(c in 1:n_attr){
    attrs[r,c] = json_data$graph$nodes[[r]]$metadata[[c]]
  }
}

data = matrix(nrow = n_nodes, ncol = n_nodes)
ecov = vector(mode='list',length = n_ecov)
for(i in 1:length((ecov))){
  ecov[[i]] = empty_matrix<-matrix(NA,n_nodes,n_nodes) 
}

for(r in 1:n_nodes){
  for(c in 1:n_nodes){
    data[r,c] = json_data$graph$edges[[n_nodes*(r-1)+c]]$metadata[[1]]
    for(e in 1:n_ecov){
      ecov[[e]][r,c] = json_data$graph$edges[[n_nodes*(r-1)+c]]$metadata[[1+e]]    
    }
  }
}


if (unfiltered == TRUE){
  RB = ECO(data, 3, FALSE)
  data = ifelse(RB > 0, 1, 0)


library(network)
bnet = network(data, directed = FALSE)
}


#### 2 ERGM 
library(ergm)
my_log <- file("output/log-computation.txt")
sink(my_log, append = TRUE, type = "output")

## 2.2 ESTIMATION 
summary(bnet ~ edges + triangle + degree(1:10))
bfit <- eval(parse(text = paste("ergm(bnet ~", my_formula, ")")))
pdf("output/mcmc-diagnostic.pdf")
mcmc.diagnostics(bfit)
dev.off()

## 2.3 GOF 
bfit.gof <- gof(bfit,control=control.gof.formula(nsim = nsim_gof)) #output simulate 
bfit.gof

pdf("output/gof.pdf")
plot(bfit.gof)
dev.off()

#### 3 OUTPUT 

my_log_est <- file("output/estimation.txt")
sink(my_log_est, append = TRUE, type = "output")
print("Estimated values and covariance matrix:")
estimate = list(bfit$coefficients,bfit$covar)
print(estimate)
print("Summary of the fitting procedure:")
summary(bfit)


closeAllConnections()



