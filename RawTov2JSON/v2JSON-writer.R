#########
# Example: start from 4 source file: 1 -> connectivity matrix [NxN] (brain0.txt) 2,3 -> edge covariate matrices [NxN] (ecov_1.txt, ecov_2.txt)  4 -> l nodal attributes [Nxl]  
#########

# Change this!
setwd("/my/path/to/files")

corr_mat = as.matrix(read.table("data/brain0.txt"),header=FALSE)
are_there_attrs = TRUE
are_there_edgecov = TRUE

n_nodes = nrow(corr_mat)

#load data/attrs and use it to initialize nodes
nodes <- list()

if(are_there_attrs==TRUE){
attrs <- as.matrix(read.table("data/attrs.txt",header=TRUE,stringsAsFactors=FALSE))
for(r in 1:n_nodes) {
    nodeid <- paste(toString(r))
    
    #let's create "metadata" object with all attributes
    metadata = list()
    for(c in 1:ncol(attrs)) {
      attrname=colnames(attrs)[c]
      value=attrs[r,c]
      metadata[attrname] = value
    }

  #JGF node uses "metadata" to store key/value for each node.
  nodes[[nodeid]] = list(
    type="node type",
    label="node label",
    metadata=metadata
  )
}
} else {
  for(r in 1:n_nodes) {
    nodeid <- paste("node-", toString(r), sep="")
    nodes[[nodeid]] = list(
      type="node type",
      label="node label"
    )
  }
}


#also.. let's initialize edges between all row/col pairs (with empty metadata)
edges <- c()
for(r in 1:n_nodes) {
  source_nodeid <- paste(toString(r))
  for(c in 1:n_nodes) {
    target_nodeid <- paste(toString(c))
    edge <- list(
      source = source_nodeid,
      target = target_nodeid,
      directed = FALSE
    )
    print(length(edges))
    edges[[length(edges)+1]] <- edge
  }
}

if(are_there_edgecov==TRUE){
#now let's load the matrices
ecov1 = as.matrix(read.table("data/ecov_1.txt"),header=FALSE)
ecov2 = as.matrix(read.table("data/ecov_2.txt"),header=FALSE)

for(r in 1:n_nodes) {
  if(r > nrow(corr_mat)) {
    print("corr_mat contains less rows than the node count in attr")
    next
  }
  for(c in 1:n_nodes) {
    if(c > ncol(corr_mat)) {
      print("corr_mat contains less cols than the node count in attr")
      next
    }
    edges[[(r-1)*n_nodes+c]]$metadata <- list(
      corr_mat=corr_mat[[r,c]],
      ecov1=ecov1[[r,c]],
      ecov2=ecov2[[r,c]]
    )
  }
}
} else {
  for(r in 1:n_nodes) {
    for(c in 1:n_nodes) {
      edges[[(r-1)*n_nodes+c]]$metadata <- list(
        corr_mat=corr_mat[[r,c]]
      )
    }
  }
}

#put everything together
data <- list(
  graph = list( 
    directed = FALSE, 
    type = 'graph type', 
    label = 'some name', 
    metadata = list( somekey = 'somevalue' ),
    nodes = nodes,
    edges = edges
  )
)

#save to file
library(rjson)
jgf <- toJSON(data, indent=4)
cat(jgf,file="network.jgf")
