
## n:number of variables
## d:degree of nodes
randCG <- function(n, d) {
  order = sample(1:n, n, replace = FALSE)
  
  amat <- matrix(0, n, n)
  prob = d / (n - 1)
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      amat[i,j] <- Rlab::rbern(1,prob = prob)
    }
  }
  amat <- amat + t(amat)
  
  k <- sample(1:n, 1)
  if (k != 1) {
    chain_cut = cut(1:n, k, 1:k)
  } else {
    chain_cut = 1:n
  }
  for (i in 1:n) {
    for (j in 1:n) {
      if (as.integer(chain_cut[which(order == i)]) > as.integer(chain_cut[which(order == j)])) {
        amat[i,j] = 0
      }
    }
  }
  rownames(amat) = c(1:n)
  colnames(amat) = c(1:n)
  
  return(amat)
}



# is.chaingraph <- function(amat) {
#   wmat <- matrix(as.integer((amat + t(amat)) > 1), nrow = nrow(amat))
#   wg <- igraph::graph.adjacency(wmat, mode = "undirected")
#   cc <- igraph::clusters(wg)
#   neworder <- order(cc$membership)
#   a <- matrix(0, nrow = length(cc$csize), ncol = length(cc$csize))
#   b <- cumsum(cc$csize)
#   wmat <- amat[neworder, neworder]
#   for(i in 1: length(cc$csize)){
#     for(j in 1: length(cc$csize)){
#       if(j != i){
#         a[i,j] <- as.integer(sum(wmat[(max(b[i-1],0)+1):b[i],
#                                       (max(b[j-1],0)+1):b[j]]) > 0)
#       }
#     }
#   }
#   rownames(a) <- colnames(a) <- as.character(1:length(b))
#   output <- ggm::isAcyclic(a)
#   for(i in 1:length(b)){
#     temp <- wmat[(max(b[i-1],0)+1):b[i], (max(b[i-1],0)+1):b[i]]
#     if (!all(temp == t(temp))) {
#       output <- FALSE
#       break
#     }
#   }
#   chainorder <- ggm::topOrder(a)
#   vertorder<-c()
#   chainsize<-c()
#   if(output == TRUE){
#     for(k in 1:length(b)){
#       vertorder <- c(vertorder, which(cc$membership == chainorder[k]-1))
#       chainsize <- c(chainsize, cc$csize[chainorder[k]])
#     }
#   }
#   return(list(result = output,
#               vert.order = vertorder,
#               chain.size = chainsize))
# }
