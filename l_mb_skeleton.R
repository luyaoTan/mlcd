  ############################# An LWF CG from the LCD package ###################################
  toy.graph <- matrix(c(0,1,1,0,0,0,0,0,0,0,0,
                        1,0,0,1,0,0,0,0,0,0,0,
                        1,0,0,0,1,0,0,0,0,0,0,
                        0,1,0,0,0,1,1,0,0,0,0,
                        0,0,0,0,0,1,0,0,1,0,0,
                        0,0,0,0,1,0,0,0,0,0,1,
                        0,0,0,0,0,0,0,1,0,0,1,
                        0,0,0,0,0,0,1,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,1,0,
                        0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0),
                      nrow = 11, byrow = TRUE)
  rownames(toy.graph) <- c("a","b","c","d","e","f","g","h","i","j","k")
  colnames(toy.graph) <- c("a","b","c","d","e","f","g","h","i","j","k")
  #############################################################################################################
  #############################################################################################################
  ############################# Plot LWF CGs ##################################################################
  #############################################################################################################
  #############################################################################################################
  # Input
  # -amat: adjacency matrix of the LWF CG
  
  ## ----------------------------------------------------------------------
  ## Author: ***************** August, 2019 (adapted from the "lcd")
  ################################################################################
  ################################################################################
  ################################################################################
  # library(igraph)
  `plotCG` <- function(amat, plain = TRUE){
    `skeletonCG` <- function(amat)
    {
      0 + (amat + t(amat) > 0)
    }
    e.total <- length(which(skeletonCG(amat) == 1))/2
    if(e.total == 0)
      return("The graph has no edge!")
    ## el <- matrix("", nrow = e.total, ncol = 2)
    el <- matrix(NA, nrow = e.total, ncol = 2)
    vset <- sort(rownames(amat))
    amat <- amat[vset,vset]
    p <- nrow(amat)
    arrow <- c()
    k <- 1
    for(i in 1:(p-1))
      for(j in (i+1):p){
        if(amat[i,j] == 1){
          ## el[k,1] <- i-1
          el[k,1] <- i
          ## el[k,2] <- j-1
          el[k,2] <- j
          k <- k + 1
          if(amat[j,i] == 0){
            arrow <- c(arrow, ">")
          } else {
            arrow <- c(arrow, "-")
          }
        } else {
          if(amat[j,i] == 1){
            ## el[k,1] <- j-1
            el[k,1] <- j
            ## el[k,2] <- i-1
            el[k,2] <- i
            k <- k + 1
            arrow <- c(arrow, ">")
          }
        }
      }
    g <- graph.empty()
    g <- add.vertices(g, p)
    V(g)$name <- vset
    g <- add.edges(g, t(el))
    if(plain){
      plot.igraph(g, #layout = layout.reingold.tilford,
                  edge.arrow.mode = arrow,
                  vertex.size = min(15, 500/p),
                  edge.arrow.size = min(1, 30/p))
    } else{
      tkplot(g, layout = layout.reingold.tilford,
             edge.arrow.mode = arrow,
             vertex.size = min(15, 500/p),
             edge.arrow.size = min(1, 30/p))
    }
    V(g)
  }
  
  
  
  
  
  
  #############################################################################################################
  #############################################################################################################
  ############################# MBC-CSP Algorithm for LWF CGs ################################################
  #############################################################################################################
  #############################################################################################################
  # -dataset: a data frame containing the variables in the model.
  
  # -cluster: an optional cluster object from package parallel.
  
  # -test: a character string, the label of the conditional independence test to be used in the
  # algorithm. If none is specified, the default test statistic is the mutual information
  # for categorical variables, the Jonckheere-Terpstra test for ordered factors and
  # the linear correlation for continuous variables.
  
  # -alpha: a numeric value, the target nominal type I error rate.
  
  # -max.sx: a positive integer, the maximum allowed size of the conditioning sets used in
  # conditional independence tests. The default is that there is no limit on size 
  
  ## ----------------------------------------------------------------------
  ## Author: ***************** August, 2019 (adapted from the "bnlearn" and "lcd")
  ################################################################################
  ################################################################################
  ################################################################################
  learn.skeleton <- function(dataset,alpha = 0.05 , test = "zf", cluster = NULL,
                      max.sx = ncol(dataset), debug = FALSE){
    data.info = bnlearn:::check.data(dataset, allow.missing = TRUE)
    complete=data.info$complete.nodes
    corr = cor(dataset)
    #####################################################
    ##################### Computing Adjs & sepset  ######
    #####################################################
    
    ##  Input: T: the target variable; alpha: significance level
    ##  Output: adj(T): the set of variables adjacent to T; sepset: separation set
    
    ## how to build a sepset with names
    ## note that dataset should be in the data frame format
    
    `adjecencies` <- function(dataset, alpha = alpha){
      #名为sepset的列表，每个元素都是一个列表，每个列表都有数据集的列名作为其内部列表的命名
      sepset <- lapply(colnames(dataset), function(.) {vec<-vector("list",ncol(dataset))
      names(vec)<-colnames(dataset)
      return(vec)})
      names(sepset)<-colnames(dataset)
      #名为Adjs的空列表
      Adjs <-vector("list",ncol(dataset))
      names(Adjs)<-colnames(dataset)
      
      for (T in colnames(dataset)){
        #message("T = ", T)
        S = character(0) # empty set
        adjT = character(0) # empty set
        for (v in setdiff(colnames(dataset),T)) {
          pval <- bnlearn:::indep.test(v, T, sx = S, test = test, data = dataset,
                                       B = 0L, alpha = alpha, complete = complete)
          # #message("v = ",v)
          # #message("pval = ",pval)
          if (abs(pval[1]) >= alpha){
            sepset[[T]][[v]] <- S
          }else{
            adjT <- append(adjT, v)
          }
        }
        # Sort adj(T ) in increasing corr(Vi, T ) value
        #相关系数升序排列
        sorted_adj <- sort(corr[,T][adjT])
        # #print(sorted_adj)
        # #print(names(sorted_adj))
        #sorted_adj 将成为一个包含排好序的数值名称的字符向量。
        sorted_adj <- names(sorted_adj)
        k = 1
        nNum = min(length(adjT),max.sx)
        if(nNum < 1){
          Adjs[[T]] <- character(0)
        }else{
          while(k <= nNum){
            # all subsets S \subset adj(T) needs to be tested whose
            # sizes do not exceed k
            for (v in sorted_adj) {
              #message("v = ", v)
              condset=setdiff(sorted_adj,v) 
              len <- length(condset)
              if(k > len){
                k <- nNum
                break
              }
              #for (j in 1:k) {
              a <- bnlearn:::allsubs.test(v, T, sx=condset, fixed = character(0), data=dataset, test=test, B = 0L,
                                          alpha = alpha, min = k, max = k, complete=complete, debug = FALSE)
              #print(a)
              if(a["p.value"] >= alpha){
                sorted_adj <- setdiff(sorted_adj,v)
                sepset[[T]][[v]] <- attr(a,"dsep.set")
              }
              #} 
            }
            k <- k+1
            nNum <- length(sorted_adj)
          }
          Adjs[[T]] <- sorted_adj
        }
      }# End of Alg 2
      # message("before symmetry correction")
      # print(Adjs)
      ## symmetry correction for removing false positives from adjV(T)
      ## i.e., if X \in adjV(Y) but Y \not\in adjV(X) then 
      ## adjV(Y) = adjV(Y)\{X} and sepset[[Y]][[X]] <- sepset[[X]][[Y]]
      for (var in colnames(dataset)) {
        for (u in Adjs[[var]]) {
          if(!(var %in% Adjs[[u]])){
            Adjs[[var]] <- setdiff(Adjs[[var]],u)
            sepset[[var]][[u]] <- sepset[[u]][[var]]
          }
        }
      }
      # message("after symmetry correction")
      # print(Adjs)
      return(list(adjs=Adjs,sepset = sepset))
    }#End adjecencies
    #############################################################################################
    ##################### MMBC-CSP Algorithm for Mb recovery ####################################
    #############################################################################################
    ## Input: T : the target variable; alpha: significant level
    ## Output: Mb(T )
    
    `findMb` <- function(T, dataset, adjs, sepset,alpha = alpha){
      adjT <- adjs[[T]]
      # bd(T) U ch(T) \subsetq Mb(T)
      mmbT <- adjT
      # complex-spouses recovery phase
      # candidates of complex-spouses
      candids <- setdiff(colnames(dataset),union(T,adjT))
      for (ch in adjT) {
        for (csp in candids) {
          pval1 <- bnlearn:::indep.test(csp, T, sx = sepset[[T]][[csp]], test = test, data = dataset,
                                        B = 0L, alpha = alpha, complete = complete)
          pval2 <- bnlearn:::indep.test(csp, T, sx = union(sepset[[T]][[csp]],ch), test = test, data = dataset,
                                        B = 0L, alpha = alpha, complete = complete)
          if(abs(pval1[1]) > alpha & abs(pval2[1]) <= alpha){
            mmbT <- union(mmbT,csp)
          }
        }
      }
      ## false positives phase (shrinking phase)
      continue <- TRUE
      delcand <- mmbT    # only those added later is allowed to be deleted
      if (length(delcand) == 0)
        continue <- FALSE
      while (continue) {       # shrink the Markov blanket
        p.val <- sapply(delcand, function(x)
          bnlearn:::indep.test(T, x, sx = setdiff(mmbT,x), test = "zf", data = dataset,
                               B = 0L, alpha = alpha, complete = complete))
        # multinom.ci.test(freq.tb, var, x, setdiff(mmbT,x))$p.value)
        ## this step could be speeded up significantly!!!
        p.val.max <- max(p.val)
        idx <- which(p.val == p.val.max)[1]
        if(p.val.max > alpha) {
          mmbT <- setdiff(mmbT, delcand[idx])
          delcand <- delcand[-idx]
        } else {
          continue <- FALSE
        }
      }
      return(mmbT)
    }#End findMb
    
    # 1. [Compute Markov Blankets]
    result <- adjecencies(dataset = dataset,alpha = alpha)
    nam <- names(dataset)
    mb = bnlearn:::smartSapply(cluster, as.list(nam), findMb, dataset = dataset, 
                               adjs = result$adjs, sepset = result$sepset, alpha=alpha)
    names(mb) = nam
    p <- ncol(dataset)
    G <- matrix(0, p, p)
    colnames(G)<-rownames(G)<-nam
    # a list similar to sepset in pcalg package
    # initialize all elements = 0
    # sepset<-apply(G, 1, as.list)
    
    # 2. [Compute Graph Structure]
    ###############################################################
    # # check symmetry in the output of the algorithm
    # message("before symmetry correction")
    # print(mb)
    ## symmetry correction for removing false positives from Mb(T)
    ## i.e., if X \in Mb(Y) but Y \not\in Mb(X) then 
    ## Mb(Y) = Mb(Y)\{X} and sepset[[Y]][[X]] <- sepset[[X]][[Y]]
    for (var in colnames(dataset)) {
      for (u in mb[[var]]) {
        if(!(var %in% mb[[u]])){
          mb[[var]] <- setdiff(mb[[var]],u)
          #sepset[[var]][[u]] <- sepset[[u]][[var]]
          G[var,mb[[var]]] <- G[mb[[var]],var] <- 1
        }
      }
    }#End check symmetry
    
    #return(list(mb,G))

    s=adjecencies(dataset = dataset,alpha = alpha)$sepset
    
     for (i in 0:p-2){
       for (u in colnames(dataset)){
        for (v in mb[[u]]){
           if(length(v)>i+1){
            a <- bnlearn:::allsubs.test(v, u, sx=setdiff(mb[[u]],v), fixed = character(0), data=dataset, test=test, B = 0L,
                                          alpha = alpha, complete=complete, debug = FALSE)
              #print(a)
            if(a["p.value"] >= alpha){
              s[[u]][[v]] <- s[[u]][[v]] <- attr(a,"dsep.set")
              G[u,v] <- G[v,u] <- 0
             }
           }
         }
       }
     }
     return(list(s=s,G=G))
  }
  
  learn.complex <- function(dataset){
    G <- learn.skeleton(dataset)$G
    s <- learn.skeleton(dataset)$s
    data.info = bnlearn:::check.data(dataset, allow.missing = TRUE)
    complete=data.info$complete.nodes
    for (u in colnames(complete)){
      for (v in setdiff(colnames(complete),u)){
        if(G[u,v] == 0){
          Suv <- s[[u]][[v]]
          
          for (w in setdiff(colnames(complete),u)){
            if(G[u,w]  == 1 && G[w,u]  == 1){
              res <- bnlearn:::indep.test(u, v, sx = unique(append(Suv, w)), test = "zf", data = dataset,
                                          B = 0L, alpha = 0.05, complete = complete)
              
              if (res < alpha) {
                G[w, u] <-  0
              }
            }
          }
          
#          nb.u <- vset[G == 1]
 #         nb.u.size <- length(nb.u)
 #         if (nb.u.size > 0) {
  #          for (j in 1:nb.u.size) {
   #           w <- nb.u[j]
    #          newsep <- unique(append(Suv, w))
     #         idx <- c(u, v, newsep)
              #res <- norm.ci.test(cov[idx, idx], n, u, v, newsep)
              
      #        res <- bnlearn:::indep.test(u, v, sx = newsep, test = "zf", data = dataset,
       #                            B = 0L, alpha = alpha, complete = complete)
              
        #      if (res < alpha) {
         #       G[w, u] <-  0
          #    }
           # }
       #   }
        }
      }
    }
    pattern(G)
    return(G)
  }


`learn.skeleton.norm` <- function(tree, cov, n, p.value, drop = TRUE)
{
  validObject(tree)
  local.ug <- c()
  vset <- rownames(cov)
  n.clique <- length(tree@cliques)
  for(i in 1:n.clique){
    idx <- tree@cliques[[i]]$vset
    #        if (length(idx) >= 10)
    new.ug <- .get.localug.pc(cov[idx, idx], n, p.value)
    #        else
    #            new.ug <- .get.localug.ic(cov[idx, idx], n, p.value)
    local.ug <- append(local.ug, new.ug)
  }
  p <- length(vset)
  amat <- matrix(0, p, p)
  rownames(amat) <- colnames(amat) <- vset
  n.clique <- length(tree@cliques)
  for(i in 1:n.clique){
    idx <- tree@cliques[[i]]$vset
    amat[idx, idx] <- 1
  }
  diag(amat) <- 0
  sep.pairs <- c()
  n.loc.sep <- length(local.ug)
  if(n.loc.sep>0)
    for(i in 1:n.loc.sep){
      u <- local.ug[[i]]@u
      v <- local.ug[[i]]@v
      if(amat[u,v] == 1){
        amat[u,v] <- amat[v,u] <- 0
        sep.pairs <- append(sep.pairs, local.ug[i])
      }
    }
  
  ## the following code is partially adapted from the "pcAlgo" function
  ## from "pcalg" package in R
  
  if (drop) {
    #ind获取可以被删除的潜在边对
    ind <- .get.exed.cand1(tree, amat)
    if (any(ind)) {
      #对ind中的数据按照第一列进行升序排序
      ind <- ind[order(ind[,1]),]
      ord <- 0
      seq_p <- 1:p
      done <- FALSE
      #获取潜在可删除边对的数量，行数
      remainingEdgeTests <- nrow(ind)
      while (!done && any(as.logical(amat))) {
        done <- TRUE
        for (i in 1:remainingEdgeTests) {
          x <- ind[i, 1]
          y <- ind[i, 2]
          if (amat[y, x]) {
            nbrsBool <- amat[, x] == 1
            nbrsBool[y] <- FALSE
            nbrs <- seq_p[nbrsBool]
            length_nbrs <- length(nbrs)
            if (length_nbrs >= ord) {
              if (length_nbrs > ord)
                done <- FALSE
              S <- seq(length = ord)
              repeat {
                p.val <- norm.ci.test(cov, n, vset[x], vset[y],
                                      vset[nbrs[S]])$p.value
                if (p.val > p.value) {
                  amat[x, y] <- amat[y, x] <- 0
                  pair <- new("sep.pair", u = vset[x],
                              v = vset[y], s = vset[nbrs[S]])
                  sep.pairs <- append(sep.pairs, pair)
                  break
                }
                else {
                  nextSet <- .getNextSet(length_nbrs, ord, S)
                  if (nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              }
            }
          }
        }
        ord <- ord + 1
      }
    }
  }
  return(list(amat=amat, sep.pairs=sep.pairs))    
}

`learn.complex.norm` <- function(skel, cov, n, p.value)
{
  wmat <- skel$amat
  vset <- rownames(wmat)
  sep.pairs <- skel$sep.pairs
  n.sep <- length(sep.pairs)
  if (n.sep == 0) return(wmat)
  for (i in 1:n.sep) {
    pair <- sep.pairs[[i]]
    for (turn in 1:2) {
      u <- if(turn == 1) pair@u else pair@v
      v <- if(turn == 1) pair@v else pair@u
      sep <- pair@s
      nb.u <- vset[skel$amat[u,] == 1]
      nb.u.size <- length(nb.u)
      if (nb.u.size > 0) {
        for (j in 1:nb.u.size) {
          w <- nb.u[j]
          newsep <- unique(append(sep, w))
          idx <- c(u, v, newsep)
          res <- norm.ci.test(cov[idx, idx], n, u, v, newsep)
          if (res$p.value < p.value &&
              (-1 - res$deviance) < wmat[w,u]) {
            wmat[w, u] <-  -1 - res$deviance
          }
        }
      }
    }
  }
  idx <- which(wmat - t(wmat) < 0)
  wmat[idx] <- 0
  wmat <- 0 + (wmat != 0)
  pattern(wmat)
}