
`mbcsp` <- function(dataset,alpha, test = "zf", cluster = NULL,
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
      }
    }
  }#End check symmetry
  
  return(mb)
}


skel_recov <- function(dataset, mb, test = "zf", alpha) {
  
  sepset <- lapply(colnames(dataset), function(.) {vec<-vector("list",ncol(dataset))
  names(vec)<-colnames(dataset)
  return(vec)})
  names(sepset)<-colnames(dataset)
  
  for (i in colnames(dataset)) {
    for (j in setdiff(colnames(dataset), i)) {
      if (!(i %in% mb[[j]]) & !(j %in% mb[[i]])) {
        if (length(mb[[i]]) < length(mb[[j]])) {
          sepset[[i]][[j]] = mb[[i]]
        } else {
          sepset[[i]][[j]] = mb[[j]]
        }
      }
    }
  }
  
  ## Skeleton Recovery
  skel = matrix(0, ncol(dataset), ncol(dataset))
  rownames(skel) = colnames(dataset)
  colnames(skel) = colnames(dataset)
  for (i in colnames(dataset)) {
    for (j in mb[[i]]) {
      skel[j,i] = 1
    }
  }
  
  for (i in 0:(ncol(dataset) - 2)) {
    for (u in colnames(dataset)) {
      for (v in setdiff(colnames(dataset), u)) {
        if (skel[u,v] == 1 & skel[v, u] == 1 & (sum(skel[, u] != 0) - 1) >= i) {
          
          adu = colnames(dataset)[which(skel[, u] == 1)]
          #names(adu) = NULL
          condset = setdiff(adu, v)
          a <- bnlearn:::allsubs.test(u, v, sx=condset, fixed = character(0), data=dataset, test=test, B = 0L,
                                      alpha = alpha, min = i, max = i, complete=TRUE, debug = FALSE)
          if(a["p.value"] >= alpha){
            # adjT <- setdiff(adjT,v)
            sepset[[u]][[v]] <- attr(a,"dsep.set")
            sepset[[v]][[u]] <- attr(a,"dsep.set")
            skel[u,v] = 0
            skel[v,u] = 0
          }
        }
      }
    }
  }
  
  return(list(skel=skel, sepset=sepset))
}



learn.vstruct <- function(dataset, mb, test = "zf",alpha) {
  skel = skel_recov(dataset, mb, test = test, alpha = alpha)
  sepset = skel$sepset
  wmat = skel$skel
  vset = colnames(dataset)
  q = length(vset)
  for (i in 1:(q-1)) {
    for (j in (i+1):q) {
      for (l in 1:q) {
        u = vset[i]
        v = vset[j]
        w = vset[l]
        if (skel$skel[u,v] == 0 && wmat[u,w] == 1 && wmat[v,w] == 1 &&
            wmat[w,u] + wmat[w,v] != 0) {
          pval1 = bnlearn:::indep.test(u, v, sx = union(sepset[[u]][[v]],w), test = test, data = dataset,
                                      B = 0L, alpha = alpha, complete = TRUE)
          pval2 = bnlearn:::indep.test(u, v, sx = sepset[[u]][[v]], test = test, data = dataset,
                                      B = 0L, alpha = alpha, complete = TRUE)
          if (pval1 < alpha && pval2 > alpha) {
            wmat[w,u] = wmat[w,v] = 0
          }
        }
      }
    }
  }
  vstruct = t(skel$skel - wmat)
  #wmat<-studeny_rules(wmat)$matrix
  #return(list(pattern = wmat, vstruct = vstruct))
  
  pattern = studeny_rules(wmat)$matrix
  # 将列名从字符串转换为数字
  new_col_names <- as.numeric(sub("X", "", colnames(pattern)))
  # 设置新的列名
  colnames(pattern) <- new_col_names
  rownames(pattern) <- new_col_names
  
  return(list(pattern=pattern, vstruct = vstruct))
  #return(list(pattern = studeny_rules(wmat)$matrix, vstruct = vstruct))
  
}



