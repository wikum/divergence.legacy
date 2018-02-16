
###
### divergence computation functions
###

### ====================================================
### quantiles
### ====================================================

quantileTransform = function(x){
  
  #if(rank_samples)
  #  x = rank(x, ties.method="average")
  
  # start with rank 0, i.e. all values with min(x) will have rank 0
  # the percentile transformation depends on the non-zero ranks
  x = x - min(x)
  
  y = x
  # separate out the positive entries and apply transformation
  z = rank(x[x > 0], ties.method="min")
  y[x > 0] = (z - min(z))/sum(x > 0)
  y
}


### apply quantiles transformation to a matrix
###     rank_samples: TRUE if ranking should be applied before converting to percentiles
getQuantileMat = function(Mat){
  
  newMat = apply(Mat, 2, function(x) quantileTransform(x))
  
  rownames(newMat) = rownames(Mat)
  colnames(newMat) = colnames(Mat)
  
  newMat
}

### ====================================================
### get expected proportion of divergent feature per sample
### ====================================================

getAlpha = function(Mat, Ranges){
  
  mean(colSums(abs(computeTernary(Mat=Mat, R=Ranges)))/nrow(Mat))
  
}

### ====================================================
### computing ranges
### ====================================================

computeSingleRange = function(x, gamma, beta, method="partial", j=NULL){
  

  if(is.null(j))
    j = max(floor(gamma * length(x)), 1)
  
  # get distance to j'th nearest neighbor from each point in x
  if(method=="fullsort"){
    xjn = sapply(1:length(x), function(i){
      sort(abs(x[-i]-x[i]))[j]
    })
  }else if(method=="partial"){
    xjn = sapply(1:length(x), function(i){
      sort(abs(x[-i]-x[i]), partial=j)[j]
    })
  }else if(method=="kdtree"){
    #xjn = sapply(1:length(x), function(i) 
    #  max(nn2(data=matrix(x[-i], ncol=1), query=c(x[i]), k=j, eps=0)$nn.dists) )
    
    w = WKNNF(x)
    xjn = sapply(1:length(x), function(i) 
      max( w$query(x[i], k=j+1, eps=0)$nn.dists ))
    
  }else if(method=="quantile_dist"){
    xjn = apply(as.matrix(dist(x, method="euclidean", diag=FALSE)), 1, function(y) 
      quantile(y, probs=gamma))
  }else if(method=="ranks"){
    xjn = sapply(1:length(x), function(i){
      y = abs(x[-i]-x[i])
      y[which(rank(y, ties.method="random") == j)]
    })
  }
  
  sel = which(xjn <= quantile(xjn, beta))
  
  x = x[sel]
  xjn = xjn[sel]
  
  r = c(
    max(x[which.min(x)] - xjn[which.min(x)], 0),
    min(x[which.max(x)] + xjn[which.max(x)], 1)
  )
  
  list(range=r, support=sel)
  
}

getRangeList = function(Mat, gamma=0.1, beta=0.95, method="partial", par=TRUE,
                         nmax=200, mmax=1000, verbose=TRUE){

	L = NULL

	if(par){

		require(parallel)

		if(nrow(Mat) > mmax && ncol(Mat) > nmax){

			# for large matrices, break-up into chunks and process
			# each setion separately with mclapply

			L = getRangesBySections(Mat=Mat, gamma=gamma, beta=beta, method=method,
    			par=par, nmax=nmax, mmax=mmax, verbose=verbose)

		}else{

			# process in parallel by each row

			L = mclapply(1:nrow(Mat), function(i){
      	tempx = NULL
      	tryCatch({
        	tempx = computeSingleRange(Mat[i, ], gamma=gamma, beta=beta, method=method, j=NULL)
      	}, error = function(e){warning(e)})
      	tempx
    	}, mc.cores=detectCores())

    	if(sum(sapply(L, is.null)) > 0 || length(L) < nrow(Mat)){

    		# did not retrurn all features; re-run without parallel
        if(verbose)
      		cat("Not all features returned; re-running without parallelization\n")

    		L = getRangeList(Mat=Mat, gamma=gamma, beta=beta, method=method,
    			par=FALSE, nmax=nmax, mmax=mmax)

    	}

		}

	}else{

		# run without any parallelization

		L = lapply(1:nrow(Mat), function(i){
      	tempx = NULL
      	tryCatch({
        	tempx = computeSingleRange(Mat[i, ], gamma=gamma, beta=beta, method=method, j=NULL)
      	}, error = function(e){warning(e)})
      	tempx
    	})

	}

	names(L) = rownames(Mat)

	L

}

getRangesBySections = function(Mat, gamma=0.1, beta=0.95, method="partial", par=TRUE,
                         nmax=200, mmax=1000, verbose=TRUE){

  dirnum = paste(sapply(1:5, function(i) sample(0:9, 1)), collapse="")
  
  dirname = sprintf("ranges_%s", dirnum)
  if(verbose)
    cat(sprintf("Creating directory %s for saving files\n", dirname))
  dir.create(dirname)
  

  # apply to each feature in the baseline data matrix

  nc = detectCores()
  #nc = floor(detectCores()/2)
  
  n = nrow(Mat)
  k = ceiling(n/mmax)
  
  slots = lapply(1:k, function(i){ 
    x = 1:mmax + ((i-1) * mmax) 
    x[x <= n]
  })
  
  #print(sapply(slots, function(x) c(x[1:2], x[(length(x)-1):length(x)], length(x)) ))
  
  for(j in 1:length(slots)){
    
    if(verbose)
    	cat(sprintf("Processing features %d to %d\n", min(slots[[j]]), max(slots[[j]]) ))
    
    partialRanges = getRangeList(Mat[slots[[j]], ], gamma=gamma, beta=beta, method=method, par=par, nmax=nmax, mmax=mmax)
    
    if(verbose)
      cat("Saving..\n")
    save(partialRanges, file=sprintf("%s/%d.rda", dirname, j))
    
    rm(partialRanges)
    gc()
  
  }
  
  # assemble
  if(verbose)
    cat("Assembling..\n")

  L = list()
  for(j in 1:length(slots)){
    
    ll = load(sprintf("%s/%d.rda", filename, j))
    
    L[[j]] = partialRanges
    
    rm(list=ll)
    gc()
    
  }

  L

}

computeRanges = function(Mat, gamma=0.1, beta=0.95, method="partial", parallel=TRUE, verbose=TRUE){
  
  stopifnot(is.matrix(Mat))

  if((! is.numeric(gamma)) || (gamma <= 0) || (gamma >= 1) ){
    warning("Invalid gamma value; using gamma=0.1\n")
    gamma = 0.1
  }

  if((! is.numeric(beta)) || (beta <= 0) || (beta >= 1) ){
    warning("Invalid beta value; using beta=0.95\n")
    beta = 0.95
  }

	if(method == "kdtree"){
		if("nabor" %in% installed.packages())
			library("nabor")
		else{
			warning("Package nabor for computing kd-trees not found; hence method=partial will be used")
			method = "partial"
		}
	}

  if(verbose){
    cat(sprintf("Computing ranges from %d reference samples for %d features\n", 
                ncol(Mat), nrow(Mat)))
    cat(sprintf("[beta=%.3f, gamma=%.3f]\n", beta, gamma))
  }

  # apply to each feature in the baseline data matrix
  L = getRangeList(Mat=Mat, gamma=gamma, beta=beta, method=method, par=parallel, verbose=verbose)
  
  R = data.frame(t(sapply(L, function(x) x$range)))
  colnames(R) = c("baseline.low", "baseline.high")
  
  S = t(sapply(L, function(x) 1*(1:ncol(Mat) %in% x$support) ))
  colnames(S) = colnames(Mat)
  
  alpha=getAlpha(Mat=Mat, Ranges=R)
  if(verbose)
    cat(sprintf("[Expected proportion of divergent features per sample=%.3f]\n", alpha))
  
  list(Ranges=R, 
       Support=S,
       gamma=gamma,
       alpha=alpha)
  
}

### ====================================================
### search for gamma
### ====================================================

findGamma = function(Mat,
	gamma=c(1:9/100, 1:9/10),
  beta=0.95, 
  alpha=0.01,
  method="partial", 
  parallel=TRUE,
  verbose=TRUE){
 
  stopifnot(is.matrix(Mat))

 if(verbose)
    cat(sprintf("Searching optimal gamma for alpha=%.5f\n", alpha))
   
  optimal_gamma = -1
  
  gamma = sort(gamma)
  names(gamma) = paste("g", 1:length(gamma), sep="")
  
  g = c()
  e = rep(NA, length(gamma))
  names(e) = names(gamma)
    
  RangesList = list()
  for(i in 1:length(gamma)){
    L = computeRanges(Mat=Mat, gamma=gamma[i], beta=beta, method=method, parallel=parallel, verbose=verbose)
    
    RangesList[[i]] = L
    e[i] = L$alpha
    
    g = c(g, gamma[i])
    if(e[i] <= alpha)
      break
  }
  names(RangesList) = names(gamma)[1:length(RangesList)]
  names(g) = names(gamma)[1:length(g)]

  optimal = FALSE
  if(length(which(e <= alpha)) < 1){
    sel = which.max(g)
  }else{
    sel = which(g == min(g[which(e <= alpha)]))
  }
  optimal_gamma = g[sel]
  R_star = RangesList[[ sel ]]
  e_star = e[sel]
  
  if(e_star <= alpha)
    optimal = TRUE
  
  temp_e = e
  names(temp_e) = g
  #print(temp_e)

  if(verbose)
    cat(sprintf("Search results for alpha=%.5f: gamma=%.5f, expectation=%.5f, optimal=%s\n", alpha, optimal_gamma, e_star, optimal))
  
  list(Ranges=R_star$Ranges, 
       Support=R_star$Support,
       gamma=optimal_gamma,
       alpha=e_star,
       optimal=optimal,
       alpha_space=data.frame(gamma=gamma, alpha=e)
       )
  
}

### ====================================================
### compute ternary form
### ====================================================

computeTernary = function(Mat, R){
  
  lower = "baseline.low"
  upper = "baseline.high"

  if(nrow(Mat) != nrow(R))
    cat(sprintf("WARNING [%d, %d]\n", nrow(Mat), nrow(R)))
  
  DMat = ((Mat < R[, lower]) * (-1)) + ((Mat > R[, upper]) * 1)
  rownames(DMat) = rownames(Mat)
  colnames(DMat) = colnames(Mat)
  DMat
  
}

### ====================================================
### compute divergences
### ====================================================

computeTernaryDigitization = function(Mat, baseMat, 
                              computeQuantiles=TRUE,
                              gamma=c(1:9/100, 1:9/10),
                              beta=0.95, 
                              alpha=0.01,
                              method="partial",
                              parallel=TRUE,
                              verbose=TRUE,
                              findGamma=TRUE, 
                              Groups=NULL,                             
                              classes=NULL){

  stopifnot(rownames(Mat) == rownames(baseMat))

  if(! is.null(Groups)){

    stopifnot(length(Groups) == ncol(Mat))

    if(! is.factor(Groups))
      Groups = factor(Groups)

    if(is.null(classes))
       classes = levels(Groups)
 
  }

  if(! is.numeric(alpha)){
    warning("Invalid alpha value\n")
    findGamma = FALSE
  }else if(alpha < 0 || alpha > 1){
    warning("Invalid alpha value\n")
    findGamma = FALSE
  }

  if(computeQuantiles){
    if(verbose)
      cat(sprintf("Computing quantiles..\n"))
    baseMat = getQuantileMat(baseMat)
    Mat = getQuantileMat(Mat)
  }

  if(findGamma){
    L = findGamma(Mat=baseMat, gamma=gamma, beta=beta, alpha=alpha, method=method, parallel=parallel, verbose=verbose)
  }
  else{
    if(verbose)
      cat(sprintf("Using gamma=.4f\n", gamma[1]))
    L = findGamma(Mat=baseMat, gamma=gamma[1], beta=beta, alpha=alpha, method=method, parallel=parallel, verbose=FALSE)
  }

  DMat_ternary = computeTernary(Mat=Mat, R=L$Ranges)
  
  baseMat_ternary = computeTernary(Mat=baseMat, R=L$Ranges)

  DMat = abs(DMat_ternary)
  
  D = rowMeans(DMat)
  N = colSums(DMat)
  
  Npos = colSums(DMat_ternary > 0)
  Nneg = colSums(DMat_ternary < 0)
  
  df = data.frame(feature=rownames(Mat), prob.div=D)

  if(! is.null(Groups)){

    classDiv = sapply(classes, function(x) rowMeans(DMat[, which(Groups == x)]))

    df = data.frame(df, classDiv)
    colnames(df) = c("feature", "prob.div", paste("prob.div.", classes, sep=""))

  }

  list(Mat.div=DMat_ternary,
      baseMat.div = baseMat_ternary,
      div = data.frame(sample=colnames(Mat), count.div=N, count.div.upper=Npos, count.div.lower=Nneg),
      features.div = df,
      Ranges = L$Ranges,
      Support = L$Support,
      gamma = L$gamma,
      alpha = L$alpha,
      optimal = L$optimal,
      alpha_space = L$alpha_space)
    
}

