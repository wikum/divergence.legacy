
###
### divergence computation functions
###

### computeSingleGenesetSupport
### ====================================================
### computing support for single gene set
### input:
###         Mat        geneset matrix  (usually in quantile), rows correspond to genes, columns correspond to samples
###         geneset    a set of gene names in geneset
###         gamma
###         beta
###         method     chose from "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
### output:
###         Centermatrix     data matrix of samples used to build support
###         Radius           distances from center sample to support boundary
###         Alpha           the expected number of divergent features per sample comuted over Mat
### ====================================================


computeSingleGenesetSupport = function(Mat, geneset, gamma, beta, method="euclidean", j=NULL){

  geneset <- intersect(rownames(Mat), geneset)
  Mat_Matrix <- Mat[geneset, ]

  if(is.null(j))
    j = max(floor(gamma * ncol(Mat_Matrix)), 1)

  # get distance to j'th nearest neighbor from each point in x

  Mat_Dist <- as.matrix(dist(t(Mat_Matrix), method = method, diag = FALSE))

  xjn <- apply(Mat_Dist, 1, function(y){
    sort(y, decreasing = FALSE)[j+1]
  })

  ### sel  indicate whether a sample is used to build the support
  sel = which(xjn <= quantile(xjn, beta))

  Radius <- xjn[sel]
  Centermatrix <- Mat_Matrix[, sel]
  Mat_Dist_m2c <- Mat_Dist[, sel]


  ### compare with radius

  v<- apply(Mat_Dist_m2c,1, function(y){
    as.numeric(all(y > Radius))
  })


  prob <- sum(v==1)/length(v)

  Alpha <- prob

  list(Radius=Radius, Centermatrix= Centermatrix, Alpha = prob)

}

### computeSingleGenesetBinaryMatrix
### ====================================================
### computing Binary vector for single gene set
### input:
###         Mat             geneset matrix  (usually in quantile), rows correspond to genes, columns correspond to samples
###         geneset         a set of gene names used to build the support (same as in computeSingleGenesetSupport )
###         Centermatrix    data matrix of samples used to build support  (data frame from output of computeSingleGenesetSupport)
###         Radius          distance from center sample to the boundary    (data frame from output of computeSingleGenesetSupport)
###         method          chose from "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
### output:
###        Binary Vector   Vector of 0/1 indicating whether each sample is in the support
### ====================================================

computeSingleGenesetBinaryMatrix = function(Mat, geneset, Centermatrix, Radius, method){

  ### get common genes
  geneset <- intersect(rownames(Mat), geneset)
  Mat_geneset <- Mat[geneset,]

  ### dimension of Mat_geneset, Centermatrix
  mn_m <- dim(Mat_geneset)
  mn_c <- dim(Centermatrix)


  Mat_Matrix <- cbind(Centermatrix, Mat_geneset)

  Mat_Dist <- as.matrix(dist(t(Mat_Matrix), method = method, diag = FALSE))


  ### get distance from Mat to Centermatrix

  Mat_Dist_m2c <-  Mat_Dist[c((mn_c[2]+1):(mn_c[2]+mn_m[2])),c(1:mn_c[2])]

  ### compare with radius
  v <- apply(Mat_Dist_m2c,1, function(y){
    as.numeric(all(y > Radius))
  })


  BinaryMatrix = as.matrix(v)

  return(BinaryMatrix)

  }


###   computeMultiGenesetSupport
### ====================================================
### computing support for multiple gene set
### input:
###         Mat        geneset matrix  (usually in quantile), rows correspond to genes, columns correspond to samples
###         genematrix    matrix of 0/1 indicating whether gene is in geneset, rows correspond to genes, columns correspond to geneset,
###         gamma
###         beta       =0.95
###         method
### output:
###        list of  Centermatrix     data matrix of samples used to build support
###        list of  Radius           distances from center sample to support boundary
###        Support                   matrix of 0/1 indicating whether sample is in the union support built by multiple geneset
###        Alpha                     the expected number of divergent features per sample computed over the baseline data matrix (if alpha=TRUE in the input, NULL otherwise)
### ====================================================

computeMultiGenesetSupport = function(Mat, genematrix, gamma, beta, method){

  Radius_list <- list()
  Centermatrix_list <- list()
  Support <- c()
  Alpha <- c()

  for(i in 1:ncol(genematrix)){


    geneset <- rownames(genematrix)[which(genematrix[,i]==1)]
    S<- computeSingleGenesetSupport(Mat,geneset,gamma, beta, method=method, j=NULL)

    Centermatrix<- S$Centermatrix
    Radius <- S$Radius


    Bv <- computeSingleGenesetBinaryMatrix(Mat,geneset,Centermatrix,Radius, method=method)

    Centermatrix_list[[i]] <- S$Centermatrix
    Radius_list[[i]] <- S$Radius
    Support <- rbind(Support, t(Bv))
    Alpha <- c(Alpha,S$Alpha)
  }

  names(Centermatrix_list)<- colnames(genematrix)
  names(Radius_list)<- colnames(genematrix)
  rownames(Support) <- colnames(genematrix)

  Alpha <- mean(Alpha)

  list(Centermatrix_list = Centermatrix_list, Radius_list = Radius_list,
       Support = Support, Alpha=Alpha)

}



### findMultiSetGammaAndSupport
### ====================================================
### find  support and gamma for multiple gene set
### input:
###         Mat        geneset matrix  (usually in quantile), rows correspond to genes, columns correspond to samples
###         genematrix    matrix of 0/1 indicating whether gene is in geneset, rows correspond to genes, columns correspond to geneset,
###         gamma      = c(0.01,0.02,...0.9)
###         beta       =0.95
###         alpha      =0.1
###         method
### output:
###         list of elements
###         list of  Centermatrix     data matrix of samples used to build support
###         list of  Radius           distances from center sample to support boundary
###         Support                   matrix of 0/1 indicating whether sample is in the union support built by multiple geneset
###         Alpha                     the expeced number of divergent features per sample computed over the baseline data matrix (if alpha=TRUE in the input, NULL otherwise)
###         optimal_gamma             minimum gamma satisfy Alpha <= alpha
### ====================================================




findMultiSetGammaAndSupport = function(Mat, genematrix, gamma, beta, alpha, method){

  #### build support for multiple gene set
  FGS <- lapply(1:length(gamma), function(i){
     MultiSetSupport <- computeMultiGenesetSupport(Mat, genematrix, gamma[i], beta, method)
     cat(sprintf("[gamma=%.3f, alpha=%.3f]\n", gamma[i], MultiSetSupport$Alpha))
     MultiSetSupport
  })

  #### find the minimal gamma satisfies: Alpha_output < alpha

  Alpha_list <- lapply(1:length(gamma), function(i){
    FGS[[i]]$Alpha
  })

  ID<- min(which(Alpha_list <= alpha) )

  if(is.null(ID)){
    return(' no optimal gamma')

  }else{

    SS <- FGS[[ID]]
    optimal_gamma <- gamma[ID]

    return(list(MultiSetSupport =SS, optimal_gamma= optimal_gamma))

    }
}


### computeMultiSetBinaryMatrix
### ====================================================
### computing Binary matrix for multiple gene set
### input:
###         Mat                  geneset matrix  (usually in quantile), rows correspond to genes, columns correspond to samples
###         genematrix           matrix of 0/1 indicating whether gene is in geneset, rows correspond to genes, columns correspond to geneset,
###         Centermatrix_list    list of data matrix of samples used to build support  (data frame from output of computeMultiGenesetSupport)
###         Radius_list          distance from center sample to the boundary    (data frame from output of ccomputeMultiGenesetSupport)
### output:
###        BinaryMatrix          Matrix of 0/1 indicating whether each sample is in the support
### ====================================================



computeMultiSetBinaryMatrix = function(Mat, genematrix, Centermatrix_list, Radius_list, method){


  TM <- lapply(1:ncol(genematrix), function(i){
    geneset <- rownames(genematrix)[which(genematrix[,i]==1)]
    Centermatrix <- Centermatrix_list[[i]]
    Radius <- Radius_list[[i]]


   BM<- computeSingleGenesetBinaryMatrix(Mat,geneset, Centermatrix, Radius, method=method)

   BinaryMatrix<- BM

  })


  BinaryMatrix <- matrix(NA, nrow=ncol(genematrix), ncol=ncol(Mat))
  for(i in 1:ncol(genematrix)){
    BinaryMatrix[i,] <- as.matrix(TM[[i]])
  }

  rownames(BinaryMatrix) <- colnames(genematrix)
  colnames(BinaryMatrix) <- colnames(Mat)

  return(BinaryMatrix)


}




### computeMultiSetDigitization
### ====================================================
###
### input:
###         Mat               geneset matrix  (usually in quantile), rows correspond to genes, columns correspond to samples
###         basemat           baseline data matrix
###         genematrix        matrix of 0/1 indicating whether gene is in geneset, rows correspond to genes, columns correspond to geneset,
###         gamma   = c(0.01,0.02,...,.1,....9)
###         beta  =0.95
###         alpha =0.01
###         method ="euclidean"
###         findGamma = TRUE   if TRUE, use findMultiSetGammaAndSupport(), otherwise use the first gamma value given
### output:
###        list with elements
###        Mat.div              Binary data matrix
###        baseMat.div          Binary baseline data matrix
###        Centermatrix_list    list of data matrix of samples used to build support  (data frame from output of computeMultiGenesetSupport)
###        Radius_list          distance from center sample to the boundary    (data frame from output of ccomputeMultiGenesetSupport)
###        Alpha                the expected number of divergent features per sample comuted over the baseline
###        optimal_gamma        if findGamma = TRUE, then optimal_gamma will be output of findMultiSetGammaAndSupport, if findGamma=FALSE, the optimal_gamma will be first element of gamma
### ====================================================

computeMultiSetDigitization= function(Mat, basemat, genematrix,gamma,beta, alpha, method="euclidean", findGamma =TRUE){

  if (findGamma ==TRUE){

    FMSGS <- findMultiSetGammaAndSupport(basemat, genematrix, gamma, beta, alpha, method = method)

    Centermatrix_list <- FMSGS$MultiSetSupport$Centermatrix_list
    Radius_list <- FMSGS$MultiSetSupport$Radius_list
    optimal_gamma <- FMSGS$optimal_gamma

  }else{
    FMSGS <- computeMultiGenesetSupport(basemat, genematrix, gamma=gamma[1], beta, method=method)
    Centermatrix_list<-FMSGS$Centermatrix_list
    Radius_list <- FMSGS$Radius_list

    optimal_gamma <- gamma[1]

  }
  cat(sprintf("Using optimal gamma %.3f\n", optimal_gamma))

  Mat.div <- computeMultiSetBinaryMatrix(Mat, genematrix, Centermatrix_list, Radius_list, method=method)
  baseMat.div <- computeMultiSetBinaryMatrix(basemat, genematrix, Centermatrix_list, Radius_list, method=method)

  Alpha <-mean(apply(baseMat.div,1, function(x){sum(x==1)/length(x)}))

  list(Mat.div = Mat.div, baseMat.div = baseMat.div,
       Centermatrix_list = Centermatrix_list, Radius_list = Radius_list,
       optimal_gamma = optimal_gamma, Alpha = Alpha)

}

geneset_list_to_matrix = function(geneset_list){

  t(sapply(unique(unlist(geneset_list)), function(x) sapply(geneset_list, function(y) x %in% y ) )) * 1

}



