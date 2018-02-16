
### ====================================================
### compute chi-squared test
### ====================================================

get_chi_sq_test = function(x, ...){
  res = list(statistic=NA, p.value=NA)
  tryCatch({
    if(ncol(x) > 1)
      res = chisq.test(x, ...)
  })
  res
}

get_chi_v = function(x, y){
  u = unique(c(x, y))
  sapply(u, function(ui) c(sum(x == ui), sum(y == ui)) )
}

get_chi_test = function(XMat, YMat, ...){
  C = sapply(1:nrow(XMat), function(j) get_chi_sq_test(get_chi_v(XMat[j,], YMat[j,]), ...), 
             simplify = FALSE)
  
  names(C) = rownames(XMat)
  C
}


chiSquaredTest = function(Mat, Groups, classes){
  
  c = get_chi_test(XMat=Mat[, which(Groups == classes[1])], YMat=Mat[, which(Groups == classes[2])])
  pvals = sapply(c, function(x) x$p.value)
  stats = sapply(c, function(x) x$statistic)

  pr = sort(pvals, decreasing=FALSE, index.return=TRUE, na.last=TRUE, method='radix')$ix

  df = data.frame(statistic=stats, pval=pvals)
  rownames(df) = rownames(Mat)

  df[pr, ]

}
