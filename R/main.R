#' Compute percentile transformations
#'
#' Function for computing the percentile transformation for one or more samples
#' supplied as columns of a matrix.
#'
#' @param Mat Matrix of data, with each column corresponding to a sample and each row corresponding to a feature.
#'
#' @return A matrix of the same dimensions as Mat with the percentile data.
#'
#' @export
#'
#' @examples
#' data(divergence)
#' percentileMat = computePercentileMatrix(breastTCGA_Mat)
#'

computePercentileMatrix <- function(Mat, rank_samples=TRUE){
	getPercentileMat(Mat=Mat, rank_samples=rank_samples)
}

# ==================================================================================================
#' Estimate the baseline support
#'
#' Function for computing the basline support for univariate features given gamma
#' and beta parameters.
#'
#' @param Mat Matrix of data, with each column corresponding to a sample and each 
#' row corresponding to a feature
#' @param gamma Parameter for selecting radius around each support point (0 < gamma < 1).
#' By default gamma = 0.1.
#' @param beta Parameter for eliminating outliers (0 < beta < 1). By default beta=0.95.
#' @param method Method of computing distance to j'th nearest neighbour. Options are 
#' "partial" (default option and usually the fastest) for using partial quicksort, 
#' "fullsort" for using full quicksort, "kdtree" for using a kd-tree (requires the package nabor), 
#' "quantile_dist" for computing quantiles over pairwise distances, "ranks" for using the rank function.
#' @param par Logical indicating whether to compute features parallelly with mclapply (defaults to TRUE).
#' @param e Logical indicating whether to return the expected number of divergent features per sample
#' estimated over the samples in Mat (defaults to TRUE).
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#'
#' @return A list with elements "Ranges" (a data frame containing the baseline interval for each feature),
#' "Support" (a binary matrix of the same dimensions as Mat indicating whether each sample was a support 
#' for a feature or not, 1=support, 0=not in the support), and "ex" (the expected number of divergent 
#' features per sample estimated over the samples in Mat if e=TRUE, NULL otherwise).
#' 
#' @keywords baseline, support
#' @export
#' 
#' @examples
#' data(divergence)
#' S = computeSupport(breastTCGA_Mat[, breastTCGA_Group=="NORMAL"])
#'

computeSupport <- function(Mat, gamma=0.1, beta=0.95, method="partial", par=TRUE, e=TRUE, verbose=TRUE){
	computeRanges(Mat=Mat, gamma=gamma, beta=beta, method=method, par=par, e=e, verbose=verbose)
}

# ==================================================================================================
#' Search for optimal gamma
#'
#' Function for searching through a range of gamma values for finding the smallest gamma that 
#' provides expected proportion of divergent features per sample less than or equal to alpha.
#'
#' @param Mat Matrix of data, with each column corresponding to a sample and each 
#' row corresponding to a feature
#' @param gamma Range of gamma values to search through. 
#' By default gamma = {0.01, 0.02, ... 0.09, 0.1, 0.2, ..., 0.9}.
#' @param beta Parameter for eliminating outliers (0 < beta < 1). By default beta=0.95.
#' @param alpha Expected proportion of divergent features per sample to be estimated
#' over the samples in Mat. By default alpha = 0.01; i.e. search for the smallest gamma that provides
#' 1% or less number of divergent features per sample. 
#' @param method Method of computing distance to j'th nearest neighbour. Options are 
#' "partial" (default option and usually the fastest) for using partial quicksort, 
#' "fullsort" for using full quicksort, "kdtree" for using a kd-tree (requires the package nabor), 
#' "quantile_dist" for computing quantiles over pairwise distances, "ranks" for using the rank function.
#' @param par Logical indicating whether to compute features parallelly with mclapply (defaults to TRUE).
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#'
#' @return A list with elements "Ranges" (a data frame containing the baseline interval for each feature for the),
#' selected gamma, "Support" (a binary matrix of the same dimensions as Mat indicating whether each sample was a support 
#' for a feature or not, 1=support, 0=not in the support), "Expectation" (the expected number of divergent 
#' features per sample estimated over the samples in Mat), "gamma" (vector of gamma values searched), "e" (vector of expected
#' proportion of divergent features per sample for each gamma value in gamma), "optimal_gamma" (selected gamma value), 
#' "optimal" (logical indicaing whether the selected gamma value provided the necessary alpha requirement), and "gamma_all" (vector
#' of gamma values provided to the function).
#' 
#' @keywords gamma
#' @export
#' 
#' @examples
#' data(divergence)
#' S = findGamma(breastTCGA_Mat[, breastTCGA_Group=="NORMAL"])
#'
findOptimalGamma <- function(Mat, gamma=c(1:9/100, 1:9/10), beta=0.95, alpha=0.01, method="partial", par=TRUE, verbose=TRUE){
     findGamma(Mat=Mat, gamma=gamma, beta=beta, alpha=alpha, method=method, par=par, verbose=verbose)
}

# ==================================================================================================
#' 
#' Compute ternary digitization
#'
#' Function for obtaining the digitized form, along with other relevant statistics and measures 
#' given a data matrix and a baseline matrix 
#'
#' @param Mat Matrix of data to be digitized (usually in quantiles), with each column corresponding to a sample and each 
#' row corresponding to a feature
#' @param baseMat Matrix of baseline data (usually in quantiles), with each column corresponding to a sample and each 
#' row corresponding to a feature
#' @param classes 
#'
#' @param gamma Range of gamma values to search through. 
#' By default gamma = {0.01, 0.02, ... 0.09, 0.1, 0.2, ..., 0.9}.
#' @param beta Parameter for eliminating outliers (0 < beta < 1). By default beta=0.95.
#' @param method Method of computing distance to j'th nearest neighbour. Options are 
#' "partial" (default option and usually the fastest) for using partial quicksort, 
#' "fullsort" for using full quicksort, "kdtree" for using a kd-tree (requires the package nabor), 
#' "quantile_dist" for computing quantiles over pairwise distances, "ranks" for using the rank function.
#' @param par Logical indicating whether to compute features parallelly with mclapply (defaults to TRUE).
#' @param alpha Expected proportion of divergent features per sample to be estimated. The optimal gamma providing
#' this level of divergence in the baseline data will be searched for.
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#'
#' @keywords digitize ternary
#' @export
#'
#' @examples
#' data(divergence)
#' D = computeDigitization(breastTCGA_Mat[, breastTCGA_Group=="NORMAL"], breastTCGA_Mat[, breastTCGA_Group!="NORMAL"])
#'

computeDigitization <- function(Mat, baseMat, 
                              classes=NULL, 
                              gamma=c(1:9/100, 1:9/10),
                              beta=0.95, 
                              method="partial",
                              par=TRUE,
                              alpha=0.01,
                              verbose=TRUE){

	computeDivergences(Mat=Mat, refMat=baseMat, 
		classes=classes, gamma=gamma, beta=beta, method=method,
		par=par, alpha=alpha, verbose=verbose)

}



