#' Compute quantile transformations
#'
#' Function for computing the quantile transformation for one or more samples
#' supplied as columns of a matrix.
#'
#' @param Mat Matrix of data, with each column corresponding to a sample and each row corresponding to a feature.
#'
#' @return A matrix of the same dimensions as Mat with the quantile data.
#'
#' @export
#'
#' @examples
#' data(divergence)
#' quantileMat = computeQuantileMatrix(breastTCGA_Mat)
#'

computeQuantileMatrix <- function(Mat){
	getQuantileMat(Mat=Mat)
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
#' @param parallel Logical indicating whether to compute features parallelly with mclapply (defaults to TRUE).
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#'
#' @return A list with elements "Ranges" (a data frame containing the baseline interval for each feature),
#' "Support" (a binary matrix of the same dimensions as Mat indicating whether each sample was a support 
#' for a feature or not, 1=support, 0=not in the support), and "alpha" (the expected number of divergent 
#' features per sample estimated over the samples).
#' 
#' @keywords baseline, support
#' @export
#' 
#' @examples
#' data(divergence)
#' S = computeUnivariateSupport(computeQuantileMatrix(breastTCGA_Mat[, breastTCGA_Group=="NORMAL"]))
#'

computeUnivariateSupport <- function(Mat, gamma=0.1, beta=0.95, method="partial", parallel=TRUE, verbose=TRUE){
	computeRanges(Mat=Mat, gamma=gamma, beta=beta, method=method, parallel=parallel, verbose=verbose)
}

# ==================================================================================================
#' Search for optimal gamma and associated support
#'
#' Function for searching through a range of gamma values for finding the smallest gamma that 
#' provides expected proportion of divergent features per sample less than or equal to alpha.
#'
#' @param Mat Matrix of data, with each column corresponding to a sample and each 
#' row corresponding to a feature
#' @param gamma Range of gamma values to search through. 
#' By default gamma = \{0.01, 0.02, ... 0.09, 0.1, 0.2, ..., 0.9\}.
#' @param beta Parameter for eliminating outliers (0 < beta < 1). By default beta=0.95.
#' @param alpha Expected proportion of divergent features per sample to be estimated
#' over the samples in Mat. By default alpha = 0.01; i.e. search for the smallest gamma that provides
#' 1\% or less number of divergent features per sample. 
#' @param method Method of computing distance to j'th nearest neighbour. Options are 
#' "partial" (default option and usually the fastest) for using partial quicksort, 
#' "fullsort" for using full quicksort, "kdtree" for using a kd-tree (requires the package nabor), 
#' "quantile_dist" for computing quantiles over pairwise distances, "ranks" for using the rank function.
#' @param parallel Logical indicating whether to compute features parallelly with mclapply (defaults to TRUE).
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#'
#' @return A list with elements "Ranges" (a data frame containing the baseline interval for each feature for the),
#' selected gamma, "Support" (a binary matrix of the same dimensions as Mat indicating whether each sample was a support 
#' for a feature or not, 1=support, 0=not in the support), "gamma" (selected gamma value), "alpha" (the expected number
#'  of divergent features per sample computed over the baseline data matrix), "optimal" (logical indicaing whether the 
#' selected gamma value provided the necessary alpha requirement), and "alpha_space" (a data frame with alpha values 
#' for each gamma searched).
#' 
#' @keywords gamma
#' @export
#' 
#' @examples
#' data(divergence)
#' S = findUnivariateGammaWithSupport(computeQuantileMatrix(breastTCGA_Mat[, breastTCGA_Group=="NORMAL"]))
#'
findUnivariateGammaWithSupport <- function(Mat, gamma=c(1:9/100, 1:9/10), beta=0.95, alpha=0.01, method="partial", parallel=TRUE, verbose=TRUE){
     findGamma(Mat=Mat, gamma=gamma, beta=beta, alpha=alpha, method=method, parallel=parallel, verbose=verbose)
}

# ==================================================================================================
#'
#' Compute the ternary matrix with digitized divergence coding
#'
#' Function for obtaining the ternary form for a matrix of data given a baseline range
#'
#' @param Mat Matrix of data to be digitized (usually in quantiles), with each column corresponding to a sample and each 
#' row corresponding to a feature
#' @param Ranges A data frame containing the baseline range of each features; this corresponds to the 'Ranges' element of 
#' the output of findUnivariateGammaWithSupport() and computeUnivariateSupport()
#'
#' @return A matrix with the same dimensions as Mat containing the ternary form data.
#'
#' @keywords ternary digitization
#' @export
#'
#' @examples
#' data(divergence)
#' S = computeUnivariateSupport(computeQuantileMatrix(breastTCGA_Mat[, breastTCGA_Group=="NORMAL"]))
#' T = computeUnivariateTernaryMatrix(Mat=computeQuantileMatrix(breastTCGA_Mat[, breastTCGA_Group!="NORMAL", ]), Ranges=S$Ranges)

computeUnivariateTernaryMatrix <- function(Mat, Ranges){
			computeTernary(Mat=Mat, R=Ranges)
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
#' @param computeQuantiles Apply quantile transformation to both data and baseline matrices (TRUE or FALSE; defaults to TRUE).
#' @param gamma Range of gamma values to search through. 
#' By default gamma = {0.01, 0.02, ... 0.09, 0.1, 0.2, ..., 0.9}.
#' @param beta Parameter for eliminating outliers (0 < beta < 1). By default beta=0.95.
#' @param alpha Expected proportion of divergent features per sample to be estimated. The optimal gamma providing
#' this level of divergence in the baseline data will be searched for.
#' @param method Method of computing distance to j'th nearest neighbour. Options are 
#' "partial" (default option and usually the fastest) for using partial quicksort, 
#' "fullsort" for using full quicksort, "kdtree" for using a kd-tree (requires the package nabor), 
#' "quantile_dist" for computing quantiles over pairwise distances, "ranks" for using the rank function.
#' @param parallel Logical indicating whether to compute features parallelly with mclapply (defaults to TRUE).
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#' @param findGamma Logical indicating whether to search for optimal gamma values through the given gamma values (defaults to 
#' TRUE). If FALSE, the first value given in gamma will be used.
#' @param Groups Factor indicating class association of samples
#' @param classes Vector of class labels
#'
#' @return A list
#'
#' @keywords digitize ternary
#' @export
#'
#' @examples
#' data(divergence)
#' D = computeDigitization(breastTCGA_Mat[, breastTCGA_Group=="NORMAL"], breastTCGA_Mat[, breastTCGA_Group!="NORMAL"])
#'

computeUnivariateDigitization <- function(Mat, baseMat, 
															computeQuantiles=TRUE,
                              gamma=c(1:9/100, 1:9/10),
                              beta=0.95, 
                              alpha=0.01,
                              method="partial",
                              parallel=TRUE,
                              verbose=TRUE,
                              findGamma=TRUE, 
                              Groups=NULL,                             
                              classes=NULL
){

	computeTernaryDigitization(Mat=Mat, baseMat=baseMat, computeQuantiles=computeQuantiles,
		gamma=gamma, beta=beta,
		alpha=alpha, method=method,
		parallel=parallel, verbose=verbose, 
		findGamma=findGamma, Groups=Groups, classes=classes)

}

# ==================================================================================================
#'
#' Compute chi-squared test
#'
#' Given a binary or ternary data matrix with class associations of samples, computes chi-squared tests
#' for each feature between two given classes
#'
#' @param Mat Matrix of digitized binary or ternary data with each column corresponding to a sample and each 
#' row corresponding to a feature
#' @param Groups Factor indicating class association of samples
#' @param classes Vector of class labels; the test will be applied between the first two classes given.
#'
#' @keywords chi-squared
#' @export
#'

computeChiSquaredTest <- function(Mat, Groups, classes){

	chiSquaredTest(Mat=Mat, Groups=Groups, classes=classes)

}



