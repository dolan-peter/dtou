#' distance to uniqueness for collection of genetic sequences.
#'
#' \code{dtou} expects that genetic sequences will be represented as character vectors using only the
#' upper-case characters ACGTNX.  It performs no type-conversions so the input
#' must already be of type character.
#'
#' For each genetic sequence in the input parameter \code{str}, the function \code{dtou} returns the numeric
#' \emph{distance to uniqueness} vector for that sequence.
#'
#' the \emph{distance to uniqueness} is defined for each base-pair in a sequence as the length of the
#' shortest unique subsequence starting at that base-pair.
#'
#' The \code{rc} parameter determines whether or not reverse complements are considered when determining what
#' counts as unique.
#'
#' Currently \code{dtou} uses recursive function calls and if the collection of genetic sequences in \code{str}
#' contains overly-long repeats then the system can crash-- including R itself.
#'
#' There are three variations of the function to help control for this problem.  First, indicating a \code{depth}
#' Will depth-limit the recursion and cause the dtou metric to return a value no greater than that depth.  A depth
#' that is too high for the amount of memory in the call-back stack used by the C++ function will still crash the
#' system.
#'
#' Secondly, the \code{optimizeForSpeed} option introduces extra over-head into the recursion (thus decreasing the
#' maximum recursion level) by keeping track of the stack size and switching to an iterative approach when the
#' copy number of a repetitive region is \strong{two}.  (see the explanatory vignette for more information)
#' If all long repetitions have copy number two, then the recursion limit is not reached.  This technique is also
#' faster than the default-- but again, at the cost of reducing the maximum level of recursion.
#'
#' The third variation combines the previous two.  It is a depth limited optimize-for-speed approach.
#'
#'
#' @param str A vector of characters.  Contains genomic sequences comprised of A,C,G,T,X, and N only.
#' @param depth A numeric scalar. Limits the max value of the dtou metric.
#' @param rc A logical scalar. Indicates whether to take reverse complements into consideration.
#' @param optimizeForSpeed A logical scalar. Introduces overhead into the recursion but my increase speed
#' @return A list of numeric vectors containing the \emph{distance to uniqueness} value for each base-pair
#' @export
#'
#' @examples
#' dtou(c("AAAAACCCGACTGGGCTCA","ACCT"),rc=TRUE)
#' dtou(c("AAAAACCCGACTGGGCTCA","ACCT"),rc=FALSE)
#' dtou(c("AAAAACCCGACTGGGCTCA","ACCT"),rc=TRUE,depth=3)
#' dtou(c("AAAAACCCGACTGGGCTCA","ACCT"),depth=3,optimizeForSpeed=TRUE)
#'
dtou=function(str,depth=NULL,rc=TRUE,optimizeForSpeed=FALSE){
	if(class(optimizeForSpeed)!="logical"){stop("Parameter <optimizeForSpeed> must be a logical")}
	if(class(str)!="character"){stop("Parameter <str> must be a character")}
	if(length(str)==0){stop("Parameter <str> must be non-empty")}
	if(!is.null(depth)){
		if(class(depth)!="numeric"&&class(depth)!="integer"){stop("Parameter <depth> must be numeric")}
		if(length(depth)==0){stop("Parameter <depth> must be non-empty")}
		if(length(depth)>1){warning("Parameter <depth> has length greater than 1... defaulting to depth[1]")}
		depth=depth[1]
		if(depth<=0){stop("Parameter <depth> must be positive")}
	}
	if(length(optimizeForSpeed)==0){
		warning("Parameter <optimizeForSpeed> is empty...defaulting to FALSE")
		optimizeForSpeed=FALSE
	}

	if(length(optimizeForSpeed)>1){warning("Parameter <optimizeForSpeed> has length greater than 1...defaulting to optimizeForSpeed[1]")}

	optimizeForSpeed=optimizeForSpeed[1]

  if(is.null(depth) && !optimizeForSpeed){return(dtou::c_dtou(str,rc))}
	if(!is.null(depth) && !optimizeForSpeed){
		return(dtou::c_dtouDepthLimit(str,rc,depth))
	}
	if(is.null(depth) && optimizeForSpeed){return(dtou::c_dtouS2(str,rc))}
	if(!is.null(depth) && optimizeForSpeed){
		return(dtou::c_dtouS2DepthLimit(str,rc,depth))
	}
}
