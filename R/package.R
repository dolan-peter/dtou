#' dtou: A package for computing the distance to uniqueness metric for a collection of genomic sequences
#'
#' The only function that should be directly invoked is \code{dtou()} which knows how to call the appropriate C functions
#'
#' @section C functions:
#' \code{c_dtou} takes a vector of characters and returns their distance to uniqueness.
#'
#' \code{c_dtouDepthLimit} is a depth limited version of the previous function.
#'
#' \code{c_dtouS2} is a variation that may improve speed.
#'
#' \code{c_dtouS2DepthLimit} is a depth limited version that may also improve speed.
#' @docType package
#' @name dtou
NULL

