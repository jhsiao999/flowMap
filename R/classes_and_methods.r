#' FR statistics and p-values generated from one single draw of the population pair comparison
#'
#' This class stores the statistics required to compute median FR statistics across random draws.
#'
#' @section FIXME:
#' ## usage
#' ## Accessors
#' getFRstats(object)
#' getPnorm(object)
#'
#' @name FRstats-class
#' @rdname frstats-class
#' @author Chiaowen Joyce Hsiao \email{joyce.hsiao1@@gmail.com}
#' @examples
#' ## see vignettes
#'
#' @export
setClass("FRstats",
         representation=representation(XX1="data.frame",
           XX2="data.frame",
           sampleMethod="character",
           sampleSize="numeric",
           ncores="numeric",
           ndraws="numeric",
           npop1="numeric",
           npop2="numeric",
           pop1Labels="numeric",
           pop2Labels="numeric",
           ww="matrix",
           runs="matrix",
           mu="matrix",
           sigma2="matrix",
           pNorm="matrix"))


setMethod("show","FRstats",function(object) { 
            cat("FR statistics object on a",object@npop1,"by",object@npop2,"two-sample comparison\n")
            cat("Based on",object@draws,"random draws\n")
            cat("with", object@sampleMethod, "sampling method\n")
            cat("of", object@sampleSize, "events from each group\n")
          })
