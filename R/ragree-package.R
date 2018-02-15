#' Rater agreement
#' 
#' Various measures for quantifying rater agreement.
#' 
#' @import ggplot2 irr stats
#' @importFrom grDevices adjustcolor
#' @importFrom graphics grid lines par plot polygon title
#' @importFrom utils write.table
#' @docType package
#' 
#' @name ragree-package
#' @aliases ragree
NULL

#' Disease severity rating of 27 patients by four ophthalmologists per person
#' 
#' Sample data in which 27 patients with a particular retinopathy had the
#' severity of their retinopathy rated by four ophthalmologists.
#' 
#' @seealso \code{\link{disease.w}}
#' @format An object of class \code{data.frame} containing 108 observations and
#' 3 variables:
#' 
#' \tabular{rll}{
#' \tab \code{patient} \tab patient number \cr
#' \tab \code{rater} \tab rater number \cr
#' \tab \code{category} \tab disease severity \cr
#' }
#' 
"disease.l"

#' Disease severity rating of 27 patients by four ophthalmologists per person
#' 
#' Sample data in which 27 patients with a particular retinopathy had the
#' severity of their retinopathy rated by four ophthalmologists.
#' 
#' @seealso \code{\link{disease.l}}
#' @format An object of class \code{data.frame} containing 27 observations and
#' 4 variables:
#' 
#' \tabular{rll}{
#' \tab \code{rater1} \tab opthalmologist 1 ratings \cr
#' \tab \code{rater2} \tab opthalmologist 2 ratings \cr
#' \tab \code{rater3} \tab opthalmologist 3 ratings \cr
#' \tab \code{rater4} \tab opthalmologist 4 ratings \cr
#' }
#' 
"disease.w"
