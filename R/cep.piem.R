#' List of species from the Piedmont Region (Italy) that can be used in the seed mixture or donor grassland composition 
#'
#' A dataframe with the 258 plant species and their name abbreviation 
#' suitable for modeling when the "data" argument of \code{\link[ResNatSeed]{trainingDB}} function 
#' is NULL,i.e. are related to Piedmont Region, Italy.
#'
#' 
#' 
#' @docType data
#' 
#' @format A dataframe with 258 rows and 2 variables:
#' \describe{
#' \item{species}{Species names according to Aeschimann et al 2004}
#' \item{cep.names}{Abbreviation of species names according to the Cornell Ecology Programs (CEP), which uses eight-letter abbreviations for species. See \code{\link[vegan]{make.cepnames}} for more details} 
#' }
#' 
#' @usage data(cep.piem)
#' 
#' @references * Aeschimann, D., Lauber, K., Moser, M. D., & Theurillat, J. P. (2004). Flora alpina (Flora of the Alps). Zanichelli, Bologna.
#' 

"cep.piem"
