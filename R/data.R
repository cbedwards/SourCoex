#' Sourdough Data
#'
#' A data frame of timeseries data from pairwise competition experiments of sourdough communities from Landis et al. 2020.
#'
#' Note from the data file of Landis et al:
#'   "CFU counts and relative abundance data from competitions, transfers 1, 3, and 6 - Pairs were transferred every 48 hours through filtered flour media (see materials and methods)."
#'
#' @format A data frame with 1080 rows and 6 variables.
#' \describe{
#'   \item{abund}{Are the values relative ("rel") or absolute ("abs") abundance?}
#'   \item{transf}{What is the transfer number? Multiply by 48 for hours since experiment start}
#'   \item{spec1}{Identifier of first competing species. See "spec.map" to get genus/species from the identifier}
#'   \item{spec2}{Identifier of second competing species.}
#'   \item{abund1}{relative or absolute abundance of species 1}
#'   \item{abund2}{relative or absolute abundance of species 2}
#' }
#' @source \url{https://elifesciences.org/articles/61644}
#' @seealso spec.map
"sour.data"


#' Map of sourdough species
#'
#' A companion to sour.data, matching species ID to the genus and species name.

#' @format A data frame with 8 rows and 2 variables.
#' \describe{
#'   \item{name.data}{Identifier used in sour.data. Note that numbers correspond to fungi, letters to bacteria.}
#'   \item{name.full}{Scientific name}
#' }
#' @source \url{https://elifesciences.org/articles/61644}
#' @seealso sour.data
"spec.map"
