#' Sourdough Data
#'
#' A data frame of timeseries data from pairwise competition experiments of sourdough communities from Landis et al. 2020.
#'
#' Note from the data file of Landis et al:
#'   "CFU counts and relative abundance data from competitions, transfers 1, 3, and 6 - Pairs were transferred every 48 hours through filtered flour media (see materials and methods)."
#'
#'   Important note in their methods: initial abundance is 2000 CFUs in 200 microliters, either 1:1 for intraspecies competition, or  entirely one species in the one-species versions.
#'     (I have added this initial condition to the data)
#'
#'   Other important note in their methods: at each transfer, 5 microliters is transfered to 195 microliters of new medium.
#'   So when modeling, the start of "transfer 4" period should be 2.5% of densities counted in the transfer 3.
#'
#' @format A data frame with 1440 rows and 7 variables.
#' \describe{
#'   \item{abund}{Are the abundance measures relative ("rel") or absolute ("abs")? In a modification from Landis 2020 data, relative abundance is no longer in percents, and instead is out of 1. Absolute measures are in CFUs/microliter}
#'   \item{spec1}{Identifier of first competing species. See "specMap" to get genus/species from the identifier}
#'   \item{spec2}{Identifier of second competing species. If same as spec1, these are single-species competition experiments}
#'   \item{rep}{Which relicate (out of 5) was this? Unique identifier for the experimental "community lineage"}
#'   \item{transf}{What is the transfer number? Each transfer is 48 hours long, and the start of each transfer dilutes to 2.5% of the end of the last one.}
#'   \item{abund1}{relative or absolute abundance of species 1}
#'   \item{abund2}{relative or absolute abundance of species 2}
#' }
#' @source \url{https://elifesciences.org/articles/61644}
#' @seealso specMap
#' @export
"sourData"


#' Map of sourdough species
#'
#' A companion to sourData, matching species ID to the genus and species name.

#' @format A data frame with 8 rows and 2 variables.
#' \describe{
#'   \item{name.data}{Identifier used in sourData. Note that numbers correspond to fungi, letters to bacteria.}
#'   \item{name.full}{Scientific name}
#' }
#' @source \url{https://elifesciences.org/articles/61644}
#' @seealso sourData
#' @export
"specMap"
