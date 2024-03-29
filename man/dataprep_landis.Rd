% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataprep_landis.R
\name{dataprep_landis}
\alias{dataprep_landis}
\title{Prep Landis et al. 2021 data (ie \code{sourData}) for model-fitting}
\usage{
dataprep_landis(specid, data.use = sourData)
}
\arguments{
\item{specid}{Vector of the first and second species for the analysis. Order matters.}

\item{data.use}{data frame of data to use. Almost assuredly sourData}
}
\value{
Returns list with several components, suitable for analysis and plotting
\itemize{
\item dat.comp - data of only the 2-species competition experiments. For \code{plotter_landis()}, joint plot.
\item dat.solo1 - data of only the single-species competition experiments for species 1. For \code{plotter_landis()}, species 1 plot.
\item dat.solo2 - data of only the single-species competition experiments for species 2. For \code{plotter_landis()}, species 2 plot.
\item dat.real.ls - combination of dat.comp, dat.solo1, dat.solo2. Used to evaluate model fit with \code{obj_helper()}.
\item x0.mat - matrix of initial conditions, with each row corresponding to one experiment. Same order as dat.real.ls, also used in \code{obj_helper()}.
}
}
\description{
This function subsets the Landis data, and returns a matrix of appropriate initial conditions, as well as various formattings of the raw data for analysis and plotting. Plays well with \code{obj_helper()}, \code{plotter_landis()}.
}
