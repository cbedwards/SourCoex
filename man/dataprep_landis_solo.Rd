% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataprep_landis_solo.R
\name{dataprep_landis_solo}
\alias{dataprep_landis_solo}
\title{Prep Landis et al. 2021 data (ie \code{sourData}) for model-fitting}
\usage{
dataprep_landis_solo(specid, data.use = sourData, aug.ls = list())
}
\arguments{
\item{specid}{Single character of the species to prep.}

\item{data.use}{Data frame of data to use. Almost assuredly sourData}

\item{aug.ls()}{Augmentation list. Currently used for setting initial conditions for Tilman experiments}
}
\value{
Returns list with several components, suitable for analysis and plotting
\itemize{
\item dat - dataframe of the data for this species, for plotting.
\item dat.real.ls - list form of the data for this species, used to evaluate model fit with \code{obj_helper()}.
\item x0.mat - matrix of initial conditions, with each row corresponding to one experiment. Same order as dat.real.ls, also used in \code{obj_helper()}.
}
}
\description{
This function subsets the Landis data for single-species model-fitting, and returns a matrix of appropriate initial conditions, as well as various formattings of the raw data for analysis and plotting. Plays well with \code{obj_helper()}, \code{plotter_landis_solo()}. Alternative to \code{dataprep_landis()} for fitting logistic growth, etc.
}
