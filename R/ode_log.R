#' Logistic Growth (1-species)
#'
#' ODE model for logistic growth. Note that is mostly an internal file: it gets used by `lsoda` in `exper_pred()`, and we don't have to interact with it directly.
#'
#' Note that the names of the parameters is provided in `parmnamesLog`, and the units of the parameters is provided in `unitsLog`. These are both used in plotting functions (to generate nice tables of parameters).
#'
#'
#' @param t time (only relevant for non-autonomous ODEs, not for us, but necessary for `deSolve` to play nice.
#' @param y current abundance of the population
#' @param parms parameters for the ODE model. In order:
#' \itemize{
#'    \item r - intrinsic growth rate
#'    \item k - carrying capacity
#'}
#' @return List of the derivate at current population density.
#'
#' @export


ode_log=function(t,y,parms) {
  # Yes, this can be solved analytically, but I want to make sure the method is working.
  #state variables:
  x=max(y[1],0);
  # Parameters:
  r=parms[1];
  k=parms[2]
  # Model:
  dx= r*x*(1 - x/k)
  dY=c(dx);
  return(list(dY));
}
# for an ode function: make vector of parameter names in order
# and make vector of units in order, using plotmath notation.
parmnamesLog=c("r","k")
unitsLog=c("growth per capita per hour",
            "individuals"
)
usethis::use_data(parmnamesLog, unitsLog, overwrite = TRUE)

#' List of parameter names for the logistic growth model
#'
#' A companion to `ode_log`.
#' @format A vector of 2 strings.

#' @seealso `ode_log`
"parmnamesLog"

#' List of parameter units for the logistic growth model
#'
#' A companion to `ode_log`.
#' @format A vector of 2 strings.

#' @seealso `ode_log`
"unitsLog"
