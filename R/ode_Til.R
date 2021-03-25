#' Tilman competition (2-species + resource)
#'
#' ODE model for Tilman competition with inert medium. Here we know that the resource isn't growing, so *r* in the classic Tilman model is 0.
#'
#' Note that I believe the amount of resources available at the start of a transfer can be set to some constant (doesn't need to be fit). This is handled in the exper_pred, but for simplicity I will set it to 195, for the mL of new resource.
#'
#' **IMPORTANT!!** For this model to work, need to pass *aug.ls$Til = 2* to `exper_pred()` (e.g. through the chain of parent functions).
#' *aug.ls$TilR* is the amount of resources provided at the start of each transfer. If undefined, defaults to 195 (e.g. Landis et al. 2021 numbersw).
#'
#' Note that the names of the parameters is provided in `parmnamesTil`, and the units of the parameters is provided in `units.LV`. These are both used in plotting functions (to generate nice tables of parameters).
#'
#'
#' @param t time (only relevant for non-autonomous ODEs, not for us, but necessary for `deSolve` to play nice.
#' @param y vector of current abundance of species 1, species 2
#' @param parms parameters for the ODE model. In order:
#' \itemize{
#'    \item r1 - intrinsic growth rate of species 1
#'    \item r2 - intrinsic growth rate of species 2
#'    \item d1 - per capita mortality rate, species 1
#'    \item d2 -  per capita mortality rate, species 2

#'}
#' @return List of derivates at current population densities.
#'
#' @export


ode_Til=function(t,y,parms) {
  # This function is based on the Tilman resouce model (e.g. R*)
  #state variables:
  x1=max(y[1], 0); x2=max(y[2], 0); R = max(y[3], 0)
  # Parameters:
  r1=parms[1];
  r2=parms[2];
  d1 = parms[3];
  d2 = parms[4];
  # Model:
  dx1 = x1 * (r1*R-d1)
  dx2 = x2 * (r2*R-d2)
  dR = -R*(r1*x1 + r2*x2)
  dY=c(dx1,dx2,dR);
  return(list(dY));
}
# for an ode function: make vector of parameter names in order
# and make vector of units in order, using plotmath notation.
#' @export
parmnamesTil=c("r[1]", "r[2]", "d[1]","d[2]", "R")
#' @export
unitsTil=c("growth per capita per unit resource",
           "growth per capita per unit resource",
           "per capita",
           "per capita",
           "resource units"
)

usethis::use_data(parmnamesTil, unitsTil, overwrite = TRUE)

#' List of parameter names for the Tilman competition (R*) model
#'
#' A companion to `ode_Til`.
#' @format A vector of 4 strings.

#' @seealso `ode_Til`
"parmnamesTil"

#' List of parameter units for the Tilman competition (R*) model
#'
#' A companion to `ode_Til`.
#' @format A vector of 4 strings.

#' @seealso `ode_Til`
"unitsTil"
