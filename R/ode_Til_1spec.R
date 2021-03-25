#' Tilman competition (1-species + resource)
#'
#' ODE model for Tilman competition with inert medium, for a single species. Here we know that the resource isn't growing, so *r* in the classic Tilman model is 0.
#'
#' Note that I believe the amount of resources available at the start of a transfer can be set to some constant (doesn't need to be fit). This is handled in the exper_pred, but for simplicity I will set it to 195, for the mL of new resource.
#'
#' **IMPORTANT!!** For this model to work, need to pass *aug.ls$Til = 1* to `exper_pred()` (e.g. through the chain of parent functions).
#' *aug.ls$TilR* is the amount of resources provided at the start of each transfer. If undefined, defaults to 195 (e.g. Landis et al. 2021 numbersw).
#'
#' Note that the names of the parameters is provided in `parmnamesTil1spec`, and the units of the parameters is provided in `unitsTil1spec`. These are both used in plotting functions (to generate nice tables of parameters).
#'
#'
#' @param t time (only relevant for non-autonomous ODEs, not for us, but necessary for `deSolve` to play nice.
#' @param y vector of current abundance of species 1, species 2
#' @param parms parameters for the ODE model. In order:
#' \itemize{
#'    \item r1 - intrinsic growth rate of species 1
#'    \item d1 - per capita mortality rate, species 1

#'}
#' @return List of derivates at current population densities.
#'
#' @export


ode_Til_1spec=function(t,y,parms) {
  # This function is based on the Tilman resouce model (e.g. R*) with 1 species
  #state variables:
  x1=max(y[1], 0); R = max(y[2], 0)
  # Parameters:
  r1=parms[1];
  d1 = parms[2];
  # Model:
  dx1 = x1 * (r1*R-d1)
  dR = -R*(r1*x1)
  dY=c(dx1,dR);
  return(list(dY));
}
# for an ode function: make vector of parameter names in order
# and make vector of units in order, using plotmath notation.
#' @export
parmnamesTil1spec=c("r[1]", "d[1]", "R")
#' @export
unitsTil1spec=c("growth per capita per unit resource",
           "per capita",
           "thousands of individuals"
)

usethis::use_data(parmnamesTil1spec, unitsTil1spec, overwrite = TRUE)

#' List of parameter names for the 1-species Tilman competition model
#'
#' A companion to `ode_Til_1spec`.
#' @format A vector of 3 strings.

#' @seealso `ode_Til_1spec`
"parmnamesTil1spec"

#' List of parameter units for the 1-species Tilman model competition model
#'
#' A companion to `ode_Til_1spec`.
#' @format A vector of 3 strings.

#' @seealso `ode_Til_1spec`
"unitsTil1spec"
