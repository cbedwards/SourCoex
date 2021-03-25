#' Lotka Volterra competition (2-species)
#'
#' ODE model for L-V competition. Note that is mostly an internal file: it gets used by `lsoda` in `exper_pred()`, and we don't have to interact with it directly.
#'
#' Note that the names of the parameters is provided in `parmnamesLV`, and the units of the parameters is provided in `unitsLV`. These are both used in plotting functions (to generate nice tables of parameters).
#'
#'
#' @param t time (only relevant for non-autonomous ODEs, not for us, but necessary for `deSolve` to play nice.
#' @param y vector of current abundance of species 1, species 2
#' @param parms parameters for the ODE model. In order:
#' \itemize{
#'    \item r1 - intrinsic growth rate of species 1
#'    \item r2 - intrinsic growth rate of species 2
#'    \item a12 - competition coefficient for species 2 on species 1
#'    \item a21 - competition coefficient for species 1 on species 2
#'    \item k1 - carrying capacity for species 1
#'    \item k2 - carrying capacity for species 2
#'}
#' @return List of derivates at current population densities.
#'
#' @export


ode_LV=function(t,y,parms) {
  # This function is based on the Lotka-Volterra competition model
  #state variables:
  x1=max(y[1],0); x2=max(y[2],0);
  # Parameters:
  r1=parms[1];
  r2=parms[2];
  a12 = parms[3]
  a21 = parms[4]
  k1=parms[5]
  k2=parms[6]
  # Model:
  dx1 = r1*x1*(1 - (x1 + a12*x2)/k1)
  dx2 = r2*x2*(1 - (x2 + a21*x1)/k2)
  dY=c(dx1,dx2);
  return(list(dY));
}
# for an ode function: make vector of parameter names in order
# and make vector of units in order, using plotmath notation.


#' @export
parmnamesLV=c("r[1]", "r[2]", "alpha[12]","alpha[21]","k[1]","k[2]")
#' @export
unitsLV=c("growth per capita per hour",
           "growth per capita per hour",
           "unitless (multiplier)",
           "unitless (multiplier)",
           "individuals",
           "individuals"
)

usethis::use_data(parmnamesLV, unitsLV, overwrite = TRUE)

#' List of parameter names for the Lotka Voltera competition model
#'
#' A companion to `ode_LV`.
#' @format A vector of 6 strings.

#' @seealso `ode_LV`
"parmnamesLV"

#' List of parameter units for the Lotka Voltera competition model
#'
#' A companion to `ode_LV`.
#' @format A vector of 6 strings.

#' @seealso `ode_LV`
"unitsLV"
