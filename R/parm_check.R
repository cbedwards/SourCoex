#' Check parameters
#'
#' This function checks to see if fitted parameters are at or near the user-defined parameters (returns informative warnings if so). This is useful for identifying when solvers may be "caught" at the boundaries, or when boundaries are too constraining.
#'
#' @param parms.cur Current estimate (from `optim()` etc)
#' @param parm.lower Lower bound that user supplied for `optim()` call
#' @param parm.upper Upper bound that user supplied
#' @param parmnames vector of the parameter names. Function uses this to identify boundary parameters
#' @param tol Tolerance parameter. Smaller values mean only estimates extremely close to the boundary are identified. Default is 0.001
#'
#' @return Nothing. Will output warnings if there are boundary parameters
#'
#' @export


parm_check = function(parms.cur, #current estimate
                      parm.lower, #lower bounds allowed
                      parm.upper, #upper bounds allowed
                      parmnames,
                      tol=0.001){ #name of parameters
  ## function to check if any parameters are at boundaries
  if(any(abs((parms.cur - parm.lower)/((parm.lower+parm.upper)/2))<tol)){
    warning(paste0("One or more parameters very near lower bounds:",
                   "\n   ",paste0((parmnames[abs((parms.cur - parm.lower)/((parm.lower+parm.upper)/2))<tol]), collapse="; ")
    ))
  }
  if(any(abs((parms.cur - parm.upper)/((parm.lower+parm.upper)/2))<tol)){
    warning(paste0("One or more parameters very near upper bounds:",
                   "\n   ",paste0((parmnames[abs((parms.cur - parm.upper)/((parm.lower+parm.upper)/2))<tol]), collapse="; ")
    ))
  }
}
