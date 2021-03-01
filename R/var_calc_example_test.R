#' Replace missing values with the mean
#'
#' @param numeric.vector A vector with at least one numeric value.
#'
#' @return The input whose NAs are replaced by the input's mean.
#'
#' @export
#'
#' @examples fill_na_mean(c(1, 2, 3, NA, 5, 3, NA, 6))
fill_na_mean <- function(numeric.vector) {
  ifelse(
    is.na(numeric.vector),
    mean(numeric.vector, na.rm = TRUE),
    numeric.vector
  )
}
