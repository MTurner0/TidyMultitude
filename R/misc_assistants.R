#' @title Scaling and Centering of Rows in Matrix-like Objects
#' 
#' @description Extends `base::scale` to center and/or scale the rows of a numeric matrix.
#' 
#' @param x a numeric matrix(like object), such as an assay.
#' @param center a logical value indicating whether rows should be centered. 
#' @param scale a logical value indicating whether rows should be scaled.
#' 
#' @export
scale_rowwise <- function(x, center = TRUE, scale = TRUE) {
  # Will only center if specified
  if (center) {
    x <- sweep(x, 1, apply(x, 1, mean))
  }
  # Will only scale if specified
  if (scale) {
    scales <- apply(x, 1, sd)
    for(i in 1:nrow(x)) { x[i, ] <- x[i, ]/scales[i] }
  }
  return(x)
}