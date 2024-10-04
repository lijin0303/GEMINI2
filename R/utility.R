#' Matrix conversion
#' @name matrix_conv
#' @param m Data.frame whose 1st column will serve as rownames
#' @return Converted matrix
#' @examples matrix_conv(init_D$Y)
#'
matrix_conv <- function(m){
  mm <- as.matrix(m[,-1])
  rownames(mm) <- m[,1]
  return(mm)}

#' Fix a level to the last in sequence for factor subject
#' @name fixed_levels
#' @param s Vector to be  factorized
#' @param f The target level to be fixed in the last
#' @return Factorized variable
#' @examples fixed_levels(s=c("a","b","d","e"),f="b")
#'
fixed_levels <- function(s,f){
  us <- unique(s)
  l <- c(us[us!=f],f)
  fs <- factor(s,levels=l)
  return(fs)
}

#' AUC calculation with series of axis values
#' @name simple_auc
#' @param y Ordered values on Y axis
#' @param x Ordered values on X axis
#' @return Calculated simple auc
#' @examples simple_auc(y = seq(1:20),x=seq(1:20))
#'
simple_auc <- function(y, x){
  # inputs already sorted, best scores first
  dy <- c(diff(y), 0)
  dx <- c(diff(x), 0)
  sum(y * dx) + sum(dy * dx)/2
}
