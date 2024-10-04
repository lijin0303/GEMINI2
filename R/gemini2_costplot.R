#' GEMINI2 cost plotting function
#'
#' @description The function visualizes a series of costs with dot plot
#' @param input A series of costs (numeric)
#' @return A ggplot2 object
#' @export
#' @import ggplot2
#' @examples
#' p <- gemini2_costplot(infered_D$cost)
#'
gemini2_costplot <- function(costs){
  df <- data.frame(x = 1:(length(costs)-1),y =costs[-1])
  costp <- ggplot(df, aes(x , y )) +
    geom_point() + geom_line()+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          axis.text = element_text(color='black'))+
    labs(x="Iteration number",y="Maximized cost function")
  return(costp)
}
