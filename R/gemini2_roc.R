#' GEMINI2 ROC curve plotting function
#'
#' @description The function plots ROC curve
#' @param confusion A confusion table obtained from \code{\link[gemini2]{gemini2_confusion}} (data.frame)
#' @param label The column name for specifying group variable
#' @return A ggplot2 object
#' @export
#' @examples
#' roc_curve <- gemini2_roc(confusion = confusion_m)
#'
gemini2_roc<- function(confusion,label=NULL,single=T){
  if(!single){
    grouplab <- enquo(label)
    roc_auc <- confusion%>%
      filter(!is.nan(TPR)&!is.nan(FPR))%>%
      group_by(!!grouplab)%>%
      summarise(roc_auc =(simple_auc(TPR,FPR)%>%round(.,4)))%>%
      as.data.frame()
    roc <- ggplot(data = confusion, aes(x = FPR, y = TPR,color = !!grouplab)) +
      geom_segment(data = data.frame(), aes(x = 0, y = 0, xend = 1, yend = 1), inherit.aes = F, color = "grey", linetype = 3) +
      scale_color_brewer(palette = "Set1", direction = 1, name = "Method") +
      geom_path() +
      xlim(0,1) +
      ylim(0,1) +
      geom_line(size=1) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(colour = "black"),
            panel.border = element_rect(color = "black", size = 1, fill = NA))}
  else{
    roc_auc <- confusion%>%
      filter(!is.nan(TPR)&!is.nan(FPR))%>%
      summarise(roc_auc =(simple_auc(TPR,FPR)%>%round(.,4)))%>%
      as.data.frame()
    roc <- ggplot(data = confusion, aes(x = FPR, y = TPR)) +
      geom_segment(data = data.frame(), aes(x = 0, y = 0, xend = 1, yend = 1), inherit.aes = F, color = "grey", linetype = 3) +
      geom_path() +
      xlim(0,1) +
      ylim(0,1) +
      geom_line(size=1) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(colour = "black"),
            panel.border = element_rect(color = "black", size = 1, fill = NA))
  }
  return(list(roc_auc=roc_auc,roc = roc))
}
