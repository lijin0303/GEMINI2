#' GEMINI2 Precision-Recall curve plotting function
#'
#' @description The function plots Precision-Recall curve
#' @param confusion A confusion table obtained from \code{\link[gemini2]{gemini2_confusion}} (data.frame)
#' @param label The column name for specifying group variable
#' @return A ggplot2 object
#' @export
#' @import RColorBrewer
#' @examples
#' PR_curve <- gemini2_PR(confusion = confusion_m)
#'
gemini2_PR<- function(confusion,label=NULL,single=T){
  if(!single){
    grouplab <- enquo(label)
    PR_auc <- confusion%>%
      filter(!is.nan(PPV)&!is.nan(TPR))%>%
      group_by(!!grouplab)%>%
      summarise(roc_auc =(simple_auc(PPV,TPR)%>%round(.,2)))%>%
      as.data.frame()
    PR <- ggplot(data = confusion, aes(x = TPR, y = PPV,color = !!grouplab)) +
      labs(color = "Method", x = "Recall", y = "Precision") +
      geom_hline(yintercept = 0.5, color = "grey", linetype = 'dashed') +
      geom_path() +
      scale_color_brewer(palette = "Set1", direction = 1, name = "Method") +
      xlim(0,1) +
      ylim(0,1) +
      geom_line(size=1) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(colour = "black"),
            panel.border = element_rect(color = "black", size = 1, fill = NA))}
  else{
    PR_auc <- confusion%>%
      filter(!is.nan(PPV)&!is.nan(TPR))%>%
      summarise(roc_auc =(simple_auc(PPV,TPR)%>%round(.,2)))%>%
      as.data.frame()
    PR <- ggplot(data = confusion, aes(x = TPR, y = PPV)) +
      labs(color = "Method", x = "Recall", y = "Precision") +
      geom_hline(yintercept = 0.5, color = "grey", linetype = 'dashed') +
      geom_path() +
      xlim(0,1) +
      ylim(0,1) +
      geom_line(size=1) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(colour = "black"),
            panel.border = element_rect(color = "black", size = 1, fill = NA))
  }
  return(list(PR_auc=PR_auc,PR = PR))
}
