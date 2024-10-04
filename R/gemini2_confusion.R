#' GEMINI2 confusion matrix generating function
#'
#' @description The function calculates the confusion matrix for llustrating the diagnostic ability of a binary classifier system as its discrimination threshold is varied.
#' @param score The list of scoring results from \code{\link[gemini2]{gemini2_score}} (list)
#' @param pos_ctrl A character vector of gene pairs serving as positive controls (character)
#' @param neg_ctrl A character vector of gene pairs serving as negative controls (character)
#' @param cut_seq A series of discrimination thresholds (numeric)
#' @param hit_cl The number of cell lines where the score of an interaction at least needs to fulfill discrimination threshold
#' in order to be counted as true positive  (numeric)
#' @param cl A character specifying a cell line to calculate confusion matrix (character)
#' @details Between parameters \code{hit_cl} and \code{cl}, only 1 need to be specified by the user
#' @return A confusion matrix
#' @export
#' @examples
#' confusion_m <- gemini2_confusion(score=scores$sensitive_lethality,
#' pos_ctrl = c("AURKB;PLK4","BAP1;UCHL5","BCL2L1;BCL2L2","BCL2L1;MCL1"),
#' neg_ctrl = c("CCL13;CCL8","CDY1;CDY1B","CDY1;CDY2A"),cut_seq=seq(10, -10, by = -0.01),hit_cl=1,cl=NULL)
#'
gemini2_confusion <- function(score,pos_ctrl,neg_ctrl,cut_seq=NULL,hit_cl=1,cl=NULL,direction="geq"){
  score%<>%replace(is.na(.), 0)
  all_ctrls <- if(length(neg_ctrl)==0){rownames(score)}else{union(pos_ctrl,neg_ctrl)}
  if(is.null(cut_seq)&direction=="geq"){
    cut_seq <- seq(ceiling(max(score)),ceiling(min(score))-1,by= -0.005)
  }else if(is.null(cut_seq)&direction=="leq"){
    cut_seq <- seq(ceiling(min(score))-1,ceiling(max(score)),by=0.005)
    }
  confusion_seq <- lapply(cut_seq,function(t){
    score_ctrls <-  if(is.null(cl)){score[all_ctrls,]}else{score[all_ctrls,cl]%>%as.data.frame()}
    score_ctrls%<>%
      apply(.,2,function(x){ x[is.na(x)] <- min(x, na.rm = T)
      return(x)})
    our_hits <- switch(direction,
      "geq" = rownames(score_ctrls)[rowSums(score_ctrls >= t, na.rm = T) >= hit_cl],
      "leq" = rownames(score_ctrls)[rowSums(score_ctrls <= t, na.rm = T) >= hit_cl])
    confusion_D <- data.frame(P = length(our_hits),
                              N = nrow(score_ctrls) - length(our_hits),
                              TP = length(intersect(pos_ctrl, our_hits)),
                              FN = length(setdiff(pos_ctrl,our_hits)),
                              cutoff = t)%>%
      mutate(FP = P - TP,
             TN = N - FN)})%>%
    bind_rows()%>%
    mutate(TPR =  TP/(TP+FN),
           TNR = TN/(TN+FP),
           FPR = FP/(FP+TN),
           PPV = TP/(TP+FP))
  return(confusion_seq)}
