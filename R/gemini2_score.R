#' GEMINI2 genetic interaction scoring function
#'
#' @description The function scores genetic interactions based on inference results
#' @param opt The list of inference results from \code{\link[gemini2]{gemini2_inference}} (list)
#' @param pc_threshold A numeric value to indicate the LFC corresponding to a positive control,
#' if no pc_gene is specified. Default is -1 (numeric)
#' @param pc_gene A character vector of any length naming genes to use as positive controls (character)
#' @param pc_weight A weighting applied to the positive control
#' (\code{pc_gene}) to filter genes whose individual phenotype is
#' more lethal than \eqn{pc_weight*y(pc_gene)}. Default is 0.5 (numeric)
#' @param nc_pairs A set of non-interacting gene pairs to define statistical significance. (character)
#' @return An object of class gemini.score containing score values for
#' strong interactions and sensitive lethality and recovery, and if \code{nc_pairs} is specified, statistical significance for each scoring type.
#' @export
#' @import mixtools
#' @examples
#' scores <- gemini2_score(opt =infered_D, pc_gene = "EEF2")
#'
gemini2_score <- function(opt,pc_threshold = -1,pc_gene=NULL,pc_weight=0.5,nc_pairs = NULL){
  S <- opt$Sf
  Y <- opt$Yf
  gene_join <- gsub("[A-Z,0-9]","",rownames(S)[1])
  S_map <- data.frame(spairs = rownames(S),stringsAsFactors = F)%>%
    mutate(g1 = gsub(paste0(gene_join,".*"),"",spairs),
           g2 = gsub(paste0(".*",gene_join),"",spairs))
  message("Calculating GEMINI scores ...")
  score_strong <- abs(S)-pmax(abs(Y[S_map$g1,]),abs(Y[S_map$g2,]))%>%
    set_rownames(S_map$spairs)
  score_sensi_lethal <- (pmin(Y[S_map$g1,],Y[S_map$g2,])-(Y[S_map$g1,]+Y[S_map$g2,]+S))%>%
    set_rownames(S_map$spairs)
  score_sensi_recov <- -score_sensi_lethal
  pcg_flag <- !is.null(pc_gene)
  pcw_flag <- !is.null(pc_threshold)
  Y1 <- Y[S_map$g1,]
  Y2 <- Y[S_map$g2,]
  cl_thresholds <- if (pcg_flag){
    if(pcw_flag){warning("Both pc_gene and pc_threshold specified - defaulting to pc_gene.")}
    apply(Y[pc_gene, ],2,median)*pc_weight
  }else if(!pcg_flag & pcw_flag){
    rep(pc_threshold,ncol(Y))
  }else{apply(Y,2,quantile,probs = 0.01)*pc_weight
  }
  remove_lethal<- Y1 < cl_thresholds[col(Y1)] | Y2 < cl_thresholds[col(Y2)]
  remove_recovery <- Y1 > cl_thresholds[col(Y1)] & Y2 > cl_thresholds[col(Y2)]
  score_sensi_lethal%<>%replace(.,remove_lethal,NA)
  score_sensi_recov%<>%replace(.,remove_recovery,NA)
  message("Calculating p values and FDR ...")
  ncp_flag <- !is.null(nc_pairs)
  if (ncp_flag) {
    if (!requireNamespace("mixtools")) {
      stop("Please install the mixtools package: install.packages('mixtools').")
    }
    # construct negetaive model using mixture of two normal distributions
    # strong
    nc_values <- as.vector(score_strong[nc_pairs, ])
    invisible(capture.output({nc_model_strong <- mixtools::normalmixEM(nc_values)}))
    # lethality
    nc_values <- as.vector(score_sensi_lethal[nc_pairs, ])%>%na.omit()
    invisible(capture.output({nc_model_lethality <- mixtools::normalmixEM(nc_values)}))
    # recovery
    nc_values <- as.vector(score_sensi_recov[nc_pairs, ])%>%na.omit()
    invisible(capture.output({nc_model_recovery <- mixtools::normalmixEM(nc_values)}))

    # calculate p-values
    # strong
    pval_strong <- nc_model_strong$lambda[1] * pnorm(
      score_strong,
      mean = nc_model_strong$mu[1],
      sd = nc_model_strong$sigma[1],
      lower.tail = FALSE
    ) +
      nc_model_strong$lambda[2] * pnorm(
        score_strong,
        mean = nc_model_strong$mu[2],
        sd = nc_model_strong$sigma[2],
        lower.tail = FALSE
      ) # null hypothesis: nc_model_strong > score_strong
    # lethality
    pval_sensi_lethal <- nc_model_lethality$lambda[1] * pnorm(
      score_sensi_lethal,
      mean = nc_model_lethality$mu[1],
      sd = nc_model_lethality$sigma[1],
      lower.tail = FALSE
    ) +
      nc_model_lethality$lambda[2] * pnorm(
        score_sensi_lethal,
        mean = nc_model_lethality$mu[2],
        sd = nc_model_lethality$sigma[2],
        lower.tail = FALSE
      ) # null hypothesis: nc_model_lethality > score_sensi_lethal
    # recovery
    pval_sensi_recov <- nc_model_recovery$lambda[1] * pnorm(
      score_sensi_recov,
      mean = nc_model_recovery$mu[1],
      sd = nc_model_recovery$sigma[1],
      lower.tail = FALSE
    ) +
      nc_model_recovery$lambda[2] * pnorm(
        score_sensi_recov,
        mean = nc_model_recovery$mu[2],
        sd = nc_model_recovery$sigma[2],
        lower.tail = FALSE
      ) # null hypothesis: nc_model_recovery > score_sensi_recov
  }

  # calculate the adjusted p-values
  if (ncp_flag) {
    fdr_strong <- apply(pval_strong, 2, function(x)
      p.adjust(x, method = "fdr"))
    fdr_sensi_lethal <- apply(pval_sensi_lethal, 2, function(x)
      p.adjust(x, method = "fdr"))
    fdr_sensi_recov <- apply(pval_sensi_recov, 2, function(x)
      p.adjust(x, method = "fdr"))
  }

  if(is.null(nc_pairs)){
    output <- list(
      score_strong = score_strong,
      score_sensitive_lethal = score_sensi_lethal,
      score_sensitive_recovery = score_sensi_recov
    )
  }else{
    output <- list(
      score_strong = score_strong,
      score_sensitive_lethal = score_sensi_lethal,
      score_sensitive_recovery = score_sensi_recov,
      nc_model_strong = nc_model_strong,
      nc_model_lethality = nc_model_lethality,
      nc_model_recovery = nc_model_recovery,
      pval_strong = pval_strong,
      pval_sensitive_lethal = pval_sensi_lethal,
      pval_sensitive_recovery = pval_sensi_recov,
      fdr_strong = fdr_strong,
      fdr_sensitive_lethal = fdr_sensi_lethal,
      fdr_sensitive_recovery = fdr_sensi_recov
    )
  }
  return(output)
}
