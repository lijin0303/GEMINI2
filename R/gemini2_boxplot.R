#' GEMINI2 boxplot visualization function
#'
#' @description A function to visualize the results of GEMINI over the raw data.
#' @param opt The list of inference results from \code{\link[gemini2]{gemini2_inference}} (list)
#' @param gene_pair The gene pair whose LFC data is used for visualization  (character)
#' @param cell_line The cell line name for extracting LFC data (character)
#' @param nc_gene A character vector naming the genes to use as a negative control,
#' to be paired with each individual g and h in a gene pair. (character)
#' @param show_inference A logical indicating whether to show the
#' inferred individual or combined values for each gene/gene pair (logical)
#' @details Raw LFC data is plotted for each gene combination (`g`-`nc_gene`, `h`-`nc_gene`, `g`-`h`) in a standard boxplot.
#' Horizontal green line segments are plotted over the box plots indicating the individual gene effects or
#' the inferred total effect of a particular gene combination. Each guide
#' pair can be colored based on the inferred sample independent effects
#' \eqn{g_i}, \eqn{h_j}, and \eqn{g_i,h_j}. Additionally, colors and shapes
#' can be used to distinguish unique guides targeting gene g and h, respectively.
#'
#' @return A ggplot2 object
#' @importFrom dplyr case_when inner_join
#' @importFrom tibble rownames_to_column
#' @export
#' @examples
#' p_box <- gemini2_boxplot(opt=infered_D,gene_pair="HDAC1/HDAC2",cell_line="MELJUSO_SKIN",nc_gene="AAVS1",show_inference=T)
#'
gemini2_boxplot <- function(opt,gene_pair,cell_line,nc_gene,show_inference=T){
  infered_S <- opt$Sf
  infered_Y <- opt$Yf
  LFC <- opt$D
  map <- opt$map
  gene_join <- gsub("[A-Z,0-9,\\-]","",rownames(opt$Sf)[1])
  g <- gsub(paste0(gene_join,".*"),"",gene_pair)
  h <- gsub(paste0(".*",gene_join),"",gene_pair)
  nc_pairs <- expand.grid(nc_gene,c(g,h),stringsAsFactors = F)%>%
    mutate(spair = pmap_chr(list(Var1,Var2),~paste0(sort(c(...)),collapse = gene_join)))%>%
    pull(spair)

  lfc_boxplot <- data.frame(lfc_cl = LFC[map%>%filter(pair %in% c(nc_pairs,gene_pair))%>%pull(rowname),cell_line])%>%
    rownames_to_column("rowname")%>%
    inner_join(map,by="rowname")%>%
    mutate(guide_g = case_when(rightS_gene==g~rightS_guide,
                               TRUE~leftS_guide),
           guide_h = case_when(leftS_gene==h ~leftS_guide,
                               TRUE~rightS_guide))%>%
    mutate(guide_g = case_when(pair==nc_pairs[2] ~"NegCtrl sgRNA",
                               TRUE~guide_g),
           guide_h = case_when(pair==nc_pairs[1]~"NegCtrl sgRNA",
                               TRUE~guide_h))%>%
    mutate(spair = case_when(pair==nc_pairs[1]~paste0(g,gene_join,"NegCtrl"),
                             pair==nc_pairs[2]~paste0(h,gene_join,"NegCtrl"),
                             TRUE~gene_pair))%>%
    select(-leftS_gene,-rightS_gene,-leftS_guide,-rightS_guide,-pair)

  colorsChoices <- c(brewer.pal(n = 9, name = "Set1")[-6],
                     brewer.pal(n = 8, name = "Dark2"))
  shapeChoices <- c(15,17,18,7,9,13,14,11,4,10,8,3,5,12,0)

  lfc_boxplot%<>%
    mutate(guide_g = fixed_levels(guide_g,f = "NegCtrl sgRNA"),
           guide_h = fixed_levels(guide_h,f = "NegCtrl sgRNA"),
           spair  =  factor(spair,levels=c(paste0(g,gene_join,"NegCtrl"),
                                           gene_pair,
                                           paste0(h,gene_join,"NegCtrl"))))
  l1 <- levels(lfc_boxplot$guide_g)%>%.[-length(.)]
  l2 <- levels(lfc_boxplot$guide_h)%>%.[-length(.)]
  n1 <- length(l1)
  n2 <- length(l2)

  p <- ggplot(lfc_boxplot,aes(x = spair,y = lfc_cl))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(aes(color = guide_g,shape = guide_h),size=3,width = 0.2)+
    scale_color_manual(values=c(colorsChoices[1:n1],"#141824"),breaks=l1)+
    scale_fill_manual(values=c(colorsChoices[1:n1],"#141824"),breaks=l1)+
    scale_shape_manual(values=c(shapeChoices[1:n2],16),breaks=l2)+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          text = element_text(color='black'),
          panel.grid = element_blank(),
          axis.text = element_text(color='black'),
          strip.background = element_rect(colour="black", fill="gray81",
                                          size=1, linetype="solid"),
          axis.title.x = element_blank(),
          legend.key = element_rect(fill = NA)
    )+
    guides(color = guide_legend(order = 1,override.aes = list(size=3.5)),
           shape = guide_legend(order = 2,override.aes = list(size=3.5)))+
    labs(color = "guide g (color)",
         shape = "guide h (shape)",
         title = paste0(gene_pair," in ",cell_line),
         y="LFC")

  inf_vals <- c(infered_Y[g,cell_line],
                infered_S[gene_pair,cell_line]+infered_Y[g,cell_line]+infered_Y[h,cell_line],
                infered_Y[h,cell_line])
  inf_df <- data.frame(
    inf_vals  = inf_vals,
    mapped_x  = c(1:3)-0.5,
    mapped_xend = c(1:3)+0.5
  )
  if (show_inference) {
    p <- p + geom_segment(
      mapping = aes(
        x = mapped_x,
        y = inf_vals,
        yend = inf_vals,
        xend = mapped_xend),
      color="Green",
      data = inf_df )
  }
  return(p)
}
