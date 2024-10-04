#' GEMINI2 Initilization function
#'
#' @description This function processes from raw screening counts to initialize
#' a list of guide effect, single/double effect, sample variance for inference
#' @param counts Raw guide pair count data obtained from NGS reads (matrix)
#' @param replicate.map Mapping between replicate name and corresponding cell line name (data.frame)
#' @param guide.annot Mapping between guide pair string,the row names of count matrix,to genes in two sides (data.frame)
#' @param nc_gene A vector or single character string for negative control genes (character)
#' @param a The noise to be added/subtracted from extreme values of guide effect. Default is 0.01 (numeric)
#' @param cores The cores to be used for code parallelization. If not specified, only 2 codes will be spared out when parallelization. (number)
#' @param gene_join The separator to join the gene pairs, default "/" (character)
#' @details
#' counts should be a matrix which rownames are guide pairs from the first column of \code{guide.annot} and colnames are sample replicate names
#' replicate.map should be a data frame with three columns in the order:
#' \itemize{
#'  \item 1st: replicate name, the names of sample replicates that has been screened, contain \code{colnames} of counts matrix
#'  \item 2nd: sample name, the cell line names corresponding to replicate names in the first column
#'  \item 3rd: label for time point, a factor (Two levels: "ETP","LTP") column specifying which replicates were screened in early time point.
#' }
#' guide.annot should be a data.frame with three columns in the order:
#' \itemize{
#'  \item 1st: guide pairs, the pattern used to join two guides here should be the same for rownames of counts
#'  \item 2nd: gene name mapping to the first/left guide in the 1st column.
#'  \item 3rd: gene name mapping to the second/right guide in the 1st column.
#' }
#'
#' @return a list of initialized variables including guide effect \code{X}, single gene effect \code{Y}, double gene effecr \code{S}
#' and replicate variance; the computed LFC, processed mapping data, the input a as well as negative control gene(s).
#' @export
#' @importFrom stats median
#' @importFrom magrittr subtract set_rownames set_colnames set_names %<>%
#' @importFrom dplyr mutate filter select bind_cols pull group_by group_map %>% distinct arrange
#' @importFrom tidyr gather separate
#' @importFrom parallel detectCores
#' @importFrom pbmcapply pbmcmapply
#' @importFrom purrr pmap_chr
#'
#' @examples
#' init_D <- gemini2_initilization(counts=COUNTS,replicate.map = repMap,guide.annot=guideAnnotation,a=0.05,nc_gene="AAVS1",cores=NULL,gene_join = "/")
#'
gemini2_initilization <- function(counts,replicate.map,guide.annot,nc_gene,a=0.01,cores=NULL,gene_join = "/"){
  message("Calculating LFC from counts and replicate.map ...")
  colnames(replicate.map) <- c("colname","samplename","TP")
  SCALE <- (sum(counts, na.rm = T) / ncol(counts))
  normed_counts <- apply(counts, 2, function(x){
    x = ((x/sum(x, na.rm = T)) * SCALE) + 32
  })
  median_normed_counts <- apply(normed_counts, 2, function(x){
    log2(x) - median(log2(x), na.rm = T)
  })
  ETP.map <- replicate.map%>%filter(TP=="ETP")
  LTP.map <- replicate.map%>%filter(TP=="LTP")
  nsample <- length(unique(LTP.map$samplename))
  stopifnot(intersect(ETP.map$samplename,LTP.map$samplename)%>%length()==0)
  n_ETP <- ETP.map%>%pull(colname)%>%unique()%>%length()
  if(n_ETP== 0){
    stop("No ETP samples identified in replicate.map")
  }else if(n_ETP == 1){
    ETP = median_normed_counts[,ETP.map$colname]
  }else if(n_ETP>1){
    ETP= median_normed_counts[,ETP.map$colname] %>%as.data.frame()%>%rowMeans(na.rm = T)
  }

  LFC <- LTP.map%>%
    group_by(samplename)%>%
    group_map(~ median_normed_counts[,.$colname]%>%
                as.data.frame(optional = T)%>%
                rowMeans(na.rm = T)%>%
                subtract(ETP))%>%
    set_names(sort(unique(LTP.map$samplename)))%>%
    bind_cols()%>%
    as.data.frame()%>%
    set_rownames(rownames(median_normed_counts))%>%
    .[,unique(LTP.map$samplename)]%>%
    as.matrix()

  message("Initializing replicate variance ...")
  sample.base <- unique(ETP.map$samplename)
  cell_lines <- unique(replicate.map$samplename)
  meanmat <- replicate.map%>%
    group_by(samplename)%>%
    group_map(~ median_normed_counts[,.$colname]%>%
                as.data.frame(optional = T)%>%
                rowMeans(na.rm = T))%>%
    set_names(sort(unique(replicate.map$samplename)))%>%
    bind_cols()%>%
    as.data.frame()%>%
    set_rownames(rownames(median_normed_counts))%>%
    .[,unique(replicate.map$samplename)]%>%as.matrix()

  sdmat <- replicate.map%>%
    group_by(samplename)%>%
    group_map(~ median_normed_counts[,.$colname]%>%
                as.data.frame(optional = T)%>%
                apply(.,1,sd,na.rm=T))%>%
    set_names(sort(unique(replicate.map$samplename)))%>%
    bind_cols()%>%
    as.data.frame()%>%
    set_rownames(rownames(median_normed_counts))%>%
    .[,unique(replicate.map$samplename)]%>%as.matrix()

  if(any(!is.na(sdmat))) {
    sdmat[, colSums(!is.na(sdmat)) == 0] = median(sdmat, na.rm = T)
    sdmat <-  apply(sdmat, 2, function(x) {
      x[is.na(x)] = median(x, na.rm = T)
      return(x)})

    beta = 0.5 * (sdmat[,!(sample.base == colnames(sdmat))] ^ 2 + sdmat[, sample.base] ^
                    2 + 2 * 1e-2) %>%
      as.matrix()  %>%
      set_colnames(cell_lines[!(sample.base == colnames(sdmat))]) %>% set_rownames(rownames(sdmat))
  }else{
    beta = matrix(
      0.2 * 0.5,
      nrow = dim(sdmat)[1],
      ncol = sum(!(sample.base == colnames(sdmat))),
      dimnames = list(rownames(sdmat), cell_lines[!(sample.base == colnames(sdmat))])
    )}

  colnames(guide.annot) <- c("rowname","leftS_gene","rightS_gene")
  mcc <- if(is.null(cores)){detectCores()-2}else{cores}
  guide_join <- gsub("[ATCG]","",guide.annot$rowname[1])
  MapInfo <- guide.annot%>%
    mutate(pair = pmap_chr(list(leftS_gene,rightS_gene),~paste0(sort(c(...)),collapse = gene_join)))%>%
    separate(col=rowname,into = c("leftS_guide","rightS_guide"),remove = F,sep="[^[:alnum:]]+")
  GG <- MapInfo%>%
    gather(key = "gene_label", value = "gene", leftS_gene, rightS_gene) %>%
    gather(key = "guide_label", value = "guide", leftS_guide, rightS_guide) %>%
    filter(gsub("\\_.*","",gene_label) == gsub("\\_.*","",guide_label))%>%
    select(gene,guide)%>%
    distinct()%>%
    arrange(gene)
  message("Initializing guide effect ...")
  Xdt <- GG%>%mutate(X = 1-a)
  message("Initializing single gene effect ...")
  Y <- pbmcmapply(function(g){
    nc.pairs <- sapply(nc_gene,function(s) paste(sort(c(g,s)),collapse = gene_join))%>%as.character()
    LFC[MapInfo%>%filter(pair %in% nc.pairs)%>%pull(rowname),]%>%as.matrix()%>%apply(.,2,median,na.rm=T)},
    unique(GG$gene),
    mc.cores = mcc)%>%as.matrix()
  if(nsample!=1){Y%<>%t()}
  Y[nc_gene,] <- 0
  message("Initializing double gene effect ...")
  S <- pbmcmapply(function(p){
    LFC[MapInfo%>%filter(pair==p)%>%pull(rowname),]%>%as.matrix()%>%apply(.,2,median,na.rm=T)},
    unique(MapInfo$pair),
    mc.cores = mcc)%>%t()-
    Y[gsub(paste0(".*\\",gene_join),"",unique(MapInfo$pair)),]-
    Y[gsub(paste0("\\",gene_join,".*"),"",unique(MapInfo$pair)),]
  if(nsample==1){S%<>%t()}
  nc_genePat <- paste(sapply(nc_gene,function(s) paste0("^",s,"\\",gene_join,
                                                        "|","\\",gene_join,s,"$")),collapse = "|")
  S[grep(nc_genePat,rownames(S),value=T),]<- 0
  if(nsample==1){
    colnames(S)=colnames(LFC)=colnames(Y)=unique(LTP.map$samplename)
  }
  return(list(
    D=LFC,X=Xdt,Y=Y,S=S,beta=beta,map = MapInfo,
    guide_join = guide_join,gene_join = gene_join,a=a,nc_gene = nc_gene))}
