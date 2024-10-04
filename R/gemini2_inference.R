#' GEMINI2 Inference function
#'
#' @description The function uses coordinate ascent method to infer a series of gene/guide effects.
#' @param input Initialized object from \code{\link[gemini2]{gemini2_initialization}} (list)
#' @param method Link function for modeling combined guide effect, "min" recommended for asymmetric combinatory screen and "prod" recommended for symmetric one. (character)
#' @param gamma 2xthe penalized factor set for double gene effect. When NULL is specificed, the function calculates it as 0.05 x (Mean guide pairs per gene pair) x quantile 50% of the initialized inverse of variance (numeric)
#' @param mu 2*the penalized factor set for single gene effect. Default is 0.1 (numeric)
#' @param nc_zero Logical flag for whether double gene effects for pairs with negative control and single gene effects for negative control genes are set to zero or not. Default is TRUE (logical)
#' @param cores The cores to be used for code parallelization. If not specified, only 2 codes will be spared out when parallelization. (numeric)
#' @param n_iterations Number of iteractions +1 to run for coordinate ascent optimization (numeric)
#' @param verbose Logical flag for whether to print all message out (logical)
#' @details
#' Link functions "max" and "sqrt" for also available for \code{method}
#'
#' @return A list of optimized guide/gene effect, inverse of variance, the calculated LFC, corresponding weights and the series of costs for optimization
#' @export
#' @examples
#' infered_D <- gemini2_inference(input=init_D,method="min",gamma=NULL,mu=0.1,nc_zero=T,cores=NULL,n_iterations=21,verbose=F)
#'
gemini2_inference <- function(input,method="min",gamma=NULL,mu=0.1,nc_zero=T,cores=NULL,n_iterations=21,verbose=F){
  #=== Initialization===
  n <- ncol(input$D)
  if(n==1){
    D <- cbind(input$D,input$D)%>%set_colnames(c(colnames(input$D),paste0("Pesudo_",colnames(input$D))))
    Xt <- input$X
    Yt <- cbind(input$Y,input$Y)%>%set_colnames(c(colnames(input$Y),paste0("Pesudo_",colnames(input$Y))))
    St <- cbind(input$S,input$S)%>%set_colnames(c(colnames(input$S),paste0("Pesudo_",colnames(input$S))))
    Tt <- (1/input$beta)*(1/2)
    Tt <- cbind(Tt,Tt)%>%set_colnames(c(colnames(Tt),paste0("Pesudo_",colnames(Tt))))
  }else{
    D <- input$D
    Xt <- input$X
    Yt <- input$Y
    St <- input$S
    Tt <- (1/input$beta)*(1/2)
  }
  guide_join <- input$guide_join
  gene_join <- input$gene_join
  a <- input$a
  nc_gene <- input$nc_gene
  xr <- c(a,1-a)
  Xmt <- Xt%>%select(-gene)%>%matrix_conv
  nc.guides <- Xt%>%filter(gene %in% nc_gene)%>%pull(guide)
  if(is.null(guide_join)){guide_join <- gsub("[ATGC]","",rownames(D)[1])}
  nc.guidePat <- paste(sapply(nc.guides,function(s) paste0(s,"\\",guide_join,"|\\",guide_join,s)),collapse = "|")
  nc_genePat <- paste(sapply(nc_gene,function(s) paste0(s,"\\",gene_join,"|\\",gene_join,s)),collapse = "|")
  Ag <- (1:nrow(Xt))[Xt$guide %in% setdiff(input$map$leftS_guide,nc.guides)]
  Pg <-  (1:nrow(Xt))[Xt$guide %in% setdiff(input$map$rightS_guide,nc.guides)]
  #=== Calculating the weights===
  ncmap <- input$map%>%
    filter(grepl(nc_genePat,pair))%>%
    mutate(Xg = gsub(nc.guidePat,"",rowname),
           gene = gsub(nc_genePat,"",pair))%>%
    select(rowname,Xg,gene)%>%
    arrange(gene)
  pairsN <- input$map%>%
    filter(!grepl(nc_genePat,pair))%>%
    group_by(pair)%>%
    summarise(n=n())%>%pull(n)%>%mean()
  if(is.null(gamma)){gamma <- 0.05*pairsN*quantile(Tt, 0.5)}
  wgil <- 1/(1+1*(abs(D[ncmap$rowname,]-Yt[ncmap$gene,]))^3)
  rownames(wgil) <- ncmap$Xg
  wgilALL <- rbind(wgil,matrix(1,nrow=length(nc.guides),ncol=ncol(wgil),dimnames=list(nc.guides,colnames(wgil))))
  gihj <- input$map
  wghl <- pmin(wgilALL[gihj$leftS_guide,],wgilALL[gihj$rightS_guide,])
  rownames(wghl) <- gihj$rowname
  #=== Implement initial cost and set indicator===
  X1 <- Xmt[input$map$leftS_guide,]
  X2 <- Xmt[input$map$rightS_guide,]
  X12 <- switch(method,"min" = pmin(X1,X2), "max" = pmax(X1,X2),"sqrt" = sqrt(X1*X2),"prod" = X1*X2)
  tosq <- D[input$map$rowname,]-Yt[input$map$leftS_gene,]*X1-Yt[input$map$rightS_gene,]*X2-St[input$map$pair,]*X12
  costll <- sum((wghl[input$map$rowname,]*0.5*log(0.5*Tt/pi))-(0.5*wghl[input$map$rowname,]*Tt*(tosq^2)))-
    sum(0.5*mu*(Yt)^2)-sum(0.5*gamma*(St)^2)
  costs <- c(costll)
  currentmax <- 1
  r <- 1
  mcc <- if(is.null(cores)){parallel::detectCores()-2}else{cores}
  ### Looping for mle implementation
  while(r < n_iterations){
    if(verbose){message("Round ", r," starts!")}
    currentmax <- costll
    if(verbose){message("Updating X ...\n")}
    gAX <- pbmcapply::pbmcmapply(function(x){
      G <- Xt$gene[x]
      GI <- Xt$guide[x]
      Map_gi <- input$map%>%filter(grepl(paste0("^",GI,"\\",guide_join,"|\\",guide_join,GI,"$"),rowname))
      Map_gi%<>%select(Dijl = rowname,Sgh = pair)%>%
        mutate(Yg = G,Xgi = GI,
               Yh = gsub(paste0(G,"\\",gene_join,"|\\",gene_join,G),"",Sgh),
               Xhj= gsub(paste0(GI,"\\",guide_join,"|\\",guide_join,GI),"",Dijl))
      Xgi <- Xmt[GI,]
      Xhj <- Xmt[Map_gi$Xhj,]
      if(method=="sqrt"){
        ### Construct for polynomial D+Ax-Cx^2-Bx^3
        za <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(D[Map_gi$Dijl,]*Yt[Map_gi$Yg,]-Yt[Map_gi$Yg,]*Yt[Map_gi$Yh,]*Xhj-0.5*(St[Map_gi$Sgh,])^2*Xhj))
        zb <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(Yt[Map_gi$Yg,])^2)
        zc <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*1.5*St[Map_gi$Sgh,]*Yt[Map_gi$Yg,]*sqrt(Xhj))
        zd <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(0.5*D[Map_gi$Dijl,]*St[Map_gi$Sgh,]*sqrt(Xhj)-0.5*St[Map_gi$Sgh,]*Yt[Map_gi$Yh,]*Xhj^(1.5)))
        ### solve for sqrt(Xgi)
        roots <- polyroot(c(zd, za, -zc,-zb))
        rr <- Re(roots)[abs(Im(roots)) < 1e-10]
        rrr <- (rr[(rr>sqrt(xr[1]) & rr<sqrt(xr[2]))])^2
        xpool <- if(length(rrr)!=0){c(0,1,rrr)}else{c(0,1)}
        xRpool <-if(length(rrr)!=0){c(xr,rrr)}else{xr}
      }else if(method=="prod"){
        A <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(D[Map_gi$Dijl,]-Yt[Map_gi$Yh,]*Xhj)*(Yt[Map_gi$Yg,]+St[Map_gi$Sgh,]*Xhj))
        B <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(Yt[Map_gi$Yg,]+St[Map_gi$Sgh,]*Xhj)^2)
        o <- A/B
        xpool <- if(o>(xr[1])&o<(xr[2])){c(0,1,o)}else{c(0,1)}
        xRpool <-if(o>(xr[1])&o<(xr[2])){c(xr,o)}else{xr}
      }else{
        SIGNeval <- switch(method,
                           "max" =0.5*(1+(sign(Xgi-Xhj)%>%ifelse(.==0,-1,.))),
                           "min" =0.5*(1-(sign(Xgi-Xhj)%>%ifelse(.==0,-1,.))))
        Xgh <- switch(method,"min" = pmin(Xgi,Xhj), "max" = pmax(Xgi,Xhj))
        A <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(D[Map_gi$Dijl,]-Yt[Map_gi$Yh,]*Xhj-St[Map_gi$Sgh,]*Xgh)*(Yt[Map_gi$Yg,]+St[Map_gi$Sgh,]*SIGNeval))
        B <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*Yt[Map_gi$Yg,]*(Yt[Map_gi$Yg,]+St[Map_gi$Sgh,]*SIGNeval))
        o <- A/B
        xpool <- if(o>(xr[1])&o<(xr[2])){c(0,1,o)}else{c(0,1)}
        xRpool <-if(o>(xr[1])&o<(xr[2])){c(xr,o)}else{xr}
      }
      C <- sapply(xpool,function(s)
      {Xgh <- switch(method,"min" = pmin(s,Xhj), "max" = pmax(s,Xhj),"sqrt" = sqrt(s*Xhj),"prod"=s*Xhj)
      ll <- sum(-0.5*wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(D[Map_gi$Dijl,]-Yt[Map_gi$Yh,]*Xhj-Yt[Map_gi$Yg,]*s-St[Map_gi$Sgh,]*Xgh)^2)
      return(ll)})
      f <- xRpool[which.max(C)]
      return(f)},
      Ag,
      mc.cores = mcc)
    gPX <- pbmcapply::pbmcmapply(function(x){
      G <- Xt$gene[x]
      GI <- Xt$guide[x]
      Map_gi <- input$map%>%filter(grepl(paste0("^",GI,"\\",guide_join,"|\\",guide_join,GI,"$"),rowname))
      Map_gi%<>%select(Dijl = rowname,Sgh = pair)%>%
        mutate(Yg = G,Xgi = GI,
               Yh = gsub(paste0(G,"\\",gene_join,"|\\",gene_join,G),"",Sgh),
               Xhj= gsub(paste0(GI,"\\",guide_join,"|\\",guide_join,GI),"",Dijl))
      Xgi <- Xmt[GI,]
      Xhj <- Xmt[Map_gi$Xhj,]
      if(method=="sqrt"){
        ### Construct for polynomial D+Ax-Cx^2-Bx^3
        za <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(D[Map_gi$Dijl,]*Yt[Map_gi$Yg,]-Yt[Map_gi$Yg,]*Yt[Map_gi$Yh,]*Xhj-0.5*(St[Map_gi$Sgh,])^2*Xhj))
        zb <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(Yt[Map_gi$Yg,])^2)
        zc <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*1.5*St[Map_gi$Sgh,]*Yt[Map_gi$Yg,]*sqrt(Xhj))
        zd <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(0.5*D[Map_gi$Dijl,]*St[Map_gi$Sgh,]*sqrt(Xhj)-0.5*St[Map_gi$Sgh,]*Yt[Map_gi$Yh,]*Xhj^(1.5)))
        ### solve for sqrt(Xgi)
        roots <- polyroot(c(zd, za, -zc,-zb))
        rr <- Re(roots)[abs(Im(roots)) < 1e-10]
        rrr <- (rr[(rr>sqrt(xr[1]) & rr<sqrt(xr[2]))])^2
        xpool <- if(length(rrr)!=0){c(0,1,rrr)}else{c(0,1)}
        xRpool <-if(length(rrr)!=0){c(xr,rrr)}else{xr}
      }else if(method=="prod"){
        A <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(D[Map_gi$Dijl,]-Yt[Map_gi$Yh,]*Xhj)*(Yt[Map_gi$Yg,]+St[Map_gi$Sgh,]*Xhj))
        B <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(Yt[Map_gi$Yg,]+St[Map_gi$Sgh,]*Xhj)^2)
        o <- A/B
        xpool <- if(o>(xr[1])&o<(xr[2])){c(0,1,o)}else{c(0,1)}
        xRpool <-if(o>(xr[1])&o<(xr[2])){c(xr,o)}else{xr}
      }else{
        SIGNeval <- switch(method,
                           "max" =0.5*(1+(sign(Xgi-Xhj)%>%ifelse(.==0,-1,.))),
                           "min" =0.5*(1-(sign(Xgi-Xhj)%>%ifelse(.==0,-1,.))))
        Xgh <- switch(method,"min" = pmin(Xgi,Xhj), "max" = pmax(Xgi,Xhj))
        A <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(D[Map_gi$Dijl,]-Yt[Map_gi$Yh,]*Xhj-St[Map_gi$Sgh,]*Xgh)*(Yt[Map_gi$Yg,]+St[Map_gi$Sgh,]*SIGNeval))
        B <- sum(wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*Yt[Map_gi$Yg,]*(Yt[Map_gi$Yg,]+St[Map_gi$Sgh,]*SIGNeval))
        o <- A/B
        xpool <- if(o>(xr[1])&o<(xr[2])){c(0,1,o)}else{c(0,1)}
        xRpool <-if(o>(xr[1])&o<(xr[2])){c(xr,o)}else{xr}
      }
      C <- sapply(xpool,function(s)
      {Xgh <- switch(method,"min" = pmin(s,Xhj), "max" = pmax(s,Xhj),"sqrt" = sqrt(s*Xhj),"prod"=s*Xhj)
      ll <- sum(-0.5*wghl[Map_gi$Dijl,]*Tt[Map_gi$Dijl,]*(D[Map_gi$Dijl,]-Yt[Map_gi$Yh,]*Xhj-Yt[Map_gi$Yg,]*s-St[Map_gi$Sgh,]*Xgh)^2)
      return(ll)})
      f <- xRpool[which.max(C)]
      return(f)},
      Pg,
      mc.cores = mcc)
    Xmt[Ag] <- gAX
    Xmt[Pg] <- gPX
    Xt$X <- Xmt
    # One round of updating Y
    if(verbose){message("Updating Y ...\n")}
    gN <- length(rownames(Yt))
    for (i in (1:gN)){
      if(verbose){if(i == round(gN/4,0)){message("25% progress completed!")}}
      if(verbose){if(i == round(gN/2,0)){message("50% progress completed!")}}
      if(verbose){if(i == round(gN*0.75,0)){message("75% progress completed!")}}
      if(verbose){if(i == gN){message("Done!")}}
      y <- rownames(Yt)[i]
      Map_yg <- input$map%>%filter(grepl(paste0("^",y,"\\",gene_join,"|\\",gene_join,y,"$"),pair))
      Map_yg%<>%
        mutate(Yg=y,
               Xgi = case_when(leftS_gene==y~leftS_guide,TRUE~rightS_guide),
               Yh = gsub(paste0(y,"\\",gene_join,"|\\",gene_join,y),"",pair),
               Xhj = case_when(leftS_gene==y~rightS_guide,TRUE~leftS_guide))%>%
        select(Dijl = rowname,Sgh = pair,Yg,Xgi,Yh,Xhj)
      Xgi <- Xmt[Map_yg$Xgi,]
      Xhj <- Xmt[Map_yg$Xhj,]
      Xgh <- switch(method,"min" = pmin(Xgi,Xhj), "max" = pmax(Xgi,Xhj),"sqrt" = sqrt(Xgi*Xhj),"prod"=Xgi*Xhj)
      A <- colSums(wghl[Map_yg$Dijl,]*Tt[Map_yg$Dijl,]*Xgi*(D[Map_yg$Dijl,]-Yt[Map_yg$Yh,]*Xhj-St[Map_yg$Sgh,]*Xgh))
      B <- colSums(wghl[Map_yg$Dijl,]*Tt[Map_yg$Dijl,]*Xgi^2)+mu
      if(nc_zero&y==nc_gene){Yt[y,] <- 0}else{Yt[y,] <- A/B}
    }
    # One round of updating S
    if(verbose){message("Updating S ...\n")}
    Srange <- if(nc_zero){rownames(St)[grep(nc_gene,rownames(St),invert=T)]}else{rownames(St)}
    gS <- pbmcapply::pbmcmapply(function(s){
      g <- strsplit(s,split=gene_join)[[1]][1]
      Map_sgh <- input$map%>%filter(pair==s)
      Map_sgh%<>%
        mutate(Yg=g,
               Xgi = case_when(leftS_gene==g~leftS_guide,TRUE~rightS_guide),
               Yh = gsub(paste0(g,"\\",gene_join,"|\\",gene_join,g),"",pair),
               Xhj = case_when(leftS_gene==g~rightS_guide,TRUE~leftS_guide))%>%
        select(Dijl = rowname,Sgh = pair,Yg,Xgi,Yh,Xhj)
      Xgi <- Xmt[Map_sgh$Xgi,]
      Xhj <- Xmt[Map_sgh$Xhj,]
      Xgh <- switch(method,"min" = pmin(Xgi,Xhj), "max" = pmax(Xgi,Xhj),"sqrt" = sqrt(Xgi*Xhj),"prod"=Xgi*Xhj)
      A <- colSums(wghl[Map_sgh$Dijl,]*Tt[Map_sgh$Dijl,]*(D[Map_sgh$Dijl,]-Yt[Map_sgh$Yg,]*Xgi-Yt[Map_sgh$Yh,]*Xhj)*Xgh)
      B <- colSums(wghl[Map_sgh$Dijl,]*Tt[Map_sgh$Dijl,]*Xgh^2)+gamma
      return(A/B)},
      Srange,
      mc.cores = mcc)%>%t()
    St[Srange,] <- gS
    # One round of updating Tau
    if(verbose){message("Updating tau ...\n")}
    gT <- pbmcapply::pbmcmapply(function(t){
      Map_tau <- input$map%>%filter(rowname==t)%>%select(Dijl=rowname,Sgh = pair,Yg=leftS_gene,Xgi=leftS_guide,
                                                         Yh=rightS_gene,Xhj=rightS_guide)
      Xgi <- Xmt[Map_tau$Xgi,]
      Xhj <- Xmt[Map_tau$Xhj,]
      ### colSums for matrix case
      Xgh <- switch(method,"min" = pmin(Xgi,Xhj), "max" = pmax(Xgi,Xhj),"sqrt" = sqrt(Xgi*Xhj),"prod"=Xgi*Xhj)
      f <- 1/((D[Map_tau$Dijl,]-Yt[Map_tau$Yg,]*Xgi-Yt[Map_tau$Yh,]*Xhj-St[Map_tau$Sgh,]*Xgh)^2+0.05)
      return(f)},
      rownames(Tt),
      mc.cores = mcc)%>%t()
    Tt[rownames(Tt),] <- gT
    # Cost evaluation
    X1 <-  Xmt[input$map$leftS_guide,]
    X2 <-  Xmt[input$map$rightS_guide,]
    X12 <- switch(method,"min" = pmin(X1,X2), "max" = pmax(X1,X2),"sqrt" = sqrt(X1*X2),"prod"=X1*X2)
    tosq <- D[input$map$rowname,]-Yt[input$map$leftS_gene,]*X1-Yt[input$map$rightS_gene,]*X2-St[input$map$pair,]*X12
    costll <- sum((wghl[input$map$rowname,]*0.5*log(0.5*Tt/pi))-(0.5*wghl[input$map$rowname,]*Tt*(tosq^2)))-
      sum(0.5*mu*(Yt)^2)-sum(0.5*gamma*(St)^2)
    message("Round ", r," has cost ",round(costll,2))
    costs <- c(costs,costll)
    r <- r+1
  }
  return(list(D=D,map = input$map,Xf=Xmt,Yf=Yt,Sf=St,Tf=Tt,cost = costs,weights = wghl))
}
