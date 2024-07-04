#'@title drawBarGraph
#'@description The enrichment results are shown in a bar chart
#'@param enrich a data frame of enrichment result
#'@param n number of bars
#'@param delta the threshold of P value
#'@author Haixiu Yang
#'@return bar graph
#'@importFrom magrittr `%>%`
#'@importFrom dplyr mutate filter select
#'@import ggplot2
#'@export
#'@examples
#'demo.data=c(1636,351,102,2932,3077,348,4137,54209,5663,5328,23621,3416,3553)
#'sample1<-doEnrich(interestGenes=demo.data)
#'drawBarGraph(enrich=sample1,n=10,delta=0.05)

drawBarGraph <- function(enrich=enrich,n=10,delta=1e-15) {
  data <- enrich %>%
    filter(p<=delta) %>%
    select(DOID, DOTerm, p,cg.len,ig.len) %>%
    mutate(log10p=-log10(p), DO=str_c(DOID, DOTerm, sep="  "))%>%
    mutate(geneRatio=as.numeric(cg.len)/as.numeric(ig.len))
  data <- data[1:n,]
  data<-na.omit(data)
  if(dim(data)[1]<n){
    cat(paste0("\033[31m","The threshold delta is too low, only ",dim(data)[1] ," nodes are less than the threshold","\n","\033[39m"))
  }

  ggplot(data, aes(x=reorder(DO, log10p), y=geneRatio, fill=log10p)) +
    geom_bar(stat="identity") +
    coord_flip() +
    scale_x_discrete(position = "top") +
    scale_fill_gradient(low='orange', high='red')

}
#'@title drawPointGraph
#'@description The enrichment results are shown in a scatter plot
#'@param enrich a data frame of enrichment result.
#'@param n number of points.
#'@param delta the threshold of P value.
#'@author Haixiu Yang
#'@return scatter graph
#'@importFrom dplyr mutate filter select
#'@importFrom ggplot2 ggplot
#'@importFrom magrittr `%>%`
#'@export
#'@examples
#'demo.data=c(1636,351,102,2932,3077,348,4137,54209,5663,5328,23621,3416,3553)
#'sample2<-doEnrich(interestGenes=demo.data)
#'drawPointGraph(enrich=sample2,n=10,delta=0.05)
drawPointGraph <- function(enrich=enrich,n=10,delta=1e-15) {
  data <- enrich %>%
    filter(p<=delta) %>%
    select(DOID, DOTerm, p, p.adjust, cg.len, gene.len, ig.len) %>%
    mutate(log10p=-log10(p), DO=str_c(DOID, DOTerm, sep="  ")) %>%
    mutate(geneRatio=as.numeric(cg.len)/as.numeric(ig.len))
  data <- data[1:n,]
  data<-na.omit(data)
  if(dim(data)[1]<n){
    cat(paste0("\033[31m","The threshold delta is too low, only ",dim(data)[1] ," nodes are less than the threshold","\n","\033[39m"))
  }

  ggplot(data) +
    geom_point(mapping=aes(x=geneRatio, y=reorder(DO, log10p), size=cg.len, color=log10p)) +
    scale_y_discrete(position = "right") +
    scale_color_gradient(low='blue', high='red')
}
#'@title writeDoTerms
#'@description Output DOterms as text
#'@param doterms a data frame of do terms.
#'@param file the address and name of the output file.
#'@author Haixiu Yang
#'@return text
#'@importFrom dplyr mutate select
#'@importFrom readr write_delim
#'@importFrom magrittr `%>%`
#'@export
#'@examples
#'demo.data=c(1636,351,102,2932,3077,348,4137,54209,5663,5328,23621,3416,3553)
#'sample3<-doEnrich(interestGenes=demo.data)
#'writeDoTerms(sample3)
writeDoTerms <- function(doterms=doterms,file="doterms.txt") {
  data <- doterms %>%
    mutate(genes=map_chr(gene.arr, str_c, collapse=",")) %>%
    mutate(parents=map_chr(parent.arr, str_c, collapse=",")) %>%
    mutate(children=map_chr(child.arr, str_c, collapse=",")) %>%
    select(DOID, DOTerm, level, genes, parents, children, gene.len, parent.len, child.len)
  write_delim(data, file=file, delim="\t", eol="\n", quote="none")
}

#'@title writeResult
#'@description Output enrichment result as text
#'@param enrich a data frame of the enrichment result.
#'@param file the address and name of the output file.
#'@param P Output only doterm information with p values less than or equal to P.
#'@param Q Output only doterm information with p.adjust values less than or equal to Q.
#'@author Haixiu Yang
#'@return text
#'@importFrom dplyr mutate  select arrange
#'@importFrom readr write_delim
#'@importFrom magrittr `%>%`
#'@export
#'@examples
#'demo.data=c(1636,351,102,2932,3077,348,4137,54209,5663,5328,23621,3416,3553)
#'sample4<-doEnrich(interestGenes=demo.data)
#'writeResult(sample4)
writeResult <- function(enrich=enrich,file="result.txt",Q=1,P=1) {
  data <- enrich %>%
    mutate(cg=map_chr(cg.arr, str_c, collapse=","))%>%
	mutate(geneRatio=paste0(cg.len,"/",ig.len)) %>%
	mutate(bgRatio=paste0(gene.len,"/",length(dotermgenes)))%>%
    select(DOID, DOTerm, p, p.adjust,geneRatio,bgRatio, cg) %>%
    arrange(p)
  data<-dplyr::filter(data,p<=P,p.adjust<=Q)
  write_delim(data, file=file, delim="\t", eol="\n", quote="none")
}
#'@title drawGraphViz
#'@description the enrichment results are shown in a tree diagram
#'@param enrich0 a data frame of the enrichment result
#'@param n the number of most significant nodes
#'@param labelfontsize the font size of nodes
#'@param numview Displays the number of intersections between the interest set and each doterm.
#'@param pview Displays the P value for each dotrem.
#'@author Haixiu Yang
#'@return tree diagram
#'@import purrr
#'@importFrom dplyr mutate filter select arrange
#'@import Rgraphviz
#'@import graph
#'@importFrom BiocGenerics unique
#'@importFrom magrittr `%>%`
#'@importFrom methods new
#'@importFrom tidyr unnest
#'@importFrom RColorBrewer brewer.pal
#'@export
#'@examples
#'demo.data=c(1636,351,102,2932,3077,348,4137,54209,5663,5328,23621,3416,3553)
#'sample5<-doEnrich(interestGenes=demo.data)
#'drawGraphViz(sample5)
#'
#'#The p-value and the number of intersections are not visible
#'drawGraphViz(sample5,numview=FALSE,pview=FALSE)
drawGraphViz <- function(enrich0=enrich, n=10,labelfontsize=14,numview=TRUE,pview=TRUE) {

  enrich0 <- enrich0 %>% arrange(p)

  data <- enrich0[1:n, ]
  nodes <- c()
  walk(data$DOID, function(DOID) {
    nodes <<- append(nodes, DOID)
    ancestors <- getAncestors(DOID, trace=TRUE)
    nodes <<- append(nodes, ancestors)
  })
  nodes <- unique(nodes)
  rEG <- new("graphNEL",nodes=nodes, edgemode="directed")

  #Add edge
  data.extends <- enrich0 %>% filter(DOID %in% nodes) %>% select(DOID, DOTerm,p, c(parent.arr),cg.len)

  edges <- data.extends %>% unnest(cols=parent.arr)
  pwalk(edges, function(DOID, parent.arr, ...){
    rEG <<- addEdge(from=parent.arr, to=DOID, graph = rEG, weights = 1)
  })

  nAttrs <- list()
  labels <- data.extends$DOTerm
  names(labels) <- data.extends$DOID
  nAttrs$label <- labels

  pvalue<-data.extends$p
  names(pvalue) <- data.extends$DOID
  nAttrs$pvalue <- pvalue

  cglen<-data.extends$cg.len
  names(cglen)<-data.extends$DOID
  nAttrs$cglen<-cglen

  labelfontsize<-rep(labelfontsize,length(pvalue))
  names(labelfontsize) <- data.extends$DOID
  nAttrs$fontsize <- labelfontsize

  #Fill color
  colors <- brewer.pal(9,'YlOrRd')
  data.extends <- enrich0 %>%
    filter(DOID %in% nodes(rEG)) %>%
    arrange(p) %>%
    mutate(log10p=-log10(p)) %>%
    mutate(fillcolor= as.character(cut(log10p, breaks=9, label=colors)) )

  fills <- data.extends$fillcolor
  names(fills) <- data.extends$DOID
  nAttrs$fillcolor <- fills

  #Fix Size
  fixedSize <- rep(FALSE, length(nodes(rEG)))
  names(fixedSize) <- nodes(rEG)
  nAttrs$fixedSize <- fixedSize

  shapes <- map(nodes(rEG), function(DOID){
    if(DOID %in% data$DOID) return('rectangle')
    else return('ellipse')
  })
  names(shapes) <- nodes(rEG)
  nAttrs$shape <- shapes

  g1layout <- agopen(rEG, name="foo",nodeAttrs = nAttrs,attrs=list(graph=list(rankdir="TB"), node=list(fixedsize=FALSE)),)
  Rgraphviz::plot(g1layout)
  if(pview==TRUE){
    for (i in 1:length(g1layout@AgNode)) {
	pval<-nAttrs[["pvalue"]][[g1layout@AgNode[[i]]@name]]
	if(pval!=1){pval<-format(pval,digit=5,scientific=TRUE)}
    text(getX(getNodeCenter(g1layout@AgNode[[i]])), getY(getNodeCenter(g1layout@AgNode[[i]])),labels=pval,pos=1, col="red",cex = 0.5)
  }
  }
  if(numview==TRUE){
    for (i in 1:length(g1layout@AgNode)) {
    text(getX(getNodeCenter(g1layout@AgNode[[i]])), getY(getNodeCenter(g1layout@AgNode[[i]])),labels=nAttrs[["cglen"]][[g1layout@AgNode[[i]]@name]],pos=3, col="dark red",cex = 0.5)
  }
  }

}

#'@title drawHeatmap
#'@description The top DOID_n nodes in the enrichment results showed the top gene_n genes with the highest weight sum.
#'@param enrich a data frame of the enrichment result
#'@param interestGenes A collection of interest genes in vector form
#'@param DOID_n There are DOID_n nodes with the highest significance in the enrichment results.
#'@param gene_n Among the selected DOID_n nodes, the top gene_n genes with the highest weight sum are selected to show.
#'@param fontsize_row Set the font size of the gene tag.
#'@param readable Logical value that controls whether the gene tag is in symbol format
#'@param ... Other parameters in the pheatmap function also apply.
#'@author Haixiu Yang
#'@return heat map
#'@import pheatmap
#'@importFrom BiocGenerics unique
#'@importFrom purrr map2
#'@importFrom BiocGenerics intersect setdiff
#'@importFrom grDevices colorRampPalette
#'@importFrom RColorBrewer brewer.pal
#'@importFrom clusterProfiler bitr
#'@export
#'@examples
#'demo.data=c(1636,351,102,2932,3077,348,4137,54209,5663,5328,23621,3416,3553)
#'sample6<-doEnrich(interestGenes=demo.data)
#'drawHeatmap(interestGenes=demo.data,enrich = sample6,gene_n = 10,fontsize_row = 8,readable=TRUE)

drawHeatmap<-function(interestGenes,enrich=enrich,DOID_n=10,gene_n=50,fontsize_row=10,readable=TRUE,...){

  n<-DOID_n
  m<-gene_n

  #数据准备
  data<-enrich[1:n,c("DOID","gene.w")]

  allnodeGene<-unique(names(unlist(data$gene.w)))
  interestgenes<-intersect(interestGenes,allnodeGene)
  diffGene<-setdiff(interestGenes,interestgenes)
  m<-min(gene_n,length(interestgenes))

  if(length(diffGene)!=0){
  note<-paste0("\033[31m","The following genes you input do not exist in the top DOID_n nodes:","\n",paste0(diffGene,collapse = " "),"\n","\033[39m")
  cat(note)
  }

  assign("weight_matrix",NULL,pos=1)
  weight_matrix<<-matrix(0,nrow = n+1,ncol = length(interestgenes),dimnames = list(c(data$DOID,"colSum"),interestgenes))

  #DOID-gene权重矩阵构建
  map2(data$gene.w,data$DOID,function(w,id){
    gene<-w[intersect(interestgenes,names(w))]
    weight_matrix[id,names(gene)]<<-as.numeric(gene)
  })
  weight_matrix[n+1,]<<-colSums(weight_matrix)
  weightmatrix<-t(weight_matrix[-(n+1),names(sort(weight_matrix[n+1,],decreasing = TRUE)[1:m])])

  if(readable==TRUE){
    #提供基因标签展示为symbol
    entrez<-row.names(weightmatrix)
    cat(paste0("\033[31m","gene symbol conversion result: ","\033[39m"))
    symbol<-clusterProfiler::bitr(entrez,fromType ="ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
    row.names(weightmatrix)<-symbol$SYMBOL

  }


  colors<-colorRampPalette(brewer.pal(11, "RdYlBu")[6:3] )(10)

  # colors<-colorRampPalette(brewer.pal(11,"RdGy")[7:2])(10)
  # colors<-colorRampPalette(brewer.pal(9,"YlGnBu")[2:7])(10)
  #colors <- colorRampPalette(c("gray", "yellow", "orange"))(25)

  pheatmap(weightmatrix,border_color = NA, cluster_cols = FALSE,color = colors,angle_col  = 45,fontsize_row = fontsize_row,...)

}

