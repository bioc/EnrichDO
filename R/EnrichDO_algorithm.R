#'@title Enrich_internal
#'@description Internal calculation of enrichment analysis
#'@author Haixiu Yang
#'@param resultDO Receives the file output by the wrireResult function, which is used to visually display the enrichment results (without running the enrichment operation again).
#'@importFrom dplyr mutate filter select
#'@importFrom stringr str_c
#'@importFrom tidyr separate
#'@import purrr
#'@import hash
#'@importFrom BiocGenerics intersect
#'@importFrom stats p.adjust fisher.test phyper pbinom chisq.test prop.test
#'@return  A \code{EnrichResult} instance.
#doterm结构输出
TermStruct <- function(resultDO){
  assign("doterms",doterms,envir=.EnrichDOenv)
  #message("The term structure is output: doterms","\n")
  enrichDOterms<-filter(doterms,doterms$DOID%in%resultDO$DOID)%>%select(DOID,parent.arr,gene.len)
  enrichResult<-separate(resultDO,col = geneRatio,sep = "/",into = c("cg.len","ig.len"))
  enrich<-merge(enrichResult,enrichDOterms,by="DOID")
  enrich<-arrange(enrich,p)
  assign("enrich",enrich,pos=1)
  message("The enrichment results you provide are stored in enrich","\n")
}
#初始化函数
init <- function(traditional) {
  if (!exists(".EnrichDOenv", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".EnrichDOenv", new.env(), envir = envir)
  }#先创建新环境
  .EnrichDOenv <- get(".EnrichDOenv", envir = .GlobalEnv)

  assign("enrich",NULL,envir = .EnrichDOenv)
  assign("doidCount",NULL,envir = .EnrichDOenv)
  assign("doterms",doterms,envir = .EnrichDOenv)

  #初始化gene权重（2024-8-22用哈希表）
  #doterms<-get("doterms",envir = .EnrichDOenv)
  if(traditional == TRUE){
    enrich <- doterms %>%
      mutate(gene.w=map2(gene.len, gene.arr, function(n, arr){ w<-rep(1, times=n); names(w)<-arr; return(w) }))
    message("\t\t -- Traditional test-- ","\n")

  }else{
    enrich <- doterms %>%
      mutate(gene.w=map2(gene.arr, weight.arr, function(g, w){ names(w)<-g; return(w);} )) #直接将weights赋值给gene.w
    message("\t\t -- Descending rights test-- ","\n")
  }

  assign("enrich", enrich,envir = .EnrichDOenv)

  enrichWeight <- hash(enrich$DOID,enrich$gene.w)
  enrichPvalue <- hash(enrich$DOID,1)
  enrichgeneArr<- hash(enrich$DOID,enrich$gene.arr)
  assign("enrichWeight",enrichWeight,envir = .EnrichDOenv)
  assign("enrichPvalue",enrichPvalue,envir = .EnrichDOenv)
  assign("enrichgeneArr",enrichgeneArr,envir = .EnrichDOenv)


  return(.EnrichDOenv)
}

#递归的相似度计算函数
computeTermSig <- function(interestGenes, level, DOID,  parents, childrenToTest,test,traditional,delta,penalize) {
  enrichWeight<-get("enrichWeight",envir = .EnrichDOenv)
  enrichPvalue<-get("enrichPvalue",envir = .EnrichDOenv)
  enrichgeneArr<-get("enrichgeneArr",envir = .EnrichDOenv)

  genes<-enrichgeneArr[[DOID]]
  weights<-enrichWeight[[DOID]]
  p<-Test(test,interestGenes, genes, weights)
  if(p==0){p<-.Machine$double.xmin}

  assign(DOID,p,envir = enrichPvalue)#对当前DOID，只在此处进行P值更新。case1中对子进行wp更新，case2中对当前和祖先节点进行w更新

  if(length(childrenToTest)==0) return(0)

  #2023.10.30 当p值大于delta时，终止后续步骤。
  if(p>=delta) return(0)

  #2024-3-20改进w的计算方式
  sigRatios <- map_dbl(childrenToTest, function(c.DOID){ return(log(enrichPvalue[[c.DOID]])/log(p)) })
  names(sigRatios) <- childrenToTest
  #添加惩罚函数
  if(penalize == TRUE){
    penal<- map_dbl(childrenToTest, function(c.DOID){
      penall<-max(1/10*log10(.Machine$double.xmin)/(log10(enrichPvalue[[c.DOID]])+log10(p)),1)
      return(penall)
    })
  }else{
    penal<-map_dbl(childrenToTest, function(c.DOID){
      return(1)
    })
  }

  names(penal)<-childrenToTest
  isigRatios <- sigRatios[sigRatios > 1]#子显著
  nsigRatios <- sigRatios[sigRatios  <=1]#父显著
  #为了方便循环操作，设置Penal_sigR存储每个节点的penalP和sigRatios
  Penal_sigR<-data.frame(ch.DOID=childrenToTest,ch.penal=penal,ch.sigRatio=sigRatios)

  if(length(isigRatios) == 0) { ## CASE1.全部子节点都没有父节点显著，将子节点相关gene的权重降低，并重新计算和记录子节点的p值

    pwalk(Penal_sigR, function(ch.DOID,ch.penal, ch.sigRatio){

      ch.weight <- enrichWeight[[ch.DOID]]
      ch.genes  <- enrichgeneArr[[ch.DOID]]

      ch.weight <- ch.weight *ch.sigRatio/ch.penal#2024-1-18#2024-4-01除以penal改成减去penal

      ch.p <- Test(test,interestGenes, ch.genes, ch.weight)

      assign(ch.DOID,ch.p,envir = enrichPvalue)
      assign(ch.DOID,ch.weight,envir = enrichWeight)

    })


  } else {
    # CASE2.部分子节点比父节点显著，则迭代每一个子节点的当前父节点，降低父节点相关gene的权重，重新计算当前父节点的每一个祖先节点权重，并对其余子节点进行分类（递归计算computeTermSig）

    Penal_sigR<-filter(Penal_sigR,sigRatios>1)#选择sigRatios>0的子节点
    pwalk(Penal_sigR, function(ch.DOID,ch.penal, ch.sigRatio){

      ch.weight <- enrichWeight[[ch.DOID]]
      ch.genes  <- enrichgeneArr[[ch.DOID]]

      #2024-3-20更新父节点的权重1
      same.genes<-intersect(ch.genes, genes)
      weights[same.genes]<-weights[same.genes]/(ch.sigRatio*ch.penal)#更新父节点中子节点基因的权重

      assign(DOID,weights,envir = enrichWeight)

      ancestors <- getAncestors(DOID)################################################需要改进

      doidCount<-get("doidCount",envir = .EnrichDOenv)
      #2024-1-18
      ancestors<-intersect(ancestors,doidCount)

      walk(ancestors, function(anc.DOID){

        anc.genes  <- enrichgeneArr[[anc.DOID]]
        anc.weight <- enrichWeight[[anc.DOID]]

        #更新父节点weight
        same.genes <- intersect(ch.genes, anc.genes)#2024-1-18
        anc.weight[same.genes] <- anc.weight[same.genes] / (ch.sigRatio*ch.penal)

        assign(anc.DOID,anc.weight,envir = enrichWeight)
      })
    })

    assign("enrichWeight",enrichWeight,envir = .EnrichDOenv)
    assign("enrichPvalue",enrichPvalue,envir = .EnrichDOenv)

    computeTermSig(interestGenes, level, DOID,parents, names(nsigRatios),test,traditional,delta,penalize)
  }
  return(0)
}

#统计模型调用
Test<-function(test,interestGenes, genes, weights){

  a <- intersect(interestGenes, genes)
  b <- intersect(setdiff(dotermgenes, interestGenes), genes)
  c <- intersect(interestGenes, setdiff(dotermgenes, genes))
  d <- intersect(setdiff(dotermgenes, interestGenes), setdiff(dotermgenes, genes))

  wa <- floor(sum(weights[a]))
  wb <- floor(sum(weights[b]))
  wc <- length(c)
  wd <- length(d)

  switch(test,
         fisherTest=fisherTest(wa,wb,wc,wd),
         hypergeomTest=hypergeomTest(wa,wb,wc,wd),
         binomTest=binomTest(wa,wb,wc),
         chisqTest=chisqTest(wa,wb,wc,wd),
         logoddTest=logoddTest(wa,wb,wc,wd))
}

#fisher检验函数
fisherTest <- function(wa,wb,wc,wd) {

  tableR<-matrix(c(wa,wb,wc,wd),nrow = 2,ncol = 2,byrow = TRUE)
  p<-fisher.test(tableR, alternative = 'greater')$p.value
  return(p)
}
###超几何分布
hypergeomTest<-function(wa,wb,wc,wd){

  tableR<-matrix(c(wa,wb,wc,wd),nrow = 2,ncol = 2,byrow = TRUE)
  p<-phyper(wa-1,wa+wb,wc+wd,wa+wc,lower.tail = FALSE)
  return(p)
}
###二项分布
binomTest<-function(wa,wb,wc){

  q<-wa-1
  size<-wa+wc
  prob<-(wa+wb)/length(dotermgenes)
  p<-pbinom(q, size, prob, lower.tail = FALSE, log.p = FALSE)
  ###q为成功的次数，size为试验次数，prob为单次试验成功的概率
  return(p)
}
###卡方分布
chisqTest<-function(wa,wb,wc,wd){

  tableR<-matrix(c(wa,wb,wc,wd),nrow = 2,ncol = 2,byrow = TRUE)
  p<-chisq.test(tableR)
  p<-p[["p.value"]]
  ###a,b,c,d四格表参数
  if(is.nan(p)){p<-1}
  return(p)
}
###对数优势比LR
logoddTest<-function(wa,wb,wc,wd){
  pval<-1
  tryCatch({
    p<-prop.test(c(wa,wc),c(wa+wb,wc+wd))
    pval<<-p[["p.value"]]
  },error = function(e){
    pval<<-1
  })
  return(pval)
}
#获取祖先节点
getAncestors <- function(DOID, trace = FALSE) {
  doterms<-get("doterms",envir = .EnrichDOenv)
  ancestors <- c()
  if(DOID=='DOID:4') return(ancestors) #root节点直接返回

  parents <- doterms[doterms$DOID==DOID,]$parent.arr[[1]]
  if(length(parents)==0) return(ancestors)
  ancestors <- append(ancestors, parents)
  walk(parents, function(p){
    p.ancestors <- getAncestors(p)
    ancestors <<- append(ancestors, p.ancestors)
  })
  #2024-1-18对ancestors去重
  ancestors<-unique(ancestors)
  #debug
  if(trace) str(ancestors)
  return(ancestors)
}






