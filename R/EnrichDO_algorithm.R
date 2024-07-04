#'@title Enrich_internal
#'@description Internal calculation of enrichment analysis
#'@author Haixiu Yang
#'@param result_do Receives the file output by the wrireResult function, which is used to visually display the enrichment results (without running the enrichment operation again).
#'@importFrom dplyr mutate filter select
#'@importFrom stringr str_c
#'@importFrom tidyr separate
#'@import purrr
#'@importFrom BiocGenerics intersect
#'@importFrom stats p.adjust fisher.test phyper pbinom chisq.test prop.test
#'@return  A \code{enrich} instance.
#doterm结构输出
TermStruct <- function(result_do){
  assign("doterms",doterms,pos=1)
  cat("The term structure is output: doterms","\n")
  enrich_doterms<-filter(doterms,doterms$DOID%in%result_do$DOID)%>%select(DOID,parent.arr,gene.len)
  enrich_result<-separate(result_do,col = geneRatio,sep = "/",into = c("cg.len","ig.len"))
  enrich<-merge(enrich_result,enrich_doterms,by="DOID")
  enrich<-arrange(enrich,p)
  assign("enrich",enrich,pos=1)
  cat("The enrichment results you provide are stored in enrich","\n")
}
#初始化函数
init <- function(traditional) {
  assign("enrich", NULL, pos = 1)
  assign("DOID_C",NULL,pos=1)
  assign("doterms",doterms,pos=1)
  #初始化gene的权重
  if(traditional == TRUE){
	enrich <<- doterms %>%
	mutate(gene.w=map2(gene.len, gene.arr, function(n, arr){ w<-rep(1, times=n); names(w)<-arr; return(w) })) %>%
	mutate(p=1.0)
	print("Traditional test")
	}else{
	  enrich <<- doterms %>%
	  mutate(gene.w=map2(gene.arr, weight.arr, function(g, w){ names(w)<-g; return(w);} )) %>% #直接将weights赋值给gene.w
	  mutate(p=1.0)
	  print("Descending rights test")
	}
}


#递归的相似度计算函数
computeTermSig <- function(interestGenes, level, DOID, genes, weights, parents, childrenToTest,test,traditional,delta,penalize) {
  p<-Test(test,interestGenes, genes, weights)
  if(p==0){p<-.Machine$double.xmin}

  enrich[enrich$DOID==DOID,]$p <<- p#对当前DOID，只在此处进行P值更新。case1中对子进行wp更新，case2中对当前和祖先节点进行w更新

  if(length(childrenToTest)==0) return(0)

  #2023.10.30 当p值大于delta时，终止后续步骤。
  if(p>=delta) return(0)

  #2024-3-20改进w的计算方式
  sigRatios <- map_dbl(childrenToTest, function(c.DOID){ return(log(enrich[enrich$DOID==c.DOID,]$p)/log(p)) })
  names(sigRatios) <- childrenToTest
  #添加惩罚函数
  if(penalize==TRUE){
    penal<- map_dbl(childrenToTest, function(c.DOID){
      penall<-max(1/10*log10(.Machine$double.xmin)/(log10(enrich[enrich$DOID==c.DOID,]$p)+log10(p)),1)
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

  if(length(isigRatios)==0) { ## CASE1.全部子节点都没有父节点显著，将子节点相关gene的权重降低，并重新计算和记录子节点的p值

    pwalk(Penal_sigR, function(ch.DOID,ch.penal, ch.sigRatio){

      ch.weight <- enrich[enrich$DOID==ch.DOID,]$gene.w[[1]]
      ch.genes  <- enrich[enrich$DOID==ch.DOID,]$gene.arr[[1]]

      ch.weight <- ch.weight *ch.sigRatio/ch.penal#2024-1-18#2024-4-01除以penal改成减去penal

      ch.p <- Test(test,interestGenes, ch.genes, ch.weight)
      enrich[enrich$DOID==ch.DOID,]$p <<- ch.p
      enrich[enrich$DOID==ch.DOID,]$gene.w <<- list(ch.weight)

    })


  } else {
    # CASE2.部分子节点比父节点显著，则迭代每一个子节点的当前父节点，降低父节点相关gene的权重，重新计算当前父节点的每一个祖先节点权重，并对其余子节点进行分类（递归计算computeTermSig）

    Penal_sigR<-filter(Penal_sigR,sigRatios>1)#选择sigRatios>0的子节点
    pwalk(Penal_sigR, function(ch.DOID,ch.penal, ch.sigRatio){

      ch.weight <- enrich[enrich$DOID==ch.DOID,]$gene.w[[1]]
      ch.genes <- enrich[enrich$DOID==ch.DOID,]$gene.arr[[1]]

    #2024-3-20更新父节点的权重1
      same.genes<-intersect(ch.genes, genes)
      weights[same.genes]<-weights[same.genes]/(ch.sigRatio*ch.penal)#更新父节点中子节点基因的权重
      enrich[enrich$DOID==DOID,]$gene.w <<- list(weights)

      ancestors <- getAncestors(DOID)

      #2024-1-18
      ancestors<-intersect(ancestors,DOID_C)

      walk(ancestors, function(anc.DOID){

        anc.genes  <- enrich[enrich$DOID==anc.DOID,]$gene.arr[[1]]
        anc.weight <- enrich[enrich$DOID==anc.DOID,]$gene.w[[1]]

        #更新父节点weight
        same.genes <- intersect(ch.genes, anc.genes)#2024-1-18
        anc.weight[same.genes] <- anc.weight[same.genes] / (ch.sigRatio*ch.penal)

        enrich[enrich$DOID==anc.DOID,]$gene.w <<- list(anc.weight)
      })
    })
    #2024-3-20更新父节点的权重2
    weights<-enrich[enrich$DOID==DOID,]$gene.w[[1]]#更新传入computeTermSig的weights

    computeTermSig(interestGenes, level, DOID, genes, weights, parents, names(nsigRatios),test,traditional,delta,penalize)
  }
  return(0)
}

#统计模型调用
Test<-function(test,interestGenes, genes, weights){
  switch(test,fisherTest=fisherTest(interestGenes, genes, weights),
         hypergeomTest=hypergeomTest(interestGenes, genes, weights),
         binomTest=binomTest(interestGenes, genes, weights),
         chisqTest=chisqTest(interestGenes, genes, weights),
         logoddTest=logoddTest(interestGenes, genes, weights))
}

#fisher检验函数
fisherTest <- function(interestGenes, genes, weights) {
  a <- intersect(interestGenes, genes)
  b <- intersect(setdiff(dotermgenes, interestGenes), genes)
  c <- intersect(interestGenes, setdiff(dotermgenes, genes))
  d <- intersect(setdiff(dotermgenes, interestGenes), setdiff(dotermgenes, genes))

  wa <- floor(sum(weights[a]))
  wb <- floor(sum(weights[b]))
  wc <- length(c)
  wd <- length(d)

  tableR<-matrix(c(wa,wb,wc,wd),nrow = 2,ncol = 2,byrow = TRUE)
  p<-fisher.test(tableR, alternative = 'greater')$p.value
  return(p)
}
###超几何分布
hypergeomTest<-function(interestGenes, genes, weights){
  a <- intersect(interestGenes, genes)
  b <- intersect(setdiff(dotermgenes, interestGenes), genes)
  c <- intersect(interestGenes, setdiff(dotermgenes, genes))
  d <- intersect(setdiff(dotermgenes, interestGenes), setdiff(dotermgenes, genes))

  wa <- floor(sum(weights[a]))
  wb <- floor(sum(weights[b]))
  wc <- length(c)
  wd <- length(d)

  tableR<-matrix(c(wa,wb,wc,wd),nrow = 2,ncol = 2,byrow = TRUE)
  p<-phyper(wa-1,wa+wb,wc+wd,wa+wc,lower.tail = FALSE)
  return(p)
}
###二项分布
binomTest<-function(interestGenes, genes, weights){
  a <- intersect(interestGenes, genes)
  b <- intersect(setdiff(dotermgenes, interestGenes), genes)
  c <- intersect(interestGenes, setdiff(dotermgenes, genes))
  d <- intersect(setdiff(dotermgenes, interestGenes), setdiff(dotermgenes, genes))

  wa <- floor(sum(weights[a]))
  wb <- floor(sum(weights[b]))
  wc <- length(c)
  wd <- length(d)
  q<-wa-1
  size<-wa+wc
  prob<-(wa+wb)/length(dotermgenes)
  p<-pbinom(q, size, prob, lower.tail = FALSE, log.p = FALSE)
  ###q为成功的次数，size为试验次数，prob为单次试验成功的概率
  return(p)
}
###卡方分布
chisqTest<-function(interestGenes, genes, weights){
  a <- length(intersect(interestGenes, genes))
  b <- length(intersect(setdiff(dotermgenes, interestGenes), genes))###要选择出dotermgenes文件中含有的geneid
  c <- length(intersect(interestGenes, setdiff(dotermgenes, genes)))
  d <- length(intersect(setdiff(dotermgenes, interestGenes), setdiff(dotermgenes, genes)))
  wa <- floor(sum(weights[a]))
  wb <- floor(sum(weights[b]))
  wc <- length(c)
  wd <- length(d)
  tableR<-matrix(c(wa,wb,wc,wd),nrow = 2,ncol = 2,byrow = TRUE)
  p<-chisq.test(tableR)
  ###a,b,c,d四格表参数
  return(p)
}
###对数优势比LR
logoddTest<-function(interestGenes, genes, weights){
  a <- length(intersect(interestGenes, genes))
  b <- length(intersect(setdiff(dotermgenes, interestGenes), genes))###要选择出dotermgenes文件中含有的geneid
  c <- length(intersect(interestGenes, setdiff(dotermgenes, genes)))
  d <- length(intersect(setdiff(dotermgenes, interestGenes), setdiff(dotermgenes, genes)))
  wa <- floor(sum(weights[a]))
  wb <- floor(sum(weights[b]))
  wc <- length(c)
  wd <- length(d)
  p<-prop.test(c(wa,wc),c(wa+wb,wc+wd))
  return(p)
}
#获取祖先节点
getAncestors <- function(DOID, trace=FALSE) {
  ancestors <- c()
  if(DOID=='DOID:4') return(ancestors) #root节点直接返回

  parents <- enrich[enrich$DOID==DOID,]$parent.arr[[1]]
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

