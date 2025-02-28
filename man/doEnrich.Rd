% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichDO_function.R
\name{doEnrich}
\alias{doEnrich}
\title{doEnrich}
\usage{
doEnrich(
  interestGenes,
  test = c("hypergeomTest", "fisherTest", "binomTest", "chisqTest", "logoddTest"),
  method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
  m = 1,
  maxGsize = 5000,
  minGsize = 5,
  traditional = FALSE,
  delta = 0.01,
  penalize = TRUE,
  allDOTerms = FALSE
)
}
\arguments{
\item{interestGenes}{a vector of gene IDs.The interest gene sets should be protein-coding genes, using the ENTREZID format from NCBI.}

\item{test}{One of 'fisherTest','hypergeomTest','binomTest','chisqTest' and 'logoddTest' statistical model. Default is hypergeomTest.}

\item{method}{One of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY','fdr' and 'none',for P value correction.}

\item{m}{Set the maximum number of ancestor layers for ontology enrichment. Default is layer 1.}

\item{maxGsize}{indicates that doterms with more annotation genes than maxGsize are ignored, and the P value of these doterms is set to 1.}

\item{minGsize}{indicates that doterms with less annotation genes than minGsize are ignored, and the P value of these doterms is set to 1.}

\item{traditional}{a logical variable, TRUE for traditional enrichment analysis, FALSE for enrichment analysis with weights. Default is FALSE.}

\item{delta}{Set the threshold of nodes, if the p value of doterm is greater than delta, the nodes are not significant, and these nodes are not weighted.Default is 0.01.}

\item{penalize}{Logical value, used to alleviate the impact of different magnitudes of p-values, default value is TRUE. When set to FALSE, the degree of reduction in weight for non-significant nodes is decreased.}

\item{allDOTerms}{Logical value, whether to store all doterms in EnrichResult, defaults is FALSE (only significant nodes are retained).}
}
\value{
A \code{EnrichResult} instance.
}
\description{
given an array of human protein-genes with NCBI ENTREZID format, this function combines topological properties of the disease ontology structure for enrichment analysis.
}
\examples{
##Input data case
#the inputdata_demo variable stores validated protein-coding genes associated with Alzheimer's disease.
Alzheimer <- read.delim(file.path(system.file('extdata', package='EnrichDO'), 'Alzheimer_curated.csv'), header = FALSE)
inputdata_demo <- Alzheimer[,1]
##doEnrich case
#The enrichment results were obtained by using demo.data
demo.data <- c(1636,351,102,2932,3077,348,4137,54209)
demo_result <- doEnrich(interestGenes=demo.data,maxGsize = 100, minGsize=10)
}
\author{
Haixiu Yang
}
