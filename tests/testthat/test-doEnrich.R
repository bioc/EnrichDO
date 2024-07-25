library(EnrichDO)

test_that("doEnrich", {
  demo.data=c(1636,351,102,2932,3077,348,4137,54209,5663,5328,23621,3416,3553)
  res <- sapply(c("fisherTest","hypergeomTest","binomTest","chisqTest","logoddTest" ), function(test){ aa<-doEnrich(interestGenes = demo.data, test = test) ;return(aa@enrich$p)})
  expect_true(all(res >= 0) && all(res <= 1))
})
