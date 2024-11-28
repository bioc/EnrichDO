# EnrichDO

***EnrichDO*** is a double weighted iterative model by integrating the DO graph topology on a global scale. It was based on the latest annotations of the human genome with DO terms, and double weighted the annotated protein-coding genes. On one hand, to reinforce the saliency of direct gene-DO annotations, different initial weights were assigned to directly annotated genes and indirectly annotated genes, respectively. On the other hand, to detect locally most significant node between the parent and its children, less significant nodes were dynamically down-weighted.

***EnrichDO*** exhibits higher accuracy that often yield more specific significant DO terms, which alleviate the over enriched problem. The input data are the protein-coding genes of the human genome, using the ENTREZID format of NCBI.

## Installation

To install this package, start R (version "4.4"), BiocManager (version "3.20") and enter:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

##Release version
BiocManager::install("EnrichDO")


## Devel version
BiocManager::install(version='devel')
BiocManager::install("EnrichDO")
```

**Note:** for other R and BiocManager versions, need to manually download the source code in the Bioconductor website for installation.

## Example

After installation, check vignettes with:

``` r
browseVigenttes("EnrichDO")
```

**Run cases** are stored in inst/scripts/EnrichDO_exampleTest.R

The **input data case** is stored at inst/extdata/Alzheimer_curated.csv

**Output example** of enrichment result is available in inst/examples/result.txt

The **thesis data** is in thesisData folder (<https://github.com/liangcheng-hrbmu/EnrichDO/tree/devel/thesisData>) and extdata_interpretation.txt explains the data source.
