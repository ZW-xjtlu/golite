golite
================

Installation
------------

Install the package with:

``` r
library(devtools)
install_github("ZhenWei10/golite")
```

Run GO enrichment analysis
--------------------------

The necessary input of GOEA should be the gene set gene ID, background gene ID, and an [orgDb object](http://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html).

``` r
library(golite)
library(magrittr)
library(org.Hs.eg.db)
```

Get example randomly sampled gene ids for background and gene set.

``` r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
all_eids_hg19 <- names(genes(txdb))

set.seed(1)
eids_bg <- sample(all_eids_hg19, 3500)
eids_set <- sample(eids_bg,300)
```

Run GO enrichment analysis given gene set and background.

``` r
GOEA(gene_set = eids_set,
     back_ground = eids_bg,
     orgDb = org.Hs.eg.db,
     interpret_term = T) %>% head(.,10) %>% knitr::kable(.,"markdown")
```

<table>
<colgroup>
<col width="65%" />
<col width="6%" />
<col width="6%" />
<col width="8%" />
<col width="7%" />
<col width="5%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">term</th>
<th align="right">freq_gs</th>
<th align="right">freq_bg</th>
<th align="right">p</th>
<th align="right">adj_BH</th>
<th align="right">OR</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">telomere cap complex</td>
<td align="right">6</td>
<td align="right">17</td>
<td align="right">0.0020701</td>
<td align="right">0.313615</td>
<td align="right">4.12</td>
</tr>
<tr class="even">
<td align="left">nuclear condensin complex</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">0.0023370</td>
<td align="right">0.313615</td>
<td align="right">8.75</td>
</tr>
<tr class="odd">
<td align="left">nucleotide-excision repair, DNA damage removal</td>
<td align="right">7</td>
<td align="right">25</td>
<td align="right">0.0039369</td>
<td align="right">0.313615</td>
<td align="right">3.27</td>
</tr>
<tr class="even">
<td align="left">Hrd1p ubiquitin ligase ERAD-M complex</td>
<td align="right">5</td>
<td align="right">14</td>
<td align="right">0.0046850</td>
<td align="right">0.313615</td>
<td align="right">4.17</td>
</tr>
<tr class="odd">
<td align="left">mitotic cytokinesis</td>
<td align="right">6</td>
<td align="right">20</td>
<td align="right">0.0052013</td>
<td align="right">0.313615</td>
<td align="right">3.50</td>
</tr>
<tr class="even">
<td align="left">transcriptional activator activity, metal ion regulated sequence-specific DNA binding</td>
<td align="right">3</td>
<td align="right">5</td>
<td align="right">0.0054719</td>
<td align="right">0.313615</td>
<td align="right">7.00</td>
</tr>
<tr class="odd">
<td align="left">alpha-1,6-mannosyltransferase activity</td>
<td align="right">6</td>
<td align="right">21</td>
<td align="right">0.0067680</td>
<td align="right">0.313615</td>
<td align="right">3.33</td>
</tr>
<tr class="even">
<td align="left">response to superoxide</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.0073245</td>
<td align="right">0.313615</td>
<td align="right">11.67</td>
</tr>
<tr class="odd">
<td align="left">RNA polymerase II transcription corepressor activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.0073245</td>
<td align="right">0.313615</td>
<td align="right">11.67</td>
</tr>
<tr class="even">
<td align="left">TFIIIC-class transcription factor binding</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.0073245</td>
<td align="right">0.313615</td>
<td align="right">11.67</td>
</tr>
</tbody>
</table>

The function can be vectorized, i.e. the input can be a `list` of multiple gene sets.

``` r
eids_sets <- lapply(1:10,function(x) sample(eids_bg,300)) 

GOEA(gene_set = eids_sets,
     back_ground = eids_bg,
     orgDb = org.Hs.eg.db) %>% summary
```

    ##       Length Class      Mode
    ##  [1,] 6      data.frame list
    ##  [2,] 6      data.frame list
    ##  [3,] 6      data.frame list
    ##  [4,] 6      data.frame list
    ##  [5,] 6      data.frame list
    ##  [6,] 6      data.frame list
    ##  [7,] 6      data.frame list
    ##  [8,] 6      data.frame list
    ##  [9,] 6      data.frame list
    ## [10,] 6      data.frame list

Run GO slim enrichment analysis
-------------------------------

GO slim is a subset of GO terms that can be defined at [here](http://geneontology.org/ontology/subsets/).

``` r
GOslimEA(gene_set = eids_set,
     back_ground = eids_bg,
       orgDb = org.Hs.eg.db,
     interpret_term = T)  %>% head(.,10) %>% knitr::kable(.,"markdown")
```

<table>
<colgroup>
<col width="61%" />
<col width="7%" />
<col width="7%" />
<col width="9%" />
<col width="9%" />
<col width="5%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">term</th>
<th align="right">freq_gs</th>
<th align="right">freq_bg</th>
<th align="right">p</th>
<th align="right">adj_BH</th>
<th align="right">OR</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">mitotic spindle elongation</td>
<td align="right">26</td>
<td align="right">196</td>
<td align="right">0.0134822</td>
<td align="right">0.5973330</td>
<td align="right">1.55</td>
</tr>
<tr class="even">
<td align="left">protein targeting to Golgi</td>
<td align="right">39</td>
<td align="right">330</td>
<td align="right">0.0178308</td>
<td align="right">0.5973330</td>
<td align="right">1.38</td>
</tr>
<tr class="odd">
<td align="left">regulation of cyclin-dependent protein serine/threonine kinase activity</td>
<td align="right">3</td>
<td align="right">11</td>
<td align="right">0.0610555</td>
<td align="right">0.7512072</td>
<td align="right">3.18</td>
</tr>
<tr class="even">
<td align="left">spindle pole body separation</td>
<td align="right">18</td>
<td align="right">144</td>
<td align="right">0.0615801</td>
<td align="right">0.7512072</td>
<td align="right">1.46</td>
</tr>
<tr class="odd">
<td align="left">alpha-1,3-mannosyltransferase activity</td>
<td align="right">31</td>
<td align="right">278</td>
<td align="right">0.0671522</td>
<td align="right">0.7512072</td>
<td align="right">1.30</td>
</tr>
<tr class="even">
<td align="left">mitotic S phase</td>
<td align="right">83</td>
<td align="right">848</td>
<td align="right">0.0672723</td>
<td align="right">0.7512072</td>
<td align="right">1.14</td>
</tr>
<tr class="odd">
<td align="left">mannosylphosphate transferase activity</td>
<td align="right">21</td>
<td align="right">185</td>
<td align="right">0.1042908</td>
<td align="right">0.9441091</td>
<td align="right">1.33</td>
</tr>
<tr class="even">
<td align="left">argininosuccinate metabolic process</td>
<td align="right">6</td>
<td align="right">40</td>
<td align="right">0.1217052</td>
<td align="right">0.9441091</td>
<td align="right">1.75</td>
</tr>
<tr class="odd">
<td align="left">maltose metabolic process</td>
<td align="right">7</td>
<td align="right">51</td>
<td align="right">0.1409450</td>
<td align="right">0.9441091</td>
<td align="right">1.60</td>
</tr>
<tr class="even">
<td align="left">acyl carrier activity</td>
<td align="right">79</td>
<td align="right">841</td>
<td align="right">0.1614878</td>
<td align="right">0.9441091</td>
<td align="right">1.10</td>
</tr>
</tbody>
</table>

you could set `EASE_score = TRUE` to get a more conservative p value.

For more information of EASE, please see [here](https://david.ncifcrf.gov/helps/functional_annotation.html#fisher).

``` r
GOslimEA(gene_set = eids_sets,
        back_ground = eids_bg, 
     orgDb = org.Hs.eg.db,
     EASE_Score = F) %>% lapply(.,function(x)x$p) %>% unlist %>% hist(main = "normal hypergeometric")
```

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
GOslimEA(gene_set = eids_sets,
         back_ground = eids_bg,
     orgDb = org.Hs.eg.db,
     EASE_Score = T) %>% lapply(.,function(x)x$p) %>% unlist %>% hist(main = "EASE score")
```

![](README_files/figure-markdown_github/unnamed-chunk-7-2.png)

with any questions, please contact <zhen.wei@xjtlu.edu.cn>.

``` r
sessionInfo()
```

    ## R version 3.4.2 (2017-09-28)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] zh_CN.UTF-8/zh_CN.UTF-8/zh_CN.UTF-8/C/zh_CN.UTF-8/zh_CN.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
    ##  [2] GenomicFeatures_1.30.3                 
    ##  [3] GenomicRanges_1.30.3                   
    ##  [4] GenomeInfoDb_1.14.0                    
    ##  [5] org.Hs.eg.db_3.5.0                     
    ##  [6] AnnotationDbi_1.40.0                   
    ##  [7] IRanges_2.12.0                         
    ##  [8] S4Vectors_0.16.0                       
    ##  [9] Biobase_2.38.0                         
    ## [10] BiocGenerics_0.24.0                    
    ## [11] magrittr_1.5                           
    ## [12] golite_1.0                             
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] SummarizedExperiment_1.8.1 progress_1.1.2            
    ##  [3] lattice_0.20-35            htmltools_0.3.6           
    ##  [5] rtracklayer_1.38.3         yaml_2.1.18               
    ##  [7] blob_1.1.1                 XML_3.98-1.10             
    ##  [9] rlang_0.2.0                pillar_1.2.1              
    ## [11] glue_1.2.0                 DBI_0.8                   
    ## [13] BiocParallel_1.12.0        bit64_0.9-7               
    ## [15] bindrcpp_0.2               matrixStats_0.53.1        
    ## [17] GenomeInfoDbData_1.0.0     bindr_0.1.1               
    ## [19] stringr_1.3.0              zlibbioc_1.24.0           
    ## [21] Biostrings_2.46.0          memoise_1.1.0             
    ## [23] evaluate_0.10.1            knitr_1.20                
    ## [25] biomaRt_2.34.2             highr_0.6                 
    ## [27] GSEABase_1.40.1            Rcpp_0.12.16              
    ## [29] xtable_1.8-2               backports_1.1.2           
    ## [31] DelayedArray_0.4.1         graph_1.56.0              
    ## [33] annotate_1.56.2            XVector_0.18.0            
    ## [35] bit_1.1-12                 Rsamtools_1.30.0          
    ## [37] RMySQL_0.10.14             digest_0.6.15             
    ## [39] stringi_1.1.7              dplyr_0.7.4               
    ## [41] grid_3.4.2                 rprojroot_1.3-2           
    ## [43] tools_3.4.2                bitops_1.0-6              
    ## [45] RCurl_1.95-4.10            tibble_1.4.2              
    ## [47] RSQLite_2.0                pkgconfig_2.0.1           
    ## [49] Matrix_1.2-12              prettyunits_1.0.2         
    ## [51] assertthat_0.2.0           rmarkdown_1.9             
    ## [53] httr_1.3.1                 R6_2.2.2                  
    ## [55] GenomicAlignments_1.14.2   compiler_3.4.2
