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

The necessary input of GOEA should be the gene set gene ID, background gene ID, and the [orgDb object](http://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html).

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

    ## Loading required package: GO.db

<table>
<colgroup>
<col width="5%" />
<col width="71%" />
<col width="4%" />
<col width="4%" />
<col width="5%" />
<col width="5%" />
<col width="2%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">term</th>
<th align="left">definition</th>
<th align="right">freq_gs</th>
<th align="right">freq_bg</th>
<th align="right">p</th>
<th align="right">adj_BH</th>
<th align="right">OR</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><a href="GO:0009967" class="uri">GO:0009967</a></td>
<td align="left">RNA polymerase III type 3 promoter transcriptional preinitiation complex assembly</td>
<td align="right">29</td>
<td align="right">191</td>
<td align="right">0.0012078</td>
<td align="right">0.5996853</td>
<td align="right">1.77</td>
</tr>
<tr class="even">
<td align="left"><a href="GO:0010647" class="uri">GO:0010647</a></td>
<td align="left">transcription factor activity, RNA polymerase II core promoter proximal region sequence-specific binding involved in preinitiation complex assembly</td>
<td align="right">32</td>
<td align="right">225</td>
<td align="right">0.0020695</td>
<td align="right">0.5996853</td>
<td align="right">1.66</td>
</tr>
<tr class="odd">
<td align="left"><a href="GO:0023056" class="uri">GO:0023056</a></td>
<td align="left">dopamine neurotransmitter receptor activity, coupled via Gi/Go</td>
<td align="right">32</td>
<td align="right">225</td>
<td align="right">0.0020695</td>
<td align="right">0.5996853</td>
<td align="right">1.66</td>
</tr>
<tr class="even">
<td align="left"><a href="GO:0006672" class="uri">GO:0006672</a></td>
<td align="left">bubble DNA binding</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">0.0023288</td>
<td align="right">0.5996853</td>
<td align="right">8.75</td>
</tr>
<tr class="odd">
<td align="left"><a href="GO:2001235" class="uri">GO:2001235</a></td>
<td align="left">chorismate mutase activity</td>
<td align="right">7</td>
<td align="right">24</td>
<td align="right">0.0030061</td>
<td align="right">0.5996853</td>
<td align="right">3.40</td>
</tr>
<tr class="even">
<td align="left"><a href="GO:0045859" class="uri">GO:0045859</a></td>
<td align="left">mucosal immune response</td>
<td align="right">17</td>
<td align="right">101</td>
<td align="right">0.0045006</td>
<td align="right">0.5996853</td>
<td align="right">1.96</td>
</tr>
<tr class="odd">
<td align="left"><a href="GO:0071900" class="uri">GO:0071900</a></td>
<td align="left">bundle of His development</td>
<td align="right">12</td>
<td align="right">61</td>
<td align="right">0.0045852</td>
<td align="right">0.5996853</td>
<td align="right">2.30</td>
</tr>
<tr class="even">
<td align="left"><a href="GO:0043408" class="uri">GO:0043408</a></td>
<td align="left">gene conversion of immunoglobulin genes</td>
<td align="right">15</td>
<td align="right">85</td>
<td align="right">0.0047575</td>
<td align="right">0.5996853</td>
<td align="right">2.06</td>
</tr>
<tr class="odd">
<td align="left"><a href="GO:1902531" class="uri">GO:1902531</a></td>
<td align="left">11-beta-hydroxysteroid dehydrogenase [NAD(P)] activity</td>
<td align="right">34</td>
<td align="right">256</td>
<td align="right">0.0047870</td>
<td align="right">0.5996853</td>
<td align="right">1.55</td>
</tr>
<tr class="even">
<td align="left"><a href="GO:0060538" class="uri">GO:0060538</a></td>
<td align="left">positive regulation of humoral immune response</td>
<td align="right">6</td>
<td align="right">20</td>
<td align="right">0.0051416</td>
<td align="right">0.5996853</td>
<td align="right">3.50</td>
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
GOEA(gene_set = eids_set,
     back_ground = eids_bg,
     orgDb = org.Hs.eg.db,
     interpret_term = T,
     GO_Slim = T)  %>% head(.,10) %>% knitr::kable(.,"markdown")
```

<table>
<colgroup>
<col width="8%" />
<col width="57%" />
<col width="6%" />
<col width="6%" />
<col width="8%" />
<col width="8%" />
<col width="4%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">term</th>
<th align="left">definition</th>
<th align="right">freq_gs</th>
<th align="right">freq_bg</th>
<th align="right">p</th>
<th align="right">adj_BH</th>
<th align="right">OR</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><a href="GO:0006629" class="uri">GO:0006629</a></td>
<td align="left">regulation of mitotic recombination</td>
<td align="right">25</td>
<td align="right">187</td>
<td align="right">0.0141012</td>
<td align="right">0.5162342</td>
<td align="right">1.56</td>
</tr>
<tr class="even">
<td align="left"><a href="GO:0008219" class="uri">GO:0008219</a></td>
<td align="left">transition metal ion transport</td>
<td align="right">39</td>
<td align="right">327</td>
<td align="right">0.0156435</td>
<td align="right">0.5162342</td>
<td align="right">1.39</td>
</tr>
<tr class="odd">
<td align="left"><a href="GO:0048856" class="uri">GO:0048856</a></td>
<td align="left">regulation of transcription involved in G1/S transition of mitotic cell cycle</td>
<td align="right">83</td>
<td align="right">848</td>
<td align="right">0.0683518</td>
<td align="right">0.9268498</td>
<td align="right">1.14</td>
</tr>
<tr class="even">
<td align="left"><a href="GO:0040007" class="uri">GO:0040007</a></td>
<td align="left">mitotic sister chromatid segregation</td>
<td align="right">17</td>
<td align="right">140</td>
<td align="right">0.0852517</td>
<td align="right">0.9268498</td>
<td align="right">1.42</td>
</tr>
<tr class="odd">
<td align="left"><a href="GO:0007165" class="uri">GO:0007165</a></td>
<td align="left">acyl binding</td>
<td align="right">77</td>
<td align="right">795</td>
<td align="right">0.0988005</td>
<td align="right">0.9268498</td>
<td align="right">1.13</td>
</tr>
<tr class="even">
<td align="left"><a href="GO:0021700" class="uri">GO:0021700</a></td>
<td align="left">citrulline metabolic process</td>
<td align="right">6</td>
<td align="right">40</td>
<td align="right">0.1221631</td>
<td align="right">0.9268498</td>
<td align="right">1.75</td>
</tr>
<tr class="odd">
<td align="left"><a href="GO:0006790" class="uri">GO:0006790</a></td>
<td align="left">mitotic spindle elongation</td>
<td align="right">7</td>
<td align="right">51</td>
<td align="right">0.1415087</td>
<td align="right">0.9268498</td>
<td align="right">1.60</td>
</tr>
<tr class="even">
<td align="left"><a href="GO:0043473" class="uri">GO:0043473</a></td>
<td align="left">DNA damage checkpoint</td>
<td align="right">2</td>
<td align="right">8</td>
<td align="right">0.1453356</td>
<td align="right">0.9268498</td>
<td align="right">2.92</td>
</tr>
<tr class="odd">
<td align="left"><a href="GO:0071554" class="uri">GO:0071554</a></td>
<td align="left">mitotic telophase</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">0.1641360</td>
<td align="right">0.9268498</td>
<td align="right">5.83</td>
</tr>
<tr class="even">
<td align="left"><a href="GO:0007049" class="uri">GO:0007049</a></td>
<td align="left">cell wall mannoprotein biosynthetic process</td>
<td align="right">26</td>
<td align="right">250</td>
<td align="right">0.1651823</td>
<td align="right">0.9268498</td>
<td align="right">1.21</td>
</tr>
</tbody>
</table>

you could set `EASE_score = TRUE` to get a more conservative p value.

For more information of EASE, please see [here](https://david.ncifcrf.gov/helps/functional_annotation.html#fisher).

``` r
GOEA(gene_set = eids_sets,
     back_ground = eids_bg, 
     orgDb = org.Hs.eg.db,
     GO_Slim = T,
     EASE_Score = F) %>% lapply(.,function(x)x$p) %>% unlist %>% hist(main = "normal hypergeometric")
```

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
GOEA(gene_set = eids_sets,
         back_ground = eids_bg,
     orgDb = org.Hs.eg.db,
     GO_Slim = T,
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
    ##  [1] GO.db_3.5.0                            
    ##  [2] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
    ##  [3] GenomicFeatures_1.30.3                 
    ##  [4] GenomicRanges_1.30.3                   
    ##  [5] GenomeInfoDb_1.14.0                    
    ##  [6] org.Hs.eg.db_3.5.0                     
    ##  [7] AnnotationDbi_1.40.0                   
    ##  [8] IRanges_2.12.0                         
    ##  [9] S4Vectors_0.16.0                       
    ## [10] Biobase_2.38.0                         
    ## [11] BiocGenerics_0.24.0                    
    ## [12] magrittr_1.5                           
    ## [13] golite_1.0                             
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.16               highr_0.6                 
    ##  [3] compiler_3.4.2             XVector_0.18.0            
    ##  [5] prettyunits_1.0.2          bitops_1.0-6              
    ##  [7] tools_3.4.2                zlibbioc_1.24.0           
    ##  [9] progress_1.1.2             biomaRt_2.34.2            
    ## [11] digest_0.6.15              bit_1.1-12                
    ## [13] lattice_0.20-35            RSQLite_2.0               
    ## [15] evaluate_0.10.1            memoise_1.1.0             
    ## [17] pkgconfig_2.0.1            Matrix_1.2-12             
    ## [19] DelayedArray_0.4.1         DBI_0.8                   
    ## [21] yaml_2.1.18                GenomeInfoDbData_1.0.0    
    ## [23] rtracklayer_1.38.3         httr_1.3.1                
    ## [25] stringr_1.3.0              knitr_1.20                
    ## [27] Biostrings_2.46.0          grid_3.4.2                
    ## [29] rprojroot_1.3-2            bit64_0.9-7               
    ## [31] R6_2.2.2                   BiocParallel_1.12.0       
    ## [33] XML_3.98-1.10              RMySQL_0.10.14            
    ## [35] rmarkdown_1.9              blob_1.1.1                
    ## [37] matrixStats_0.53.1         GenomicAlignments_1.14.2  
    ## [39] Rsamtools_1.30.0           backports_1.1.2           
    ## [41] htmltools_0.3.6            SummarizedExperiment_1.8.1
    ## [43] assertthat_0.2.0           stringi_1.1.7             
    ## [45] RCurl_1.95-4.10
