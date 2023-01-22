R Notebook
================

``` bash
wget -i url_data
```

``` r
#charger les librairies 
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(phyloseq)
library(DECIPHER)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

``` r
library(phangorn)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

``` r
library(ggplot2)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
library(shiny)
library(miniUI)
#library(caret)
library(pls)
```

    ## 
    ## Attaching package: 'pls'

    ## The following object is masked from 'package:ape':
    ## 
    ##     mvr

    ## The following object is masked from 'package:stats':
    ## 
    ##     loadings

``` r
library(e1071)
library(ggplot2)
library(randomForest)
```

    ## randomForest 4.7-1.1

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:randomForest':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggrepel)
#library(nlme)
library(devtools)
```

    ## Loading required package: usethis

``` r
library(reshape2)
library(PMA)
#library(structSSI)
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
library(ggnetwork)
library(intergraph)
library(scales)
library(phyloseqGraphTest)
library(genefilter)
library(impute)
```

``` r
#la variable (path) indique le chemin qui permettra d’accéder aux objets dont nous allons avoir besoin.
path <- "/home/rstudio/CC3/data"
```

``` r
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz.1", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz.1", full.names = TRUE))
#créer une variable qui extrait le nom de cet échantillon.
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
```

``` r
#Cette étape permet de visualiser la qualité des séquences grâce au Q score associé à chaque nucléotide
plotQualityProfile(fnFs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](data_import_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](data_import_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
# créer des objets (filtFs et filtRs) pour stoquer les séquences filtrées.
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

``` r
#Dans cette étapes Nous utiliserons les paramètres de filtrage standard
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(180,200),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)
```

    ##                           reads.in reads.out
    ## SRR19667519_R1.fastq.gz.1   382439    209875
    ## SRR19667520_R1.fastq.gz.1   247892    138532
    ## SRR19667521_R1.fastq.gz.1   203290    165327
    ## SRR19667522_R1.fastq.gz.1   343027    284831
    ## SRR19667523_R1.fastq.gz.1   282482    236135
    ## SRR19667524_R1.fastq.gz.1   249004    207388

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 143741700 total bases in 798565 reads from 4 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 102746800 total bases in 513734 reads from 3 samples will be used for learning the error rates.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](data_import_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 209875 reads in 97418 unique sequences.
    ## Sample 2 - 138532 reads in 69297 unique sequences.
    ## Sample 3 - 165327 reads in 65260 unique sequences.
    ## Sample 4 - 284831 reads in 95316 unique sequences.
    ## Sample 5 - 236135 reads in 90206 unique sequences.
    ## Sample 6 - 207388 reads in 81713 unique sequences.
    ## Sample 7 - 148563 reads in 52344 unique sequences.
    ## Sample 8 - 131851 reads in 50188 unique sequences.
    ## Sample 9 - 166764 reads in 58162 unique sequences.
    ## Sample 10 - 487805 reads in 152369 unique sequences.
    ## Sample 11 - 70423 reads in 30333 unique sequences.
    ## Sample 12 - 56999 reads in 25143 unique sequences.
    ## Sample 13 - 757795 reads in 243645 unique sequences.
    ## Sample 14 - 185854 reads in 66696 unique sequences.
    ## Sample 15 - 372732 reads in 147419 unique sequences.
    ## Sample 16 - 135186 reads in 65061 unique sequences.
    ## Sample 17 - 168638 reads in 83190 unique sequences.
    ## Sample 18 - 105820 reads in 63907 unique sequences.
    ## Sample 19 - 194056 reads in 90553 unique sequences.
    ## Sample 20 - 139264 reads in 65046 unique sequences.
    ## Sample 21 - 119050 reads in 58593 unique sequences.
    ## Sample 22 - 215072 reads in 91219 unique sequences.
    ## Sample 23 - 171858 reads in 64026 unique sequences.
    ## Sample 24 - 93218 reads in 33147 unique sequences.
    ## Sample 25 - 112729 reads in 55818 unique sequences.
    ## Sample 26 - 1346687 reads in 461515 unique sequences.
    ## Sample 27 - 220389 reads in 96892 unique sequences.
    ## Sample 28 - 121758 reads in 62060 unique sequences.
    ## Sample 29 - 227891 reads in 110651 unique sequences.
    ## Sample 30 - 241995 reads in 105903 unique sequences.
    ## Sample 31 - 362218 reads in 145968 unique sequences.
    ## Sample 32 - 107614 reads in 57454 unique sequences.
    ## Sample 33 - 93996 reads in 43235 unique sequences.
    ## Sample 34 - 230923 reads in 106832 unique sequences.
    ## Sample 35 - 100718 reads in 52124 unique sequences.
    ## Sample 36 - 186241 reads in 66490 unique sequences.
    ## Sample 37 - 196783 reads in 87431 unique sequences.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 209875 reads in 119169 unique sequences.
    ## Sample 2 - 138532 reads in 81487 unique sequences.
    ## Sample 3 - 165327 reads in 81137 unique sequences.
    ## Sample 4 - 284831 reads in 135391 unique sequences.
    ## Sample 5 - 236135 reads in 122132 unique sequences.
    ## Sample 6 - 207388 reads in 100022 unique sequences.
    ## Sample 7 - 148563 reads in 66167 unique sequences.
    ## Sample 8 - 131851 reads in 61549 unique sequences.
    ## Sample 9 - 166764 reads in 81360 unique sequences.
    ## Sample 10 - 487805 reads in 189690 unique sequences.
    ## Sample 11 - 70423 reads in 38185 unique sequences.
    ## Sample 12 - 56999 reads in 30876 unique sequences.
    ## Sample 13 - 757795 reads in 334839 unique sequences.
    ## Sample 14 - 185854 reads in 79520 unique sequences.
    ## Sample 15 - 372732 reads in 204045 unique sequences.
    ## Sample 16 - 135186 reads in 82762 unique sequences.
    ## Sample 17 - 168638 reads in 104623 unique sequences.
    ## Sample 18 - 105820 reads in 77020 unique sequences.
    ## Sample 19 - 194056 reads in 111066 unique sequences.
    ## Sample 20 - 139264 reads in 79613 unique sequences.
    ## Sample 21 - 119050 reads in 75948 unique sequences.
    ## Sample 22 - 215072 reads in 113746 unique sequences.
    ## Sample 23 - 171858 reads in 82506 unique sequences.
    ## Sample 24 - 93218 reads in 42416 unique sequences.
    ## Sample 25 - 112729 reads in 67304 unique sequences.
    ## Sample 26 - 1346687 reads in 789223 unique sequences.
    ## Sample 27 - 220389 reads in 159097 unique sequences.
    ## Sample 28 - 121758 reads in 105561 unique sequences.
    ## Sample 29 - 227891 reads in 162636 unique sequences.
    ## Sample 30 - 241995 reads in 171721 unique sequences.
    ## Sample 31 - 362218 reads in 279195 unique sequences.
    ## Sample 32 - 107614 reads in 69331 unique sequences.
    ## Sample 33 - 93996 reads in 54555 unique sequences.
    ## Sample 34 - 230923 reads in 134447 unique sequences.
    ## Sample 35 - 100718 reads in 61237 unique sequences.
    ## Sample 36 - 186241 reads in 90245 unique sequences.
    ## Sample 37 - 196783 reads in 99328 unique sequences.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1073 sequence variants were inferred from 97418 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
dadaRs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1080 sequence variants were inferred from 119169 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 2545 paired-reads (in 278 unique pairings) successfully merged out of 202090 (in 62355 pairings) input.

    ## 1868 paired-reads (in 212 unique pairings) successfully merged out of 132388 (in 43801 pairings) input.

    ## 1158 paired-reads (in 200 unique pairings) successfully merged out of 157085 (in 62022 pairings) input.

    ## 860 paired-reads (in 248 unique pairings) successfully merged out of 273650 (in 92741 pairings) input.

    ## 696 paired-reads (in 217 unique pairings) successfully merged out of 224401 (in 93431 pairings) input.

    ## 228 paired-reads (in 74 unique pairings) successfully merged out of 195908 (in 62064 pairings) input.

    ## 377 paired-reads (in 116 unique pairings) successfully merged out of 142186 (in 32096 pairings) input.

    ## 133 paired-reads (in 47 unique pairings) successfully merged out of 123550 (in 35881 pairings) input.

    ## 331 paired-reads (in 122 unique pairings) successfully merged out of 160197 (in 43749 pairings) input.

    ## 828 paired-reads (in 263 unique pairings) successfully merged out of 468936 (in 124325 pairings) input.

    ## 36 paired-reads (in 19 unique pairings) successfully merged out of 66796 (in 17112 pairings) input.

    ## 82 paired-reads (in 23 unique pairings) successfully merged out of 53575 (in 15531 pairings) input.

    ## 1189 paired-reads (in 373 unique pairings) successfully merged out of 730424 (in 211399 pairings) input.

    ## 46 paired-reads (in 21 unique pairings) successfully merged out of 177471 (in 53967 pairings) input.

    ## 1009 paired-reads (in 220 unique pairings) successfully merged out of 353302 (in 137536 pairings) input.

    ## 769 paired-reads (in 96 unique pairings) successfully merged out of 124325 (in 58349 pairings) input.

    ## 541 paired-reads (in 121 unique pairings) successfully merged out of 157775 (in 74694 pairings) input.

    ## 1159 paired-reads (in 140 unique pairings) successfully merged out of 96041 (in 45322 pairings) input.

    ## 834 paired-reads (in 193 unique pairings) successfully merged out of 180823 (in 81295 pairings) input.

    ## 2055 paired-reads (in 224 unique pairings) successfully merged out of 129970 (in 66822 pairings) input.

    ## 1716 paired-reads (in 215 unique pairings) successfully merged out of 109662 (in 45800 pairings) input.

    ## 3058 paired-reads (in 445 unique pairings) successfully merged out of 202195 (in 82078 pairings) input.

    ## 287 paired-reads (in 73 unique pairings) successfully merged out of 163918 (in 48228 pairings) input.

    ## 97 paired-reads (in 32 unique pairings) successfully merged out of 88487 (in 21687 pairings) input.

    ## 121 paired-reads (in 26 unique pairings) successfully merged out of 102986 (in 44142 pairings) input.

    ## 6369 paired-reads (in 1096 unique pairings) successfully merged out of 1307613 (in 420723 pairings) input.

    ## 10639 paired-reads (in 205 unique pairings) successfully merged out of 211256 (in 60985 pairings) input.

    ## 548 paired-reads (in 63 unique pairings) successfully merged out of 113154 (in 38112 pairings) input.

    ## 1309 paired-reads (in 163 unique pairings) successfully merged out of 215772 (in 73572 pairings) input.

    ## 882 paired-reads (in 143 unique pairings) successfully merged out of 231900 (in 66633 pairings) input.

    ## 1757 paired-reads (in 212 unique pairings) successfully merged out of 348553 (in 104248 pairings) input.

    ## 7751 paired-reads (in 313 unique pairings) successfully merged out of 102503 (in 29738 pairings) input.

    ## 1293 paired-reads (in 123 unique pairings) successfully merged out of 89787 (in 22755 pairings) input.

    ## 2683 paired-reads (in 299 unique pairings) successfully merged out of 221674 (in 67300 pairings) input.

    ## 927 paired-reads (in 130 unique pairings) successfully merged out of 95471 (in 32240 pairings) input.

    ## 122 paired-reads (in 21 unique pairings) successfully merged out of 175912 (in 59160 pairings) input.

    ## 296 paired-reads (in 47 unique pairings) successfully merged out of 182654 (in 74336 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                                                                                                   sequence
    ## 147                                                                                                                                               TGTATAAGAGACAGCCTACGGGTGGCTGCATTCCAAACTGGAAGTCTAGAGTGCAGGAGAGGAGAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGAGATTAGGAAGAACACCAGTGGCGAAGGCGACTCTCTGGACTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTAGTAGTC
    ## 177 CCTACGGGTGGCTGCAGTGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAAGGTTAAGTCAGTGGTCAAACTGCGGAGCTCAACTCCGTATCGCCATTGAAACTGGTCTTCTTGAGTGAGTGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCACTCAACTGACGCTGAAGCACGAAAGCGTGGGTATCGAACAGGATTAGATACCCTAGTAGTC
    ## 291                                                                                                                    CCTACGGGTGGCTGCAAAGGTCTGCGGTGAAAGACCGAAGCTAAACTTCGGTAAGCCGTGGAAACCATGCAGCTGGAGTGCAGCAGAGGATCGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGACGATCTGGGCTGTAACTGACGCTCAGTCCCGAAAGCGTGGGGAGCAAATAGGATTAGATACCCTAGTAGTC
    ## 374                                             CCTACGGGTGGCTGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTTCCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGATAGTGGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAACACCAGTGGCGAAGGCGACTATCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTAGTAGTC
    ## 399                                             CCTACGGGTGGCTGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTAGTAGTC
    ## 407                                             CCTACGGGTGGCTGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGATTTACTGGGCGTAAAGCGCGCGTAGGTGGCCAATTAAGTCAAATGTGAAATCCCCGAGCTTAACTTGGGAATTGCATTCGATACTGGTTGGCTAGAGTATGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAATACTGACACTGAGGTGCGAAAGCATGGGGAGCAAACAGGATTAGATACCCTAGTAGTC
    ##     abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 147       125     302     796    180         0      0      1   TRUE
    ## 177       111      97      47     38         0      0      2   TRUE
    ## 291        78     282     208    153         0      0      2   TRUE
    ## 374        66     148       3     82         0      0      2   TRUE
    ## 399        63     185       1     82         0      0      2   TRUE
    ## 407        62     184       2     82         0      0      2   TRUE

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]   37 4570

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  200  203  204  205  207  211  212  217  220  221  224  225  226  227  228  229 
    ##   36    1    2    1    1    5    1    1    1    2    8    2   36   29   23   96 
    ##  231  232  233  234  235  236  237  238  240  243  245  247  249  251  252  253 
    ##   70    2    2    2    3   10   22    1    2    1    1    2    1    1    1    3 
    ##  259  265  266  267  268  269  270  271  272  274  275  276  280  281  284  290 
    ##    1    2    1   10   34    3    3    2    6   45    7    2    2    1    3    1 
    ##  293  294  295  296  297  298  299  300  301  302  304  305  308  310  313  314 
    ##    1    8    8   12  924 2562  137   58   11    1    2    6    2    1    7    1 
    ##  315  317  318  319  320  321  322  323  324  325  327  328  330  331  333  334 
    ##    1    1    1    1   10    2    3    3    1    2    3    2    2    3    7   12 
    ##  336  337  338  339  340  341  342  343  344  346  347  350  351  352  353  354 
    ##    3    1   21    4    1    2   46   52    1    1    1    2    2    4   21   11 
    ##  355  356  357  358  359  361  362  363  364  365  366  367  368 
    ##   22   40    2    1    3    9    2    8    2    7    7    4    3

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 2179 bimeras out of 4570 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]   37 2391

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.8055619

``` r
#Tableau de suivi,voir combien de séquences ont été éliminées à chaque étape
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##              input filtered denoisedF denoisedR merged nonchim
    ## SRR19667519 382439   209875    207143    204509   2545    1966
    ## SRR19667520 247892   138532    136352    134304   1868    1510
    ## SRR19667521 203290   165327    161622    160528   1158     854
    ## SRR19667522 343027   284831    279719    278433    860     633
    ## SRR19667523 282482   236135    230959    229288    696     489
    ## SRR19667524 249004   207388    203130    199717    228     201

``` bash
#attribuer une taxonomie aux variants de séquence
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1
```

    ## --2023-01-18 10:07:06--  https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1
    ## Resolving zenodo.org (zenodo.org)... 188.185.124.72
    ## Connecting to zenodo.org (zenodo.org)|188.185.124.72|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 137283333 (131M) [application/octet-stream]
    ## Saving to: ‘silva_nr99_v138.1_train_set.fa.gz?download=1.6’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 11.7M 11s
    ##     50K .......... .......... .......... .......... ..........  0% 8.51M 13s
    ##    100K .......... .......... .......... .......... ..........  0% 12.4M 12s
    ##    150K .......... .......... .......... .......... ..........  0% 12.4M 12s
    ##    200K .......... .......... .......... .......... ..........  0% 15.6M 11s
    ##    250K .......... .......... .......... .......... ..........  0% 12.5M 11s
    ##    300K .......... .......... .......... .......... ..........  0% 16.0M 11s
    ##    350K .......... .......... .......... .......... ..........  0% 8.61M 11s
    ##    400K .......... .......... .......... .......... ..........  0% 15.7M 11s
    ##    450K .......... .......... .......... .......... ..........  0% 14.8M 11s
    ##    500K .......... .......... .......... .......... ..........  0% 16.0M 10s
    ##    550K .......... .......... .......... .......... ..........  0% 10.5M 11s
    ##    600K .......... .......... .......... .......... ..........  0% 15.2M 10s
    ##    650K .......... .......... .......... .......... ..........  0% 16.4M 10s
    ##    700K .......... .......... .......... .......... ..........  0% 13.0M 10s
    ##    750K .......... .......... .......... .......... ..........  0% 13.6M 10s
    ##    800K .......... .......... .......... .......... ..........  0% 14.5M 10s
    ##    850K .......... .......... .......... .......... ..........  0% 16.3M 10s
    ##    900K .......... .......... .......... .......... ..........  0% 12.0M 10s
    ##    950K .......... .......... .......... .......... ..........  0% 13.8M 10s
    ##   1000K .......... .......... .......... .......... ..........  0% 15.8M 10s
    ##   1050K .......... .......... .......... .......... ..........  0% 14.2M 10s
    ##   1100K .......... .......... .......... .......... ..........  0% 15.4M 10s
    ##   1150K .......... .......... .......... .......... ..........  0% 14.7M 10s
    ##   1200K .......... .......... .......... .......... ..........  0% 11.7M 10s
    ##   1250K .......... .......... .......... .......... ..........  0% 18.7M 10s
    ##   1300K .......... .......... .......... .......... ..........  1% 12.7M 10s
    ##   1350K .......... .......... .......... .......... ..........  1% 14.3M 10s
    ##   1400K .......... .......... .......... .......... ..........  1% 17.3M 10s
    ##   1450K .......... .......... .......... .......... ..........  1% 18.8M 10s
    ##   1500K .......... .......... .......... .......... ..........  1% 16.7M 9s
    ##   1550K .......... .......... .......... .......... ..........  1% 12.3M 9s
    ##   1600K .......... .......... .......... .......... ..........  1% 14.6M 9s
    ##   1650K .......... .......... .......... .......... ..........  1% 10.9M 10s
    ##   1700K .......... .......... .......... .......... ..........  1% 14.9M 10s
    ##   1750K .......... .......... .......... .......... ..........  1% 16.5M 9s
    ##   1800K .......... .......... .......... .......... ..........  1% 19.2M 9s
    ##   1850K .......... .......... .......... .......... ..........  1% 14.7M 9s
    ##   1900K .......... .......... .......... .......... ..........  1% 16.3M 9s
    ##   1950K .......... .......... .......... .......... ..........  1% 11.4M 9s
    ##   2000K .......... .......... .......... .......... ..........  1% 12.5M 9s
    ##   2050K .......... .......... .......... .......... ..........  1% 13.5M 9s
    ##   2100K .......... .......... .......... .......... ..........  1% 18.3M 9s
    ##   2150K .......... .......... .......... .......... ..........  1% 15.3M 9s
    ##   2200K .......... .......... .......... .......... ..........  1% 19.4M 9s
    ##   2250K .......... .......... .......... .......... ..........  1% 15.9M 9s
    ##   2300K .......... .......... .......... .......... ..........  1% 16.0M 9s
    ##   2350K .......... .......... .......... .......... ..........  1% 13.7M 9s
    ##   2400K .......... .......... .......... .......... ..........  1% 14.2M 9s
    ##   2450K .......... .......... .......... .......... ..........  1% 19.4M 9s
    ##   2500K .......... .......... .......... .......... ..........  1% 13.4M 9s
    ##   2550K .......... .......... .......... .......... ..........  1% 15.3M 9s
    ##   2600K .......... .......... .......... .......... ..........  1% 13.6M 9s
    ##   2650K .......... .......... .......... .......... ..........  2% 17.3M 9s
    ##   2700K .......... .......... .......... .......... ..........  2% 15.2M 9s
    ##   2750K .......... .......... .......... .......... ..........  2% 13.1M 9s
    ##   2800K .......... .......... .......... .......... ..........  2% 18.0M 9s
    ##   2850K .......... .......... .......... .......... ..........  2% 18.6M 9s
    ##   2900K .......... .......... .......... .......... ..........  2% 14.5M 9s
    ##   2950K .......... .......... .......... .......... ..........  2% 13.5M 9s
    ##   3000K .......... .......... .......... .......... ..........  2% 16.2M 9s
    ##   3050K .......... .......... .......... .......... ..........  2% 15.1M 9s
    ##   3100K .......... .......... .......... .......... ..........  2% 15.7M 9s
    ##   3150K .......... .......... .......... .......... ..........  2% 16.5M 9s
    ##   3200K .......... .......... .......... .......... ..........  2% 17.0M 9s
    ##   3250K .......... .......... .......... .......... ..........  2% 18.9M 9s
    ##   3300K .......... .......... .......... .......... ..........  2% 14.7M 9s
    ##   3350K .......... .......... .......... .......... ..........  2% 24.2M 9s
    ##   3400K .......... .......... .......... .......... ..........  2% 15.3M 9s
    ##   3450K .......... .......... .......... .......... ..........  2% 19.6M 9s
    ##   3500K .......... .......... .......... .......... ..........  2% 10.7M 9s
    ##   3550K .......... .......... .......... .......... ..........  2% 16.3M 9s
    ##   3600K .......... .......... .......... .......... ..........  2% 44.6M 9s
    ##   3650K .......... .......... .......... .......... ..........  2% 15.3M 9s
    ##   3700K .......... .......... .......... .......... ..........  2% 19.6M 9s
    ##   3750K .......... .......... .......... .......... ..........  2% 15.5M 9s
    ##   3800K .......... .......... .......... .......... ..........  2% 28.5M 9s
    ##   3850K .......... .......... .......... .......... ..........  2% 17.8M 9s
    ##   3900K .......... .......... .......... .......... ..........  2% 15.4M 9s
    ##   3950K .......... .......... .......... .......... ..........  2% 10.9M 9s
    ##   4000K .......... .......... .......... .......... ..........  3% 15.5M 9s
    ##   4050K .......... .......... .......... .......... ..........  3% 17.5M 9s
    ##   4100K .......... .......... .......... .......... ..........  3% 16.8M 9s
    ##   4150K .......... .......... .......... .......... ..........  3% 14.7M 9s
    ##   4200K .......... .......... .......... .......... ..........  3% 31.1M 8s
    ##   4250K .......... .......... .......... .......... ..........  3% 13.8M 8s
    ##   4300K .......... .......... .......... .......... ..........  3% 15.2M 8s
    ##   4350K .......... .......... .......... .......... ..........  3% 15.4M 8s
    ##   4400K .......... .......... .......... .......... ..........  3% 21.8M 8s
    ##   4450K .......... .......... .......... .......... ..........  3% 15.4M 8s
    ##   4500K .......... .......... .......... .......... ..........  3% 27.6M 8s
    ##   4550K .......... .......... .......... .......... ..........  3% 14.1M 8s
    ##   4600K .......... .......... .......... .......... ..........  3% 17.2M 8s
    ##   4650K .......... .......... .......... .......... ..........  3% 16.7M 8s
    ##   4700K .......... .......... .......... .......... ..........  3% 25.0M 8s
    ##   4750K .......... .......... .......... .......... ..........  3% 11.1M 8s
    ##   4800K .......... .......... .......... .......... ..........  3% 36.3M 8s
    ##   4850K .......... .......... .......... .......... ..........  3% 14.6M 8s
    ##   4900K .......... .......... .......... .......... ..........  3% 16.1M 8s
    ##   4950K .......... .......... .......... .......... ..........  3% 20.0M 8s
    ##   5000K .......... .......... .......... .......... ..........  3% 24.8M 8s
    ##   5050K .......... .......... .......... .......... ..........  3% 16.3M 8s
    ##   5100K .......... .......... .......... .......... ..........  3% 15.8M 8s
    ##   5150K .......... .......... .......... .......... ..........  3% 16.5M 8s
    ##   5200K .......... .......... .......... .......... ..........  3% 17.3M 8s
    ##   5250K .......... .......... .......... .......... ..........  3% 13.8M 8s
    ##   5300K .......... .......... .......... .......... ..........  3% 45.9M 8s
    ##   5350K .......... .......... .......... .......... ..........  4% 14.1M 8s
    ##   5400K .......... .......... .......... .......... ..........  4% 15.4M 8s
    ##   5450K .......... .......... .......... .......... ..........  4% 12.4M 8s
    ##   5500K .......... .......... .......... .......... ..........  4% 42.3M 8s
    ##   5550K .......... .......... .......... .......... ..........  4% 13.6M 8s
    ##   5600K .......... .......... .......... .......... ..........  4% 14.8M 8s
    ##   5650K .......... .......... .......... .......... ..........  4% 15.5M 8s
    ##   5700K .......... .......... .......... .......... ..........  4% 34.3M 8s
    ##   5750K .......... .......... .......... .......... ..........  4% 18.1M 8s
    ##   5800K .......... .......... .......... .......... ..........  4% 15.9M 8s
    ##   5850K .......... .......... .......... .......... ..........  4% 15.2M 8s
    ##   5900K .......... .......... .......... .......... ..........  4% 16.4M 8s
    ##   5950K .......... .......... .......... .......... ..........  4% 15.6M 8s
    ##   6000K .......... .......... .......... .......... ..........  4% 14.3M 8s
    ##   6050K .......... .......... .......... .......... ..........  4% 40.9M 8s
    ##   6100K .......... .......... .......... .......... ..........  4% 18.4M 8s
    ##   6150K .......... .......... .......... .......... ..........  4% 18.7M 8s
    ##   6200K .......... .......... .......... .......... ..........  4% 29.7M 8s
    ##   6250K .......... .......... .......... .......... ..........  4% 14.0M 8s
    ##   6300K .......... .......... .......... .......... ..........  4% 16.4M 8s
    ##   6350K .......... .......... .......... .......... ..........  4% 15.4M 8s
    ##   6400K .......... .......... .......... .......... ..........  4% 26.1M 8s
    ##   6450K .......... .......... .......... .......... ..........  4% 17.9M 8s
    ##   6500K .......... .......... .......... .......... ..........  4% 17.2M 8s
    ##   6550K .......... .......... .......... .......... ..........  4% 15.6M 8s
    ##   6600K .......... .......... .......... .......... ..........  4% 20.0M 8s
    ##   6650K .......... .......... .......... .......... ..........  4% 15.3M 8s
    ##   6700K .......... .......... .......... .......... ..........  5% 41.4M 8s
    ##   6750K .......... .......... .......... .......... ..........  5% 14.4M 8s
    ##   6800K .......... .......... .......... .......... ..........  5% 21.4M 8s
    ##   6850K .......... .......... .......... .......... ..........  5% 13.0M 8s
    ##   6900K .......... .......... .......... .......... ..........  5% 57.5M 8s
    ##   6950K .......... .......... .......... .......... ..........  5% 14.2M 8s
    ##   7000K .......... .......... .......... .......... ..........  5% 15.1M 8s
    ##   7050K .......... .......... .......... .......... ..........  5% 16.3M 8s
    ##   7100K .......... .......... .......... .......... ..........  5% 16.3M 8s
    ##   7150K .......... .......... .......... .......... ..........  5% 15.5M 8s
    ##   7200K .......... .......... .......... .......... ..........  5% 15.5M 8s
    ##   7250K .......... .......... .......... .......... ..........  5% 24.2M 8s
    ##   7300K .......... .......... .......... .......... ..........  5% 16.3M 8s
    ##   7350K .......... .......... .......... .......... ..........  5% 23.4M 8s
    ##   7400K .......... .......... .......... .......... ..........  5% 16.8M 8s
    ##   7450K .......... .......... .......... .......... ..........  5% 23.8M 8s
    ##   7500K .......... .......... .......... .......... ..........  5% 15.2M 8s
    ##   7550K .......... .......... .......... .......... ..........  5% 15.2M 8s
    ##   7600K .......... .......... .......... .......... ..........  5% 10.7M 8s
    ##   7650K .......... .......... .......... .......... ..........  5% 56.1M 8s
    ##   7700K .......... .......... .......... .......... ..........  5% 16.3M 8s
    ##   7750K .......... .......... .......... .......... ..........  5% 17.2M 8s
    ##   7800K .......... .......... .......... .......... ..........  5% 19.3M 8s
    ##   7850K .......... .......... .......... .......... ..........  5% 29.6M 8s
    ##   7900K .......... .......... .......... .......... ..........  5% 16.7M 8s
    ##   7950K .......... .......... .......... .......... ..........  5% 16.3M 8s
    ##   8000K .......... .......... .......... .......... ..........  6% 28.8M 8s
    ##   8050K .......... .......... .......... .......... ..........  6% 16.7M 8s
    ##   8100K .......... .......... .......... .......... ..........  6% 14.9M 8s
    ##   8150K .......... .......... .......... .......... ..........  6% 14.4M 8s
    ##   8200K .......... .......... .......... .......... ..........  6% 43.1M 8s
    ##   8250K .......... .......... .......... .......... ..........  6% 19.5M 8s
    ##   8300K .......... .......... .......... .......... ..........  6% 14.3M 8s
    ##   8350K .......... .......... .......... .......... ..........  6% 37.9M 8s
    ##   8400K .......... .......... .......... .......... ..........  6% 13.3M 8s
    ##   8450K .......... .......... .......... .......... ..........  6% 16.6M 8s
    ##   8500K .......... .......... .......... .......... ..........  6% 19.3M 8s
    ##   8550K .......... .......... .......... .......... ..........  6% 26.6M 8s
    ##   8600K .......... .......... .......... .......... ..........  6% 14.4M 8s
    ##   8650K .......... .......... .......... .......... ..........  6% 30.5M 8s
    ##   8700K .......... .......... .......... .......... ..........  6% 16.1M 8s
    ##   8750K .......... .......... .......... .......... ..........  6% 15.3M 8s
    ##   8800K .......... .......... .......... .......... ..........  6% 17.9M 7s
    ##   8850K .......... .......... .......... .......... ..........  6% 31.5M 7s
    ##   8900K .......... .......... .......... .......... ..........  6% 13.8M 7s
    ##   8950K .......... .......... .......... .......... ..........  6% 22.8M 7s
    ##   9000K .......... .......... .......... .......... ..........  6% 23.7M 7s
    ##   9050K .......... .......... .......... .......... ..........  6% 18.8M 7s
    ##   9100K .......... .......... .......... .......... ..........  6% 20.0M 7s
    ##   9150K .......... .......... .......... .......... ..........  6% 13.6M 7s
    ##   9200K .......... .......... .......... .......... ..........  6% 15.9M 7s
    ##   9250K .......... .......... .......... .......... ..........  6% 13.8M 7s
    ##   9300K .......... .......... .......... .......... ..........  6% 42.0M 7s
    ##   9350K .......... .......... .......... .......... ..........  7% 21.5M 7s
    ##   9400K .......... .......... .......... .......... ..........  7% 13.8M 7s
    ##   9450K .......... .......... .......... .......... ..........  7% 41.6M 7s
    ##   9500K .......... .......... .......... .......... ..........  7% 15.7M 7s
    ##   9550K .......... .......... .......... .......... ..........  7% 16.5M 7s
    ##   9600K .......... .......... .......... .......... ..........  7% 17.4M 7s
    ##   9650K .......... .......... .......... .......... ..........  7% 20.3M 7s
    ##   9700K .......... .......... .......... .......... ..........  7% 14.5M 7s
    ##   9750K .......... .......... .......... .......... ..........  7% 17.4M 7s
    ##   9800K .......... .......... .......... .......... ..........  7% 16.6M 7s
    ##   9850K .......... .......... .......... .......... ..........  7% 37.4M 7s
    ##   9900K .......... .......... .......... .......... ..........  7% 17.5M 7s
    ##   9950K .......... .......... .......... .......... ..........  7% 9.20M 7s
    ##  10000K .......... .......... .......... .......... ..........  7% 43.9M 7s
    ##  10050K .......... .......... .......... .......... ..........  7% 16.6M 7s
    ##  10100K .......... .......... .......... .......... ..........  7% 14.8M 7s
    ##  10150K .......... .......... .......... .......... ..........  7% 43.8M 7s
    ##  10200K .......... .......... .......... .......... ..........  7% 14.2M 7s
    ##  10250K .......... .......... .......... .......... ..........  7% 18.1M 7s
    ##  10300K .......... .......... .......... .......... ..........  7% 49.7M 7s
    ##  10350K .......... .......... .......... .......... ..........  7% 13.9M 7s
    ##  10400K .......... .......... .......... .......... ..........  7% 13.1M 7s
    ##  10450K .......... .......... .......... .......... ..........  7% 16.1M 7s
    ##  10500K .......... .......... .......... .......... ..........  7% 40.4M 7s
    ##  10550K .......... .......... .......... .......... ..........  7% 13.0M 7s
    ##  10600K .......... .......... .......... .......... ..........  7% 16.5M 7s
    ##  10650K .......... .......... .......... .......... ..........  7% 14.3M 7s
    ##  10700K .......... .......... .......... .......... ..........  8% 39.6M 7s
    ##  10750K .......... .......... .......... .......... ..........  8% 15.2M 7s
    ##  10800K .......... .......... .......... .......... ..........  8% 22.4M 7s
    ##  10850K .......... .......... .......... .......... ..........  8% 38.4M 7s
    ##  10900K .......... .......... .......... .......... ..........  8% 16.1M 7s
    ##  10950K .......... .......... .......... .......... ..........  8% 16.5M 7s
    ##  11000K .......... .......... .......... .......... ..........  8% 24.1M 7s
    ##  11050K .......... .......... .......... .......... ..........  8% 20.2M 7s
    ##  11100K .......... .......... .......... .......... ..........  8% 13.7M 7s
    ##  11150K .......... .......... .......... .......... ..........  8% 16.8M 7s
    ##  11200K .......... .......... .......... .......... ..........  8% 37.0M 7s
    ##  11250K .......... .......... .......... .......... ..........  8% 14.7M 7s
    ##  11300K .......... .......... .......... .......... ..........  8% 20.9M 7s
    ##  11350K .......... .......... .......... .......... ..........  8% 16.0M 7s
    ##  11400K .......... .......... .......... .......... ..........  8% 28.4M 7s
    ##  11450K .......... .......... .......... .......... ..........  8% 21.4M 7s
    ##  11500K .......... .......... .......... .......... ..........  8% 13.5M 7s
    ##  11550K .......... .......... .......... .......... ..........  8% 15.4M 7s
    ##  11600K .......... .......... .......... .......... ..........  8% 45.7M 7s
    ##  11650K .......... .......... .......... .......... ..........  8% 15.6M 7s
    ##  11700K .......... .......... .......... .......... ..........  8% 15.6M 7s
    ##  11750K .......... .......... .......... .......... ..........  8% 17.8M 7s
    ##  11800K .......... .......... .......... .......... ..........  8% 27.2M 7s
    ##  11850K .......... .......... .......... .......... ..........  8% 16.0M 7s
    ##  11900K .......... .......... .......... .......... ..........  8% 21.8M 7s
    ##  11950K .......... .......... .......... .......... ..........  8% 15.8M 7s
    ##  12000K .......... .......... .......... .......... ..........  8% 18.6M 7s
    ##  12050K .......... .......... .......... .......... ..........  9% 16.2M 7s
    ##  12100K .......... .......... .......... .......... ..........  9% 33.5M 7s
    ##  12150K .......... .......... .......... .......... ..........  9% 15.2M 7s
    ##  12200K .......... .......... .......... .......... ..........  9% 14.4M 7s
    ##  12250K .......... .......... .......... .......... ..........  9% 42.3M 7s
    ##  12300K .......... .......... .......... .......... ..........  9% 14.6M 7s
    ##  12350K .......... .......... .......... .......... ..........  9% 14.2M 7s
    ##  12400K .......... .......... .......... .......... ..........  9% 15.7M 7s
    ##  12450K .......... .......... .......... .......... ..........  9% 18.0M 7s
    ##  12500K .......... .......... .......... .......... ..........  9% 37.8M 7s
    ##  12550K .......... .......... .......... .......... ..........  9% 14.2M 7s
    ##  12600K .......... .......... .......... .......... ..........  9% 14.6M 7s
    ##  12650K .......... .......... .......... .......... ..........  9% 25.7M 7s
    ##  12700K .......... .......... .......... .......... ..........  9% 30.2M 7s
    ##  12750K .......... .......... .......... .......... ..........  9% 16.0M 7s
    ##  12800K .......... .......... .......... .......... ..........  9% 15.7M 7s
    ##  12850K .......... .......... .......... .......... ..........  9% 17.3M 7s
    ##  12900K .......... .......... .......... .......... ..........  9% 58.3M 7s
    ##  12950K .......... .......... .......... .......... ..........  9% 13.8M 7s
    ##  13000K .......... .......... .......... .......... ..........  9% 18.5M 7s
    ##  13050K .......... .......... .......... .......... ..........  9% 40.8M 7s
    ##  13100K .......... .......... .......... .......... ..........  9% 14.9M 7s
    ##  13150K .......... .......... .......... .......... ..........  9% 16.9M 7s
    ##  13200K .......... .......... .......... .......... ..........  9% 15.9M 7s
    ##  13250K .......... .......... .......... .......... ..........  9% 47.5M 7s
    ##  13300K .......... .......... .......... .......... ..........  9% 14.6M 7s
    ##  13350K .......... .......... .......... .......... ..........  9% 18.4M 7s
    ##  13400K .......... .......... .......... .......... .......... 10% 33.6M 7s
    ##  13450K .......... .......... .......... .......... .......... 10% 13.4M 7s
    ##  13500K .......... .......... .......... .......... .......... 10% 15.5M 7s
    ##  13550K .......... .......... .......... .......... .......... 10% 44.8M 7s
    ##  13600K .......... .......... .......... .......... .......... 10% 16.0M 7s
    ##  13650K .......... .......... .......... .......... .......... 10% 15.8M 7s
    ##  13700K .......... .......... .......... .......... .......... 10% 17.0M 7s
    ##  13750K .......... .......... .......... .......... .......... 10% 36.6M 7s
    ##  13800K .......... .......... .......... .......... .......... 10% 14.5M 7s
    ##  13850K .......... .......... .......... .......... .......... 10% 47.9M 7s
    ##  13900K .......... .......... .......... .......... .......... 10% 18.1M 7s
    ##  13950K .......... .......... .......... .......... .......... 10% 9.22M 7s
    ##  14000K .......... .......... .......... .......... .......... 10% 14.8M 7s
    ##  14050K .......... .......... .......... .......... .......... 10% 16.6M 7s
    ##  14100K .......... .......... .......... .......... .......... 10% 19.0M 7s
    ##  14150K .......... .......... .......... .......... .......... 10% 15.8M 7s
    ##  14200K .......... .......... .......... .......... .......... 10% 18.5M 7s
    ##  14250K .......... .......... .......... .......... .......... 10% 17.2M 7s
    ##  14300K .......... .......... .......... .......... .......... 10% 17.6M 7s
    ##  14350K .......... .......... .......... .......... .......... 10% 15.6M 7s
    ##  14400K .......... .......... .......... .......... .......... 10% 18.8M 7s
    ##  14450K .......... .......... .......... .......... .......... 10% 17.8M 7s
    ##  14500K .......... .......... .......... .......... .......... 10% 17.6M 7s
    ##  14550K .......... .......... .......... .......... .......... 10% 14.7M 7s
    ##  14600K .......... .......... .......... .......... .......... 10% 18.8M 7s
    ##  14650K .......... .......... .......... .......... .......... 10% 17.8M 7s
    ##  14700K .......... .......... .......... .......... .......... 11% 17.3M 7s
    ##  14750K .......... .......... .......... .......... .......... 11% 16.0M 7s
    ##  14800K .......... .......... .......... .......... .......... 11% 17.6M 7s
    ##  14850K .......... .......... .......... .......... .......... 11% 18.2M 7s
    ##  14900K .......... .......... .......... .......... .......... 11% 16.3M 7s
    ##  14950K .......... .......... .......... .......... .......... 11% 14.9M 7s
    ##  15000K .......... .......... .......... .......... .......... 11% 19.4M 7s
    ##  15050K .......... .......... .......... .......... .......... 11% 18.0M 7s
    ##  15100K .......... .......... .......... .......... .......... 11% 17.1M 7s
    ##  15150K .......... .......... .......... .......... .......... 11% 14.1M 7s
    ##  15200K .......... .......... .......... .......... .......... 11% 18.4M 7s
    ##  15250K .......... .......... .......... .......... .......... 11% 18.0M 7s
    ##  15300K .......... .......... .......... .......... .......... 11% 17.9M 7s
    ##  15350K .......... .......... .......... .......... .......... 11% 15.6M 7s
    ##  15400K .......... .......... .......... .......... .......... 11% 18.0M 7s
    ##  15450K .......... .......... .......... .......... .......... 11% 19.1M 7s
    ##  15500K .......... .......... .......... .......... .......... 11% 18.5M 7s
    ##  15550K .......... .......... .......... .......... .......... 11% 15.8M 7s
    ##  15600K .......... .......... .......... .......... .......... 11% 17.6M 7s
    ##  15650K .......... .......... .......... .......... .......... 11% 19.1M 7s
    ##  15700K .......... .......... .......... .......... .......... 11% 18.7M 7s
    ##  15750K .......... .......... .......... .......... .......... 11% 16.2M 7s
    ##  15800K .......... .......... .......... .......... .......... 11% 16.4M 7s
    ##  15850K .......... .......... .......... .......... .......... 11% 18.1M 7s
    ##  15900K .......... .......... .......... .......... .......... 11% 18.1M 7s
    ##  15950K .......... .......... .......... .......... .......... 11% 16.9M 7s
    ##  16000K .......... .......... .......... .......... .......... 11% 18.8M 7s
    ##  16050K .......... .......... .......... .......... .......... 12% 17.7M 7s
    ##  16100K .......... .......... .......... .......... .......... 12% 19.2M 7s
    ##  16150K .......... .......... .......... .......... .......... 12% 15.6M 7s
    ##  16200K .......... .......... .......... .......... .......... 12% 17.6M 7s
    ##  16250K .......... .......... .......... .......... .......... 12% 19.2M 7s
    ##  16300K .......... .......... .......... .......... .......... 12% 17.9M 7s
    ##  16350K .......... .......... .......... .......... .......... 12% 15.3M 7s
    ##  16400K .......... .......... .......... .......... .......... 12% 17.5M 7s
    ##  16450K .......... .......... .......... .......... .......... 12% 17.8M 7s
    ##  16500K .......... .......... .......... .......... .......... 12% 18.7M 7s
    ##  16550K .......... .......... .......... .......... .......... 12% 16.8M 7s
    ##  16600K .......... .......... .......... .......... .......... 12% 17.9M 7s
    ##  16650K .......... .......... .......... .......... .......... 12% 17.7M 7s
    ##  16700K .......... .......... .......... .......... .......... 12% 18.6M 7s
    ##  16750K .......... .......... .......... .......... .......... 12% 14.6M 7s
    ##  16800K .......... .......... .......... .......... .......... 12% 16.6M 7s
    ##  16850K .......... .......... .......... .......... .......... 12% 16.9M 7s
    ##  16900K .......... .......... .......... .......... .......... 12% 17.3M 7s
    ##  16950K .......... .......... .......... .......... .......... 12% 15.8M 7s
    ##  17000K .......... .......... .......... .......... .......... 12% 17.5M 7s
    ##  17050K .......... .......... .......... .......... .......... 12% 17.7M 7s
    ##  17100K .......... .......... .......... .......... .......... 12% 18.6M 7s
    ##  17150K .......... .......... .......... .......... .......... 12% 15.5M 7s
    ##  17200K .......... .......... .......... .......... .......... 12% 17.5M 7s
    ##  17250K .......... .......... .......... .......... .......... 12% 18.6M 7s
    ##  17300K .......... .......... .......... .......... .......... 12% 19.0M 7s
    ##  17350K .......... .......... .......... .......... .......... 12% 15.7M 7s
    ##  17400K .......... .......... .......... .......... .......... 13% 18.2M 7s
    ##  17450K .......... .......... .......... .......... .......... 13% 18.6M 7s
    ##  17500K .......... .......... .......... .......... .......... 13% 17.6M 7s
    ##  17550K .......... .......... .......... .......... .......... 13% 15.4M 7s
    ##  17600K .......... .......... .......... .......... .......... 13% 18.4M 7s
    ##  17650K .......... .......... .......... .......... .......... 13% 18.3M 7s
    ##  17700K .......... .......... .......... .......... .......... 13% 17.7M 7s
    ##  17750K .......... .......... .......... .......... .......... 13% 16.1M 7s
    ##  17800K .......... .......... .......... .......... .......... 13% 18.5M 7s
    ##  17850K .......... .......... .......... .......... .......... 13% 18.2M 7s
    ##  17900K .......... .......... .......... .......... .......... 13% 18.0M 7s
    ##  17950K .......... .......... .......... .......... .......... 13% 16.0M 7s
    ##  18000K .......... .......... .......... .......... .......... 13% 18.6M 7s
    ##  18050K .......... .......... .......... .......... .......... 13% 18.2M 7s
    ##  18100K .......... .......... .......... .......... .......... 13% 16.8M 7s
    ##  18150K .......... .......... .......... .......... .......... 13% 16.0M 7s
    ##  18200K .......... .......... .......... .......... .......... 13% 18.6M 7s
    ##  18250K .......... .......... .......... .......... .......... 13% 18.3M 7s
    ##  18300K .......... .......... .......... .......... .......... 13% 17.8M 7s
    ##  18350K .......... .......... .......... .......... .......... 13% 14.9M 7s
    ##  18400K .......... .......... .......... .......... .......... 13% 17.8M 7s
    ##  18450K .......... .......... .......... .......... .......... 13% 18.7M 7s
    ##  18500K .......... .......... .......... .......... .......... 13% 17.7M 7s
    ##  18550K .......... .......... .......... .......... .......... 13% 16.0M 7s
    ##  18600K .......... .......... .......... .......... .......... 13% 18.2M 7s
    ##  18650K .......... .......... .......... .......... .......... 13% 17.8M 7s
    ##  18700K .......... .......... .......... .......... .......... 13% 15.7M 7s
    ##  18750K .......... .......... .......... .......... .......... 14% 15.3M 7s
    ##  18800K .......... .......... .......... .......... .......... 14% 17.9M 7s
    ##  18850K .......... .......... .......... .......... .......... 14% 18.7M 7s
    ##  18900K .......... .......... .......... .......... .......... 14% 17.1M 7s
    ##  18950K .......... .......... .......... .......... .......... 14% 16.2M 7s
    ##  19000K .......... .......... .......... .......... .......... 14% 16.3M 7s
    ##  19050K .......... .......... .......... .......... .......... 14% 18.3M 7s
    ##  19100K .......... .......... .......... .......... .......... 14% 12.6M 7s
    ##  19150K .......... .......... .......... .......... .......... 14% 18.8M 7s
    ##  19200K .......... .......... .......... .......... .......... 14% 20.0M 7s
    ##  19250K .......... .......... .......... .......... .......... 14% 22.1M 7s
    ##  19300K .......... .......... .......... .......... .......... 14% 20.7M 7s
    ##  19350K .......... .......... .......... .......... .......... 14% 16.2M 7s
    ##  19400K .......... .......... .......... .......... .......... 14% 21.5M 7s
    ##  19450K .......... .......... .......... .......... .......... 14% 21.5M 7s
    ##  19500K .......... .......... .......... .......... .......... 14% 20.7M 7s
    ##  19550K .......... .......... .......... .......... .......... 14% 18.3M 7s
    ##  19600K .......... .......... .......... .......... .......... 14% 20.5M 7s
    ##  19650K .......... .......... .......... .......... .......... 14% 21.0M 6s
    ##  19700K .......... .......... .......... .......... .......... 14% 22.3M 6s
    ##  19750K .......... .......... .......... .......... .......... 14% 18.5M 6s
    ##  19800K .......... .......... .......... .......... .......... 14% 21.5M 6s
    ##  19850K .......... .......... .......... .......... .......... 14% 23.1M 6s
    ##  19900K .......... .......... .......... .......... .......... 14% 22.0M 6s
    ##  19950K .......... .......... .......... .......... .......... 14% 17.8M 6s
    ##  20000K .......... .......... .......... .......... .......... 14% 21.0M 6s
    ##  20050K .......... .......... .......... .......... .......... 14% 21.1M 6s
    ##  20100K .......... .......... .......... .......... .......... 15% 21.4M 6s
    ##  20150K .......... .......... .......... .......... .......... 15% 19.0M 6s
    ##  20200K .......... .......... .......... .......... .......... 15% 22.4M 6s
    ##  20250K .......... .......... .......... .......... .......... 15% 23.5M 6s
    ##  20300K .......... .......... .......... .......... .......... 15% 21.8M 6s
    ##  20350K .......... .......... .......... .......... .......... 15% 19.0M 6s
    ##  20400K .......... .......... .......... .......... .......... 15% 22.3M 6s
    ##  20450K .......... .......... .......... .......... .......... 15% 21.5M 6s
    ##  20500K .......... .......... .......... .......... .......... 15% 21.7M 6s
    ##  20550K .......... .......... .......... .......... .......... 15% 19.6M 6s
    ##  20600K .......... .......... .......... .......... .......... 15% 22.3M 6s
    ##  20650K .......... .......... .......... .......... .......... 15% 21.8M 6s
    ##  20700K .......... .......... .......... .......... .......... 15% 21.2M 6s
    ##  20750K .......... .......... .......... .......... .......... 15% 17.5M 6s
    ##  20800K .......... .......... .......... .......... .......... 15% 22.2M 6s
    ##  20850K .......... .......... .......... .......... .......... 15% 22.2M 6s
    ##  20900K .......... .......... .......... .......... .......... 15% 22.5M 6s
    ##  20950K .......... .......... .......... .......... .......... 15% 19.8M 6s
    ##  21000K .......... .......... .......... .......... .......... 15% 21.9M 6s
    ##  21050K .......... .......... .......... .......... .......... 15% 21.0M 6s
    ##  21100K .......... .......... .......... .......... .......... 15% 20.6M 6s
    ##  21150K .......... .......... .......... .......... .......... 15% 18.8M 6s
    ##  21200K .......... .......... .......... .......... .......... 15% 22.6M 6s
    ##  21250K .......... .......... .......... .......... .......... 15% 23.1M 6s
    ##  21300K .......... .......... .......... .......... .......... 15% 21.5M 6s
    ##  21350K .......... .......... .......... .......... .......... 15% 20.1M 6s
    ##  21400K .......... .......... .......... .......... .......... 15% 21.9M 6s
    ##  21450K .......... .......... .......... .......... .......... 16% 21.6M 6s
    ##  21500K .......... .......... .......... .......... .......... 16% 22.4M 6s
    ##  21550K .......... .......... .......... .......... .......... 16% 18.2M 6s
    ##  21600K .......... .......... .......... .......... .......... 16% 19.7M 6s
    ##  21650K .......... .......... .......... .......... .......... 16% 20.9M 6s
    ##  21700K .......... .......... .......... .......... .......... 16% 17.9M 6s
    ##  21750K .......... .......... .......... .......... .......... 16% 17.4M 6s
    ##  21800K .......... .......... .......... .......... .......... 16% 20.8M 6s
    ##  21850K .......... .......... .......... .......... .......... 16% 22.2M 6s
    ##  21900K .......... .......... .......... .......... .......... 16% 21.1M 6s
    ##  21950K .......... .......... .......... .......... .......... 16% 15.3M 6s
    ##  22000K .......... .......... .......... .......... .......... 16% 15.4M 6s
    ##  22050K .......... .......... .......... .......... .......... 16% 21.6M 6s
    ##  22100K .......... .......... .......... .......... .......... 16% 21.6M 6s
    ##  22150K .......... .......... .......... .......... .......... 16% 15.1M 6s
    ##  22200K .......... .......... .......... .......... .......... 16% 22.5M 6s
    ##  22250K .......... .......... .......... .......... .......... 16% 19.5M 6s
    ##  22300K .......... .......... .......... .......... .......... 16% 20.5M 6s
    ##  22350K .......... .......... .......... .......... .......... 16% 13.5M 6s
    ##  22400K .......... .......... .......... .......... .......... 16% 20.6M 6s
    ##  22450K .......... .......... .......... .......... .......... 16% 21.0M 6s
    ##  22500K .......... .......... .......... .......... .......... 16% 22.4M 6s
    ##  22550K .......... .......... .......... .......... .......... 16% 19.0M 6s
    ##  22600K .......... .......... .......... .......... .......... 16% 23.2M 6s
    ##  22650K .......... .......... .......... .......... .......... 16% 21.4M 6s
    ##  22700K .......... .......... .......... .......... .......... 16% 22.1M 6s
    ##  22750K .......... .......... .......... .......... .......... 17% 11.9M 6s
    ##  22800K .......... .......... .......... .......... .......... 17% 20.5M 6s
    ##  22850K .......... .......... .......... .......... .......... 17% 21.6M 6s
    ##  22900K .......... .......... .......... .......... .......... 17% 21.5M 6s
    ##  22950K .......... .......... .......... .......... .......... 17% 18.1M 6s
    ##  23000K .......... .......... .......... .......... .......... 17% 22.7M 6s
    ##  23050K .......... .......... .......... .......... .......... 17% 20.3M 6s
    ##  23100K .......... .......... .......... .......... .......... 17% 22.7M 6s
    ##  23150K .......... .......... .......... .......... .......... 17% 19.3M 6s
    ##  23200K .......... .......... .......... .......... .......... 17% 21.1M 6s
    ##  23250K .......... .......... .......... .......... .......... 17% 21.1M 6s
    ##  23300K .......... .......... .......... .......... .......... 17% 14.5M 6s
    ##  23350K .......... .......... .......... .......... .......... 17% 18.3M 6s
    ##  23400K .......... .......... .......... .......... .......... 17% 22.0M 6s
    ##  23450K .......... .......... .......... .......... .......... 17% 16.7M 6s
    ##  23500K .......... .......... .......... .......... .......... 17% 39.0M 6s
    ##  23550K .......... .......... .......... .......... .......... 17% 16.0M 6s
    ##  23600K .......... .......... .......... .......... .......... 17% 26.2M 6s
    ##  23650K .......... .......... .......... .......... .......... 17% 23.9M 6s
    ##  23700K .......... .......... .......... .......... .......... 17% 17.1M 6s
    ##  23750K .......... .......... .......... .......... .......... 17% 16.9M 6s
    ##  23800K .......... .......... .......... .......... .......... 17% 43.2M 6s
    ##  23850K .......... .......... .......... .......... .......... 17% 19.3M 6s
    ##  23900K .......... .......... .......... .......... .......... 17% 15.0M 6s
    ##  23950K .......... .......... .......... .......... .......... 17% 32.4M 6s
    ##  24000K .......... .......... .......... .......... .......... 17% 15.4M 6s
    ##  24050K .......... .......... .......... .......... .......... 17% 15.2M 6s
    ##  24100K .......... .......... .......... .......... .......... 18% 34.7M 6s
    ##  24150K .......... .......... .......... .......... .......... 18% 13.9M 6s
    ##  24200K .......... .......... .......... .......... .......... 18% 51.3M 6s
    ##  24250K .......... .......... .......... .......... .......... 18% 14.3M 6s
    ##  24300K .......... .......... .......... .......... .......... 18% 21.9M 6s
    ##  24350K .......... .......... .......... .......... .......... 18% 26.9M 6s
    ##  24400K .......... .......... .......... .......... .......... 18% 16.3M 6s
    ##  24450K .......... .......... .......... .......... .......... 18% 47.9M 6s
    ##  24500K .......... .......... .......... .......... .......... 18% 14.4M 6s
    ##  24550K .......... .......... .......... .......... .......... 18% 47.0M 6s
    ##  24600K .......... .......... .......... .......... .......... 18% 15.5M 6s
    ##  24650K .......... .......... .......... .......... .......... 18% 12.9M 6s
    ##  24700K .......... .......... .......... .......... .......... 18% 40.2M 6s
    ##  24750K .......... .......... .......... .......... .......... 18% 18.5M 6s
    ##  24800K .......... .......... .......... .......... .......... 18% 15.5M 6s
    ##  24850K .......... .......... .......... .......... .......... 18% 46.3M 6s
    ##  24900K .......... .......... .......... .......... .......... 18% 16.7M 6s
    ##  24950K .......... .......... .......... .......... .......... 18% 17.9M 6s
    ##  25000K .......... .......... .......... .......... .......... 18% 43.8M 6s
    ##  25050K .......... .......... .......... .......... .......... 18% 24.0M 6s
    ##  25100K .......... .......... .......... .......... .......... 18% 33.2M 6s
    ##  25150K .......... .......... .......... .......... .......... 18% 12.9M 6s
    ##  25200K .......... .......... .......... .......... .......... 18% 24.2M 6s
    ##  25250K .......... .......... .......... .......... .......... 18% 23.8M 6s
    ##  25300K .......... .......... .......... .......... .......... 18% 13.6M 6s
    ##  25350K .......... .......... .......... .......... .......... 18% 60.7M 6s
    ##  25400K .......... .......... .......... .......... .......... 18% 14.2M 6s
    ##  25450K .......... .......... .......... .......... .......... 19% 21.4M 6s
    ##  25500K .......... .......... .......... .......... .......... 19% 33.2M 6s
    ##  25550K .......... .......... .......... .......... .......... 19% 16.2M 6s
    ##  25600K .......... .......... .......... .......... .......... 19% 44.1M 6s
    ##  25650K .......... .......... .......... .......... .......... 19% 15.2M 6s
    ##  25700K .......... .......... .......... .......... .......... 19% 12.2M 6s
    ##  25750K .......... .......... .......... .......... .......... 19% 78.1M 6s
    ##  25800K .......... .......... .......... .......... .......... 19% 12.5M 6s
    ##  25850K .......... .......... .......... .......... .......... 19% 75.2M 6s
    ##  25900K .......... .......... .......... .......... .......... 19% 14.7M 6s
    ##  25950K .......... .......... .......... .......... .......... 19% 27.8M 6s
    ##  26000K .......... .......... .......... .......... .......... 19% 16.1M 6s
    ##  26050K .......... .......... .......... .......... .......... 19% 43.3M 6s
    ##  26100K .......... .......... .......... .......... .......... 19% 20.6M 6s
    ##  26150K .......... .......... .......... .......... .......... 19% 19.0M 6s
    ##  26200K .......... .......... .......... .......... .......... 19% 25.3M 6s
    ##  26250K .......... .......... .......... .......... .......... 19% 18.1M 6s
    ##  26300K .......... .......... .......... .......... .......... 19% 37.7M 6s
    ##  26350K .......... .......... .......... .......... .......... 19% 17.0M 6s
    ##  26400K .......... .......... .......... .......... .......... 19% 18.0M 6s
    ##  26450K .......... .......... .......... .......... .......... 19% 25.6M 6s
    ##  26500K .......... .......... .......... .......... .......... 19% 15.7M 6s
    ##  26550K .......... .......... .......... .......... .......... 19% 33.2M 6s
    ##  26600K .......... .......... .......... .......... .......... 19% 40.3M 6s
    ##  26650K .......... .......... .......... .......... .......... 19% 18.7M 6s
    ##  26700K .......... .......... .......... .......... .......... 19% 17.4M 6s
    ##  26750K .......... .......... .......... .......... .......... 19% 15.3M 6s
    ##  26800K .......... .......... .......... .......... .......... 20% 43.0M 6s
    ##  26850K .......... .......... .......... .......... .......... 20% 12.6M 6s
    ##  26900K .......... .......... .......... .......... .......... 20% 31.6M 6s
    ##  26950K .......... .......... .......... .......... .......... 20% 18.2M 6s
    ##  27000K .......... .......... .......... .......... .......... 20% 58.6M 6s
    ##  27050K .......... .......... .......... .......... .......... 20% 16.1M 6s
    ##  27100K .......... .......... .......... .......... .......... 20% 37.5M 6s
    ##  27150K .......... .......... .......... .......... .......... 20% 13.5M 6s
    ##  27200K .......... .......... .......... .......... .......... 20% 56.6M 6s
    ##  27250K .......... .......... .......... .......... .......... 20% 15.8M 6s
    ##  27300K .......... .......... .......... .......... .......... 20% 43.8M 6s
    ##  27350K .......... .......... .......... .......... .......... 20% 22.5M 6s
    ##  27400K .......... .......... .......... .......... .......... 20% 17.7M 6s
    ##  27450K .......... .......... .......... .......... .......... 20% 42.5M 6s
    ##  27500K .......... .......... .......... .......... .......... 20% 18.4M 6s
    ##  27550K .......... .......... .......... .......... .......... 20% 19.7M 6s
    ##  27600K .......... .......... .......... .......... .......... 20% 22.5M 6s
    ##  27650K .......... .......... .......... .......... .......... 20% 13.7M 6s
    ##  27700K .......... .......... .......... .......... .......... 20% 41.3M 6s
    ##  27750K .......... .......... .......... .......... .......... 20% 12.9M 6s
    ##  27800K .......... .......... .......... .......... .......... 20% 22.2M 6s
    ##  27850K .......... .......... .......... .......... .......... 20% 53.5M 6s
    ##  27900K .......... .......... .......... .......... .......... 20% 16.8M 6s
    ##  27950K .......... .......... .......... .......... .......... 20% 16.9M 6s
    ##  28000K .......... .......... .......... .......... .......... 20% 34.7M 6s
    ##  28050K .......... .......... .......... .......... .......... 20% 40.1M 6s
    ##  28100K .......... .......... .......... .......... .......... 20% 20.8M 6s
    ##  28150K .......... .......... .......... .......... .......... 21% 15.5M 6s
    ##  28200K .......... .......... .......... .......... .......... 21% 37.2M 6s
    ##  28250K .......... .......... .......... .......... .......... 21% 7.78M 6s
    ##  28300K .......... .......... .......... .......... .......... 21% 43.6M 6s
    ##  28350K .......... .......... .......... .......... .......... 21% 16.4M 6s
    ##  28400K .......... .......... .......... .......... .......... 21% 23.8M 6s
    ##  28450K .......... .......... .......... .......... .......... 21% 42.1M 6s
    ##  28500K .......... .......... .......... .......... .......... 21% 22.9M 6s
    ##  28550K .......... .......... .......... .......... .......... 21% 14.6M 6s
    ##  28600K .......... .......... .......... .......... .......... 21% 46.3M 6s
    ##  28650K .......... .......... .......... .......... .......... 21% 15.7M 6s
    ##  28700K .......... .......... .......... .......... .......... 21% 41.3M 6s
    ##  28750K .......... .......... .......... .......... .......... 21% 15.0M 6s
    ##  28800K .......... .......... .......... .......... .......... 21% 16.2M 6s
    ##  28850K .......... .......... .......... .......... .......... 21% 19.2M 6s
    ##  28900K .......... .......... .......... .......... .......... 21% 29.6M 6s
    ##  28950K .......... .......... .......... .......... .......... 21% 29.4M 6s
    ##  29000K .......... .......... .......... .......... .......... 21% 23.3M 6s
    ##  29050K .......... .......... .......... .......... .......... 21% 34.9M 6s
    ##  29100K .......... .......... .......... .......... .......... 21% 12.5M 6s
    ##  29150K .......... .......... .......... .......... .......... 21% 23.4M 6s
    ##  29200K .......... .......... .......... .......... .......... 21% 57.3M 6s
    ##  29250K .......... .......... .......... .......... .......... 21% 12.2M 6s
    ##  29300K .......... .......... .......... .......... .......... 21% 53.2M 6s
    ##  29350K .......... .......... .......... .......... .......... 21% 9.74M 6s
    ##  29400K .......... .......... .......... .......... .......... 21% 25.3M 6s
    ##  29450K .......... .......... .......... .......... .......... 22% 51.4M 6s
    ##  29500K .......... .......... .......... .......... .......... 22% 13.8M 6s
    ##  29550K .......... .......... .......... .......... .......... 22% 16.1M 6s
    ##  29600K .......... .......... .......... .......... .......... 22% 21.9M 6s
    ##  29650K .......... .......... .......... .......... .......... 22% 15.6M 6s
    ##  29700K .......... .......... .......... .......... .......... 22% 32.7M 6s
    ##  29750K .......... .......... .......... .......... .......... 22% 14.3M 6s
    ##  29800K .......... .......... .......... .......... .......... 22% 24.0M 6s
    ##  29850K .......... .......... .......... .......... .......... 22% 19.5M 6s
    ##  29900K .......... .......... .......... .......... .......... 22% 45.7M 6s
    ##  29950K .......... .......... .......... .......... .......... 22% 11.0M 6s
    ##  30000K .......... .......... .......... .......... .......... 22% 16.4M 6s
    ##  30050K .......... .......... .......... .......... .......... 22% 24.6M 6s
    ##  30100K .......... .......... .......... .......... .......... 22% 49.6M 6s
    ##  30150K .......... .......... .......... .......... .......... 22% 17.7M 6s
    ##  30200K .......... .......... .......... .......... .......... 22% 39.6M 6s
    ##  30250K .......... .......... .......... .......... .......... 22% 11.4M 6s
    ##  30300K .......... .......... .......... .......... .......... 22% 42.6M 6s
    ##  30350K .......... .......... .......... .......... .......... 22% 21.9M 6s
    ##  30400K .......... .......... .......... .......... .......... 22% 15.4M 6s
    ##  30450K .......... .......... .......... .......... .......... 22% 28.6M 6s
    ##  30500K .......... .......... .......... .......... .......... 22% 18.1M 6s
    ##  30550K .......... .......... .......... .......... .......... 22% 19.3M 6s
    ##  30600K .......... .......... .......... .......... .......... 22% 27.8M 6s
    ##  30650K .......... .......... .......... .......... .......... 22% 11.8M 6s
    ##  30700K .......... .......... .......... .......... .......... 22% 77.4M 6s
    ##  30750K .......... .......... .......... .......... .......... 22% 14.6M 6s
    ##  30800K .......... .......... .......... .......... .......... 23% 52.3M 6s
    ##  30850K .......... .......... .......... .......... .......... 23% 14.1M 6s
    ##  30900K .......... .......... .......... .......... .......... 23% 15.1M 6s
    ##  30950K .......... .......... .......... .......... .......... 23% 35.3M 6s
    ##  31000K .......... .......... .......... .......... .......... 23% 42.2M 5s
    ##  31050K .......... .......... .......... .......... .......... 23% 30.0M 5s
    ##  31100K .......... .......... .......... .......... .......... 23% 11.0M 5s
    ##  31150K .......... .......... .......... .......... .......... 23% 17.3M 5s
    ##  31200K .......... .......... .......... .......... .......... 23% 33.1M 5s
    ##  31250K .......... .......... .......... .......... .......... 23% 18.8M 5s
    ##  31300K .......... .......... .......... .......... .......... 23% 19.3M 5s
    ##  31350K .......... .......... .......... .......... .......... 23% 28.2M 5s
    ##  31400K .......... .......... .......... .......... .......... 23% 16.0M 5s
    ##  31450K .......... .......... .......... .......... .......... 23% 31.6M 5s
    ##  31500K .......... .......... .......... .......... .......... 23% 84.1M 5s
    ##  31550K .......... .......... .......... .......... .......... 23% 17.3M 5s
    ##  31600K .......... .......... .......... .......... .......... 23% 14.4M 5s
    ##  31650K .......... .......... .......... .......... .......... 23% 17.6M 5s
    ##  31700K .......... .......... .......... .......... .......... 23% 32.6M 5s
    ##  31750K .......... .......... .......... .......... .......... 23% 17.0M 5s
    ##  31800K .......... .......... .......... .......... .......... 23% 70.0M 5s
    ##  31850K .......... .......... .......... .......... .......... 23% 13.4M 5s
    ##  31900K .......... .......... .......... .......... .......... 23% 13.6M 5s
    ##  31950K .......... .......... .......... .......... .......... 23% 49.2M 5s
    ##  32000K .......... .......... .......... .......... .......... 23% 16.0M 5s
    ##  32050K .......... .......... .......... .......... .......... 23% 36.2M 5s
    ##  32100K .......... .......... .......... .......... .......... 23% 27.2M 5s
    ##  32150K .......... .......... .......... .......... .......... 24% 16.0M 5s
    ##  32200K .......... .......... .......... .......... .......... 24% 39.3M 5s
    ##  32250K .......... .......... .......... .......... .......... 24% 13.6M 5s
    ##  32300K .......... .......... .......... .......... .......... 24% 67.9M 5s
    ##  32350K .......... .......... .......... .......... .......... 24% 16.6M 5s
    ##  32400K .......... .......... .......... .......... .......... 24% 18.0M 5s
    ##  32450K .......... .......... .......... .......... .......... 24% 30.0M 5s
    ##  32500K .......... .......... .......... .......... .......... 24% 13.6M 5s
    ##  32550K .......... .......... .......... .......... .......... 24% 67.9M 5s
    ##  32600K .......... .......... .......... .......... .......... 24% 18.4M 5s
    ##  32650K .......... .......... .......... .......... .......... 24% 30.2M 5s
    ##  32700K .......... .......... .......... .......... .......... 24% 12.7M 5s
    ##  32750K .......... .......... .......... .......... .......... 24% 35.4M 5s
    ##  32800K .......... .......... .......... .......... .......... 24% 16.1M 5s
    ##  32850K .......... .......... .......... .......... .......... 24% 19.1M 5s
    ##  32900K .......... .......... .......... .......... .......... 24% 36.1M 5s
    ##  32950K .......... .......... .......... .......... .......... 24% 12.6M 5s
    ##  33000K .......... .......... .......... .......... .......... 24% 77.1M 5s
    ##  33050K .......... .......... .......... .......... .......... 24% 17.0M 5s
    ##  33100K .......... .......... .......... .......... .......... 24% 30.2M 5s
    ##  33150K .......... .......... .......... .......... .......... 24% 12.7M 5s
    ##  33200K .......... .......... .......... .......... .......... 24% 42.2M 5s
    ##  33250K .......... .......... .......... .......... .......... 24% 20.1M 5s
    ##  33300K .......... .......... .......... .......... .......... 24% 22.5M 5s
    ##  33350K .......... .......... .......... .......... .......... 24% 30.9M 5s
    ##  33400K .......... .......... .......... .......... .......... 24% 15.6M 5s
    ##  33450K .......... .......... .......... .......... .......... 24% 29.8M 5s
    ##  33500K .......... .......... .......... .......... .......... 25% 72.5M 5s
    ##  33550K .......... .......... .......... .......... .......... 25% 18.9M 5s
    ##  33600K .......... .......... .......... .......... .......... 25% 19.7M 5s
    ##  33650K .......... .......... .......... .......... .......... 25% 17.9M 5s
    ##  33700K .......... .......... .......... .......... .......... 25% 22.1M 5s
    ##  33750K .......... .......... .......... .......... .......... 25% 43.8M 5s
    ##  33800K .......... .......... .......... .......... .......... 25% 22.7M 5s
    ##  33850K .......... .......... .......... .......... .......... 25% 19.3M 5s
    ##  33900K .......... .......... .......... .......... .......... 25% 79.9M 5s
    ##  33950K .......... .......... .......... .......... .......... 25% 17.2M 5s
    ##  34000K .......... .......... .......... .......... .......... 25% 18.3M 5s
    ##  34050K .......... .......... .......... .......... .......... 25% 32.5M 5s
    ##  34100K .......... .......... .......... .......... .......... 25% 34.7M 5s
    ##  34150K .......... .......... .......... .......... .......... 25% 21.4M 5s
    ##  34200K .......... .......... .......... .......... .......... 25% 27.2M 5s
    ##  34250K .......... .......... .......... .......... .......... 25% 19.9M 5s
    ##  34300K .......... .......... .......... .......... .......... 25% 59.6M 5s
    ##  34350K .......... .......... .......... .......... .......... 25% 17.8M 5s
    ##  34400K .......... .......... .......... .......... .......... 25% 16.0M 5s
    ##  34450K .......... .......... .......... .......... .......... 25% 39.4M 5s
    ##  34500K .......... .......... .......... .......... .......... 25% 38.7M 5s
    ##  34550K .......... .......... .......... .......... .......... 25% 19.9M 5s
    ##  34600K .......... .......... .......... .......... .......... 25% 39.0M 5s
    ##  34650K .......... .......... .......... .......... .......... 25% 20.1M 5s
    ##  34700K .......... .......... .......... .......... .......... 25% 19.8M 5s
    ##  34750K .......... .......... .......... .......... .......... 25% 19.2M 5s
    ##  34800K .......... .......... .......... .......... .......... 25% 30.6M 5s
    ##  34850K .......... .......... .......... .......... .......... 26% 31.5M 5s
    ##  34900K .......... .......... .......... .......... .......... 26% 25.8M 5s
    ##  34950K .......... .......... .......... .......... .......... 26% 20.4M 5s
    ##  35000K .......... .......... .......... .......... .......... 26% 24.9M 5s
    ##  35050K .......... .......... .......... .......... .......... 26% 17.9M 5s
    ##  35100K .......... .......... .......... .......... .......... 26% 36.1M 5s
    ##  35150K .......... .......... .......... .......... .......... 26% 25.0M 5s
    ##  35200K .......... .......... .......... .......... .......... 26% 19.6M 5s
    ##  35250K .......... .......... .......... .......... .......... 26% 16.0M 5s
    ##  35300K .......... .......... .......... .......... .......... 26% 44.6M 5s
    ##  35350K .......... .......... .......... .......... .......... 26% 10.2M 5s
    ##  35400K .......... .......... .......... .......... .......... 26% 62.9M 5s
    ##  35450K .......... .......... .......... .......... .......... 26% 20.9M 5s
    ##  35500K .......... .......... .......... .......... .......... 26% 46.3M 5s
    ##  35550K .......... .......... .......... .......... .......... 26% 14.4M 5s
    ##  35600K .......... .......... .......... .......... .......... 26% 17.8M 5s
    ##  35650K .......... .......... .......... .......... .......... 26% 40.6M 5s
    ##  35700K .......... .......... .......... .......... .......... 26% 29.4M 5s
    ##  35750K .......... .......... .......... .......... .......... 26% 30.0M 5s
    ##  35800K .......... .......... .......... .......... .......... 26% 20.0M 5s
    ##  35850K .......... .......... .......... .......... .......... 26% 14.7M 5s
    ##  35900K .......... .......... .......... .......... .......... 26% 50.0M 5s
    ##  35950K .......... .......... .......... .......... .......... 26% 16.1M 5s
    ##  36000K .......... .......... .......... .......... .......... 26% 59.8M 5s
    ##  36050K .......... .......... .......... .......... .......... 26% 14.7M 5s
    ##  36100K .......... .......... .......... .......... .......... 26% 50.1M 5s
    ##  36150K .......... .......... .......... .......... .......... 27% 14.4M 5s
    ##  36200K .......... .......... .......... .......... .......... 27% 82.1M 5s
    ##  36250K .......... .......... .......... .......... .......... 27% 20.1M 5s
    ##  36300K .......... .......... .......... .......... .......... 27% 17.2M 5s
    ##  36350K .......... .......... .......... .......... .......... 27% 13.5M 5s
    ##  36400K .......... .......... .......... .......... .......... 27% 32.2M 5s
    ##  36450K .......... .......... .......... .......... .......... 27% 89.4M 5s
    ##  36500K .......... .......... .......... .......... .......... 27% 15.8M 5s
    ##  36550K .......... .......... .......... .......... .......... 27% 29.6M 5s
    ##  36600K .......... .......... .......... .......... .......... 27% 21.3M 5s
    ##  36650K .......... .......... .......... .......... .......... 27% 23.8M 5s
    ##  36700K .......... .......... .......... .......... .......... 27% 56.1M 5s
    ##  36750K .......... .......... .......... .......... .......... 27% 14.6M 5s
    ##  36800K .......... .......... .......... .......... .......... 27% 26.2M 5s
    ##  36850K .......... .......... .......... .......... .......... 27% 23.3M 5s
    ##  36900K .......... .......... .......... .......... .......... 27% 35.9M 5s
    ##  36950K .......... .......... .......... .......... .......... 27% 19.8M 5s
    ##  37000K .......... .......... .......... .......... .......... 27% 34.3M 5s
    ##  37050K .......... .......... .......... .......... .......... 27% 19.8M 5s
    ##  37100K .......... .......... .......... .......... .......... 27% 21.1M 5s
    ##  37150K .......... .......... .......... .......... .......... 27% 23.6M 5s
    ##  37200K .......... .......... .......... .......... .......... 27% 25.6M 5s
    ##  37250K .......... .......... .......... .......... .......... 27% 18.5M 5s
    ##  37300K .......... .......... .......... .......... .......... 27% 20.8M 5s
    ##  37350K .......... .......... .......... .......... .......... 27% 48.6M 5s
    ##  37400K .......... .......... .......... .......... .......... 27% 20.8M 5s
    ##  37450K .......... .......... .......... .......... .......... 27% 38.8M 5s
    ##  37500K .......... .......... .......... .......... .......... 28% 19.2M 5s
    ##  37550K .......... .......... .......... .......... .......... 28% 22.9M 5s
    ##  37600K .......... .......... .......... .......... .......... 28% 24.9M 5s
    ##  37650K .......... .......... .......... .......... .......... 28% 18.6M 5s
    ##  37700K .......... .......... .......... .......... .......... 28% 34.9M 5s
    ##  37750K .......... .......... .......... .......... .......... 28% 14.3M 5s
    ##  37800K .......... .......... .......... .......... .......... 28% 43.2M 5s
    ##  37850K .......... .......... .......... .......... .......... 28% 19.4M 5s
    ##  37900K .......... .......... .......... .......... .......... 28% 36.0M 5s
    ##  37950K .......... .......... .......... .......... .......... 28% 18.5M 5s
    ##  38000K .......... .......... .......... .......... .......... 28% 35.9M 5s
    ##  38050K .......... .......... .......... .......... .......... 28% 38.7M 5s
    ##  38100K .......... .......... .......... .......... .......... 28% 18.3M 5s
    ##  38150K .......... .......... .......... .......... .......... 28% 25.2M 5s
    ##  38200K .......... .......... .......... .......... .......... 28% 23.2M 5s
    ##  38250K .......... .......... .......... .......... .......... 28% 18.3M 5s
    ##  38300K .......... .......... .......... .......... .......... 28% 18.8M 5s
    ##  38350K .......... .......... .......... .......... .......... 28% 36.0M 5s
    ##  38400K .......... .......... .......... .......... .......... 28% 16.6M 5s
    ##  38450K .......... .......... .......... .......... .......... 28% 33.8M 5s
    ##  38500K .......... .......... .......... .......... .......... 28% 61.9M 5s
    ##  38550K .......... .......... .......... .......... .......... 28% 23.5M 5s
    ##  38600K .......... .......... .......... .......... .......... 28% 23.7M 5s
    ##  38650K .......... .......... .......... .......... .......... 28% 25.4M 5s
    ##  38700K .......... .......... .......... .......... .......... 28% 21.9M 5s
    ##  38750K .......... .......... .......... .......... .......... 28% 29.8M 5s
    ##  38800K .......... .......... .......... .......... .......... 28% 13.8M 5s
    ##  38850K .......... .......... .......... .......... .......... 29% 17.8M 5s
    ##  38900K .......... .......... .......... .......... .......... 29% 33.8M 5s
    ##  38950K .......... .......... .......... .......... .......... 29% 48.8M 5s
    ##  39000K .......... .......... .......... .......... .......... 29% 18.0M 5s
    ##  39050K .......... .......... .......... .......... .......... 29% 14.9M 5s
    ##  39100K .......... .......... .......... .......... .......... 29% 41.2M 5s
    ##  39150K .......... .......... .......... .......... .......... 29% 17.4M 5s
    ##  39200K .......... .......... .......... .......... .......... 29% 41.1M 5s
    ##  39250K .......... .......... .......... .......... .......... 29% 19.9M 5s
    ##  39300K .......... .......... .......... .......... .......... 29% 26.1M 5s
    ##  39350K .......... .......... .......... .......... .......... 29% 37.1M 5s
    ##  39400K .......... .......... .......... .......... .......... 29% 29.7M 5s
    ##  39450K .......... .......... .......... .......... .......... 29% 14.2M 5s
    ##  39500K .......... .......... .......... .......... .......... 29% 24.3M 5s
    ##  39550K .......... .......... .......... .......... .......... 29% 21.5M 5s
    ##  39600K .......... .......... .......... .......... .......... 29% 32.2M 5s
    ##  39650K .......... .......... .......... .......... .......... 29% 16.2M 5s
    ##  39700K .......... .......... .......... .......... .......... 29% 24.3M 5s
    ##  39750K .......... .......... .......... .......... .......... 29% 24.9M 5s
    ##  39800K .......... .......... .......... .......... .......... 29% 26.8M 5s
    ##  39850K .......... .......... .......... .......... .......... 29% 23.9M 5s
    ##  39900K .......... .......... .......... .......... .......... 29% 12.3M 5s
    ##  39950K .......... .......... .......... .......... .......... 29% 38.5M 5s
    ##  40000K .......... .......... .......... .......... .......... 29% 15.1M 5s
    ##  40050K .......... .......... .......... .......... .......... 29% 43.2M 5s
    ##  40100K .......... .......... .......... .......... .......... 29% 19.3M 5s
    ##  40150K .......... .......... .......... .......... .......... 29% 30.6M 5s
    ##  40200K .......... .......... .......... .......... .......... 30% 14.5M 5s
    ##  40250K .......... .......... .......... .......... .......... 30% 37.4M 5s
    ##  40300K .......... .......... .......... .......... .......... 30% 20.8M 5s
    ##  40350K .......... .......... .......... .......... .......... 30% 30.7M 5s
    ##  40400K .......... .......... .......... .......... .......... 30% 12.9M 5s
    ##  40450K .......... .......... .......... .......... .......... 30% 49.3M 5s
    ##  40500K .......... .......... .......... .......... .......... 30% 14.5M 5s
    ##  40550K .......... .......... .......... .......... .......... 30% 38.7M 5s
    ##  40600K .......... .......... .......... .......... .......... 30% 18.0M 5s
    ##  40650K .......... .......... .......... .......... .......... 30% 30.3M 5s
    ##  40700K .......... .......... .......... .......... .......... 30% 19.4M 5s
    ##  40750K .......... .......... .......... .......... .......... 30% 25.6M 5s
    ##  40800K .......... .......... .......... .......... .......... 30% 26.7M 5s
    ##  40850K .......... .......... .......... .......... .......... 30% 15.2M 5s
    ##  40900K .......... .......... .......... .......... .......... 30% 73.2M 5s
    ##  40950K .......... .......... .......... .......... .......... 30% 40.2M 5s
    ##  41000K .......... .......... .......... .......... .......... 30% 15.3M 5s
    ##  41050K .......... .......... .......... .......... .......... 30% 16.1M 5s
    ##  41100K .......... .......... .......... .......... .......... 30% 38.9M 5s
    ##  41150K .......... .......... .......... .......... .......... 30% 12.8M 5s
    ##  41200K .......... .......... .......... .......... .......... 30% 38.1M 5s
    ##  41250K .......... .......... .......... .......... .......... 30% 14.9M 5s
    ##  41300K .......... .......... .......... .......... .......... 30% 37.8M 5s
    ##  41350K .......... .......... .......... .......... .......... 30% 72.0M 5s
    ##  41400K .......... .......... .......... .......... .......... 30% 16.0M 5s
    ##  41450K .......... .......... .......... .......... .......... 30% 19.8M 5s
    ##  41500K .......... .......... .......... .......... .......... 30% 38.2M 5s
    ##  41550K .......... .......... .......... .......... .......... 31% 15.1M 5s
    ##  41600K .......... .......... .......... .......... .......... 31% 53.9M 5s
    ##  41650K .......... .......... .......... .......... .......... 31% 15.6M 5s
    ##  41700K .......... .......... .......... .......... .......... 31% 14.4M 5s
    ##  41750K .......... .......... .......... .......... .......... 31% 38.2M 5s
    ##  41800K .......... .......... .......... .......... .......... 31% 21.3M 5s
    ##  41850K .......... .......... .......... .......... .......... 31% 55.3M 5s
    ##  41900K .......... .......... .......... .......... .......... 31% 14.5M 5s
    ##  41950K .......... .......... .......... .......... .......... 31% 12.4M 5s
    ##  42000K .......... .......... .......... .......... .......... 31% 39.8M 5s
    ##  42050K .......... .......... .......... .......... .......... 31% 25.3M 5s
    ##  42100K .......... .......... .......... .......... .......... 31% 30.1M 5s
    ##  42150K .......... .......... .......... .......... .......... 31% 17.6M 5s
    ##  42200K .......... .......... .......... .......... .......... 31% 26.0M 5s
    ##  42250K .......... .......... .......... .......... .......... 31% 62.7M 5s
    ##  42300K .......... .......... .......... .......... .......... 31% 27.6M 5s
    ##  42350K .......... .......... .......... .......... .......... 31% 23.1M 5s
    ##  42400K .......... .......... .......... .......... .......... 31% 13.4M 5s
    ##  42450K .......... .......... .......... .......... .......... 31% 28.6M 5s
    ##  42500K .......... .......... .......... .......... .......... 31% 92.5M 5s
    ##  42550K .......... .......... .......... .......... .......... 31% 17.7M 5s
    ##  42600K .......... .......... .......... .......... .......... 31% 45.9M 5s
    ##  42650K .......... .......... .......... .......... .......... 31% 14.5M 5s
    ##  42700K .......... .......... .......... .......... .......... 31% 41.1M 5s
    ##  42750K .......... .......... .......... .......... .......... 31% 12.7M 5s
    ##  42800K .......... .......... .......... .......... .......... 31% 60.0M 5s
    ##  42850K .......... .......... .......... .......... .......... 31% 14.2M 5s
    ##  42900K .......... .......... .......... .......... .......... 32% 40.8M 5s
    ##  42950K .......... .......... .......... .......... .......... 32% 16.0M 5s
    ##  43000K .......... .......... .......... .......... .......... 32% 42.1M 5s
    ##  43050K .......... .......... .......... .......... .......... 32% 15.2M 5s
    ##  43100K .......... .......... .......... .......... .......... 32% 31.7M 5s
    ##  43150K .......... .......... .......... .......... .......... 32% 16.6M 5s
    ##  43200K .......... .......... .......... .......... .......... 32% 38.8M 5s
    ##  43250K .......... .......... .......... .......... .......... 32% 16.8M 5s
    ##  43300K .......... .......... .......... .......... .......... 32% 33.6M 5s
    ##  43350K .......... .......... .......... .......... .......... 32% 65.2M 5s
    ##  43400K .......... .......... .......... .......... .......... 32% 19.1M 5s
    ##  43450K .......... .......... .......... .......... .......... 32% 22.6M 5s
    ##  43500K .......... .......... .......... .......... .......... 32% 25.2M 5s
    ##  43550K .......... .......... .......... .......... .......... 32% 20.3M 5s
    ##  43600K .......... .......... .......... .......... .......... 32% 50.6M 5s
    ##  43650K .......... .......... .......... .......... .......... 32% 15.1M 5s
    ##  43700K .......... .......... .......... .......... .......... 32% 22.1M 5s
    ##  43750K .......... .......... .......... .......... .......... 32% 29.0M 5s
    ##  43800K .......... .......... .......... .......... .......... 32% 23.6M 5s
    ##  43850K .......... .......... .......... .......... .......... 32% 18.1M 5s
    ##  43900K .......... .......... .......... .......... .......... 32% 33.0M 5s
    ##  43950K .......... .......... .......... .......... .......... 32% 18.3M 5s
    ##  44000K .......... .......... .......... .......... .......... 32% 37.0M 5s
    ##  44050K .......... .......... .......... .......... .......... 32% 17.3M 5s
    ##  44100K .......... .......... .......... .......... .......... 32% 14.7M 5s
    ##  44150K .......... .......... .......... .......... .......... 32% 27.3M 5s
    ##  44200K .......... .......... .......... .......... .......... 33% 14.4M 5s
    ##  44250K .......... .......... .......... .......... .......... 33% 26.2M 5s
    ##  44300K .......... .......... .......... .......... .......... 33% 75.2M 4s
    ##  44350K .......... .......... .......... .......... .......... 33% 16.2M 4s
    ##  44400K .......... .......... .......... .......... .......... 33% 13.1M 4s
    ##  44450K .......... .......... .......... .......... .......... 33% 80.5M 4s
    ##  44500K .......... .......... .......... .......... .......... 33% 16.3M 4s
    ##  44550K .......... .......... .......... .......... .......... 33% 37.7M 4s
    ##  44600K .......... .......... .......... .......... .......... 33% 15.8M 4s
    ##  44650K .......... .......... .......... .......... .......... 33% 8.08M 4s
    ##  44700K .......... .......... .......... .......... .......... 33% 93.6M 4s
    ##  44750K .......... .......... .......... .......... .......... 33% 14.7M 4s
    ##  44800K .......... .......... .......... .......... .......... 33% 19.9M 4s
    ##  44850K .......... .......... .......... .......... .......... 33% 29.2M 4s
    ##  44900K .......... .......... .......... .......... .......... 33% 20.5M 4s
    ##  44950K .......... .......... .......... .......... .......... 33% 26.0M 4s
    ##  45000K .......... .......... .......... .......... .......... 33% 59.7M 4s
    ##  45050K .......... .......... .......... .......... .......... 33% 17.5M 4s
    ##  45100K .......... .......... .......... .......... .......... 33% 34.9M 4s
    ##  45150K .......... .......... .......... .......... .......... 33% 17.2M 4s
    ##  45200K .......... .......... .......... .......... .......... 33% 26.2M 4s
    ##  45250K .......... .......... .......... .......... .......... 33% 22.6M 4s
    ##  45300K .......... .......... .......... .......... .......... 33% 15.5M 4s
    ##  45350K .......... .......... .......... .......... .......... 33% 35.9M 4s
    ##  45400K .......... .......... .......... .......... .......... 33% 20.8M 4s
    ##  45450K .......... .......... .......... .......... .......... 33% 20.6M 4s
    ##  45500K .......... .......... .......... .......... .......... 33% 75.2M 4s
    ##  45550K .......... .......... .......... .......... .......... 34% 14.5M 4s
    ##  45600K .......... .......... .......... .......... .......... 34% 23.1M 4s
    ##  45650K .......... .......... .......... .......... .......... 34% 26.4M 4s
    ##  45700K .......... .......... .......... .......... .......... 34% 35.7M 4s
    ##  45750K .......... .......... .......... .......... .......... 34% 18.6M 4s
    ##  45800K .......... .......... .......... .......... .......... 34% 22.2M 4s
    ##  45850K .......... .......... .......... .......... .......... 34% 38.3M 4s
    ##  45900K .......... .......... .......... .......... .......... 34% 16.1M 4s
    ##  45950K .......... .......... .......... .......... .......... 34% 22.4M 4s
    ##  46000K .......... .......... .......... .......... .......... 34% 15.0M 4s
    ##  46050K .......... .......... .......... .......... .......... 34% 51.9M 4s
    ##  46100K .......... .......... .......... .......... .......... 34% 25.6M 4s
    ##  46150K .......... .......... .......... .......... .......... 34% 20.7M 4s
    ##  46200K .......... .......... .......... .......... .......... 34% 13.3M 4s
    ##  46250K .......... .......... .......... .......... .......... 34% 16.0M 4s
    ##  46300K .......... .......... .......... .......... .......... 34% 32.0M 4s
    ##  46350K .......... .......... .......... .......... .......... 34% 27.5M 4s
    ##  46400K .......... .......... .......... .......... .......... 34% 30.5M 4s
    ##  46450K .......... .......... .......... .......... .......... 34% 16.1M 4s
    ##  46500K .......... .......... .......... .......... .......... 34% 17.2M 4s
    ##  46550K .......... .......... .......... .......... .......... 34% 31.2M 4s
    ##  46600K .......... .......... .......... .......... .......... 34% 16.0M 4s
    ##  46650K .......... .......... .......... .......... .......... 34% 40.3M 4s
    ##  46700K .......... .......... .......... .......... .......... 34% 13.1M 4s
    ##  46750K .......... .......... .......... .......... .......... 34% 46.3M 4s
    ##  46800K .......... .......... .......... .......... .......... 34% 16.6M 4s
    ##  46850K .......... .......... .......... .......... .......... 34% 31.7M 4s
    ##  46900K .......... .......... .......... .......... .......... 35% 20.4M 4s
    ##  46950K .......... .......... .......... .......... .......... 35% 18.2M 4s
    ##  47000K .......... .......... .......... .......... .......... 35% 24.0M 4s
    ##  47050K .......... .......... .......... .......... .......... 35% 22.8M 4s
    ##  47100K .......... .......... .......... .......... .......... 35% 27.2M 4s
    ##  47150K .......... .......... .......... .......... .......... 35% 17.5M 4s
    ##  47200K .......... .......... .......... .......... .......... 35% 24.8M 4s
    ##  47250K .......... .......... .......... .......... .......... 35% 18.6M 4s
    ##  47300K .......... .......... .......... .......... .......... 35% 16.5M 4s
    ##  47350K .......... .......... .......... .......... .......... 35% 40.9M 4s
    ##  47400K .......... .......... .......... .......... .......... 35% 12.6M 4s
    ##  47450K .......... .......... .......... .......... .......... 35% 29.8M 4s
    ##  47500K .......... .......... .......... .......... .......... 35% 28.9M 4s
    ##  47550K .......... .......... .......... .......... .......... 35% 17.1M 4s
    ##  47600K .......... .......... .......... .......... .......... 35% 17.1M 4s
    ##  47650K .......... .......... .......... .......... .......... 35% 31.2M 4s
    ##  47700K .......... .......... .......... .......... .......... 35% 34.6M 4s
    ##  47750K .......... .......... .......... .......... .......... 35% 14.1M 4s
    ##  47800K .......... .......... .......... .......... .......... 35% 43.2M 4s
    ##  47850K .......... .......... .......... .......... .......... 35% 18.1M 4s
    ##  47900K .......... .......... .......... .......... .......... 35% 38.7M 4s
    ##  47950K .......... .......... .......... .......... .......... 35% 23.2M 4s
    ##  48000K .......... .......... .......... .......... .......... 35% 22.2M 4s
    ##  48050K .......... .......... .......... .......... .......... 35% 19.5M 4s
    ##  48100K .......... .......... .......... .......... .......... 35% 35.3M 4s
    ##  48150K .......... .......... .......... .......... .......... 35% 15.4M 4s
    ##  48200K .......... .......... .......... .......... .......... 35% 20.7M 4s
    ##  48250K .......... .......... .......... .......... .......... 36% 23.7M 4s
    ##  48300K .......... .......... .......... .......... .......... 36% 61.0M 4s
    ##  48350K .......... .......... .......... .......... .......... 36% 22.6M 4s
    ##  48400K .......... .......... .......... .......... .......... 36% 29.0M 4s
    ##  48450K .......... .......... .......... .......... .......... 36% 19.0M 4s
    ##  48500K .......... .......... .......... .......... .......... 36% 31.9M 4s
    ##  48550K .......... .......... .......... .......... .......... 36% 23.3M 4s
    ##  48600K .......... .......... .......... .......... .......... 36% 15.6M 4s
    ##  48650K .......... .......... .......... .......... .......... 36% 35.9M 4s
    ##  48700K .......... .......... .......... .......... .......... 36% 13.8M 4s
    ##  48750K .......... .......... .......... .......... .......... 36% 31.5M 4s
    ##  48800K .......... .......... .......... .......... .......... 36% 17.7M 4s
    ##  48850K .......... .......... .......... .......... .......... 36% 32.4M 4s
    ##  48900K .......... .......... .......... .......... .......... 36% 22.4M 4s
    ##  48950K .......... .......... .......... .......... .......... 36% 19.9M 4s
    ##  49000K .......... .......... .......... .......... .......... 36% 42.2M 4s
    ##  49050K .......... .......... .......... .......... .......... 36% 21.7M 4s
    ##  49100K .......... .......... .......... .......... .......... 36% 13.3M 4s
    ##  49150K .......... .......... .......... .......... .......... 36% 20.4M 4s
    ##  49200K .......... .......... .......... .......... .......... 36% 21.7M 4s
    ##  49250K .......... .......... .......... .......... .......... 36% 20.1M 4s
    ##  49300K .......... .......... .......... .......... .......... 36% 41.3M 4s
    ##  49350K .......... .......... .......... .......... .......... 36% 13.6M 4s
    ##  49400K .......... .......... .......... .......... .......... 36% 22.9M 4s
    ##  49450K .......... .......... .......... .......... .......... 36% 29.9M 4s
    ##  49500K .......... .......... .......... .......... .......... 36% 19.4M 4s
    ##  49550K .......... .......... .......... .......... .......... 36% 14.9M 4s
    ##  49600K .......... .......... .......... .......... .......... 37% 41.6M 4s
    ##  49650K .......... .......... .......... .......... .......... 37% 14.9M 4s
    ##  49700K .......... .......... .......... .......... .......... 37% 39.3M 4s
    ##  49750K .......... .......... .......... .......... .......... 37% 24.0M 4s
    ##  49800K .......... .......... .......... .......... .......... 37% 19.3M 4s
    ##  49850K .......... .......... .......... .......... .......... 37% 16.8M 4s
    ##  49900K .......... .......... .......... .......... .......... 37% 39.5M 4s
    ##  49950K .......... .......... .......... .......... .......... 37% 12.1M 4s
    ##  50000K .......... .......... .......... .......... .......... 37% 86.6M 4s
    ##  50050K .......... .......... .......... .......... .......... 37% 13.0M 4s
    ##  50100K .......... .......... .......... .......... .......... 37% 61.2M 4s
    ##  50150K .......... .......... .......... .......... .......... 37% 13.4M 4s
    ##  50200K .......... .......... .......... .......... .......... 37% 44.8M 4s
    ##  50250K .......... .......... .......... .......... .......... 37% 16.8M 4s
    ##  50300K .......... .......... .......... .......... .......... 37% 39.3M 4s
    ##  50350K .......... .......... .......... .......... .......... 37% 13.3M 4s
    ##  50400K .......... .......... .......... .......... .......... 37% 42.8M 4s
    ##  50450K .......... .......... .......... .......... .......... 37% 14.6M 4s
    ##  50500K .......... .......... .......... .......... .......... 37% 38.7M 4s
    ##  50550K .......... .......... .......... .......... .......... 37% 15.6M 4s
    ##  50600K .......... .......... .......... .......... .......... 37% 31.8M 4s
    ##  50650K .......... .......... .......... .......... .......... 37% 58.4M 4s
    ##  50700K .......... .......... .......... .......... .......... 37% 20.4M 4s
    ##  50750K .......... .......... .......... .......... .......... 37% 12.9M 4s
    ##  50800K .......... .......... .......... .......... .......... 37% 29.3M 4s
    ##  50850K .......... .......... .......... .......... .......... 37% 57.9M 4s
    ##  50900K .......... .......... .......... .......... .......... 38% 30.6M 4s
    ##  50950K .......... .......... .......... .......... .......... 38% 13.8M 4s
    ##  51000K .......... .......... .......... .......... .......... 38% 49.7M 4s
    ##  51050K .......... .......... .......... .......... .......... 38% 14.0M 4s
    ##  51100K .......... .......... .......... .......... .......... 38% 49.5M 4s
    ##  51150K .......... .......... .......... .......... .......... 38% 18.9M 4s
    ##  51200K .......... .......... .......... .......... .......... 38% 14.8M 4s
    ##  51250K .......... .......... .......... .......... .......... 38% 16.9M 4s
    ##  51300K .......... .......... .......... .......... .......... 38% 21.6M 4s
    ##  51350K .......... .......... .......... .......... .......... 38% 52.5M 4s
    ##  51400K .......... .......... .......... .......... .......... 38% 15.4M 4s
    ##  51450K .......... .......... .......... .......... .......... 38% 58.0M 4s
    ##  51500K .......... .......... .......... .......... .......... 38% 14.0M 4s
    ##  51550K .......... .......... .......... .......... .......... 38% 14.1M 4s
    ##  51600K .......... .......... .......... .......... .......... 38% 19.1M 4s
    ##  51650K .......... .......... .......... .......... .......... 38% 21.6M 4s
    ##  51700K .......... .......... .......... .......... .......... 38% 56.3M 4s
    ##  51750K .......... .......... .......... .......... .......... 38% 19.5M 4s
    ##  51800K .......... .......... .......... .......... .......... 38% 40.1M 4s
    ##  51850K .......... .......... .......... .......... .......... 38% 16.1M 4s
    ##  51900K .......... .......... .......... .......... .......... 38% 36.9M 4s
    ##  51950K .......... .......... .......... .......... .......... 38% 21.3M 4s
    ##  52000K .......... .......... .......... .......... .......... 38% 15.0M 4s
    ##  52050K .......... .......... .......... .......... .......... 38% 36.6M 4s
    ##  52100K .......... .......... .......... .......... .......... 38% 17.7M 4s
    ##  52150K .......... .......... .......... .......... .......... 38% 27.2M 4s
    ##  52200K .......... .......... .......... .......... .......... 38% 80.6M 4s
    ##  52250K .......... .......... .......... .......... .......... 39% 18.1M 4s
    ##  52300K .......... .......... .......... .......... .......... 39% 15.3M 4s
    ##  52350K .......... .......... .......... .......... .......... 39% 15.2M 4s
    ##  52400K .......... .......... .......... .......... .......... 39% 32.3M 4s
    ##  52450K .......... .......... .......... .......... .......... 39% 59.6M 4s
    ##  52500K .......... .......... .......... .......... .......... 39% 21.4M 4s
    ##  52550K .......... .......... .......... .......... .......... 39% 15.5M 4s
    ##  52600K .......... .......... .......... .......... .......... 39% 11.7M 4s
    ##  52650K .......... .......... .......... .......... .......... 39% 35.8M 4s
    ##  52700K .......... .......... .......... .......... .......... 39% 84.0M 4s
    ##  52750K .......... .......... .......... .......... .......... 39% 17.2M 4s
    ##  52800K .......... .......... .......... .......... .......... 39% 31.1M 4s
    ##  52850K .......... .......... .......... .......... .......... 39% 21.3M 4s
    ##  52900K .......... .......... .......... .......... .......... 39% 24.2M 4s
    ##  52950K .......... .......... .......... .......... .......... 39% 50.1M 4s
    ##  53000K .......... .......... .......... .......... .......... 39% 15.0M 4s
    ##  53050K .......... .......... .......... .......... .......... 39% 52.9M 4s
    ##  53100K .......... .......... .......... .......... .......... 39% 14.1M 4s
    ##  53150K .......... .......... .......... .......... .......... 39% 14.2M 4s
    ##  53200K .......... .......... .......... .......... .......... 39% 34.8M 4s
    ##  53250K .......... .......... .......... .......... .......... 39% 51.8M 4s
    ##  53300K .......... .......... .......... .......... .......... 39% 21.6M 4s
    ##  53350K .......... .......... .......... .......... .......... 39% 25.1M 4s
    ##  53400K .......... .......... .......... .......... .......... 39% 25.7M 4s
    ##  53450K .......... .......... .......... .......... .......... 39% 19.3M 4s
    ##  53500K .......... .......... .......... .......... .......... 39% 50.1M 4s
    ##  53550K .......... .......... .......... .......... .......... 39% 16.2M 4s
    ##  53600K .......... .......... .......... .......... .......... 40% 16.9M 4s
    ##  53650K .......... .......... .......... .......... .......... 40% 27.3M 4s
    ##  53700K .......... .......... .......... .......... .......... 40% 64.1M 4s
    ##  53750K .......... .......... .......... .......... .......... 40% 13.0M 4s
    ##  53800K .......... .......... .......... .......... .......... 40% 97.7M 4s
    ##  53850K .......... .......... .......... .......... .......... 40% 15.4M 4s
    ##  53900K .......... .......... .......... .......... .......... 40% 35.1M 4s
    ##  53950K .......... .......... .......... .......... .......... 40% 24.1M 4s
    ##  54000K .......... .......... .......... .......... .......... 40% 18.6M 4s
    ##  54050K .......... .......... .......... .......... .......... 40% 47.4M 4s
    ##  54100K .......... .......... .......... .......... .......... 40% 12.5M 4s
    ##  54150K .......... .......... .......... .......... .......... 40% 28.4M 4s
    ##  54200K .......... .......... .......... .......... .......... 40% 28.7M 4s
    ##  54250K .......... .......... .......... .......... .......... 40% 25.4M 4s
    ##  54300K .......... .......... .......... .......... .......... 40% 20.7M 4s
    ##  54350K .......... .......... .......... .......... .......... 40% 26.4M 4s
    ##  54400K .......... .......... .......... .......... .......... 40% 28.4M 4s
    ##  54450K .......... .......... .......... .......... .......... 40% 15.3M 4s
    ##  54500K .......... .......... .......... .......... .......... 40% 24.1M 4s
    ##  54550K .......... .......... .......... .......... .......... 40% 39.2M 4s
    ##  54600K .......... .......... .......... .......... .......... 40% 26.0M 4s
    ##  54650K .......... .......... .......... .......... .......... 40% 25.2M 4s
    ##  54700K .......... .......... .......... .......... .......... 40% 24.9M 4s
    ##  54750K .......... .......... .......... .......... .......... 40% 19.9M 4s
    ##  54800K .......... .......... .......... .......... .......... 40% 37.8M 4s
    ##  54850K .......... .......... .......... .......... .......... 40% 16.3M 4s
    ##  54900K .......... .......... .......... .......... .......... 40% 22.5M 4s
    ##  54950K .......... .......... .......... .......... .......... 41% 36.9M 4s
    ##  55000K .......... .......... .......... .......... .......... 41% 28.9M 4s
    ##  55050K .......... .......... .......... .......... .......... 41% 45.6M 4s
    ##  55100K .......... .......... .......... .......... .......... 41% 18.6M 4s
    ##  55150K .......... .......... .......... .......... .......... 41% 17.2M 4s
    ##  55200K .......... .......... .......... .......... .......... 41% 38.8M 4s
    ##  55250K .......... .......... .......... .......... .......... 41% 92.2M 4s
    ##  55300K .......... .......... .......... .......... .......... 41% 14.1M 4s
    ##  55350K .......... .......... .......... .......... .......... 41% 78.8M 4s
    ##  55400K .......... .......... .......... .......... .......... 41% 13.9M 4s
    ##  55450K .......... .......... .......... .......... .......... 41% 90.6M 4s
    ##  55500K .......... .......... .......... .......... .......... 41% 23.3M 4s
    ##  55550K .......... .......... .......... .......... .......... 41% 16.5M 4s
    ##  55600K .......... .......... .......... .......... .......... 41% 6.20M 4s
    ##  55650K .......... .......... .......... .......... .......... 41% 91.4M 4s
    ##  55700K .......... .......... .......... .......... .......... 41% 16.3M 4s
    ##  55750K .......... .......... .......... .......... .......... 41% 81.2M 4s
    ##  55800K .......... .......... .......... .......... .......... 41% 14.3M 4s
    ##  55850K .......... .......... .......... .......... .......... 41% 95.9M 4s
    ##  55900K .......... .......... .......... .......... .......... 41% 8.72M 4s
    ##  55950K .......... .......... .......... .......... .......... 41% 15.5M 4s
    ##  56000K .......... .......... .......... .......... .......... 41% 47.1M 4s
    ##  56050K .......... .......... .......... .......... .......... 41% 17.7M 4s
    ##  56100K .......... .......... .......... .......... .......... 41% 33.6M 4s
    ##  56150K .......... .......... .......... .......... .......... 41% 39.8M 4s
    ##  56200K .......... .......... .......... .......... .......... 41% 23.8M 4s
    ##  56250K .......... .......... .......... .......... .......... 41% 23.9M 4s
    ##  56300K .......... .......... .......... .......... .......... 42% 15.6M 4s
    ##  56350K .......... .......... .......... .......... .......... 42% 32.0M 4s
    ##  56400K .......... .......... .......... .......... .......... 42% 31.4M 4s
    ##  56450K .......... .......... .......... .......... .......... 42% 17.1M 4s
    ##  56500K .......... .......... .......... .......... .......... 42% 28.7M 4s
    ##  56550K .......... .......... .......... .......... .......... 42% 16.6M 4s
    ##  56600K .......... .......... .......... .......... .......... 42% 58.1M 4s
    ##  56650K .......... .......... .......... .......... .......... 42% 19.9M 4s
    ##  56700K .......... .......... .......... .......... .......... 42% 38.5M 4s
    ##  56750K .......... .......... .......... .......... .......... 42% 17.5M 4s
    ##  56800K .......... .......... .......... .......... .......... 42% 42.4M 4s
    ##  56850K .......... .......... .......... .......... .......... 42% 30.9M 4s
    ##  56900K .......... .......... .......... .......... .......... 42% 14.7M 4s
    ##  56950K .......... .......... .......... .......... .......... 42% 60.4M 4s
    ##  57000K .......... .......... .......... .......... .......... 42% 15.4M 4s
    ##  57050K .......... .......... .......... .......... .......... 42% 14.8M 4s
    ##  57100K .......... .......... .......... .......... .......... 42% 17.4M 4s
    ##  57150K .......... .......... .......... .......... .......... 42% 30.4M 4s
    ##  57200K .......... .......... .......... .......... .......... 42% 24.1M 4s
    ##  57250K .......... .......... .......... .......... .......... 42% 45.3M 4s
    ##  57300K .......... .......... .......... .......... .......... 42% 16.8M 4s
    ##  57350K .......... .......... .......... .......... .......... 42% 45.2M 4s
    ##  57400K .......... .......... .......... .......... .......... 42% 10.7M 4s
    ##  57450K .......... .......... .......... .......... .......... 42% 92.0M 4s
    ##  57500K .......... .......... .......... .......... .......... 42% 19.4M 4s
    ##  57550K .......... .......... .......... .......... .......... 42% 15.5M 4s
    ##  57600K .......... .......... .......... .......... .......... 43% 48.9M 4s
    ##  57650K .......... .......... .......... .......... .......... 43% 14.7M 4s
    ##  57700K .......... .......... .......... .......... .......... 43% 48.0M 4s
    ##  57750K .......... .......... .......... .......... .......... 43% 21.4M 4s
    ##  57800K .......... .......... .......... .......... .......... 43% 35.3M 4s
    ##  57850K .......... .......... .......... .......... .......... 43% 12.8M 4s
    ##  57900K .......... .......... .......... .......... .......... 43% 80.0M 4s
    ##  57950K .......... .......... .......... .......... .......... 43% 24.4M 4s
    ##  58000K .......... .......... .......... .......... .......... 43% 46.2M 4s
    ##  58050K .......... .......... .......... .......... .......... 43% 16.6M 4s
    ##  58100K .......... .......... .......... .......... .......... 43% 45.7M 4s
    ##  58150K .......... .......... .......... .......... .......... 43% 16.7M 4s
    ##  58200K .......... .......... .......... .......... .......... 43% 40.0M 4s
    ##  58250K .......... .......... .......... .......... .......... 43% 16.1M 4s
    ##  58300K .......... .......... .......... .......... .......... 43% 24.6M 4s
    ##  58350K .......... .......... .......... .......... .......... 43% 25.4M 4s
    ##  58400K .......... .......... .......... .......... .......... 43% 19.0M 4s
    ##  58450K .......... .......... .......... .......... .......... 43% 50.3M 4s
    ##  58500K .......... .......... .......... .......... .......... 43% 23.3M 4s
    ##  58550K .......... .......... .......... .......... .......... 43% 25.5M 4s
    ##  58600K .......... .......... .......... .......... .......... 43% 23.4M 4s
    ##  58650K .......... .......... .......... .......... .......... 43% 19.1M 4s
    ##  58700K .......... .......... .......... .......... .......... 43% 20.9M 4s
    ##  58750K .......... .......... .......... .......... .......... 43% 17.7M 4s
    ##  58800K .......... .......... .......... .......... .......... 43% 37.5M 4s
    ##  58850K .......... .......... .......... .......... .......... 43% 18.6M 4s
    ##  58900K .......... .......... .......... .......... .......... 43% 25.9M 4s
    ##  58950K .......... .......... .......... .......... .......... 44% 20.8M 4s
    ##  59000K .......... .......... .......... .......... .......... 44% 65.3M 4s
    ##  59050K .......... .......... .......... .......... .......... 44% 14.2M 4s
    ##  59100K .......... .......... .......... .......... .......... 44% 10.1M 4s
    ##  59150K .......... .......... .......... .......... .......... 44% 73.2M 4s
    ##  59200K .......... .......... .......... .......... .......... 44% 17.8M 4s
    ##  59250K .......... .......... .......... .......... .......... 44% 52.6M 4s
    ##  59300K .......... .......... .......... .......... .......... 44% 18.3M 4s
    ##  59350K .......... .......... .......... .......... .......... 44% 32.2M 4s
    ##  59400K .......... .......... .......... .......... .......... 44% 50.5M 4s
    ##  59450K .......... .......... .......... .......... .......... 44% 23.1M 4s
    ##  59500K .......... .......... .......... .......... .......... 44% 13.5M 4s
    ##  59550K .......... .......... .......... .......... .......... 44% 18.4M 4s
    ##  59600K .......... .......... .......... .......... .......... 44% 45.2M 4s
    ##  59650K .......... .......... .......... .......... .......... 44% 23.4M 4s
    ##  59700K .......... .......... .......... .......... .......... 44% 32.3M 4s
    ##  59750K .......... .......... .......... .......... .......... 44% 16.8M 4s
    ##  59800K .......... .......... .......... .......... .......... 44% 18.1M 4s
    ##  59850K .......... .......... .......... .......... .......... 44% 40.0M 4s
    ##  59900K .......... .......... .......... .......... .......... 44% 14.0M 4s
    ##  59950K .......... .......... .......... .......... .......... 44% 53.8M 4s
    ##  60000K .......... .......... .......... .......... .......... 44% 16.3M 4s
    ##  60050K .......... .......... .......... .......... .......... 44% 17.7M 4s
    ##  60100K .......... .......... .......... .......... .......... 44% 44.0M 4s
    ##  60150K .......... .......... .......... .......... .......... 44% 43.2M 4s
    ##  60200K .......... .......... .......... .......... .......... 44% 14.3M 4s
    ##  60250K .......... .......... .......... .......... .......... 44% 51.6M 4s
    ##  60300K .......... .......... .......... .......... .......... 45% 14.2M 4s
    ##  60350K .......... .......... .......... .......... .......... 45% 28.9M 4s
    ##  60400K .......... .......... .......... .......... .......... 45% 26.1M 4s
    ##  60450K .......... .......... .......... .......... .......... 45% 19.5M 4s
    ##  60500K .......... .......... .......... .......... .......... 45% 12.4M 4s
    ##  60550K .......... .......... .......... .......... .......... 45% 31.6M 4s
    ##  60600K .......... .......... .......... .......... .......... 45% 17.7M 4s
    ##  60650K .......... .......... .......... .......... .......... 45% 45.1M 4s
    ##  60700K .......... .......... .......... .......... .......... 45%  105M 4s
    ##  60750K .......... .......... .......... .......... .......... 45% 13.8M 4s
    ##  60800K .......... .......... .......... .......... .......... 45% 51.6M 4s
    ##  60850K .......... .......... .......... .......... .......... 45% 24.1M 4s
    ##  60900K .......... .......... .......... .......... .......... 45% 15.1M 4s
    ##  60950K .......... .......... .......... .......... .......... 45% 19.8M 4s
    ##  61000K .......... .......... .......... .......... .......... 45% 31.6M 4s
    ##  61050K .......... .......... .......... .......... .......... 45% 57.1M 4s
    ##  61100K .......... .......... .......... .......... .......... 45% 12.2M 4s
    ##  61150K .......... .......... .......... .......... .......... 45% 20.1M 4s
    ##  61200K .......... .......... .......... .......... .......... 45% 30.2M 4s
    ##  61250K .......... .......... .......... .......... .......... 45% 40.1M 4s
    ##  61300K .......... .......... .......... .......... .......... 45% 20.1M 4s
    ##  61350K .......... .......... .......... .......... .......... 45% 16.5M 4s
    ##  61400K .......... .......... .......... .......... .......... 45% 51.9M 4s
    ##  61450K .......... .......... .......... .......... .......... 45% 14.2M 4s
    ##  61500K .......... .......... .......... .......... .......... 45% 18.0M 4s
    ##  61550K .......... .......... .......... .......... .......... 45% 32.6M 4s
    ##  61600K .......... .......... .......... .......... .......... 45% 16.2M 4s
    ##  61650K .......... .......... .......... .......... .......... 46% 16.4M 3s
    ##  61700K .......... .......... .......... .......... .......... 46% 53.7M 3s
    ##  61750K .......... .......... .......... .......... .......... 46% 12.1M 3s
    ##  61800K .......... .......... .......... .......... .......... 46% 73.1M 3s
    ##  61850K .......... .......... .......... .......... .......... 46% 16.6M 3s
    ##  61900K .......... .......... .......... .......... .......... 46% 30.8M 3s
    ##  61950K .......... .......... .......... .......... .......... 46% 13.5M 3s
    ##  62000K .......... .......... .......... .......... .......... 46% 12.1M 3s
    ##  62050K .......... .......... .......... .......... .......... 46% 40.5M 3s
    ##  62100K .......... .......... .......... .......... .......... 46% 23.7M 3s
    ##  62150K .......... .......... .......... .......... .......... 46% 20.0M 3s
    ##  62200K .......... .......... .......... .......... .......... 46% 16.0M 3s
    ##  62250K .......... .......... .......... .......... .......... 46% 44.4M 3s
    ##  62300K .......... .......... .......... .......... .......... 46% 16.5M 3s
    ##  62350K .......... .......... .......... .......... .......... 46% 11.5M 3s
    ##  62400K .......... .......... .......... .......... .......... 46% 81.7M 3s
    ##  62450K .......... .......... .......... .......... .......... 46% 27.8M 3s
    ##  62500K .......... .......... .......... .......... .......... 46% 15.3M 3s
    ##  62550K .......... .......... .......... .......... .......... 46% 31.8M 3s
    ##  62600K .......... .......... .......... .......... .......... 46% 17.0M 3s
    ##  62650K .......... .......... .......... .......... .......... 46% 39.2M 3s
    ##  62700K .......... .......... .......... .......... .......... 46% 20.8M 3s
    ##  62750K .......... .......... .......... .......... .......... 46% 32.1M 3s
    ##  62800K .......... .......... .......... .......... .......... 46% 18.5M 3s
    ##  62850K .......... .......... .......... .......... .......... 46% 18.6M 3s
    ##  62900K .......... .......... .......... .......... .......... 46% 20.4M 3s
    ##  62950K .......... .......... .......... .......... .......... 46% 40.5M 3s
    ##  63000K .......... .......... .......... .......... .......... 47% 18.7M 3s
    ##  63050K .......... .......... .......... .......... .......... 47% 48.3M 3s
    ##  63100K .......... .......... .......... .......... .......... 47% 17.6M 3s
    ##  63150K .......... .......... .......... .......... .......... 47% 31.2M 3s
    ##  63200K .......... .......... .......... .......... .......... 47% 19.4M 3s
    ##  63250K .......... .......... .......... .......... .......... 47% 27.3M 3s
    ##  63300K .......... .......... .......... .......... .......... 47% 24.7M 3s
    ##  63350K .......... .......... .......... .......... .......... 47% 19.7M 3s
    ##  63400K .......... .......... .......... .......... .......... 47% 72.0M 3s
    ##  63450K .......... .......... .......... .......... .......... 47% 20.0M 3s
    ##  63500K .......... .......... .......... .......... .......... 47% 26.8M 3s
    ##  63550K .......... .......... .......... .......... .......... 47% 18.1M 3s
    ##  63600K .......... .......... .......... .......... .......... 47% 32.3M 3s
    ##  63650K .......... .......... .......... .......... .......... 47% 23.6M 3s
    ##  63700K .......... .......... .......... .......... .......... 47% 28.4M 3s
    ##  63750K .......... .......... .......... .......... .......... 47% 16.4M 3s
    ##  63800K .......... .......... .......... .......... .......... 47% 37.9M 3s
    ##  63850K .......... .......... .......... .......... .......... 47% 44.7M 3s
    ##  63900K .......... .......... .......... .......... .......... 47% 17.2M 3s
    ##  63950K .......... .......... .......... .......... .......... 47% 39.3M 3s
    ##  64000K .......... .......... .......... .......... .......... 47% 17.2M 3s
    ##  64050K .......... .......... .......... .......... .......... 47% 34.7M 3s
    ##  64100K .......... .......... .......... .......... .......... 47% 18.2M 3s
    ##  64150K .......... .......... .......... .......... .......... 47% 27.6M 3s
    ##  64200K .......... .......... .......... .......... .......... 47% 17.5M 3s
    ##  64250K .......... .......... .......... .......... .......... 47% 54.9M 3s
    ##  64300K .......... .......... .......... .......... .......... 47% 33.4M 3s
    ##  64350K .......... .......... .......... .......... .......... 48% 20.4M 3s
    ##  64400K .......... .......... .......... .......... .......... 48% 15.7M 3s
    ##  64450K .......... .......... .......... .......... .......... 48% 52.5M 3s
    ##  64500K .......... .......... .......... .......... .......... 48% 16.7M 3s
    ##  64550K .......... .......... .......... .......... .......... 48% 31.9M 3s
    ##  64600K .......... .......... .......... .......... .......... 48% 17.8M 3s
    ##  64650K .......... .......... .......... .......... .......... 48% 13.5M 3s
    ##  64700K .......... .......... .......... .......... .......... 48% 77.9M 3s
    ##  64750K .......... .......... .......... .......... .......... 48% 13.7M 3s
    ##  64800K .......... .......... .......... .......... .......... 48% 13.8M 3s
    ##  64850K .......... .......... .......... .......... .......... 48% 44.1M 3s
    ##  64900K .......... .......... .......... .......... .......... 48% 14.1M 3s
    ##  64950K .......... .......... .......... .......... .......... 48% 47.5M 3s
    ##  65000K .......... .......... .......... .......... .......... 48% 10.3M 3s
    ##  65050K .......... .......... .......... .......... .......... 48% 62.3M 3s
    ##  65100K .......... .......... .......... .......... .......... 48% 14.1M 3s
    ##  65150K .......... .......... .......... .......... .......... 48% 39.2M 3s
    ##  65200K .......... .......... .......... .......... .......... 48% 17.2M 3s
    ##  65250K .......... .......... .......... .......... .......... 48% 38.0M 3s
    ##  65300K .......... .......... .......... .......... .......... 48% 54.4M 3s
    ##  65350K .......... .......... .......... .......... .......... 48% 19.7M 3s
    ##  65400K .......... .......... .......... .......... .......... 48% 12.2M 3s
    ##  65450K .......... .......... .......... .......... .......... 48%  106M 3s
    ##  65500K .......... .......... .......... .......... .......... 48% 9.68M 3s
    ##  65550K .......... .......... .......... .......... .......... 48% 66.6M 3s
    ##  65600K .......... .......... .......... .......... .......... 48% 16.2M 3s
    ##  65650K .......... .......... .......... .......... .......... 49% 52.3M 3s
    ##  65700K .......... .......... .......... .......... .......... 49% 13.1M 3s
    ##  65750K .......... .......... .......... .......... .......... 49% 18.3M 3s
    ##  65800K .......... .......... .......... .......... .......... 49% 26.6M 3s
    ##  65850K .......... .......... .......... .......... .......... 49% 99.9M 3s
    ##  65900K .......... .......... .......... .......... .......... 49% 11.2M 3s
    ##  65950K .......... .......... .......... .......... .......... 49% 15.4M 3s
    ##  66000K .......... .......... .......... .......... .......... 49% 22.0M 3s
    ##  66050K .......... .......... .......... .......... .......... 49% 96.7M 3s
    ##  66100K .......... .......... .......... .......... .......... 49% 19.8M 3s
    ##  66150K .......... .......... .......... .......... .......... 49% 27.4M 3s
    ##  66200K .......... .......... .......... .......... .......... 49% 20.8M 3s
    ##  66250K .......... .......... .......... .......... .......... 49% 21.6M 3s
    ##  66300K .......... .......... .......... .......... .......... 49% 58.6M 3s
    ##  66350K .......... .......... .......... .......... .......... 49% 16.2M 3s
    ##  66400K .......... .......... .......... .......... .......... 49% 16.1M 3s
    ##  66450K .......... .......... .......... .......... .......... 49% 41.8M 3s
    ##  66500K .......... .......... .......... .......... .......... 49% 50.9M 3s
    ##  66550K .......... .......... .......... .......... .......... 49% 24.3M 3s
    ##  66600K .......... .......... .......... .......... .......... 49% 17.9M 3s
    ##  66650K .......... .......... .......... .......... .......... 49% 49.4M 3s
    ##  66700K .......... .......... .......... .......... .......... 49% 12.8M 3s
    ##  66750K .......... .......... .......... .......... .......... 49% 23.8M 3s
    ##  66800K .......... .......... .......... .......... .......... 49% 13.7M 3s
    ##  66850K .......... .......... .......... .......... .......... 49% 30.9M 3s
    ##  66900K .......... .......... .......... .......... .......... 49% 17.9M 3s
    ##  66950K .......... .......... .......... .......... .......... 49% 31.4M 3s
    ##  67000K .......... .......... .......... .......... .......... 50% 16.5M 3s
    ##  67050K .......... .......... .......... .......... .......... 50% 79.4M 3s
    ##  67100K .......... .......... .......... .......... .......... 50% 11.2M 3s
    ##  67150K .......... .......... .......... .......... .......... 50% 12.3M 3s
    ##  67200K .......... .......... .......... .......... .......... 50% 30.3M 3s
    ##  67250K .......... .......... .......... .......... .......... 50% 72.2M 3s
    ##  67300K .......... .......... .......... .......... .......... 50% 18.2M 3s
    ##  67350K .......... .......... .......... .......... .......... 50% 42.6M 3s
    ##  67400K .......... .......... .......... .......... .......... 50% 18.3M 3s
    ##  67450K .......... .......... .......... .......... .......... 50% 23.8M 3s
    ##  67500K .......... .......... .......... .......... .......... 50% 31.0M 3s
    ##  67550K .......... .......... .......... .......... .......... 50% 20.7M 3s
    ##  67600K .......... .......... .......... .......... .......... 50% 33.7M 3s
    ##  67650K .......... .......... .......... .......... .......... 50% 19.5M 3s
    ##  67700K .......... .......... .......... .......... .......... 50% 32.5M 3s
    ##  67750K .......... .......... .......... .......... .......... 50% 18.5M 3s
    ##  67800K .......... .......... .......... .......... .......... 50% 23.5M 3s
    ##  67850K .......... .......... .......... .......... .......... 50% 29.0M 3s
    ##  67900K .......... .......... .......... .......... .......... 50% 15.5M 3s
    ##  67950K .......... .......... .......... .......... .......... 50% 25.2M 3s
    ##  68000K .......... .......... .......... .......... .......... 50% 43.2M 3s
    ##  68050K .......... .......... .......... .......... .......... 50% 24.2M 3s
    ##  68100K .......... .......... .......... .......... .......... 50% 23.8M 3s
    ##  68150K .......... .......... .......... .......... .......... 50% 22.9M 3s
    ##  68200K .......... .......... .......... .......... .......... 50% 50.8M 3s
    ##  68250K .......... .......... .......... .......... .......... 50% 15.3M 3s
    ##  68300K .......... .......... .......... .......... .......... 50% 25.5M 3s
    ##  68350K .......... .......... .......... .......... .......... 51% 15.3M 3s
    ##  68400K .......... .......... .......... .......... .......... 51% 85.4M 3s
    ##  68450K .......... .......... .......... .......... .......... 51% 16.4M 3s
    ##  68500K .......... .......... .......... .......... .......... 51% 50.9M 3s
    ##  68550K .......... .......... .......... .......... .......... 51% 12.9M 3s
    ##  68600K .......... .......... .......... .......... .......... 51% 87.3M 3s
    ##  68650K .......... .......... .......... .......... .......... 51% 14.4M 3s
    ##  68700K .......... .......... .......... .......... .......... 51% 15.8M 3s
    ##  68750K .......... .......... .......... .......... .......... 51% 37.2M 3s
    ##  68800K .......... .......... .......... .......... .......... 51% 19.6M 3s
    ##  68850K .......... .......... .......... .......... .......... 51% 61.1M 3s
    ##  68900K .......... .......... .......... .......... .......... 51% 38.1M 3s
    ##  68950K .......... .......... .......... .......... .......... 51% 19.6M 3s
    ##  69000K .......... .......... .......... .......... .......... 51% 27.8M 3s
    ##  69050K .......... .......... .......... .......... .......... 51% 20.0M 3s
    ##  69100K .......... .......... .......... .......... .......... 51% 31.6M 3s
    ##  69150K .......... .......... .......... .......... .......... 51% 20.3M 3s
    ##  69200K .......... .......... .......... .......... .......... 51% 28.0M 3s
    ##  69250K .......... .......... .......... .......... .......... 51% 23.2M 3s
    ##  69300K .......... .......... .......... .......... .......... 51% 28.0M 3s
    ##  69350K .......... .......... .......... .......... .......... 51% 16.5M 3s
    ##  69400K .......... .......... .......... .......... .......... 51% 15.1M 3s
    ##  69450K .......... .......... .......... .......... .......... 51% 22.0M 3s
    ##  69500K .......... .......... .......... .......... .......... 51% 97.9M 3s
    ##  69550K .......... .......... .......... .......... .......... 51% 20.0M 3s
    ##  69600K .......... .......... .......... .......... .......... 51% 25.0M 3s
    ##  69650K .......... .......... .......... .......... .......... 51% 24.3M 3s
    ##  69700K .......... .......... .......... .......... .......... 52% 18.4M 3s
    ##  69750K .......... .......... .......... .......... .......... 52% 26.7M 3s
    ##  69800K .......... .......... .......... .......... .......... 52% 32.7M 3s
    ##  69850K .......... .......... .......... .......... .......... 52% 23.0M 3s
    ##  69900K .......... .......... .......... .......... .......... 52% 22.2M 3s
    ##  69950K .......... .......... .......... .......... .......... 52% 22.2M 3s
    ##  70000K .......... .......... .......... .......... .......... 52% 20.9M 3s
    ##  70050K .......... .......... .......... .......... .......... 52% 35.0M 3s
    ##  70100K .......... .......... .......... .......... .......... 52% 26.4M 3s
    ##  70150K .......... .......... .......... .......... .......... 52% 15.8M 3s
    ##  70200K .......... .......... .......... .......... .......... 52% 33.5M 3s
    ##  70250K .......... .......... .......... .......... .......... 52% 18.8M 3s
    ##  70300K .......... .......... .......... .......... .......... 52% 31.4M 3s
    ##  70350K .......... .......... .......... .......... .......... 52% 19.8M 3s
    ##  70400K .......... .......... .......... .......... .......... 52% 24.4M 3s
    ##  70450K .......... .......... .......... .......... .......... 52% 31.7M 3s
    ##  70500K .......... .......... .......... .......... .......... 52% 19.4M 3s
    ##  70550K .......... .......... .......... .......... .......... 52% 20.0M 3s
    ##  70600K .......... .......... .......... .......... .......... 52% 47.2M 3s
    ##  70650K .......... .......... .......... .......... .......... 52% 16.0M 3s
    ##  70700K .......... .......... .......... .......... .......... 52% 46.3M 3s
    ##  70750K .......... .......... .......... .......... .......... 52% 15.2M 3s
    ##  70800K .......... .......... .......... .......... .......... 52% 42.5M 3s
    ##  70850K .......... .......... .......... .......... .......... 52% 7.15M 3s
    ##  70900K .......... .......... .......... .......... .......... 52% 35.3M 3s
    ##  70950K .......... .......... .......... .......... .......... 52% 12.4M 3s
    ##  71000K .......... .......... .......... .......... .......... 52%  101M 3s
    ##  71050K .......... .......... .......... .......... .......... 53% 25.8M 3s
    ##  71100K .......... .......... .......... .......... .......... 53% 12.1M 3s
    ##  71150K .......... .......... .......... .......... .......... 53% 18.1M 3s
    ##  71200K .......... .......... .......... .......... .......... 53% 23.2M 3s
    ##  71250K .......... .......... .......... .......... .......... 53% 96.6M 3s
    ##  71300K .......... .......... .......... .......... .......... 53% 15.0M 3s
    ##  71350K .......... .......... .......... .......... .......... 53% 15.1M 3s
    ##  71400K .......... .......... .......... .......... .......... 53% 44.7M 3s
    ##  71450K .......... .......... .......... .......... .......... 53% 15.8M 3s
    ##  71500K .......... .......... .......... .......... .......... 53% 41.9M 3s
    ##  71550K .......... .......... .......... .......... .......... 53% 16.3M 3s
    ##  71600K .......... .......... .......... .......... .......... 53% 31.5M 3s
    ##  71650K .......... .......... .......... .......... .......... 53% 12.6M 3s
    ##  71700K .......... .......... .......... .......... .......... 53% 19.1M 3s
    ##  71750K .......... .......... .......... .......... .......... 53% 30.1M 3s
    ##  71800K .......... .......... .......... .......... .......... 53% 17.0M 3s
    ##  71850K .......... .......... .......... .......... .......... 53% 33.7M 3s
    ##  71900K .......... .......... .......... .......... .......... 53% 46.1M 3s
    ##  71950K .......... .......... .......... .......... .......... 53% 13.7M 3s
    ##  72000K .......... .......... .......... .......... .......... 53% 24.6M 3s
    ##  72050K .......... .......... .......... .......... .......... 53% 17.9M 3s
    ##  72100K .......... .......... .......... .......... .......... 53% 69.3M 3s
    ##  72150K .......... .......... .......... .......... .......... 53% 11.5M 3s
    ##  72200K .......... .......... .......... .......... .......... 53% 46.1M 3s
    ##  72250K .......... .......... .......... .......... .......... 53% 22.1M 3s
    ##  72300K .......... .......... .......... .......... .......... 53% 17.1M 3s
    ##  72350K .......... .......... .......... .......... .......... 54% 14.5M 3s
    ##  72400K .......... .......... .......... .......... .......... 54% 34.4M 3s
    ##  72450K .......... .......... .......... .......... .......... 54% 42.8M 3s
    ##  72500K .......... .......... .......... .......... .......... 54% 19.4M 3s
    ##  72550K .......... .......... .......... .......... .......... 54% 25.5M 3s
    ##  72600K .......... .......... .......... .......... .......... 54% 24.5M 3s
    ##  72650K .......... .......... .......... .......... .......... 54% 22.8M 3s
    ##  72700K .......... .......... .......... .......... .......... 54% 57.0M 3s
    ##  72750K .......... .......... .......... .......... .......... 54% 16.4M 3s
    ##  72800K .......... .......... .......... .......... .......... 54% 20.3M 3s
    ##  72850K .......... .......... .......... .......... .......... 54% 33.0M 3s
    ##  72900K .......... .......... .......... .......... .......... 54% 34.7M 3s
    ##  72950K .......... .......... .......... .......... .......... 54% 29.0M 3s
    ##  73000K .......... .......... .......... .......... .......... 54% 14.4M 3s
    ##  73050K .......... .......... .......... .......... .......... 54% 31.8M 3s
    ##  73100K .......... .......... .......... .......... .......... 54% 15.2M 3s
    ##  73150K .......... .......... .......... .......... .......... 54% 39.9M 3s
    ##  73200K .......... .......... .......... .......... .......... 54% 14.9M 3s
    ##  73250K .......... .......... .......... .......... .......... 54% 51.1M 3s
    ##  73300K .......... .......... .......... .......... .......... 54% 17.5M 3s
    ##  73350K .......... .......... .......... .......... .......... 54% 34.2M 3s
    ##  73400K .......... .......... .......... .......... .......... 54% 46.8M 3s
    ##  73450K .......... .......... .......... .......... .......... 54% 20.8M 3s
    ##  73500K .......... .......... .......... .......... .......... 54% 24.5M 3s
    ##  73550K .......... .......... .......... .......... .......... 54% 15.1M 3s
    ##  73600K .......... .......... .......... .......... .......... 54% 33.4M 3s
    ##  73650K .......... .......... .......... .......... .......... 54% 42.6M 3s
    ##  73700K .......... .......... .......... .......... .......... 55% 14.7M 3s
    ##  73750K .......... .......... .......... .......... .......... 55% 19.7M 3s
    ##  73800K .......... .......... .......... .......... .......... 55% 36.8M 3s
    ##  73850K .......... .......... .......... .......... .......... 55% 21.5M 3s
    ##  73900K .......... .......... .......... .......... .......... 55% 24.0M 3s
    ##  73950K .......... .......... .......... .......... .......... 55% 25.4M 3s
    ##  74000K .......... .......... .......... .......... .......... 55% 13.0M 3s
    ##  74050K .......... .......... .......... .......... .......... 55% 46.2M 3s
    ##  74100K .......... .......... .......... .......... .......... 55% 18.5M 3s
    ##  74150K .......... .......... .......... .......... .......... 55% 16.6M 3s
    ##  74200K .......... .......... .......... .......... .......... 55% 36.8M 3s
    ##  74250K .......... .......... .......... .......... .......... 55% 14.1M 3s
    ##  74300K .......... .......... .......... .......... .......... 55% 87.3M 3s
    ##  74350K .......... .......... .......... .......... .......... 55% 12.4M 3s
    ##  74400K .......... .......... .......... .......... .......... 55% 11.4M 3s
    ##  74450K .......... .......... .......... .......... .......... 55% 51.4M 3s
    ##  74500K .......... .......... .......... .......... .......... 55% 20.2M 3s
    ##  74550K .......... .......... .......... .......... .......... 55% 23.6M 3s
    ##  74600K .......... .......... .......... .......... .......... 55% 28.1M 3s
    ##  74650K .......... .......... .......... .......... .......... 55% 27.2M 3s
    ##  74700K .......... .......... .......... .......... .......... 55% 38.2M 3s
    ##  74750K .......... .......... .......... .......... .......... 55% 18.2M 3s
    ##  74800K .......... .......... .......... .......... .......... 55% 27.5M 3s
    ##  74850K .......... .......... .......... .......... .......... 55% 20.2M 3s
    ##  74900K .......... .......... .......... .......... .......... 55% 25.7M 3s
    ##  74950K .......... .......... .......... .......... .......... 55% 23.2M 3s
    ##  75000K .......... .......... .......... .......... .......... 55% 16.4M 3s
    ##  75050K .......... .......... .......... .......... .......... 56% 49.7M 3s
    ##  75100K .......... .......... .......... .......... .......... 56% 17.3M 3s
    ##  75150K .......... .......... .......... .......... .......... 56% 34.5M 3s
    ##  75200K .......... .......... .......... .......... .......... 56% 24.0M 3s
    ##  75250K .......... .......... .......... .......... .......... 56% 15.3M 3s
    ##  75300K .......... .......... .......... .......... .......... 56% 37.3M 3s
    ##  75350K .......... .......... .......... .......... .......... 56% 14.9M 3s
    ##  75400K .......... .......... .......... .......... .......... 56% 44.5M 3s
    ##  75450K .......... .......... .......... .......... .......... 56% 25.9M 3s
    ##  75500K .......... .......... .......... .......... .......... 56% 33.0M 3s
    ##  75550K .......... .......... .......... .......... .......... 56% 16.4M 3s
    ##  75600K .......... .......... .......... .......... .......... 56% 31.9M 3s
    ##  75650K .......... .......... .......... .......... .......... 56% 29.0M 3s
    ##  75700K .......... .......... .......... .......... .......... 56% 16.2M 3s
    ##  75750K .......... .......... .......... .......... .......... 56% 37.2M 3s
    ##  75800K .......... .......... .......... .......... .......... 56% 11.6M 3s
    ##  75850K .......... .......... .......... .......... .......... 56% 43.2M 3s
    ##  75900K .......... .......... .......... .......... .......... 56% 24.7M 3s
    ##  75950K .......... .......... .......... .......... .......... 56% 18.9M 3s
    ##  76000K .......... .......... .......... .......... .......... 56% 23.7M 3s
    ##  76050K .......... .......... .......... .......... .......... 56% 21.9M 3s
    ##  76100K .......... .......... .......... .......... .......... 56% 37.5M 3s
    ##  76150K .......... .......... .......... .......... .......... 56% 28.1M 3s
    ##  76200K .......... .......... .......... .......... .......... 56% 20.4M 3s
    ##  76250K .......... .......... .......... .......... .......... 56% 16.1M 3s
    ##  76300K .......... .......... .......... .......... .......... 56% 48.7M 3s
    ##  76350K .......... .......... .......... .......... .......... 56% 11.5M 3s
    ##  76400K .......... .......... .......... .......... .......... 57% 81.6M 3s
    ##  76450K .......... .......... .......... .......... .......... 57% 23.9M 3s
    ##  76500K .......... .......... .......... .......... .......... 57% 14.9M 3s
    ##  76550K .......... .......... .......... .......... .......... 57% 45.8M 3s
    ##  76600K .......... .......... .......... .......... .......... 57% 10.8M 3s
    ##  76650K .......... .......... .......... .......... .......... 57% 78.8M 3s
    ##  76700K .......... .......... .......... .......... .......... 57% 15.3M 3s
    ##  76750K .......... .......... .......... .......... .......... 57% 36.3M 3s
    ##  76800K .......... .......... .......... .......... .......... 57% 17.8M 3s
    ##  76850K .......... .......... .......... .......... .......... 57% 45.2M 3s
    ##  76900K .......... .......... .......... .......... .......... 57% 13.4M 3s
    ##  76950K .......... .......... .......... .......... .......... 57% 22.6M 3s
    ##  77000K .......... .......... .......... .......... .......... 57% 48.2M 3s
    ##  77050K .......... .......... .......... .......... .......... 57% 15.6M 3s
    ##  77100K .......... .......... .......... .......... .......... 57% 14.8M 3s
    ##  77150K .......... .......... .......... .......... .......... 57% 38.0M 3s
    ##  77200K .......... .......... .......... .......... .......... 57% 22.5M 3s
    ##  77250K .......... .......... .......... .......... .......... 57% 14.5M 3s
    ##  77300K .......... .......... .......... .......... .......... 57% 50.6M 3s
    ##  77350K .......... .......... .......... .......... .......... 57% 13.6M 3s
    ##  77400K .......... .......... .......... .......... .......... 57% 13.0M 3s
    ##  77450K .......... .......... .......... .......... .......... 57% 42.7M 3s
    ##  77500K .......... .......... .......... .......... .......... 57% 22.6M 3s
    ##  77550K .......... .......... .......... .......... .......... 57% 16.3M 3s
    ##  77600K .......... .......... .......... .......... .......... 57% 31.0M 3s
    ##  77650K .......... .......... .......... .......... .......... 57% 12.6M 3s
    ##  77700K .......... .......... .......... .......... .......... 57% 86.5M 3s
    ##  77750K .......... .......... .......... .......... .......... 58% 25.5M 3s
    ##  77800K .......... .......... .......... .......... .......... 58% 15.7M 3s
    ##  77850K .......... .......... .......... .......... .......... 58% 45.7M 3s
    ##  77900K .......... .......... .......... .......... .......... 58% 13.0M 3s
    ##  77950K .......... .......... .......... .......... .......... 58% 22.0M 3s
    ##  78000K .......... .......... .......... .......... .......... 58% 18.2M 3s
    ##  78050K .......... .......... .......... .......... .......... 58% 27.3M 3s
    ##  78100K .......... .......... .......... .......... .......... 58% 21.0M 3s
    ##  78150K .......... .......... .......... .......... .......... 58% 33.6M 3s
    ##  78200K .......... .......... .......... .......... .......... 58% 52.2M 3s
    ##  78250K .......... .......... .......... .......... .......... 58% 11.8M 3s
    ##  78300K .......... .......... .......... .......... .......... 58% 76.5M 3s
    ##  78350K .......... .......... .......... .......... .......... 58% 17.5M 3s
    ##  78400K .......... .......... .......... .......... .......... 58% 13.5M 3s
    ##  78450K .......... .......... .......... .......... .......... 58% 45.8M 3s
    ##  78500K .......... .......... .......... .......... .......... 58% 19.5M 3s
    ##  78550K .......... .......... .......... .......... .......... 58% 45.1M 3s
    ##  78600K .......... .......... .......... .......... .......... 58% 16.0M 3s
    ##  78650K .......... .......... .......... .......... .......... 58% 29.5M 3s
    ##  78700K .......... .......... .......... .......... .......... 58% 54.8M 3s
    ##  78750K .......... .......... .......... .......... .......... 58% 22.6M 3s
    ##  78800K .......... .......... .......... .......... .......... 58% 26.5M 3s
    ##  78850K .......... .......... .......... .......... .......... 58% 21.8M 3s
    ##  78900K .......... .......... .......... .......... .......... 58% 21.7M 3s
    ##  78950K .......... .......... .......... .......... .......... 58% 38.3M 3s
    ##  79000K .......... .......... .......... .......... .......... 58% 18.0M 3s
    ##  79050K .......... .......... .......... .......... .......... 59% 40.3M 3s
    ##  79100K .......... .......... .......... .......... .......... 59% 16.0M 3s
    ##  79150K .......... .......... .......... .......... .......... 59% 34.8M 3s
    ##  79200K .......... .......... .......... .......... .......... 59% 15.3M 3s
    ##  79250K .......... .......... .......... .......... .......... 59% 51.4M 3s
    ##  79300K .......... .......... .......... .......... .......... 59% 14.3M 3s
    ##  79350K .......... .......... .......... .......... .......... 59% 37.7M 3s
    ##  79400K .......... .......... .......... .......... .......... 59% 15.4M 3s
    ##  79450K .......... .......... .......... .......... .......... 59% 44.0M 3s
    ##  79500K .......... .......... .......... .......... .......... 59% 15.6M 3s
    ##  79550K .......... .......... .......... .......... .......... 59% 17.0M 3s
    ##  79600K .......... .......... .......... .......... .......... 59% 45.3M 3s
    ##  79650K .......... .......... .......... .......... .......... 59% 16.1M 3s
    ##  79700K .......... .......... .......... .......... .......... 59% 39.6M 3s
    ##  79750K .......... .......... .......... .......... .......... 59% 16.0M 3s
    ##  79800K .......... .......... .......... .......... .......... 59% 42.8M 3s
    ##  79850K .......... .......... .......... .......... .......... 59% 12.8M 3s
    ##  79900K .......... .......... .......... .......... .......... 59% 64.3M 3s
    ##  79950K .......... .......... .......... .......... .......... 59% 16.6M 3s
    ##  80000K .......... .......... .......... .......... .......... 59% 31.6M 3s
    ##  80050K .......... .......... .......... .......... .......... 59% 16.2M 3s
    ##  80100K .......... .......... .......... .......... .......... 59% 36.3M 3s
    ##  80150K .......... .......... .......... .......... .......... 59% 18.4M 3s
    ##  80200K .......... .......... .......... .......... .......... 59% 31.1M 3s
    ##  80250K .......... .......... .......... .......... .......... 59% 43.0M 3s
    ##  80300K .......... .......... .......... .......... .......... 59% 27.6M 3s
    ##  80350K .......... .......... .......... .......... .......... 59% 14.6M 3s
    ##  80400K .......... .......... .......... .......... .......... 60% 14.9M 3s
    ##  80450K .......... .......... .......... .......... .......... 60% 37.1M 3s
    ##  80500K .......... .......... .......... .......... .......... 60% 20.1M 3s
    ##  80550K .......... .......... .......... .......... .......... 60% 17.0M 3s
    ##  80600K .......... .......... .......... .......... .......... 60% 16.7M 3s
    ##  80650K .......... .......... .......... .......... .......... 60% 14.2M 3s
    ##  80700K .......... .......... .......... .......... .......... 60% 22.3M 3s
    ##  80750K .......... .......... .......... .......... .......... 60% 23.6M 3s
    ##  80800K .......... .......... .......... .......... .......... 60% 47.6M 3s
    ##  80850K .......... .......... .......... .......... .......... 60% 14.4M 3s
    ##  80900K .......... .......... .......... .......... .......... 60% 44.2M 3s
    ##  80950K .......... .......... .......... .......... .......... 60% 13.9M 3s
    ##  81000K .......... .......... .......... .......... .......... 60% 42.8M 3s
    ##  81050K .......... .......... .......... .......... .......... 60% 18.9M 3s
    ##  81100K .......... .......... .......... .......... .......... 60% 36.3M 2s
    ##  81150K .......... .......... .......... .......... .......... 60% 17.3M 2s
    ##  81200K .......... .......... .......... .......... .......... 60% 24.2M 2s
    ##  81250K .......... .......... .......... .......... .......... 60% 43.5M 2s
    ##  81300K .......... .......... .......... .......... .......... 60% 22.4M 2s
    ##  81350K .......... .......... .......... .......... .......... 60% 17.9M 2s
    ##  81400K .......... .......... .......... .......... .......... 60% 22.5M 2s
    ##  81450K .......... .......... .......... .......... .......... 60% 18.4M 2s
    ##  81500K .......... .......... .......... .......... .......... 60% 24.6M 2s
    ##  81550K .......... .......... .......... .......... .......... 60% 16.0M 2s
    ##  81600K .......... .......... .......... .......... .......... 60% 44.4M 2s
    ##  81650K .......... .......... .......... .......... .......... 60% 15.8M 2s
    ##  81700K .......... .......... .......... .......... .......... 60% 32.4M 2s
    ##  81750K .......... .......... .......... .......... .......... 61% 20.9M 2s
    ##  81800K .......... .......... .......... .......... .......... 61% 29.5M 2s
    ##  81850K .......... .......... .......... .......... .......... 61% 23.2M 2s
    ##  81900K .......... .......... .......... .......... .......... 61% 15.3M 2s
    ##  81950K .......... .......... .......... .......... .......... 61% 44.9M 2s
    ##  82000K .......... .......... .......... .......... .......... 61% 19.8M 2s
    ##  82050K .......... .......... .......... .......... .......... 61% 21.6M 2s
    ##  82100K .......... .......... .......... .......... .......... 61% 26.4M 2s
    ##  82150K .......... .......... .......... .......... .......... 61% 16.4M 2s
    ##  82200K .......... .......... .......... .......... .......... 61% 51.9M 2s
    ##  82250K .......... .......... .......... .......... .......... 61% 37.9M 2s
    ##  82300K .......... .......... .......... .......... .......... 61% 20.9M 2s
    ##  82350K .......... .......... .......... .......... .......... 61% 13.9M 2s
    ##  82400K .......... .......... .......... .......... .......... 61% 42.0M 2s
    ##  82450K .......... .......... .......... .......... .......... 61% 28.5M 2s
    ##  82500K .......... .......... .......... .......... .......... 61% 25.7M 2s
    ##  82550K .......... .......... .......... .......... .......... 61% 25.0M 2s
    ##  82600K .......... .......... .......... .......... .......... 61% 14.8M 2s
    ##  82650K .......... .......... .......... .......... .......... 61% 81.0M 2s
    ##  82700K .......... .......... .......... .......... .......... 61% 20.9M 2s
    ##  82750K .......... .......... .......... .......... .......... 61% 15.8M 2s
    ##  82800K .......... .......... .......... .......... .......... 61% 24.4M 2s
    ##  82850K .......... .......... .......... .......... .......... 61% 21.5M 2s
    ##  82900K .......... .......... .......... .......... .......... 61% 33.9M 2s
    ##  82950K .......... .......... .......... .......... .......... 61% 15.1M 2s
    ##  83000K .......... .......... .......... .......... .......... 61% 39.9M 2s
    ##  83050K .......... .......... .......... .......... .......... 61% 15.5M 2s
    ##  83100K .......... .......... .......... .......... .......... 62% 40.6M 2s
    ##  83150K .......... .......... .......... .......... .......... 62% 17.1M 2s
    ##  83200K .......... .......... .......... .......... .......... 62% 32.9M 2s
    ##  83250K .......... .......... .......... .......... .......... 62% 14.0M 2s
    ##  83300K .......... .......... .......... .......... .......... 62% 36.1M 2s
    ##  83350K .......... .......... .......... .......... .......... 62% 46.2M 2s
    ##  83400K .......... .......... .......... .......... .......... 62% 23.9M 2s
    ##  83450K .......... .......... .......... .......... .......... 62% 24.8M 2s
    ##  83500K .......... .......... .......... .......... .......... 62% 23.7M 2s
    ##  83550K .......... .......... .......... .......... .......... 62% 16.8M 2s
    ##  83600K .......... .......... .......... .......... .......... 62% 61.0M 2s
    ##  83650K .......... .......... .......... .......... .......... 62% 16.2M 2s
    ##  83700K .......... .......... .......... .......... .......... 62% 14.4M 2s
    ##  83750K .......... .......... .......... .......... .......... 62% 64.2M 2s
    ##  83800K .......... .......... .......... .......... .......... 62% 11.4M 2s
    ##  83850K .......... .......... .......... .......... .......... 62% 69.0M 2s
    ##  83900K .......... .......... .......... .......... .......... 62% 22.2M 2s
    ##  83950K .......... .......... .......... .......... .......... 62% 24.0M 2s
    ##  84000K .......... .......... .......... .......... .......... 62% 20.8M 2s
    ##  84050K .......... .......... .......... .......... .......... 62% 35.5M 2s
    ##  84100K .......... .......... .......... .......... .......... 62% 38.6M 2s
    ##  84150K .......... .......... .......... .......... .......... 62% 18.8M 2s
    ##  84200K .......... .......... .......... .......... .......... 62% 30.6M 2s
    ##  84250K .......... .......... .......... .......... .......... 62% 23.7M 2s
    ##  84300K .......... .......... .......... .......... .......... 62% 21.5M 2s
    ##  84350K .......... .......... .......... .......... .......... 62% 27.7M 2s
    ##  84400K .......... .......... .......... .......... .......... 62% 16.4M 2s
    ##  84450K .......... .......... .......... .......... .......... 63% 40.7M 2s
    ##  84500K .......... .......... .......... .......... .......... 63% 16.7M 2s
    ##  84550K .......... .......... .......... .......... .......... 63% 26.8M 2s
    ##  84600K .......... .......... .......... .......... .......... 63% 24.3M 2s
    ##  84650K .......... .......... .......... .......... .......... 63% 19.0M 2s
    ##  84700K .......... .......... .......... .......... .......... 63% 21.2M 2s
    ##  84750K .......... .......... .......... .......... .......... 63% 25.5M 2s
    ##  84800K .......... .......... .......... .......... .......... 63% 39.1M 2s
    ##  84850K .......... .......... .......... .......... .......... 63% 24.3M 2s
    ##  84900K .......... .......... .......... .......... .......... 63% 20.0M 2s
    ##  84950K .......... .......... .......... .......... .......... 63% 16.7M 2s
    ##  85000K .......... .......... .......... .......... .......... 63% 30.5M 2s
    ##  85050K .......... .......... .......... .......... .......... 63% 72.1M 2s
    ##  85100K .......... .......... .......... .......... .......... 63% 19.7M 2s
    ##  85150K .......... .......... .......... .......... .......... 63% 15.2M 2s
    ##  85200K .......... .......... .......... .......... .......... 63% 17.8M 2s
    ##  85250K .......... .......... .......... .......... .......... 63% 17.8M 2s
    ##  85300K .......... .......... .......... .......... .......... 63% 39.7M 2s
    ##  85350K .......... .......... .......... .......... .......... 63% 32.3M 2s
    ##  85400K .......... .......... .......... .......... .......... 63% 24.4M 2s
    ##  85450K .......... .......... .......... .......... .......... 63% 23.2M 2s
    ##  85500K .......... .......... .......... .......... .......... 63% 28.9M 2s
    ##  85550K .......... .......... .......... .......... .......... 63% 16.1M 2s
    ##  85600K .......... .......... .......... .......... .......... 63% 15.3M 2s
    ##  85650K .......... .......... .......... .......... .......... 63% 43.8M 2s
    ##  85700K .......... .......... .......... .......... .......... 63% 19.3M 2s
    ##  85750K .......... .......... .......... .......... .......... 63% 22.9M 2s
    ##  85800K .......... .......... .......... .......... .......... 64% 27.2M 2s
    ##  85850K .......... .......... .......... .......... .......... 64% 21.0M 2s
    ##  85900K .......... .......... .......... .......... .......... 64% 21.9M 2s
    ##  85950K .......... .......... .......... .......... .......... 64% 21.5M 2s
    ##  86000K .......... .......... .......... .......... .......... 64% 30.2M 2s
    ##  86050K .......... .......... .......... .......... .......... 64% 23.8M 2s
    ##  86100K .......... .......... .......... .......... .......... 64% 15.8M 2s
    ##  86150K .......... .......... .......... .......... .......... 64% 63.7M 2s
    ##  86200K .......... .......... .......... .......... .......... 64% 23.4M 2s
    ##  86250K .......... .......... .......... .......... .......... 64% 34.6M 2s
    ##  86300K .......... .......... .......... .......... .......... 64% 16.6M 2s
    ##  86350K .......... .......... .......... .......... .......... 64% 33.5M 2s
    ##  86400K .......... .......... .......... .......... .......... 64% 19.3M 2s
    ##  86450K .......... .......... .......... .......... .......... 64% 28.6M 2s
    ##  86500K .......... .......... .......... .......... .......... 64% 22.1M 2s
    ##  86550K .......... .......... .......... .......... .......... 64% 41.9M 2s
    ##  86600K .......... .......... .......... .......... .......... 64% 13.9M 2s
    ##  86650K .......... .......... .......... .......... .......... 64% 20.1M 2s
    ##  86700K .......... .......... .......... .......... .......... 64% 32.9M 2s
    ##  86750K .......... .......... .......... .......... .......... 64% 28.3M 2s
    ##  86800K .......... .......... .......... .......... .......... 64% 15.7M 2s
    ##  86850K .......... .......... .......... .......... .......... 64% 18.1M 2s
    ##  86900K .......... .......... .......... .......... .......... 64% 33.2M 2s
    ##  86950K .......... .......... .......... .......... .......... 64% 37.3M 2s
    ##  87000K .......... .......... .......... .......... .......... 64% 20.9M 2s
    ##  87050K .......... .......... .......... .......... .......... 64% 14.3M 2s
    ##  87100K .......... .......... .......... .......... .......... 65% 74.2M 2s
    ##  87150K .......... .......... .......... .......... .......... 65% 15.3M 2s
    ##  87200K .......... .......... .......... .......... .......... 65% 38.7M 2s
    ##  87250K .......... .......... .......... .......... .......... 65% 13.7M 2s
    ##  87300K .......... .......... .......... .......... .......... 65% 44.6M 2s
    ##  87350K .......... .......... .......... .......... .......... 65% 16.4M 2s
    ##  87400K .......... .......... .......... .......... .......... 65% 42.1M 2s
    ##  87450K .......... .......... .......... .......... .......... 65% 17.7M 2s
    ##  87500K .......... .......... .......... .......... .......... 65% 59.0M 2s
    ##  87550K .......... .......... .......... .......... .......... 65% 11.8M 2s
    ##  87600K .......... .......... .......... .......... .......... 65% 18.9M 2s
    ##  87650K .......... .......... .......... .......... .......... 65% 26.7M 2s
    ##  87700K .......... .......... .......... .......... .......... 65% 71.1M 2s
    ##  87750K .......... .......... .......... .......... .......... 65% 16.6M 2s
    ##  87800K .......... .......... .......... .......... .......... 65% 29.2M 2s
    ##  87850K .......... .......... .......... .......... .......... 65% 20.3M 2s
    ##  87900K .......... .......... .......... .......... .......... 65% 20.6M 2s
    ##  87950K .......... .......... .......... .......... .......... 65% 18.9M 2s
    ##  88000K .......... .......... .......... .......... .......... 65% 49.5M 2s
    ##  88050K .......... .......... .......... .......... .......... 65% 17.6M 2s
    ##  88100K .......... .......... .......... .......... .......... 65% 39.0M 2s
    ##  88150K .......... .......... .......... .......... .......... 65% 37.6M 2s
    ##  88200K .......... .......... .......... .......... .......... 65% 17.1M 2s
    ##  88250K .......... .......... .......... .......... .......... 65% 29.6M 2s
    ##  88300K .......... .......... .......... .......... .......... 65% 22.3M 2s
    ##  88350K .......... .......... .......... .......... .......... 65% 17.2M 2s
    ##  88400K .......... .......... .......... .......... .......... 65% 42.1M 2s
    ##  88450K .......... .......... .......... .......... .......... 66% 25.1M 2s
    ##  88500K .......... .......... .......... .......... .......... 66% 24.9M 2s
    ##  88550K .......... .......... .......... .......... .......... 66% 24.1M 2s
    ##  88600K .......... .......... .......... .......... .......... 66% 47.1M 2s
    ##  88650K .......... .......... .......... .......... .......... 66% 35.9M 2s
    ##  88700K .......... .......... .......... .......... .......... 66% 14.2M 2s
    ##  88750K .......... .......... .......... .......... .......... 66% 20.2M 2s
    ##  88800K .......... .......... .......... .......... .......... 66% 22.1M 2s
    ##  88850K .......... .......... .......... .......... .......... 66% 62.4M 2s
    ##  88900K .......... .......... .......... .......... .......... 66% 16.2M 2s
    ##  88950K .......... .......... .......... .......... .......... 66% 17.0M 2s
    ##  89000K .......... .......... .......... .......... .......... 66% 17.0M 2s
    ##  89050K .......... .......... .......... .......... .......... 66% 46.9M 2s
    ##  89100K .......... .......... .......... .......... .......... 66% 38.9M 2s
    ##  89150K .......... .......... .......... .......... .......... 66% 17.8M 2s
    ##  89200K .......... .......... .......... .......... .......... 66% 24.0M 2s
    ##  89250K .......... .......... .......... .......... .......... 66% 26.3M 2s
    ##  89300K .......... .......... .......... .......... .......... 66% 18.2M 2s
    ##  89350K .......... .......... .......... .......... .......... 66% 25.6M 2s
    ##  89400K .......... .......... .......... .......... .......... 66% 34.0M 2s
    ##  89450K .......... .......... .......... .......... .......... 66% 18.5M 2s
    ##  89500K .......... .......... .......... .......... .......... 66% 50.7M 2s
    ##  89550K .......... .......... .......... .......... .......... 66% 26.7M 2s
    ##  89600K .......... .......... .......... .......... .......... 66% 24.5M 2s
    ##  89650K .......... .......... .......... .......... .......... 66% 25.1M 2s
    ##  89700K .......... .......... .......... .......... .......... 66% 24.3M 2s
    ##  89750K .......... .......... .......... .......... .......... 66% 21.6M 2s
    ##  89800K .......... .......... .......... .......... .......... 67% 36.1M 2s
    ##  89850K .......... .......... .......... .......... .......... 67% 17.2M 2s
    ##  89900K .......... .......... .......... .......... .......... 67% 35.4M 2s
    ##  89950K .......... .......... .......... .......... .......... 67% 11.8M 2s
    ##  90000K .......... .......... .......... .......... .......... 67% 33.6M 2s
    ##  90050K .......... .......... .......... .......... .......... 67% 38.9M 2s
    ##  90100K .......... .......... .......... .......... .......... 67% 17.2M 2s
    ##  90150K .......... .......... .......... .......... .......... 67% 18.1M 2s
    ##  90200K .......... .......... .......... .......... .......... 67% 22.4M 2s
    ##  90250K .......... .......... .......... .......... .......... 67% 50.0M 2s
    ##  90300K .......... .......... .......... .......... .......... 67% 15.2M 2s
    ##  90350K .......... .......... .......... .......... .......... 67% 18.2M 2s
    ##  90400K .......... .......... .......... .......... .......... 67% 42.3M 2s
    ##  90450K .......... .......... .......... .......... .......... 67% 87.2M 2s
    ##  90500K .......... .......... .......... .......... .......... 67% 23.5M 2s
    ##  90550K .......... .......... .......... .......... .......... 67% 20.1M 2s
    ##  90600K .......... .......... .......... .......... .......... 67% 37.0M 2s
    ##  90650K .......... .......... .......... .......... .......... 67% 18.6M 2s
    ##  90700K .......... .......... .......... .......... .......... 67% 69.0M 2s
    ##  90750K .......... .......... .......... .......... .......... 67% 16.9M 2s
    ##  90800K .......... .......... .......... .......... .......... 67% 23.7M 2s
    ##  90850K .......... .......... .......... .......... .......... 67% 41.7M 2s
    ##  90900K .......... .......... .......... .......... .......... 67% 47.7M 2s
    ##  90950K .......... .......... .......... .......... .......... 67% 20.5M 2s
    ##  91000K .......... .......... .......... .......... .......... 67% 22.6M 2s
    ##  91050K .......... .......... .......... .......... .......... 67% 33.5M 2s
    ##  91100K .......... .......... .......... .......... .......... 67% 25.2M 2s
    ##  91150K .......... .......... .......... .......... .......... 68% 13.0M 2s
    ##  91200K .......... .......... .......... .......... .......... 68% 16.9M 2s
    ##  91250K .......... .......... .......... .......... .......... 68% 35.2M 2s
    ##  91300K .......... .......... .......... .......... .......... 68% 99.5M 2s
    ##  91350K .......... .......... .......... .......... .......... 68% 16.0M 2s
    ##  91400K .......... .......... .......... .......... .......... 68% 16.6M 2s
    ##  91450K .......... .......... .......... .......... .......... 68% 48.2M 2s
    ##  91500K .......... .......... .......... .......... .......... 68% 15.4M 2s
    ##  91550K .......... .......... .......... .......... .......... 68% 17.0M 2s
    ##  91600K .......... .......... .......... .......... .......... 68% 28.3M 2s
    ##  91650K .......... .......... .......... .......... .......... 68% 23.6M 2s
    ##  91700K .......... .......... .......... .......... .......... 68% 20.5M 2s
    ##  91750K .......... .......... .......... .......... .......... 68% 19.0M 2s
    ##  91800K .......... .......... .......... .......... .......... 68% 26.1M 2s
    ##  91850K .......... .......... .......... .......... .......... 68% 95.0M 2s
    ##  91900K .......... .......... .......... .......... .......... 68% 13.3M 2s
    ##  91950K .......... .......... .......... .......... .......... 68% 15.0M 2s
    ##  92000K .......... .......... .......... .......... .......... 68% 27.3M 2s
    ##  92050K .......... .......... .......... .......... .......... 68% 95.8M 2s
    ##  92100K .......... .......... .......... .......... .......... 68% 15.9M 2s
    ##  92150K .......... .......... .......... .......... .......... 68% 36.4M 2s
    ##  92200K .......... .......... .......... .......... .......... 68% 17.7M 2s
    ##  92250K .......... .......... .......... .......... .......... 68% 37.4M 2s
    ##  92300K .......... .......... .......... .......... .......... 68% 54.2M 2s
    ##  92350K .......... .......... .......... .......... .......... 68% 11.0M 2s
    ##  92400K .......... .......... .......... .......... .......... 68% 74.2M 2s
    ##  92450K .......... .......... .......... .......... .......... 68% 13.9M 2s
    ##  92500K .......... .......... .......... .......... .......... 69% 48.9M 2s
    ##  92550K .......... .......... .......... .......... .......... 69% 16.1M 2s
    ##  92600K .......... .......... .......... .......... .......... 69% 32.7M 2s
    ##  92650K .......... .......... .......... .......... .......... 69% 17.3M 2s
    ##  92700K .......... .......... .......... .......... .......... 69% 39.0M 2s
    ##  92750K .......... .......... .......... .......... .......... 69% 20.5M 2s
    ##  92800K .......... .......... .......... .......... .......... 69% 19.3M 2s
    ##  92850K .......... .......... .......... .......... .......... 69% 35.7M 2s
    ##  92900K .......... .......... .......... .......... .......... 69% 15.8M 2s
    ##  92950K .......... .......... .......... .......... .......... 69% 19.7M 2s
    ##  93000K .......... .......... .......... .......... .......... 69% 30.6M 2s
    ##  93050K .......... .......... .......... .......... .......... 69% 72.8M 2s
    ##  93100K .......... .......... .......... .......... .......... 69% 13.0M 2s
    ##  93150K .......... .......... .......... .......... .......... 69% 16.4M 2s
    ##  93200K .......... .......... .......... .......... .......... 69% 30.4M 2s
    ##  93250K .......... .......... .......... .......... .......... 69% 89.7M 2s
    ##  93300K .......... .......... .......... .......... .......... 69% 21.9M 2s
    ##  93350K .......... .......... .......... .......... .......... 69% 23.9M 2s
    ##  93400K .......... .......... .......... .......... .......... 69% 27.2M 2s
    ##  93450K .......... .......... .......... .......... .......... 69% 22.4M 2s
    ##  93500K .......... .......... .......... .......... .......... 69% 49.5M 2s
    ##  93550K .......... .......... .......... .......... .......... 69% 14.3M 2s
    ##  93600K .......... .......... .......... .......... .......... 69% 35.1M 2s
    ##  93650K .......... .......... .......... .......... .......... 69% 14.5M 2s
    ##  93700K .......... .......... .......... .......... .......... 69% 18.9M 2s
    ##  93750K .......... .......... .......... .......... .......... 69% 24.5M 2s
    ##  93800K .......... .......... .......... .......... .......... 70% 93.5M 2s
    ##  93850K .......... .......... .......... .......... .......... 70% 17.6M 2s
    ##  93900K .......... .......... .......... .......... .......... 70% 48.5M 2s
    ##  93950K .......... .......... .......... .......... .......... 70% 13.5M 2s
    ##  94000K .......... .......... .......... .......... .......... 70% 54.7M 2s
    ##  94050K .......... .......... .......... .......... .......... 70% 19.5M 2s
    ##  94100K .......... .......... .......... .......... .......... 70% 18.7M 2s
    ##  94150K .......... .......... .......... .......... .......... 70% 39.7M 2s
    ##  94200K .......... .......... .......... .......... .......... 70% 15.9M 2s
    ##  94250K .......... .......... .......... .......... .......... 70% 83.4M 2s
    ##  94300K .......... .......... .......... .......... .......... 70% 14.1M 2s
    ##  94350K .......... .......... .......... .......... .......... 70% 13.7M 2s
    ##  94400K .......... .......... .......... .......... .......... 70% 54.2M 2s
    ##  94450K .......... .......... .......... .......... .......... 70% 16.7M 2s
    ##  94500K .......... .......... .......... .......... .......... 70% 20.3M 2s
    ##  94550K .......... .......... .......... .......... .......... 70% 32.9M 2s
    ##  94600K .......... .......... .......... .......... .......... 70% 19.7M 2s
    ##  94650K .......... .......... .......... .......... .......... 70% 30.3M 2s
    ##  94700K .......... .......... .......... .......... .......... 70% 97.4M 2s
    ##  94750K .......... .......... .......... .......... .......... 70% 13.0M 2s
    ##  94800K .......... .......... .......... .......... .......... 70% 16.5M 2s
    ##  94850K .......... .......... .......... .......... .......... 70% 25.2M 2s
    ##  94900K .......... .......... .......... .......... .......... 70% 78.2M 2s
    ##  94950K .......... .......... .......... .......... .......... 70% 22.2M 2s
    ##  95000K .......... .......... .......... .......... .......... 70% 14.8M 2s
    ##  95050K .......... .......... .......... .......... .......... 70% 16.2M 2s
    ##  95100K .......... .......... .......... .......... .......... 70% 29.2M 2s
    ##  95150K .......... .......... .......... .......... .......... 71% 50.9M 2s
    ##  95200K .......... .......... .......... .......... .......... 71% 22.0M 2s
    ##  95250K .......... .......... .......... .......... .......... 71% 23.6M 2s
    ##  95300K .......... .......... .......... .......... .......... 71% 25.7M 2s
    ##  95350K .......... .......... .......... .......... .......... 71% 19.2M 2s
    ##  95400K .......... .......... .......... .......... .......... 71% 40.7M 2s
    ##  95450K .......... .......... .......... .......... .......... 71% 18.2M 2s
    ##  95500K .......... .......... .......... .......... .......... 71% 42.3M 2s
    ##  95550K .......... .......... .......... .......... .......... 71% 12.5M 2s
    ##  95600K .......... .......... .......... .......... .......... 71% 27.4M 2s
    ##  95650K .......... .......... .......... .......... .......... 71% 38.8M 2s
    ##  95700K .......... .......... .......... .......... .......... 71% 16.2M 2s
    ##  95750K .......... .......... .......... .......... .......... 71% 18.5M 2s
    ##  95800K .......... .......... .......... .......... .......... 71% 29.2M 2s
    ##  95850K .......... .......... .......... .......... .......... 71% 35.8M 2s
    ##  95900K .......... .......... .......... .......... .......... 71% 27.5M 2s
    ##  95950K .......... .......... .......... .......... .......... 71% 18.6M 2s
    ##  96000K .......... .......... .......... .......... .......... 71% 11.4M 2s
    ##  96050K .......... .......... .......... .......... .......... 71% 83.1M 2s
    ##  96100K .......... .......... .......... .......... .......... 71% 30.2M 2s
    ##  96150K .......... .......... .......... .......... .......... 71% 14.9M 2s
    ##  96200K .......... .......... .......... .......... .......... 71% 27.1M 2s
    ##  96250K .......... .......... .......... .......... .......... 71% 15.7M 2s
    ##  96300K .......... .......... .......... .......... .......... 71% 45.2M 2s
    ##  96350K .......... .......... .......... .......... .......... 71% 19.7M 2s
    ##  96400K .......... .......... .......... .......... .......... 71% 36.1M 2s
    ##  96450K .......... .......... .......... .......... .......... 71% 17.8M 2s
    ##  96500K .......... .......... .......... .......... .......... 72% 38.5M 2s
    ##  96550K .......... .......... .......... .......... .......... 72% 20.1M 2s
    ##  96600K .......... .......... .......... .......... .......... 72% 44.5M 2s
    ##  96650K .......... .......... .......... .......... .......... 72% 11.6M 2s
    ##  96700K .......... .......... .......... .......... .......... 72% 78.5M 2s
    ##  96750K .......... .......... .......... .......... .......... 72% 14.3M 2s
    ##  96800K .......... .......... .......... .......... .......... 72% 44.0M 2s
    ##  96850K .......... .......... .......... .......... .......... 72% 18.5M 2s
    ##  96900K .......... .......... .......... .......... .......... 72% 16.1M 2s
    ##  96950K .......... .......... .......... .......... .......... 72% 19.4M 2s
    ##  97000K .......... .......... .......... .......... .......... 72% 36.1M 2s
    ##  97050K .......... .......... .......... .......... .......... 72% 46.7M 2s
    ##  97100K .......... .......... .......... .......... .......... 72% 17.9M 2s
    ##  97150K .......... .......... .......... .......... .......... 72% 28.0M 2s
    ##  97200K .......... .......... .......... .......... .......... 72% 22.9M 2s
    ##  97250K .......... .......... .......... .......... .......... 72% 21.5M 2s
    ##  97300K .......... .......... .......... .......... .......... 72% 43.3M 2s
    ##  97350K .......... .......... .......... .......... .......... 72% 18.2M 2s
    ##  97400K .......... .......... .......... .......... .......... 72% 43.0M 2s
    ##  97450K .......... .......... .......... .......... .......... 72% 18.2M 2s
    ##  97500K .......... .......... .......... .......... .......... 72% 43.1M 2s
    ##  97550K .......... .......... .......... .......... .......... 72% 19.0M 2s
    ##  97600K .......... .......... .......... .......... .......... 72% 17.3M 2s
    ##  97650K .......... .......... .......... .......... .......... 72% 15.0M 2s
    ##  97700K .......... .......... .......... .......... .......... 72% 39.9M 2s
    ##  97750K .......... .......... .......... .......... .......... 72% 25.3M 2s
    ##  97800K .......... .......... .......... .......... .......... 72% 26.8M 2s
    ##  97850K .......... .......... .......... .......... .......... 73% 14.2M 2s
    ##  97900K .......... .......... .......... .......... .......... 73% 14.8M 2s
    ##  97950K .......... .......... .......... .......... .......... 73% 44.3M 2s
    ##  98000K .......... .......... .......... .......... .......... 73% 25.7M 2s
    ##  98050K .......... .......... .......... .......... .......... 73% 24.1M 2s
    ##  98100K .......... .......... .......... .......... .......... 73% 22.1M 2s
    ##  98150K .......... .......... .......... .......... .......... 73% 25.6M 2s
    ##  98200K .......... .......... .......... .......... .......... 73% 20.6M 2s
    ##  98250K .......... .......... .......... .......... .......... 73% 38.8M 2s
    ##  98300K .......... .......... .......... .......... .......... 73% 19.9M 2s
    ##  98350K .......... .......... .......... .......... .......... 73% 23.4M 2s
    ##  98400K .......... .......... .......... .......... .......... 73% 67.8M 2s
    ##  98450K .......... .......... .......... .......... .......... 73% 20.9M 2s
    ##  98500K .......... .......... .......... .......... .......... 73% 36.3M 2s
    ##  98550K .......... .......... .......... .......... .......... 73% 14.1M 2s
    ##  98600K .......... .......... .......... .......... .......... 73% 30.0M 2s
    ##  98650K .......... .......... .......... .......... .......... 73% 22.7M 2s
    ##  98700K .......... .......... .......... .......... .......... 73% 41.1M 2s
    ##  98750K .......... .......... .......... .......... .......... 73% 17.0M 2s
    ##  98800K .......... .......... .......... .......... .......... 73% 36.8M 2s
    ##  98850K .......... .......... .......... .......... .......... 73% 31.1M 2s
    ##  98900K .......... .......... .......... .......... .......... 73% 34.6M 2s
    ##  98950K .......... .......... .......... .......... .......... 73% 16.9M 2s
    ##  99000K .......... .......... .......... .......... .......... 73% 16.1M 2s
    ##  99050K .......... .......... .......... .......... .......... 73% 46.2M 2s
    ##  99100K .......... .......... .......... .......... .......... 73% 23.5M 2s
    ##  99150K .......... .......... .......... .......... .......... 73% 16.3M 2s
    ##  99200K .......... .......... .......... .......... .......... 74% 45.6M 2s
    ##  99250K .......... .......... .......... .......... .......... 74% 12.0M 2s
    ##  99300K .......... .......... .......... .......... .......... 74% 51.0M 2s
    ##  99350K .......... .......... .......... .......... .......... 74% 22.1M 2s
    ##  99400K .......... .......... .......... .......... .......... 74% 26.6M 2s
    ##  99450K .......... .......... .......... .......... .......... 74% 16.1M 2s
    ##  99500K .......... .......... .......... .......... .......... 74% 36.1M 2s
    ##  99550K .......... .......... .......... .......... .......... 74% 33.9M 2s
    ##  99600K .......... .......... .......... .......... .......... 74% 28.6M 2s
    ##  99650K .......... .......... .......... .......... .......... 74% 21.7M 2s
    ##  99700K .......... .......... .......... .......... .......... 74% 26.4M 2s
    ##  99750K .......... .......... .......... .......... .......... 74% 17.3M 2s
    ##  99800K .......... .......... .......... .......... .......... 74% 46.2M 2s
    ##  99850K .......... .......... .......... .......... .......... 74% 19.4M 2s
    ##  99900K .......... .......... .......... .......... .......... 74% 35.4M 2s
    ##  99950K .......... .......... .......... .......... .......... 74% 13.5M 2s
    ## 100000K .......... .......... .......... .......... .......... 74% 38.7M 2s
    ## 100050K .......... .......... .......... .......... .......... 74% 18.4M 2s
    ## 100100K .......... .......... .......... .......... .......... 74% 19.0M 2s
    ## 100150K .......... .......... .......... .......... .......... 74% 26.4M 2s
    ## 100200K .......... .......... .......... .......... .......... 74% 10.6M 2s
    ## 100250K .......... .......... .......... .......... .......... 74% 62.7M 2s
    ## 100300K .......... .......... .......... .......... .......... 74% 62.2M 2s
    ## 100350K .......... .......... .......... .......... .......... 74% 13.6M 2s
    ## 100400K .......... .......... .......... .......... .......... 74% 17.0M 2s
    ## 100450K .......... .......... .......... .......... .......... 74% 52.2M 2s
    ## 100500K .......... .......... .......... .......... .......... 75% 13.8M 2s
    ## 100550K .......... .......... .......... .......... .......... 75% 18.8M 2s
    ## 100600K .......... .......... .......... .......... .......... 75% 56.2M 2s
    ## 100650K .......... .......... .......... .......... .......... 75% 15.0M 2s
    ## 100700K .......... .......... .......... .......... .......... 75% 61.6M 2s
    ## 100750K .......... .......... .......... .......... .......... 75% 15.6M 2s
    ## 100800K .......... .......... .......... .......... .......... 75% 27.8M 2s
    ## 100850K .......... .......... .......... .......... .......... 75% 22.2M 2s
    ## 100900K .......... .......... .......... .......... .......... 75% 30.9M 2s
    ## 100950K .......... .......... .......... .......... .......... 75% 18.2M 2s
    ## 101000K .......... .......... .......... .......... .......... 75% 43.6M 2s
    ## 101050K .......... .......... .......... .......... .......... 75% 16.7M 2s
    ## 101100K .......... .......... .......... .......... .......... 75% 31.9M 2s
    ## 101150K .......... .......... .......... .......... .......... 75% 18.1M 2s
    ## 101200K .......... .......... .......... .......... .......... 75% 16.9M 2s
    ## 101250K .......... .......... .......... .......... .......... 75% 51.8M 2s
    ## 101300K .......... .......... .......... .......... .......... 75% 12.5M 2s
    ## 101350K .......... .......... .......... .......... .......... 75% 13.4M 2s
    ## 101400K .......... .......... .......... .......... .......... 75% 35.8M 2s
    ## 101450K .......... .......... .......... .......... .......... 75% 28.8M 2s
    ## 101500K .......... .......... .......... .......... .......... 75% 48.1M 2s
    ## 101550K .......... .......... .......... .......... .......... 75% 13.6M 1s
    ## 101600K .......... .......... .......... .......... .......... 75% 18.2M 1s
    ## 101650K .......... .......... .......... .......... .......... 75% 35.5M 1s
    ## 101700K .......... .......... .......... .......... .......... 75% 46.0M 1s
    ## 101750K .......... .......... .......... .......... .......... 75% 12.1M 1s
    ## 101800K .......... .......... .......... .......... .......... 75% 81.7M 1s
    ## 101850K .......... .......... .......... .......... .......... 76% 13.6M 1s
    ## 101900K .......... .......... .......... .......... .......... 76% 17.1M 1s
    ## 101950K .......... .......... .......... .......... .......... 76% 31.3M 1s
    ## 102000K .......... .......... .......... .......... .......... 76% 25.1M 1s
    ## 102050K .......... .......... .......... .......... .......... 76% 17.4M 1s
    ## 102100K .......... .......... .......... .......... .......... 76% 40.5M 1s
    ## 102150K .......... .......... .......... .......... .......... 76% 12.3M 1s
    ## 102200K .......... .......... .......... .......... .......... 76% 89.0M 1s
    ## 102250K .......... .......... .......... .......... .......... 76% 23.2M 1s
    ## 102300K .......... .......... .......... .......... .......... 76% 15.6M 1s
    ## 102350K .......... .......... .......... .......... .......... 76% 39.8M 1s
    ## 102400K .......... .......... .......... .......... .......... 76% 14.6M 1s
    ## 102450K .......... .......... .......... .......... .......... 76% 47.0M 1s
    ## 102500K .......... .......... .......... .......... .......... 76% 16.8M 1s
    ## 102550K .......... .......... .......... .......... .......... 76% 40.1M 1s
    ## 102600K .......... .......... .......... .......... .......... 76% 17.5M 1s
    ## 102650K .......... .......... .......... .......... .......... 76% 24.5M 1s
    ## 102700K .......... .......... .......... .......... .......... 76% 27.5M 1s
    ## 102750K .......... .......... .......... .......... .......... 76% 16.3M 1s
    ## 102800K .......... .......... .......... .......... .......... 76% 33.7M 1s
    ## 102850K .......... .......... .......... .......... .......... 76% 12.1M 1s
    ## 102900K .......... .......... .......... .......... .......... 76% 92.1M 1s
    ## 102950K .......... .......... .......... .......... .......... 76% 16.7M 1s
    ## 103000K .......... .......... .......... .......... .......... 76% 16.0M 1s
    ## 103050K .......... .......... .......... .......... .......... 76% 50.4M 1s
    ## 103100K .......... .......... .......... .......... .......... 76% 15.6M 1s
    ## 103150K .......... .......... .......... .......... .......... 76% 35.5M 1s
    ## 103200K .......... .......... .......... .......... .......... 77% 20.8M 1s
    ## 103250K .......... .......... .......... .......... .......... 77% 29.8M 1s
    ## 103300K .......... .......... .......... .......... .......... 77% 24.7M 1s
    ## 103350K .......... .......... .......... .......... .......... 77% 21.8M 1s
    ## 103400K .......... .......... .......... .......... .......... 77% 19.2M 1s
    ## 103450K .......... .......... .......... .......... .......... 77% 25.3M 1s
    ## 103500K .......... .......... .......... .......... .......... 77% 29.6M 1s
    ## 103550K .......... .......... .......... .......... .......... 77% 17.5M 1s
    ## 103600K .......... .......... .......... .......... .......... 77% 51.2M 1s
    ## 103650K .......... .......... .......... .......... .......... 77% 10.6M 1s
    ## 103700K .......... .......... .......... .......... .......... 77% 46.0M 1s
    ## 103750K .......... .......... .......... .......... .......... 77% 19.2M 1s
    ## 103800K .......... .......... .......... .......... .......... 77% 40.9M 1s
    ## 103850K .......... .......... .......... .......... .......... 77% 19.8M 1s
    ## 103900K .......... .......... .......... .......... .......... 77% 27.2M 1s
    ## 103950K .......... .......... .......... .......... .......... 77% 24.2M 1s
    ## 104000K .......... .......... .......... .......... .......... 77% 24.1M 1s
    ## 104050K .......... .......... .......... .......... .......... 77% 20.3M 1s
    ## 104100K .......... .......... .......... .......... .......... 77% 19.5M 1s
    ## 104150K .......... .......... .......... .......... .......... 77% 41.2M 1s
    ## 104200K .......... .......... .......... .......... .......... 77% 35.9M 1s
    ## 104250K .......... .......... .......... .......... .......... 77% 20.7M 1s
    ## 104300K .......... .......... .......... .......... .......... 77% 23.3M 1s
    ## 104350K .......... .......... .......... .......... .......... 77% 25.8M 1s
    ## 104400K .......... .......... .......... .......... .......... 77% 16.6M 1s
    ## 104450K .......... .......... .......... .......... .......... 77% 73.7M 1s
    ## 104500K .......... .......... .......... .......... .......... 77% 18.7M 1s
    ## 104550K .......... .......... .......... .......... .......... 78% 19.8M 1s
    ## 104600K .......... .......... .......... .......... .......... 78% 36.7M 1s
    ## 104650K .......... .......... .......... .......... .......... 78% 35.9M 1s
    ## 104700K .......... .......... .......... .......... .......... 78% 27.7M 1s
    ## 104750K .......... .......... .......... .......... .......... 78% 12.8M 1s
    ## 104800K .......... .......... .......... .......... .......... 78% 17.6M 1s
    ## 104850K .......... .......... .......... .......... .......... 78% 41.9M 1s
    ## 104900K .......... .......... .......... .......... .......... 78% 11.5M 1s
    ## 104950K .......... .......... .......... .......... .......... 78% 32.9M 1s
    ## 105000K .......... .......... .......... .......... .......... 78% 22.5M 1s
    ## 105050K .......... .......... .......... .......... .......... 78% 18.0M 1s
    ## 105100K .......... .......... .......... .......... .......... 78% 45.9M 1s
    ## 105150K .......... .......... .......... .......... .......... 78% 10.6M 1s
    ## 105200K .......... .......... .......... .......... .......... 78% 52.1M 1s
    ## 105250K .......... .......... .......... .......... .......... 78% 21.3M 1s
    ## 105300K .......... .......... .......... .......... .......... 78% 17.0M 1s
    ## 105350K .......... .......... .......... .......... .......... 78% 44.8M 1s
    ## 105400K .......... .......... .......... .......... .......... 78% 15.4M 1s
    ## 105450K .......... .......... .......... .......... .......... 78% 18.8M 1s
    ## 105500K .......... .......... .......... .......... .......... 78% 46.5M 1s
    ## 105550K .......... .......... .......... .......... .......... 78% 16.4M 1s
    ## 105600K .......... .......... .......... .......... .......... 78% 44.7M 1s
    ## 105650K .......... .......... .......... .......... .......... 78% 16.2M 1s
    ## 105700K .......... .......... .......... .......... .......... 78% 32.1M 1s
    ## 105750K .......... .......... .......... .......... .......... 78% 45.1M 1s
    ## 105800K .......... .......... .......... .......... .......... 78% 22.7M 1s
    ## 105850K .......... .......... .......... .......... .......... 78% 14.4M 1s
    ## 105900K .......... .......... .......... .......... .......... 79% 72.5M 1s
    ## 105950K .......... .......... .......... .......... .......... 79% 10.9M 1s
    ## 106000K .......... .......... .......... .......... .......... 79% 78.1M 1s
    ## 106050K .......... .......... .......... .......... .......... 79% 20.9M 1s
    ## 106100K .......... .......... .......... .......... .......... 79% 39.7M 1s
    ## 106150K .......... .......... .......... .......... .......... 79% 15.2M 1s
    ## 106200K .......... .......... .......... .......... .......... 79% 46.6M 1s
    ## 106250K .......... .......... .......... .......... .......... 79% 22.5M 1s
    ## 106300K .......... .......... .......... .......... .......... 79% 58.4M 1s
    ## 106350K .......... .......... .......... .......... .......... 79% 13.2M 1s
    ## 106400K .......... .......... .......... .......... .......... 79% 16.5M 1s
    ## 106450K .......... .......... .......... .......... .......... 79% 45.0M 1s
    ## 106500K .......... .......... .......... .......... .......... 79% 17.8M 1s
    ## 106550K .......... .......... .......... .......... .......... 79% 18.0M 1s
    ## 106600K .......... .......... .......... .......... .......... 79% 29.6M 1s
    ## 106650K .......... .......... .......... .......... .......... 79% 21.3M 1s
    ## 106700K .......... .......... .......... .......... .......... 79% 35.5M 1s
    ## 106750K .......... .......... .......... .......... .......... 79% 35.3M 1s
    ## 106800K .......... .......... .......... .......... .......... 79% 21.5M 1s
    ## 106850K .......... .......... .......... .......... .......... 79% 29.1M 1s
    ## 106900K .......... .......... .......... .......... .......... 79% 14.5M 1s
    ## 106950K .......... .......... .......... .......... .......... 79% 35.3M 1s
    ## 107000K .......... .......... .......... .......... .......... 79% 19.9M 1s
    ## 107050K .......... .......... .......... .......... .......... 79% 16.7M 1s
    ## 107100K .......... .......... .......... .......... .......... 79% 31.0M 1s
    ## 107150K .......... .......... .......... .......... .......... 79% 14.7M 1s
    ## 107200K .......... .......... .......... .......... .......... 79% 35.0M 1s
    ## 107250K .......... .......... .......... .......... .......... 80% 28.0M 1s
    ## 107300K .......... .......... .......... .......... .......... 80% 23.6M 1s
    ## 107350K .......... .......... .......... .......... .......... 80% 15.6M 1s
    ## 107400K .......... .......... .......... .......... .......... 80% 34.7M 1s
    ## 107450K .......... .......... .......... .......... .......... 80% 16.8M 1s
    ## 107500K .......... .......... .......... .......... .......... 80% 51.9M 1s
    ## 107550K .......... .......... .......... .......... .......... 80% 15.0M 1s
    ## 107600K .......... .......... .......... .......... .......... 80% 43.7M 1s
    ## 107650K .......... .......... .......... .......... .......... 80% 16.6M 1s
    ## 107700K .......... .......... .......... .......... .......... 80% 50.4M 1s
    ## 107750K .......... .......... .......... .......... .......... 80% 11.1M 1s
    ## 107800K .......... .......... .......... .......... .......... 80% 36.0M 1s
    ## 107850K .......... .......... .......... .......... .......... 80% 12.9M 1s
    ## 107900K .......... .......... .......... .......... .......... 80% 86.8M 1s
    ## 107950K .......... .......... .......... .......... .......... 80% 10.2M 1s
    ## 108000K .......... .......... .......... .......... .......... 80% 99.1M 1s
    ## 108050K .......... .......... .......... .......... .......... 80% 10.9M 1s
    ## 108100K .......... .......... .......... .......... .......... 80% 90.1M 1s
    ## 108150K .......... .......... .......... .......... .......... 80% 11.6M 1s
    ## 108200K .......... .......... .......... .......... .......... 80% 98.3M 1s
    ## 108250K .......... .......... .......... .......... .......... 80% 17.3M 1s
    ## 108300K .......... .......... .......... .......... .......... 80% 33.4M 1s
    ## 108350K .......... .......... .......... .......... .......... 80% 16.1M 1s
    ## 108400K .......... .......... .......... .......... .......... 80% 17.3M 1s
    ## 108450K .......... .......... .......... .......... .......... 80% 34.3M 1s
    ## 108500K .......... .......... .......... .......... .......... 80% 20.1M 1s
    ## 108550K .......... .......... .......... .......... .......... 81% 32.6M 1s
    ## 108600K .......... .......... .......... .......... .......... 81% 8.33M 1s
    ## 108650K .......... .......... .......... .......... .......... 81% 81.8M 1s
    ## 108700K .......... .......... .......... .......... .......... 81% 16.8M 1s
    ## 108750K .......... .......... .......... .......... .......... 81% 17.0M 1s
    ## 108800K .......... .......... .......... .......... .......... 81% 15.4M 1s
    ## 108850K .......... .......... .......... .......... .......... 81% 17.2M 1s
    ## 108900K .......... .......... .......... .......... .......... 81% 32.6M 1s
    ## 108950K .......... .......... .......... .......... .......... 81% 15.8M 1s
    ## 109000K .......... .......... .......... .......... .......... 81% 46.5M 1s
    ## 109050K .......... .......... .......... .......... .......... 81% 17.1M 1s
    ## 109100K .......... .......... .......... .......... .......... 81% 49.6M 1s
    ## 109150K .......... .......... .......... .......... .......... 81% 15.1M 1s
    ## 109200K .......... .......... .......... .......... .......... 81% 47.0M 1s
    ## 109250K .......... .......... .......... .......... .......... 81% 16.1M 1s
    ## 109300K .......... .......... .......... .......... .......... 81% 44.8M 1s
    ## 109350K .......... .......... .......... .......... .......... 81% 15.5M 1s
    ## 109400K .......... .......... .......... .......... .......... 81% 55.8M 1s
    ## 109450K .......... .......... .......... .......... .......... 81% 16.1M 1s
    ## 109500K .......... .......... .......... .......... .......... 81% 47.5M 1s
    ## 109550K .......... .......... .......... .......... .......... 81% 15.5M 1s
    ## 109600K .......... .......... .......... .......... .......... 81% 41.2M 1s
    ## 109650K .......... .......... .......... .......... .......... 81% 17.0M 1s
    ## 109700K .......... .......... .......... .......... .......... 81% 47.9M 1s
    ## 109750K .......... .......... .......... .......... .......... 81% 16.7M 1s
    ## 109800K .......... .......... .......... .......... .......... 81% 45.5M 1s
    ## 109850K .......... .......... .......... .......... .......... 81% 15.0M 1s
    ## 109900K .......... .......... .......... .......... .......... 82% 75.6M 1s
    ## 109950K .......... .......... .......... .......... .......... 82% 14.5M 1s
    ## 110000K .......... .......... .......... .......... .......... 82% 16.4M 1s
    ## 110050K .......... .......... .......... .......... .......... 82% 46.3M 1s
    ## 110100K .......... .......... .......... .......... .......... 82% 19.5M 1s
    ## 110150K .......... .......... .......... .......... .......... 82% 32.2M 1s
    ## 110200K .......... .......... .......... .......... .......... 82% 17.6M 1s
    ## 110250K .......... .......... .......... .......... .......... 82% 40.0M 1s
    ## 110300K .......... .......... .......... .......... .......... 82% 17.0M 1s
    ## 110350K .......... .......... .......... .......... .......... 82% 39.6M 1s
    ## 110400K .......... .......... .......... .......... .......... 82% 17.4M 1s
    ## 110450K .......... .......... .......... .......... .......... 82% 47.1M 1s
    ## 110500K .......... .......... .......... .......... .......... 82% 17.9M 1s
    ## 110550K .......... .......... .......... .......... .......... 82% 43.1M 1s
    ## 110600K .......... .......... .......... .......... .......... 82% 17.7M 1s
    ## 110650K .......... .......... .......... .......... .......... 82% 47.5M 1s
    ## 110700K .......... .......... .......... .......... .......... 82% 11.7M 1s
    ## 110750K .......... .......... .......... .......... .......... 82% 44.0M 1s
    ## 110800K .......... .......... .......... .......... .......... 82% 17.0M 1s
    ## 110850K .......... .......... .......... .......... .......... 82% 47.2M 1s
    ## 110900K .......... .......... .......... .......... .......... 82% 11.3M 1s
    ## 110950K .......... .......... .......... .......... .......... 82% 80.2M 1s
    ## 111000K .......... .......... .......... .......... .......... 82% 15.0M 1s
    ## 111050K .......... .......... .......... .......... .......... 82% 57.0M 1s
    ## 111100K .......... .......... .......... .......... .......... 82% 10.5M 1s
    ## 111150K .......... .......... .......... .......... .......... 82% 16.1M 1s
    ## 111200K .......... .......... .......... .......... .......... 82% 30.1M 1s
    ## 111250K .......... .......... .......... .......... .......... 83% 18.6M 1s
    ## 111300K .......... .......... .......... .......... .......... 83% 40.5M 1s
    ## 111350K .......... .......... .......... .......... .......... 83% 14.5M 1s
    ## 111400K .......... .......... .......... .......... .......... 83% 48.0M 1s
    ## 111450K .......... .......... .......... .......... .......... 83% 16.6M 1s
    ## 111500K .......... .......... .......... .......... .......... 83% 45.9M 1s
    ## 111550K .......... .......... .......... .......... .......... 83% 15.9M 1s
    ## 111600K .......... .......... .......... .......... .......... 83% 39.1M 1s
    ## 111650K .......... .......... .......... .......... .......... 83% 16.8M 1s
    ## 111700K .......... .......... .......... .......... .......... 83% 43.5M 1s
    ## 111750K .......... .......... .......... .......... .......... 83% 16.6M 1s
    ## 111800K .......... .......... .......... .......... .......... 83% 54.3M 1s
    ## 111850K .......... .......... .......... .......... .......... 83% 17.4M 1s
    ## 111900K .......... .......... .......... .......... .......... 83% 44.0M 1s
    ## 111950K .......... .......... .......... .......... .......... 83% 13.9M 1s
    ## 112000K .......... .......... .......... .......... .......... 83% 15.8M 1s
    ## 112050K .......... .......... .......... .......... .......... 83% 48.3M 1s
    ## 112100K .......... .......... .......... .......... .......... 83% 17.2M 1s
    ## 112150K .......... .......... .......... .......... .......... 83% 43.8M 1s
    ## 112200K .......... .......... .......... .......... .......... 83% 17.5M 1s
    ## 112250K .......... .......... .......... .......... .......... 83% 46.9M 1s
    ## 112300K .......... .......... .......... .......... .......... 83% 14.9M 1s
    ## 112350K .......... .......... .......... .......... .......... 83% 49.1M 1s
    ## 112400K .......... .......... .......... .......... .......... 83% 16.4M 1s
    ## 112450K .......... .......... .......... .......... .......... 83% 51.1M 1s
    ## 112500K .......... .......... .......... .......... .......... 83% 15.1M 1s
    ## 112550K .......... .......... .......... .......... .......... 83% 55.5M 1s
    ## 112600K .......... .......... .......... .......... .......... 84% 15.1M 1s
    ## 112650K .......... .......... .......... .......... .......... 84% 52.9M 1s
    ## 112700K .......... .......... .......... .......... .......... 84% 17.2M 1s
    ## 112750K .......... .......... .......... .......... .......... 84% 14.8M 1s
    ## 112800K .......... .......... .......... .......... .......... 84% 39.4M 1s
    ## 112850K .......... .......... .......... .......... .......... 84% 17.6M 1s
    ## 112900K .......... .......... .......... .......... .......... 84% 41.8M 1s
    ## 112950K .......... .......... .......... .......... .......... 84% 19.7M 1s
    ## 113000K .......... .......... .......... .......... .......... 84% 29.9M 1s
    ## 113050K .......... .......... .......... .......... .......... 84% 14.4M 1s
    ## 113100K .......... .......... .......... .......... .......... 84% 64.1M 1s
    ## 113150K .......... .......... .......... .......... .......... 84% 14.4M 1s
    ## 113200K .......... .......... .......... .......... .......... 84% 15.4M 1s
    ## 113250K .......... .......... .......... .......... .......... 84% 42.6M 1s
    ## 113300K .......... .......... .......... .......... .......... 84% 17.0M 1s
    ## 113350K .......... .......... .......... .......... .......... 84% 39.5M 1s
    ## 113400K .......... .......... .......... .......... .......... 84% 16.3M 1s
    ## 113450K .......... .......... .......... .......... .......... 84% 50.0M 1s
    ## 113500K .......... .......... .......... .......... .......... 84% 15.9M 1s
    ## 113550K .......... .......... .......... .......... .......... 84% 46.8M 1s
    ## 113600K .......... .......... .......... .......... .......... 84% 14.8M 1s
    ## 113650K .......... .......... .......... .......... .......... 84% 16.9M 1s
    ## 113700K .......... .......... .......... .......... .......... 84% 43.2M 1s
    ## 113750K .......... .......... .......... .......... .......... 84% 12.3M 1s
    ## 113800K .......... .......... .......... .......... .......... 84% 85.3M 1s
    ## 113850K .......... .......... .......... .......... .......... 84% 13.4M 1s
    ## 113900K .......... .......... .......... .......... .......... 84% 92.3M 1s
    ## 113950K .......... .......... .......... .......... .......... 85% 13.7M 1s
    ## 114000K .......... .......... .......... .......... .......... 85% 20.7M 1s
    ## 114050K .......... .......... .......... .......... .......... 85% 30.8M 1s
    ## 114100K .......... .......... .......... .......... .......... 85% 17.4M 1s
    ## 114150K .......... .......... .......... .......... .......... 85% 35.2M 1s
    ## 114200K .......... .......... .......... .......... .......... 85% 17.3M 1s
    ## 114250K .......... .......... .......... .......... .......... 85% 45.1M 1s
    ## 114300K .......... .......... .......... .......... .......... 85% 16.7M 1s
    ## 114350K .......... .......... .......... .......... .......... 85% 38.4M 1s
    ## 114400K .......... .......... .......... .......... .......... 85% 16.1M 1s
    ## 114450K .......... .......... .......... .......... .......... 85% 45.9M 1s
    ## 114500K .......... .......... .......... .......... .......... 85% 18.4M 1s
    ## 114550K .......... .......... .......... .......... .......... 85% 38.8M 1s
    ## 114600K .......... .......... .......... .......... .......... 85% 17.2M 1s
    ## 114650K .......... .......... .......... .......... .......... 85% 44.1M 1s
    ## 114700K .......... .......... .......... .......... .......... 85% 16.3M 1s
    ## 114750K .......... .......... .......... .......... .......... 85% 35.9M 1s
    ## 114800K .......... .......... .......... .......... .......... 85% 19.5M 1s
    ## 114850K .......... .......... .......... .......... .......... 85% 33.3M 1s
    ## 114900K .......... .......... .......... .......... .......... 85% 15.8M 1s
    ## 114950K .......... .......... .......... .......... .......... 85% 43.7M 1s
    ## 115000K .......... .......... .......... .......... .......... 85% 15.4M 1s
    ## 115050K .......... .......... .......... .......... .......... 85% 61.5M 1s
    ## 115100K .......... .......... .......... .......... .......... 85% 13.4M 1s
    ## 115150K .......... .......... .......... .......... .......... 85% 17.9M 1s
    ## 115200K .......... .......... .......... .......... .......... 85% 36.6M 1s
    ## 115250K .......... .......... .......... .......... .......... 86% 19.9M 1s
    ## 115300K .......... .......... .......... .......... .......... 86% 29.0M 1s
    ## 115350K .......... .......... .......... .......... .......... 86% 20.9M 1s
    ## 115400K .......... .......... .......... .......... .......... 86% 31.2M 1s
    ## 115450K .......... .......... .......... .......... .......... 86% 17.7M 1s
    ## 115500K .......... .......... .......... .......... .......... 86% 41.9M 1s
    ## 115550K .......... .......... .......... .......... .......... 86% 15.7M 1s
    ## 115600K .......... .......... .......... .......... .......... 86% 52.7M 1s
    ## 115650K .......... .......... .......... .......... .......... 86% 14.9M 1s
    ## 115700K .......... .......... .......... .......... .......... 86% 54.8M 1s
    ## 115750K .......... .......... .......... .......... .......... 86% 15.2M 1s
    ## 115800K .......... .......... .......... .......... .......... 86% 44.8M 1s
    ## 115850K .......... .......... .......... .......... .......... 86% 14.2M 1s
    ## 115900K .......... .......... .......... .......... .......... 86% 85.6M 1s
    ## 115950K .......... .......... .......... .......... .......... 86% 12.9M 1s
    ## 116000K .......... .......... .......... .......... .......... 86% 13.7M 1s
    ## 116050K .......... .......... .......... .......... .......... 86% 36.6M 1s
    ## 116100K .......... .......... .......... .......... .......... 86% 21.6M 1s
    ## 116150K .......... .......... .......... .......... .......... 86% 32.2M 1s
    ## 116200K .......... .......... .......... .......... .......... 86% 13.5M 1s
    ## 116250K .......... .......... .......... .......... .......... 86% 93.9M 1s
    ## 116300K .......... .......... .......... .......... .......... 86% 9.69M 1s
    ## 116350K .......... .......... .......... .......... .......... 86% 18.1M 1s
    ## 116400K .......... .......... .......... .......... .......... 86% 42.9M 1s
    ## 116450K .......... .......... .......... .......... .......... 86% 16.2M 1s
    ## 116500K .......... .......... .......... .......... .......... 86% 45.9M 1s
    ## 116550K .......... .......... .......... .......... .......... 86% 16.6M 1s
    ## 116600K .......... .......... .......... .......... .......... 87% 45.4M 1s
    ## 116650K .......... .......... .......... .......... .......... 87% 17.9M 1s
    ## 116700K .......... .......... .......... .......... .......... 87% 44.9M 1s
    ## 116750K .......... .......... .......... .......... .......... 87% 16.0M 1s
    ## 116800K .......... .......... .......... .......... .......... 87% 45.9M 1s
    ## 116850K .......... .......... .......... .......... .......... 87% 18.7M 1s
    ## 116900K .......... .......... .......... .......... .......... 87% 42.7M 1s
    ## 116950K .......... .......... .......... .......... .......... 87% 16.5M 1s
    ## 117000K .......... .......... .......... .......... .......... 87% 44.9M 1s
    ## 117050K .......... .......... .......... .......... .......... 87% 13.5M 1s
    ## 117100K .......... .......... .......... .......... .......... 87% 82.8M 1s
    ## 117150K .......... .......... .......... .......... .......... 87% 14.4M 1s
    ## 117200K .......... .......... .......... .......... .......... 87% 12.0M 1s
    ## 117250K .......... .......... .......... .......... .......... 87% 59.9M 1s
    ## 117300K .......... .......... .......... .......... .......... 87% 15.5M 1s
    ## 117350K .......... .......... .......... .......... .......... 87% 25.3M 1s
    ## 117400K .......... .......... .......... .......... .......... 87% 23.7M 1s
    ## 117450K .......... .......... .......... .......... .......... 87% 26.2M 1s
    ## 117500K .......... .......... .......... .......... .......... 87% 25.3M 1s
    ## 117550K .......... .......... .......... .......... .......... 87% 15.1M 1s
    ## 117600K .......... .......... .......... .......... .......... 87% 23.9M 1s
    ## 117650K .......... .......... .......... .......... .......... 87% 20.9M 1s
    ## 117700K .......... .......... .......... .......... .......... 87% 26.0M 1s
    ## 117750K .......... .......... .......... .......... .......... 87% 26.2M 1s
    ## 117800K .......... .......... .......... .......... .......... 87% 18.8M 1s
    ## 117850K .......... .......... .......... .......... .......... 87% 33.6M 1s
    ## 117900K .......... .......... .......... .......... .......... 87% 18.7M 1s
    ## 117950K .......... .......... .......... .......... .......... 88% 26.0M 1s
    ## 118000K .......... .......... .......... .......... .......... 88% 21.4M 1s
    ## 118050K .......... .......... .......... .......... .......... 88% 30.5M 1s
    ## 118100K .......... .......... .......... .......... .......... 88% 20.6M 1s
    ## 118150K .......... .......... .......... .......... .......... 88% 25.7M 1s
    ## 118200K .......... .......... .......... .......... .......... 88% 23.5M 1s
    ## 118250K .......... .......... .......... .......... .......... 88% 26.1M 1s
    ## 118300K .......... .......... .......... .......... .......... 88% 22.3M 1s
    ## 118350K .......... .......... .......... .......... .......... 88% 24.1M 1s
    ## 118400K .......... .......... .......... .......... .......... 88% 23.0M 1s
    ## 118450K .......... .......... .......... .......... .......... 88% 24.3M 1s
    ## 118500K .......... .......... .......... .......... .......... 88% 23.5M 1s
    ## 118550K .......... .......... .......... .......... .......... 88% 22.5M 1s
    ## 118600K .......... .......... .......... .......... .......... 88% 25.5M 1s
    ## 118650K .......... .......... .......... .......... .......... 88% 24.5M 1s
    ## 118700K .......... .......... .......... .......... .......... 88% 21.6M 1s
    ## 118750K .......... .......... .......... .......... .......... 88% 13.9M 1s
    ## 118800K .......... .......... .......... .......... .......... 88% 35.8M 1s
    ## 118850K .......... .......... .......... .......... .......... 88% 17.7M 1s
    ## 118900K .......... .......... .......... .......... .......... 88% 46.0M 1s
    ## 118950K .......... .......... .......... .......... .......... 88% 19.2M 1s
    ## 119000K .......... .......... .......... .......... .......... 88% 33.6M 1s
    ## 119050K .......... .......... .......... .......... .......... 88% 14.4M 1s
    ## 119100K .......... .......... .......... .......... .......... 88% 16.2M 1s
    ## 119150K .......... .......... .......... .......... .......... 88% 39.5M 1s
    ## 119200K .......... .......... .......... .......... .......... 88% 12.5M 1s
    ## 119250K .......... .......... .......... .......... .......... 88% 21.3M 1s
    ## 119300K .......... .......... .......... .......... .......... 89% 38.5M 1s
    ## 119350K .......... .......... .......... .......... .......... 89% 16.5M 1s
    ## 119400K .......... .......... .......... .......... .......... 89% 36.7M 1s
    ## 119450K .......... .......... .......... .......... .......... 89% 18.9M 1s
    ## 119500K .......... .......... .......... .......... .......... 89% 36.0M 1s
    ## 119550K .......... .......... .......... .......... .......... 89% 17.9M 1s
    ## 119600K .......... .......... .......... .......... .......... 89% 40.3M 1s
    ## 119650K .......... .......... .......... .......... .......... 89% 16.4M 1s
    ## 119700K .......... .......... .......... .......... .......... 89% 45.5M 1s
    ## 119750K .......... .......... .......... .......... .......... 89% 17.3M 1s
    ## 119800K .......... .......... .......... .......... .......... 89% 41.7M 1s
    ## 119850K .......... .......... .......... .......... .......... 89% 15.9M 1s
    ## 119900K .......... .......... .......... .......... .......... 89% 56.6M 1s
    ## 119950K .......... .......... .......... .......... .......... 89% 14.6M 1s
    ## 120000K .......... .......... .......... .......... .......... 89% 69.8M 1s
    ## 120050K .......... .......... .......... .......... .......... 89% 17.9M 1s
    ## 120100K .......... .......... .......... .......... .......... 89% 38.6M 1s
    ## 120150K .......... .......... .......... .......... .......... 89% 16.6M 1s
    ## 120200K .......... .......... .......... .......... .......... 89% 25.3M 1s
    ## 120250K .......... .......... .......... .......... .......... 89% 19.5M 1s
    ## 120300K .......... .......... .......... .......... .......... 89% 33.0M 1s
    ## 120350K .......... .......... .......... .......... .......... 89% 16.1M 1s
    ## 120400K .......... .......... .......... .......... .......... 89% 16.8M 1s
    ## 120450K .......... .......... .......... .......... .......... 89% 51.9M 1s
    ## 120500K .......... .......... .......... .......... .......... 89% 13.6M 1s
    ## 120550K .......... .......... .......... .......... .......... 89% 28.2M 1s
    ## 120600K .......... .......... .......... .......... .......... 89% 17.4M 1s
    ## 120650K .......... .......... .......... .......... .......... 90% 42.2M 1s
    ## 120700K .......... .......... .......... .......... .......... 90% 16.2M 1s
    ## 120750K .......... .......... .......... .......... .......... 90% 39.0M 1s
    ## 120800K .......... .......... .......... .......... .......... 90% 16.7M 1s
    ## 120850K .......... .......... .......... .......... .......... 90% 43.6M 1s
    ## 120900K .......... .......... .......... .......... .......... 90% 16.6M 1s
    ## 120950K .......... .......... .......... .......... .......... 90% 50.1M 1s
    ## 121000K .......... .......... .......... .......... .......... 90% 16.0M 1s
    ## 121050K .......... .......... .......... .......... .......... 90% 53.8M 1s
    ## 121100K .......... .......... .......... .......... .......... 90% 12.9M 1s
    ## 121150K .......... .......... .......... .......... .......... 90% 32.9M 1s
    ## 121200K .......... .......... .......... .......... .......... 90% 15.4M 1s
    ## 121250K .......... .......... .......... .......... .......... 90% 45.7M 1s
    ## 121300K .......... .......... .......... .......... .......... 90% 16.7M 1s
    ## 121350K .......... .......... .......... .......... .......... 90% 36.3M 1s
    ## 121400K .......... .......... .......... .......... .......... 90% 16.0M 1s
    ## 121450K .......... .......... .......... .......... .......... 90% 51.0M 1s
    ## 121500K .......... .......... .......... .......... .......... 90% 15.7M 1s
    ## 121550K .......... .......... .......... .......... .......... 90% 16.6M 1s
    ## 121600K .......... .......... .......... .......... .......... 90% 29.6M 1s
    ## 121650K .......... .......... .......... .......... .......... 90% 18.4M 1s
    ## 121700K .......... .......... .......... .......... .......... 90% 49.1M 1s
    ## 121750K .......... .......... .......... .......... .......... 90% 17.0M 1s
    ## 121800K .......... .......... .......... .......... .......... 90% 37.2M 1s
    ## 121850K .......... .......... .......... .......... .......... 90% 16.4M 1s
    ## 121900K .......... .......... .......... .......... .......... 90% 45.9M 1s
    ## 121950K .......... .......... .......... .......... .......... 91% 15.2M 1s
    ## 122000K .......... .......... .......... .......... .......... 91% 50.5M 1s
    ## 122050K .......... .......... .......... .......... .......... 91% 16.4M 1s
    ## 122100K .......... .......... .......... .......... .......... 91% 45.2M 1s
    ## 122150K .......... .......... .......... .......... .......... 91% 17.6M 1s
    ## 122200K .......... .......... .......... .......... .......... 91% 47.1M 1s
    ## 122250K .......... .......... .......... .......... .......... 91% 18.0M 1s
    ## 122300K .......... .......... .......... .......... .......... 91% 39.1M 1s
    ## 122350K .......... .......... .......... .......... .......... 91% 17.5M 1s
    ## 122400K .......... .......... .......... .......... .......... 91% 43.1M 1s
    ## 122450K .......... .......... .......... .......... .......... 91% 13.8M 1s
    ## 122500K .......... .......... .......... .......... .......... 91% 66.9M 1s
    ## 122550K .......... .......... .......... .......... .......... 91% 15.4M 1s
    ## 122600K .......... .......... .......... .......... .......... 91% 17.3M 1s
    ## 122650K .......... .......... .......... .......... .......... 91% 41.9M 1s
    ## 122700K .......... .......... .......... .......... .......... 91% 18.9M 1s
    ## 122750K .......... .......... .......... .......... .......... 91% 40.7M 1s
    ## 122800K .......... .......... .......... .......... .......... 91% 18.9M 1s
    ## 122850K .......... .......... .......... .......... .......... 91% 43.8M 1s
    ## 122900K .......... .......... .......... .......... .......... 91% 19.2M 1s
    ## 122950K .......... .......... .......... .......... .......... 91% 39.3M 1s
    ## 123000K .......... .......... .......... .......... .......... 91% 18.6M 1s
    ## 123050K .......... .......... .......... .......... .......... 91% 45.3M 1s
    ## 123100K .......... .......... .......... .......... .......... 91% 18.4M 0s
    ## 123150K .......... .......... .......... .......... .......... 91% 41.9M 0s
    ## 123200K .......... .......... .......... .......... .......... 91% 15.2M 0s
    ## 123250K .......... .......... .......... .......... .......... 91% 82.0M 0s
    ## 123300K .......... .......... .......... .......... .......... 92% 17.4M 0s
    ## 123350K .......... .......... .......... .......... .......... 92% 54.2M 0s
    ## 123400K .......... .......... .......... .......... .......... 92% 14.8M 0s
    ## 123450K .......... .......... .......... .......... .......... 92%  106M 0s
    ## 123500K .......... .......... .......... .......... .......... 92% 14.5M 0s
    ## 123550K .......... .......... .......... .......... .......... 92% 19.8M 0s
    ## 123600K .......... .......... .......... .......... .......... 92% 34.3M 0s
    ## 123650K .......... .......... .......... .......... .......... 92% 18.6M 0s
    ## 123700K .......... .......... .......... .......... .......... 92% 38.2M 0s
    ## 123750K .......... .......... .......... .......... .......... 92% 16.0M 0s
    ## 123800K .......... .......... .......... .......... .......... 92% 47.1M 0s
    ## 123850K .......... .......... .......... .......... .......... 92% 19.1M 0s
    ## 123900K .......... .......... .......... .......... .......... 92% 46.2M 0s
    ## 123950K .......... .......... .......... .......... .......... 92% 16.0M 0s
    ## 124000K .......... .......... .......... .......... .......... 92% 53.9M 0s
    ## 124050K .......... .......... .......... .......... .......... 92% 17.1M 0s
    ## 124100K .......... .......... .......... .......... .......... 92% 48.0M 0s
    ## 124150K .......... .......... .......... .......... .......... 92% 17.2M 0s
    ## 124200K .......... .......... .......... .......... .......... 92% 67.2M 0s
    ## 124250K .......... .......... .......... .......... .......... 92% 12.3M 0s
    ## 124300K .......... .......... .......... .......... .......... 92% 86.4M 0s
    ## 124350K .......... .......... .......... .......... .......... 92% 16.3M 0s
    ## 124400K .......... .......... .......... .......... .......... 92% 49.6M 0s
    ## 124450K .......... .......... .......... .......... .......... 92% 13.7M 0s
    ## 124500K .......... .......... .......... .......... .......... 92% 65.1M 0s
    ## 124550K .......... .......... .......... .......... .......... 92% 13.2M 0s
    ## 124600K .......... .......... .......... .......... .......... 92% 52.0M 0s
    ## 124650K .......... .......... .......... .......... .......... 93% 12.1M 0s
    ## 124700K .......... .......... .......... .......... .......... 93% 89.2M 0s
    ## 124750K .......... .......... .......... .......... .......... 93% 13.1M 0s
    ## 124800K .......... .......... .......... .......... .......... 93% 15.1M 0s
    ## 124850K .......... .......... .......... .......... .......... 93% 53.3M 0s
    ## 124900K .......... .......... .......... .......... .......... 93% 15.7M 0s
    ## 124950K .......... .......... .......... .......... .......... 93% 45.6M 0s
    ## 125000K .......... .......... .......... .......... .......... 93% 14.8M 0s
    ## 125050K .......... .......... .......... .......... .......... 93% 49.8M 0s
    ## 125100K .......... .......... .......... .......... .......... 93% 14.3M 0s
    ## 125150K .......... .......... .......... .......... .......... 93% 16.9M 0s
    ## 125200K .......... .......... .......... .......... .......... 93% 33.5M 0s
    ## 125250K .......... .......... .......... .......... .......... 93% 15.8M 0s
    ## 125300K .......... .......... .......... .......... .......... 93% 58.3M 0s
    ## 125350K .......... .......... .......... .......... .......... 93% 14.0M 0s
    ## 125400K .......... .......... .......... .......... .......... 93% 36.5M 0s
    ## 125450K .......... .......... .......... .......... .......... 93% 20.9M 0s
    ## 125500K .......... .......... .......... .......... .......... 93% 30.2M 0s
    ## 125550K .......... .......... .......... .......... .......... 93% 14.6M 0s
    ## 125600K .......... .......... .......... .......... .......... 93% 21.0M 0s
    ## 125650K .......... .......... .......... .......... .......... 93% 26.8M 0s
    ## 125700K .......... .......... .......... .......... .......... 93% 15.0M 0s
    ## 125750K .......... .......... .......... .......... .......... 93% 34.2M 0s
    ## 125800K .......... .......... .......... .......... .......... 93% 13.4M 0s
    ## 125850K .......... .......... .......... .......... .......... 93% 43.2M 0s
    ## 125900K .......... .......... .......... .......... .......... 93% 13.5M 0s
    ## 125950K .......... .......... .......... .......... .......... 93% 22.3M 0s
    ## 126000K .......... .......... .......... .......... .......... 94% 19.3M 0s
    ## 126050K .......... .......... .......... .......... .......... 94% 26.9M 0s
    ## 126100K .......... .......... .......... .......... .......... 94% 20.1M 0s
    ## 126150K .......... .......... .......... .......... .......... 94% 30.6M 0s
    ## 126200K .......... .......... .......... .......... .......... 94% 13.1M 0s
    ## 126250K .......... .......... .......... .......... .......... 94%  103M 0s
    ## 126300K .......... .......... .......... .......... .......... 94% 18.4M 0s
    ## 126350K .......... .......... .......... .......... .......... 94% 34.6M 0s
    ## 126400K .......... .......... .......... .......... .......... 94% 11.5M 0s
    ## 126450K .......... .......... .......... .......... .......... 94% 43.8M 0s
    ## 126500K .......... .......... .......... .......... .......... 94% 15.7M 0s
    ## 126550K .......... .......... .......... .......... .......... 94% 38.7M 0s
    ## 126600K .......... .......... .......... .......... .......... 94% 11.0M 0s
    ## 126650K .......... .......... .......... .......... .......... 94% 80.1M 0s
    ## 126700K .......... .......... .......... .......... .......... 94% 13.7M 0s
    ## 126750K .......... .......... .......... .......... .......... 94% 18.7M 0s
    ## 126800K .......... .......... .......... .......... .......... 94% 25.9M 0s
    ## 126850K .......... .......... .......... .......... .......... 94% 21.4M 0s
    ## 126900K .......... .......... .......... .......... .......... 94% 29.5M 0s
    ## 126950K .......... .......... .......... .......... .......... 94% 19.1M 0s
    ## 127000K .......... .......... .......... .......... .......... 94% 27.3M 0s
    ## 127050K .......... .......... .......... .......... .......... 94% 25.6M 0s
    ## 127100K .......... .......... .......... .......... .......... 94% 27.8M 0s
    ## 127150K .......... .......... .......... .......... .......... 94% 17.8M 0s
    ## 127200K .......... .......... .......... .......... .......... 94% 35.8M 0s
    ## 127250K .......... .......... .......... .......... .......... 94% 18.8M 0s
    ## 127300K .......... .......... .......... .......... .......... 94% 30.5M 0s
    ## 127350K .......... .......... .......... .......... .......... 95% 18.3M 0s
    ## 127400K .......... .......... .......... .......... .......... 95% 41.3M 0s
    ## 127450K .......... .......... .......... .......... .......... 95% 14.3M 0s
    ## 127500K .......... .......... .......... .......... .......... 95% 21.4M 0s
    ## 127550K .......... .......... .......... .......... .......... 95% 25.4M 0s
    ## 127600K .......... .......... .......... .......... .......... 95% 18.6M 0s
    ## 127650K .......... .......... .......... .......... .......... 95% 42.1M 0s
    ## 127700K .......... .......... .......... .......... .......... 95% 16.9M 0s
    ## 127750K .......... .......... .......... .......... .......... 95% 49.7M 0s
    ## 127800K .......... .......... .......... .......... .......... 95% 15.0M 0s
    ## 127850K .......... .......... .......... .......... .......... 95% 34.0M 0s
    ## 127900K .......... .......... .......... .......... .......... 95% 19.3M 0s
    ## 127950K .......... .......... .......... .......... .......... 95% 32.1M 0s
    ## 128000K .......... .......... .......... .......... .......... 95% 18.4M 0s
    ## 128050K .......... .......... .......... .......... .......... 95% 39.8M 0s
    ## 128100K .......... .......... .......... .......... .......... 95% 17.9M 0s
    ## 128150K .......... .......... .......... .......... .......... 95% 36.7M 0s
    ## 128200K .......... .......... .......... .......... .......... 95% 16.7M 0s
    ## 128250K .......... .......... .......... .......... .......... 95% 46.7M 0s
    ## 128300K .......... .......... .......... .......... .......... 95% 16.6M 0s
    ## 128350K .......... .......... .......... .......... .......... 95% 17.7M 0s
    ## 128400K .......... .......... .......... .......... .......... 95% 37.5M 0s
    ## 128450K .......... .......... .......... .......... .......... 95% 15.0M 0s
    ## 128500K .......... .......... .......... .......... .......... 95% 41.3M 0s
    ## 128550K .......... .......... .......... .......... .......... 95% 19.2M 0s
    ## 128600K .......... .......... .......... .......... .......... 95% 33.3M 0s
    ## 128650K .......... .......... .......... .......... .......... 95% 20.7M 0s
    ## 128700K .......... .......... .......... .......... .......... 96% 31.1M 0s
    ## 128750K .......... .......... .......... .......... .......... 96% 18.4M 0s
    ## 128800K .......... .......... .......... .......... .......... 96% 15.0M 0s
    ## 128850K .......... .......... .......... .......... .......... 96% 11.8M 0s
    ## 128900K .......... .......... .......... .......... .......... 96% 46.3M 0s
    ## 128950K .......... .......... .......... .......... .......... 96% 17.1M 0s
    ## 129000K .......... .......... .......... .......... .......... 96% 43.0M 0s
    ## 129050K .......... .......... .......... .......... .......... 96% 16.3M 0s
    ## 129100K .......... .......... .......... .......... .......... 96% 46.4M 0s
    ## 129150K .......... .......... .......... .......... .......... 96% 12.8M 0s
    ## 129200K .......... .......... .......... .......... .......... 96% 78.9M 0s
    ## 129250K .......... .......... .......... .......... .......... 96% 14.6M 0s
    ## 129300K .......... .......... .......... .......... .......... 96% 15.8M 0s
    ## 129350K .......... .......... .......... .......... .......... 96% 16.7M 0s
    ## 129400K .......... .......... .......... .......... .......... 96% 34.4M 0s
    ## 129450K .......... .......... .......... .......... .......... 96% 16.9M 0s
    ## 129500K .......... .......... .......... .......... .......... 96% 42.4M 0s
    ## 129550K .......... .......... .......... .......... .......... 96% 15.7M 0s
    ## 129600K .......... .......... .......... .......... .......... 96% 44.6M 0s
    ## 129650K .......... .......... .......... .......... .......... 96% 16.8M 0s
    ## 129700K .......... .......... .......... .......... .......... 96% 44.6M 0s
    ## 129750K .......... .......... .......... .......... .......... 96% 16.6M 0s
    ## 129800K .......... .......... .......... .......... .......... 96% 54.4M 0s
    ## 129850K .......... .......... .......... .......... .......... 96% 14.5M 0s
    ## 129900K .......... .......... .......... .......... .......... 96% 17.5M 0s
    ## 129950K .......... .......... .......... .......... .......... 96% 16.6M 0s
    ## 130000K .......... .......... .......... .......... .......... 97% 15.9M 0s
    ## 130050K .......... .......... .......... .......... .......... 97% 45.7M 0s
    ## 130100K .......... .......... .......... .......... .......... 97% 16.6M 0s
    ## 130150K .......... .......... .......... .......... .......... 97% 43.4M 0s
    ## 130200K .......... .......... .......... .......... .......... 97% 16.0M 0s
    ## 130250K .......... .......... .......... .......... .......... 97% 52.4M 0s
    ## 130300K .......... .......... .......... .......... .......... 97% 16.7M 0s
    ## 130350K .......... .......... .......... .......... .......... 97% 15.4M 0s
    ## 130400K .......... .......... .......... .......... .......... 97% 38.1M 0s
    ## 130450K .......... .......... .......... .......... .......... 97% 18.4M 0s
    ## 130500K .......... .......... .......... .......... .......... 97% 35.5M 0s
    ## 130550K .......... .......... .......... .......... .......... 97% 18.6M 0s
    ## 130600K .......... .......... .......... .......... .......... 97% 39.3M 0s
    ## 130650K .......... .......... .......... .......... .......... 97% 16.5M 0s
    ## 130700K .......... .......... .......... .......... .......... 97% 43.0M 0s
    ## 130750K .......... .......... .......... .......... .......... 97% 11.8M 0s
    ## 130800K .......... .......... .......... .......... .......... 97% 85.3M 0s
    ## 130850K .......... .......... .......... .......... .......... 97% 18.0M 0s
    ## 130900K .......... .......... .......... .......... .......... 97% 52.7M 0s
    ## 130950K .......... .......... .......... .......... .......... 97% 16.9M 0s
    ## 131000K .......... .......... .......... .......... .......... 97% 49.7M 0s
    ## 131050K .......... .......... .......... .......... .......... 97% 17.8M 0s
    ## 131100K .......... .......... .......... .......... .......... 97% 41.0M 0s
    ## 131150K .......... .......... .......... .......... .......... 97% 14.3M 0s
    ## 131200K .......... .......... .......... .......... .......... 97% 60.9M 0s
    ## 131250K .......... .......... .......... .......... .......... 97% 17.8M 0s
    ## 131300K .......... .......... .......... .......... .......... 97% 46.9M 0s
    ## 131350K .......... .......... .......... .......... .......... 98% 16.0M 0s
    ## 131400K .......... .......... .......... .......... .......... 98% 47.3M 0s
    ## 131450K .......... .......... .......... .......... .......... 98% 17.4M 0s
    ## 131500K .......... .......... .......... .......... .......... 98% 35.3M 0s
    ## 131550K .......... .......... .......... .......... .......... 98% 13.7M 0s
    ## 131600K .......... .......... .......... .......... .......... 98% 19.6M 0s
    ## 131650K .......... .......... .......... .......... .......... 98% 37.2M 0s
    ## 131700K .......... .......... .......... .......... .......... 98% 19.9M 0s
    ## 131750K .......... .......... .......... .......... .......... 98% 30.7M 0s
    ## 131800K .......... .......... .......... .......... .......... 98% 17.0M 0s
    ## 131850K .......... .......... .......... .......... .......... 98% 40.6M 0s
    ## 131900K .......... .......... .......... .......... .......... 98% 16.6M 0s
    ## 131950K .......... .......... .......... .......... .......... 98% 48.7M 0s
    ## 132000K .......... .......... .......... .......... .......... 98% 14.4M 0s
    ## 132050K .......... .......... .......... .......... .......... 98% 53.5M 0s
    ## 132100K .......... .......... .......... .......... .......... 98% 16.2M 0s
    ## 132150K .......... .......... .......... .......... .......... 98% 43.3M 0s
    ## 132200K .......... .......... .......... .......... .......... 98% 16.8M 0s
    ## 132250K .......... .......... .......... .......... .......... 98% 34.6M 0s
    ## 132300K .......... .......... .......... .......... .......... 98% 20.1M 0s
    ## 132350K .......... .......... .......... .......... .......... 98% 15.5M 0s
    ## 132400K .......... .......... .......... .......... .......... 98% 22.9M 0s
    ## 132450K .......... .......... .......... .......... .......... 98% 27.1M 0s
    ## 132500K .......... .......... .......... .......... .......... 98% 26.8M 0s
    ## 132550K .......... .......... .......... .......... .......... 98% 23.0M 0s
    ## 132600K .......... .......... .......... .......... .......... 98% 21.7M 0s
    ## 132650K .......... .......... .......... .......... .......... 98% 25.6M 0s
    ## 132700K .......... .......... .......... .......... .......... 99% 21.6M 0s
    ## 132750K .......... .......... .......... .......... .......... 99% 24.6M 0s
    ## 132800K .......... .......... .......... .......... .......... 99% 23.6M 0s
    ## 132850K .......... .......... .......... .......... .......... 99% 15.7M 0s
    ## 132900K .......... .......... .......... .......... .......... 99% 38.6M 0s
    ## 132950K .......... .......... .......... .......... .......... 99% 14.7M 0s
    ## 133000K .......... .......... .......... .......... .......... 99% 15.2M 0s
    ## 133050K .......... .......... .......... .......... .......... 99% 41.0M 0s
    ## 133100K .......... .......... .......... .......... .......... 99% 18.3M 0s
    ## 133150K .......... .......... .......... .......... .......... 99% 31.0M 0s
    ## 133200K .......... .......... .......... .......... .......... 99% 16.1M 0s
    ## 133250K .......... .......... .......... .......... .......... 99% 44.0M 0s
    ## 133300K .......... .......... .......... .......... .......... 99% 17.6M 0s
    ## 133350K .......... .......... .......... .......... .......... 99% 14.7M 0s
    ## 133400K .......... .......... .......... .......... .......... 99% 53.1M 0s
    ## 133450K .......... .......... .......... .......... .......... 99% 13.9M 0s
    ## 133500K .......... .......... .......... .......... .......... 99% 14.7M 0s
    ## 133550K .......... .......... .......... .......... .......... 99% 48.6M 0s
    ## 133600K .......... .......... .......... .......... .......... 99% 15.0M 0s
    ## 133650K .......... .......... .......... .......... .......... 99% 11.0M 0s
    ## 133700K .......... .......... .......... .......... .......... 99% 51.2M 0s
    ## 133750K .......... .......... .......... .......... .......... 99% 16.0M 0s
    ## 133800K .......... .......... .......... .......... .......... 99% 46.1M 0s
    ## 133850K .......... .......... .......... .......... .......... 99% 10.5M 0s
    ## 133900K .......... .......... .......... .......... .......... 99% 67.4M 0s
    ## 133950K .......... .......... .......... .......... .......... 99% 11.0M 0s
    ## 134000K .......... .......... .......... .......... .......... 99% 92.1M 0s
    ## 134050K .......... .....                                      100% 5.78M=6.1s
    ## 
    ## 2023-01-18 10:07:14 (21.4 MB/s) - ‘silva_nr99_v138.1_train_set.fa.gz?download=1.6’ saved [137283333/137283333]

``` r
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz?download=1", multithread=TRUE)
```

``` r
#Inspecter les affectations taxonomiques
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class        
    ## [1,] "Bacteria" NA             NA           
    ## [2,] "Bacteria" "Firmicutes"   "Bacilli"    
    ## [3,] "Bacteria" NA             NA           
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia"
    ## [5,] "Bacteria" "Firmicutes"   "Clostridia" 
    ## [6,] "Bacteria" "Firmicutes"   "Clostridia" 
    ##      Order                                 Family                     
    ## [1,] NA                                    NA                         
    ## [2,] "Erysipelotrichales"                  "Erysipelatoclostridiaceae"
    ## [3,] NA                                    NA                         
    ## [4,] "Bacteroidales"                       "Muribaculaceae"           
    ## [5,] "Clostridiales"                       "Clostridiaceae"           
    ## [6,] "Peptostreptococcales-Tissierellales" "Anaerovoracaceae"         
    ##      Genus                        
    ## [1,] NA                           
    ## [2,] "Coprobacillus"              
    ## [3,] NA                           
    ## [4,] NA                           
    ## [5,] "Clostridium sensu stricto 1"
    ## [6,] "Mogibacterium"

``` r
library(phyloseq); packageVersion("phyloseq")
```

    ## [1] '1.40.0'

``` r
library(Biostrings); packageVersion("Biostrings")
```

    ## [1] '2.64.1'

``` r
library(ggplot2); packageVersion("ggplot2")
```

    ## [1] '3.3.6'

``` r
theme_set(theme_bw())
```

# Les séquences d’ADN mise en ligne par l’auteur étaient dans un seul dossier (Ne sont pas séparées selon les régions d’étude), cela m’a empéché d’établir un graphique qui montre la taxonomie des microorganismes selon la région.
