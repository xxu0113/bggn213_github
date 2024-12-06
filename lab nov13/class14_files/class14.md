# class14:rnaseq mini


## import data

counts metadata what DESeq calls colDATA as it describes the columns in
the counts \## data cleanup

``` r
counts <- read.csv("GSE37704_featurecounts.csv", row.names=1)
metadata <-  read.csv("GSE37704_metadata.csv")
head(counts)
```

                    length SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ENSG00000186092    918         0         0         0         0         0
    ENSG00000279928    718         0         0         0         0         0
    ENSG00000279457   1982        23        28        29        29        28
    ENSG00000278566    939         0         0         0         0         0
    ENSG00000273547    939         0         0         0         0         0
    ENSG00000187634   3214       124       123       205       207       212
                    SRR493371
    ENSG00000186092         0
    ENSG00000279928         0
    ENSG00000279457        46
    ENSG00000278566         0
    ENSG00000273547         0
    ENSG00000187634       258

``` r
head(metadata)
```

             id     condition
    1 SRR493366 control_sirna
    2 SRR493367 control_sirna
    3 SRR493368 control_sirna
    4 SRR493369      hoxa1_kd
    5 SRR493370      hoxa1_kd
    6 SRR493371      hoxa1_kd

we want the columns in counts to match the rows in metadata

``` r
colnames(counts)
```

    [1] "length"    "SRR493366" "SRR493367" "SRR493368" "SRR493369" "SRR493370"
    [7] "SRR493371"

``` r
metadata$id
```

    [1] "SRR493366" "SRR493367" "SRR493368" "SRR493369" "SRR493370" "SRR493371"

we can get rid of the first column in counts

``` r
countData <- counts[,-1]
head(countData)
```

                    SRR493366 SRR493367 SRR493368 SRR493369 SRR493370 SRR493371
    ENSG00000186092         0         0         0         0         0         0
    ENSG00000279928         0         0         0         0         0         0
    ENSG00000279457        23        28        29        29        28        46
    ENSG00000278566         0         0         0         0         0         0
    ENSG00000273547         0         0         0         0         0         0
    ENSG00000187634       124       123       205       207       212       258

``` r
colnames(countData)
```

    [1] "SRR493366" "SRR493367" "SRR493368" "SRR493369" "SRR493370" "SRR493371"

``` r
colnames(countData) == metadata$id
```

    [1] TRUE TRUE TRUE TRUE TRUE TRUE

``` r
all(colnames(countData) == metadata$id)
```

    [1] TRUE

``` r
x <- c(T,F,T,T)
if(all(x)){
  cat("Me happy")
}  else{
  cat("Me no happy")
}
```

    Me no happy

\##filter outzero counts it is practice to remove any genes that we have
no data for.

``` r
tp.keep.inds <- rowSums(countData) > 0
cleanCounts <- countData[tp.keep.inds,]
head(cleanCounts)
```

                    SRR493366 SRR493367 SRR493368 SRR493369 SRR493370 SRR493371
    ENSG00000279457        23        28        29        29        28        46
    ENSG00000187634       124       123       205       207       212       258
    ENSG00000188976      1637      1831      2383      1226      1326      1504
    ENSG00000187961       120       153       180       236       255       357
    ENSG00000187583        24        48        65        44        48        64
    ENSG00000187642         4         9        16        14        16        16

## setup fro DESeup

``` r
#/ message: false
library(DESeq2)
```

    Loading required package: S4Vectors

    Loading required package: stats4

    Loading required package: BiocGenerics


    Attaching package: 'BiocGenerics'

    The following objects are masked from 'package:stats':

        IQR, mad, sd, var, xtabs

    The following objects are masked from 'package:base':

        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
        table, tapply, union, unique, unsplit, which.max, which.min


    Attaching package: 'S4Vectors'

    The following object is masked from 'package:utils':

        findMatches

    The following objects are masked from 'package:base':

        expand.grid, I, unname

    Loading required package: IRanges

    Loading required package: GenomicRanges

    Loading required package: GenomeInfoDb

    Loading required package: SummarizedExperiment

    Loading required package: MatrixGenerics

    Loading required package: matrixStats


    Attaching package: 'MatrixGenerics'

    The following objects are masked from 'package:matrixStats':

        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars

    Loading required package: Biobase

    Welcome to Bioconductor

        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.


    Attaching package: 'Biobase'

    The following object is masked from 'package:MatrixGenerics':

        rowMedians

    The following objects are masked from 'package:matrixStats':

        anyMissing, rowMedians

``` r
dds <- DESeqDataSetFromMatrix(countData=cleanCounts,
                       colData=metadata,
                       design= ~condition)
```

    Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    design formula are characters, converting to factors

\##DESeq

``` r
dds <- DESeq(dds)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

``` r
res <-results(dds)
```

\##Inspect Restuls

``` r
head(res)
```

    log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    Wald test p-value: condition hoxa1 kd vs control sirna 
    DataFrame with 6 rows and 6 columns
                     baseMean log2FoldChange     lfcSE       stat      pvalue
                    <numeric>      <numeric> <numeric>  <numeric>   <numeric>
    ENSG00000279457   29.9136      0.1792571 0.3248216   0.551863 5.81042e-01
    ENSG00000187634  183.2296      0.4264571 0.1402658   3.040350 2.36304e-03
    ENSG00000188976 1651.1881     -0.6927205 0.0548465 -12.630158 1.43989e-36
    ENSG00000187961  209.6379      0.7297556 0.1318599   5.534326 3.12428e-08
    ENSG00000187583   47.2551      0.0405765 0.2718928   0.149237 8.81366e-01
    ENSG00000187642   11.9798      0.5428105 0.5215599   1.040744 2.97994e-01
                           padj
                      <numeric>
    ENSG00000279457 6.86555e-01
    ENSG00000187634 5.15718e-03
    ENSG00000188976 1.76549e-35
    ENSG00000187961 1.13413e-07
    ENSG00000187583 9.19031e-01
    ENSG00000187642 4.03379e-01

\##data viz

``` r
plot(res$log2FoldChange,-log(res$padj))
```

![](class14_files/figure-commonmark/unnamed-chunk-14-1.png)

\##pathway analysis

``` r
head(res)
```

    log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    Wald test p-value: condition hoxa1 kd vs control sirna 
    DataFrame with 6 rows and 6 columns
                     baseMean log2FoldChange     lfcSE       stat      pvalue
                    <numeric>      <numeric> <numeric>  <numeric>   <numeric>
    ENSG00000279457   29.9136      0.1792571 0.3248216   0.551863 5.81042e-01
    ENSG00000187634  183.2296      0.4264571 0.1402658   3.040350 2.36304e-03
    ENSG00000188976 1651.1881     -0.6927205 0.0548465 -12.630158 1.43989e-36
    ENSG00000187961  209.6379      0.7297556 0.1318599   5.534326 3.12428e-08
    ENSG00000187583   47.2551      0.0405765 0.2718928   0.149237 8.81366e-01
    ENSG00000187642   11.9798      0.5428105 0.5215599   1.040744 2.97994e-01
                           padj
                      <numeric>
    ENSG00000279457 6.86555e-01
    ENSG00000187634 5.15718e-03
    ENSG00000188976 1.76549e-35
    ENSG00000187961 1.13413e-07
    ENSG00000187583 9.19031e-01
    ENSG00000187642 4.03379e-01

\##anotation of genes

first I need to translate my ensemble id in my res object. to entrez and
gene symbol formats

for this I will use the AnnotationDbi package and it is mapids()
function,

``` r
library(AnnotationDbi)
library(org.Hs.eg.db)
```

``` r
columns(org.Hs.eg.db)
```

     [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
     [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
    [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
    [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
    [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    [26] "UNIPROT"     

lets map it to symbol, entrezid, genename

``` r
res$genename <- mapIds(org.Hs.eg.db,keys=rownames(res),keytype="ENSEMBL",column="GENENAME")
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$entrez_id <- mapIds(org.Hs.eg.db,keys=rownames(res),keytype="ENSEMBL",column="ENTREZID")
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$symbol <-mapIds(org.Hs.eg.db,keys=rownames(res),keytype="ENSEMBL",column="SYMBOL")
```

    'select()' returned 1:many mapping between keys and columns

``` r
head(res)
```

    log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    Wald test p-value: condition hoxa1 kd vs control sirna 
    DataFrame with 6 rows and 9 columns
                     baseMean log2FoldChange     lfcSE       stat      pvalue
                    <numeric>      <numeric> <numeric>  <numeric>   <numeric>
    ENSG00000279457   29.9136      0.1792571 0.3248216   0.551863 5.81042e-01
    ENSG00000187634  183.2296      0.4264571 0.1402658   3.040350 2.36304e-03
    ENSG00000188976 1651.1881     -0.6927205 0.0548465 -12.630158 1.43989e-36
    ENSG00000187961  209.6379      0.7297556 0.1318599   5.534326 3.12428e-08
    ENSG00000187583   47.2551      0.0405765 0.2718928   0.149237 8.81366e-01
    ENSG00000187642   11.9798      0.5428105 0.5215599   1.040744 2.97994e-01
                           padj               genename   entrez_id      symbol
                      <numeric>            <character> <character> <character>
    ENSG00000279457 6.86555e-01                     NA          NA          NA
    ENSG00000187634 5.15718e-03 sterile alpha motif ..      148398      SAMD11
    ENSG00000188976 1.76549e-35 NOC2 like nucleolar ..       26155       NOC2L
    ENSG00000187961 1.13413e-07 kelch like family me..      339451      KLHL17
    ENSG00000187583 9.19031e-01 pleckstrin homology ..       84069     PLEKHN1
    ENSG00000187642 4.03379e-01 PPARGC1 and ESRR ind..       84808       PERM1

Before going further, lets focus in on a subset of top genes, we can use
as a starting point log2FC of +2/-2 and a adjusted p value of 0.05.

``` r
top.inds <- (abs(res$log2FoldChange) >2) & (abs(res$padj) < 0.05)
top.inds[is.na(top.inds)] <- FALSE
```

lets save our top genes to file

``` r
top.genes <- res[top.inds, ]
write.csv(top.genes, file = "top_genes.csv", row.names = FALSE)
```

``` r
library(gage)
```

``` r
library(gageData)
library(pathview)
```

    ##############################################################################
    Pathview is an open source software package distributed under GNU General
    Public License version 3 (GPLv3). Details of GPLv3 is available at
    http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    formally cite the original Pathview paper (not just mention it) in publications
    or products. For details, do citation("pathview") within R.

    The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ##############################################################################

``` r
data(kegg.sets.hs)
data(sigmet.idx.hs)

# Focus on signaling and metabolic pathways only
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

# Examine the first 3 pathways
head(kegg.sets.hs, 3)
```

    $`hsa00232 Caffeine metabolism`
    [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   

    $`hsa00983 Drug metabolism - other enzymes`
     [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"   "1551"  
     [9] "1553"   "1576"   "1577"   "1806"   "1807"   "1890"   "221223" "2990"  
    [17] "3251"   "3614"   "3615"   "3704"   "51733"  "54490"  "54575"  "54576" 
    [25] "54577"  "54578"  "54579"  "54600"  "54657"  "54658"  "54659"  "54963" 
    [33] "574537" "64816"  "7083"   "7084"   "7172"   "7363"   "7364"   "7365"  
    [41] "7366"   "7367"   "7371"   "7372"   "7378"   "7498"   "79799"  "83549" 
    [49] "8824"   "8833"   "9"      "978"   

    $`hsa00230 Purine metabolism`
      [1] "100"    "10201"  "10606"  "10621"  "10622"  "10623"  "107"    "10714" 
      [9] "108"    "10846"  "109"    "111"    "11128"  "11164"  "112"    "113"   
     [17] "114"    "115"    "122481" "122622" "124583" "132"    "158"    "159"   
     [25] "1633"   "171568" "1716"   "196883" "203"    "204"    "205"    "221823"
     [33] "2272"   "22978"  "23649"  "246721" "25885"  "2618"   "26289"  "270"   
     [41] "271"    "27115"  "272"    "2766"   "2977"   "2982"   "2983"   "2984"  
     [49] "2986"   "2987"   "29922"  "3000"   "30833"  "30834"  "318"    "3251"  
     [57] "353"    "3614"   "3615"   "3704"   "377841" "471"    "4830"   "4831"  
     [65] "4832"   "4833"   "4860"   "4881"   "4882"   "4907"   "50484"  "50940" 
     [73] "51082"  "51251"  "51292"  "5136"   "5137"   "5138"   "5139"   "5140"  
     [81] "5141"   "5142"   "5143"   "5144"   "5145"   "5146"   "5147"   "5148"  
     [89] "5149"   "5150"   "5151"   "5152"   "5153"   "5158"   "5167"   "5169"  
     [97] "51728"  "5198"   "5236"   "5313"   "5315"   "53343"  "54107"  "5422"  
    [105] "5424"   "5425"   "5426"   "5427"   "5430"   "5431"   "5432"   "5433"  
    [113] "5434"   "5435"   "5436"   "5437"   "5438"   "5439"   "5440"   "5441"  
    [121] "5471"   "548644" "55276"  "5557"   "5558"   "55703"  "55811"  "55821" 
    [129] "5631"   "5634"   "56655"  "56953"  "56985"  "57804"  "58497"  "6240"  
    [137] "6241"   "64425"  "646625" "654364" "661"    "7498"   "8382"   "84172" 
    [145] "84265"  "84284"  "84618"  "8622"   "8654"   "87178"  "8833"   "9060"  
    [153] "9061"   "93034"  "953"    "9533"   "954"    "955"    "956"    "957"   
    [161] "9583"   "9615"  

the gage function wants a vector of importance as input with gene names
as labels kegg speaks entrez

``` r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez_id
head(foldchanges)
```

           <NA>      148398       26155      339451       84069       84808 
     0.17925708  0.42645712 -0.69272046  0.72975561  0.04057653  0.54281049 

run gage with these values

``` r
keggres = gage(foldchanges, gsets=kegg.sets.hs)
```

``` r
head(keggres$less)
```

                                             p.geomean stat.mean        p.val
    hsa04110 Cell cycle                   8.995727e-06 -4.378644 8.995727e-06
    hsa03030 DNA replication              9.424076e-05 -3.951803 9.424076e-05
    hsa03013 RNA transport                1.246882e-03 -3.059466 1.246882e-03
    hsa03440 Homologous recombination     3.066756e-03 -2.852899 3.066756e-03
    hsa04114 Oocyte meiosis               3.784520e-03 -2.698128 3.784520e-03
    hsa00010 Glycolysis / Gluconeogenesis 8.961413e-03 -2.405398 8.961413e-03
                                                q.val set.size         exp1
    hsa04110 Cell cycle                   0.001448312      121 8.995727e-06
    hsa03030 DNA replication              0.007586381       36 9.424076e-05
    hsa03013 RNA transport                0.066915974      144 1.246882e-03
    hsa03440 Homologous recombination     0.121861535       28 3.066756e-03
    hsa04114 Oocyte meiosis               0.121861535      102 3.784520e-03
    hsa00010 Glycolysis / Gluconeogenesis 0.212222694       53 8.961413e-03

``` r
pathview(foldchanges,pathway.id="hsa04110")
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/xiaowen/Desktop/lab nov13

    Info: Writing image file hsa04110.pathview.png

``` r
library(gage)
library(gageData)

data(go.sets.hs)
data(go.subs.hs)

gobpsets=go.sets.hs[go.subs.hs$BP]
gores <- gage(foldchanges,gsets=gobpsets)
head(gores$less)
```

                                                p.geomean stat.mean        p.val
    GO:0048285 organelle fission             1.536227e-15 -8.063910 1.536227e-15
    GO:0000280 nuclear division              4.286961e-15 -7.939217 4.286961e-15
    GO:0007067 mitosis                       4.286961e-15 -7.939217 4.286961e-15
    GO:0000087 M phase of mitotic cell cycle 1.169934e-14 -7.797496 1.169934e-14
    GO:0007059 chromosome segregation        2.028624e-11 -6.878340 2.028624e-11
    GO:0000236 mitotic prometaphase          1.729553e-10 -6.695966 1.729553e-10
                                                    q.val set.size         exp1
    GO:0048285 organelle fission             5.841698e-12      376 1.536227e-15
    GO:0000280 nuclear division              5.841698e-12      352 4.286961e-15
    GO:0007067 mitosis                       5.841698e-12      352 4.286961e-15
    GO:0000087 M phase of mitotic cell cycle 1.195672e-11      362 1.169934e-14
    GO:0007059 chromosome segregation        1.658603e-08      142 2.028624e-11
    GO:0000236 mitotic prometaphase          1.178402e-07       84 1.729553e-10

to run reactome

``` r
sig_genes <- res[res$padj <= 0.05 & !is.na(res$padj), "symbol"]
print(paste("Total number of significant genes:", length(sig_genes)))
```

    [1] "Total number of significant genes: 8147"

``` r
sig_genes
```

    ENSG00000187634 ENSG00000188976 ENSG00000187961 ENSG00000188290 ENSG00000187608 
           "SAMD11"         "NOC2L"        "KLHL17"          "HES4"         "ISG15" 
    ENSG00000188157 ENSG00000078808 ENSG00000176022 ENSG00000160087 ENSG00000162572 
             "AGRN"          "SDF4"       "B3GALT6"        "UBE2J2"        "SCNN1D" 
    ENSG00000131584 ENSG00000127054 ENSG00000107404 ENSG00000162576 ENSG00000221978 
            "ACAP3"        "INTS11"          "DVL1"         "MXRA8"         "CCNL2" 
    ENSG00000160072 ENSG00000197785 ENSG00000228594 ENSG00000197530 ENSG00000189339 
           "ATAD3B"        "ATAD3A"        "FNDC10"          "MIB2"      "SLC35E2B" 
    ENSG00000215790 ENSG00000162585 ENSG00000157933 ENSG00000157873 ENSG00000157870 
                 NA        "FAAP20"           "SKI"      "TNFRSF14"        "PRXL2B" 
    ENSG00000158109 ENSG00000130764 ENSG00000116198 ENSG00000131697 ENSG00000116237 
           "TPRG1L"        "LRRC47"        "CEP104"         "NPHP4"          "ICMT" 
    ENSG00000097021 ENSG00000069812 ENSG00000215788 ENSG00000171680 ENSG00000162408 
            "ACOT7"          "HES2"      "TNFRSF25"       "PLEKHG5"          "NOL9" 
    ENSG00000162413 ENSG00000041988 ENSG00000007923 ENSG00000049245 ENSG00000116288 
           "KLHL21"         "THAP3"       "DNAJC11"         "VAMP3"         "PARK7" 
    ENSG00000116285 ENSG00000142599 ENSG00000074800 ENSG00000180758 ENSG00000049239 
           "ERRFI1"          "RERE"          "ENO1"        "GPR157"          "H6PD" 
    ENSG00000171621 ENSG00000171608 ENSG00000171603 ENSG00000178585 ENSG00000162441 
            "SPSB1"        "PIK3CD"        "CLSTN1"      "CTNNBIP1"          "LZIC" 
    ENSG00000173614 ENSG00000130939 ENSG00000054523 ENSG00000142657 ENSG00000160049 
           "NMNAT1"         "UBE4B"         "KIF1B"           "PGD"          "DFFA" 
    ENSG00000120948 ENSG00000116649 ENSG00000171824 ENSG00000120942 ENSG00000116661 
           "TARDBP"           "SRM"       "EXOSC10"        "UBIAD1"         "FBXO2" 
    ENSG00000116670 ENSG00000162490 ENSG00000177674 ENSG00000177000 ENSG00000011021 
           "MAD2L2"        "DRAXIN"        "AGTRAP"         "MTHFR"         "CLCN6" 
    ENSG00000116688 ENSG00000116691 ENSG00000028137 ENSG00000048707 ENSG00000162496 
             "MFN2"          "MIIP"      "TNFRSF1B"        "VPS13D"         "DHRS3" 
    ENSG00000162493 ENSG00000116731 ENSG00000189337 ENSG00000171729 ENSG00000142634 
             "PDPN"         "PRDM2"          "KAZN"        "TMEM51"         "EFHD2" 
    ENSG00000116138 ENSG00000197312 ENSG00000116786 ENSG00000162458 ENSG00000173641 
          "DNAJC16"          "DDI2"       "PLEKHM2"        "FBLIM1"         "HSPB7" 
    ENSG00000142627 ENSG00000132881 ENSG00000055070 ENSG00000058453 ENSG00000117122 
            "EPHA2"       "CPLANE2"         "SZRD1"         "CROCC"         "MFAP2" 
    ENSG00000159363 ENSG00000117118 ENSG00000179051 ENSG00000074964 ENSG00000159423 
          "ATP13A2"          "SDHB"          "RCC2"     "ARHGEF10L"       "ALDH4A1" 
    ENSG00000127463 ENSG00000053372 ENSG00000040487 ENSG00000077549 ENSG00000158747 
             "EMC1"         "MRTO4"       "SLC66A1"         "CAPZB"          "NBL1" 
    ENSG00000162542 ENSG00000162545 ENSG00000158825 ENSG00000158828 ENSG00000244038 
            "TMCO4"       "CAMK2N1"           "CDA"         "PINK1"         "DDOST" 
    ENSG00000117245 ENSG00000189410 ENSG00000127483 ENSG00000075151 ENSG00000117298 
            "KIF17"         "SH2D5"        "HP1BP3"        "EIF4G3"          "ECE1" 
    ENSG00000162551 ENSG00000076864 ENSG00000142798 ENSG00000184677 ENSG00000004487 
             "ALPL"       "RAP1GAP"         "HSPG2"        "ZBTB40"         "KDM1A" 
    ENSG00000169641 ENSG00000125944 ENSG00000125945 ENSG00000204219 ENSG00000007968 
            "LUZP1"        "HNRNPR"        "ZNF436"         "TCEA3"          "E2F2" 
    ENSG00000117318 ENSG00000142676 ENSG00000011007 ENSG00000057757 ENSG00000011009 
              "ID3"         "RPL11"          "ELOA"        "PITHD1"        "LYPLA2" 
    ENSG00000117308 ENSG00000117305 ENSG00000179163 ENSG00000189266 ENSG00000188529 
             "GALE"         "HMGCL"         "FUCA1"         "PNRC2"        "SRSF10" 
    ENSG00000001460 ENSG00000001461 ENSG00000117602 ENSG00000133226 ENSG00000169504 
            "STPG1"        "NIPAL3"         "RCAN3"         "SRRM1"         "CLIC4" 
    ENSG00000020633 ENSG00000117614 ENSG00000117616 ENSG00000183726 ENSG00000204178 
            "RUNX3"          "SYF2"         "RSRP1"       "TMEM50A"         "MACO1" 
    ENSG00000117643 ENSG00000162430 ENSG00000127423 ENSG00000182749 ENSG00000117632 
           "MAN1C1"       "SELENON"         "AUNIP"         "PAQR7"         "STMN1" 
    ENSG00000158006 ENSG00000158008 ENSG00000130695 ENSG00000142669 ENSG00000117682 
           "PAFAH2"         "EXTL1"         "CEP85"      "SH3BGRL3"         "DHDDS" 
    ENSG00000198830 ENSG00000117676 ENSG00000117713 ENSG00000060642 ENSG00000204160 
            "HMGN2"       "RPS6KA1"        "ARID1A"          "PIGV"       "ZDHHC18" 
    ENSG00000198746 ENSG00000090273 ENSG00000158246 ENSG00000142784 ENSG00000186501 
          "GPATCH3"          "NUDC"        "TENT5B"         "WDTC1"       "TMEM222" 
    ENSG00000142733 ENSG00000181773 ENSG00000158195 ENSG00000126705 ENSG00000117751 
           "MAP3K6"          "GPR3"         "WASF2"         "AHDC1"        "PPP1R8" 
    ENSG00000130775 ENSG00000117748 ENSG00000158156 ENSG00000158161 ENSG00000126698 
          "THEMIS2"          "RPA2"          "XKR8"          "EYA3"        "DNAJC8" 
    ENSG00000130772 ENSG00000204138 ENSG00000120656 ENSG00000162419 ENSG00000198492 
            "MED18"       "PHACTR4"         "TAF12"         "GMEB1"        "YTHDF2" 
    ENSG00000159023 ENSG00000253304 ENSG00000116350 ENSG00000060656 ENSG00000162512 
            "EPB41"      "TMEM200B"         "SRSF4"         "PTPRU"          "SDC3" 
    ENSG00000134644 ENSG00000060688 ENSG00000121766 ENSG00000168528 ENSG00000162517 
             "PUM1"       "SNRNP40"       "ZCCHC17"       "SERINC2"          "PEF1" 
    ENSG00000084636 ENSG00000134668 ENSG00000121774 ENSG00000025800 ENSG00000084652 
          "COL16A1"        "SPOCD1"       "KHDRBS1"         "KPNA6"         "TXLNA" 
    ENSG00000084623 ENSG00000116478 ENSG00000175130 ENSG00000225828 ENSG00000160058 
            "EIF3I"         "HDAC1"      "MARCKSL1"       "FAM229A"         "BSDC1" 
    ENSG00000176261 ENSG00000162521 ENSG00000162520 ENSG00000162522 ENSG00000134684 
          "ZBTB8OS"         "RBBP4"          "SYNC"         "NHSL3"         "YARS1" 
    ENSG00000116497 ENSG00000121900 ENSG00000116514 ENSG00000004455 ENSG00000134686 
          "S100PBP"        "TMEM54"        "RNF19B"           "AK2"          "PHC2" 
    ENSG00000163866 ENSG00000116560 ENSG00000146463 ENSG00000142687 ENSG00000020129 
           "SMIM12"          "SFPQ"         "ZMYM4"     "KIAA0319L"          "NCDN" 
    ENSG00000142686 ENSG00000092853 ENSG00000134698 ENSG00000171812 ENSG00000214193 
         "C1orf216"         "CLSPN"          "AGO4"        "COL8A2"        "SH3D21" 
    ENSG00000142694 ENSG00000196182 ENSG00000116885 ENSG00000163874 ENSG00000163875 
            "EVA1B"         "STK40"         "OSCP1"       "ZC3H12A"         "MEAF6" 
    ENSG00000163877 ENSG00000134697 ENSG00000134690 ENSG00000196449 ENSG00000197982 
            "SNIP1"          "GNL2"         "CDCA8"          "YRDC"      "C1orf122" 
    ENSG00000188786 ENSG00000204084 ENSG00000183431 ENSG00000183386 ENSG00000116954 
             "MTF1"        "INPP5B"         "SF3A3"          "FHL3"         "RRAGC" 
    ENSG00000214114 ENSG00000158315 ENSG00000174574 ENSG00000168653 ENSG00000090621 
            "MYCBP"        "RHBDL2"       "AKIRIN1"        "NDUFS5"        "PABPC4" 
    ENSG00000163909 ENSG00000084072 ENSG00000043514 ENSG00000168389 ENSG00000131236 
             "HEYL"          "PPIE"         "TRIT1"        "MFSD2A"          "CAP1" 
    ENSG00000131238 ENSG00000187801 ENSG00000187815 ENSG00000179862 ENSG00000171793 
             "PPT1"        "ZFP69B"         "ZFP69"        "CITED4"         "CTPS1" 
    ENSG00000010803 ENSG00000127124 ENSG00000127125 ENSG00000171960 ENSG00000065978 
            "SCMH1"        "HIVEP3"          "PPCS"          "PPIH"          "YBX1" 
    ENSG00000117385 ENSG00000177868 ENSG00000117394 ENSG00000117395 ENSG00000066056 
             "P3H1"          "SVBP"        "SLC2A1"      "EBNA1BP2"          "TIE1" 
    ENSG00000117399 ENSG00000066322 ENSG00000159479 ENSG00000142949 ENSG00000126091 
            "CDC20"        "ELOVL1"          "MED8"         "PTPRF"       "ST3GAL3" 
    ENSG00000132768 ENSG00000117410 ENSG00000117411 ENSG00000196517 ENSG00000187147 
             "DPH2"       "ATP6V0B"       "B4GALT2"        "SLC6A9"        "RNF220" 
    ENSG00000142945 ENSG00000142937 ENSG00000173846 ENSG00000188396 ENSG00000222009 
            "KIF2C"          "RPS8"          "PLK3"        "DYNLT4"        "BTBD19" 
    ENSG00000117425 ENSG00000126107 ENSG00000126088 ENSG00000162415 ENSG00000070759 
            "PTCH2"        "HECTD3"          "UROD"        "ZSWIM5"         "TESK2" 
    ENSG00000132763 ENSG00000117450 ENSG00000132780 ENSG00000159592 ENSG00000159596 
           "MMACHC"         "PRDX1"          "NASP"       "GPBP1L1"        "TMEM69" 
    ENSG00000197429 ENSG00000086015 ENSG00000117461 ENSG00000085998 ENSG00000085999 
              "IPP"         "MAST2"        "PIK3R3"       "POMGNT1"        "RAD54L" 
    ENSG00000117481 ENSG00000123472 ENSG00000159658 ENSG00000123473 ENSG00000132122 
            "NSUN4"        "ATPAF1"       "EFCAB14"          "STIL"        "SPATA6" 
    ENSG00000185104 ENSG00000123080 ENSG00000123091 ENSG00000085832 ENSG00000117859 
             "FAF1"        "CDKN2C"         "RNF11"         "EPS15"        "OSBPL9" 
    ENSG00000169213 ENSG00000117862 ENSG00000198841 ENSG00000157077 ENSG00000154222 
            "RAB3B"       "TXNDC12"         "KTI12"        "ZFYVE9"        "CC2D1B" 
    ENSG00000085840 ENSG00000162377 ENSG00000203995 ENSG00000116171 ENSG00000174348 
             "ORC1"          "COA7"        "ZYG11A"          "SCP2"          "PODN" 
    ENSG00000157184 ENSG00000162384 ENSG00000162385 ENSG00000157193 ENSG00000058804 
             "CPT2"          "CZIB"         "MAGOH"          "LRP8"          "NDC1" 
    ENSG00000058799 ENSG00000116212 ENSG00000116209 ENSG00000215883 ENSG00000116221 
            "YIPF1"        "LRRC42"        "TMEM59"        "CYB5RL"        "MRPL37" 
    ENSG00000162390 ENSG00000243725 ENSG00000162396 ENSG00000116133 ENSG00000162402 
           "ACOT11"          "TTC4"         "PARS2"        "DHCR24"         "USP24" 
    ENSG00000162407 ENSG00000162409 ENSG00000173406 ENSG00000162601 ENSG00000177606 
            "PLPP3"        "PRKAA2"          "DAB1"         "MYSM1"           "JUN" 
    ENSG00000162599 ENSG00000162604 ENSG00000132849 ENSG00000162607 ENSG00000116641 
             "NFIA"         "TM2D1"          "PATJ"          "USP1"         "DOCK7" 
    ENSG00000088035 ENSG00000142856 ENSG00000079739 ENSG00000185483 ENSG00000162437 
             "ALG6"       "ITGB3BP"          "PGM1"          "ROR1"        "RAVER2" 
    ENSG00000162434 ENSG00000162433 ENSG00000213625 ENSG00000184588 ENSG00000118473 
             "JAK1"           "AK4"        "LEPROT"         "PDE4B"         "SGIP1" 
    ENSG00000152763 ENSG00000198160 ENSG00000116704 ENSG00000142864 ENSG00000116717 
            "DNAI4"         "MIER1"       "SLC35D1"        "SERBP1"       "GADD45A" 
    ENSG00000172380 ENSG00000162595 ENSG00000116729 ENSG00000024526 ENSG00000116754 
            "GNG12"        "DIRAS3"           "WLS"        "DEPDC1"        "SRSF11" 
    ENSG00000132485 ENSG00000172260 ENSG00000254685 ENSG00000116791 ENSG00000162623 
           "ZRANB2"         "NEGR1"          "FPGT"          "CRYZ"          "TYW3" 
    ENSG00000117054 ENSG00000117069 ENSG00000142892 ENSG00000154027 ENSG00000077254 
            "ACADM"    "ST6GALNAC5"          "PIGK"           "AK5"         "USP33" 
    ENSG00000180488 ENSG00000162614 ENSG00000162613 ENSG00000162616 ENSG00000137960 
            "MIGA1"          "NEXN"         "FUBP1"        "DNAJB4"         "GIPC2" 
    ENSG00000122420 ENSG00000137959 ENSG00000162618 ENSG00000117114 ENSG00000137941 
            "PTGFR"        "IFI44L"        "ADGRL4"        "ADGRL2"         "TTLL7" 
    ENSG00000142875 ENSG00000117133 ENSG00000174021 ENSG00000117151 ENSG00000117155 
           "PRKACB"          "RPF1"          "GNG5"          "CTBS"        "SSX2IP" 
    ENSG00000171517 ENSG00000162643 ENSG00000097096 ENSG00000162642 ENSG00000142867 
            "LPAR3"         "DNAI3"         "SYDE2"       "C1orf52"         "BCL10" 
    ENSG00000153904 ENSG00000142871 ENSG00000171502 ENSG00000122417 ENSG00000137975 
            "DDAH1"          "CCN1"       "COL24A1"         "ODF2L"         "CLCA2" 
    ENSG00000097033 ENSG00000183291 ENSG00000143013 ENSG00000065243 ENSG00000137947 
          "SH3GLB1"       "SELENOF"          "LMO4"          "PKN2"         "GTF2B" 
    ENSG00000137944 ENSG00000213516 ENSG00000117228 ENSG00000162645 ENSG00000171488 
            "KYAT3"        "RBMXL1"          "GBP1"          "GBP2"        "LRRC8C" 
    ENSG00000171492 ENSG00000122482 ENSG00000097046 ENSG00000069702 ENSG00000172031 
           "LRRC8D"        "ZNF644"          "CDC7"        "TGFBR3"         "EPHX4" 
    ENSG00000069712 ENSG00000067208 ENSG00000122406 ENSG00000143033 ENSG00000117500 
                 NA          "EVI5"          "RPL5"          "MTF2"         "TMED5" 
    ENSG00000122483 ENSG00000117505 ENSG00000137936 ENSG00000067334 ENSG00000023909 
           "CCDC18"           "DR1"         "BCAR3"       "DNTTIP2"          "GCLM" 
    ENSG00000137962 ENSG00000117525 ENSG00000117519 ENSG00000172339 ENSG00000117569 
         "ARHGAP29"            "F3"          "CNN3"         "ALG14"         "PTBP2" 
    ENSG00000188641 ENSG00000162627 ENSG00000117598 ENSG00000117600 ENSG00000156869 
             "DPYD"          "SNX7"        "PLPPR5"        "PLPPR4"         "FRRS1" 
    ENSG00000117620 ENSG00000156875 ENSG00000156876 ENSG00000122435 ENSG00000137992 
          "SLC35A3"       "MFSD14A"         "SASS6"        "TRMT13"           "DBT" 
    ENSG00000079335 ENSG00000162692 ENSG00000162694 ENSG00000170989 ENSG00000162631 
           "CDC14A"         "VCAM1"         "EXTL2"         "S1PR1"         "NTNG1" 
    ENSG00000085491 ENSG00000162636 ENSG00000121957 ENSG00000121940 ENSG00000085433 
         "SLC25A24"         "EEIG2"         "GPSM2"         "CLCC1"         "WDR47" 
    ENSG00000197780 ENSG00000215717 ENSG00000116299 ENSG00000031698 ENSG00000134222 
            "TAF13"      "TMEM167B"       "ELAPOR1"         "SARS1"         "PSRC1" 
    ENSG00000134243 ENSG00000143028 ENSG00000174151 ENSG00000181754 ENSG00000065135 
            "SORT1"         "SYPL2"      "CYB561D1"        "AMIGO1"         "GNAI3" 
    ENSG00000213366 ENSG00000184371 ENSG00000134248 ENSG00000177272 ENSG00000156171 
            "GSTM2"          "CSF1"       "LAMTOR5"         "KCNA3"         "DRAM2" 
    ENSG00000162777 ENSG00000173947 ENSG00000085465 ENSG00000116455 ENSG00000197852 
          "DENND2D"        "CIMAP3"         "OVGP1"         "WDR77"         "INKA2" 
    ENSG00000171385 ENSG00000143079 ENSG00000116489 ENSG00000155363 ENSG00000155366 
            "KCND3"     "CTTNBP2NL"        "CAPZA1"         "MOV10"          "RHOC" 
    ENSG00000155367 ENSG00000184599 ENSG00000155380 ENSG00000081026 ENSG00000116793 
            "PPM1J"         "TAFA3"       "SLC16A1"         "MAGI3"         "PHTF1" 
    ENSG00000081019 ENSG00000134262 ENSG00000118655 ENSG00000163349 ENSG00000197323 
            "RSBN1"         "AP4B1"       "DCLRE1B"         "HIPK1"        "TRIM33" 
    ENSG00000116752 ENSG00000175984 ENSG00000052723 ENSG00000134259 ENSG00000173218 
            "BCAS2"       "DENND2C"         "SIKE1"           "NGF"        "VANGL1" 
    ENSG00000163399 ENSG00000116815 ENSG00000134247 ENSG00000116830 ENSG00000134253 
           "ATP1A1"          "CD58"        "PTGFRN"          "TTF2"        "TRIM45" 
    ENSG00000198162 ENSG00000183508 ENSG00000196505 ENSG00000143067 ENSG00000092621 
           "MAN1A2"        "TENT5C"         "GDAP2"        "ZNF697"         "PHGDH" 
    ENSG00000134250 ENSG00000273136 ENSG00000188610 ENSG00000171943 ENSG00000263513 
           "NOTCH2"        "NBPF26"        "FAM72B"       "SRGAP2C"        "FAM72C" 
    ENSG00000266338 ENSG00000215784 ENSG00000174827 ENSG00000186141 ENSG00000131788 
           "NBPF15"        "FAM72D"         "PDZK1"        "POLR3C"         "PIAS3" 
    ENSG00000198483 ENSG00000143127 ENSG00000131779 ENSG00000265241 ENSG00000271601 
          "ANKRD35"        "ITGA10"        "PEX11B"         "RBM8A"         "LIX1L" 
    ENSG00000272031 ENSG00000265972 ENSG00000264343 ENSG00000268043 ENSG00000131791 
         "ANKRD34A"         "TXNIP"     "NOTCH2NLA"        "NBPF12"        "PRKAB2" 
    ENSG00000131778 ENSG00000270629 ENSG00000178104 ENSG00000150337 ENSG00000184678 
            "CHD1L"        "NBPF14"       "PDE4DIP"        "FCGR1A"        "H2BC21" 
    ENSG00000159164 ENSG00000023902 ENSG00000143401 ENSG00000117362 ENSG00000118292 
             "SV2A"       "PLEKHO1"        "ANP32E"         "APH1A"       "C1orf54" 
    ENSG00000266472 ENSG00000143369 ENSG00000143384 ENSG00000163131 ENSG00000143387 
           "MRPS21"          "ECM1"          "MCL1"          "CTSS"          "CTSK" 
    ENSG00000143409 ENSG00000143363 ENSG00000197622 ENSG00000213190 ENSG00000143434 
           "MINDY1"        "PRUNE1"      "CDC42SE1"        "MLLT11"        "SEMA6C" 
    ENSG00000163155 ENSG00000159352 ENSG00000143393 ENSG00000143390 ENSG00000143416 
           "LYSMD1"         "PSMD4"         "PI4KB"          "RFX5"      "SELENBP1" 
    ENSG00000143367 ENSG00000143376 ENSG00000143436 ENSG00000197747 ENSG00000163191 
            "TUFT1"         "SNX27"         "MRPL9"       "S100A10"       "S100A11" 
    ENSG00000197956 ENSG00000196154 ENSG00000196754 ENSG00000188643 ENSG00000189171 
           "S100A6"        "S100A4"        "S100A2"       "S100A16"       "S100A13" 
    ENSG00000143621 ENSG00000143624 ENSG00000143554 ENSG00000143614 ENSG00000160741 
             "ILF2"         "INTS3"       "SLC27A3"       "GATAD2B"         "CRTC2" 
    ENSG00000143570 ENSG00000143545 ENSG00000143549 ENSG00000143612 ENSG00000143569 
          "SLC39A1"         "RAB13"          "TPM3"       "C1orf43"        "UBAP2L" 
    ENSG00000143575 ENSG00000143515 ENSG00000160714 ENSG00000160710 ENSG00000163344 
             "HAX1"        "ATP8B2"        "UBE2Q1"          "ADAR"          "PMVK" 
    ENSG00000163346 ENSG00000163348 ENSG00000160691 ENSG00000173207 ENSG00000160688 
           "PBXIP1"         "PYGO2"          "SHC1"         "CKS1B"         "FLAD1" 
    ENSG00000160685 ENSG00000243364 ENSG00000169241 ENSG00000177628 ENSG00000160767 
           "ZBTB7B"         "EFNA4"       "SLC50A1"          "GBA1"       "ENTREP3" 
    ENSG00000116521 ENSG00000176444 ENSG00000143630 ENSG00000160752 ENSG00000160753 
           "SCAMP3"          "CLK2"          "HCN3"          "FDPS"         "RUSC1" 
    ENSG00000125459 ENSG00000163374 ENSG00000132676 ENSG00000143622 ENSG00000132680 
            "MSTO1"        "YY1AP1"          "DAP3"          "RIT1"         "KHDC4" 
    ENSG00000116584 ENSG00000116586 ENSG00000254726 ENSG00000160789 ENSG00000198952 
          "ARHGEF2"       "LAMTOR2"         "MEX3A"          "LMNA"          "SMG5" 
    ENSG00000198715 ENSG00000163468 ENSG00000116604 ENSG00000183856 ENSG00000160818 
             "GLMP"          "CCT3"         "MEF2D"        "IQGAP3"       "GPATCH4" 
    ENSG00000132688 ENSG00000143320 ENSG00000143319 ENSG00000143303 ENSG00000143314 
              "NES"        "CRABP2"       "ISG20L2"      "METTL25B"        "MRPL24" 
    ENSG00000187800 ENSG00000132694 ENSG00000117036 ENSG00000183853 ENSG00000158473 
            "PEAR1"      "ARHGEF11"          "ETV3"       "KIRREL1"          "CD1D" 
    ENSG00000163565 ENSG00000213085 ENSG00000158710 ENSG00000143315 ENSG00000162729 
            "IFI16"        "CFAP45"        "TAGLN2"          "PIGM"         "IGSF8" 
    ENSG00000162734 ENSG00000132716 ENSG00000122218 ENSG00000162736 ENSG00000162738 
            "PEA15"         "DCAF8"          "COPA"         "NCSTN"        "VANGL2" 
    ENSG00000158773 ENSG00000158796 ENSG00000143222 ENSG00000143258 ENSG00000158850 
             "USF1"          "DEDD"          "UFC1"         "USP21"       "B4GALT3" 
    ENSG00000158882 ENSG00000158887 ENSG00000143252 ENSG00000081721 ENSG00000118217 
          "TOMM40L"           "MPZ"          "SDHC"        "DUSP12"          "ATF6" 
    ENSG00000198929 ENSG00000239887 ENSG00000152332 ENSG00000117143 ENSG00000162733 
           "NOS1AP"      "C1orf226"         "UHMK1"          "UAP1"          "DDR2" 
    ENSG00000132196 ENSG00000117152 ENSG00000143248 ENSG00000143228 ENSG00000185630 
          "HSD17B7"          "RGS4"          "RGS5"          "NUF2"          "PBX1" 
    ENSG00000143149 ENSG00000143183 ENSG00000143157 ENSG00000152382 ENSG00000143190 
          "ALDH9A1"         "TMCO1"          "POGK"         "TADA1"        "POU2F1" 
    ENSG00000143162 ENSG00000197965 ENSG00000143164 ENSG00000143147 ENSG00000143155 
            "CREG1"         "MPZL1"         "DCAF6"        "GPR161"         "TIPRL" 
    ENSG00000213064 ENSG00000143153 ENSG00000117475 ENSG00000117479 ENSG00000000460 
           "SFT2D2"        "ATP1B1"         "BLZF1"       "SLC19A2"         "FIRRM" 
    ENSG00000000457 ENSG00000075945 ENSG00000116132 ENSG00000117523 ENSG00000117533 
            "SCYL3"        "KIFAP3"         "PRRX1"        "PRRC2C"         "VAMP4" 
    ENSG00000010165 ENSG00000197959 ENSG00000094975 ENSG00000120337 ENSG00000117586 
          "METTL13"          "DNM3"          "SUCO"       "TNFSF18"        "TNFSF4" 
    ENSG00000120334 ENSG00000117593 ENSG00000185278 ENSG00000135870 ENSG00000152061 
            "CENPL"         "DARS2"        "ZBTB37"         "RC3H1"      "RABGAP1L" 
    ENSG00000116161 ENSG00000120333 ENSG00000143207 ENSG00000116183 ENSG00000075391 
           "CACYBP"        "MRPS14"          "COP1"        "PAPPA2"        "RASAL2" 
    ENSG00000116191 ENSG00000116199 ENSG00000186283 ENSG00000169905 ENSG00000135837 
          "RALGPS2"        "FAM20B"         "TOR3A"      "TOR1AIP2"        "CEP350" 
    ENSG00000135823 ENSG00000153029 ENSG00000162783 ENSG00000143333 ENSG00000135829 
             "STX6"           "MR1"          "IER5"         "RGS16"          "DHX9" 
    ENSG00000135862 ENSG00000116698 ENSG00000116701 ENSG00000162704 ENSG00000143344 
            "LAMC1"          "SMG7"          "NCF2"         "ARPC5"          "RGL1" 
    ENSG00000198756 ENSG00000198860 ENSG00000116667 ENSG00000135842 ENSG00000116668 
         "COLGALT2"        "TSEN15"       "C1orf21"        "NIBAN1"          "SWT1" 
    ENSG00000116679 ENSG00000143341 ENSG00000116690 ENSG00000047410 ENSG00000157181 
         "IVNS1ABP"         "HMCN1"          "PRG4"           "TPR"          "ODR4" 
    ENSG00000073756 ENSG00000116711 ENSG00000162670 ENSG00000150681 ENSG00000116741 
            "PTGS2"       "PLA2G4A"        "BRINP3"         "RGS18"          "RGS2" 
    ENSG00000116750 ENSG00000134371 ENSG00000162630 ENSG00000162687 ENSG00000066279 
            "UCHL5"         "CDC73"       "B3GALT2"         "KCNT2"          "ASPM" 
    ENSG00000177888 ENSG00000213047 ENSG00000151414 ENSG00000116833 ENSG00000162702 
           "ZBTB41"       "DENND1B"          "NEK7"         "NR5A2"        "ZNF281" 
    ENSG00000118193 ENSG00000118197 ENSG00000118200 ENSG00000116852 ENSG00000174307 
            "KIF14"         "DDX59"       "CAMSAP2"        "KIF21B"        "PHLDA3" 
    ENSG00000159176 ENSG00000134369 ENSG00000163431 ENSG00000134375 ENSG00000143862 
            "CSRP1"          "NAV1"         "LMOD1"       "TIMM17A"         "ARL8A" 
    ENSG00000077152 ENSG00000077157 ENSG00000117139 ENSG00000183155 ENSG00000117153 
            "UBE2T"      "PPP1R12B"         "KDM5B"         "RABIF"        "KLHL12" 
    ENSG00000159346 ENSG00000159348 ENSG00000143847 ENSG00000163485 ENSG00000133048 
          "ADIPOR1"        "CYB5R1"        "PPFIA4"        "ADORA1"        "CHI3L1" 
    ENSG00000159388 ENSG00000122176 ENSG00000058668 ENSG00000058673 ENSG00000257315 
             "BTG2"          "FMOD"        "ATP2B4"       "ZC3H11A"         "ZBED6" 
    ENSG00000182004 ENSG00000143845 ENSG00000170498 ENSG00000143850 ENSG00000158615 
            "SNRPE"         "ETNK2"         "KISS1"       "PLEKHA6"      "PPP1R15B" 
    ENSG00000133056 ENSG00000163531 ENSG00000133069 ENSG00000163545 ENSG00000158715 
          "PIK3C2B"         "NFASC"         "TMCC2"         "NUAK2"       "SLC45A3" 
    ENSG00000069275 ENSG00000117280 ENSG00000276600 ENSG00000196550 ENSG00000266028 
           "NUCKS1"         "RAB29"         "RAB7B"        "FAM72A"        "SRGAP2" 
    ENSG00000263528 ENSG00000266094 ENSG00000143486 ENSG00000143479 ENSG00000162889 
            "IKBKE"        "RASSF5"         "EIF2D"         "DYRK3"      "MAPKAPK2" 
    ENSG00000162892 ENSG00000162894 ENSG00000123836 ENSG00000196352 ENSG00000117335 
             "IL24"          "FCMR"        "PFKFB2"          "CD55"          "CD46" 
    ENSG00000008118 ENSG00000196878 ENSG00000123689 ENSG00000117594 ENSG00000162757 
           "CAMK1G"         "LAMB3"          "G0S2"       "HSD11B1"       "C1orf74" 
    ENSG00000117597 ENSG00000143469 ENSG00000054392 ENSG00000143473 ENSG00000082512 
            "UTP25"         "SYT14"          "HHAT"         "KCNH1"         "TRAF5" 
    ENSG00000170385 ENSG00000117650 ENSG00000123684 ENSG00000143493 ENSG00000143476 
          "SLC30A1"          "NEK2"        "LPGAT1"         "INTS7"           "DTL" 
    ENSG00000065600 ENSG00000162772 ENSG00000123685 ENSG00000117697 ENSG00000203705 
            "PACC1"          "ATF3"         "BATF3"          "NSL1"        "TATDN3" 
    ENSG00000143494 ENSG00000174606 ENSG00000136643 ENSG00000117707 ENSG00000143499 
            "VASH2"        "ANGEL2"       "RPS6KC1"         "PROX1"         "SMYD2" 
    ENSG00000152104 ENSG00000117724 ENSG00000082482 ENSG00000136636 ENSG00000196482 
           "PTPN14"         "CENPF"         "KCNK2"         "KCTD3"         "ESRRG" 
    ENSG00000067533 ENSG00000092969 ENSG00000136628 ENSG00000118873 ENSG00000117791 
            "RRP15"         "TGFB2"         "EPRS1"      "RAB3GAP2"        "MTARC2" 
    ENSG00000143507 ENSG00000154305 ENSG00000186063 ENSG00000143502 ENSG00000162909 
           "DUSP10"          "MIA3"          "AIDA"         "SUSD4"         "CAPN2" 
    ENSG00000143756 ENSG00000143753 ENSG00000143748 ENSG00000162923 ENSG00000143786 
           "FBXO28"         "DEGS1"           "NVL"         "WDR26"         "CNIH3" 
    ENSG00000143815 ENSG00000154380 ENSG00000143819 ENSG00000143811 ENSG00000143751 
              "LBR"          "ENAH"         "EPHX1"         "PYCR2"          "SDE2" 
    ENSG00000163041 ENSG00000183814 ENSG00000143799 ENSG00000143776 ENSG00000181450 
            "H3-3A"          "LIN9"         "PARP1"      "CDC42BPA"        "ZNF678" 
    ENSG00000143740 ENSG00000143761 ENSG00000181873 ENSG00000162913 ENSG00000154358 
           "SNAP47"          "ARF1"         "IBA57"     "OBSCN-AS1"         "OBSCN" 
    ENSG00000154370 ENSG00000181218 ENSG00000168159 ENSG00000135776 ENSG00000135801 
           "TRIM11"        "H2AC25"        "RNF187"        "ABCB10"         "TAF5L" 
    ENSG00000143641 ENSG00000135744 ENSG00000119280 ENSG00000143643 ENSG00000182118 
           "GALNT2"           "AGT"      "C1orf198"         "TTC13"        "FAM89A" 
    ENSG00000116903 ENSG00000010072 ENSG00000116918 ENSG00000116991 ENSG00000135749 
            "EXOC8"         "SPRTN"         "TSNAX"       "SIPA1L2"         "PCNX2" 
    ENSG00000183780 ENSG00000059588 ENSG00000173726 ENSG00000116957 ENSG00000143669 
          "SLC35F3"        "TARBP1"        "TOMM20"              NA          "LYST" 
    ENSG00000116962 ENSG00000077585 ENSG00000086619 ENSG00000116977 ENSG00000119285 
             "NID1"       "GPR137B"         "ERO1B"        "LGALS8"        "HEATR1" 
    ENSG00000155816 ENSG00000091483 ENSG00000203668 ENSG00000174371 ENSG00000143702 
             "FMN2"            "FH"          "CHML"          "EXO1"        "CEP170" 
    ENSG00000054282 ENSG00000117020 ENSG00000179456 ENSG00000035687 ENSG00000121644 
          "SDCCAG8"          "AKT3"        "ZBTB18"         "ADSS2"         "DESI2" 
    ENSG00000185420 ENSG00000162851 ENSG00000162852 ENSG00000143653 ENSG00000197472 
            "SMYD3"         "TFB2M"          "CNST"        "SCCPDH"        "ZNF695" 
    ENSG00000188295 ENSG00000196418 ENSG00000162722 ENSG00000175137 ENSG00000185220 
           "ZNF669"        "ZNF124"        "TRIM58"       "SH3BP5L"         "PGBD2" 
    ENSG00000143727 ENSG00000151353 ENSG00000130508 ENSG00000182551 ENSG00000171863 
             "ACP1"        "TMEM18"          "PXDN"          "ADI1"          "RPS7" 
    ENSG00000176887 ENSG00000115738 ENSG00000134313 ENSG00000143797 ENSG00000151693 
            "SOX11"           "ID2"     "KIDINS220"        "MBOAT2"         "ASAP2" 
    ENSG00000119185 ENSG00000134330 ENSG00000134308 ENSG00000172059 ENSG00000171848 
         "ITGB1BP1"          "IAH1"         "YWHAQ"         "KLF11"          "RRM2" 
    ENSG00000115756 ENSG00000115758 ENSG00000115761 ENSG00000162975 ENSG00000134318 
           "HPCAL1"          "ODC1"         "NOL10"         "KCNF1"         "ROCK2" 
    ENSG00000169016 ENSG00000196208 ENSG00000071575 ENSG00000162981 ENSG00000151779 
             "E2F6"         "GREB1"         "TRIB2"        "LRATD1"          "NBAS" 
    ENSG00000079785 ENSG00000197872 ENSG00000143867 ENSG00000118965 ENSG00000068697 
             "DDX1"         "CYRIA"          "OSR1"         "WDR35"       "LAPTM4A" 
    ENSG00000115884 ENSG00000055917 ENSG00000143878 ENSG00000118960 ENSG00000119771 
             "SDC1"          "PUM2"          "RHOB"        "HS1BP3"        "KLHL29" 
    ENSG00000119778 ENSG00000173960 ENSG00000115128 ENSG00000115129 ENSG00000198399 
           "ATAD2B"        "UBXN2A"         "SF3B6"        "TP53I3"         "ITSN2" 
    ENSG00000138092 ENSG00000138031 ENSG00000143970 ENSG00000084731 ENSG00000157833 
            "CENPO"         "ADCY3"         "ASXL2"         "KIF3C"        "GAREM2" 
    ENSG00000138029 ENSG00000138018 ENSG00000213699 ENSG00000115163 ENSG00000119777 
            "HADHB"       "SELENOI"       "SLC35F6"         "CENPA"       "TMEM214" 
    ENSG00000084693 ENSG00000228474 ENSG00000138080 ENSG00000138073 ENSG00000138074 
            "AGBL5"          "OST4"       "EMILIN1"          "PREB"        "SLC5A6" 
    ENSG00000138085 ENSG00000084774 ENSG00000115207 ENSG00000163795 ENSG00000115241 
           "ATRAID"           "CAD"        "GTF3C2"        "ZNF513"         "PPM1G" 
    ENSG00000115216 ENSG00000138002 ENSG00000243943 ENSG00000176714 ENSG00000198522 
            "NRBP1"        "IFT172"        "ZNF512"       "CCDC121"          "GPN1" 
    ENSG00000119760 ENSG00000158019 ENSG00000075426 ENSG00000213639 ENSG00000163811 
           "SUPT7L"        "BABAM2"         "FOSL2"        "PPP1CB"         "WDR43" 
    ENSG00000119801 ENSG00000213626 ENSG00000172954 ENSG00000158125 ENSG00000152683 
            "YPEL5"           "LBH"        "LCLAT1"           "XDH"       "SLC30A6" 
    ENSG00000119820 ENSG00000115760 ENSG00000018699 ENSG00000049323 ENSG00000152689 
            "YIPF4"         "BIRC6"         "TTC27"         "LTBP1"       "RASGRP3" 
    ENSG00000119812 ENSG00000150938 ENSG00000171055 ENSG00000115808 ENSG00000008869 
           "FAM98A"         "CRIM1"          "FEZ2"          "STRN"       "HEATR5B" 
    ENSG00000055332 ENSG00000218739 ENSG00000115825 ENSG00000163171 ENSG00000115841 
          "EIF2AK2"       "CEBPZOS"         "PRKD3"      "CDC42EP3"         "RMDN2" 
    ENSG00000138061 ENSG00000119787 ENSG00000115875 ENSG00000163214 ENSG00000115904 
           "CYP1B1"          "ATL2"         "SRSF7"         "DHX57"          "SOS1" 
    ENSG00000011566 ENSG00000183023 ENSG00000162878 ENSG00000143924 ENSG00000115944 
           "MAP4K3"        "SLC8A1"         "PKDCC"          "EML4"       "COX7A2L" 
    ENSG00000057935 ENSG00000152518 ENSG00000152527 ENSG00000138095 ENSG00000138032 
             "MTA3"       "ZFP36L2"       "PLEKHH2"        "LRPPRC"         "PPM1B" 
    ENSG00000138078 ENSG00000143919 ENSG00000068784 ENSG00000171132 ENSG00000116016 
            "PREPL"        "CAMKMT"         "SRBD1"         "PRKCE"         "EPAS1" 
    ENSG00000250565 ENSG00000151665 ENSG00000171150 ENSG00000180398 ENSG00000068724 
         "ATP6V1E2"          "PIGF"         "SOCS5"         "MCFD2"         "TTC7A" 
    ENSG00000143933 ENSG00000095002 ENSG00000116062 ENSG00000138081 ENSG00000170802 
            "CALM2"          "MSH2"          "MSH6"        "FBXO11"         "FOXN2" 
    ENSG00000243244 ENSG00000143942 ENSG00000068912 ENSG00000115306 ENSG00000214595 
            "STON1"         "CHAC2"        "ERLEC1"        "SPTBN1"          "EML6" 
    ENSG00000115310 ENSG00000143947 ENSG00000115355 ENSG00000163001 ENSG00000275052 
             "RTN4"        "RPS27A"       "CCDC88A"        "CFAP36"       "PPP4R3B" 
    ENSG00000115380 ENSG00000055813 ENSG00000028116 ENSG00000115421 ENSG00000162929 
           "EFEMP1"       "CCDC85A"          "VRK2"        "PAPOLG"         "SANBR" 
    ENSG00000173209 ENSG00000115464 ENSG00000082898 ENSG00000115484 ENSG00000173163 
                 NA         "USP34"          "XPO1"          "CCT4"        "COMMD1" 
    ENSG00000170340 ENSG00000186889 ENSG00000014641 ENSG00000169764 ENSG00000143952 
           "B3GNT2"        "TMEM17"          "MDH1"          "UGP2"         "VPS54" 
    ENSG00000197329 ENSG00000179833 ENSG00000115902 ENSG00000011523 ENSG00000138069 
            "PELI1"       "SERTAD2"        "SLC1A4"         "CEP68"         "RAB1A" 
    ENSG00000138071 ENSG00000198369 ENSG00000115946 ENSG00000221823 ENSG00000119865 
            "ACTR2"        "SPRED2"          "PNO1"        "PPP3R1"        "CNRIP1" 
    ENSG00000169604 ENSG00000198380 ENSG00000115977 ENSG00000196975 ENSG00000087338 
           "ANTXR1"         "GFPT1"          "AAK1"         "ANXA4"         "GMCL1" 
    ENSG00000169564 ENSG00000116001 ENSG00000116005 ENSG00000035141 ENSG00000124357 
            "PCBP1"          "TIA1"        "PCYOX1"       "FAM136A"          "NAGK" 
    ENSG00000124374 ENSG00000075292 ENSG00000135636 ENSG00000003137 ENSG00000144040 
           "PAIP2B"        "ZNF638"          "DYSF"       "CYP26B1"         "SFXN5" 
    ENSG00000135631 ENSG00000135632 ENSG00000135617 ENSG00000135624 ENSG00000163013 
        "RAB11FIP5"         "SMYD5"        "PRADC1"          "CCT7"        "FBXO41" 
    ENSG00000124356 ENSG00000163017 ENSG00000187605 ENSG00000163170 ENSG00000114978 
           "STAMBP"         "ACTG2"          "TET3"         "BOLA3"         "MOB1A" 
    ENSG00000065911 ENSG00000204843 ENSG00000239779 ENSG00000115289 ENSG00000115307 
           "MTHFD2"         "DCTN1"          "WBP1"         "PCGF1"          "AUP1" 
    ENSG00000115325 ENSG00000135622 ENSG00000159399 ENSG00000115353 ENSG00000115363 
             "DOK1"        "SEMA4F"           "HK2"         "TACR1"         "EVA1A" 
    ENSG00000115364 ENSG00000005436 ENSG00000163541 ENSG00000186854 ENSG00000034510 
           "MRPL19"         "GCFC2"        "SUCLG1"       "TRABD2A"        "TMSB10" 
    ENSG00000176407 ENSG00000152291 ENSG00000042445 ENSG00000042493 ENSG00000168906 
            "KCMF1"        "TGOLN2"        "RETSAT"          "CAPG"         "MAT2A" 
    ENSG00000115486 ENSG00000168899 ENSG00000168890 ENSG00000168883 ENSG00000168874 
             "GGCX"         "VAMP5"      "TMEM150A"         "USP39"         "ATOH8" 
    ENSG00000132305 ENSG00000115548 ENSG00000239305 ENSG00000153561 ENSG00000172071 
             "IMMT"         "KDM3A"        "RNF103"        "RMND5A"       "EIF2AK3" 
    ENSG00000153574 ENSG00000233757 ENSG00000115042 ENSG00000174501 ENSG00000084090 
             "RPIA"        "ZNF892"        "FAHD2A"      "ANKRD36C"        "STARD7" 
    ENSG00000135956 ENSG00000144021 ENSG00000198885 ENSG00000121152 ENSG00000114982 
          "TMEM127"         "CIAO1"      "ITPRIPL1"         "NCAPH"        "KANSL3" 
    ENSG00000114988 ENSG00000158158 ENSG00000168763 ENSG00000213337 ENSG00000168758 
           "LMAN2L"         "CNNM4"         "CNNM3"       "ANKRD39"        "SEMA4C" 
    ENSG00000135976 ENSG00000075568 ENSG00000183513 ENSG00000115446 ENSG00000071073 
          "ANKRD36"       "TMEM131"          "COA5"         "UNC50"        "MGAT4A" 
    ENSG00000185414 ENSG00000144218 ENSG00000170485 ENSG00000071082 ENSG00000204634 
           "MRPL30"          "AFF3"         "NPAS2"         "RPL31"        "TBC1D8" 
    ENSG00000158435 ENSG00000175874 ENSG00000071054 ENSG00000115594 ENSG00000115598 
           "CNOT11"         "CREG2"        "MAP4K4"         "IL1R1"        "IL1RL2" 
    ENSG00000115602 ENSG00000115604 ENSG00000135953 ENSG00000135972 ENSG00000115641 
           "IL1RL1"        "IL18R1"         "MFSD9"         "MRPS9"          "FHL2" 
    ENSG00000071051 ENSG00000144057 ENSG00000135968 ENSG00000169756 ENSG00000153201 
             "NCK2"       "ST6GAL2"          "GCC2"         "LIMS1"        "RANBP2" 
    ENSG00000163006 ENSG00000172985 ENSG00000186522 ENSG00000198142 ENSG00000144063 
          "CCDC138"        "SH3RF3"      "SEPTIN10"        "SOWAHC"          "MALL" 
    ENSG00000169679 ENSG00000153094 ENSG00000153107 ENSG00000153214 ENSG00000125630 
             "BUB1"       "BCL2L11"        "ANAPC1"       "TMEM87B"        "POLR1B" 
    ENSG00000169607 ENSG00000115008 ENSG00000125538 ENSG00000144134 ENSG00000115084 
           "CKAP2L"          "IL1A"          "IL1B"        "RABL2A"       "SLC35F5" 
    ENSG00000115091 ENSG00000115107 ENSG00000155368 ENSG00000144120 ENSG00000088179 
            "ACTR3"        "STEAP3"           "DBI"       "TMEM177"         "PTPN4" 
    ENSG00000115109 ENSG00000144118 ENSG00000074047 ENSG00000074054 ENSG00000211460 
          "EPB41L5"          "RALB"          "GLI2"        "CLASP1"           "TSN" 
    ENSG00000163161 ENSG00000169967 ENSG00000072163 ENSG00000144231 ENSG00000144233 
            "ERCC3"        "MAP3K2"         "LIMS2"        "POLR2D"      "AMMECR1L" 
    ENSG00000136731 ENSG00000136720 ENSG00000196604 ENSG00000136699 ENSG00000152082 
            "UGGT1"        "HS6ST1"         "POTEF"         "SMPD4"         "MZT2B" 
    ENSG00000136710 ENSG00000136718 ENSG00000072135 ENSG00000152102 ENSG00000115762 
          "CCDC115"          "IMP4"        "PTPN18"       "FAM168B"       "PLEKHB2" 
    ENSG00000173272 ENSG00000183840 ENSG00000150551 ENSG00000176771 ENSG00000152127 
            "MZT2A"         "GPR39"         "LYPD1"        "NCKAP5"         "MGAT5" 
    ENSG00000082258 ENSG00000115839 ENSG00000121988 ENSG00000048991 ENSG00000144224 
            "CCNT2"      "RAB3GAP1"        "ZRANB3"        "R3HDM1"         "UBXN4" 
    ENSG00000076003 ENSG00000115866 ENSG00000150540 ENSG00000144228 ENSG00000168702 
             "MCM6"         "DARS1"          "HNMT"         "SPOPL"         "LRP1B" 
    ENSG00000115919 ENSG00000121964 ENSG00000169554 ENSG00000115947 ENSG00000204406 
             "KYNU"         "GTDC1"          "ZEB2"          "ORC4"          "MBD5" 
    ENSG00000135999 ENSG00000187123 ENSG00000168288 ENSG00000115963 ENSG00000123609 
             "EPC2"         "LYPD6"        "MMADHC"          "RND3"           "NMI" 
    ENSG00000123610 ENSG00000080345 ENSG00000162980 ENSG00000115145 ENSG00000157827 
          "TNFAIP6"          "RIF1"         "ARL5A"         "STAM2"         "FMNL2" 
    ENSG00000153234 ENSG00000115159 ENSG00000115170 ENSG00000144283 ENSG00000115183 
            "NR4A2"          "GPD2"         "ACVR1"          "PKP4"         "TANC1" 
    ENSG00000123636 ENSG00000153250 ENSG00000136560 ENSG00000115233 ENSG00000197635 
            "BAZ2B"         "RBMS1"          "TANK"        "PSMD14"          "DPP4" 
    ENSG00000078098 ENSG00000182263 ENSG00000082438 ENSG00000136531 ENSG00000115339 
              "FAP"          "FIGN"        "COBLL1"         "SCN2A"        "GALNT3" 
    ENSG00000169432 ENSG00000136546 ENSG00000172318 ENSG00000163072 ENSG00000152253 
            "SCN9A"         "SCN7A"       "B3GALT1"       "NOSTRIN"         "SPC25" 
    ENSG00000138382 ENSG00000144357 ENSG00000115806 ENSG00000123600 ENSG00000115827 
           "METTL5"          "UBR3"       "GORASP2"        "METTL8"        "DCAF17" 
    ENSG00000071967 ENSG00000077380 ENSG00000115840 ENSG00000128708 ENSG00000172878 
           "CYBRD1"       "DYNC1I2"      "SLC25A12"          "HAT1"       "METAP1D" 
    ENSG00000115844 ENSG00000091409 ENSG00000152256 ENSG00000091428 ENSG00000091436 
             "DLX2"         "ITGA6"          "PDK1"       "RAPGEF4"       "MAP3K20" 
    ENSG00000144354 ENSG00000172845 ENSG00000138430 ENSG00000138433 ENSG00000144306 
            "CDCA7"           "SP3"          "OLA1"          "CIR1"         "SCRN3" 
    ENSG00000163328 ENSG00000115935 ENSG00000128656 ENSG00000115966 ENSG00000154518 
           "GPR155"         "WIPF1"          "CHN1"          "ATF2"       "ATP5MC3" 
    ENSG00000144320 ENSG00000170144 ENSG00000116044 ENSG00000018510 ENSG00000128655 
             "LNPK"       "HNRNPA3"        "NFE2L2"          "AGPS"        "PDE11A" 
    ENSG00000079156 ENSG00000079150 ENSG00000116095 ENSG00000187231 ENSG00000170035 
           "OSBPL6"         "FKBP7"       "PLEKHA3"        "SESTD1"        "UBE2E3" 
    ENSG00000115232 ENSG00000188452 ENSG00000115252 ENSG00000077232 ENSG00000162998 
            "ITGA4"         "CERKL"         "PDE1A"       "DNAJC10"          "FRZB" 
    ENSG00000061676 ENSG00000163002 ENSG00000170396 ENSG00000138448 ENSG00000144369 
           "NCKAP1"         "NUP35"       "ZNF804A"         "ITGAV"       "FAM171B" 
    ENSG00000064989 ENSG00000003436 ENSG00000144366 ENSG00000168542 ENSG00000204262 
           "CALCRL"          "TFPI"         "GULP1"        "COL3A1"        "COL5A2" 
    ENSG00000138449 ENSG00000138381 ENSG00000128699 ENSG00000151689 ENSG00000189362 
          "SLC40A1"        "ASNSD1"        "ORMDL1"         "INPP1"         "NEMP2" 
    ENSG00000115419 ENSG00000115415 ENSG00000128641 ENSG00000168497 ENSG00000144339 
              "GLS"         "STAT1"         "MYO1B"        "CAVIN2"        "TMEFF2" 
    ENSG00000196950 ENSG00000081320 ENSG00000138411 ENSG00000144395 ENSG00000197121 
         "SLC39A10"        "STK17B"         "HECW2"       "CCDC150"         "PGAP1" 
    ENSG00000065413 ENSG00000115520 ENSG00000144381 ENSG00000115541 ENSG00000115540 
          "ANKRD44"        "COQ10B"         "HSPD1"         "HSPE1"          "MOB4" 
    ENSG00000162944 ENSG00000247626 ENSG00000119042 ENSG00000178074 ENSG00000196141 
            "RFTN2"         "MARS2"         "SATB2"       "C2orf69"       "SPATS2L" 
    ENSG00000163535 ENSG00000138356 ENSG00000082153 ENSG00000013441 ENSG00000240344 
             "SGO2"          "AOX1"          "BZW1"          "CLK1"         "PPIL3" 
    ENSG00000155744 ENSG00000003402 ENSG00000064012 ENSG00000082146 ENSG00000155755 
            "HYCC2"         "CFLAR"         "CASP8"        "STRADB"       "TMEM237" 
    ENSG00000082126 ENSG00000138395 ENSG00000155760 ENSG00000182329 ENSG00000116030 
             "MPP4"         "CDK15"          "FZD7"      "KIAA2012"         "SUMO1" 
    ENSG00000055044 ENSG00000138439 ENSG00000138442 ENSG00000144426 ENSG00000119004 
            "NOP58"       "FAM117B"         "WDR12"        "NBEAL1"       "CYP20A1" 
    ENSG00000138443 ENSG00000173166 ENSG00000116117 ENSG00000114933 ENSG00000114942 
             "ABI2"         "RAPH1"        "PARD3B"        "INO80D"        "EEF1B2" 
    ENSG00000183671 ENSG00000204186 ENSG00000114948 ENSG00000118246 ENSG00000118263 
           "CMKLR2"         "ZDBF2"        "ADAM23"       "FASTKD2"          "KLF7" 
    ENSG00000163249 ENSG00000138413 ENSG00000115020 ENSG00000197713 ENSG00000144445 
           "CCNYL1"          "IDH1"       "PIKFYVE"           "RPE"       "KANSL1L" 
    ENSG00000115365 ENSG00000030419 ENSG00000144451 ENSG00000138376 ENSG00000138363 
           "LANCL1"         "IKZF2"        "SPAG16"         "BARD1"          "ATIC" 
    ENSG00000115414 ENSG00000118242 ENSG00000115425 ENSG00000144583 ENSG00000138375 
              "FN1"          "MREG"          "PECR"       "MARCHF4"      "SMARCAL1" 
    ENSG00000197756 ENSG00000115457 ENSG00000115461 ENSG00000231672 ENSG00000079308 
           "RPL37A"        "IGFBP2"        "IGFBP5"         "DIRC3"          "TNS1" 
    ENSG00000163466 ENSG00000127837 ENSG00000127838 ENSG00000135926 ENSG00000144579 
            "ARPC2"          "AAMP"          "PNKD"        "TMBIM1"        "CTDSP1" 
    ENSG00000135913 ENSG00000144580 ENSG00000074582 ENSG00000163482 ENSG00000135912 
            "USP37"         "CNOT9"         "BCS1L"         "STK36"         "TTLL4" 
    ENSG00000135929 ENSG00000144567 ENSG00000198925 ENSG00000163516 ENSG00000127824 
          "CYP27A1"       "RETREG2"         "ATG9A"        "ANKZF1"        "TUBA4A" 
    ENSG00000135924 ENSG00000123992 ENSG00000175084 ENSG00000072195 ENSG00000144591 
           "DNAJB2"         "DNPEP"           "DES"          "SPEG"         "GMPPA" 
    ENSG00000123989 ENSG00000188760 ENSG00000124006 ENSG00000144589 ENSG00000116106 
             "CHPF"       "TMEM198"         "OBSL1"       "STK11IP"         "EPHA4" 
    ENSG00000116120 ENSG00000123983 ENSG00000152049 ENSG00000171951 ENSG00000152056 
            "FARSB"         "ACSL3"         "KCNE4"          "SCG2"         "AP1S3" 
    ENSG00000085449 ENSG00000135919 ENSG00000135905 ENSG00000144468 ENSG00000168958 
            "WDFY1"      "SERPINE2"        "DOCK10"        "RHBDD1"           "MFF" 
    ENSG00000123977 ENSG00000187957 ENSG00000163053 ENSG00000185404 ENSG00000067066 
             "DAW1"          "DNER"      "SLC16A14"        "SP140L"         "SP100" 
    ENSG00000135932 ENSG00000135916 ENSG00000173692 ENSG00000135931 ENSG00000115053 
            "CAB39"         "ITM2C"         "PSMD1"         "ARMC9"           "NCL" 
    ENSG00000187514 ENSG00000156973 ENSG00000163273 ENSG00000135930 ENSG00000204120 
             "PTMA"         "PDE6D"          "NPPC"        "EIF4E2"        "GIGYF2" 
    ENSG00000168918 ENSG00000077044 ENSG00000085982 ENSG00000123485 ENSG00000188042 
           "INPP5D"          "DGKD"         "USP40"         "HJURP"         "ARL4C" 
    ENSG00000157985 ENSG00000198612 ENSG00000115648 ENSG00000124831 ENSG00000132329 
            "AGAP1"         "COPS8"          "MLPH"       "LRRFIP1"         "RAMP1" 
    ENSG00000184182 ENSG00000144488 ENSG00000178752 ENSG00000144485 ENSG00000132326 
            "UBE2F"         "ESPNL"          "ERFE"          "HES6"          "PER2" 
    ENSG00000204104 ENSG00000065802 ENSG00000130414 ENSG00000063660 ENSG00000142330 
         "TRAF3IP1"          "ASB1"       "NDUFA10"          "GPC1"        "CAPN10" 
    ENSG00000162804 ENSG00000115687 ENSG00000146205 ENSG00000115677 ENSG00000168385 
            "SNED1"          "PASK"          "ANO7"         "HDLBP"       "SEPTIN2" 
    ENSG00000006607 ENSG00000115694 ENSG00000176720 ENSG00000176946 ENSG00000168397 
            "FARP2"         "STK25"           "BOK"         "THAP4"         "ATG4B" 
    ENSG00000168393 ENSG00000180902 ENSG00000134121 ENSG00000144619 ENSG00000113851 
            "DTYMK"        "D2HGDH"          "CHL1"         "CNTN4"          "CRBN" 
    ENSG00000144455 ENSG00000175928 ENSG00000134108 ENSG00000134109 ENSG00000071282 
            "SUMF1"         "LRRN1"         "ARL8B"         "EDEM1"         "LMCD1" 
    ENSG00000180914 ENSG00000070950 ENSG00000196220 ENSG00000168137 ENSG00000163719 
             "OXTR"         "RAD18"        "SRGAP3"         "SETD5"        "MTMR14" 
    ENSG00000156983 ENSG00000114026 ENSG00000134072 ENSG00000171148 ENSG00000214021 
            "BRPF1"          "OGG1"         "CAMK1"         "TADA3"         "TTLL3" 
    ENSG00000156990 ENSG00000163701 ENSG00000163703 ENSG00000163704 ENSG00000125037 
           "RPUSD3"        "IL17RE"        "CRELD1"         "PRRT3"          "EMC3" 
    ENSG00000144554 ENSG00000254999 ENSG00000134086 ENSG00000134070 ENSG00000157014 
           "FANCD2"          "BRK1"           "VHL"         "IRAK2"        "TATDN2" 
    ENSG00000157020 ENSG00000196639 ENSG00000197548 ENSG00000132170 ENSG00000075975 
            "SEC13"          "HRH1"          "ATG7"         "PPARG"         "MKRN2" 
    ENSG00000132155 ENSG00000144712 ENSG00000144713 ENSG00000144711 ENSG00000170876 
             "RAF1"         "CAND2"         "RPL32"        "IQSEC1"        "TMEM43" 
    ENSG00000170860 ENSG00000177463 ENSG00000131368 ENSG00000144597 ENSG00000131373 
             "LSM3"         "NR2C2"        "MRPS25"          "EAF1"         "HACL1" 
    ENSG00000206560 ENSG00000131386 ENSG00000154813 ENSG00000154814 ENSG00000131378 
          "ANKRD28"       "GALNT15"          "DPH3"        "OXNAD1"         "RFTN1" 
    ENSG00000182568 ENSG00000144566 ENSG00000114166 ENSG00000129810 ENSG00000151789 
            "SATB1"         "RAB5A"         "KAT2B"          "SGO1"       "ZNF385D" 
    ENSG00000182247 ENSG00000170142 ENSG00000174738 ENSG00000163491 ENSG00000033867 
           "UBE2E2"        "UBE2E1"         "NR1D2"         "NEK10"        "SLC4A7" 
    ENSG00000144642 ENSG00000163513 ENSG00000163527 ENSG00000152642 ENSG00000153551 
            "RBMS3"        "TGFBR2"         "STT3B"         "GPD1L"         "CMTM7" 
    ENSG00000091317 ENSG00000144635 ENSG00000182973 ENSG00000170266 ENSG00000188167 
            "CMTM6"      "DYNC1LI1"        "CNOT10"          "GLB1"         "TMPPE" 
    ENSG00000170275 ENSG00000173705 ENSG00000153558 ENSG00000153560 ENSG00000163539 
            "CRTAP"         "SUSD5"         "FBXL2"          "UBP1"        "CLASP2" 
    ENSG00000178567 ENSG00000076242 ENSG00000093167 ENSG00000144674 ENSG00000172936 
         "EPM2AIP1"          "MLH1"       "LRRFIP2"        "GOLGA4"         "MYD88" 
    ENSG00000172939 ENSG00000093217 ENSG00000114739 ENSG00000157036 ENSG00000114742 
            "OXSR1"          "XYLB"        "ACVR2B"          "EXOG"         "WDR48" 
    ENSG00000114745 ENSG00000168026 ENSG00000144655 ENSG00000168334 ENSG00000144659 
          "GORASP1"        "TTC21A"        "CSRNP1"         "XIRP1"      "SLC25A38" 
    ENSG00000168028 ENSG00000114784 ENSG00000188846 ENSG00000168036 ENSG00000182606 
             "RPSA"         "EIF1B"         "RPL14"        "CTNNB1"         "TRAK1" 
    ENSG00000093183 ENSG00000114857 ENSG00000181061 ENSG00000144647 ENSG00000163788 
           "SEC22C"          "NKTR"        "HIGD1A"       "POMGNT2"          "SNRK" 
    ENSG00000011198 ENSG00000179152 ENSG00000163808 ENSG00000169964 ENSG00000163812 
            "ABHD5"         "TCAIM"         "KIF15"        "TMEM42"        "ZDHHC3" 
    ENSG00000075914 ENSG00000163814 ENSG00000249992 ENSG00000144791 ENSG00000163827 
           "EXOSC7"         "CDCP1"       "TMEM158"         "LIMD1"         "LRRC2" 
    ENSG00000178038 ENSG00000160796 ENSG00000181555 ENSG00000114648 ENSG00000076201 
           "ALS2CL"        "NBEAL2"         "SETD2"        "KLHL18"        "PTPN23" 
    ENSG00000114650 ENSG00000163832 ENSG00000173473 ENSG00000132153 ENSG00000047849 
             "SCAP"          "ELP6"       "SMARCC1"         "DHX30"          "MAP4" 
    ENSG00000164045 ENSG00000172113 ENSG00000164050 ENSG00000114270 ENSG00000225697 
           "CDC25A"          "NME6"        "PLXNB1"        "COL7A1"       "SLC26A6" 
    ENSG00000008300 ENSG00000213672 ENSG00000068745 ENSG00000178467 ENSG00000178252 
           "CELSR3"       "NCKIPSD"         "IP6K2"         "P4HTM"          "WDR6" 
    ENSG00000178057 ENSG00000178035 ENSG00000172053 ENSG00000172046 ENSG00000177352 
          "NDUFAF3"        "IMPDH2"         "QARS1"         "USP19"        "CCDC71" 
    ENSG00000185909 ENSG00000188315 ENSG00000114316 ENSG00000067560 ENSG00000173402 
          "KLHDC8B"       "C3orf62"          "USP4"          "RHOA"          "DAG1" 
    ENSG00000173531 ENSG00000164068 ENSG00000176095 ENSG00000182179 ENSG00000183763 
             "MST1"        "RNF123"         "IP6K1"          "UBA7"         "TRAIP" 
    ENSG00000004534 ENSG00000003756 ENSG00000001617 ENSG00000114353 ENSG00000012171 
             "RBM6"          "RBM5"        "SEMA3F"         "GNAI2"        "SEMA3B" 
    ENSG00000214706 ENSG00000243477 ENSG00000068001 ENSG00000068028 ENSG00000126062 
            "IFRD2"         "NAA80"         "HYAL2"        "RASSF1"       "TMEM115" 
    ENSG00000114737 ENSG00000114738 ENSG00000145050 ENSG00000259956 ENSG00000164080 
             "CISH"      "MAPKAPK3"          "MANF"        "RBM15B"       "RAD54L2" 
    ENSG00000164081 ENSG00000114767 ENSG00000041880 ENSG00000090097 ENSG00000114779 
           "TEX264"          "RRP9"         "PARP3"         "PCBP4"       "ABHD14B" 
    ENSG00000248487 ENSG00000162244 ENSG00000164086 ENSG00000164087 ENSG00000023330 
          "ABHD14A"         "RPL29"         "DUSP7"         "POC1A"         "ALAS1" 
    ENSG00000164088 ENSG00000114841 ENSG00000163930 ENSG00000010322 ENSG00000168268 
            "PPM1M"         "DNAH1"          "BAP1"         "NISCH"        "NT5DC2" 
    ENSG00000163939 ENSG00000114904 ENSG00000163935 ENSG00000163933 ENSG00000163932 
            "PBRM1"          "NEK4"        "SFMBT1"          "RFT1"         "PRKCD" 
    ENSG00000163931 ENSG00000272886 ENSG00000113812 ENSG00000113811 ENSG00000157445 
              "TKT"         "DCP1A"         "ACTR8"       "SELENOK"      "CACNA2D3" 
    ENSG00000114251 ENSG00000163946 ENSG00000163947 ENSG00000144730 ENSG00000157500 
            "WNT5A"         "TASOR"       "ARHGEF3"        "IL17RD"         "APPL1" 
    ENSG00000174840 ENSG00000174839 ENSG00000163681 ENSG00000136068 ENSG00000168297 
            "PDE12"       "DENND6A"         "SLMAP"          "FLNB"           "PXK" 
    ENSG00000168291 ENSG00000168301 ENSG00000168306 ENSG00000163689 ENSG00000144724 
             "PDHB"         "KCTD6"         "ACOX2"      "CFAP20DC"         "PTPRG" 
    ENSG00000163634 ENSG00000163635 ENSG00000163636 ENSG00000163637 ENSG00000163638 
            "THOC7"         "ATXN7"         "PSMD6"      "PRICKLE2"       "ADAMTS9" 
    ENSG00000151276 ENSG00000144741 ENSG00000144749 ENSG00000163376 ENSG00000172340 
            "MAGI1"      "SLC25A26"         "LRIG1"        "KBTBD8"        "SUCLG2" 
    ENSG00000163378 ENSG00000144747 ENSG00000144746 ENSG00000114541 ENSG00000187098 
             "EOGT"          "TMF1"       "ARL6IP5"        "FRMD4B"          "MITF" 
    ENSG00000163602 ENSG00000144736 ENSG00000172986 ENSG00000121440 ENSG00000185008 
             "RYBP"          "SHQ1"        "GXYLT2"        "PDZRN3"         "ROBO2" 
    ENSG00000169855 ENSG00000114480 ENSG00000206538 ENSG00000083937 ENSG00000163320 
            "ROBO1"          "GBE1"         "VGLL3"        "CHMP2B"        "CGGBP1" 
    ENSG00000175105 ENSG00000179021 ENSG00000044524 ENSG00000184500 ENSG00000169379 
           "ZNF654"       "C3orf38"         "EPHA3"         "PROS1"        "ARL13B" 
    ENSG00000178700 ENSG00000178694 ENSG00000170854 ENSG00000080822 ENSG00000080819 
            "DHFR2"         "NSUN3"         "RIOX2"        "CLDND1"          "CPOX" 
    ENSG00000057019 ENSG00000144810 ENSG00000184220 ENSG00000168386 ENSG00000036054 
           "DCBLD2"        "COL8A1"         "CMSS1"       "FILIP1L"       "TBC1D23" 
    ENSG00000154174 ENSG00000206535 ENSG00000181458 ENSG00000114354 ENSG00000154175 
           "TOMM70"          "LNP1"       "TMEM45A"           "TFG"        "ABI3BP" 
    ENSG00000174173 ENSG00000114391 ENSG00000182504 ENSG00000144815 ENSG00000144802 
          "TRMT10C"         "RPL24"         "CEP97"         "NXPE3"        "NFKBIZ" 
    ENSG00000170017 ENSG00000114423 ENSG00000114439 ENSG00000196776 ENSG00000144821 
            "ALCAM"          "CBLB"           "BBX"          "CD47"         "MYH15" 
    ENSG00000163507 ENSG00000198919 ENSG00000177707 ENSG00000240891 ENSG00000144824 
            "CIP2A"         "DZIP3"       "NECTIN3"        "PLCXD2"        "PHLDB2" 
    ENSG00000144827 ENSG00000144834 ENSG00000091986 ENSG00000163607 ENSG00000144857 
           "ABHD10"        "TAGLN3"        "CCDC80"        "GTPBP8"           "BOC" 
    ENSG00000206530 ENSG00000163611 ENSG00000114573 ENSG00000178075 ENSG00000172020 
           "CFAP44"        "SPICE1"       "ATP6V1A"       "GRAMD1C"         "GAP43" 
    ENSG00000185565 ENSG00000121578 ENSG00000113845 ENSG00000082701 ENSG00000163428 
            "LSAMP"       "B4GALT4"       "TIMMDC1"         "GSK3B"        "LRRC58" 
    ENSG00000065518 ENSG00000144840 ENSG00000153767 ENSG00000051341 ENSG00000160124 
           "NDUFB4"         "RABL3"        "GTF2E1"          "POLQ"         "MIX23" 
    ENSG00000114030 ENSG00000138496 ENSG00000163840 ENSG00000173193 ENSG00000169087 
            "KPNA1"         "PARP9"         "DTX3L"        "PARP14"       "HSPBAP1" 
    ENSG00000138463 ENSG00000121542 ENSG00000206527 ENSG00000065534 ENSG00000175455 
          "SLC49A4"        "SEC22A"         "HACD2"          "MYLK"        "CCDC14" 
    ENSG00000160145 ENSG00000082781 ENSG00000173706 ENSG00000221955 ENSG00000163848 
            "KALRN"         "ITGB5"          "HEG1"       "SLC12A8"        "ZNF148" 
    ENSG00000114544 ENSG00000070476 ENSG00000159685 ENSG00000114554 ENSG00000163870 
          "SLC41A3"          "ZXDC"        "CHCHD6"        "PLXNA1"         "TPRA1" 
    ENSG00000073111 ENSG00000114626 ENSG00000074416 ENSG00000058262 ENSG00000175792 
             "MCM2"         "ABTB1"          "MGLL"       "SEC61A1"        "RUVBL1" 
    ENSG00000132394 ENSG00000179348 ENSG00000163902 ENSG00000075785 ENSG00000177646 
           "EEFSEC"         "GATA2"          "RPN1"         "RAB7A"         "ACAD9" 
    ENSG00000114654 ENSG00000172780 ENSG00000181789 ENSG00000129071 ENSG00000172765 
            "EFCC1"         "RAB43"         "COPG1"          "MBD4"         "TMCC1" 
    ENSG00000196455 ENSG00000017260 ENSG00000034533 ENSG00000114686 ENSG00000138246 
           "PIK3R4"        "ATP2C1"         "ASTE1"         "MRPL3"       "DNAJC13" 
    ENSG00000129048 ENSG00000081307 ENSG00000091527 ENSG00000163781 ENSG00000154917 
            "ACKR4"          "UBA5"          "CDV3"        "TOPBP1"         "RAB6B" 
    ENSG00000174640 ENSG00000114019 ENSG00000154928 ENSG00000073711 ENSG00000174579 
          "SLCO2A1"        "AMOTL2"         "EPHB1"       "PPP2R3A"          "MSL2" 
    ENSG00000114054 ENSG00000158092 ENSG00000158163 ENSG00000114098 ENSG00000158186 
             "PCCB"          "NCK1"        "DZIP1L"         "ARMC8"          "MRAS" 
    ENSG00000114107 ENSG00000051382 ENSG00000183770 ENSG00000206262 ENSG00000114115 
            "CEP70"        "PIK3CB"         "FOXL2"       "FOXL2NB"          "RBP1" 
    ENSG00000114120 ENSG00000155893 ENSG00000177311 ENSG00000155903 ENSG00000069849 
         "SLC25A36"        "PXYLP1"        "ZBTB38"         "RASA2"        "ATP1B3" 
    ENSG00000114126 ENSG00000114127 ENSG00000163710 ENSG00000163714 ENSG00000175040 
            "TFDP2"          "XRN1"       "PCOLCE2"        "U2SURP"         "CHST2" 
    ENSG00000181744 ENSG00000152952 ENSG00000114698 ENSG00000188313 ENSG00000163754 
           "DIPK2A"         "PLOD2"        "PLSCR4"        "PLSCR1"          "GYG1" 
    ENSG00000071794 ENSG00000163755 ENSG00000169908 ENSG00000114744 ENSG00000082996 
             "HLTF"          "HPS3"        "TM4SF1"        "COMMD2"         "RNF13" 
    ENSG00000196428 ENSG00000144895 ENSG00000198843 ENSG00000144893 ENSG00000152580 
          "TSC22D2"         "EIF2A"       "SELENOT"        "MED12L"        "IGSF10" 
    ENSG00000114771 ENSG00000152601 ENSG00000169860 ENSG00000114790 ENSG00000174953 
            "AADAC"         "MBNL1"         "P2RY1"      "ARHGEF26"         "DHX36" 
    ENSG00000196549 ENSG00000174928 ENSG00000163655 ENSG00000169282 ENSG00000114850 
              "MME"       "C3orf33"          "GMPS"        "KCNAB1"          "SSR3" 
    ENSG00000163659 ENSG00000163660 ENSG00000197415 ENSG00000163661 ENSG00000079257 
           "TIPARP"         "CCNL1"         "VEPH1"          "PTX3"           "LXN" 
    ENSG00000118855 ENSG00000250588 ENSG00000168811 ENSG00000180044 ENSG00000113810 
            "MFSD1"              NA         "IL12A"       "C3orf80"          "SMC4" 
    ENSG00000213186 ENSG00000186432 ENSG00000163590 ENSG00000169255 ENSG00000169251 
           "TRIM59"         "KPNA4"         "PPM1L"      "B3GALNT1"          "NMD3" 
    ENSG00000114209 ENSG00000173905 ENSG00000085276 ENSG00000008952 ENSG00000163558 
           "PDCD10"        "GOLIM4"         "MECOM"         "SEC62"         "PRKCI" 
    ENSG00000136603 ENSG00000013297 ENSG00000013293 ENSG00000163584 ENSG00000163577 
             "SKIL"        "CLDN11"       "SLC7A14"       "RPL22L1"        "EIF5A2" 
    ENSG00000154310 ENSG00000075651 ENSG00000075420 ENSG00000144959 ENSG00000114346 
             "TNIK"          "PLD1"        "FNDC3B"         "NCEH1"          "ECT2" 
    ENSG00000177694 ENSG00000177565 ENSG00000172667 ENSG00000121879 ENSG00000121864 
         "NAALADL2"       "TBL1XR1"         "ZMAT3"        "PIK3CA"        "ZNF639" 
    ENSG00000171109 ENSG00000136518 ENSG00000058056 ENSG00000163728 ENSG00000205981 
             "MFN1"        "ACTL6A"         "USP13"         "TTC14"       "DNAJC19" 
    ENSG00000043093 ENSG00000114796 ENSG00000180834 ENSG00000114770 ENSG00000161202 
          "DCUN1D1"        "KLHL24"        "MAP6D1"         "ABCC5"          "DVL3" 
    ENSG00000161203 ENSG00000214160 ENSG00000145194 ENSG00000175166 ENSG00000114867 
            "AP2M1"          "ALG3"          "ECE2"         "PSMD2"        "EIF4G1" 
    ENSG00000114859 ENSG00000163882 ENSG00000182580 ENSG00000177383 ENSG00000156931 
            "CLCN2"        "POLR2H"         "EPHB3"        "MAGEF1"          "VPS8" 
    ENSG00000187068 ENSG00000113790 ENSG00000163898 ENSG00000163904 ENSG00000073792 
          "C3orf70"        "EHHADH"          "LIPH"         "SENP2"       "IGF2BP2" 
    ENSG00000136527 ENSG00000244405 ENSG00000058866 ENSG00000113838 ENSG00000090520 
            "TRA2B"          "ETV5"          "DGKG"        "TBCCD1"       "DNAJB11" 
    ENSG00000156976 ENSG00000163918 ENSG00000073849 ENSG00000163923 ENSG00000127241 
           "EIF4A2"          "RFC4"       "ST6GAL1"        "RPL39L"         "MASP1" 
    ENSG00000145012 ENSG00000090530 ENSG00000163347 ENSG00000196083 ENSG00000152492 
              "LPP"          "P3H2"         "CLDN1"        "IL1RAP"        "CCDC50" 
    ENSG00000180611 ENSG00000198836 ENSG00000114315 ENSG00000172061 ENSG00000133657 
           "MB21D2"          "OPA1"          "HES1"        "LRRC15"       "ATP13A3" 
    ENSG00000145014 ENSG00000185112 ENSG00000114331 ENSG00000184203 ENSG00000061938 
           "TMEM44"        "FAM43A"         "ACAP2"        "PPP1R2"          "TNK2" 
    ENSG00000161217 ENSG00000163961 ENSG00000174004 ENSG00000163964 ENSG00000174007 
           "PCYT1A"        "RNF168"         "NRROS"          "PIGX"         "CEP19" 
    ENSG00000180370 ENSG00000119231 ENSG00000114503 ENSG00000270170 ENSG00000075711 
             "PAK2"         "SENP5"         "NCBP2"      "NCBP2AS2"          "DLG1" 
    ENSG00000122068 ENSG00000182903 ENSG00000174227 ENSG00000185619 ENSG00000178950 
           "FYTTD1"        "ZNF721"          "PIGG"         "PCGF3"           "GAK" 
    ENSG00000145214 ENSG00000127415 ENSG00000127418 ENSG00000159674 ENSG00000159692 
             "DGKQ"          "IDUA"        "FGFRL1"         "SPON2"         "CTBP1" 
    ENSG00000179979 ENSG00000163950 ENSG00000168936 ENSG00000013810 ENSG00000109685 
                 NA          "SLBP"       "TMEM129"         "TACC3"          "NSD2" 
    ENSG00000185049 ENSG00000243449 ENSG00000185818 ENSG00000123933 ENSG00000063978 
            "NELFA"        "NICOL1"         "NAT8L"          "MXD4"          "RNF4" 
    ENSG00000168884 ENSG00000197386 ENSG00000159788 ENSG00000163956 ENSG00000145220 
            "TNIP2"           "HTT"         "RGS12"        "LRPAP1"          "LYAR" 
    ENSG00000168824 ENSG00000168818 ENSG00000170891 ENSG00000152953 ENSG00000072840 
             "NSG1"         "STX18"         "CYTL1"        "STK32B"           "EVC" 
    ENSG00000013288 ENSG00000179010 ENSG00000178988 ENSG00000186222 ENSG00000170871 
           "MAN2B2"        "MRFAP1"      "MRFAP1L1"       "BLOC1S4"      "KIAA0232" 
    ENSG00000196526 ENSG00000125089 ENSG00000087008 ENSG00000071127 ENSG00000178163 
            "AFAP1"        "SH3TC1"         "ACOX3"          "WDR1"       "ZNF518B" 
    ENSG00000038219 ENSG00000137449 ENSG00000118564 ENSG00000237765 ENSG00000169762 
           "BOD1L1"         "CPEB2"         "FBXL5"       "FAM200B"         "TAPT1" 
    ENSG00000169744 ENSG00000151552 ENSG00000002549 ENSG00000118579 ENSG00000163257 
             "LDB2"          "QDPR"          "LAP3"         "MED28"        "DCAF16" 
    ENSG00000109805 ENSG00000178177 ENSG00000145147 ENSG00000152990 ENSG00000109606 
            "NCAPG"         "LCORL"         "SLIT2"        "ADGRA3"         "DHX15" 
    ENSG00000181982 ENSG00000038210 ENSG00000168228 ENSG00000091490 ENSG00000168214 
          "CCDC149"        "PI4K2B"        "ZCCHC4"        "SEL1L3"          "RBPJ" 
    ENSG00000109680 ENSG00000109689 ENSG00000169851 ENSG00000181826 ENSG00000169299 
          "TBC1D19"         "STIM2"         "PCDH7"         "RELL1"          "PGM2" 
    ENSG00000065882 ENSG00000109787 ENSG00000174130 ENSG00000197712 ENSG00000121895 
           "TBC1D1"          "KLF3"          "TLR6"      "FAM114A1"       "TMEM156" 
    ENSG00000109790 ENSG00000157796 ENSG00000035928 ENSG00000163682 ENSG00000163683 
            "KLHL5"         "WDR19"          "RFC1"          "RPL9"        "SMIM14" 
    ENSG00000078140 ENSG00000078177 ENSG00000163697 ENSG00000154277 ENSG00000064042 
            "UBE2K"         "N4BP2"         "APBB2"         "UCHL1"        "LIMCH1" 
    ENSG00000109133 ENSG00000182308 ENSG00000124406 ENSG00000183783 ENSG00000163281 
           "TMEM33"       "DCAF4L1"        "ATP8A1"         "KCTD8"        "GNPDA2" 
    ENSG00000151834 ENSG00000109158 ENSG00000163288 ENSG00000145246 ENSG00000145244 
           "GABRA2"        "GABRA4"        "GABRB1"        "ATP10D"         "CORIN" 
    ENSG00000170448 ENSG00000163293 ENSG00000109180 ENSG00000188993 ENSG00000163071 
            "NFXL1"        "NIPAL1"        "OCIAD1"        "LRRC66"       "SPATA18" 
    ENSG00000109189 ENSG00000226887 ENSG00000128045 ENSG00000145216 ENSG00000072201 
            "USP46"    "ERVMER34-1"       "RASL11B"        "FIP1L1"          "LNX1" 
    ENSG00000109220 ENSG00000134853 ENSG00000157404 ENSG00000128052 ENSG00000134851 
            "CHIC2"        "PDGFRA"           "KIT"           "KDR"       "TMEM165" 
    ENSG00000134852 ENSG00000174799 ENSG00000109265 ENSG00000128059 ENSG00000128050 
            "CLOCK"        "CEP135"         "CRACD"          "PPAT"         "PAICS" 
    ENSG00000174780 ENSG00000171476 ENSG00000084093 ENSG00000084092 ENSG00000047315 
            "SRP72"          "HOPX"          "REST"          "NOA1"        "POLR2B" 
    ENSG00000163453 ENSG00000150471 ENSG00000145242 ENSG00000132463 ENSG00000173542 
           "IGFBP7"        "ADGRL3"         "EPHA5"         "GRSF1"         "MOB1B" 
    ENSG00000080493 ENSG00000132466 ENSG00000081051 ENSG00000169429 ENSG00000124875 
           "SLC4A4"       "ANKRD17"           "AFP"         "CXCL8"         "CXCL6" 
    ENSG00000163739 ENSG00000163735 ENSG00000163734 ENSG00000081041 ENSG00000124882 
            "CXCL1"         "CXCL5"         "CXCL3"         "CXCL2"          "EREG" 
    ENSG00000109321 ENSG00000169116 ENSG00000138757 ENSG00000138768 ENSG00000138744 
             "AREG"         "PARM1"         "G3BP2"          "USO1"          "NAAA" 
    ENSG00000198301 ENSG00000138750 ENSG00000138760 ENSG00000138771 ENSG00000138764 
            "SDAD1"         "NUP54"        "SCARB2"       "SHROOM3"         "CCNG2" 
    ENSG00000138767 ENSG00000169288 ENSG00000138772 ENSG00000163291 ENSG00000152784 
           "CNOT6L"         "MRPL1"         "ANXA3"         "PAQR3"         "PRDM8" 
    ENSG00000138675 ENSG00000138669 ENSG00000138668 ENSG00000145293 ENSG00000145284 
             "FGF5"         "PRKG2"        "HNRNPD"        "ENOPH1"          "SCD5" 
    ENSG00000138674 ENSG00000138663 ENSG00000138678 ENSG00000163625 ENSG00000109339 
           "SEC31A"         "COPS4"         "GPAT3"         "WDFY3"        "MAPK10" 
    ENSG00000163629 ENSG00000172493 ENSG00000198189 ENSG00000170502 ENSG00000118785 
           "PTPN13"          "AFF1"      "HSD17B11"         "NUDT9"          "SPP1" 
    ENSG00000118762 ENSG00000118777 ENSG00000138642 ENSG00000138641 ENSG00000177432 
             "PKD2"         "ABCG2"         "HERC6"         "HERC3"        "NAP1L5" 
    ENSG00000138640 ENSG00000180346 ENSG00000185477 ENSG00000145335 ENSG00000163104 
           "FAM13A"         "TIGD2"        "GPRIN3"          "SNCA"      "SMARCAD1" 
    ENSG00000163110 ENSG00000138696 ENSG00000182168 ENSG00000138698 ENSG00000151247 
           "PDLIM5"        "BMPR1B"         "UNC5C"      "RAP1GDS1"         "EIF4E" 
    ENSG00000164024 ENSG00000197894 ENSG00000109270 ENSG00000164031 ENSG00000164032 
           "METAP1"          "ADH5"       "LAMTOR3"       "DNAJB14"         "H2AZ1" 
    ENSG00000145358 ENSG00000138814 ENSG00000109320 ENSG00000109323 ENSG00000109332 
           "DDIT4L"        "PPP3CA"         "NFKB1"         "MANBA"        "UBE2D3" 
    ENSG00000145354 ENSG00000164038 ENSG00000164039 ENSG00000138778 ENSG00000168769 
            "CISD2"        "SLC9B2"          "BDH2"         "CENPE"          "TET2" 
    ENSG00000138777 ENSG00000138780 ENSG00000145348 ENSG00000164022 ENSG00000155011 
             "PPA2"         "GSTCD"          "TBCK"         "AIMP1"          "DKK2" 
    ENSG00000138801 ENSG00000164023 ENSG00000155016 ENSG00000138796 ENSG00000138795 
           "PAPSS1"         "SGMS2"        "CYP2U1"          "HADH"          "LEF1" 
    ENSG00000198856 ENSG00000138802 ENSG00000005059 ENSG00000138798 ENSG00000170522 
             "OSTC"        "SEC24B"          "MCUB"           "EGF"        "ELOVL6" 
    ENSG00000164093 ENSG00000174749 ENSG00000138660 ENSG00000145365 ENSG00000073331 
            "PITX2"       "FAM241A"         "AP1AR"          "TIFA"         "ALPK1" 
    ENSG00000138658 ENSG00000145362 ENSG00000145349 ENSG00000180801 ENSG00000164099 
            "ZGRF1"          "ANK2"        "CAMK2D"          "ARSJ"        "PRSS12" 
    ENSG00000150961 ENSG00000172403 ENSG00000145390 ENSG00000164096 ENSG00000164109 
           "SEC24D"        "SYNPO2"         "USP53"        "C4orf3"        "MAD2L1" 
    ENSG00000138738 ENSG00000164111 ENSG00000123737 ENSG00000145386 ENSG00000138686 
            "PRDM5"         "ANXA5"        "EXOSC9"         "CCNA2"          "BBS7" 
    ENSG00000138741 ENSG00000138688 ENSG00000181004 ENSG00000138685 ENSG00000164056 
            "TRPC3"         "BLTP1"         "BBS12"          "FGF2"         "SPRY1" 
    ENSG00000151458 ENSG00000196159 ENSG00000164066 ENSG00000164070 ENSG00000142731 
          "ANKRD50"          "FAT4"          "INTU"        "HSPA4L"          "PLK4" 
    ENSG00000164073 ENSG00000164040 ENSG00000077684 ENSG00000151466 ENSG00000151470 
            "MFSD8"        "PGRMC2"         "JADE1"         "SCLT1"       "C4orf33" 
    ENSG00000138650 ENSG00000254535 ENSG00000189184 ENSG00000151012 ENSG00000151014 
           "PCDH10"       "PABPC4L"        "PCDH18"       "SLC7A11"          "NOCT" 
    ENSG00000137463 ENSG00000109390 ENSG00000164134 ENSG00000145391 ENSG00000196782 
            "MGARP"        "NDUFC1"         "NAA15"         "SETD7"         "MAML3" 
    ENSG00000153132 ENSG00000179387 ENSG00000170153 ENSG00000109445 ENSG00000164136 
             "CLGN"        "ELMOD2"        "RNF150"        "ZNF330"          "IL15" 
    ENSG00000109452 ENSG00000170185 ENSG00000109458 ENSG00000153147 ENSG00000164161 
           "INPP4B"         "USP38"          "GAB1"       "SMARCA5"          "HHIP" 
    ENSG00000164162 ENSG00000164163 ENSG00000170365 ENSG00000151612 ENSG00000137473 
          "ANAPC10"         "ABCE1"         "SMAD1"        "ZNF827"         "TTC29" 
    ENSG00000151617 ENSG00000071205 ENSG00000151623 ENSG00000170390 ENSG00000198589 
            "EDNRA"      "ARHGAP10"         "NR3C2"         "DCLK2"          "LRBA" 
    ENSG00000181541 ENSG00000145425 ENSG00000109686 ENSG00000059691 ENSG00000109670 
          "MAB21L2"         "RPS3A"        "SH3D19"          "GATB"         "FBXW7" 
    ENSG00000170006 ENSG00000164144 ENSG00000109654 ENSG00000121211 ENSG00000171566 
          "TMEM154"        "ARFIP1"         "TRIM2"          "MND1"         "PLRG1" 
    ENSG00000061918 ENSG00000256043 ENSG00000145431 ENSG00000205208 ENSG00000171503 
          "GUCY1B1"          "CTSO"         "PDGFC"       "C4orf46"         "ETFDH" 
    ENSG00000171497 ENSG00000052795 ENSG00000170088 ENSG00000052802 ENSG00000109472 
             "PPID"         "FNIP2"       "TMEM192"         "MSMO1"           "CPE" 
    ENSG00000196104 ENSG00000109511 ENSG00000137628 ENSG00000154447 ENSG00000137601 
           "SPOCK3"        "ANXA10"         "DDX60"        "SH3RF1"          "NEK1" 
    ENSG00000109572 ENSG00000056050 ENSG00000109576 ENSG00000109586 ENSG00000164104 
            "CLCN3"          "HPF1"         "AADAT"        "GALNT7"         "HMGB2" 
    ENSG00000150628 ENSG00000164122 ENSG00000150630 ENSG00000109674 ENSG00000218336 
           "SPATA4"          "ASB5"         "VEGFC"         "NEIL3"         "TENM3" 
    ENSG00000151718 ENSG00000168556 ENSG00000182552 ENSG00000173320 ENSG00000164305 
             "WWC2"          "ING2"         "RWDD4"         "STOX2"         "CASP3" 
    ENSG00000164306 ENSG00000151725 ENSG00000151726 ENSG00000151729 ENSG00000164323 
          "PRIMPOL"         "CENPU"         "ACSL1"       "SLC25A4"        "CFAP97" 
    ENSG00000109775 ENSG00000168491 ENSG00000154553 ENSG00000164342 ENSG00000145476 
            "UFSP2"       "CCDC110"        "PDLIM3"          "TLR3"        "CYP4V2" 
    ENSG00000083857 ENSG00000073578 ENSG00000249915 ENSG00000063438 ENSG00000066230 
             "FAT1"          "SDHA"         "PDCD6"          "AHRR"        "SLC9A3" 
    ENSG00000028310 ENSG00000071539 ENSG00000145506 ENSG00000113504 ENSG00000049656 
             "BRD9"        "TRIP13"          "NKD2"       "SLC12A7"       "CLPTM1L" 
    ENSG00000153395 ENSG00000171421 ENSG00000186493 ENSG00000164151 ENSG00000133398 
           "LPCAT1"        "MRPL36"              NA          "ICE1"         "MED10" 
    ENSG00000037474 ENSG00000145545 ENSG00000112941 ENSG00000124275 ENSG00000112902 
            "NSUN2"        "SRD5A1"        "TENT4A"          "MTRR"        "SEMA5A" 
    ENSG00000150753 ENSG00000164237 ENSG00000112977 ENSG00000039139 ENSG00000038382 
             "CCT5"          "CMBL"           "DAP"         "DNAH5"          "TRIO" 
    ENSG00000145569 ENSG00000154122 ENSG00000173545 ENSG00000145555 ENSG00000113361 
          "OTULINL"          "ANKH"        "ZNF622"         "MYO10"          "CDH6" 
    ENSG00000113360 ENSG00000133401 ENSG00000113384 ENSG00000056097 ENSG00000113387 
           "DROSHA"         "PDZD2"        "GOLPH3"           "ZFR"          "SUB1" 
    ENSG00000113389 ENSG00000113407 ENSG00000151388 ENSG00000039560 ENSG00000168724 
             "NPR3"         "TARS1"      "ADAMTS12"         "RAI14"       "DNAJC21" 
    ENSG00000113494 ENSG00000152582 ENSG00000168685 ENSG00000164187 ENSG00000145604 
             "PRLR"         "SPEF2"          "IL7R"        "LMBRD2"          "SKP2" 
    ENSG00000164190 ENSG00000197603 ENSG00000113569 ENSG00000168621 ENSG00000113594 
            "NIPBL"       "CPLANE1"        "NUP155"          "GDNF"          "LIFR" 
    ENSG00000145623 ENSG00000164327 ENSG00000153071 ENSG00000171522 ENSG00000132356 
             "OSMR"        "RICTOR"          "DAB2"        "PTGER4"        "PRKAA1" 
    ENSG00000145592 ENSG00000112936 ENSG00000083720 ENSG00000205765 ENSG00000112964 
            "RPL37"            "C7"         "OXCT1"        "RIMOC1"           "GHR" 
    ENSG00000177453 ENSG00000112972 ENSG00000151881 ENSG00000172244 ENSG00000172239 
            "NIM1K"        "HMGCS1"       "TMEM267"       "C5orf34"         "PAIP1" 
    ENSG00000112992 ENSG00000170571 ENSG00000213949 ENSG00000164171 ENSG00000134363 
              "NNT"           "EMB"         "ITGA1"         "ITGA2"           "FST" 
    ENSG00000185305 ENSG00000178996 ENSG00000164283 ENSG00000164294 ENSG00000067248 
            "ARL15"         "SNX18"          "ESM1"          "GPX8"         "DHX29" 
    ENSG00000067113 ENSG00000177058 ENSG00000164509 ENSG00000134352 ENSG00000164512 
            "PLPP1"       "SLC38A9"        "IL31RA"         "IL6ST"       "ANKRD55" 
    ENSG00000155545 ENSG00000062194 ENSG00000169067 ENSG00000145632 ENSG00000035499 
            "MIER3"         "GPBP1"        "ACTBL2"          "PLK2"       "DEPDC1B" 
    ENSG00000049167 ENSG00000164182 ENSG00000188725 ENSG00000130449 ENSG00000068796 
            "ERCC8"       "NDUFAF2"        "SMIM15"        "ZSWIM6"         "KIF2A" 
    ENSG00000086189 ENSG00000086200 ENSG00000153006 ENSG00000049192 ENSG00000123219 
            "DIMT1"         "IPO11"      "SREK1IP1"       "ADAMTS6"         "CENPK" 
    ENSG00000113595 ENSG00000197860 ENSG00000123213 ENSG00000069020 ENSG00000145675 
           "TRIM23"          "SGTB"           "NLN"         "MAST4"        "PIK3R1" 
    ENSG00000145740 ENSG00000134057 ENSG00000153044 ENSG00000134056 ENSG00000273841 
          "SLC30A5"         "CCNB1"         "CENPH"          "KGD4"          "TAF9" 
    ENSG00000152939 ENSG00000205572 ENSG00000172062 ENSG00000145734 ENSG00000131711 
         "MARVELD2"        "SERF1A"          "SMN1"          "BDP1"         "MAP1B" 
    ENSG00000113048 ENSG00000049883 ENSG00000083312 ENSG00000157107 ENSG00000157111 
           "MRPS27"         "PTCD2"         "TNPO1"         "FCHO2"       "TMEM171" 
    ENSG00000164331 ENSG00000214944 ENSG00000171617 ENSG00000049860 ENSG00000164347 
           "ANKRA2"      "ARHGEF28"          "ENC1"          "HEXB"          "GFM2" 
    ENSG00000113161 ENSG00000113163 ENSG00000145703 ENSG00000164220 ENSG00000181104 
            "HMGCR"         "CERT1"        "IQGAP2"         "F2RL2"           "F2R" 
    ENSG00000164251 ENSG00000132846 ENSG00000171530 ENSG00000132842 ENSG00000085365 
            "F2RL1"         "ZBED3"          "TBCA"         "AP3B1"        "SCAMP1" 
    ENSG00000145685 ENSG00000152409 ENSG00000152413 ENSG00000164329 ENSG00000177034 
           "LHFPL2"           "JMY"        "HOMER1"         "TENT2"          "MTX3" 
    ENSG00000164300 ENSG00000039319 ENSG00000228716 ENSG00000113318 ENSG00000113319 
          "SERINC5"       "ZFYVE16"          "DHFR"          "MSH3"       "RASGRF2" 
    ENSG00000145687 ENSG00000152348 ENSG00000174695 ENSG00000038427 ENSG00000164176 
            "SSBP2"         "ATG10"      "TMEM167A"          "VCAN"         "EDIL3" 
    ENSG00000145715 ENSG00000164180 ENSG00000081189 ENSG00000153140 ENSG00000176055 
            "RASA1"      "TMEM161B"         "MEF2C"         "CETN3"        "MBLAC2" 
    ENSG00000113356 ENSG00000176018 ENSG00000113369 ENSG00000175745 ENSG00000133302 
           "POLR3G"        "LYSMD3"        "ARRDC3"         "NR2F1"          "SLF1" 
    ENSG00000175471 ENSG00000164291 ENSG00000164292 ENSG00000173221 ENSG00000118985 
            "MCTP1"          "ARSK"       "RHOBTB3"          "GLRX"          "ELL2" 
    ENSG00000153113 ENSG00000113441 ENSG00000058729 ENSG00000174136 ENSG00000153922 
             "CAST"         "LNPEP"         "RIOK2"          "RGMB"          "CHD1" 
    ENSG00000174132 ENSG00000145730 ENSG00000181751 ENSG00000112874 ENSG00000184349 
          "FAM174A"           "PAM"         "MACIR"        "NUDT12"         "EFNA5" 
    ENSG00000145743 ENSG00000198961 ENSG00000112893 ENSG00000164209 ENSG00000145777 
           "FBXL17"          "PJA2"        "MAN2A1"      "SLC25A46"          "TSLP" 
    ENSG00000134987 ENSG00000152495 ENSG00000164211 ENSG00000134986 ENSG00000134982 
            "WDR36"         "CAMK4"        "STARD4"          "NREP"           "APC" 
    ENSG00000129625 ENSG00000172795 ENSG00000171444 ENSG00000164219 ENSG00000145780 
            "REEP5"          "DCP2"           "MCC"        "PGGT1B"         "FEM1C" 
    ENSG00000134970 ENSG00000145782 ENSG00000177879 ENSG00000145781 ENSG00000092421 
            "TMED7"         "ATG12"         "AP3S1"       "COMMD10"        "SEMA6A" 
    ENSG00000145779 ENSG00000133835 ENSG00000184838 ENSG00000113083 ENSG00000205302 
          "TNFAIP8"       "HSD17B4"         "PRR16"           "LOX"          "SNX2" 
    ENSG00000064652 ENSG00000061455 ENSG00000168944 ENSG00000151292 ENSG00000168916 
            "SNX24"         "PRDM6"        "CEP120"       "CSNK1G3"        "ZNF608" 
    ENSG00000155324 ENSG00000164904 ENSG00000164902 ENSG00000113368 ENSG00000173926 
          "GRAMD2B"       "ALDH7A1"          "PHAX"         "LMNB1"       "MARCHF3" 
    ENSG00000164244 ENSG00000064651 ENSG00000138829 ENSG00000113396 ENSG00000145808 
            "PRRC1"       "SLC12A2"          "FBN2"       "SLC27A6"      "ADAMTS19" 
    ENSG00000198108 ENSG00000169567 ENSG00000186687 ENSG00000158985 ENSG00000217128 
            "CHSY3"         "HINT1"         "LYRM7"      "CDC42SE2"         "FNIP1" 
    ENSG00000131435 ENSG00000197208 ENSG00000197375 ENSG00000197536 ENSG00000125347 
           "PDLIM4"       "SLC22A4"       "SLC22A5"              NA          "IRF1" 
    ENSG00000113522 ENSG00000131437 ENSG00000164402 ENSG00000155329 ENSG00000170606 
            "RAD50"         "KIF3A"       "SEPTIN8"       "ZCCHC10"         "HSPA4" 
    ENSG00000113583 ENSG00000213585 ENSG00000081059 ENSG00000113558 ENSG00000113575 
          "C5orf15"         "VDAC1"          "TCF7"          "SKP1"        "PPP2CA" 
    ENSG00000119048 ENSG00000237190 ENSG00000164615 ENSG00000145833 ENSG00000181904 
            "UBE2B"    "CDKN2AIPNL"         "CAMLG"         "DDX46"       "C5orf24" 
    ENSG00000113621 ENSG00000069011 ENSG00000113648 ENSG00000120708 ENSG00000152377 
          "TXNDC15"         "PITX1"     "MACROH2A1"         "TGFBI"        "SPOCK1" 
    ENSG00000177733 ENSG00000031003 ENSG00000112984 ENSG00000094880 ENSG00000158402 
          "HNRNPA0"        "FAM13B"        "KIF20A"         "CDC23"        "CDC25C" 
    ENSG00000120709 ENSG00000132563 ENSG00000120738 ENSG00000120705 ENSG00000113013 
           "FAM53C"         "REEP2"          "EGR1"          "ETF1"         "HSPA9" 
    ENSG00000044115 ENSG00000170464 ENSG00000184584 ENSG00000171604 ENSG00000185129 
           "CTNNA1"       "DNAJC18"        "STING1"         "CXXC5"          "PURA" 
    ENSG00000182700 ENSG00000120306 ENSG00000113068 ENSG00000113070 ENSG00000176087 
             "IGIP"        "CYSTM1"         "PFDN1"         "HBEGF"       "SLC35A4" 
    ENSG00000113119 ENSG00000113141 ENSG00000170445 ENSG00000146007 ENSG00000204967 
            "TMCO6"            "IK"         "HARS1"         "ZMAT2"        "PCDHA4" 
    ENSG00000204963 ENSG00000112852 ENSG00000113205 ENSG00000113209 ENSG00000113212 
           "PCDHA7"        "PCDHB2"        "PCDHB3"        "PCDHB5"        "PCDHB7" 
    ENSG00000120322 ENSG00000272674 ENSG00000177839 ENSG00000120324 ENSG00000197479 
           "PCDHB8"       "PCDHB16"        "PCDHB9"       "PCDHB10"       "PCDHB11" 
    ENSG00000120328 ENSG00000187372 ENSG00000120327 ENSG00000113248 ENSG00000178913 
          "PCDHB12"       "PCDHB13"       "PCDHB14"       "PCDHB15"          "TAF7" 
    ENSG00000204956 ENSG00000081853 ENSG00000254245 ENSG00000254221 ENSG00000262576 
          "PCDHGA1"       "PCDHGA2"       "PCDHGA3"       "PCDHGB1"       "PCDHGA4" 
    ENSG00000253910 ENSG00000253485 ENSG00000262209 ENSG00000253731 ENSG00000253537 
          "PCDHGB2"       "PCDHGA5"       "PCDHGB3"       "PCDHGA6"       "PCDHGA7" 
    ENSG00000253953 ENSG00000276547 ENSG00000253305 ENSG00000253846 ENSG00000254122 
          "PCDHGB4"       "PCDHGB5"       "PCDHGB6"      "PCDHGA10"       "PCDHGB7" 
    ENSG00000240184 ENSG00000240764 ENSG00000131504 ENSG00000197948 ENSG00000013561 
          "PCDHGC3"       "PCDHGC5"        "DIAPH1"        "FCHSD1"         "RNF14" 
    ENSG00000113552 ENSG00000131507 ENSG00000187678 ENSG00000113578 ENSG00000113580 
           "GNPDA1"        "NDFIP1"         "SPRY4"          "FGF1"         "NR3C1" 
    ENSG00000183775 ENSG00000186314 ENSG00000156463 ENSG00000133706 ENSG00000091009 
           "KCTD16"       "PRELID2"        "SH3RF2"         "LARS1"         "RBM27" 
    ENSG00000113649 ENSG00000113657 ENSG00000145868 ENSG00000169252 ENSG00000169247 
           "TCERG1"        "DPYSL3"        "FBXO38"         "ADRB2"        "SH3TC2" 
    ENSG00000173210 ENSG00000157510 ENSG00000164284 ENSG00000145882 ENSG00000113712 
           "ABLIM3"       "AFAP1L1"        "GRPEL2"       "PCYOX1L"       "CSNK1A1" 
    ENSG00000155846 ENSG00000155850 ENSG00000113716 ENSG00000113721 ENSG00000183876 
         "PPARGC1B"       "SLC26A2"        "HMGXB3"        "PDGFRB"          "ARSI" 
    ENSG00000070814 ENSG00000164587 ENSG00000070614 ENSG00000171992 ENSG00000132912 
            "TCOF1"         "RPS14"         "NDST1"         "SYNPO"         "DCTN4" 
    ENSG00000256235 ENSG00000211445 ENSG00000145901 ENSG00000197043 ENSG00000198624 
            "SMIM3"          "GPX3"         "TNIP1"         "ANXA6"        "CCDC69" 
    ENSG00000123643 ENSG00000177556 ENSG00000145907 ENSG00000155511 ENSG00000055147 
          "SLC36A1"         "ATOX1"         "G3BP1"         "GRIA1"      "FAM114A2" 
    ENSG00000037749 ENSG00000164576 ENSG00000170271 ENSG00000155508 ENSG00000082516 
            "MFAP3"        "SAP30L"        "FAXDC2"         "CNOT8"        "GEMIN5" 
    ENSG00000170624 ENSG00000055163 ENSG00000135074 ENSG00000172548 ENSG00000155858 
             "SGCD"        "CYFIP2"        "ADAM19"        "NIPAL4"         "LSM11" 
    ENSG00000113282 ENSG00000164330 ENSG00000145860 ENSG00000164332 ENSG00000170214 
           "CLINT1"          "EBF1"        "RNF145"        "UBLCP1"        "ADRA1B" 
    ENSG00000113312 ENSG00000170234 ENSG00000135083 ENSG00000145861 ENSG00000221886 
             "TTC1"        "PWWP2A"         "CCNJL"       "C1QTNF2"       "FAM200C" 
    ENSG00000164609 ENSG00000164611 ENSG00000113328 ENSG00000072571 ENSG00000038274 
             "SLU7"         "PTTG1"         "CCNG1"          "HMMR"         "MAT2B" 
    ENSG00000113645 ENSG00000113643 ENSG00000120137 ENSG00000184347 ENSG00000040275 
             "WWC1"         "RARS1"         "PANK3"         "SLIT3"         "SPDL1" 
    ENSG00000204767 ENSG00000181163 ENSG00000072803 ENSG00000072786 ENSG00000214357 
          "INSYN2B"          "NPM1"        "FBXW11"         "STK10"       "NEURL1B" 
    ENSG00000120129 ENSG00000113719 ENSG00000113732 ENSG00000164463 ENSG00000113734 
            "DUSP1"        "ERGIC1"      "ATP6V0E1"        "CREBRF"         "BNIP1" 
    ENSG00000145919 ENSG00000113742 ENSG00000120149 ENSG00000184845 ENSG00000164466 
             "BOD1"         "CPEB4"          "MSX2"          "DRD1"         "SFXN1" 
    ENSG00000051596 ENSG00000170085 ENSG00000122203 ENSG00000175414 ENSG00000146066 
            "THOC3"         "SIMC1"      "KIAA1191"         "ARL10"        "HIGD2A" 
    ENSG00000175416 ENSG00000113194 ENSG00000146083 ENSG00000074317 ENSG00000087206 
             "CLTB"          "FAF2"         "RNF44"          "SNCB"         "UIMC1" 
    ENSG00000113761 ENSG00000165671 ENSG00000169228 ENSG00000213347 ENSG00000169230 
           "ZNF346"          "NSD1"         "RAB24"          "MXD3"       "PRELID1" 
    ENSG00000169223 ENSG00000198055 ENSG00000196923 ENSG00000027847 ENSG00000145916 
            "LMAN2"          "GRK6"        "PDLIM7"       "B4GALT7"        "RMND5B" 
    ENSG00000197451 ENSG00000113240 ENSG00000169131 ENSG00000178338 ENSG00000178187 
          "HNRNPAB"          "CLK4"       "ZNF354A"       "ZNF354B"        "ZNF454" 
    ENSG00000177932 ENSG00000087116 ENSG00000176783 ENSG00000169045 ENSG00000127022 
          "ZNF354C"       "ADAMTS2"         "RUFY1"       "HNRNPH1"          "CANX" 
    ENSG00000161021 ENSG00000161013 ENSG00000161011 ENSG00000161010 ENSG00000197226 
            "MAML1"        "MGAT4B"        "SQSTM1"         "MRNIP"       "TBC1D9B" 
    ENSG00000131459 ENSG00000113300 ENSG00000146063 ENSG00000204628 ENSG00000112685 
            "GFPT2"         "CNOT6"        "TRIM41"         "RACK1"         "EXOC2" 
    ENSG00000164379 ENSG00000137273 ENSG00000054598 ENSG00000112699 ENSG00000124535 
            "FOXQ1"         "FOXF2"         "FOXC1"          "GMDS"        "WRNIP1" 
    ENSG00000021355 ENSG00000124588 ENSG00000137274 ENSG00000137267 ENSG00000137285 
         "SERPINB1"          "NQO2"          "BPHL"        "TUBB2A"        "TUBB2B" 
    ENSG00000137266 ENSG00000168994 ENSG00000145945 ENSG00000112739 ENSG00000198721 
         "SLC22A23"         "PXDC1"        "FAM50B"         "PRP4K"          "ECI2" 
    ENSG00000153046 ENSG00000214113 ENSG00000124783 ENSG00000096696 ENSG00000239264 
             "CDYL"         "LYRM4"          "SSR1"           "DSP"        "TXNDC5" 
    ENSG00000124786 ENSG00000137203 ENSG00000111846 ENSG00000111843 ENSG00000137210 
          "SLC35B3"        "TFAP2A"         "GCNT2"       "TMEM14C"       "TMEM14B" 
    ENSG00000197977 ENSG00000224531 ENSG00000111859 ENSG00000145979 ENSG00000145990 
           "ELOVL2"        "SMIM13"         "NEDD9"        "TBC1D7"         "GFOD1" 
    ENSG00000124523 ENSG00000225921 ENSG00000010017 ENSG00000050393 ENSG00000180537 
            "SIRT5"          "NOL7"        "RANBP9"         "MCUR1"        "RNF182" 
    ENSG00000112149 ENSG00000008083 ENSG00000137198 ENSG00000124788 ENSG00000112186 
             "CD83"        "JARID2"          "GMPR"         "ATXN1"          "CAP2" 
    ENSG00000137414 ENSG00000137364 ENSG00000165097 ENSG00000124795 ENSG00000137393 
           "FAM8A1"          "TPMT"         "KDM1B"           "DEK"       "RNF144B" 
    ENSG00000172201 ENSG00000172197 ENSG00000145996 ENSG00000124766 ENSG00000124532 
              "ID4"        "MBOAT1"        "CDKAL1"          "SOX4"          "MRS2" 
    ENSG00000112294 ENSG00000111802 ENSG00000112308 ENSG00000112312 ENSG00000112343 
          "ALDH5A1"          "TDP2"       "C6orf62"          "GMNN"        "TRIM38" 
    ENSG00000010704 ENSG00000180573 ENSG00000158373 ENSG00000273802 ENSG00000275713 
              "HFE"         "H2AC6"         "H2BC5"         "H2BC8"         "H2BC9" 
    ENSG00000158406 ENSG00000186470 ENSG00000124508 ENSG00000026950 ENSG00000111801 
             "H4C8"        "BTN3A2"        "BTN2A2"        "BTN3A1"        "BTN3A3" 
    ENSG00000112763 ENSG00000182952 ENSG00000181315 ENSG00000124635 ENSG00000096654 
           "BTN2A1"         "HMGN4"        "ZNF322"        "H2BC11"        "ZNF184" 
    ENSG00000198315 ENSG00000137185 ENSG00000137338 ENSG00000158691 ENSG00000204713 
          "ZKSCAN8"        "ZSCAN9"         "PGBD1"       "ZSCAN12"        "TRIM27" 
    ENSG00000206503 ENSG00000066379 ENSG00000204619 ENSG00000204592 ENSG00000204569 
            "HLA-A"        "POLR1H"       "PPP1R11"         "HLA-E"       "PPP1R10" 
    ENSG00000137343 ENSG00000146112 ENSG00000137404 ENSG00000137337 ENSG00000196230 
            "ATAT1"       "PPP1R18"           "NRM"          "MDC1"          "TUBB" 
    ENSG00000137312 ENSG00000137331 ENSG00000204580 ENSG00000204536 ENSG00000137310 
            "FLOT1"          "IER3"          "DDR1"        "CCHCR1"         "TCF19" 
    ENSG00000204525 ENSG00000234745 ENSG00000204520 ENSG00000204516 ENSG00000198563 
            "HLA-C"         "HLA-B"          "MICA"          "MICB"        "DDX39B" 
    ENSG00000204469 ENSG00000204463 ENSG00000204439 ENSG00000204438 ENSG00000204435 
           "PRRC2A"          "BAG6"       "C6orf47"        "GPANK1"        "CSNK2B" 
    ENSG00000204428 ENSG00000204427 ENSG00000213722 ENSG00000213719 ENSG00000204394 
           "LY6G5C"       "ABHD16A"         "DDAH2"         "CLIC1"         "VARS1" 
    ENSG00000204389 ENSG00000204388 ENSG00000204387 ENSG00000204386 ENSG00000204371 
           "HSPA1A"        "HSPA1B"        "SNHG32"          "NEU1"         "EHMT2" 
    ENSG00000204356 ENSG00000204344 ENSG00000213676 ENSG00000221988 ENSG00000204310 
            "NELFE"         "STK19"         "ATF6B"          "PPT2"        "AGPAT1" 
    ENSG00000204308 ENSG00000204304 ENSG00000204267 ENSG00000204264 ENSG00000168394 
             "RNF5"          "PBX2"          "TAP2"         "PSMB8"          "TAP1" 
    ENSG00000242574 ENSG00000204257 ENSG00000204256 ENSG00000204231 ENSG00000112473 
          "HLA-DMB"       "HLA-DMA"          "BRD2"          "RXRB"       "SLC39A7" 
    ENSG00000204227 ENSG00000231500 ENSG00000227057 ENSG00000237441 ENSG00000231925 
            "RING1"         "RPS18"         "WDR46"          "RGL2"         "TAPBP" 
    ENSG00000236104 ENSG00000204209 ENSG00000237649 ENSG00000197283 ENSG00000030110 
           "ZBTB22"          "DAXX"         "KIFC1"       "SYNGAP1"          "BAK1" 
    ENSG00000096433 ENSG00000137288 ENSG00000161904 ENSG00000124493 ENSG00000137309 
            "ITPR3"         "UQCC2"         "LEMD2"          "GRM4"         "HMGA1" 
    ENSG00000186577 ENSG00000196821 ENSG00000124562 ENSG00000065060 ENSG00000064995 
           "SMIM29"         "ILRUN"         "SNRPC"        "BLTP3A"         "TAF11" 
    ENSG00000146197 ENSG00000023892 ENSG00000112039 ENSG00000198755 ENSG00000007866 
           "SCUBE3"          "DEF6"         "FANCE"        "RPL10A"         "TEAD3" 
    ENSG00000096060 ENSG00000156711 ENSG00000112078 ENSG00000112081 ENSG00000124762 
            "FKBP5"        "MAPK13"        "KCTD20"         "SRSF3"        "CDKN1A" 
    ENSG00000137168 ENSG00000137409 ENSG00000137193 ENSG00000137200 ENSG00000198937 
            "PPIL1"         "MTCH1"          "PIM1"         "CMTR1"       "CCDC167" 
    ENSG00000112139 ENSG00000156639 ENSG00000183826 ENSG00000124767 ENSG00000112167 
            "MDGA1"        "ZFAND3"         "BTBD9"          "GLO1"        "SAYSD1" 
    ENSG00000146122 ENSG00000124615 ENSG00000124596 ENSG00000112561 ENSG00000188112 
            "DAAM2"         "MOCS1"         "OARD1"          "TFEB"      "C6orf132" 
    ENSG00000124496 ENSG00000024048 ENSG00000124659 ENSG00000146223 ENSG00000221821 
           "TRERF1"          "UBR2"          "TBCC"        "RPL7L1"      "C6orf226" 
    ENSG00000137161 ENSG00000124587 ENSG00000112640 ENSG00000124702 ENSG00000044090 
            "CNPY3"          "PEX6"       "PPP2R5D"        "KLHDC3"          "CUL7" 
    ENSG00000137171 ENSG00000112651 ENSG00000112655 ENSG00000112658 ENSG00000112659 
             "KLC4"         "MRPL2"          "PTK7"           "SRF"          "CUL9" 
    ENSG00000171467 ENSG00000124574 ENSG00000204052 ENSG00000137207 ENSG00000124571 
           "ZNF318"        "ABCC10"        "LRRC73"         "YIPF3"          "XPO5" 
    ENSG00000170734 ENSG00000172432 ENSG00000096080 ENSG00000112715 ENSG00000180992 
             "POLH"        "GTPBP2"       "MRPS18A"         "VEGFA"        "MRPL14" 
    ENSG00000137216 ENSG00000112759 ENSG00000096384 ENSG00000157593 ENSG00000124608 
          "TMEM63B"       "SLC29A1"      "HSP90AB1"       "SLC35B2"         "AARS2" 
    ENSG00000124813 ENSG00000172348 ENSG00000146072 ENSG00000198087 ENSG00000244694 
            "RUNX2"         "RCAN2"      "TNFRSF21"         "CD2AP"        "PTCHD4" 
    ENSG00000031691 ENSG00000112118 ENSG00000096093 ENSG00000065308 ENSG00000096092 
            "CENPQ"          "MCM3"         "EFHC1"         "TRAM2"       "TMEM14A" 
    ENSG00000112144 ENSG00000012660 ENSG00000001084 ENSG00000137269 ENSG00000151917 
            "CILK1"        "ELOVL5"          "GCLC"         "LRRC1"         "BEND6" 
    ENSG00000112200 ENSG00000112208 ENSG00000146143 ENSG00000112245 ENSG00000135298 
           "ZNF451"          "BAG2"         "PRIM2"        "PTP4A1"        "ADGRB3" 
    ENSG00000168216 ENSG00000082269 ENSG00000112305 ENSG00000119900 ENSG00000079841 
           "LMBRD1"       "FAM135A"         "SMAP1"        "OGFRL1"         "RIMS1" 
    ENSG00000185760 ENSG00000135297 ENSG00000156508 ENSG00000119899 ENSG00000156535 
            "KCNQ5"          "MTO1"        "EEF1A1"       "SLC17A5"         "CD109" 
    ENSG00000111799 ENSG00000112695 ENSG00000112697 ENSG00000196586 ENSG00000146243 
          "COL12A1"        "COX7A2"       "TMEM30A"          "MYO6"      "IRAK1BP1" 
    ENSG00000146247 ENSG00000118418 ENSG00000135338 ENSG00000118402 ENSG00000112742 
             "PHIP"         "HMGN3"          "LCA5"        "ELOVL4"           "TTK" 
    ENSG00000083123 ENSG00000112773 ENSG00000146242 ENSG00000083097 ENSG00000013375 
           "BCKDHB"        "TENT5A"          "TPBG"         "DOP1A"          "PGM3" 
    ENSG00000013392 ENSG00000065833 ENSG00000146250 ENSG00000065615 ENSG00000135315 
           "RWDD2A"           "ME1"        "PRSS35"        "CYB5R4"        "CEP162" 
    ENSG00000135318 ENSG00000135317 ENSG00000188994 ENSG00000146282 ENSG00000135336 
             "NT5E"         "SNX14"        "ZNF292"         "RARS2"          "ORC3" 
    ENSG00000135334 ENSG00000111880 ENSG00000146278 ENSG00000146281 ENSG00000198833 
          "AKIRIN2"         "RNGTT"         "PNRC1"        "PM20D2"        "UBE2J1" 
    ENSG00000025039 ENSG00000112159 ENSG00000135333 ENSG00000014123 ENSG00000123545 
            "RRAGD"          "MDN1"         "EPHA7"          "UFL1"       "NDUFAF4" 
    ENSG00000146263 ENSG00000112234 ENSG00000132423 ENSG00000132424 ENSG00000123552 
           "MMS22L"         "FBXL4"          "COQ3"         "PNISR"         "USP45" 
    ENSG00000112237 ENSG00000112249 ENSG00000085382 ENSG00000112276 ENSG00000057657 
             "CCNC"         "ASCC3"         "HACE1"        "POPDC1"         "PRDM1" 
    ENSG00000112297 ENSG00000130347 ENSG00000130348 ENSG00000178409 ENSG00000112320 
           "CRYBG1"       "RTN4IP1"         "QRSL1"         "BEND3"          "SOBP" 
    ENSG00000025796 ENSG00000112335 ENSG00000118689 ENSG00000080546 ENSG00000183137 
            "SEC63"          "SNX3"         "FOXO3"         "SESN1"       "CEP57L1" 
    ENSG00000112367 ENSG00000168438 ENSG00000155111 ENSG00000155115 ENSG00000197498 
             "FIG4"         "CDC40"         "CDK19"        "GTF3C6"          "RPF2" 
    ENSG00000173214 ENSG00000009413 ENSG00000056972 ENSG00000010810 ENSG00000074935 
           "MFSD4B"         "REV3L"      "TRAF3IP2"           "FYN"         "TUBE1" 
    ENSG00000112769 ENSG00000277443 ENSG00000196591 ENSG00000111817 ENSG00000189241 
            "LAMA4"        "MARCKS"         "HDAC2"           "DSE"        "TSPYL1" 
    ENSG00000178033 ENSG00000111832 ENSG00000153975 ENSG00000196911 ENSG00000183807 
           "CALHM5"         "RWDD1"          "ZUP1"         "KPNA5"       "FAM162B" 
    ENSG00000164465 ENSG00000047932 ENSG00000153989 ENSG00000196376 ENSG00000111860 
           "DCBLD1"          "GOPC"          "NUS1"       "SLC35F1"        "CEP85L" 
    ENSG00000198523 ENSG00000111875 ENSG00000111885 ENSG00000152661 ENSG00000111897 
              "PLN"         "ASF1A"        "MAN1A1"          "GJA1"       "SERINC1" 
    ENSG00000172594 ENSG00000146373 ENSG00000111907 ENSG00000111906 ENSG00000111912 
          "SMPDL3A"        "RNF217"       "TPD52L1"         "HDDC2"         "NCOA7" 
    ENSG00000111911 ENSG00000203760 ENSG00000146374 ENSG00000118518 ENSG00000093144 
            "HINT3"         "CENPW"         "RSPO3"        "RNF146"        "ECHDC1" 
    ENSG00000152894 ENSG00000196569 ENSG00000146376 ENSG00000198945 ENSG00000164484 
            "PTPRK"         "LAMA2"      "ARHGAP18"       "L3MBTL3"      "TMEM200A" 
    ENSG00000079819 ENSG00000112282 ENSG00000197594 ENSG00000118523 ENSG00000079931 
          "EPB41L2"         "MED23"         "ENPP1"          "CCN2"         "MOXD1" 
    ENSG00000112306 ENSG00000112319 ENSG00000118526 ENSG00000028839 ENSG00000146411 
            "RPS12"          "EYA4"         "TCF21"         "TBPL1"       "SLC2A12" 
    ENSG00000118515 ENSG00000118514 ENSG00000112339 ENSG00000135541 ENSG00000146410 
             "SGK1"       "ALDH8A1"         "HBS1L"          "AHI1"         "MTFR2" 
    ENSG00000029363 ENSG00000135525 ENSG00000197442 ENSG00000118503 ENSG00000112378 
           "BCLAF1"          "MAP7"        "MAP3K5"       "TNFAIP3"          "PERP" 
    ENSG00000024862 ENSG00000146386 ENSG00000112406 ENSG00000009844 ENSG00000112414 
          "CCDC28A"        "ABRACL"          "HECA"          "VTA1"        "ADGRG6" 
    ENSG00000010818 ENSG00000146416 ENSG00000034693 ENSG00000001036 ENSG00000112419 
           "HIVEP2"          "AIG1"          "PEX3"         "FUCA2"       "PHACTR2" 
    ENSG00000118495 ENSG00000152818 ENSG00000118496 ENSG00000146414 ENSG00000164506 
           "PLAGL1"          "UTRN"        "FBXO30"         "SHPRH"        "STXBP5" 
    ENSG00000203727 ENSG00000055208 ENSG00000186625 ENSG00000120253 ENSG00000120265 
            "SAMD5"          "TAB2"        "KATNA1"         "NUP43"         "PCMT1" 
    ENSG00000120256 ENSG00000203722 ENSG00000131015 ENSG00000111981 ENSG00000120278 
            "LRP11"        "RAET1G"         "ULBP2"         "ULBP1"       "PLEKHG1" 
    ENSG00000120254 ENSG00000131016 ENSG00000181472 ENSG00000146476 ENSG00000131018 
          "MTHFD1L"        "AKAP12"         "ZBTB2"         "ARMT1"         "SYNE1" 
    ENSG00000120279 ENSG00000112029 ENSG00000112031 ENSG00000091844 ENSG00000153721 
            "MYCT1"         "FBXO5"        "MTRF1L"         "RGS17"        "CNKSR3" 
    ENSG00000213079 ENSG00000029639 ENSG00000049618 ENSG00000215712 ENSG00000130340 
            "SCAF8"         "TFB1M"        "ARID1B"       "TMEM242"          "SNX9" 
    ENSG00000078269 ENSG00000122335 ENSG00000146433 ENSG00000146425 ENSG00000092820 
            "SYNJ2"        "SERAC1"       "TMEM181"        "DYNLT1"           "EZR" 
    ENSG00000112096 ENSG00000120437 ENSG00000120438 ENSG00000197081 ENSG00000085511 
                 NA         "ACAT2"          "TCP1"         "IGF2R"        "MAP3K4" 
    ENSG00000026652 ENSG00000112531 ENSG00000112541 ENSG00000198818 ENSG00000071242 
           "AGPAT4"           "QKI"        "PDE10A"        "SFT2D1"       "RPS6KA2" 
    ENSG00000130396 ENSG00000184465 ENSG00000185127 ENSG00000130024 ENSG00000198719 
             "AFDN"         "WDR27"      "C6orf120"         "PHF10"          "DLL1" 
    ENSG00000112584 ENSG00000177706 ENSG00000197461 ENSG00000188191 ENSG00000164828 
          "FAM120B"        "FAM20C"         "PDGFA"       "PRKAR1B"          "SUN1" 
    ENSG00000105963 ENSG00000146540 ENSG00000164850 ENSG00000178381 ENSG00000164877 
            "ADAP1"       "C7orf50"         "GPER1"       "ZFAND2A"       "MICALL2" 
    ENSG00000198517 ENSG00000225968 ENSG00000106268 ENSG00000106266 ENSG00000106263 
             "MAFK"         "ELFN1"         "NUDT1"          "SNX8"         "EIF3B" 
    ENSG00000136213 ENSG00000106009 ENSG00000106012 ENSG00000136295 ENSG00000174945 
           "CHST12"         "BRAT1"          "IQCE"         "TTYH3"          "AMZ1" 
    ENSG00000146535 ENSG00000198286 ENSG00000146555 ENSG00000242802 ENSG00000146587 
            "GNA12"        "CARD11"          "SDK1"         "AP5Z1"          "RBAK" 
    ENSG00000157954 ENSG00000164638 ENSG00000182095 ENSG00000075624 ENSG00000075618 
            "WIPI2"       "SLC29A4"        "TNRC18"          "ACTB"         "FSCN1" 
    ENSG00000011275 ENSG00000122512 ENSG00000106305 ENSG00000086232 ENSG00000178397 
           "RNF216"          "PMS2"         "AIMP2"       "EIF2AK1"       "FAM220A" 
    ENSG00000136238 ENSG00000164535 ENSG00000136247 ENSG00000146576 ENSG00000205903 
             "RAC1"         "DAGLB"        "ZDHHC4"        "INTS15"        "ZNF316" 
    ENSG00000164631 ENSG00000106392 ENSG00000164654 ENSG00000106399 ENSG00000106415 
            "ZNF12"       "C1GALT1"          "MIOS"          "RPA3"        "GLCCI1" 
    ENSG00000189043 ENSG00000106443 ENSG00000005108 ENSG00000106460 ENSG00000122644 
           "NDUFA4"         "PHF14"        "THSD7A"      "TMEM106B"         "ARL4A" 
    ENSG00000006468 ENSG00000106524 ENSG00000106537 ENSG00000106546 ENSG00000071189 
             "ETV1"        "ANKMY2"       "TSPAN13"           "AHR"         "SNX13" 
    ENSG00000048052 ENSG00000122691 ENSG00000105849 ENSG00000105855 ENSG00000164649 
            "HDAC9"        "TWIST1"        "POLR1F"         "ITGB8"        "CDCA7L" 
    ENSG00000105889 ENSG00000136244 ENSG00000122591 ENSG00000122550 ENSG00000156928 
          "STEAP1B"           "IL6"         "HYCC1"         "KLHL7"        "MALSU1" 
    ENSG00000136231 ENSG00000164548 ENSG00000169193 ENSG00000105926 ENSG00000105928 
          "IGF2BP3"         "TRA2A"       "CCDC126"         "PALS2"         "GSDME" 
    ENSG00000172115 ENSG00000050344 ENSG00000122566 ENSG00000122565 ENSG00000005020 
             "CYCS"        "NFE2L3"     "HNRNPA2B1"          "CBX3"         "SKAP2" 
    ENSG00000105991 ENSG00000105996 ENSG00000105997 ENSG00000197576 ENSG00000106004 
            "HOXA1"         "HOXA2"         "HOXA3"         "HOXA4"         "HOXA5" 
    ENSG00000106052 ENSG00000153814 ENSG00000146592 ENSG00000136193 ENSG00000106080 
          "TAX1BP1"         "JAZF1"         "CREB5"         "SCRN1"        "FKBP14" 
    ENSG00000106086 ENSG00000180354 ENSG00000106100 ENSG00000106105 ENSG00000240583 
          "PLEKHA8"         "MTURN"          "NOD1"         "GARS1"          "AQP1" 
    ENSG00000106355 ENSG00000105778 ENSG00000170852 ENSG00000122642 ENSG00000122643 
             "LSM5"          "AVL9"        "KBTBD2"         "FKBP9"        "NT5C3A" 
    ENSG00000122507 ENSG00000164619 ENSG00000173852 ENSG00000122557 ENSG00000122547 
             "BBS9"         "BMPER"       "DPY19L1"       "HERPUD2"         "EEPD1" 
    ENSG00000011426 ENSG00000010270 ENSG00000078053 ENSG00000006715 ENSG00000006451 
             "ANLN"      "STARD3NL"          "AMPH"         "VPS41"          "RALA" 
    ENSG00000168303 ENSG00000175600 ENSG00000122641 ENSG00000164543 ENSG00000106603 
           "MPLKIP"         "SUGCT"         "INHBA"        "STK17A"          "COA1" 
    ENSG00000106605 ENSG00000136279 ENSG00000106624 ENSG00000106628 ENSG00000015520 
            "BLVRA"          "DBNL"         "AEBP1"         "POLD2"        "NPC1L1" 
    ENSG00000136271 ENSG00000158604 ENSG00000105953 ENSG00000122515 ENSG00000196262 
            "DDX56"         "TMED4"          "OGDH"         "ZMIZ2"          "PPIA" 
    ENSG00000105968 ENSG00000146676 ENSG00000136274 ENSG00000146674 ENSG00000136205 
            "H2AZ2"          "PURB"         "NACAD"        "IGFBP3"          "TNS3" 
    ENSG00000136273 ENSG00000158683 ENSG00000183696 ENSG00000132436 ENSG00000106070 
             "HUS1"        "PKD1L1"          "UPP1"        "FIGNL1"         "GRB10" 
    ENSG00000146648 ENSG00000178665 ENSG00000146729 ENSG00000239789 ENSG00000146733 
             "EGFR"        "ZNF713"      "NIPSNAP2"        "MRPS17"          "PSPH" 
    ENSG00000146731 ENSG00000129103 ENSG00000106153 ENSG00000234444 ENSG00000173041 
            "CCT6A"         "SUMF2"        "CHCHD2"        "ZNF736"        "ZNF680" 
    ENSG00000152926 ENSG00000213462 ENSG00000146757 ENSG00000196715 ENSG00000126522 
           "ZNF117"        "ERV3-1"         "ZNF92"      "VKORC1L1"           "ASL" 
    ENSG00000241258 ENSG00000169902 ENSG00000126524 ENSG00000196313 ENSG00000009954 
             "CRCP"         "TPST1"          "SBDS"        "POM121"         "BAZ1B" 
    ENSG00000071462 ENSG00000106077 ENSG00000189143 ENSG00000165171 ENSG00000049540 
            "BUD23"        "ABHD11"         "CLDN4"       "METTL27"           "ELN" 
    ENSG00000049541 ENSG00000106665 ENSG00000263001 ENSG00000272391 ENSG00000127946 
             "RFC2"         "CLIP2"         "GTF2I"       "POM121C"          "HIP1" 
    ENSG00000005486 ENSG00000127948 ENSG00000127952 ENSG00000146701 ENSG00000177679 
           "RHBDD2"           "POR"        "STYXL1"          "MDH2"         "SRRM3" 
    ENSG00000106211 ENSG00000170027 ENSG00000188372 ENSG00000091073 ENSG00000135205 
            "HSPB1"         "YWHAG"           "ZP3"          "DTX2"       "CCDC146" 
    ENSG00000127951 ENSG00000186088 ENSG00000127947 ENSG00000135211 ENSG00000187391 
             "FGL2"          "GSAP"        "PTPN12"        "TMEM60"         "MAGI2" 
    ENSG00000127955 ENSG00000135218 ENSG00000075223 ENSG00000019991 ENSG00000075213 
            "GNAI1"          "CD36"        "SEMA3C"           "HGF"        "SEMA3A" 
    ENSG00000153993 ENSG00000135185 ENSG00000005469 ENSG00000085563 ENSG00000006634 
           "SEMA3D"       "TMEM243"          "CROT"         "ABCB1"          "DBF4" 
    ENSG00000075142 ENSG00000164647 ENSG00000105792 ENSG00000157224 ENSG00000157240 
              "SRI"        "STEAP1"        "CFAP69"        "CLDN12"          "FZD1" 
    ENSG00000127914 ENSG00000127980 ENSG00000127993 ENSG00000105810 ENSG00000004766 
            "AKAP9"          "PEX1"         "RBM48"          "CDK6"         "VPS50" 
    ENSG00000105825 ENSG00000127920 ENSG00000105829 ENSG00000164692 ENSG00000127995 
            "TFPI2"         "GNG11"          "BET1"        "COL1A2"         "CASD1" 
    ENSG00000127990 ENSG00000242265 ENSG00000105854 ENSG00000004799 ENSG00000158560 
             "SGCE"         "PEG10"          "PON2"          "PDK4"       "DYNC1I1" 
    ENSG00000004864 ENSG00000127922 ENSG00000070669 ENSG00000164715 ENSG00000205356 
         "SLC25A13"          "SEM1"          "ASNS"         "LMTK2"        "TECPR1" 
    ENSG00000164713 ENSG00000006453 ENSG00000166448 ENSG00000196367 ENSG00000198742 
             "BRI3"      "BAIAP2L1"       "TMEM130"         "TRRAP"        "SMURF1" 
    ENSG00000241685 ENSG00000106244 ENSG00000160917 ENSG00000241468 ENSG00000160908 
           "ARPC1A"         "PDAP1"         "CPSF4"        "ATP5MF"        "ZNF394" 
    ENSG00000221909 ENSG00000197037 ENSG00000146833 ENSG00000166526 ENSG00000168090 
          "FAM200A"       "ZSCAN25"         "TRIM4"          "ZNF3"         "COPS6" 
    ENSG00000166508 ENSG00000166997 ENSG00000214309 ENSG00000146826 ENSG00000213420 
             "MCM7"         "CNPY4"        "MBLAC1"      "TRAPPC14"          "GPC2" 
    ENSG00000146834 ENSG00000166925 ENSG00000106351 ENSG00000077454 ENSG00000106333 
            "MEPCE"       "TSC22D4"         "AGFG2"         "LRCH4"        "PCOLCE" 
    ENSG00000106330 ENSG00000146830 ENSG00000196411 ENSG00000087077 ENSG00000087087 
           "MOSPD3"        "GIGYF1"         "EPHB4"         "TRIP6"          "SRRT" 
    ENSG00000087085 ENSG00000169871 ENSG00000106366 ENSG00000128564 ENSG00000106397 
             "ACHE"        "TRIM56"      "SERPINE1"           "VGF"         "PLOD3" 
    ENSG00000106400 ENSG00000257923 ENSG00000160999 ENSG00000160991 ENSG00000161036 
           "ZNHIT1"          "CUX1"         "SH2B2"         "ORAI2"         "LRWD1" 
    ENSG00000128606 ENSG00000170632 ENSG00000161048 ENSG00000161057 ENSG00000189056 
           "LRRC17"        "ARMC10"       "NAPEPLD"         "PSMC2"          "RELN" 
    ENSG00000164815 ENSG00000005483 ENSG00000135250 ENSG00000135249 ENSG00000146776 
             "ORC5"         "KMT2E"         "SRPK2"         "RINT1"       "ATXN7L1" 
    ENSG00000008282 ENSG00000105835 ENSG00000253276 ENSG00000005249 ENSG00000105856 
            "SYPL1"         "NAMPT"       "CCDC71L"       "PRKAR2B"          "HBP1" 
    ENSG00000164597 ENSG00000105879 ENSG00000091140 ENSG00000091136 ENSG00000135241 
             "COG5"         "CBLL1"           "DLD"         "LAMB1"        "PNPLA8" 
    ENSG00000177683 ENSG00000128590 ENSG00000173114 ENSG00000128512 ENSG00000198839 
            "THAP5"        "DNAJB9"         "LRRN3"         "DOCK4"        "ZNF277" 
    ENSG00000146802 ENSG00000164603 ENSG00000135272 ENSG00000135269 ENSG00000105971 
          "TMEM168"        "SAMTOR"         "MDFIC"           "TES"          "CAV2" 
    ENSG00000105974 ENSG00000105976 ENSG00000198898 ENSG00000105989 ENSG00000077063 
             "CAV1"           "MET"        "CAPZA2"          "WNT2"       "CTTNBP2" 
    ENSG00000128534 ENSG00000184408 ENSG00000106025 ENSG00000106034 ENSG00000002745 
             "LSM8"         "KCND2"       "TSPAN12"         "CPED1"         "WNT16" 
    ENSG00000128609 ENSG00000106299 ENSG00000170775 ENSG00000128513 ENSG00000004059 
           "NDUFA5"          "WASL"         "GPR37"          "POT1"          "ARF5" 
    ENSG00000197157 ENSG00000106348 ENSG00000165055 ENSG00000128595 ENSG00000128596 
             "SND1"        "IMPDH1"       "METTL2B"          "CALU"       "CCDC136" 
    ENSG00000128591 ENSG00000128524 ENSG00000128602 ENSG00000158467 ENSG00000128578 
             "FLNC"       "ATP6V1F"           "SMO"        "AHCYL2"        "STRIP2" 
    ENSG00000186591 ENSG00000091732 ENSG00000128607 ENSG00000128510 ENSG00000106477 
            "UBE2H"        "ZC3HC1"       "KLHDC10"          "CPA4"         "CEP41" 
    ENSG00000106484 ENSG00000128585 ENSG00000106554 ENSG00000131558 ENSG00000205060 
             "MEST"         "MKLN1"        "CHCHD3"         "EXOC4"       "SLC35B4" 
    ENSG00000085662 ENSG00000198074 ENSG00000172331 ENSG00000122786 ENSG00000146859 
           "AKR1B1"       "AKR1B10"          "BPGM"         "CALD1"       "TMEM140" 
    ENSG00000155561 ENSG00000189320 ENSG00000105887 ENSG00000181072 ENSG00000105894 
           "NUP205"       "FAM180A"          "MTPN"         "CHRM2"           "PTN" 
    ENSG00000157680 ENSG00000182158 ENSG00000122779 ENSG00000122778 ENSG00000146858 
             "DGKI"       "CREB3L2"        "TRIM24"      "KIAA1549"      "ZC3HAV1L" 
    ENSG00000105939 ENSG00000105948 ENSG00000157741 ENSG00000146963 ENSG00000064393 
          "ZC3HAV1"         "IFT56"          "UBN2"        "LUC7L2"         "HIPK2" 
    ENSG00000059377 ENSG00000059378 ENSG00000006459 ENSG00000133606 ENSG00000133597 
           "TBXAS1"        "PARP12"         "KDM7A"         "MKRN1"         "ADCK2" 
    ENSG00000006530 ENSG00000106028 ENSG00000159784 ENSG00000159840 ENSG00000170379 
              "AGK"         "SSBP1"       "FAM131B"           "ZYX"         "TCAF2" 
    ENSG00000055130 ENSG00000106462 ENSG00000155660 ENSG00000170265 ENSG00000170260 
             "CUL1"          "EZH2"         "PDIA4"        "ZNF282"        "ZNF212" 
    ENSG00000204946 ENSG00000196453 ENSG00000106479 ENSG00000171130 ENSG00000214022 
           "ZNF783"        "ZNF777"        "ZNF862"      "ATP6V0E2"        "REPIN1" 
    ENSG00000055118 ENSG00000197150 ENSG00000164885 ENSG00000164889 ENSG00000164896 
            "KCNH2"         "ABCB8"          "CDK5"        "SLC4A2"         "FASTK" 
    ENSG00000133612 ENSG00000033050 ENSG00000033100 ENSG00000082014 ENSG00000013374 
            "AGAP3"         "ABCF2"         "CHPF2"       "SMARCD3"          "NUB1" 
    ENSG00000106617 ENSG00000178234 ENSG00000055609 ENSG00000196584 ENSG00000133627 
           "PRKAG2"       "GALNT11"         "KMT2C"         "XRCC2"        "ACTR3B" 
    ENSG00000186480 ENSG00000164690 ENSG00000009335 ENSG00000105993 ENSG00000146918 
           "INSIG1"           "SHH"         "UBE3C"        "DNAJB6"        "NCAPG2" 
    ENSG00000117868 ENSG00000182378 ENSG00000185291 ENSG00000169100 ENSG00000169084 
            "ESYT2"        "PLCXD1"         "IL3RA"       "SLC25A6"         "DHRSX" 
    ENSG00000214717 ENSG00000006756 ENSG00000157399 ENSG00000101825 ENSG00000183943 
            "ZBED1"          "ARSD"          "ARSL"         "MXRA5"          "PRKX" 
    ENSG00000130021 ENSG00000101846 ENSG00000006757 ENSG00000101849 ENSG00000146950 
             "PUDP"           "STS"        "PNPLA4"         "TBL1X"       "SHROOM2" 
    ENSG00000047644 ENSG00000101871 ENSG00000047648 ENSG00000101911 ENSG00000205542 
             "WWC3"          "MID1"       "ARHGAP6"         "PRPS2"        "TMSB4X" 
    ENSG00000123595 ENSG00000196459 ENSG00000046653 ENSG00000181544 ENSG00000130150 
            "RAB9A"       "TRAPPC2"         "GPM6B"         "FANCB"        "MOSPD2" 
    ENSG00000169239 ENSG00000182287 ENSG00000126010 ENSG00000047230 ENSG00000102054 
             "CA5B"         "AP1S2"          "GRPR"         "CTPS2"         "RBBP7" 
    ENSG00000102098 ENSG00000008086 ENSG00000044446 ENSG00000131828 ENSG00000147010 
            "SCML2"         "CDKL5"         "PHKA2"         "PDHA1"       "SH3KBP1" 
    ENSG00000173674 ENSG00000177189 ENSG00000012174 ENSG00000102172 ENSG00000102174 
           "EIF1AX"       "RPS6KA3"        "MBTPS2"           "SMS"          "PHEX" 
    ENSG00000123130 ENSG00000130066 ENSG00000184831 ENSG00000130741 ENSG00000101868 
            "ACOT9"          "SAT1"          "APOO"        "EIF2S3"         "POLA1" 
    ENSG00000169306 ENSG00000198814 ENSG00000198947 ENSG00000147027 ENSG00000130962 
         "IL1RAPL1"            "GK"           "DMD"        "TMEM47"         "PRRG1" 
    ENSG00000165169 ENSG00000147041 ENSG00000101955 ENSG00000156313 ENSG00000165175 
           "DYNLT3"         "SYTL5"          "SRPX"          "RPGR"       "MID1IP1" 
    ENSG00000183337 ENSG00000182220 ENSG00000185753 ENSG00000180182 ENSG00000147044 
             "BCOR"       "ATP6AP2"       "CXorf38"         "MED14"          "CASK" 
    ENSG00000189221 ENSG00000069509 ENSG00000147050 ENSG00000065923 ENSG00000102221 
             "MAOA"        "FUNDC1"         "KDM6A"        "SLC9A7"         "JADE3" 
    ENSG00000182872 ENSG00000130985 ENSG00000102225 ENSG00000102226 ENSG00000197779 
            "RBM10"          "UBA1"         "CDK16"         "USP11"         "ZNF81" 
    ENSG00000017483 ENSG00000068438 ENSG00000102317 ENSG00000101940 ENSG00000101945 
          "SLC38A5"         "FTSJ1"          "RBM3"         "WDR13"       "SUV39H1" 
    ENSG00000094631 ENSG00000102103 ENSG00000102100 ENSG00000068308 ENSG00000102057 
            "HDAC6"         "PQBP1"       "SLC35A2"         "OTUD5"         "KCND1" 
    ENSG00000147144 ENSG00000243279 ENSG00000196998 ENSG00000068394 ENSG00000102007 
          "CCDC120"         "PRAF2"         "WDR45"         "GPKOW"          "PLP2" 
    ENSG00000012211 ENSG00000171365 ENSG00000158352 ENSG00000189369 ENSG00000184194 
         "PRICKLE3"         "CLCN5"       "SHROOM4"         "GSPT2"        "GPR173" 
    ENSG00000126012 ENSG00000072501 ENSG00000072506 ENSG00000184083 ENSG00000130119 
            "KDM5C"         "SMC1A"      "HSD17B10"       "FAM120C"         "GNL3L" 
    ENSG00000187601 ENSG00000083750 ENSG00000188021 ENSG00000126970 ENSG00000147065 
           "MAGEH1"         "RRAGB"        "UBQLN2"         "ZC4H2"           "MSN" 
    ENSG00000089472 ENSG00000131080 ENSG00000130052 ENSG00000090776 ENSG00000181191 
             "HEPH"         "EDA2R"        "STARD8"         "EFNB1"          "PJA1" 
    ENSG00000120509 ENSG00000090889 ENSG00000184481 ENSG00000184634 ENSG00000147130 
           "PDZD11"         "KIF4A"         "FOXO4"         "MED12"         "ZMYM3" 
    ENSG00000147140 ENSG00000147166 ENSG00000147162 ENSG00000147174 ENSG00000204131 
             "NONO"      "ITGB1BP2"           "OGT"          "GCNA"         "NHSL2" 
    ENSG00000186871 ENSG00000067177 ENSG00000147100 ENSG00000131269 ENSG00000094841 
           "ERCC6L"         "PHKA1"       "SLC16A2"         "ABCB7"          "UPRT" 
    ENSG00000085224 ENSG00000102158 ENSG00000165240 ENSG00000102144 ENSG00000187325 
             "ATRX"         "MAGT1"         "ATP7A"          "PGK1"         "TAF9B" 
    ENSG00000179300 ENSG00000131171 ENSG00000072133 ENSG00000155008 ENSG00000102271 
             "RTL3"       "SH3BGRL"       "RPS6KA6"         "APOOL"         "KLHL4" 
    ENSG00000147202 ENSG00000165194 ENSG00000000003 ENSG00000102359 ENSG00000101811 
           "DIAPH2"        "PCDH19"        "TSPAN6"         "SRPX2"         "CSTF2" 
    ENSG00000182489 ENSG00000188917 ENSG00000126950 ENSG00000102384 ENSG00000102393 
             "XKRX"        "TRMT2B"       "TMEM35A"         "CENPI"           "GLA" 
    ENSG00000126945 ENSG00000126947 ENSG00000198960 ENSG00000184867 ENSG00000158164 
          "HNRNPH2"        "ARMCX1"        "ARMCX6"        "ARMCX2"       "TMSB15A" 
    ENSG00000158301 ENSG00000180964 ENSG00000182916 ENSG00000185222 ENSG00000166681 
          "GPRASP2"        "TCEAL8"        "TCEAL7"        "TCEAL9"          "BEX3" 
    ENSG00000133142 ENSG00000196507 ENSG00000123562 ENSG00000123575 ENSG00000123572 
           "TCEAL4"        "TCEAL3"       "MORF4L2"       "FAM199X"           "NRK" 
    ENSG00000157502 ENSG00000133138 ENSG00000133131 ENSG00000089682 ENSG00000147234 
           "PWWP3B"       "TBC1D8B"         "MORC4"         "RBM41"        "FRMPD3" 
    ENSG00000147224 ENSG00000157514 ENSG00000101843 ENSG00000101844 ENSG00000197565 
            "PRPS1"       "TSC22D3"        "PSMD10"         "ATG4A"        "COL4A6" 
    ENSG00000188153 ENSG00000068366 ENSG00000157600 ENSG00000101935 ENSG00000077264 
           "COL4A5"         "ACSL4"       "TMEM164"       "AMMECR1"          "PAK3" 
    ENSG00000072315 ENSG00000126016 ENSG00000123496 ENSG00000130224 ENSG00000102024 
            "TRPC5"          "AMOT"       "IL13RA2"         "LRCH2"          "PLS3" 
    ENSG00000003096 ENSG00000131724 ENSG00000174460 ENSG00000101856 ENSG00000077713 
           "KLHL13"       "IL13RA1"       "ZCCHC12"        "PGRMC1"      "SLC25A43" 
    ENSG00000005022 ENSG00000077721 ENSG00000186416 ENSG00000125354 ENSG00000198918 
          "SLC25A5"         "UBE2A"          "NKRF"       "SEPTIN6"         "RPL39" 
    ENSG00000125351 ENSG00000101882 ENSG00000177485 ENSG00000005893 ENSG00000171155 
            "UPF3B"          "NKAP"        "ZBTB33"         "LAMP2"     "C1GALT1C1" 
    ENSG00000125676 ENSG00000101966 ENSG00000101972 ENSG00000009694 ENSG00000188706 
            "THOC2"          "XIAP"         "STAG2"         "TENM1"        "ZDHHC9" 
    ENSG00000085185 ENSG00000102034 ENSG00000056277 ENSG00000102078 ENSG00000123728 
           "BCORL1"          "ELF4"       "ZNF280C"      "SLC25A14"         "RAP2C" 
    ENSG00000076716 ENSG00000156531 ENSG00000165704 ENSG00000156504 ENSG00000101928 
             "GPC4"          "PHF6"         "HPRT1"        "PABIR2"        "MOSPD1" 
    ENSG00000184785 ENSG00000212747 ENSG00000134590 ENSG00000203950 ENSG00000186376 
           "SMIM10"         "RTL8B"         "RTL8C"         "RTL8A"        "ZNF75D" 
    ENSG00000198689 ENSG00000022267 ENSG00000129680 ENSG00000147274 ENSG00000182195 
           "SLC9A6"          "FHL1"        "MAP7D3"          "RBMX"         "LDOC1" 
    ENSG00000102081 ENSG00000010404 ENSG00000269556 ENSG00000013619 ENSG00000063601 
             "FMR1"           "IDS"      "TMEM185A"        "MAMLD1"         "MTMR1" 
    ENSG00000160131 ENSG00000147383 ENSG00000182492 ENSG00000262919 ENSG00000130821 
            "VMA21"         "NSDHL"           "BGN"          "CCNQ"        "SLC6A8" 
    ENSG00000185825 ENSG00000180879 ENSG00000198910 ENSG00000102032 ENSG00000184216 
           "BCAP31"          "SSR4"         "L1CAM"         "RENBP"         "IRAK1" 
    ENSG00000169057 ENSG00000196924 ENSG00000147403 ENSG00000071553 ENSG00000203879 
            "MECP2"          "FLNA"         "RPL10"       "ATP6AP1"          "GDI1" 
    ENSG00000071859 ENSG00000130827 ENSG00000196976 ENSG00000102178 ENSG00000126903 
           "FAM50A"        "PLXNA3"         "LAGE3"         "UBL4A"       "SLC10A3" 
    ENSG00000071889 ENSG00000160211 ENSG00000160219 ENSG00000130826 ENSG00000130830 
            "FAM3A"          "G6PD"          "GAB3"          "DKC1"          "MPP1" 
    ENSG00000185010 ENSG00000165775 ENSG00000155961 ENSG00000185973 ENSG00000147364 
               "F8"        "FUNDC2"        "RAB39B"         "TMLHE"        "FBXO25" 
    ENSG00000180190 ENSG00000182372 ENSG00000104728 ENSG00000091879 ENSG00000275591 
             "TDRP"          "CLN8"      "ARHGEF10"        "ANGPT2"          "XKR5" 
    ENSG00000275342 ENSG00000253958 ENSG00000147324 ENSG00000104626 ENSG00000173273 
            "PRAG1"        "CLDN23"        "MFHAS1"          "ERI1"          "TNKS" 
    ENSG00000104643 ENSG00000154319 ENSG00000154328 ENSG00000079459 ENSG00000164733 
            "MTMR9"       "FAM167A"         "NEIL2"         "FDFT1"          "CTSB" 
    ENSG00000186523 ENSG00000154359 ENSG00000250305 ENSG00000164741 ENSG00000104723 
          "FAM86B1"        "LONRF1"        "TRMT9B"          "DLC1"         "TUSC3" 
    ENSG00000155970 ENSG00000104219 ENSG00000155975 ENSG00000003989 ENSG00000104213 
            "MICU3"        "ZDHHC2"        "VPS37A"        "SLC7A2"        "PDGFRL" 
    ENSG00000129422 ENSG00000078674 ENSG00000104763 ENSG00000156011 ENSG00000104611 
            "MTUS1"          "PCM1"         "ASAH1"          "PSD3"        "SH2D4A" 
    ENSG00000147408 ENSG00000147416 ENSG00000130227 ENSG00000158863 ENSG00000168453 
       "CSGALNACT1"      "ATP6V1B2"          "XPO7"        "FHIP2B"            "HR" 
    ENSG00000168487 ENSG00000168495 ENSG00000104635 ENSG00000120896 ENSG00000120913 
             "BMP1"        "POLR3D"      "SLC39A14"        "SORBS3"        "PDLIM2" 
    ENSG00000179388 ENSG00000008853 ENSG00000120889 ENSG00000173535 ENSG00000173530 
             "EGR3"       "RHOBTB2"     "TNFRSF10B"     "TNFRSF10C"     "TNFRSF10D" 
    ENSG00000104689 ENSG00000147457 ENSG00000104679 ENSG00000134013 ENSG00000197217 
        "TNFRSF10A"         "CHMP7"        "R3HCC1"         "LOXL2"        "ENTPD4" 
    ENSG00000147454 ENSG00000167034 ENSG00000159167 ENSG00000277586 ENSG00000147459 
         "SLC25A37"        "NKX3-1"          "STC1"          "NEFL"         "DOCK5" 
    ENSG00000184661 ENSG00000221914 ENSG00000104765 ENSG00000092964 ENSG00000120899 
            "CDCA2"       "PPP2R2A"        "BNIP3L"        "DPYSL2"         "PTK2B" 
    ENSG00000120885 ENSG00000168077 ENSG00000147419 ENSG00000171320 ENSG00000168078 
              "CLU"        "SCARA3"        "CCDC25"         "ESCO2"           "PBK" 
    ENSG00000186918 ENSG00000104290 ENSG00000012232 ENSG00000197892 ENSG00000120875 
           "ZNF395"          "FZD3"         "EXTL3"        "KIF13B"         "DUSP4" 
    ENSG00000133872 ENSG00000104660 ENSG00000104671 ENSG00000157110 ENSG00000104687 
            "SARAF"      "LEPROTL1"         "DCTN6"         "RBPMS"           "GSR" 
    ENSG00000104695 ENSG00000172733 ENSG00000157168 ENSG00000198042 ENSG00000133874 
           "PPP2CB"          "PURG"          "NRG1"         "MAK16"        "RNF122" 
    ENSG00000183779 ENSG00000147475 ENSG00000147471 ENSG00000020181 ENSG00000156675 
           "ZNF703"        "ERLIN2"         "PLPBP"        "ADGRA2"     "RAB11FIP1" 
    ENSG00000187840 ENSG00000085788 ENSG00000147535 ENSG00000147548 ENSG00000077782 
         "EIF4EBP1"         "DDHD2"         "PLPP5"          "NSD3"         "FGFR1" 
    ENSG00000147526 ENSG00000169499 ENSG00000168615 ENSG00000197140 ENSG00000176907 
            "TACC1"       "PLEKHA2"         "ADAM9"        "ADAM32"          "TCIM" 
    ENSG00000104332 ENSG00000147533 ENSG00000147536 ENSG00000158669 ENSG00000070718 
            "SFRP1"        "GOLGA7"         "GINS4"         "GPAT4"         "AP3M2" 
    ENSG00000104368 ENSG00000104365 ENSG00000070501 ENSG00000078668 ENSG00000168575 
             "PLAT"         "IKBKB"          "POLB"         "VDAC3"       "SLC20A2" 
    ENSG00000176209 ENSG00000168522 ENSG00000185900 ENSG00000165102 ENSG00000164808 
           "SMIM19"          "FNTA"          "POMK"        "HGSNAT"         "SPIDR" 
    ENSG00000221869 ENSG00000253729 ENSG00000104738 ENSG00000019549 ENSG00000168300 
            "CEBPD"         "PRKDC"          "MCM4"         "SNAI2"        "PCMTD1" 
    ENSG00000196711 ENSG00000023287 ENSG00000047249 ENSG00000187735 ENSG00000120992 
           "ALKAL1"        "RB1CC1"       "ATP6V1H"         "TCEA1"        "LYPLA1" 
    ENSG00000167904 ENSG00000137574 ENSG00000008988 ENSG00000181690 ENSG00000170791 
           "TMEM68"          "TGS1"         "RPS20"         "PLAG1"        "CHCHD7" 
    ENSG00000181195 ENSG00000104331 ENSG00000215114 ENSG00000137575 ENSG00000035681 
             "PENK"         "BPNT2"        "UBXN2B"         "SDCBP"         "NSMAF" 
    ENSG00000198846 ENSG00000198363 ENSG00000137563 ENSG00000185728 ENSG00000104442 
              "TOX"          "ASPH"           "GGH"        "YTHDF3"         "ARMC1" 
    ENSG00000066855 ENSG00000205268 ENSG00000147573 ENSG00000147571 ENSG00000185697 
            "MTFR1"         "PDE7A"        "TRIM55"           "CRH"         "MYBL1" 
    ENSG00000104218 ENSG00000137573 ENSG00000140396 ENSG00000067167 ENSG00000178860 
            "CSPP1"         "SULF1"         "NCOA2"         "TRAM1"           "MSC" 
    ENSG00000104321 ENSG00000147601 ENSG00000147604 ENSG00000121039 ENSG00000040341 
            "TRPA1"         "TERF1"          "RPL7"         "RDH10"         "STAU2" 
    ENSG00000104343 ENSG00000175606 ENSG00000104381 ENSG00000091656 ENSG00000164751 
            "UBE2W"        "TMEM70"         "GDAP1"         "ZFHX4"          "PEX2" 
    ENSG00000171033 ENSG00000104427 ENSG00000164683 ENSG00000205189 ENSG00000076641 
             "PKIA"       "ZC2HC1A"          "HEY1"        "ZBTB10"          "PAG1" 
    ENSG00000164687 ENSG00000170323 ENSG00000133731 ENSG00000104231 ENSG00000104497 
            "FABP5"         "FABP4"         "IMPA1"        "ZFAND1"         "SNX16" 
    ENSG00000147614 ENSG00000085719 ENSG00000156103 ENSG00000104312 ENSG00000164823 
         "ATP6V0D2"         "CPNE3"         "MMP16"         "RIPK2"        "OSGIN2" 
    ENSG00000104325 ENSG00000180694 ENSG00000123119 ENSG00000253250 ENSG00000155100 
            "DECR1"        "TMEM64"        "NECAB1"       "C8orf88"        "OTUD6B" 
    ENSG00000188343 ENSG00000164953 ENSG00000164949 ENSG00000197275 ENSG00000164944 
           "CIBAR1"        "TMEM67"           "GEM"        "RAD54B"         "VIRMA" 
    ENSG00000164941 ENSG00000175305 ENSG00000164938 ENSG00000175895 ENSG00000156172 
            "INTS8"         "CCNE2"      "TP53INP1"       "PLEKHF2"       "CFAP418" 
    ENSG00000156466 ENSG00000156471 ENSG00000169439 ENSG00000104324 ENSG00000132561 
             "GDF6"        "PTDSS1"          "SDC2"           "CPQ"         "MATN2" 
    ENSG00000156482 ENSG00000104356 ENSG00000104375 ENSG00000156486 ENSG00000164920 
            "RPL30"          "POP1"          "STK3"         "KCNS2"          "OSR2" 
    ENSG00000132549 ENSG00000164919 ENSG00000156509 ENSG00000147669 ENSG00000104450 
           "VPS13B"         "COX6C"        "FBXO43"        "POLR2K"         "SPAG1" 
    ENSG00000186106 ENSG00000070756 ENSG00000164924 ENSG00000120963 ENSG00000048392 
          "ANKRD46"        "PABPC1"         "YWHAZ"        "ZNF706"         "RRM2B" 
    ENSG00000104517 ENSG00000155090 ENSG00000155096 ENSG00000155097 ENSG00000164929 
             "UBR5"         "KLF10"         "AZIN1"      "ATP6V1C1"         "BAALC" 
    ENSG00000164932 ENSG00000164933 ENSG00000164934 ENSG00000147650 ENSG00000169946 
           "CTHRC1"      "SLC25A32"        "DCAF13"         "LRP12"         "ZFPM2" 
    ENSG00000164830 ENSG00000154188 ENSG00000104408 ENSG00000120526 ENSG00000120533 
             "OXR1"        "ANGPT1"         "EIF3E"        "NUDCD1"          "ENY2" 
    ENSG00000147654 ENSG00000147642 ENSG00000164796 ENSG00000104447 ENSG00000147677 
            "EBAG9"          "SYBU"         "CSMD3"         "TRPS1"         "EIF3H" 
    ENSG00000147679 ENSG00000164754 ENSG00000182197 ENSG00000177570 ENSG00000164761 
            "UTP23"         "RAD21"          "EXT1"        "SAMD12"     "TNFRSF11B" 
    ENSG00000184374 ENSG00000147676 ENSG00000136999 ENSG00000136960 ENSG00000136982 
          "COLEC10"          "MAL2"          "CCN3"         "ENPP2"         "DSCC1" 
    ENSG00000187955 ENSG00000172172 ENSG00000172167 ENSG00000172164 ENSG00000170961 
          "COL14A1"        "MRPL13"          "MTBP"         "SNTB1"          "HAS2" 
    ENSG00000156787 ENSG00000165156 ENSG00000156802 ENSG00000156795 ENSG00000156804 
          "TBC1D31"          "ZHX1"         "ATAD2"         "NTAQ1"        "FBXO32" 
    ENSG00000176853 ENSG00000214814 ENSG00000164983 ENSG00000170881 ENSG00000170873 
          "FAM91A1"        "FER1L6"        "TMEM65"        "RNF139"         "MTSS1" 
    ENSG00000104549 ENSG00000164961 ENSG00000156831 ENSG00000173334 ENSG00000136997 
             "SQLE"        "WASHC5"        "NSMCE2"         "TRIB1"           "MYC" 
    ENSG00000132294 ENSG00000165071 ENSG00000129292 ENSG00000042832 ENSG00000155926 
            "EFR3A"        "TMEM71"       "PHF20L1"            "TG"           "SLA" 
    ENSG00000104415 ENSG00000104419 ENSG00000008513 ENSG00000131773 ENSG00000167632 
             "CCN4"         "NDRG1"       "ST3GAL1"       "KHDRBS3"       "TRAPPC9" 
    ENSG00000104472 ENSG00000123908 ENSG00000169398 ENSG00000105339 ENSG00000204882 
           "CHRAC1"          "AGO2"          "PTK2"        "DENND3"         "GPR20" 
    ENSG00000171045 ENSG00000180155 ENSG00000160932 ENSG00000181638 ENSG00000250571 
          "TSNARE1"         "LYNX1"          "LY6E"         "ZFP41"          "GLI4" 
    ENSG00000014164 ENSG00000204839 ENSG00000147813 ENSG00000104529 ENSG00000179886 
            "ZC3H3"         "MROH6"         "NAPRT"         "EEF1D"         "TIGD5" 
    ENSG00000104524 ENSG00000104522 ENSG00000181135 ENSG00000180921 ENSG00000180900 
            "PYCR3"          "GFUS"        "ZNF707"        "FAM83H"         "SCRIB" 
    ENSG00000261150 ENSG00000178209 ENSG00000178719 ENSG00000178896 ENSG00000179091 
            "EPPK1"          "PLEC"         "GRINA"        "EXOSC4"          "CYC1" 
    ENSG00000179632 ENSG00000235173 ENSG00000179832 ENSG00000261236 ENSG00000185122 
             "MAF1"          "HGH1"         "MROH1"          "BOP1"          "HSF1" 
    ENSG00000185000 ENSG00000185803 ENSG00000071894 ENSG00000160949 ENSG00000187954 
            "DGAT1"       "SLC52A2"         "CPSF1"         "TONSL"       "ZFTRAF1" 
    ENSG00000167702 ENSG00000160957 ENSG00000213563 ENSG00000147799 ENSG00000198169 
            "KIFC2"        "RECQL4"       "C8orf82"      "ARHGAP39"        "ZNF251" 
    ENSG00000161016 ENSG00000197363 ENSG00000147789 ENSG00000170619 ENSG00000196150 
             "RPL8"        "ZNF517"          "ZNF7"        "COMMD5"        "ZNF250" 
    ENSG00000182307 ENSG00000107104 ENSG00000080503 ENSG00000147852 ENSG00000080298 
          "C8orf33"         "KANK1"       "SMARCA2"         "VLDLR"          "RFX3" 
    ENSG00000106688 ENSG00000205808 ENSG00000106993 ENSG00000147853 ENSG00000096968 
           "SLC1A1"         "PLPP6"       "CDC37L1"           "AK3"          "JAK2" 
    ENSG00000120217 ENSG00000197646 ENSG00000107036 ENSG00000099219 ENSG00000183354 
            "CD274"      "PDCD1LG2"          "RIC1"         "ERMP1"         "BRD10" 
    ENSG00000137040 ENSG00000137033 ENSG00000137038 ENSG00000107165 ENSG00000153714 
           "RANBP6"          "IL33"         "DMAC1"         "TYRP1"       "LURAP1L" 
    ENSG00000147862 ENSG00000175893 ENSG00000164975 ENSG00000164985 ENSG00000173068 
             "NFIB"       "ZDHHC21"        "SNAPC3"         "PSIP1"          "BNC2" 
    ENSG00000044459 ENSG00000178031 ENSG00000147874 ENSG00000137145 ENSG00000137154 
            "CNTLN"      "ADAMTSL1"         "HAUS6"       "DENND4C"          "RPS6" 
    ENSG00000177076 ENSG00000188352 ENSG00000188921 ENSG00000198642 ENSG00000184995 
            "ACER2"         "FOCAD"         "HACD4"         "KLHL9"          "IFNE" 
    ENSG00000099810 ENSG00000147883 ENSG00000176399 ENSG00000198680 ENSG00000120159 
             "MTAP"        "CDKN2B"        "DMRTA1"         "TUSC1"         "CAAP1" 
    ENSG00000137055 ENSG00000120156 ENSG00000147894 ENSG00000122729 ENSG00000107201 
             "PLAA"           "TEK"       "C9orf72"          "ACO1"          "RIGI" 
    ENSG00000197579 ENSG00000086061 ENSG00000086062 ENSG00000107262 ENSG00000086065 
           "TOPORS"        "DNAJA1"       "B4GALT1"          "BAG1"         "CHMP5" 
    ENSG00000086102 ENSG00000165272 ENSG00000165271 ENSG00000230453 ENSG00000107341 
             "NFX1"          "AQP3"          "NOL6"      "ANKRD18B"        "UBE2R2" 
    ENSG00000137073 ENSG00000186638 ENSG00000164976 ENSG00000164970 ENSG00000137100 
            "UBAP2"         "KIF24"         "MYORG"       "FAM219A"         "DCTN3" 
    ENSG00000147955 ENSG00000137094 ENSG00000165280 ENSG00000221829 ENSG00000165283 
          "SIGMAR1"        "DNAJB5"           "VCP"         "FANCG"        "STOML2" 
    ENSG00000005238 ENSG00000198722 ENSG00000198853 ENSG00000137135 ENSG00000198467 
            "ATOSB"        "UNC13B"         "RUSC2"      "ARHGEF39"          "TPM2" 
    ENSG00000137076 ENSG00000070610 ENSG00000107185 ENSG00000137133 ENSG00000137103 
             "TLN1"          "GBA2"          "RGP1"         "HINT2"        "TMEM8B" 
    ENSG00000196196 ENSG00000122707 ENSG00000122694 ENSG00000185972 ENSG00000122705 
            "HRCT1"          "RECK"        "GLIPR2"          "CCIN"          "CLTA" 
    ENSG00000137075 ENSG00000165304 ENSG00000137106 ENSG00000137054 ENSG00000147912 
            "RNF38"          "MELK"         "GRHPR"        "POLR1E"        "FBXO10" 
    ENSG00000122696 ENSG00000107338 ENSG00000137124 ENSG00000154330 ENSG00000187866 
         "SLC25A51"           "SHB"       "ALDH1B1"          "PGM5"        "PABIR1" 
    ENSG00000165060 ENSG00000119139 ENSG00000165072 ENSG00000198887 ENSG00000119138 
              "FXN"          "TJP2"        "MAMDC2"          "SMC5"          "KLF9" 
    ENSG00000135048 ENSG00000107362 ENSG00000107372 ENSG00000135046 ENSG00000135045 
           "CEMIP2"       "ABHD17B"        "ZFAND5"         "ANXA1"       "C9orf40" 
    ENSG00000156017 ENSG00000099139 ENSG00000135002 ENSG00000187210 ENSG00000106772 
          "CARNMT1"         "PCSK5"           "RFK"         "GCNT1"        "PRUNE2" 
    ENSG00000197969 ENSG00000156052 ENSG00000148019 ENSG00000135069 ENSG00000106829 
           "VPS13A"          "GNAQ"         "CEP78"         "PSAT1"          "TLE4" 
    ENSG00000148057 ENSG00000135018 ENSG00000165115 ENSG00000165119 ENSG00000178966 
             "IDNK"        "UBQLN1"         "KIF27"        "HNRNPK"          "RMI1" 
    ENSG00000197506 ENSG00000135049 ENSG00000135070 ENSG00000083223 ENSG00000180447 
          "SLC28A3"       "AGTPBP1"         "ISCA1"          "TUT7"          "GAS1" 
    ENSG00000135047 ENSG00000130045 ENSG00000148082 ENSG00000123975 ENSG00000187742 
             "CTSL"         "NXNL2"          "SHC3"          "CKS2"      "SECISBP2" 
    ENSG00000187764 ENSG00000148090 ENSG00000165030 ENSG00000090054 ENSG00000196305 
           "SEMA4D"           "AUH"         "NFIL3"        "SPTLC1"         "IARS1" 
    ENSG00000188312 ENSG00000165233 ENSG00000131669 ENSG00000048828 ENSG00000197724 
            "CENPP"        "CARD19"         "NINJ1"       "FAM120A"          "PHF2" 
    ENSG00000148110 ENSG00000148120 ENSG00000185920 ENSG00000130958 ENSG00000165244 
          "MFSD14B"         "AOPEP"         "PTCH1"       "SLC35D2"        "ZNF367" 
    ENSG00000081377 ENSG00000158122 ENSG00000136943 ENSG00000136842 ENSG00000136925 
           "CDC14B"        "PRXL2C"          "CTSV"         "TMOD1"         "TSTD2" 
    ENSG00000136936 ENSG00000136932 ENSG00000136938 ENSG00000095380 ENSG00000106785 
              "XPA"          "TRMO"        "ANP32B"          "NANS"        "TRIM14" 
    ENSG00000106789 ENSG00000095383 ENSG00000136928 ENSG00000119514 ENSG00000204291 
           "CORO2A"        "TBC1D2"        "GABBR2"       "GALNT12"       "COL15A1" 
    ENSG00000106799 ENSG00000119508 ENSG00000023318 ENSG00000119509 ENSG00000136891 
           "TGFBR1"         "NR4A3"         "ERP44"          "INVS"         "TEX10" 
    ENSG00000066697 ENSG00000170681 ENSG00000136870 ENSG00000165152 ENSG00000155827 
          "MSANTD3"        "CAVIN4"        "ZNF189"         "PGAP4"         "RNF20" 
    ENSG00000136824 ENSG00000136783 ENSG00000165029 ENSG00000070214 ENSG00000106701 
             "SMC2"     "NIPSNAP3A"         "ABCA1"       "SLC44A1"         "FSD1L" 
    ENSG00000106692 ENSG00000095209 ENSG00000119318 ENSG00000136826 ENSG00000070061 
             "FKTN"       "TMEM38B"        "RAD23B"          "KLF4"          "ELP1" 
    ENSG00000119328 ENSG00000119326 ENSG00000106771 ENSG00000070159 ENSG00000243444 
          "ABITRAM"       "CTNNAL1"       "TMEM245"         "PTPN3"              NA 
    ENSG00000136810 ENSG00000165124 ENSG00000136813 ENSG00000106853 ENSG00000059769 
              "TXN"         "SVEP1"         "ECPAS"         "PTGR1"       "DNAJC25" 
    ENSG00000148154 ENSG00000106868 ENSG00000119314 ENSG00000119471 ENSG00000165185 
             "UGCG"         "SUSD1"         "PTBP3"         "HSDL2"      "KIAA1958" 
    ENSG00000148153 ENSG00000119321 ENSG00000136868 ENSG00000148218 ENSG00000148229 
             "INIP"        "FKBP15"       "SLC31A1"          "ALAD"         "POLE3" 
    ENSG00000157657 ENSG00000196739 ENSG00000106948 ENSG00000095397 ENSG00000136888 
           "ZNF618"       "COL27A1"          "AKNA"          "WHRN"      "ATP6V1G1" 
    ENSG00000181634 ENSG00000041982 ENSG00000182752 ENSG00000148219 ENSG00000119401 
          "TNFSF15"           "TNC"         "PAPPA"         "ASTN2"        "TRIM32" 
    ENSG00000136861 ENSG00000106780 ENSG00000095261 ENSG00000119403 ENSG00000056558 
         "CDK5RAP2"         "MEGF9"         "PSMD5"         "PHF19"         "TRAF1" 
    ENSG00000106804 ENSG00000119397 ENSG00000148180 ENSG00000148175 ENSG00000136848 
               "C5"         "CNTRL"           "GSN"          "STOM"        "DAB2IP" 
    ENSG00000175764 ENSG00000119421 ENSG00000119446 ENSG00000095303 ENSG00000197233 
           "TTLL11"        "NDUFA8"         "RBM18"         "PTGS1"         "OR1J2" 
    ENSG00000239590 ENSG00000165202 ENSG00000136940 ENSG00000056586 ENSG00000186130 
            "OR1J4"         "OR1Q1"          "PDCL"         "RC3H2"         "ZBTB6" 
    ENSG00000119408 ENSG00000136930 ENSG00000136942 ENSG00000173611 ENSG00000119414 
             "NEK6"         "PSMB7"         "RPL35"          "SCAI"         "PPP6C" 
    ENSG00000136933 ENSG00000044574 ENSG00000196814 ENSG00000177125 ENSG00000136859 
           "RABEPK"         "HSPA5"        "MVB12B"        "ZBTB34"       "ANGPTL2" 
    ENSG00000136856 ENSG00000196152 ENSG00000197958 ENSG00000136854 ENSG00000160401 
           "SLC2A8"         "ZNF79"         "RPL12"        "STXBP1"       "CFAP157" 
    ENSG00000136807 ENSG00000106991 ENSG00000136840 ENSG00000136908 ENSG00000167106 
             "CDK9"           "ENG"    "ST6GALNAC4"          "DPM2"         "EEIG1" 
    ENSG00000148339 ENSG00000148334 ENSG00000171159 ENSG00000106976 ENSG00000167113 
         "SLC25A25"        "PTGES2"          "BBLN"          "DNM1"          "COQ4" 
    ENSG00000167123 ENSG00000136811 ENSG00000119392 ENSG00000197694 ENSG00000119333 
           "CERCAM"          "ODF2"          "GLE1"        "SPTAN1"       "DYNC2I2" 
    ENSG00000160447 ENSG00000160445 ENSG00000171097 ENSG00000136802 ENSG00000175287 
             "PKN3"          "ZER1"         "KYAT1"        "LRRC8A"        "PHYHD1" 
    ENSG00000095319 ENSG00000148343 ENSG00000167130 ENSG00000095321 ENSG00000119383 
           "NUP188"         "MIGA2"        "DOLPP1"          "CRAT"          "PTPA" 
    ENSG00000188483 ENSG00000148335 ENSG00000148331 ENSG00000167157 ENSG00000136816 
            "IER5L"         "NTMT1"          "ASB6"         "PRRX2"         "TOR1B" 
    ENSG00000136827 ENSG00000136819 ENSG00000148358 ENSG00000107130 ENSG00000130707 
            "TOR1A"       "C9orf78"        "GPR107"          "NCS1"          "ASS1" 
    ENSG00000107164 ENSG00000130713 ENSG00000097007 ENSG00000130720 ENSG00000126883 
            "FUBP3"        "EXOSC2"          "ABL1"        "FIBCD1"        "NUP214" 
    ENSG00000160539 ENSG00000130723 ENSG00000130714 ENSG00000107263 ENSG00000160563 
            "PLPP7"              NA         "POMT1"       "RAPGEF1"         "MED27" 
    ENSG00000196358 ENSG00000107290 ENSG00000125484 ENSG00000165698 ENSG00000148308 
            "NTNG2"          "SETX"        "GTF3C4"        "SPACA9"        "GTF3C5" 
    ENSG00000160271 ENSG00000148297 ENSG00000148303 ENSG00000148290 ENSG00000148300 
           "RALGDS"         "MED22"         "RPL7A"         "SURF1"         "REXO4" 
    ENSG00000160325 ENSG00000123453 ENSG00000160293 ENSG00000169925 ENSG00000196363 
           "CACFD1"         "SARDH"          "VAV2"          "BRD3"          "WDR5" 
    ENSG00000186350 ENSG00000130635 ENSG00000196422 ENSG00000148411 ENSG00000238227 
             "RXRA"        "COL5A1"       "PPP1R26"         "NACC2"       "TMEM250" 
    ENSG00000165661 ENSG00000165684 ENSG00000165688 ENSG00000148384 ENSG00000148400 
            "QSOX2"        "SNAPC4"         "PMPCA"        "INPP5E"        "NOTCH1" 
    ENSG00000172889 ENSG00000169692 ENSG00000196642 ENSG00000054148 ENSG00000107223 
            "EGFL7"        "AGPAT2"         "RABL6"         "PHPT1"          "EDF1" 
    ENSG00000127191 ENSG00000159069 ENSG00000107317 ENSG00000148362 ENSG00000107331 
            "TRAF2"         "FBXW5"         "PTGDS"          "PAXX"         "ABCA2" 
    ENSG00000107281 ENSG00000186193 ENSG00000197355 ENSG00000176978 ENSG00000176248 
            "NPDC1"        "SAPCD2"        "UAP1L1"          "DPP7"        "ANAPC2" 
    ENSG00000176101 ENSG00000187713 ENSG00000188229 ENSG00000188986 ENSG00000165802 
            "SSNA1"       "TMEM203"        "TUBB4B"         "NELFB"          "NSMF" 
    ENSG00000130653 ENSG00000197070 ENSG00000177963 ENSG00000185627 ENSG00000185885 
           "PNPLA7"        "ARRDC1"         "RIC8A"        "PSMD13"        "IFITM1" 
    ENSG00000142089 ENSG00000184363 ENSG00000174915 ENSG00000023191 ENSG00000161328 
           "IFITM3"          "PKP3"        "PTDSS2"          "RNH1"        "LRRC56" 
    ENSG00000099849 ENSG00000070047 ENSG00000185507 ENSG00000177030 ENSG00000177106 
           "RASSF7"         "PHRF1"          "IRF7"         "DEAF1"        "EPS8L2" 
    ENSG00000184524 ENSG00000177542 ENSG00000177595 ENSG00000177600 ENSG00000177666 
            "CEND1"      "SLC25A22"         "PIDD1"         "RPLP2"        "PNPLA2" 
    ENSG00000177700 ENSG00000214063 ENSG00000177830 ENSG00000078902 ENSG00000182208 
           "POLR2L"        "TSPAN4"         "CHID1"        "TOLLIP"          "MOB2" 
    ENSG00000184545 ENSG00000244242 ENSG00000117984 ENSG00000214026 ENSG00000167244 
            "DUSP8"       "IFITM10"          "CTSD"        "MRPL23"          "IGF2" 
    ENSG00000110651 ENSG00000129757 ENSG00000110628 ENSG00000205531 ENSG00000110619 
             "CD81"        "CDKN1C"      "SLC22A18"        "NAP1L4"         "CARS1" 
    ENSG00000005801 ENSG00000148985 ENSG00000177105 ENSG00000167323 ENSG00000167325 
           "ZNF195"         "PGAP2"          "RHOG"         "STIM1"          "RRM1" 
    ENSG00000167333 ENSG00000167332 ENSG00000121236 ENSG00000132256 ENSG00000132274 
           "TRIM68"        "OR51E2"         "TRIM6"         "TRIM5"        "TRIM22" 
    ENSG00000170955 ENSG00000166311 ENSG00000166313 ENSG00000132275 ENSG00000166340 
           "CAVIN3"         "SMPD1"         "APBB1"          "RRP8"          "TPP1" 
    ENSG00000158042 ENSG00000182261 ENSG00000175390 ENSG00000166402 ENSG00000166441 
           "MRPL17"        "NLRP10"         "EIF3F"           "TUB"        "RPL27A" 
    ENSG00000166444 ENSG00000166452 ENSG00000175348 ENSG00000175352 ENSG00000184014 
          "DENND2B"         "AKIP1"        "TMEM9B"         "NRIP3"       "DENND5A" 
    ENSG00000166483 ENSG00000148926 ENSG00000133805 ENSG00000110315 ENSG00000072952 
             "WEE1"           "ADM"         "AMPD3"        "RNF141"         "IRAG1" 
    ENSG00000110321 ENSG00000110328 ENSG00000050165 ENSG00000133816 ENSG00000197702 
           "EIF4G2"       "GALNT18"          "DKK3"        "MICAL2"         "PARVA" 
    ENSG00000187079 ENSG00000133794 ENSG00000148925 ENSG00000197601 ENSG00000133818 
            "TEAD1"         "BMAL1"        "BTBD10"          "FAR1"         "RRAS2" 
    ENSG00000129083 ENSG00000152270 ENSG00000175868 ENSG00000188487 ENSG00000110693 
            "COPB1"         "PDE3B"         "CALCB"          "INSC"          "SOX6" 
    ENSG00000110700 ENSG00000011405 ENSG00000070081 ENSG00000188211 ENSG00000166788 
            "RPS13"       "PIK3C2A"         "NUCB2"       "NCR3LG1"         "SAAL1" 
    ENSG00000110756 ENSG00000110768 ENSG00000134333 ENSG00000074319 ENSG00000151116 
             "HPS5"        "GTF2H1"          "LDHA"        "TSG101"         "UEVLD" 
    ENSG00000151117 ENSG00000129173 ENSG00000166833 ENSG00000187398 ENSG00000109881 
          "TMEM86A"          "E2F8"          "NAV2"         "LUZP2"        "CCDC34" 
    ENSG00000205213 ENSG00000148943 ENSG00000176697 ENSG00000121621 ENSG00000169519 
             "LGR4"         "LIN7C"          "BDNF"        "KIF18A"       "METTL15" 
    ENSG00000109911 ENSG00000049449 ENSG00000060749 ENSG00000176148 ENSG00000280550 
             "ELP4"          "RCN1"         "QSER1"       "TCP11L1"              NA 
    ENSG00000110422 ENSG00000110427 ENSG00000085063 ENSG00000135363 ENSG00000135387 
            "HIPK3"     "KIAA1549L"          "CD59"          "LMO2"       "CAPRIN1" 
    ENSG00000135372 ENSG00000166016 ENSG00000121691 ENSG00000135373 ENSG00000149089 
            "NAT10"         "ABTB2"           "CAT"           "EHF"          "APIP" 
    ENSG00000110435 ENSG00000026508 ENSG00000110436 ENSG00000149090 ENSG00000179431 
             "PDHX"          "CD44"        "SLC1A2"         "PAMR1"          "FJX1" 
    ENSG00000166326 ENSG00000179241 ENSG00000175104 ENSG00000166349 ENSG00000166181 
           "TRIM44"       "LDLRAD3"         "TRAF6"          "RAG1"          "API5" 
    ENSG00000052841 ENSG00000166199 ENSG00000187479 ENSG00000085117 ENSG00000157570 
            "TTC17"        "ALKBH3"      "C11orf96"          "CD82"       "TSPAN18" 
    ENSG00000175274 ENSG00000181830 ENSG00000121671 ENSG00000234776 ENSG00000121680 
          "TP53I11"       "SLC35C1"          "CRY2"         "FREY1"         "PEX16" 
    ENSG00000110492 ENSG00000110497 ENSG00000175224 ENSG00000175220 ENSG00000175216 
              "MDK"        "AMBRA1"         "ATG13"       "ARHGAP1"         "CKAP5" 
    ENSG00000134569 ENSG00000149179 ENSG00000149182 ENSG00000165912 ENSG00000134574 
             "LRP4"        "CSTPP1"       "ARFGAP2"       "PACSIN3"          "DDB2" 
    ENSG00000134575 ENSG00000110514 ENSG00000165915 ENSG00000165916 ENSG00000149187 
             "ACP2"          "MADD"      "SLC39A13"         "PSMC3"         "CELF1" 
    ENSG00000109919 ENSG00000109920 ENSG00000149177 ENSG00000149136 ENSG00000149150 
            "MTCH2"         "FNBP4"         "PTPRJ"         "SSRP1"       "SLC43A1" 
    ENSG00000134809 ENSG00000156587 ENSG00000149131 ENSG00000156599 ENSG00000156603 
           "TIMM10"        "UBE2L6"      "SERPING1"        "ZDHHC5"         "MED19" 
    ENSG00000211450 ENSG00000110031 ENSG00000186660 ENSG00000189057 ENSG00000166801 
          "SELENOH"          "LPXN"         "ZFP91"       "FAM111B"       "FAM111A" 
    ENSG00000110042 ENSG00000110048 ENSG00000166889 ENSG00000166900 ENSG00000166902 
             "DTX4"          "OSBP"         "PATL1"          "STX3"        "MRPL16" 
    ENSG00000183134 ENSG00000110107 ENSG00000110108 ENSG00000006118 ENSG00000167992 
           "PTGDR2"        "PRPF19"       "TMEM109"      "TMEM132A"          "VWCE" 
    ENSG00000167986 ENSG00000149476 ENSG00000162144 ENSG00000134780 ENSG00000124920 
             "DDB1"          "TKFC"      "CYB561A3"         "DAGLA"          "MYRF" 
    ENSG00000134825 ENSG00000168496 ENSG00000134824 ENSG00000149485 ENSG00000221968 
          "TMEM258"          "FEN1"         "FADS2"         "FADS1"         "FADS3" 
    ENSG00000167996 ENSG00000149503 ENSG00000124942 ENSG00000149499 ENSG00000149541 
             "FTH1"        "INCENP"         "AHNAK"          "EML3"        "B3GAT3" 
    ENSG00000185085 ENSG00000204922 ENSG00000162191 ENSG00000214753 ENSG00000162222 
            "INTS5"         "UQCC3"         "UBXN1"      "HNRNPUL2"         "TTC9C" 
    ENSG00000168002 ENSG00000162227 ENSG00000185475 ENSG00000162231 ENSG00000168003 
           "POLR2G"         "TAF6L"      "TMEM179B"          "NXF1"        "SLC3A2" 
    ENSG00000184743 ENSG00000133318 ENSG00000188070 ENSG00000168005 ENSG00000110583 
             "ATL3"          "RTN3"          "ZFTA"       "SPINDOC"         "NAA40" 
    ENSG00000168439 ENSG00000149761 ENSG00000110011 ENSG00000173486 ENSG00000173264 
            "STIP1"        "NUDT22"        "DNAJC4"         "FKBP2"        "GPR137" 
    ENSG00000173153 ENSG00000126432 ENSG00000162302 ENSG00000168066 ENSG00000168067 
            "ESRRA"         "PRDX5"       "RPS6KA4"           "SF1"        "MAP4K2" 
    ENSG00000133895 ENSG00000171219 ENSG00000110047 ENSG00000110046 ENSG00000146670 
             "MEN1"      "CDC42BPG"          "EHD1"         "ATG2A"         "CDCA5" 
    ENSG00000149809 ENSG00000162298 ENSG00000014216 ENSG00000014138 ENSG00000133884 
           "TM7SF2"         "SYVN1"         "CAPN1"         "POLA2"          "DPF2" 
    ENSG00000126391 ENSG00000168056 ENSG00000176973 ENSG00000173442 ENSG00000197136 
            "FRMD8"         "LTBP3"        "FAM89B"       "EHBP1L1"         "PCNX3" 
    ENSG00000173039 ENSG00000172922 ENSG00000172757 ENSG00000172500 ENSG00000175602 
             "RELA"      "RNASEH2C"          "CFL1"          "FIBP"       "CCDC85B" 
    ENSG00000175592 ENSG00000175573 ENSG00000175550 ENSG00000175467 ENSG00000175376 
            "FOSL1"      "C11orf68"         "DRAP1"         "SART1"        "EIF1AD" 
    ENSG00000175334 ENSG00000087365 ENSG00000175115 ENSG00000174996 ENSG00000174807 
            "BANF1"         "SF3B2"         "PACS1"          "KLC2"         "CD248" 
    ENSG00000174791 ENSG00000174576 ENSG00000174547 ENSG00000174080 ENSG00000239306 
             "RIN1"         "NPAS4"        "MRPL11"          "CTSF"         "RBM14" 
    ENSG00000173914 ENSG00000173715 ENSG00000173599 ENSG00000173227 ENSG00000173020 
            "RBM4B"        "TOP6BL"            "PC"         "SYT12"          "GRK2" 
    ENSG00000172932 ENSG00000175482 ENSG00000175505 ENSG00000172531 ENSG00000175634 
         "ANKRD13D"         "POLD4"         "CLCF1"        "PPP1CA"       "RPS6KB2" 
    ENSG00000172725 ENSG00000172663 ENSG00000110697 ENSG00000167797 ENSG00000084207 
           "CORO1B"       "TMEM134"       "PITPNM1"       "CDK2AP2"         "GSTP1" 
    ENSG00000110057 ENSG00000006534 ENSG00000110717 ENSG00000110719 ENSG00000110721 
          "UNC93B1"       "ALDH3B1"        "NDUFS8"        "TCIRG1"          "CHKA" 
    ENSG00000110066 ENSG00000171067 ENSG00000110075 ENSG00000069482 ENSG00000132749 
            "KMT5B"      "C11orf24"        "PPP6R3"           "GAL"        "TESMIN" 
    ENSG00000110090 ENSG00000172935 ENSG00000162341 ENSG00000172927 ENSG00000110092 
            "CPT1A"        "MRGPRF"         "TPCN2"         "MYEOV"         "CCND1" 
    ENSG00000131620 ENSG00000168040 ENSG00000085733 ENSG00000162105 ENSG00000172893 
             "ANO1"          "FADD"          "CTTN"        "SHANK2"         "DHCR7" 
    ENSG00000158483 ENSG00000137497 ENSG00000149357 ENSG00000165458 ENSG00000162129 
         "FAM86C1P"         "NUMA1"       "LAMTOR1"        "INPPL1"          "CLPB" 
    ENSG00000186635 ENSG00000214530 ENSG00000110237 ENSG00000054967 ENSG00000054965 
            "ARAP1"       "STARD10"      "ARHGEF17"          "RELT"       "FAM168A" 
    ENSG00000175582 ENSG00000181924 ENSG00000175567 ENSG00000168014 ENSG00000149380 
            "RAB6A"          "COA4"          "UCP2"         "C2CD3"         "P4HA3" 
    ENSG00000165434 ENSG00000077514 ENSG00000166439 ENSG00000166435 ENSG00000162139 
           "PGM2L1"         "POLD3"        "RNF169"         "XRRA1"          "NEU3" 
    ENSG00000137486 ENSG00000149273 ENSG00000158555 ENSG00000149257 ENSG00000198382 
            "ARRB1"          "RPS3"         "GDPD5"      "SERPINH1"         "UVRAG" 
    ENSG00000137492 ENSG00000182704 ENSG00000078124 ENSG00000149260 ENSG00000149269 
           "THAP12"          "TSKU"         "ACER3"         "CAPN5"          "PAK1" 
    ENSG00000074201 ENSG00000048649 ENSG00000087884 ENSG00000149262 ENSG00000188997 
           "CLNS1A"          "RSF1"         "AAMDC"         "INTS4"        "KCTD21" 
    ENSG00000033327 ENSG00000149256 ENSG00000165490 ENSG00000137502 ENSG00000137494 
             "GAB2"         "TENM4"         "DDIAS"         "RAB30"       "ANKRD42" 
    ENSG00000137500 ENSG00000171204 ENSG00000171202 ENSG00000137504 ENSG00000179071 
          "CCDC90B"      "TMEM126B"      "TMEM126A"        "CREBZF"        "CCDC89" 
    ENSG00000137501 ENSG00000073921 ENSG00000074266 ENSG00000149201 ENSG00000150687 
            "SYTL2"        "PICALM"           "EED"        "CCDC81"        "PRSS23" 
    ENSG00000123892 ENSG00000109861 ENSG00000086991 ENSG00000110172 ENSG00000180773 
            "RAB38"          "CTSC"          "NOX4"       "CHORDC1"       "SLC36A4" 
    ENSG00000166002 ENSG00000166012 ENSG00000042429 ENSG00000214376 ENSG00000110218 
            "SMCO4"         "TAF1D"         "MED17"         "VSTM5"         "PANX1" 
    ENSG00000020922 ENSG00000263465 ENSG00000149218 ENSG00000149212 ENSG00000077458 
            "MRE11"         "SRSF8"        "ENDOD1"         "SESN3"        "FAM76B" 
    ENSG00000087053 ENSG00000149231 ENSG00000165895 ENSG00000170647 ENSG00000137672 
            "MTMR2"        "CCDC82"      "ARHGAP42"              NA         "TRPC6" 
    ENSG00000023445 ENSG00000110330 ENSG00000152558 ENSG00000166670 ENSG00000196611 
            "BIRC3"         "BIRC2"       "TMEM123"         "MMP10"          "MMP1" 
    ENSG00000149968 ENSG00000137745 ENSG00000137692 ENSG00000187240 ENSG00000170962 
             "MMP3"         "MMP13"       "DCUN1D5"       "DYNC2H1"         "PDGFD" 
    ENSG00000196954 ENSG00000137752 ENSG00000204397 ENSG00000152578 ENSG00000170903 
            "CASP4"         "CASP1"        "CARD16"         "GRIA4"       "MSANTD4" 
    ENSG00000149313 ENSG00000152402 ENSG00000137760 ENSG00000075239 ENSG00000149308 
         "AASDHPPT"       "GUCY1A2"        "ALKBH8"         "ACAT1"          "NPAT" 
    ENSG00000149311 ENSG00000178202 ENSG00000178105 ENSG00000149289 ENSG00000137710 
              "ATM"       "POGLUT3"         "DDX10"       "ZC3H12C"           "RDX" 
    ENSG00000137714 ENSG00000137727 ENSG00000204381 ENSG00000170145 ENSG00000137713 
             "FDX1"      "ARHGAP20"          "LAYN"          "SIK2"       "PPP2R1B" 
    ENSG00000109846 ENSG00000150764 ENSG00000150768 ENSG00000150779 ENSG00000150782 
            "CRYAB"        "DIXDC1"          "DLAT"        "TIMM8B"          "IL18" 
    ENSG00000150787 ENSG00000149292 ENSG00000048028 ENSG00000166741 ENSG00000076053 
              "PTS"         "TTC12"         "USP28"          "NNMT"          "RBM7" 
    ENSG00000076043 ENSG00000182985 ENSG00000160584 ENSG00000168092 ENSG00000149591 
            "REXO2"         "CADM1"          "SIK3"      "PAFAH1B2"         "TAGLN" 
    ENSG00000186318 ENSG00000110274 ENSG00000177098 ENSG00000167283 ENSG00000118058 
            "BACE1"        "CEP164"         "SCN4B"        "ATP5MG"         "KMT2A" 
    ENSG00000118096 ENSG00000019144 ENSG00000110367 ENSG00000118181 ENSG00000137700 
            "IFT46"        "PHLDB1"          "DDX6"         "RPS25"       "SLC37A4" 
    ENSG00000149428 ENSG00000256269 ENSG00000188486 ENSG00000172269 ENSG00000172375 
            "HYOU1"          "HMBS"          "H2AX"        "DPAGT1"        "C2CD2L" 
    ENSG00000172273 ENSG00000160703 ENSG00000110395 ENSG00000076706 ENSG00000173456 
            "HINFP"         "NLRX1"           "CBL"          "MCAM"         "RNF26" 
    ENSG00000154096 ENSG00000110400 ENSG00000184232 ENSG00000181264 ENSG00000196914 
             "THY1"       "NECTIN1"           "OAF"         "TLCD5"      "ARHGEF12" 
    ENSG00000149403 ENSG00000154114 ENSG00000109929 ENSG00000137642 ENSG00000154127 
            "GRIK4"         "TBCEL"          "SC5D"         "SORL1"       "UBASH3B" 
    ENSG00000109944 ENSG00000109971 ENSG00000023171 ENSG00000110002 ENSG00000154144 
              "JHY"         "HSPA8"       "GRAMD1B"         "VWA5A"         "TBRG1" 
    ENSG00000110013 ENSG00000154146 ENSG00000120458 ENSG00000154134 ENSG00000154133 
             "SIAE"          "NRGN"       "MSANTD2"         "ROBO3"         "ROBO4" 
    ENSG00000149548 ENSG00000134955 ENSG00000165495 ENSG00000149547 ENSG00000134910 
           "CCDC15"       "SLC37A2"        "PKNOX2"          "EI24"         "STT3A" 
    ENSG00000149554 ENSG00000198331 ENSG00000110060 ENSG00000064309 ENSG00000197798 
            "CHEK1"         "HYLS1"          "PUS3"          "CDON"       "FAM118B" 
    ENSG00000182934 ENSG00000110074 ENSG00000110063 ENSG00000149571 ENSG00000134954 
            "SRPRA"       "FOXRED1"          "DCPS"       "KIRREL3"          "ETS1" 
    ENSG00000151702 ENSG00000120457 ENSG00000174370 ENSG00000134909 ENSG00000084234 
             "FLI1"         "KCNJ5"     "KCNJ5-AS1"      "ARHGAP32"         "APLP2" 
    ENSG00000196323 ENSG00000134917 ENSG00000120451 ENSG00000182667 ENSG00000166086 
           "ZBTB44"       "ADAMTS8"         "SNX19"           "NTM"          "JAM3" 
    ENSG00000151503 ENSG00000151502 ENSG00000151498 ENSG00000151240 ENSG00000180525 
           "NCAPD3"        "VPS26B"         "ACAD8"         "DIP2C"     "DIP2C-AS1" 
    ENSG00000107929 ENSG00000107937 ENSG00000067064 ENSG00000047056 ENSG00000067057 
           "LARP4B"        "GTPBP4"          "IDI1"         "WDR37"          "PFKP" 
    ENSG00000107959 ENSG00000067082 ENSG00000151632 ENSG00000187134 ENSG00000196139 
           "PITRM1"          "KLF6"        "AKR1C2"        "AKR1C1"        "AKR1C3" 
    ENSG00000173848 ENSG00000196372 ENSG00000108021 ENSG00000057608 ENSG00000134452 
             "NET1"         "ASB13"        "TASOR2"          "GDI2"          "FBH1" 
    ENSG00000134470 ENSG00000170525 ENSG00000198879 ENSG00000165629 ENSG00000107485 
           "IL15RA"        "PFKFB3"        "SFMBT2"       "ATP5F1C"         "GATA3" 
    ENSG00000048740 ENSG00000148429 ENSG00000134463 ENSG00000148426 ENSG00000151465 
            "CELF2"        "USP6NL"        "ECHDC3"       "PROSER2"        "CDC123" 
    ENSG00000183049 ENSG00000123240 ENSG00000065328 ENSG00000107537 ENSG00000086475 
           "CAMK1D"          "OPTN"         "MCM10"          "PHYH"        "SEPHS1" 
    ENSG00000165626 ENSG00000151474 ENSG00000065809 ENSG00000187522 ENSG00000152455 
            "BEND7"        "FRMD4A"       "FAM107B"        "HSPA14"       "SUV39H2" 
    ENSG00000152464 ENSG00000152465 ENSG00000077943 ENSG00000148481 ENSG00000148484 
            "RPP38"          "NMT2"         "ITGA8"        "MINDY3"          "RSU1" 
    ENSG00000026025 ENSG00000165996 ENSG00000136738 ENSG00000165997 ENSG00000120594 
              "VIM"         "HACD1"          "STAM"         "ARL5B"        "PLXDC2" 
    ENSG00000204682 ENSG00000180592 ENSG00000078403 ENSG00000148444 ENSG00000150867 
        "MIR1915HG"        "SKIDA1"        "MLLT10"        "COMMD3"       "PIP4K2A" 
    ENSG00000148450 ENSG00000165312 ENSG00000107863 ENSG00000099256 ENSG00000148459 
            "MSRB2"         "OTUD1"      "ARHGAP21"       "PRTFDC1"         "PDSS1" 
    ENSG00000107890 ENSG00000120539 ENSG00000107897 ENSG00000099246 ENSG00000169126 
          "ANKRD26"         "MASTL"         "ACBD5"         "RAB18"         "ODAD2" 
    ENSG00000150054 ENSG00000095739 ENSG00000197321 ENSG00000165757 ENSG00000107968 
             "MPP7"         "BAMBI"          "SVIL"          "JCAD"        "MAP3K8" 
    ENSG00000183621 ENSG00000148516 ENSG00000170759 ENSG00000150093 ENSG00000099250 
           "ZNF438"          "ZEB1"         "KIF5B"         "ITGB1"          "NRP1" 
    ENSG00000148498 ENSG00000108094 ENSG00000095794 ENSG00000108100 ENSG00000177283 
            "PARD3"          "CUL2"          "CREM"          "CCNY"          "FZD8" 
    ENSG00000198105 ENSG00000189180 ENSG00000196693 ENSG00000165733 ENSG00000165731 
           "ZNF248"        "ZNF33A"        "ZNF33B"          "BMS1"           "RET" 
    ENSG00000169826 ENSG00000169813 ENSG00000196793 ENSG00000169740 ENSG00000107562 
       "CSGALNACT2"        "HNRNPF"        "ZNF239"         "ZNF32"        "CXCL12" 
    ENSG00000107551 ENSG00000165507 ENSG00000165512 ENSG00000165406 ENSG00000172661 
           "RASSF4"         "DEPP1"         "ZNF22"       "MARCHF8"       "WASHC2C" 
    ENSG00000266412 ENSG00000107643 ENSG00000128805 ENSG00000165633 ENSG00000099290 
            "NCOA4"         "MAPK8"      "ARHGAP22"         "VSTM4"       "WASHC2A" 
    ENSG00000188611 ENSG00000198964 ENSG00000204147 ENSG00000185532 ENSG00000107984 
            "ASAH2"         "SGMS1"        "ASAH2B"         "PRKG1"          "DKK1" 
    ENSG00000122952 ENSG00000151151 ENSG00000072401 ENSG00000108064 ENSG00000122870 
            "ZWINT"          "IPMK"        "UBE2D1"          "TFAM"         "BICC1" 
    ENSG00000108091 ENSG00000151150 ENSG00000170312 ENSG00000196932 ENSG00000150347 
            "CCDC6"          "ANK3"          "CDK1"        "TMEM26"        "ARID5B" 
    ENSG00000182010 ENSG00000138311 ENSG00000181915 ENSG00000122877 ENSG00000171988 
            "RTKN2"        "ZNF365"           "ADO"          "EGR2"        "JMJD1C" 
    ENSG00000108176 ENSG00000096717 ENSG00000148634 ENSG00000138347 ENSG00000108187 
          "DNAJC12"         "SIRT1"         "HERC4"          "MYPN"          "PBLD" 
    ENSG00000096746 ENSG00000204130 ENSG00000122912 ENSG00000060339 ENSG00000107625 
          "HNRNPH3"         "RUFY2"      "SLC25A16"         "CCAR1"         "DDX50" 
    ENSG00000165732 ENSG00000198954 ENSG00000122862 ENSG00000122958 ENSG00000156502 
            "DDX21"         "KIFBP"          "SRGN"        "VPS26A"       "SUPV3L1" 
    ENSG00000156515 ENSG00000197467 ENSG00000042286 ENSG00000079332 ENSG00000180817 
              "HK1"       "COL13A1"         "AIFM2"         "SAR1A"          "PPA1" 
    ENSG00000172731 ENSG00000148730 ENSG00000166228 ENSG00000107731 ENSG00000198246 
           "LRRC20"      "EIF4EBP2"         "PCBD1"         "UNC5B"       "SLC29A3" 
    ENSG00000107738 ENSG00000197746 ENSG00000122863 ENSG00000138303 ENSG00000166295 
             "VSIR"          "PSAP"         "CHST3"         "ASCC1"       "ANAPC16" 
    ENSG00000107745 ENSG00000122884 ENSG00000138286 ENSG00000213551 ENSG00000182180 
            "MICU1"         "P4HA1"      "FAM149B1"        "DNAJC9"        "MRPS16" 
    ENSG00000138279 ENSG00000166343 ENSG00000107758 ENSG00000166317 ENSG00000176986 
            "ANXA7"         "MSS51"        "PPP3CB"       "SYNPO2L"        "SEC24C" 
    ENSG00000196968 ENSG00000214655 ENSG00000122861 ENSG00000222047 ENSG00000035403 
            "FUT11"        "ZSWIM8"          "PLAU"      "C10orf55"           "VCL" 
    ENSG00000156110 ENSG00000156650 ENSG00000156671 ENSG00000165655 ENSG00000151208 
              "ADK"         "KAT6B"         "SAMD8"        "ZNF503"          "DLG5" 
    ENSG00000148606 ENSG00000138326 ENSG00000253626 ENSG00000189129 ENSG00000122359 
           "POLR3A"         "RPS24"       "EIF5AL1"         "PLAC9"        "ANXA11" 
    ENSG00000122378 ENSG00000108219 ENSG00000107771 ENSG00000122367 ENSG00000107779 
           "PRXL2A"       "TSPAN14"        "CCSER2"          "LDB3"        "BMPR1A" 
    ENSG00000173267 ENSG00000122376 ENSG00000107789 ENSG00000198682 ENSG00000138138 
             "SNCG"         "SHLD2"        "MINPP1"        "PAPSS2"         "ATAD1" 
    ENSG00000171862 ENSG00000138134 ENSG00000107796 ENSG00000026103 ENSG00000119922 
             "PTEN"      "STAMBPL1"         "ACTA2"           "FAS"         "IFIT2" 
    ENSG00000185745 ENSG00000152778 ENSG00000152779 ENSG00000152782 ENSG00000138182 
            "IFIT1"         "IFIT5"      "SLC16A12"         "PANK1"        "KIF20B" 
    ENSG00000148680 ENSG00000148677 ENSG00000180628 ENSG00000165338 ENSG00000119938 
             "HTR7"        "ANKRD1"         "PCGF5"        "HECTD2"       "PPP1R3C" 
    ENSG00000107854 ENSG00000095564 ENSG00000198060 ENSG00000119912 ENSG00000138160 
            "TNKS2"         "BTAF1"       "MARCHF5"           "IDE"         "KIF11" 
    ENSG00000152804 ENSG00000138190 ENSG00000138119 ENSG00000138180 ENSG00000148690 
             "HHEX"         "EXOC6"          "MYOF"         "CEP55"      "FRA10AC1" 
    ENSG00000138193 ENSG00000173145 ENSG00000119969 ENSG00000107438 ENSG00000059573 
            "PLCE1"         "NOC3L"         "HELLS"        "PDLIM1"      "ALDH18A1" 
    ENSG00000119977 ENSG00000138185 ENSG00000177853 ENSG00000077147 ENSG00000196233 
            "TCTN3"        "ENTPD1"       "ZNF518A"        "TM9SF3"          "LCOR" 
    ENSG00000155640 ENSG00000213390 ENSG00000165879 ENSG00000181274 ENSG00000052749 
                 NA      "ARHGAP19"         "FRAT1"         "FRAT2"         "RRP12" 
    ENSG00000171314 ENSG00000171307 ENSG00000155229 ENSG00000165886 ENSG00000171160 
            "PGAM1"       "ZDHHC16"         "MMS19"         "UBTD1"         "MORN4" 
    ENSG00000155252 ENSG00000119986 ENSG00000155254 ENSG00000155256 ENSG00000166024 
           "PI4K2A"         "AVPI1"      "MARVELD1"       "ZFYVE27"       "R3HCC1L" 
    ENSG00000138131 ENSG00000119943 ENSG00000107521 ENSG00000119946 ENSG00000120053 
            "LOXL4"       "PYROXD2"          "HPS1"         "CNNM1"          "GOT1" 
    ENSG00000155287 ENSG00000198018 ENSG00000014919 ENSG00000107554 ENSG00000107566 
         "SLC25A28"        "ENTPD7"         "COX15"         "DNMBP"        "ERLIN1" 
    ENSG00000213341 ENSG00000196072 ENSG00000099194 ENSG00000166135 ENSG00000119906 
             "CHUK"       "BLOC1S2"           "SCD"        "HIF1AN"          "SLF2" 
    ENSG00000055950 ENSG00000107816 ENSG00000107819 ENSG00000107821 ENSG00000166167 
           "MRPL43"         "LZTS2"         "SFXN3"       "KAZALD1"          "BTRC" 
    ENSG00000107829 ENSG00000107833 ENSG00000198408 ENSG00000120029 ENSG00000166189 
            "FBXW4"          "NPM3"           "OGA"         "ARMH3"          "HPS6" 
    ENSG00000198728 ENSG00000148840 ENSG00000166197 ENSG00000119915 ENSG00000107862 
             "LDB1"         "PPRC1"         "NOLC1"        "ELOVL3"          "GBF1" 
    ENSG00000077150 ENSG00000107874 ENSG00000171206 ENSG00000138175 ENSG00000156398 
            "NFKB2"        "CUEDC2"         "TRIM8"          "ARL3"         "SFXN2" 
    ENSG00000166275 ENSG00000148842 ENSG00000076685 ENSG00000148798 ENSG00000156374 
           "BORCS7"         "CNNM2"         "NT5C2"           "INA"         "PCGF6" 
    ENSG00000148843 ENSG00000107957 ENSG00000107960 ENSG00000065618 ENSG00000148834 
           "PDCD11"      "SH3PXD2A"          "STN1"       "COL17A1"         "GSTO1" 
    ENSG00000065621 ENSG00000148841 ENSG00000148700 ENSG00000119950 ENSG00000138166 
            "GSTO2"        "ITPRIP"          "ADD3"          "MXI1"         "DUSP5" 
    ENSG00000108055 ENSG00000108061 ENSG00000119927 ENSG00000197142 ENSG00000165806 
             "SMC3"         "SHOC2"          "GPAM"         "ACSL5"         "CASP7" 
    ENSG00000198924 ENSG00000169129 ENSG00000099204 ENSG00000151553 ENSG00000107518 
          "DCLRE1A"       "AFAP1L2"        "ABLIM1"        "FHIP2A"        "ATRNL1" 
    ENSG00000151892 ENSG00000165868 ENSG00000187164 ENSG00000165650 ENSG00000151893 
            "GFRA1"       "HSPA12A"         "SHTN1"         "PDZD8"        "CACUL1" 
    ENSG00000107581 ENSG00000119979 ENSG00000183605 ENSG00000165672 ENSG00000198873 
            "EIF3A"       "DENND10"         "SFXN4"         "PRDX3"          "GRK5" 
    ENSG00000148908 ENSG00000151929 ENSG00000197771 ENSG00000120008 ENSG00000138162 
            "RGS10"          "BAG3"         "MCMBP"         "WDR11"         "TACC2" 
    ENSG00000107679 ENSG00000166033 ENSG00000154473 ENSG00000189319 ENSG00000165660 
          "PLEKHA1"         "HTRA1"          "BUB3"        "FAM53B"      "ABRAXAS2" 
    ENSG00000019995 ENSG00000175029 ENSG00000188690 ENSG00000107949 ENSG00000089876 
           "ZRANB1"         "CTBP2"          "UROS"         "BCCIP"         "DHX32" 
    ENSG00000203780 ENSG00000148848 ENSG00000150760 ENSG00000132334 ENSG00000148773 
            "FANK1"        "ADAM12"         "DOCK1"         "PTPRE"         "MKI67" 
    ENSG00000151640 ENSG00000148814 ENSG00000130640 ENSG00000203772 ENSG00000120645 
           "DPYSL4"        "LRRC27"       "TUBGCP2"          "SPRN"        "IQSEC3" 
    ENSG00000073614 ENSG00000120647 ENSG00000060237 ENSG00000082805 ENSG00000111186 
            "KDM5A"        "CCDC77"          "WNK1"          "ERC1"         "WNT5B" 
    ENSG00000171823 ENSG00000006831 ENSG00000151062 ENSG00000151067 ENSG00000004478 
           "FBXL14"       "ADIPOR2"      "CACNA2D4"       "CACNA1C"         "FKBP4" 
    ENSG00000111203 ENSG00000111206 ENSG00000171792 ENSG00000078246 ENSG00000197905 
            "ITFG2"         "FOXM1"         "RHNO1"         "TULP3"         "TEAD4" 
    ENSG00000011105 ENSG00000111224 ENSG00000078237 ENSG00000047621 ENSG00000111247 
           "TSPAN9"        "PARP11"         "TIGAR"        "FERRY3"      "RAD51AP1" 
    ENSG00000185652 ENSG00000010278 ENSG00000067182 ENSG00000111321 ENSG00000111639 
             "NTF3"           "CD9"      "TNFRSF1A"          "LTBR"        "MRPL51" 
    ENSG00000010292 ENSG00000111640 ENSG00000010295 ENSG00000111641 ENSG00000111642 
           "NCAPD2"         "GAPDH"         "IFFO1"          "NOP2"          "CHD4" 
    ENSG00000111653 ENSG00000089693 ENSG00000159335 ENSG00000089692 ENSG00000250510 
             "ING4"          "MLF2"          "PTMS"          "LAG3"        "GPR162" 
    ENSG00000111665 ENSG00000111667 ENSG00000111669 ENSG00000111674 ENSG00000111678 
            "CDCA3"          "USP5"          "TPI1"          "ENO2"      "C12orf57" 
    ENSG00000126749 ENSG00000111684 ENSG00000182326 ENSG00000139182 ENSG00000139197 
             "EMG1"        "LPCAT3"           "C1S"        "CLSTN3"          "PEX5" 
    ENSG00000177675 ENSG00000173262 ENSG00000059804 ENSG00000065970 ENSG00000166532 
          "CD163L1"       "SLC2A14"        "SLC2A3"         "FOXJ2"        "RIMKLB" 
    ENSG00000111752 ENSG00000003056 ENSG00000175899 ENSG00000069493 ENSG00000110848 
             "PHC1"          "M6PR"           "A2M"        "CLEC2D"          "CD69" 
    ENSG00000110852 ENSG00000173391 ENSG00000139112 ENSG00000111196 ENSG00000060140 
           "CLEC2B"          "OLR1"     "GABARAPL1"        "MAGOHB"         "STYK1" 
    ENSG00000060138 ENSG00000139083 ENSG00000070018 ENSG00000111266 ENSG00000111269 
             "YBX3"          "ETV6"          "LRP6"        "DUSP16"        "CREBL2" 
    ENSG00000111276 ENSG00000178878 ENSG00000013588 ENSG00000013583 ENSG00000111305 
           "CDKN1B"        "APOLD1"        "GPRC5A"         "HEBP1"          "GSG1" 
    ENSG00000246705 ENSG00000111348 ENSG00000151491 ENSG00000023734 ENSG00000023697 
             "H2AJ"       "ARHGDIB"          "EPS8"         "STRAP"          "DERA" 
    ENSG00000008394 ENSG00000052126 ENSG00000139154 ENSG00000172572 ENSG00000111711 
            "MGST1"       "PLEKHA5"         "AEBP2"         "PDE3A"        "GOLT1B" 
    ENSG00000111716 ENSG00000121361 ENSG00000069431 ENSG00000111726 ENSG00000111728 
             "LDHB"         "KCNJ8"         "ABCC9"          "CMAS"       "ST8SIA1" 
    ENSG00000111731 ENSG00000139163 ENSG00000134532 ENSG00000060982 ENSG00000133703 
            "C2CD5"         "ETNK1"          "SOX5"         "BCAT1"          "KRAS" 
    ENSG00000152936 ENSG00000123094 ENSG00000123096 ENSG00000064102 ENSG00000064115 
           "LMNTD1"        "RASSF8"          "SSPN"        "INTS13"        "TM7SF3" 
    ENSG00000152944 ENSG00000211455 ENSG00000029153 ENSG00000110841 ENSG00000061794 
            "MED21"        "STK38L"         "BMAL2"       "PPFIBP1"        "MRPS35" 
    ENSG00000087448 ENSG00000087494 ENSG00000123106 ENSG00000087502 ENSG00000187950 
           "KLHL42"         "PTHLH"        "CCDC91"        "ERGIC2"         "OVCH1" 
    ENSG00000133687 ENSG00000110888 ENSG00000110900 ENSG00000013573 ENSG00000139146 
            "TMTC1"       "CAPRIN2"       "TSPAN11"         "DDX11"       "SINHCAF" 
    ENSG00000170456 ENSG00000151743 ENSG00000174718 ENSG00000139132 ENSG00000057294 
          "DENND5B"          "AMN1"         "RESF1"          "FGD4"          "PKP2" 
    ENSG00000175548 ENSG00000151229 ENSG00000018236 ENSG00000151233 ENSG00000015153 
           "ALG10B"       "SLC2A13"         "CNTN1"        "GXYLT1"          "YAF2" 
    ENSG00000139174 ENSG00000129317 ENSG00000198001 ENSG00000151239 ENSG00000139173 
         "PRICKLE1"         "PUS7L"         "IRAK4"          "TWF1"       "TMEM117" 
    ENSG00000177119 ENSG00000189079 ENSG00000111371 ENSG00000134294 ENSG00000139211 
             "ANO6"         "ARID2"       "SLC38A1"       "SLC38A2"        "AMIGO2" 
    ENSG00000079337 ENSG00000211584 ENSG00000061273 ENSG00000111424 ENSG00000134291 
          "RAPGEF3"       "SLC48A1"         "HDAC7"           "VDR"      "TMEM106C" 
    ENSG00000079387 ENSG00000152556 ENSG00000177981 ENSG00000167528 ENSG00000139620 
            "SENP1"          "PFKM"          "ASB8"        "ZNF641"        "KANSL2" 
    ENSG00000129315 ENSG00000174233 ENSG00000167535 ENSG00000174243 ENSG00000167548 
            "CCNT1"         "ADCY6"        "CACNB3"         "DDX23"         "KMT2D" 
    ENSG00000139636 ENSG00000123416 ENSG00000167552 ENSG00000167553 ENSG00000135451 
           "LMBR1L"        "TUBA1B"        "TUBA1A"        "TUBA1C"         "TROAP" 
    ENSG00000186897 ENSG00000178401 ENSG00000123352 ENSG00000110844 ENSG00000161791 
            "C1QL4"       "DNAJC22"        "SPATS2"       "PRPF40B"         "FMNL3" 
    ENSG00000139644 ENSG00000167566 ENSG00000161800 ENSG00000110881 ENSG00000178449 
           "TMBIM6"       "NCKAP5L"       "RACGAP1"         "ASIC1"         "COX14" 
    ENSG00000139624 ENSG00000050405 ENSG00000161813 ENSG00000066084 ENSG00000123268 
            "CERS5"         "LIMA1"         "LARP4"         "DIP2B"          "ATF1" 
    ENSG00000185432 ENSG00000110911 ENSG00000050426 ENSG00000110925 ENSG00000135457 
            "TMT1A"       "SLC11A2"        "LETMD1"        "CSRNP2"         "TFCP2" 
    ENSG00000184271 ENSG00000183283 ENSG00000170545 ENSG00000139629 ENSG00000050438 
           "POU6F1"        "DAZAP2"         "SMAGP"        "GALNT6"        "SLC4A8" 
    ENSG00000196876 ENSG00000139567 ENSG00000135503 ENSG00000161835 ENSG00000123358 
            "SCN8A"        "ACVRL1"        "ACVR1B"       "TAMALIN"         "NR4A1" 
    ENSG00000123395 ENSG00000167767 ENSG00000186442 ENSG00000170421 ENSG00000111057 
           "ATG101"         "KRT80"          "KRT3"          "KRT8"         "KRT18" 
    ENSG00000063046 ENSG00000111077 ENSG00000167778 ENSG00000139651 ENSG00000182544 
            "EIF4B"          "TNS2"        "SPRYD3"        "ZNF740"         "MFSD5" 
    ENSG00000135476 ENSG00000123349 ENSG00000139637 ENSG00000094914 ENSG00000185591 
            "ESPL1"         "PFDN5"          "MYG1"          "AAAS"           "SP1" 
    ENSG00000139625 ENSG00000170653 ENSG00000012822 ENSG00000037965 ENSG00000094916 
          "MAP3K12"          "ATF7"      "CALCOCO1"         "HOXC8"          "CBX5" 
    ENSG00000135486 ENSG00000111481 ENSG00000161638 ENSG00000170439 ENSG00000135424 
          "HNRNPA1"         "COPZ1"         "ITGA5"         "TMT1B"         "ITGA7" 
    ENSG00000135437 ENSG00000135404 ENSG00000135414 ENSG00000170473 ENSG00000065357 
             "RDH5"          "CD63"         "GDF11"          "PYM1"          "DGKA" 
    ENSG00000185664 ENSG00000123374 ENSG00000197728 ENSG00000170515 ENSG00000229117 
             "PMEL"          "CDK2"         "RPS26"         "PA2G4"         "RPL41" 
    ENSG00000135482 ENSG00000139641 ENSG00000196465 ENSG00000092841 ENSG00000139645 
           "ZC3H10"         "ESYT1"         "MYL6B"          "MYL6"       "ANKRD52" 
    ENSG00000135469 ENSG00000062485 ENSG00000257727 ENSG00000170581 ENSG00000111602 
           "COQ10A"            "CS"         "CNPY2"         "STAT2"      "TIMELESS" 
    ENSG00000076108 ENSG00000198056 ENSG00000166860 ENSG00000166886 ENSG00000166888 
            "BAZ2A"         "PRIM1"        "ZBTB39"          "NAB2"         "STAT6" 
    ENSG00000123384 ENSG00000182199 ENSG00000139269 ENSG00000111087 ENSG00000166986 
             "LRP1"         "SHMT2"         "INHBE"          "GLI1"         "MARS1" 
    ENSG00000166987 ENSG00000175203 ENSG00000155980 ENSG00000166908 ENSG00000178498 
             "MBD6"         "DCTN2"         "KIF5A"       "PIP4K2C"          "DTX3" 
    ENSG00000135506 ENSG00000135446 ENSG00000139266 ENSG00000037897 ENSG00000175215 
              "OS9"          "CDK4"       "MARCHF9"        "METTL1"        "CTDSP2" 
    ENSG00000135655 ENSG00000111110 ENSG00000177990 ENSG00000196935 ENSG00000184575 
            "USP15"         "PPM1H"       "DPY19L2"        "SRGAP1"          "XPOT" 
    ENSG00000183735 ENSG00000153179 ENSG00000135677 ENSG00000174106 ENSG00000174099 
             "TBK1"        "RASSF3"           "GNS"         "LEMD3"         "MSRB3" 
    ENSG00000149948 ENSG00000139233 ENSG00000155957 ENSG00000111530 ENSG00000111581 
            "HMGA2"          "LLPH"        "TMBIM4"         "CAND1"        "NUP107" 
    ENSG00000135679 ENSG00000135678 ENSG00000166225 ENSG00000111596 ENSG00000127329 
             "MDM2"           "CPM"          "FRS2"         "CNOT2"         "PTPRB" 
    ENSG00000133858 ENSG00000139291 ENSG00000080371 ENSG00000121749 ENSG00000072657 
           "ZFC3H1"        "TMEM19"         "RAB21"       "TBC1D15"         "TRHDE" 
    ENSG00000253719 ENSG00000111615 ENSG00000139289 ENSG00000179941 ENSG00000091039 
         "ATXN7L3B"          "KRR1"        "PHLDA1"         "BBS10"        "OSBPL8" 
    ENSG00000186908 ENSG00000165891 ENSG00000067798 ENSG00000067715 ENSG00000177425 
          "ZDHHC17"          "E2F7"          "NAV3"          "SYT1"          "PAWR" 
    ENSG00000058272 ENSG00000139304 ENSG00000111052 ENSG00000127720 ENSG00000179104 
         "PPP1R12A"         "PTPRQ"         "LIN7A"       "METTL25"         "TMTC2" 
    ENSG00000072041 ENSG00000231738 ENSG00000133640 ENSG00000139324 ENSG00000049130 
          "SLC6A15"       "TSPAN19"        "LRRIQ1"         "TMTC3"         "KITLG" 
    ENSG00000139318 ENSG00000139323 ENSG00000070961 ENSG00000139329 ENSG00000011465 
            "DUSP6"         "POC1B"        "ATP2B1"           "LUM"           "DCN" 
    ENSG00000133639 ENSG00000102189 ENSG00000173598 ENSG00000177889 ENSG00000198015 
             "BTG1"          "EEA1"         "NUDT4"         "UBE2N"        "MRPL42" 
    ENSG00000120833 ENSG00000169372 ENSG00000136040 ENSG00000173588 ENSG00000184752 
            "SOCS2"         "CRADD"        "PLXNC1"         "CEP83"       "NDUFA12" 
    ENSG00000180263 ENSG00000028203 ENSG00000136014 ENSG00000074527 ENSG00000139343 
             "FGD6"          "VEZT"         "USP44"          "NTN4"         "SNRPF" 
    ENSG00000139350 ENSG00000120802 ENSG00000075415 ENSG00000185046 ENSG00000111647 
            "NEDD1"          "TMPO"       "SLC25A3"        "ANKS1B"        "BLTP3B" 
    ENSG00000075089 ENSG00000136021 ENSG00000139354 ENSG00000151572 ENSG00000120800 
            "ACTR6"         "SCYL2"        "GAS2L3"          "ANO4"         "UTP20" 
    ENSG00000111666 ENSG00000136048 ENSG00000120860 ENSG00000075188 ENSG00000185480 
            "CHPT1"         "DRAM1"        "WASHC3"         "NUP37"        "PARPBP" 
    ENSG00000111696 ENSG00000166598 ENSG00000204954 ENSG00000120820 ENSG00000111727 
           "NT5DC3"       "HSP90B1"         "UQCC6"        "GLT8D2"         "HCFC2" 
    ENSG00000120837 ENSG00000198431 ENSG00000171310 ENSG00000136052 ENSG00000136010 
             "NFYB"        "TXNRD1"        "CHST11"       "SLC41A2"       "ALDH1L2" 
    ENSG00000136051 ENSG00000136044 ENSG00000235162 ENSG00000074590 ENSG00000136026 
           "WASHC4"         "APPL2"      "C12orf75"         "NUAK1"         "CKAP4" 
    ENSG00000166046 ENSG00000013503 ENSG00000151135 ENSG00000008405 ENSG00000136045 
          "TCP11L2"        "POLR3B"       "TMEM263"          "CRY1"          "PWP1" 
    ENSG00000198855 ENSG00000075856 ENSG00000136003 ENSG00000183160 ENSG00000110876 
             "FICD"         "SART3"          "ISCU"       "TMEM119"        "SELPLG" 
    ENSG00000110880 ENSG00000076248 ENSG00000076555 ENSG00000110906 ENSG00000139438 
           "CORO1C"           "UNG"         "ACACB"        "KCTD10"       "FAM222A" 
    ENSG00000111199 ENSG00000139433 ENSG00000139437 ENSG00000139436 ENSG00000076513 
            "TRPV4"          "GLTP"          "TCHP"          "GIT2"      "ANKRD13A" 
    ENSG00000174456 ENSG00000174437 ENSG00000196510 ENSG00000204856 ENSG00000111237 
         "C12orf76"        "ATP2A2"        "ANAPC7"       "FAM216A"         "VPS29" 
    ENSG00000196850 ENSG00000186298 ENSG00000198324 ENSG00000111252 ENSG00000204842 
            "PPTC7"        "PPP1CC"        "PHETA1"         "SH2B3"         "ATXN2" 
    ENSG00000111275 ENSG00000089022 ENSG00000089248 ENSG00000111300 ENSG00000135148 
            "ALDH2"      "MAPKAPK5"         "ERP29"         "NAA25"        "TRAFD1" 
    ENSG00000173064 ENSG00000089009 ENSG00000179295 ENSG00000111331 ENSG00000123064 
           "HECTD4"          "RPL6"        "PTPN11"          "OAS3"         "DDX54" 
    ENSG00000139405 ENSG00000186815 ENSG00000089060 ENSG00000151176 ENSG00000122965 
            "RITA1"         "TPCN1"        "SLC8B1"         "PLBD2"         "RBM19" 
    ENSG00000089225 ENSG00000135111 ENSG00000111412 ENSG00000174989 ENSG00000111445 
             "TBX5"          "TBX3"       "SPRING1"         "FBXW8"          "RFC5" 
    ENSG00000176871 ENSG00000176834 ENSG00000089220 ENSG00000135090 ENSG00000111707 
             "WSB2"        "VSIG10"         "PEBP1"         "TAOK3"         "SUDS3" 
    ENSG00000152137 ENSG00000111725 ENSG00000122966 ENSG00000111737 ENSG00000089154 
            "HSPB8"        "PRKAB1"           "CIT"         "RAB35"          "GCN1" 
    ENSG00000089157 ENSG00000089159 ENSG00000170855 ENSG00000257218 ENSG00000088986 
            "RPLP0"           "PXN"        "TRIAP1"          "GATC"        "DYNLL1" 
    ENSG00000110871 ENSG00000022840 ENSG00000110917 ENSG00000175970 ENSG00000122971 
             "COQ5"         "RNF10"          "MLEC"       "UNC119B"         "ACADS" 
    ENSG00000157837 ENSG00000157895 ENSG00000089041 ENSG00000110931 ENSG00000089094 
            "SPPL3"      "C12orf43"         "P2RX7"        "CAMKK2"         "KDM2B" 
    ENSG00000188735 ENSG00000139718 ENSG00000158023 ENSG00000110987 ENSG00000175727 
         "TMEM120B"        "SETD1B"       "CFAP251"         "BCL7A"         "MLXIP" 
    ENSG00000139719 ENSG00000130779 ENSG00000033030 ENSG00000111011 ENSG00000184445 
           "VPS33A"         "CLIP1"        "ZCCHC8"         "RSRC2"         "KNTC1" 
    ENSG00000130787 ENSG00000139722 ENSG00000090975 ENSG00000111328 ENSG00000139697 
            "HIP1R"        "VPS37B"       "PITPNM2"       "CDK2AP1"         "SBNO1" 
    ENSG00000183955 ENSG00000150977 ENSG00000086598 ENSG00000111361 ENSG00000111358 
            "KMT5A"        "RILPL2"         "TMED2"        "EIF2B1"        "GTF2H3" 
    ENSG00000168778 ENSG00000185344 ENSG00000119242 ENSG00000250091 ENSG00000179195 
            "TCTN2"      "ATP6V0A2"        "CCDC92"      "DNAH10OS"        "ZNF664" 
    ENSG00000196498 ENSG00000073060 ENSG00000184992 ENSG00000081760 ENSG00000139364 
            "NCOR2"        "SCARB1"        "BRI3BP"          "AACS"      "TMEM132B" 
    ENSG00000139370 ENSG00000111450 ENSG00000132341 ENSG00000061936 ENSG00000198598 
          "SLC15A4"          "STX2"           "RAN"        "SFSWAP"         "MMP17" 
    ENSG00000183495 ENSG00000185163 ENSG00000184967 ENSG00000112787 ENSG00000177084 
            "EP400"         "DDX51"         "NOC4L"        "FBRSL1"          "POLE" 
    ENSG00000176894 ENSG00000072609 ENSG00000196458 ENSG00000198393 ENSG00000198040 
            "PXMP2"          "CHFR"        "ZNF605"         "ZNF26"         "ZNF84" 
    ENSG00000196387 ENSG00000256223 ENSG00000196199 ENSG00000121741 ENSG00000121743 
           "ZNF140"         "ZNF10"      "MPHOSPH8"         "ZMYM2"          "GJA3" 
    ENSG00000132953 ENSG00000150457 ENSG00000165480 ENSG00000151835 ENSG00000127863 
             "XPO4"         "LATS2"          "SKA3"          "SACS"      "TNFRSF19" 
    ENSG00000027001 ENSG00000182957 ENSG00000151846 ENSG00000180730 ENSG00000127870 
            "MIPEP"       "SPATA13"        "PABPC3"        "SHISA2"          "RNF6" 
    ENSG00000132970 ENSG00000152484 ENSG00000122026 ENSG00000122034 ENSG00000139517 
            "WASF3"         "USP12"         "RPL21"         "GTF3A"          "LNX2" 
    ENSG00000186184 ENSG00000152520 ENSG00000102755 ENSG00000132963 ENSG00000139514 
           "POLR1D"          "PAN3"          "FLT1"          "POMP"        "SLC7A1" 
    ENSG00000122042 ENSG00000102781 ENSG00000189403 ENSG00000132952 ENSG00000102802 
             "UBL3"       "KATNAL1"         "HMGB1"         "USPL1"         "MEDAG" 
    ENSG00000120694 ENSG00000187676 ENSG00000073910 ENSG00000139618 ENSG00000139597 
            "HSPH1"        "B3GLCT"           "FRY"         "BRCA2"       "N4BP2L1" 
    ENSG00000083642 ENSG00000133121 ENSG00000133119 ENSG00000172915 ENSG00000133104 
            "PDS5B"       "STARD13"          "RFC3"          "NBEA"         "SPART" 
    ENSG00000133111 ENSG00000120697 ENSG00000120699 ENSG00000133110 ENSG00000133107 
            "RFXAP"          "ALG5"        "EXOSC8"         "POSTN"         "TRPC4" 
    ENSG00000120686 ENSG00000120685 ENSG00000183722 ENSG00000150907 ENSG00000102743 
             "UFM1"       "PROSER1"        "LHFPL6"         "FOXO1"      "SLC25A15" 
    ENSG00000120690 ENSG00000165572 ENSG00000120696 ENSG00000102780 ENSG00000023516 
             "ELF1"        "KBTBD6"        "KBTBD7"          "DGKH"        "AKAP11" 
    ENSG00000133106 ENSG00000120675 ENSG00000120658 ENSG00000179630 ENSG00000102804 
           "EPSTI1"       "DNAJC15"         "ENOX1"         "LACC1"       "TSC22D1" 
    ENSG00000180332 ENSG00000174032 ENSG00000123200 ENSG00000136141 ENSG00000139684 
            "KCTD4"      "SLC25A30"        "ZC3H13"         "LRCH1"           "ESD" 
    ENSG00000136143 ENSG00000136159 ENSG00000136156 ENSG00000139687 ENSG00000139679 
           "SUCLA2"        "NUDT15"         "ITM2B"           "RB1"         "LPAR6" 
    ENSG00000136161 ENSG00000102531 ENSG00000136147 ENSG00000102753 ENSG00000204977 
           "RCBTB2"        "FNDC3A"         "PHF11"         "KPNA3"        "TRIM13" 
    ENSG00000176124 ENSG00000136104 ENSG00000139668 ENSG00000253710 ENSG00000253797 
            "DLEU1"      "RNASEH2B"         "WDFY2"         "ALG11"        "UTP14C" 
    ENSG00000136114 ENSG00000136108 ENSG00000139675 ENSG00000165416 ENSG00000139734 
            "THSD1"         "CKAP2"     "HNRNPA1L2"         "SUGT1"        "DIAPH3" 
    ENSG00000184226 ENSG00000276644 ENSG00000204899 ENSG00000136122 ENSG00000083520 
            "PCDH9"         "DACH1"          "MZT1"          "BORA"          "DIS3" 
    ENSG00000083535 ENSG00000102554 ENSG00000118922 ENSG00000136111 ENSG00000118939 
            "PIBF1"          "KLF5"         "KLF12"        "TBC1D4"         "UCHL3" 
    ENSG00000136153 ENSG00000178695 ENSG00000005812 ENSG00000136160 ENSG00000136158 
             "LMO7"        "KCTD12"         "FBXL3"         "EDNRB"         "SPRY2" 
    ENSG00000184564 ENSG00000183098 ENSG00000152749 ENSG00000125257 ENSG00000134874 
          "SLITRK6"          "GPC6"        "GPR180"         "ABCC4"         "DZIP1" 
    ENSG00000102580 ENSG00000102595 ENSG00000139793 ENSG00000125249 ENSG00000065150 
           "DNAJC3"         "UGGT2"         "MBNL2"         "RAP2A"          "IPO5" 
    ENSG00000102572 ENSG00000134882 ENSG00000125304 ENSG00000175198 ENSG00000134864 
            "STK24"         "UBAC2"        "TM9SF2"          "PCCA"         "GGACT" 
    ENSG00000102452 ENSG00000198542 ENSG00000134900 ENSG00000134901 ENSG00000204442 
            "NALCN"        "ITGBL1"          "TPP2"       "POGLUT2"         "NALF1" 
    ENSG00000174405 ENSG00000139826 ENSG00000102524 ENSG00000185950 ENSG00000187498 
             "LIG4"        "ABHD13"      "TNFSF13B"          "IRS2"        "COL4A1" 
    ENSG00000134871 ENSG00000134905 ENSG00000153487 ENSG00000088448 ENSG00000102606 
           "COL4A2"         "CARS2"          "ING1"       "ANKRD10"       "ARHGEF7" 
    ENSG00000126216 ENSG00000126218 ENSG00000126226 ENSG00000185896 ENSG00000150403 
          "TUBGCP3"           "F10"         "PCID2"         "LAMP1"         "TMCO3" 
    ENSG00000198176 ENSG00000184497 ENSG00000183087 ENSG00000130177 ENSG00000169062 
            "TFDP1"      "TMEM255B"          "GAS6"         "CDC16"         "UPF3A" 
    ENSG00000198824 ENSG00000100814 ENSG00000129484 ENSG00000129566 ENSG00000100823 
           "CHAMP1"      "CCNB1IP1"         "PARP2"          "TEP1"         "APEX1" 
    ENSG00000165782 ENSG00000165801 ENSG00000165804 ENSG00000092199 ENSG00000129472 
           "PIP4P1"      "ARHGEF40"        "ZNF219"        "HNRNPC"         "RAB2B" 
    ENSG00000092203 ENSG00000165821 ENSG00000129562 ENSG00000100439 ENSG00000155463 
             "TOX4"         "SALL2"          "DAD1"         "ABHD4"         "OXA1L" 
    ENSG00000172590 ENSG00000157227 ENSG00000197324 ENSG00000100461 ENSG00000100462 
           "MRPL52"         "MMP14"         "LRP10"         "RBM23"         "PRMT5" 
    ENSG00000129474 ENSG00000139880 ENSG00000179933 ENSG00000092068 ENSG00000215271 
            "AJUBA"         "CDH24"     "C14orf119"        "SLC7A8"              NA 
    ENSG00000129473 ENSG00000100836 ENSG00000092096 ENSG00000136367 ENSG00000100889 
           "BCL2L2"        "PABPN1"      "SLC22A17"         "ZFHX2"          "PCK2" 
    ENSG00000100897 ENSG00000092010 ENSG00000100908 ENSG00000092098 ENSG00000100926 
           "DCAF11"         "PSME1"          "EMC9"         "RNF31"        "TM9SF1" 
    ENSG00000092330 ENSG00000100949 ENSG00000157379 ENSG00000196943 ENSG00000213903 
            "TINF2"       "RABGGTA"         "DHRS1"          "NOP9"         "LTB4R" 
    ENSG00000100968 ENSG00000100441 ENSG00000168952 ENSG00000129480 ENSG00000151413 
           "NFATC4"         "KHNYN"        "STXBP6"          "DTD2"         "NUBPL" 
    ENSG00000151320 ENSG00000151322 ENSG00000165389 ENSG00000129518 ENSG00000165410 
            "AKAP6"         "NPAS3"        "SPTSSA"          "EAPP"          "CFL2" 
    ENSG00000198604 ENSG00000151327 ENSG00000092020 ENSG00000174373 ENSG00000151332 
            "BAZ1A"      "FAM177A1"       "PPP2R3C"      "RALGAPA1"          "MBIP" 
    ENSG00000183032 ENSG00000139874 ENSG00000176435 ENSG00000182400 ENSG00000100941 
         "SLC25A21"         "SSTR1"       "CLEC14A"      "TRAPPC6B"           "PNN" 
    ENSG00000150527 ENSG00000165355 ENSG00000165379 ENSG00000179454 ENSG00000198718 
             "MIA2"        "FBXO33"         "LRFN5"        "KLHL28"      "TOGARAM1" 
    ENSG00000185246 ENSG00000100442 ENSG00000187790 ENSG00000129534 ENSG00000213741 
           "PRPF39"         "FKBP3"         "FANCM"      "MIS18BP1"         "RPS29" 
    ENSG00000165501 ENSG00000165502 ENSG00000165506 ENSG00000100479 ENSG00000165525 
             "LRR1"       "RPL36AL"        "DNAAF2"         "POLE2"          "NEMF" 
    ENSG00000012983 ENSG00000198513 ENSG00000151748 ENSG00000100503 ENSG00000100504 
           "MAP4K5"          "ATL1"          "SAV1"           "NIN"          "PYGL" 
    ENSG00000139921 ENSG00000139926 ENSG00000186469 ENSG00000087303 ENSG00000125384 
             "TMX1"         "FRMD6"          "GNG2"          "NID2"        "PTGER2" 
    ENSG00000087301 ENSG00000197930 ENSG00000198252 ENSG00000100522 ENSG00000073712 
          "TXNDC16"         "ERO1A"          "STYX"       "GNPNAT1"        "FERMT2" 
    ENSG00000100523 ENSG00000125378 ENSG00000100526 ENSG00000100528 ENSG00000197045 
            "DDHD1"          "BMP4"         "CDKN3"         "CNIH1"          "GMFB" 
    ENSG00000198554 ENSG00000168175 ENSG00000126787 ENSG00000178974 ENSG00000126775 
            "WDHD1"     "MAPK1IP1L"        "DLGAP5"        "FBXO34"         "ATG14" 
    ENSG00000126777 ENSG00000139946 ENSG00000070269 ENSG00000070367 ENSG00000139977 
             "KTN1"         "PELI2"       "TMEM260"         "EXOC5"         "NAA30" 
    ENSG00000165617 ENSG00000100592 ENSG00000126790 ENSG00000050130 ENSG00000126773 
            "DACT1"         "DAAM1"       "L3HYPDH"         "JKAMP"         "PCNX4" 
    ENSG00000100612 ENSG00000100614 ENSG00000126778 ENSG00000139974 ENSG00000100644 
            "DHRS7"         "PPM1A"          "SIX1"       "SLC38A6"         "HIF1A" 
    ENSG00000023608 ENSG00000140015 ENSG00000126785 ENSG00000126821 ENSG00000054654 
           "SNAPC1"         "KCNH5"          "RHOJ"         "SGPP1"         "SYNE2" 
    ENSG00000100714 ENSG00000089775 ENSG00000179841 ENSG00000126804 ENSG00000126803 
           "MTHFD1"        "ZBTB25"         "AKAP5"         "ZBTB1"         "HSPA2" 
    ENSG00000258289 ENSG00000125952 ENSG00000171723 ENSG00000072415 ENSG00000134001 
           "CHURC1"           "MAX"          "GPHN"         "PALS1"        "EIF2S1" 
    ENSG00000100568 ENSG00000072042 ENSG00000072121 ENSG00000185650 ENSG00000072110 
            "VTI1B"         "RDH11"       "ZFYVE26"       "ZFP36L1"         "ACTN1" 
    ENSG00000139990 ENSG00000081177 ENSG00000100626 ENSG00000100632 ENSG00000029364 
            "DCAF5"          "EXD2"       "GALNT16"           "ERH"       "SLC39A9" 
    ENSG00000100647 ENSG00000100650 ENSG00000198732 ENSG00000100731 ENSG00000197555 
            "SUSD6"         "SRSF5"         "SMOC1"         "PCNX1"       "SIPA1L1" 
    ENSG00000205683 ENSG00000119599 ENSG00000165861 ENSG00000080815 ENSG00000170468 
             "DPF3"         "DCAF4"        "ZFYVE1"         "PSEN1"         "RIOX1" 
    ENSG00000177465 ENSG00000176903 ENSG00000156030 ENSG00000119636 ENSG00000119711 
            "ACOT4"         "PNMA1"        "MIDEAS"         "BBOF1"       "ALDH6A1" 
    ENSG00000205659 ENSG00000119688 ENSG00000119655 ENSG00000119681 ENSG00000119682 
            "LIN52"         "ABCD4"          "NPC2"         "LTBP2"         "AREL1" 
    ENSG00000119689 ENSG00000119630 ENSG00000170348 ENSG00000170345 ENSG00000140044 
             "DLST"           "PGF"        "TMED10"           "FOS"          "JDP2" 
    ENSG00000119686 ENSG00000119685 ENSG00000133935 ENSG00000119699 ENSG00000089916 
           "FLVCR2"         "TTLL5"         "ERG28"         "TGFB3"      "GPATCH2L" 
    ENSG00000071246 ENSG00000119669 ENSG00000198894 ENSG00000009830 ENSG00000100580 
            "VASH1"       "IRF2BPL"          "CIPC"         "POMT2"         "TMED8" 
    ENSG00000100596 ENSG00000119705 ENSG00000100603 ENSG00000021645 ENSG00000211448 
           "SPTLC2"         "SLIRP"          "SNW1"         "NRXN3"          "DIO2" 
    ENSG00000100629 ENSG00000165417 ENSG00000071537 ENSG00000185070 ENSG00000054983 
           "CEP128"        "GTF2A1"         "SEL1L"         "FLRT2"          "GALC" 
    ENSG00000042317 ENSG00000070778 ENSG00000100722 ENSG00000165521 ENSG00000053254 
           "SPATA7"        "PTPN21"        "ZC3H14"          "EML5"         "FOXN3" 
    ENSG00000042088 ENSG00000100764 ENSG00000198668 ENSG00000165914 ENSG00000133943 
             "TDP1"         "PSMC1"         "CALM1"         "TTC7B"        "DGLUCY" 
    ENSG00000119714 ENSG00000100796 ENSG00000140092 ENSG00000066427 ENSG00000165934 
            "GPR68"       "PPP4R3A"         "FBLN5"         "ATXN3"         "CPSF2" 
    ENSG00000100600 ENSG00000066455 ENSG00000100605 ENSG00000165943 ENSG00000153485 
             "LGMN"        "GOLGA5"         "ITPK1"         "MOAP1"         "LYSET" 
    ENSG00000170270 ENSG00000012963 ENSG00000011114 ENSG00000100628 ENSG00000089737 
             "GON7"          "UBR7"         "BTBD7"          "ASB2"         "DDX24" 
    ENSG00000119632 ENSG00000119698 ENSG00000100697 ENSG00000176438 ENSG00000182512 
          "IFI27L2"        "PPP4R4"        "DICER1"         "SYNE3"         "GLRX5" 
    ENSG00000168398 ENSG00000066739 ENSG00000100744 ENSG00000090060 ENSG00000100749 
           "BDKRB2"         "ATG2B"         "GSKIP"        "PAPOLA"          "VRK1" 
    ENSG00000183576 ENSG00000090061 ENSG00000205476 ENSG00000066629 ENSG00000100811 
            "SETD3"          "CCNK"       "CCDC85C"          "EML1"           "YY1" 
    ENSG00000197119 ENSG00000140105 ENSG00000185559 ENSG00000197102 ENSG00000080824 
         "SLC25A29"         "WARS1"          "DLK1"       "DYNC1H1"      "HSP90AA1" 
    ENSG00000080823 ENSG00000100865 ENSG00000196663 ENSG00000198752 ENSG00000100664 
              "MOK"          "CINP"        "TECPR2"      "CDC42BPB"          "EIF5" 
    ENSG00000166170 ENSG00000126215 ENSG00000100711 ENSG00000088808 ENSG00000203485 
             "BAG5"         "XRCC3"       "ZFYVE21"      "PPP1R13B"          "INF2" 
    ENSG00000179627 ENSG00000099814 ENSG00000185567 ENSG00000170779 ENSG00000183828 
           "ZBTB42"       "CEP170B"        "AHNAK2"         "CDCA4"        "NUDT14" 
    ENSG00000184887 ENSG00000179364 ENSG00000182809 ENSG00000185347 ENSG00000170113 
            "BTBD6"         "PACS2"         "CRIP2"         "TEDC1"         "NIPA1" 
    ENSG00000273749 ENSG00000275835 ENSG00000254585 ENSG00000182636 ENSG00000128739 
           "CYFIP1"       "TUBGCP5"        "MAGEL2"           "NDN"         "SNRPN" 
    ENSG00000206190 ENSG00000186297 ENSG00000128731 ENSG00000034053 ENSG00000185115 
           "ATP10A"        "GABRA5"         "HERC2"         "APBA2"        "NSMCE3" 
    ENSG00000104067 ENSG00000187951 ENSG00000166912 ENSG00000169926 ENSG00000198826 
             "TJP1"  "LOC100288637"        "MTMR10"         "KLF13"     "ARHGAP11A" 
    ENSG00000166922 ENSG00000166923 ENSG00000248905 ENSG00000169857 ENSG00000134153 
             "SCG5"         "GREM1"          "FMN1"          "AVEN"          "EMC7" 
    ENSG00000134152 ENSG00000128463 ENSG00000140199 ENSG00000182117 ENSG00000176454 
          "KATNBL1"          "EMC4"       "SLC12A6"         "NOP10"        "LPCAT4" 
    ENSG00000175265 ENSG00000215252 ENSG00000198146 ENSG00000186073 ENSG00000134138 
          "GOLGA8A"       "GOLGA8B"        "ZNF770"         "CDIN1"         "MEIS2" 
    ENSG00000166068 ENSG00000171262 ENSG00000172575 ENSG00000137801 ENSG00000166073 
           "SPRED1"        "FAM98B"       "RASGRP1"         "THBS1"        "GPR176" 
    ENSG00000128829 ENSG00000140319 ENSG00000104081 ENSG00000156970 ENSG00000188549 
          "EIF2AK4"         "SRP14"           "BMF"         "BUB1B"        "CCDC9B" 
    ENSG00000128944 ENSG00000128928 ENSG00000169105 ENSG00000137812 ENSG00000051180 
           "KNSTRN"           "IVD"        "CHST14"          "KNL1"         "RAD51" 
    ENSG00000104129 ENSG00000166140 ENSG00000166145 ENSG00000104142 ENSG00000128917 
          "DNAJC17"       "ZFYVE19"        "SPINT1"         "VPS18"          "DLL4" 
    ENSG00000128965 ENSG00000128908 ENSG00000187446 ENSG00000104147 ENSG00000137804 
            "CHAC1"         "INO80"          "CHP1"          "OIP5"        "NUSAP1" 
    ENSG00000137806 ENSG00000137815 ENSG00000092445 ENSG00000174197 ENSG00000103966 
          "NDUFAF1"          "RTF1"         "TYRO3"           "MGA"          "EHD4" 
    ENSG00000166887 ENSG00000103978 ENSG00000103994 ENSG00000180979 ENSG00000137814 
            "VPS39"       "TMEM87A"        "ZNF106"        "LRRC57"         "HAUS2" 
    ENSG00000140326 ENSG00000159459 ENSG00000140265 ENSG00000137822 ENSG00000067369 
            "CDAN1"          "UBR1"       "ZSCAN29"       "TUBGCP4"       "TP53BP1" 
    ENSG00000166963 ENSG00000168781 ENSG00000167004 ENSG00000140259 ENSG00000092470 
            "MAP1A"       "PPIP5K1"         "PDIA3"         "MFAP1"         "WDR76" 
    ENSG00000171877 ENSG00000137770 ENSG00000104131 ENSG00000166710 ENSG00000185880 
            "FRMD5"       "CTDSPL2"         "EIF3J"           "B2M"        "TRIM69" 
    ENSG00000140263 ENSG00000138606 ENSG00000171763 ENSG00000104154 ENSG00000104164 
             "SORD"           "SHF"         "AFG2B"       "SLC30A4"       "BLOC1S6" 
    ENSG00000137767 ENSG00000137872 ENSG00000104177 ENSG00000128951 ENSG00000166147 
             "SQOR"        "SEMA6D"         "MYEF2"           "DUT"          "FBN1" 
    ENSG00000103995 ENSG00000255302 ENSG00000138593 ENSG00000166200 ENSG00000140285 
           "CEP152"          "EID1"     "SECISBP2L"         "COPS2"          "FGF7" 
    ENSG00000104047 ENSG00000104043 ENSG00000092439 ENSG00000081014 ENSG00000104093 
            "DTWD1"        "ATP8B4"         "TRPM7"         "AP4E1"         "DMXL2" 
    ENSG00000140280 ENSG00000128872 ENSG00000138594 ENSG00000069956 ENSG00000197535 
           "LYSMD2"         "TMOD2"         "TMOD3"         "MAPK6"         "MYO5A" 
    ENSG00000128989 ENSG00000047346 ENSG00000137766 ENSG00000171016 ENSG00000138587 
           "ARPP19"         "ATOSA"        "UNC13C"         "PYGO1"          "MNS1" 
    ENSG00000140262 ENSG00000255529 ENSG00000137845 ENSG00000128923 ENSG00000157450 
            "TCF12"        "POLR2M"        "ADAM10"        "MINDY2"        "RNF111" 
    ENSG00000157456 ENSG00000157483 ENSG00000140299 ENSG00000182718 ENSG00000069667 
            "CCNB2"         "MYO1E"         "BNIP2"         "ANXA2"          "RORA" 
    ENSG00000129003 ENSG00000140416 ENSG00000103642 ENSG00000185088 ENSG00000166128 
           "VPS13C"          "TPM1"         "LACTB"        "RPS27L"         "RAB8B" 
    ENSG00000138613 ENSG00000074410 ENSG00000166797 ENSG00000166794 ENSG00000169118 
            "APH1B"          "CA12"        "CIAO2A"          "PPIB"       "CSNK1G1" 
    ENSG00000166803 ENSG00000103671 ENSG00000180357 ENSG00000180304 ENSG00000166831 
            "PCLAF"         "TRIP4"        "ZNF609"          "OAZ2"        "RBPMS2" 
    ENSG00000140451 ENSG00000241839 ENSG00000090487 ENSG00000103707 ENSG00000138617 
             "PIF1"       "PLEKHO2"         "SPG21"         "MTFMT"        "PARP16" 
    ENSG00000074603 ENSG00000074696 ENSG00000138614 ENSG00000074621 ENSG00000166938 
             "DPP8"         "HACD3"        "INTS14"       "SLC24A1"         "DIS3L" 
    ENSG00000075131 ENSG00000169032 ENSG00000174446 ENSG00000174444 ENSG00000174442 
            "TIPIN"        "MAP2K1"        "SNAPC5"          "RPL4"        "ZWILCH" 
    ENSG00000188501 ENSG00000166949 ENSG00000103591 ENSG00000137764 ENSG00000033800 
             "LCTL"         "SMAD3"         "AAGAB"        "MAP2K5"         "PIAS1" 
    ENSG00000128973 ENSG00000169018 ENSG00000137809 ENSG00000103647 ENSG00000140350 
             "CLN6"         "FEM1B"        "ITGA11"        "CORO2B"        "ANP32A" 
    ENSG00000138604 ENSG00000137819 ENSG00000137807 ENSG00000137818 ENSG00000140332 
             "GLCE"         "PAQR5"         "KIF23"         "RPLP1"          "TLE3" 
    ENSG00000137831 ENSG00000166173 ENSG00000137821 ENSG00000187720 ENSG00000067225 
             "UACA"         "LARP6"        "LRRC49"         "THSD4"           "PKM" 
    ENSG00000137817 ENSG00000213614 ENSG00000159322 ENSG00000067141 ENSG00000103855 
            "PARP6"          "HEXA"         "ADPGK"          "NEO1"         "CD276" 
    ENSG00000129038 ENSG00000067221 ENSG00000140464 ENSG00000138623 ENSG00000179335 
            "LOXL1"        "STOML1"           "PML"        "SEMA7A"          "CLK3" 
    ENSG00000179151 ENSG00000140465 ENSG00000103653 ENSG00000140474 ENSG00000140497 
             "EDC3"        "CYP1A1"           "CSK"          "ULK3"        "SCAMP2" 
    ENSG00000178802 ENSG00000178761 ENSG00000198794 ENSG00000169410 ENSG00000169371 
              "MPI"       "FAM219B"        "SCAMP5"         "PTPN9"         "SNUPN" 
    ENSG00000177971 ENSG00000173548 ENSG00000173546 ENSG00000140367 ENSG00000167196 
             "IMP3"         "SNX33"         "CSPG4"        "UBE2Q2"        "FBXO22" 
    ENSG00000140374 ENSG00000117906 ENSG00000140391 ENSG00000173517 ENSG00000167202 
             "ETFA"          "RCN2"        "TSPAN3"         "PEAK1"       "TBC1D2B" 
    ENSG00000166411 ENSG00000136381 ENSG00000188266 ENSG00000041357 ENSG00000136378 
            "IDH3A"         "IREB2"          "HYKK"         "PSMA4"       "ADAMTS7" 
    ENSG00000185787 ENSG00000086666 ENSG00000103876 ENSG00000172379 ENSG00000136379 
          "MORF4L1"        "ZFAND6"           "FAH"         "ARNT2"       "ABHD17C" 
    ENSG00000140406 ENSG00000172345 ENSG00000182774 ENSG00000103942 ENSG00000169612 
           "TLNRD1"        "STARD5"         "RPS17"        "HOMER2"         "RAMAC" 
    ENSG00000064726 ENSG00000136404 ENSG00000176371 ENSG00000140612 ENSG00000073417 
            "BTBD1"        "TM6SF1"        "ZSCAN2"        "SEC11A"         "PDE8A" 
    ENSG00000259494 ENSG00000181026 ENSG00000157766 ENSG00000140511 ENSG00000140545 
           "MRPL46"           "AEN"          "ACAN"        "HAPLN3"         "MFGE8" 
    ENSG00000140525 ENSG00000140521 ENSG00000140534 ENSG00000166813 ENSG00000166821 
            "FANCI"          "POLG"         "TICRR"          "KIF7"        "PEX11A" 
    ENSG00000157823 ENSG00000242498 ENSG00000140548 ENSG00000182054 ENSG00000185033 
            "AP3S2"         "ARPIN"        "ZNF710"          "IDH2"        "SEMA4B" 
    ENSG00000185043 ENSG00000182768 ENSG00000140575 ENSG00000140577 ENSG00000197299 
             "CIB1"          "NGRN"        "IQGAP1"         "CRTC3"           "BLM" 
    ENSG00000140564 ENSG00000182511 ENSG00000196547 ENSG00000184508 ENSG00000140553 
            "FURIN"           "FES"        "MAN2A2"         "HDDC3"        "UNC45A" 
    ENSG00000166965 ENSG00000198901 ENSG00000176463 ENSG00000140563 ENSG00000185551 
            "RCCD1"          "PRC1"       "SLCO3A1"         "MCTP2"         "NR2F2" 
    ENSG00000140450 ENSG00000140443 ENSG00000182253 ENSG00000103852 ENSG00000068305 
           "ARRDC4"         "IGF1R"          "SYNM"         "TTC23"         "MEF2A" 
    ENSG00000183060 ENSG00000184254 ENSG00000131871 ENSG00000131876 ENSG00000184277 
           "LYSMD4"       "ALDH1A3"       "SELENOS"        "SNRPA1"         "TM2D3" 
    ENSG00000185418 ENSG00000161981 ENSG00000007384 ENSG00000103152 ENSG00000103148 
            "TARS3"       "SNRNP25"        "RHBDF1"           "MPG"         "NPRL3" 
    ENSG00000007392 ENSG00000167930 ENSG00000076344 ENSG00000103126 ENSG00000086504 
            "LUC7L"       "FAM234A"         "RGS11"         "AXIN1"        "MRPL28" 
    ENSG00000129925 ENSG00000103202 ENSG00000242612 ENSG00000090565 ENSG00000103326 
            "PGAP6"          "NME4"         "DECR2"     "RAB11FIP3"        "CAPN15" 
    ENSG00000007541 ENSG00000257108 ENSG00000197562 ENSG00000103266 ENSG00000127580 
             "PIGQ"        "NHLRC4"        "RAB40C"         "STUB1"         "WDR24" 
    ENSG00000127585 ENSG00000103260 ENSG00000103253 ENSG00000007376 ENSG00000127586 
           "FBXL16"         "METRN"         "HAGHL"        "RPUSD1"        "CHTF18" 
    ENSG00000103227 ENSG00000196557 ENSG00000103275 ENSG00000007516 ENSG00000090581 
             "LMF1"       "CACNA1H"         "UBE2I"        "BAIAP3"         "GNPTG" 
    ENSG00000059145 ENSG00000103249 ENSG00000100726 ENSG00000187535 ENSG00000131634 
             "UNKL"         "CLCN7"         "TELO2"        "IFT140"       "TMEM204" 
    ENSG00000206053 ENSG00000138834 ENSG00000074071 ENSG00000197774 ENSG00000095906 
             "JPT2"      "MAPK8IP3"        "MRPS34"          "EME2"         "NUBP2" 
    ENSG00000140990 ENSG00000140988 ENSG00000183751 ENSG00000167962 ENSG00000065054 
          "NDUFB10"          "RPS2"          "TBL3"        "ZNF598"        "NHERF2" 
    ENSG00000103197 ENSG00000008710 ENSG00000167965 ENSG00000184207 ENSG00000205937 
             "TSC2"          "PKD1"         "MLST8"           "PGP"         "RNPS1" 
    ENSG00000167972 ENSG00000162063 ENSG00000162062 ENSG00000162068 ENSG00000185883 
            "ABCA3"          "CCNF"         "TEDC2"          "NTN3"       "ATP6V0C" 
    ENSG00000162066 ENSG00000140992 ENSG00000167977 ENSG00000005001 ENSG00000162076 
           "AMDHD2"         "PDPK1"         "KCTD5"        "PRSS22"       "FLYWCH2" 
    ENSG00000059122 ENSG00000127564 ENSG00000006327 ENSG00000008517 ENSG00000085644 
          "FLYWCH1"        "PKMYT1"     "TNFRSF12A"          "IL32"        "ZNF213" 
    ENSG00000006194 ENSG00000140993 ENSG00000167981 ENSG00000122390 ENSG00000167984 
           "ZNF263"         "TIGD7"        "ZNF597"         "NAA60"         "NLRC3" 
    ENSG00000126602 ENSG00000005339 ENSG00000162104 ENSG00000126603 ENSG00000168140 
            "TRAP1"        "CREBBP"         "ADCY9"         "GLIS2"          "VASN" 
    ENSG00000103423 ENSG00000103415 ENSG00000153443 ENSG00000102858 ENSG00000168096 
           "DNAJA3"         "HMOX2"        "UBALD1"         "MGRN1"         "ANKS3" 
    ENSG00000103199 ENSG00000067836 ENSG00000140632 ENSG00000118900 ENSG00000118894 
           "ZNF500"         "ROGDI"         "GLYR1"          "UBN1"       "EEF2KMT" 
    ENSG00000183044 ENSG00000187555 ENSG00000182831 ENSG00000183454 ENSG00000213853 
             "ABAT"          "USP7"       "HAPSTR1"        "GRIN2A"          "EMP2" 
    ENSG00000175643 ENSG00000185338 ENSG00000189067 ENSG00000184602 ENSG00000153066 
             "RMI2"         "SOCS1"         "LITAF"           "SNN"       "TXNDC11" 
    ENSG00000171490 ENSG00000103342 ENSG00000048471 ENSG00000103381 ENSG00000186260 
           "RSL1D1"         "GSPT1"         "SNX29"        "CPPED1"         "MRTFB" 
    ENSG00000140694 ENSG00000103512 ENSG00000179889 ENSG00000166780 ENSG00000072864 
             "PARN"         "NOMO1"        "PDXDC1"        "BMERB1"          "NDE1" 
    ENSG00000133392 ENSG00000103222 ENSG00000103489 ENSG00000170540 ENSG00000157106 
            "MYH11"         "ABCC1"         "XYLT1"       "ARL6IP1"          "SMG1" 
    ENSG00000170537 ENSG00000205730 ENSG00000103534 ENSG00000006007 ENSG00000103544 
             "TMC7"      "ITPRIPL2"          "TMC5"          "GDE1"        "VPS35L" 
    ENSG00000174628 ENSG00000167191 ENSG00000066654 ENSG00000196678 ENSG00000005189 
             "IQCK"        "GPRC5B"       "THUMPD1"          "ERI2"         "REXO5" 
    ENSG00000188215 ENSG00000102897 ENSG00000011638 ENSG00000197006 ENSG00000140740 
          "DCUN1D3"         "LYRM1"         "LDAF1"        "METTL9"        "UQCRC2" 
    ENSG00000185716 ENSG00000103319 ENSG00000140743 ENSG00000103404 ENSG00000103365 
            "MOSMO"         "EEF2K"          "CDR2"         "USP31"          "GGA2" 
    ENSG00000103356 ENSG00000103353 ENSG00000083093 ENSG00000166847 ENSG00000166851 
            "EARS2"         "UBFD1"         "PALB2"         "DCTN5"          "PLK1" 
    ENSG00000122257 ENSG00000090905 ENSG00000140750 ENSG00000155592 ENSG00000077238 
            "RBBP6"        "TNRC6A"      "ARHGAP17"       "ZKSCAN2"          "IL4R" 
    ENSG00000103522 ENSG00000169180 ENSG00000188603 ENSG00000176046 ENSG00000184110 
            "IL21R"          "XPO6"          "CLN3"         "NUPR1"         "EIF3C" 
    ENSG00000176953 ENSG00000254206 ENSG00000079616 ENSG00000167371 ENSG00000013364 
         "NFATC2IP"       "NPIPB11"         "KIF22"         "PRRT2"           "MVP" 
    ENSG00000103502 ENSG00000174938 ENSG00000174939 ENSG00000174943 ENSG00000149930 
            "CDIPT"        "SEZ6L2"        "ASPHD1"        "KCTD13"         "TAOK2" 
    ENSG00000149929 ENSG00000149925 ENSG00000149923 ENSG00000090238 ENSG00000102882 
           "HIRIP3"         "ALDOA"         "PPP4C"         "YPEL3"         "MAPK3" 
    ENSG00000102879 ENSG00000180035 ENSG00000179958 ENSG00000179918 ENSG00000169951 
           "CORO1A"         "ZNF48"        "DCTPP1"        "SEPHS2"        "ZNF764" 
    ENSG00000229809 ENSG00000197162 ENSG00000156853 ENSG00000080603 ENSG00000103549 
           "ZNF688"        "ZNF785"        "ZNF689"         "SRCAP"         "RNF40" 
    ENSG00000102870 ENSG00000099385 ENSG00000150281 ENSG00000099364 ENSG00000175938 
           "ZNF629"         "BCL7C"          "CTF1"        "FBXL19"         "ORAI3" 
    ENSG00000099381 ENSG00000099365 ENSG00000103496 ENSG00000167395 ENSG00000103507 
           "SETD1A"         "STX1B"          "STX4"        "ZNF646"         "BCKDK" 
    ENSG00000089280 ENSG00000103490 ENSG00000140691 ENSG00000140682 ENSG00000197302 
              "FUS"        "PYCARD"         "ARMC5"       "TGFB1I1"        "KRBOX5" 
    ENSG00000185947 ENSG00000171241 ENSG00000069329 ENSG00000091651 ENSG00000155330 
           "ZNF267"        "SHCBP1"         "VPS35"          "ORC6"      "C16orf87" 
    ENSG00000166123 ENSG00000069345 ENSG00000129636 ENSG00000102893 ENSG00000102935 
             "GPT2"        "DNAJA2"         "ITFG1"          "PHKB"        "ZNF423" 
    ENSG00000155393 ENSG00000121281 ENSG00000140807 ENSG00000103449 ENSG00000177200 
           "HEATR3"         "ADCY7"          "NKD1"         "SALL1"          "CHD9" 
    ENSG00000103479 ENSG00000166971 ENSG00000103494 ENSG00000140718 ENSG00000087245 
             "RBL2"         "AKTIP"      "RPGRIP1L"           "FTO"          "MMP2" 
    ENSG00000087253 ENSG00000159461 ENSG00000167005 ENSG00000087263 ENSG00000125124 
           "LPCAT2"          "AMFR"        "NUDT21"        "OGFOD1"          "BBS2" 
    ENSG00000125148 ENSG00000169715 ENSG00000198417 ENSG00000187193 ENSG00000102900 
             "MT2A"          "MT1E"          "MT1F"          "MT1X"         "NUP93" 
    ENSG00000051108 ENSG00000140848 ENSG00000159579 ENSG00000102931 ENSG00000102934 
          "HERPUD1"         "CPNE2"        "RSPRY1"        "ARL2BP"          "PLLP" 
    ENSG00000006210 ENSG00000088682 ENSG00000135736 ENSG00000205336 ENSG00000140854 
           "CX3CL1"          "COQ9"      "CCDC102A"        "ADGRG1"        "KATNB1" 
    ENSG00000103005 ENSG00000102996 ENSG00000070761 ENSG00000070770 ENSG00000181938 
             "USB1"         "MMP15"        "CFAP20"       "CSNK2A2"         "GINS3" 
    ENSG00000103034 ENSG00000125107 ENSG00000103042 ENSG00000140937 ENSG00000166548 
            "NDRG4"         "CNOT1"       "SLC38A7"         "CDH11"           "TK2" 
    ENSG00000140931 ENSG00000172840 ENSG00000172831 ENSG00000125149 ENSG00000237172 
            "CMTM3"          "PDP2"          "CES2"         "PHAF1"        "B3GNT9" 
    ENSG00000102871 ENSG00000135722 ENSG00000196123 ENSG00000205250 ENSG00000168701 
            "TRADD"         "FBXL8"       "MATCAP1"          "E2F4"       "TMEM208" 
    ENSG00000135740 ENSG00000196155 ENSG00000159714 ENSG00000102974 ENSG00000124074 
           "SLC9A5"       "PLEKHG4"        "ZDHHC1"          "CTCF"         "ENKD1" 
    ENSG00000141098 ENSG00000141084 ENSG00000102901 ENSG00000038358 ENSG00000159792 
            "GFOD2"       "RANBP10"         "CENPT"          "EDC4"         "PSKH1" 
    ENSG00000205220 ENSG00000124067 ENSG00000072736 ENSG00000103066 ENSG00000103064 
           "PSMB10"       "SLC12A4"        "NFATC3"       "PLA2G15"        "SLC7A6" 
    ENSG00000103061 ENSG00000132600 ENSG00000103047 ENSG00000141076 ENSG00000168807 
         "SLC7A6OS"         "PRMT7"        "TANGO6"          "UTP4"         "SNTB2" 
    ENSG00000132612 ENSG00000132603 ENSG00000132604 ENSG00000103018 ENSG00000102908 
            "VPS4A"          "NIP7"         "TERF2"         "CYB5B"         "NFAT5" 
    ENSG00000181019 ENSG00000198373 ENSG00000090861 ENSG00000157350 ENSG00000103051 
             "NQO1"          "WWP2"         "AARS1"       "ST3GAL2"          "COG4" 
    ENSG00000189091 ENSG00000157368 ENSG00000132613 ENSG00000103043 ENSG00000180917 
            "SF3B3"          "IL34"         "MTSS2"         "VAC14"         "CMTR2" 
    ENSG00000172137 ENSG00000040199 ENSG00000166747 ENSG00000224470 ENSG00000102967 
            "CALB2"        "PHLPP2"         "AP1G1"        "ATXN1L"         "DHODH" 
    ENSG00000140836 ENSG00000168411 ENSG00000184517 ENSG00000050820 ENSG00000153774 
            "ZFHX3"         "RFWD3"          "ZFP1"         "BCAR1"         "CFDP1" 
    ENSG00000166822 ENSG00000205084 ENSG00000034713 ENSG00000065427 ENSG00000103111 
         "TMEM170A"       "TMEM231"     "GABARAPL2"         "KARS1"         "MON1B" 
    ENSG00000171724 ENSG00000166446 ENSG00000103121 ENSG00000166451 ENSG00000140905 
            "VAT1L"         "CDYL2"          "CMC2"         "CENPN"          "GCSH" 
    ENSG00000153815 ENSG00000086696 ENSG00000135698 ENSG00000230989 ENSG00000260300 
             "CMIP"       "HSD17B2"      "MPHOSPH6"         "HSBP1"              NA 
    ENSG00000140961 ENSG00000140943 ENSG00000103160 ENSG00000103168 ENSG00000103175 
           "OSGIN1"        "MBTPS1"         "HSDL1"         "TAF1C"         "WFDC1" 
    ENSG00000140950 ENSG00000103187 ENSG00000135686 ENSG00000103194 ENSG00000103196 
            "MEAK7"         "COTL1"        "KLHL36"         "USP10"      "CRISPLD2" 
    ENSG00000153786 ENSG00000135709 ENSG00000131149 ENSG00000131153 ENSG00000131148 
           "ZDHHC7"      "KIAA0513"          "GSE1"         "GINS2"          "EMC8" 
    ENSG00000131143 ENSG00000103241 ENSG00000176692 ENSG00000176678 ENSG00000103264 
           "COX4I1"         "FOXF1"         "FOXC2"         "FOXL1"        "FBXO31" 
    ENSG00000140941 ENSG00000140948 ENSG00000104731 ENSG00000103257 ENSG00000225614 
         "MAP1LC3B"       "ZCCHC14"        "KLHDC4"        "SLC7A5"        "ZNF469" 
    ENSG00000179588 ENSG00000051523 ENSG00000174177 ENSG00000103335 ENSG00000167513 
            "ZFPM1"          "CYBA"          "CTU2"        "PIEZO1"          "CDT1" 
    ENSG00000198931 ENSG00000167526 ENSG00000167523 ENSG00000185324 ENSG00000158792 
             "APRT"         "RPL13"       "SPATA33"         "CDK10"       "SPATA2L" 
    ENSG00000075399 ENSG00000158805 ENSG00000187741 ENSG00000204991 ENSG00000258947 
           "VPS9D1"        "ZNF276"         "FANCA"        "SPIRE2"         "TUBB3" 
    ENSG00000140995 ENSG00000177946 ENSG00000003249 ENSG00000141013 ENSG00000183688 
             "DEF8"              NA        "DBNDD1"          "GAS8"         "RFLNB" 
    ENSG00000141252 ENSG00000167695 ENSG00000179409 ENSG00000167693 ENSG00000159842 
            "VPS53"        "TLCD3A"        "GEMIN4"           "NXN"           "ABR" 
    ENSG00000108953 ENSG00000167193 ENSG00000197879 ENSG00000132376 ENSG00000167705 
            "YWHAE"           "CRK"         "MYO1C"        "INPP5K"          "RILP" 
    ENSG00000174231 ENSG00000185561 ENSG00000167716 ENSG00000132383 ENSG00000177374 
            "PRPF8"         "TLCD2"         "WDR81"          "RPA1"          "HIC1" 
    ENSG00000070366 ENSG00000167720 ENSG00000167721 ENSG00000070444 ENSG00000040531 
             "SMG6"           "SRR"          "TSR1"           "MNT"          "CTNS" 
    ENSG00000213977 ENSG00000177602 ENSG00000074356 ENSG00000004660 ENSG00000074370 
          "TAX1BP3"        "HASPIN"         "NCBP3"        "CAMKK1"        "ATP2A3" 
    ENSG00000167740 ENSG00000132382 ENSG00000141456 ENSG00000141480 ENSG00000142507 
           "CYB5D2"       "MYBBP1A"         "PELP1"         "ARRB2"         "PSMB6" 
    ENSG00000129219 ENSG00000141503 ENSG00000108523 ENSG00000108518 ENSG00000108509 
             "PLD2"         "MINK1"        "RNF167"          "PFN1"        "CAMTA2" 
    ENSG00000129250 ENSG00000180787 ENSG00000167840 ENSG00000029725 ENSG00000108559 
            "KIF1C"          "ZFP3"        "ZNF232"        "RABEP1"         "NUP88" 
    ENSG00000129197 ENSG00000108561 ENSG00000072849 ENSG00000167842 ENSG00000129195 
            "RPAIN"         "C1QBP"         "DERL2"         "MIS12"        "PIMREG" 
    ENSG00000198920 ENSG00000108590 ENSG00000108839 ENSG00000141505 ENSG00000132535 
         "KIAA0753"         "MED31"        "ALOX12"         "ASGR1"          "DLG4" 
    ENSG00000072778 ENSG00000004975 ENSG00000040633 ENSG00000132507 ENSG00000215041 
           "ACADVL"          "DVL2"         "PHF23"         "EIF5A"        "NEURL4" 
    ENSG00000174292 ENSG00000169992 ENSG00000174282 ENSG00000259224 ENSG00000181222 
             "TNK1"         "NLGN2"         "ZBTB4"       "SLC35G6"        "POLR2A" 
    ENSG00000239697 ENSG00000129255 ENSG00000129245 ENSG00000141504 ENSG00000141510 
          "TNFSF12"         "MPDU1"          "FXR2"          "SAT2"          "TP53" 
    ENSG00000108947 ENSG00000132510 ENSG00000167874 ENSG00000170004 ENSG00000170037 
            "EFNB3"         "KDM6B"        "TMEM88"          "CHD3"        "CNTROB" 
    ENSG00000179111 ENSG00000179094 ENSG00000220205 ENSG00000196544 ENSG00000178999 
             "HES7"          "PER1"         "VAMP2"        "BORCS6"         "AURKB" 
    ENSG00000178971 ENSG00000178921 ENSG00000161970 ENSG00000166579 ENSG00000133026 
             "CTC1"          "PFAS"         "RPL26"         "NDEL1"         "MYH10" 
    ENSG00000065320 ENSG00000133028 ENSG00000170222 ENSG00000065559 ENSG00000141052 
             "NTN1"          "SCO1"         "ADPRM"        "MAP2K4"         "MYOCD" 
    ENSG00000006744 ENSG00000153976 ENSG00000109099 ENSG00000125409 ENSG00000221926 
            "ELAC2"      "HS3ST3A1"         "PMP22"         "TEKT3"        "TRIM16" 
    ENSG00000187607 ENSG00000170425 ENSG00000011295 ENSG00000141027 ENSG00000166582 
          "ZNF286A"       "ADORA2B"         "TTC19"         "NCOR1"         "CENPV" 
    ENSG00000187688 ENSG00000154803 ENSG00000141030 ENSG00000072310 ENSG00000175662 
            "TRPV2"          "FLCN"         "COPS3"        "SREBF1"        "TOM1L2" 
    ENSG00000141034 ENSG00000177731 ENSG00000176994 ENSG00000176974 ENSG00000249459 
             "GID4"          "FLII"         "SMCR8"         "SHMT1"              NA 
    ENSG00000141127 ENSG00000188522 ENSG00000108641 ENSG00000166482 ENSG00000128482 
          "PRPSAP2"        "FAM83G"          "B9D1"         "MFAP4"        "RNF112" 
    ENSG00000083290 ENSG00000124422 ENSG00000109016 ENSG00000274180 ENSG00000034152 
             "ULK2"         "USP22"        "DHRS7B"         "NATD1"        "MAP2K3" 
    ENSG00000212719 ENSG00000168961 ENSG00000109083 ENSG00000109079 ENSG00000004142 
        "LINC02693"        "LGALS9"         "IFT20"       "TNFAIP1"       "POLDIP2" 
    ENSG00000244045 ENSG00000004139 ENSG00000076351 ENSG00000087111 ENSG00000076382 
          "TMEM199"         "SARM1"       "SLC46A1"          "PIGS"         "SPAG5" 
    ENSG00000132581 ENSG00000109113 ENSG00000198242 ENSG00000160602 ENSG00000076604 
             "SDF2"         "RAB34"        "RPL23A"          "NEK8"         "TRAF4" 
    ENSG00000173065 ENSG00000132591 ENSG00000132589 ENSG00000109118 ENSG00000196535 
          "FAM222B"         "ERAL1"         "FLOT2"         "PHF12"        "MYO18A" 
    ENSG00000108256 ENSG00000160551 ENSG00000167543 ENSG00000108262 ENSG00000198720 
           "NUFIP2"         "TAOK1"       "TP53I13"          "GIT1"      "ANKRD13B" 
    ENSG00000167549 ENSG00000141298 ENSG00000108578 ENSG00000108582 ENSG00000108587 
            "CORO6"          "SSH2"          "BLMH"           "CPD"         "GOSR1" 
    ENSG00000176390 ENSG00000176208 ENSG00000184060 ENSG00000181481 ENSG00000196712 
            "CRLF3"         "ATAD5"         "ADAP2"        "RNF135"           "NF1" 
    ENSG00000126860 ENSG00000131242 ENSG00000172301 ENSG00000178691 ENSG00000126858 
            "EVI2A"     "RAB11FIP4"         "COPRS"         "SUZ12"         "RHOT1" 
    ENSG00000010244 ENSG00000108671 ENSG00000176749 ENSG00000176658 ENSG00000006042 
           "ZNF207"        "PSMD11"        "CDK5R1"         "MYO1D"        "TMEM98" 
    ENSG00000198783 ENSG00000005156 ENSG00000092871 ENSG00000185379 ENSG00000073536 
           "ZNF830"          "LIG3"          "RFFL"        "RAD51D"          "NLE1" 
    ENSG00000166750 ENSG00000172716 ENSG00000154760 ENSG00000108733 ENSG00000006125 
            "SLFN5"        "SLFN11"        "SLFN13"         "PEX12"         "AP2B1" 
    ENSG00000270885 ENSG00000273611 ENSG00000278259 ENSG00000278540 ENSG00000276234 
          "RASL10B"        "ZNHIT3"         "MYO19"         "ACACA"        "TADA2A" 
    ENSG00000276023 ENSG00000275066 ENSG00000278053 ENSG00000278845 ENSG00000274211 
           "DUSP14"         "SYNRG"         "DDX52"        "MRPL45"         "SOCS7" 
    ENSG00000275832 ENSG00000275023 ENSG00000277258 ENSG00000002834 ENSG00000067191 
         "ARHGAP23"         "MLLT6"         "PCGF2"         "LASP1"        "CACNB1" 
    ENSG00000108298 ENSG00000108306 ENSG00000125686 ENSG00000167258 ENSG00000131748 
            "RPL19"        "FBXL20"          "MED1"         "CDK12"        "STARD3" 
    ENSG00000141736 ENSG00000141741 ENSG00000108344 ENSG00000108342 ENSG00000126351 
            "ERBB2"         "MIEN1"         "PSMD3"          "CSF3"          "THRA" 
    ENSG00000126368 ENSG00000108349 ENSG00000171475 ENSG00000094804 ENSG00000131759 
            "NR1D1"         "CASC3"         "WIPF2"          "CDC6"          "RARA" 
    ENSG00000131747 ENSG00000141753 ENSG00000073584 ENSG00000221852 ENSG00000188581 
            "TOP2A"        "IGFBP4"       "SMARCE1"      "KRTAP1-5"      "KRTAP1-1" 
    ENSG00000212724 ENSG00000213416 ENSG00000131737 ENSG00000171346 ENSG00000171345 
         "KRTAP2-3"     "KRTAP4-12"         "KRT34"         "KRT15"         "KRT19" 
    ENSG00000173812 ENSG00000173805 ENSG00000173801 ENSG00000141756 ENSG00000131473 
             "EIF1"          "HAP1"           "JUP"        "FKBP10"          "ACLY" 
    ENSG00000173786 ENSG00000168259 ENSG00000168256 ENSG00000108773 ENSG00000108774 
              "CNP"        "DNAJC7"       "NKIRAS2"         "KAT2A"         "RAB5C" 
    ENSG00000173757 ENSG00000126561 ENSG00000168610 ENSG00000177469 ENSG00000033627 
           "STAT5B"        "STAT5A"         "STAT3"        "CAVIN1"      "ATP6V0A1" 
    ENSG00000108784 ENSG00000108786 ENSG00000068120 ENSG00000108788 ENSG00000131470 
            "NAGLU"       "HSD17B1"         "COASY"           "MLX"       "PSMC3IP" 
    ENSG00000141699 ENSG00000131462 ENSG00000037042 ENSG00000108797 ENSG00000126562 
          "RETREG3"         "TUBG1"         "TUBG2"       "CNTNAP1"          "WNK4" 
    ENSG00000126581 ENSG00000131480 ENSG00000267060 ENSG00000131469 ENSG00000108828 
            "BECN1"          "AOC2"       "PTGES3L"         "RPL27"          "VAT1" 
    ENSG00000108830 ENSG00000012048 ENSG00000188554 ENSG00000184988 ENSG00000175906 
             "RND2"         "BRCA1"          "NBR1"      "TMEM106A"         "ARL4D" 
    ENSG00000067596 ENSG00000175832 ENSG00000108861 ENSG00000108852 ENSG00000091947 
             "DHX8"          "ETV4"         "DUSP3"          "MPP2"       "TMEM101" 
    ENSG00000141349 ENSG00000108840 ENSG00000125319 ENSG00000087152 ENSG00000013306 
            "G6PC3"         "HDAC5"          "HROB"       "ATXN7L3"      "SLC25A39" 
    ENSG00000030582 ENSG00000161682 ENSG00000180340 ENSG00000161692 ENSG00000182963 
              "GRN"      "FAM171A2"          "FZD2"         "DBF4B"          "GJC1" 
    ENSG00000186185 ENSG00000136448 ENSG00000186834 ENSG00000006062 ENSG00000159314 
           "KIF18B"          "NMT1"        "HEXIM1"       "MAP3K14"      "ARHGAP27" 
    ENSG00000225190 ENSG00000073969 ENSG00000004897 ENSG00000141279 ENSG00000108424 
          "PLEKHM1"           "NSF"         "CDC27"        "NPEPPS"         "KPNB1" 
    ENSG00000006025 ENSG00000189120 ENSG00000167182 ENSG00000108439 ENSG00000005243 
           "OSBPL7"           "SP6"           "SP2"          "PNPO"         "COPZ2" 
    ENSG00000082641 ENSG00000108468 ENSG00000173917 ENSG00000120093 ENSG00000120075 
           "NFE2L1"          "CBX1"         "HOXB2"         "HOXB3"         "HOXB5" 
    ENSG00000260027 ENSG00000136436 ENSG00000159199 ENSG00000159202 ENSG00000159217 
            "HOXB7"      "CALCOCO2"       "ATP5MC1"         "UBE2Z"       "IGF2BP1" 
    ENSG00000167085 ENSG00000121067 ENSG00000121073 ENSG00000005882 ENSG00000167100 
             "PHB1"          "SPOP"       "SLC35B1"          "PDK2"        "SAMD14" 
    ENSG00000108819 ENSG00000108821 ENSG00000154920 ENSG00000108829 ENSG00000167107 
          "PPP1R9B"        "COL1A1"          "EME1"        "LRRC59"         "ACSF2" 
    ENSG00000108846 ENSG00000154945 ENSG00000108848 ENSG00000141232 ENSG00000239672 
            "ABCC3"       "ANKRD40"        "LUC7L3"          "TOB1"          "NME1" 
    ENSG00000011258 ENSG00000011260 ENSG00000141198 ENSG00000108960 ENSG00000141179 
            "MBTD1"         "UTP18"        "TOM1L1"           "MMD"          "PCTP" 
    ENSG00000183691 ENSG00000153933 ENSG00000121060 ENSG00000121064 ENSG00000181610 
              "NOG"          "DGKE"        "TRIM25"        "SCPEP1"        "MRPS23" 
    ENSG00000180891 ENSG00000136451 ENSG00000136450 ENSG00000264364 ENSG00000011143 
           "CUEDC1"         "VEZF1"         "SRSF1"        "DYNLL2"          "MKS1" 
    ENSG00000108389 ENSG00000108384 ENSG00000108395 ENSG00000182628 ENSG00000068489 
            "MTMR4"        "RAD51C"        "TRIM37"          "SKA2"         "PRR11" 
    ENSG00000153982 ENSG00000175155 ENSG00000062716 ENSG00000108443 ENSG00000068097 
            "GDPD1"         "YPEL2"          "VMP1"       "RPS6KB1"        "HEATR6" 
    ENSG00000170832 ENSG00000062725 ENSG00000170836 ENSG00000141376 ENSG00000121068 
            "USP32"        "APPBP2"         "PPM1D"         "BCAS3"          "TBX2" 
    ENSG00000136492 ENSG00000108510 ENSG00000087995 ENSG00000146872 ENSG00000011028 
            "BRIP1"         "MED13"       "METTL2A"          "TLK2"          "MRC2" 
    ENSG00000170921 ENSG00000008283 ENSG00000159640 ENSG00000173826 ENSG00000136485 
            "TANC2"        "CYB561"           "ACE"         "KCNH6"         "DCAF7" 
    ENSG00000198909 ENSG00000136490 ENSG00000108604 ENSG00000178607 ENSG00000136478 
           "MAP3K3"         "LIMD2"       "SMARCD2"          "ERN1"          "TEX2" 
    ENSG00000108654 ENSG00000108854 ENSG00000120063 ENSG00000154240 ENSG00000154229 
             "DDX5"        "SMURF2"         "GNA13"        "CEP112"         "PRKCA" 
    ENSG00000198265 ENSG00000154217 ENSG00000130935 ENSG00000182481 ENSG00000108932 
             "HELZ"       "PITPNC1"         "NOL11"         "KPNA2"       "SLC16A6" 
    ENSG00000070540 ENSG00000108946 ENSG00000154262 ENSG00000125398 ENSG00000133195 
            "WIPI1"       "PRKAR1A"         "ABCA6"          "SOX9"      "SLC39A11" 
    ENSG00000180616 ENSG00000166685 ENSG00000133193 ENSG00000179604 ENSG00000172809 
            "SSTR2"          "COG1"          "VCF1"      "CDC42EP4"         "RPL38" 
    ENSG00000141540 ENSG00000170412 ENSG00000161513 ENSG00000109089 ENSG00000180901 
            "TTYH2"        "GPRC5C"          "FDXR"         "CDR2L"         "KCTD2" 
    ENSG00000189159 ENSG00000125450 ENSG00000125447 ENSG00000125445 ENSG00000177885 
             "JPT1"         "NUP85"          "GGA3"         "MRPS7"          "GRB2" 
    ENSG00000177728 ENSG00000266714 ENSG00000108479 ENSG00000132475 ENSG00000132478 
           "TMEM94"        "MYO15B"         "GALK1"         "H3-3B"           "UNK" 
    ENSG00000132471 ENSG00000132481 ENSG00000141569 ENSG00000188878 ENSG00000161533 
             "WBP2"        "TRIM47"        "TRIM65"          "FBF1"         "ACOX1" 
    ENSG00000167881 ENSG00000186919 ENSG00000182473 ENSG00000185262 ENSG00000161542 
            "SRP68"          "ZACN"         "EXOC7"        "UBALD2"       "PRPSAP1" 
    ENSG00000175931 ENSG00000129667 ENSG00000161544 ENSG00000182534 ENSG00000181038 
            "UBE2O"        "RHBDF2"          "CYGB"         "MXRA7"       "METTL23" 
    ENSG00000161547 ENSG00000092931 ENSG00000129657 ENSG00000184640 ENSG00000078687 
            "SRSF2"        "MFSD11"       "SEC14L1"       "SEPTIN9"        "TNRC6C" 
    ENSG00000108639 ENSG00000167900 ENSG00000183077 ENSG00000089685 ENSG00000184557 
           "SYNGR2"           "TK1"         "AFMID"         "BIRC5"         "SOCS3" 
    ENSG00000087157 ENSG00000108669 ENSG00000055483 ENSG00000108679 ENSG00000171302 
             "PGS1"         "CYTH1"         "USP36"      "LGALS3BP"         "CANT1" 
    ENSG00000173918 ENSG00000167291 ENSG00000141519 ENSG00000171298 ENSG00000181523 
          "C1QTNF1"       "TBC1D16"        "CCDC40"           "GAA"          "SGSH" 
    ENSG00000181045 ENSG00000173818 ENSG00000171246 ENSG00000141564 ENSG00000141577 
         "SLC26A11"         "ENDOV"         "NPTX1"         "RPTOR"        "CEP131" 
    ENSG00000157637 ENSG00000266074 ENSG00000184009 ENSG00000185504 ENSG00000204237 
         "SLC38A10"        "BAHCC1"         "ACTG1"       "FAAP100"         "OXLD1" 
    ENSG00000183048 ENSG00000225663 ENSG00000185624 ENSG00000183684 ENSG00000185813 
         "SLC25A10"        "MCRIP1"          "P4HB"        "ALYREF"         "PCYT2" 
    ENSG00000197063 ENSG00000183010 ENSG00000169750 ENSG00000169727 ENSG00000169710 
             "MAFG"         "PYCR1"          "RAC3"          "GPS1"          "FASN" 
    ENSG00000176155 ENSG00000141526 ENSG00000141551 ENSG00000141574 ENSG00000181396 
           "CCDC57"       "SLC16A3"        "CSNK1D"        "SECTM1"        "OGFOD3" 
    ENSG00000169660 ENSG00000178927 ENSG00000141562 ENSG00000141560 ENSG00000175711 
             "HEXD"         "CYBC1"          "NARF"        "FN3KRP"       "B3GNTL1" 
    ENSG00000176845 ENSG00000101557 ENSG00000079134 ENSG00000158270 ENSG00000176890 
           "METRNL"         "USP14"         "THOC1"       "COLEC12"          "TYMS" 
    ENSG00000176105 ENSG00000080986 ENSG00000101596 ENSG00000132205 ENSG00000101577 
             "YES1"         "NDC80"        "SMCHD1"       "EMILIN2"         "LPIN2" 
    ENSG00000101608 ENSG00000118680 ENSG00000177426 ENSG00000082397 ENSG00000206432 
           "MYL12A"        "MYL12B"         "TGIF1"       "EPB41L3"      "TMEM200C" 
    ENSG00000088756 ENSG00000101680 ENSG00000173482 ENSG00000206418 ENSG00000168502 
         "ARHGAP28"         "LAMA1"         "PTPRM"         "RAB12"         "MTCL1" 
    ENSG00000101745 ENSG00000128791 ENSG00000017797 ENSG00000154845 ENSG00000168461 
          "ANKRD12"         "TWSG1"        "RALBP1"        "PPP4R1"         "RAB31" 
    ENSG00000101558 ENSG00000154856 ENSG00000134265 ENSG00000154864 ENSG00000141404 
             "VAPA"        "APCDD1"          "NAPG"        "PIEZO2"          "GNAL" 
    ENSG00000255112 ENSG00000154889 ENSG00000141401 ENSG00000176014 ENSG00000141391 
           "CHMP1B"         "MPPE1"         "IMPA2"         "TUBB6"      "PRELID3A" 
    ENSG00000134278 ENSG00000101624 ENSG00000175354 ENSG00000085415 ENSG00000101639 
           "SPIRE1"         "CEP76"         "PTPN2"         "SEH1L"        "CEP192" 
    ENSG00000101654 ENSG00000175322 ENSG00000067900 ENSG00000141446 ENSG00000101752 
             "RNMT"        "ZNF519"         "ROCK1"         "ESCO1"          "MIB1" 
    ENSG00000141448 ENSG00000101773 ENSG00000134508 ENSG00000101782 ENSG00000141458 
            "GATA6"         "RBBP8"       "CABLES1"         "RIOK3"          "NPC1" 
    ENSG00000154065 ENSG00000198795 ENSG00000141380 ENSG00000046604 ENSG00000153339 
          "ANKRD29"        "ZNF521"          "SS18"          "DSG2"       "TRAPPC8" 
    ENSG00000134758 ENSG00000141441 ENSG00000141431 ENSG00000134769 ENSG00000166974 
           "RNF138"        "GAREM1"         "ASXL3"          "DTNA"        "MAPRE2" 
    ENSG00000186812 ENSG00000186814 ENSG00000172466 ENSG00000153391 ENSG00000141429 
           "ZNF397"       "ZSCAN30"         "ZNF24"        "INO80C"        "GALNT1" 
    ENSG00000141425 ENSG00000141424 ENSG00000134759 ENSG00000075643 ENSG00000134775 
           "RPRD1A"       "SLC39A6"          "ELP2"         "MOCOS"         "FHOD3" 
    ENSG00000134779 ENSG00000150477 ENSG00000078142 ENSG00000141469 ENSG00000152223 
            "TPGS2"      "KIAA1328"        "PIK3C3"       "SLC14A1"          "EPG5" 
    ENSG00000152229 ENSG00000152234 ENSG00000101638 ENSG00000167220 ENSG00000134049 
          "PSTPIP2"       "ATP5F1A"       "ST8SIA5"         "HDHD2"       "IER3IP1" 
    ENSG00000134030 ENSG00000101665 ENSG00000141627 ENSG00000177576 ENSG00000101670 
             "CTIF"         "SMAD7"           "DYM"      "C18orf32"          "LIPG" 
    ENSG00000167315 ENSG00000154832 ENSG00000154839 ENSG00000082212 ENSG00000176624 
            "ACAA2"         "CXXC1"          "SKA1"           "ME2"         "MEX3C" 
    ENSG00000187323 ENSG00000134046 ENSG00000166845 ENSG00000041353 ENSG00000166510 
              "DCC"          "MBD2"      "C18orf54"        "RAB27B"        "CCDC68" 
    ENSG00000196628 ENSG00000091164 ENSG00000091157 ENSG00000134440 ENSG00000081923 
             "TCF4"         "TXNL1"          "WDR7"         "NARS1"        "ATP8B1" 
    ENSG00000049759 ENSG00000198796 ENSG00000074657 ENSG00000074695 ENSG00000183287 
           "NEDD4L"         "ALPK2"        "ZNF532"         "LMAN1"         "CCBE1" 
    ENSG00000176641 ENSG00000134444 ENSG00000141664 ENSG00000081913 ENSG00000171791 
           "RNF152"         "RELCH"        "ZCCHC2"        "PHLPP1"          "BCL2" 
    ENSG00000119537 ENSG00000119541 ENSG00000166396 ENSG00000166401 ENSG00000171451 
             "KDSR"         "VPS4B"      "SERPINB7"      "SERPINB8"          "DSEL" 
    ENSG00000166479 ENSG00000150636 ENSG00000206052 ENSG00000150637 ENSG00000170677 
             "TMX3"      "CCDC102B"          "DOK6"         "CD226"         "SOCS6" 
    ENSG00000141668 ENSG00000166342 ENSG00000141665 ENSG00000166347 ENSG00000133313 
            "CBLN2"         "NETO1"        "FBXO15"         "CYB5A"         "CNDP2" 
    ENSG00000180011 ENSG00000179981 ENSG00000101493 ENSG00000131196 ENSG00000060069 
            "PTGR3"         "TSHZ1"        "ZNF516"        "NFATC1"         "CTDP1" 
    ENSG00000122490 ENSG00000226742 ENSG00000101546 ENSG00000178184 ENSG00000196476 
          "SLC66A2"       "HSBP1L1"          "RBFA"        "PARD6G"      "C20orf96" 
    ENSG00000247315 ENSG00000177732 ENSG00000125841 ENSG00000101255 ENSG00000125826 
           "ZCCHC3"         "SOX12"         "NRSN2"         "TRIB3"         "RBCK1" 
    ENSG00000125898 ENSG00000125818 ENSG00000101298 ENSG00000088832 ENSG00000088833 
          "FAM110A"         "PSMF1"          "SNPH"        "FKBP1A"        "NSFL1C" 
    ENSG00000198053 ENSG00000125835 ENSG00000215251 ENSG00000088899 ENSG00000198171 
            "SIRPA"         "SNRPB"       "FASTKD5"         "LZTS3"        "DDRGK1" 
    ENSG00000125877 ENSG00000088836 ENSG00000088854 ENSG00000088812 ENSG00000101220 
             "ITPA"       "SLC4A11"        "DNAAF9"          "ATRN"        "ADISSP" 
    ENSG00000101224 ENSG00000088888 ENSG00000125779 ENSG00000101236 ENSG00000088826 
           "CDC25B"          "MAVS"         "PANK2"         "RNF24"          "SMOX" 
    ENSG00000171873 ENSG00000171867 ENSG00000101265 ENSG00000089057 ENSG00000132646 
           "ADRA1D"          "PRNP"        "RASSF2"       "SLC23A2"          "PCNA" 
    ENSG00000101290 ENSG00000125772 ENSG00000125885 ENSG00000125845 ENSG00000132623 
             "CDS2"        "GPCPD1"          "MCM8"          "BMP2"        "ANKEF1" 
    ENSG00000132639 ENSG00000101384 ENSG00000132640 ENSG00000172296 ENSG00000089048 
           "SNAP25"          "JAG1"         "BTBD3"        "SPTLC3"          "ESF1" 
    ENSG00000125848 ENSG00000089177 ENSG00000125868 ENSG00000125844 ENSG00000089006 
            "FLRT3"        "KIF16B"          "DSTN"         "RRBP1"          "SNX5" 
    ENSG00000089091 ENSG00000132664 ENSG00000101310 ENSG00000232388 ENSG00000125821 
           "DZANK1"        "POLR3F"        "SEC23B"        "SMIM26"          "DTD1" 
    ENSG00000132669 ENSG00000173418 ENSG00000101343 ENSG00000088930 ENSG00000125812 
             "RIN2"         "NAA20"        "CRNKL1"          "XRN2"          "GZF1" 
    ENSG00000125814 ENSG00000101439 ENSG00000154930 ENSG00000197586 ENSG00000100994 
             "NAPB"          "CST3"         "ACSS1"        "ENTPD6"          "PYGB" 
    ENSG00000100997 ENSG00000101003 ENSG00000101004 ENSG00000149531 ENSG00000101294 
           "ABHD12"         "GINS1"          "NINL"              NA          "HM13" 
    ENSG00000125968 ENSG00000171552 ENSG00000088325 ENSG00000101306 ENSG00000088356 
              "ID1"        "BCL2L1"          "TPX2"         "MYLK2"         "PDRG1" 
    ENSG00000101337 ENSG00000101346 ENSG00000101350 ENSG00000171456 ENSG00000149600 
           "TM9SF4"        "POFUT1"         "KIF3B"         "ASXL1"        "COMMD7" 
    ENSG00000088305 ENSG00000101367 ENSG00000101391 ENSG00000078699 ENSG00000125967 
           "DNMT3B"        "MAPRE1"      "CDK5RAP1"       "CBFA2T2"        "NECAB3" 
    ENSG00000101412 ENSG00000101417 ENSG00000101421 ENSG00000125970 ENSG00000125977 
             "E2F1"         "PXMP4"        "CHMP4B"          "RALY"        "EIF2S2" 
    ENSG00000101444 ENSG00000078747 ENSG00000125971 ENSG00000101460 ENSG00000198646 
             "AHCY"          "ITCH"       "DYNLRB1"      "MAP1LC3A"         "NCOA6" 
    ENSG00000078804 ENSG00000131069 ENSG00000100983 ENSG00000100991 ENSG00000088298 
         "TP53INP2"         "ACSS2"           "GSS"       "TRPC4AP"         "EDEM2" 
    ENSG00000101000 ENSG00000125966 ENSG00000242372 ENSG00000204183 ENSG00000125965 
            "PROCR"         "MMP24"          "EIF6"              NA          "GDF5" 
    ENSG00000125991 ENSG00000244005 ENSG00000125995 ENSG00000131051 ENSG00000025293 
           "ERGIC3"          "NFS1"         "ROMO1"         "RBM39"         "PHF20" 
    ENSG00000171222 ENSG00000088367 ENSG00000101335 ENSG00000118707 ENSG00000101084 
           "SCAND1"       "EPB41L1"          "MYL9"         "TGIF2"        "RAB5IF" 
    ENSG00000101079 ENSG00000149636 ENSG00000101347 ENSG00000080839 ENSG00000101363 
            "NDRG3"          "DSN1"        "SAMHD1"          "RBL1"        "MANBAL" 
    ENSG00000197122 ENSG00000166619 ENSG00000053438 ENSG00000132792 ENSG00000132821 
              "SRC"         "BLCAP"          "NNAT"       "CTNNBL1"        "VSTM2L" 
    ENSG00000101407 ENSG00000198959 ENSG00000149633 ENSG00000170471 ENSG00000101442 
             "TTI1"          "TGM2"      "KIAA1755"       "RALGAPB"         "ACTR5" 
    ENSG00000101447 ENSG00000204103 ENSG00000198900 ENSG00000124181 ENSG00000174306 
           "FAM83D"          "MAFB"          "TOP1"         "PLCG1"          "ZHX3" 
    ENSG00000124177 ENSG00000124193 ENSG00000101057 ENSG00000124191 ENSG00000149596 
             "CHD6"         "SRSF6"         "MYBL2"          "TOX2"          "JPH2" 
    ENSG00000132823 ENSG00000132824 ENSG00000168734 ENSG00000101109 ENSG00000124134 
            "OSER1"       "SERINC3"          "PKIG"          "STK4"         "KCNS1" 
    ENSG00000124145 ENSG00000204070 ENSG00000124155 ENSG00000101457 ENSG00000175063 
             "SDC4"          "SYS1"          "PIGT"       "DNTTIP1"         "UBE2C" 
    ENSG00000101473 ENSG00000168612 ENSG00000064601 ENSG00000100979 ENSG00000124160 
            "ACOT8"        "ZSWIM1"          "CTSA"          "PLTP"         "NCOA5" 
    ENSG00000080189 ENSG00000172315 ENSG00000197496 ENSG00000064655 ENSG00000124151 
          "SLC35C2"        "TP53RK"       "SLC2A10"          "EYA2"         "NCOA3" 
    ENSG00000196562 ENSG00000124198 ENSG00000124207 ENSG00000124214 ENSG00000124201 
            "SULF2"       "ARFGEF2"         "CSE1L"         "STAU1"         "ZNFX1" 
    ENSG00000124212 ENSG00000158470 ENSG00000197818 ENSG00000158480 ENSG00000124226 
            "PTGIS"       "B4GALT5"        "SLC9A8"        "SPATA2"        "RNF114" 
    ENSG00000124216 ENSG00000240849 ENSG00000172216 ENSG00000196396 ENSG00000042062 
            "SNAI1"         "PEDS1"         "CEBPB"         "PTPN1"        "RIPOR3" 
    ENSG00000124171 ENSG00000124243 ENSG00000101126 ENSG00000101096 ENSG00000054793 
           "PARD6B"         "BCAS4"          "ADNP"        "NFATC2"         "ATP9A" 
    ENSG00000020256 ENSG00000182463 ENSG00000171940 ENSG00000101134 ENSG00000124098 
            "ZFP64"         "TSHZ2"        "ZNF217"          "DOK5"       "FAM210B" 
    ENSG00000087586 ENSG00000101138 ENSG00000087589 ENSG00000022277 ENSG00000087510 
            "AURKA"         "CSTF1"         "CASS4"          "RTF2"        "TFAP2C" 
    ENSG00000124225 ENSG00000124209 ENSG00000124164 ENSG00000198768 ENSG00000124222 
           "PMEPA1"        "RAB22A"          "VAPB"       "APCDD1L"         "STX16" 
    ENSG00000087460 ENSG00000101158 ENSG00000101166 ENSG00000196227 ENSG00000179242 
             "GNAS"        "NELFCD"      "PRELID3B"       "FAM217B"          "CDH4" 
    ENSG00000149657 ENSG00000184402 ENSG00000130702 ENSG00000171858 ENSG00000149679 
           "LSM14B"        "SS18L1"         "LAMA5"         "RPS21"       "CABLES2" 
    ENSG00000101187 ENSG00000060491 ENSG00000101191 ENSG00000101193 ENSG00000101194 
          "SLCO4A1"          "OGFR"         "DIDO1"          "GID8"       "SLC17A9" 
    ENSG00000101199 ENSG00000125534 ENSG00000130589 ENSG00000101216 ENSG00000197457 
          "ARFGAP1"         "PPDPF"         "HELZ2"         "GMEB2"         "STMN3" 
    ENSG00000125520 ENSG00000130584 ENSG00000101152 ENSG00000198276 ENSG00000130590 
         "SLC2A4RG"        "ZBTB46"        "DNAJC5"         "UCKL1"        "SAMD10" 
    ENSG00000101161 ENSG00000171703 ENSG00000171700 ENSG00000203880 ENSG00000141934 
            "PRPF6"         "TCEA2"         "RGS19"        "PCMTD2"         "PLPP2" 
    ENSG00000129946 ENSG00000141933 ENSG00000099804 ENSG00000172270 ENSG00000099822 
             "SHC2"         "TPGS1"         "CDC34"           "BSG"          "HCN2" 
    ENSG00000099821 ENSG00000070423 ENSG00000070404 ENSG00000099864 ENSG00000011304 
           "POLRMT"        "RNF126"         "FSTL3"          "PALM"         "PTBP1" 
    ENSG00000129951 ENSG00000175221 ENSG00000116017 ENSG00000065268 ENSG00000182087 
           "PLPPR3"         "MED16"        "ARID3A"         "WDR18"       "TMEM259" 
    ENSG00000064666 ENSG00000064687 ENSG00000180448 ENSG00000099817 ENSG00000167468 
             "CNN2"         "ABCA7"      "ARHGAP45"        "POLR2E"          "GPX4" 
    ENSG00000064932 ENSG00000118046 ENSG00000167470 ENSG00000228300 ENSG00000160953 
            "SBNO2"         "STK11"          "MIDN"       "FAM174C"        "PWWP3A" 
    ENSG00000130005 ENSG00000115268 ENSG00000119559 ENSG00000185761 ENSG00000181588 
             "GAMT"         "RPS15"      "C19orf25"      "ADAMTSL5"         "MEX3D" 
    ENSG00000071564 ENSG00000129968 ENSG00000227500 ENSG00000213638 ENSG00000133275 
             "TCF3"       "ABHD17A"        "SCAMP4"         "ADAT3"       "CSNK1G2" 
    ENSG00000133243 ENSG00000099875 ENSG00000172081 ENSG00000065000 ENSG00000104885 
            "BTBD2"         "MKNK2"         "MOB3A"         "AP3D1"         "DOT1L" 
    ENSG00000104904 ENSG00000005206 ENSG00000176619 ENSG00000099860 ENSG00000176490 
             "OAZ1"        "SPPL2B"         "LMNB2"       "GADD45B"        "DIRAS1" 
    ENSG00000141873 ENSG00000172009 ENSG00000186300 ENSG00000104964 ENSG00000088256 
          "SLC39A3"         "THOP1"        "ZNF555"          "TLE5"         "GNA11" 
    ENSG00000141905 ENSG00000105325 ENSG00000161091 ENSG00000183397 ENSG00000064961 
             "NFIC"          "FZR1"        "MFSD12"       "TEKTIP1"        "HMG20B" 
    ENSG00000179855 ENSG00000006638 ENSG00000186111 ENSG00000167654 ENSG00000167657 
            "GIPC3"        "TBXA2R"       "PIP5K1C"         "ATCAY"         "DAPK3" 
    ENSG00000105229 ENSG00000178951 ENSG00000126934 ENSG00000105255 ENSG00000167670 
            "PIAS4"        "ZBTB7A"        "MAP2K2"          "FSD1"        "CHAF1A" 
    ENSG00000167671 ENSG00000167674 ENSG00000167680 ENSG00000185361 ENSG00000142002 
            "UBXN6"        "HDGFL2"        "SEMA6B"     "TNFAIP8L1"          "DPP9" 
    ENSG00000141965 ENSG00000127666 ENSG00000105355 ENSG00000276043 ENSG00000105426 
            "FEM1A"        "TICAM1"         "PLIN3"         "UHRF1"         "PTPRS" 
    ENSG00000130255 ENSG00000174917 ENSG00000167733 ENSG00000196365 ENSG00000031823 
            "RPL36"       "MICOS13"      "HSD11B1L"         "LONP1"        "RANBP3" 
    ENSG00000087903 ENSG00000125656 ENSG00000088247 ENSG00000125648 ENSG00000125657 
             "RFX2"          "CLPP"         "KHSRP"      "SLC25A23"        "TNFSF9" 
    ENSG00000125734 ENSG00000130544 ENSG00000171105 ENSG00000198816 ENSG00000090674 
           "GPR108"        "ZNF557"          "INSR"        "ZNF358"        "MCOLN1" 
    ENSG00000032444 ENSG00000181029 ENSG00000142459 ENSG00000178531 ENSG00000104980 
           "PNPLA6"       "TRAPPC5"         "EVI5L"         "CTXN1"        "TIMM44" 
    ENSG00000066044 ENSG00000186994 ENSG00000167772 ENSG00000185236 ENSG00000099783 
           "ELAVL1"         "KANK3"       "ANGPTL4"        "RAB11B"        "HNRNPM" 
    ENSG00000142347 ENSG00000167785 ENSG00000130803 ENSG00000196110 ENSG00000188321 
            "MYO1F"        "ZNF558"        "ZNF317"        "ZNF699"        "ZNF559" 
    ENSG00000130818 ENSG00000197961 ENSG00000171469 ENSG00000171466 ENSG00000198258 
           "ZNF426"        "ZNF121"        "ZNF561"        "ZNF562"          "UBL5" 
    ENSG00000127445 ENSG00000105088 ENSG00000080573 ENSG00000130813 ENSG00000244165 
             "PIN1"         "OLFM2"        "COL5A3"          "SHFL"        "P2RY11" 
    ENSG00000130816 ENSG00000267534 ENSG00000105364 ENSG00000090339 ENSG00000105376 
            "DNMT1"         "S1PR2"         "MRPL4"         "ICAM1"         "ICAM5" 
    ENSG00000105401 ENSG00000065989 ENSG00000079999 ENSG00000180739 ENSG00000130734 
            "CDC37"         "PDE4A"         "KEAP1"         "S1PR5"         "ATG4D" 
    ENSG00000129347 ENSG00000129353 ENSG00000129351 ENSG00000079805 ENSG00000099203 
             "KRI1"       "SLC44A2"          "ILF3"          "DNM2"         "TMED1" 
    ENSG00000130733 ENSG00000130164 ENSG00000161888 ENSG00000130158 ENSG00000105514 
            "YIPF2"          "LDLR"         "SPC24"         "DOCK6"         "RAB3D" 
    ENSG00000105518 ENSG00000130175 ENSG00000130159 ENSG00000130176 ENSG00000198551 
          "TMEM205"        "PRKCSH"         "ECSIT"          "CNN1"        "ZNF627" 
    ENSG00000102575 ENSG00000197044 ENSG00000171295 ENSG00000198429 ENSG00000197647 
             "ACP5"        "ZNF441"        "ZNF440"         "ZNF69"        "ZNF433" 
    ENSG00000223547 ENSG00000196646 ENSG00000197857 ENSG00000173875 ENSG00000104774 
           "ZNF844"        "ZNF136"         "ZNF44"        "ZNF791"        "MAN2B1" 
    ENSG00000105583 ENSG00000095059 ENSG00000132004 ENSG00000105576 ENSG00000198356 
          "WDR83OS"          "DHPS"         "FBXW9"         "TNPO2"          "GET3" 
    ENSG00000171223 ENSG00000167815 ENSG00000104889 ENSG00000105612 ENSG00000179115 
             "JUNB"         "PRDX2"      "RNASEH2A"        "DNASE2"         "FARSA" 
    ENSG00000179218 ENSG00000179262 ENSG00000179271 ENSG00000008441 ENSG00000104903 
             "CALR"        "RAD23A"    "GADD45GIP1"          "NFIX"          "LYL1" 
    ENSG00000160888 ENSG00000141837 ENSG00000104957 ENSG00000037757 ENSG00000132003 
             "IER2"       "CACNA1A"         "YJU2B"          "MRI1"        "ZSWIM4" 
    ENSG00000132024 ENSG00000132017 ENSG00000104998 ENSG00000187867 ENSG00000141858 
           "CC2D1A"        "DCAF15"        "IL27RA"         "PALM3"         "SAMD1" 
    ENSG00000072062 ENSG00000105011 ENSG00000072071 ENSG00000123146 ENSG00000123136 
           "PRKACA"         "ASF1B"        "ADGRL1"        "ADGRE5"        "DDX39A" 
    ENSG00000123143 ENSG00000123159 ENSG00000099797 ENSG00000127507 ENSG00000141867 
             "PKN1"         "GIPC1"          "TECR"        "ADGRE2"          "BRD4" 
    ENSG00000105127 ENSG00000011451 ENSG00000167460 ENSG00000167461 ENSG00000105058 
            "AKAP8"           "WIZ"          "TPM4"         "RAB8A"        "FAM32A" 
    ENSG00000072958 ENSG00000127528 ENSG00000127527 ENSG00000127526 ENSG00000214046 
            "AP1M1"          "KLF2"       "EPS15L1"       "SLC35E1"         "SMIM7" 
    ENSG00000131351 ENSG00000099331 ENSG00000099330 ENSG00000160113 ENSG00000130312 
            "HAUS8"         "MYO9B"         "OCEL1"         "NR2F6"        "MRPL34" 
    ENSG00000130311 ENSG00000074855 ENSG00000130299 ENSG00000130304 ENSG00000167483 
             "DDA1"          "ANO8"        "GTPBP3"       "SLC27A1"        "NIBAN3" 
    ENSG00000130477 ENSG00000130479 ENSG00000105640 ENSG00000007080 ENSG00000105643 
           "UNC13A"         "MAP1S"        "RPL18A"       "CCDC124"        "ARRDC2" 
    ENSG00000254858 ENSG00000105649 ENSG00000105650 ENSG00000130522 ENSG00000130520 
          "MPV17L2"         "RAB3A"         "PDE4C"          "JUND"          "LSM4" 
    ENSG00000130517 ENSG00000130513 ENSG00000130511 ENSG00000105655 ENSG00000105656 
           "PGPEP1"         "GDF15"         "SSBP4"        "ISYNA1"           "ELL" 
    ENSG00000105701 ENSG00000006016 ENSG00000105662 ENSG00000105664 ENSG00000223802 
            "FKBP8"         "CRLF1"         "CRTC1"          "COMP"         "CERS1" 
    ENSG00000105676 ENSG00000254901 ENSG00000064490 ENSG00000129933 ENSG00000167491 
            "ARMC6"        "BORCS8"        "RFXANK"          "MAU2"       "GATAD2A" 
    ENSG00000160161 ENSG00000105717 ENSG00000064547 ENSG00000181896 ENSG00000256771 
            "CILP2"          "PBX4"         "LPAR2"        "ZNF101"        "ZNF253" 
    ENSG00000184635 ENSG00000256229 ENSG00000237440 ENSG00000188171 ENSG00000160352 
            "ZNF93"        "ZNF486"        "ZNF737"        "ZNF626"        "ZNF714" 
    ENSG00000196705 ENSG00000182141 ENSG00000172687 ENSG00000196268 ENSG00000197020 
           "ZNF431"        "ZNF708"        "ZNF738"        "ZNF493"        "ZNF100" 
    ENSG00000198521 ENSG00000160321 ENSG00000167232 ENSG00000213096 ENSG00000105171 
            "ZNF43"        "ZNF208"         "ZNF91"        "ZNF254"          "POP4" 
    ENSG00000131943 ENSG00000105176 ENSG00000168813 ENSG00000178904 ENSG00000105185 
         "C19orf12"          "URI1"        "ZNF507"       "DPY19L3"         "PDCD5" 
    ENSG00000105186 ENSG00000213965 ENSG00000131944 ENSG00000131941 ENSG00000076650 
          "ANKRD27"        "NUDT19"        "FAAP24"         "RHPN2"       "GPATCH1" 
    ENSG00000130881 ENSG00000153879 ENSG00000124299 ENSG00000153885 ENSG00000166398 
             "LRP3"         "CEBPG"          "PEPD"        "KCTD15"        "GARRE1" 
    ENSG00000105220 ENSG00000126261 ENSG00000089335 ENSG00000180884 ENSG00000089351 
              "GPI"          "UBA2"        "ZNF302"        "ZNF792"       "GRAMD1A" 
    ENSG00000089327 ENSG00000105698 ENSG00000189001 ENSG00000105677 ENSG00000249115 
            "FXYD5"          "USF2"          "SBSN"       "TMEM147"         "HAUS5" 
    ENSG00000126254 ENSG00000126267 ENSG00000205155 ENSG00000004776 ENSG00000004777 
            "RBM42"        "COX6B1"        "PSENEN"         "HSPB6"      "ARHGAP33" 
    ENSG00000105290 ENSG00000105270 ENSG00000075702 ENSG00000126247 ENSG00000167635 
            "APLP1"         "CLIP3"         "WDR62"        "CAPNS1"        "ZNF146" 
    ENSG00000142065 ENSG00000181007 ENSG00000186017 ENSG00000254004 ENSG00000186020 
            "ZFP14"         "ZFP82"        "ZNF566"        "ZNF260"        "ZNF529" 
    ENSG00000267041 ENSG00000197863 ENSG00000251247 ENSG00000198453 ENSG00000196967 
           "ZNF850"        "ZNF790"        "ZNF345"        "ZNF568"       "ZNF585A" 
    ENSG00000245680 ENSG00000188283 ENSG00000181666 ENSG00000196437 ENSG00000171827 
          "ZNF585B"        "ZNF383"        "ZNF875"        "ZNF569"        "ZNF570" 
    ENSG00000188227 ENSG00000120784 ENSG00000011332 ENSG00000167645 ENSG00000099341 
           "ZNF793"         "ZFP30"          "DPF1"         "YIF1B"         "PSMD8" 
    ENSG00000130402 ENSG00000104823 ENSG00000104824 ENSG00000068903 ENSG00000104825 
            "ACTN4"          "ECH1"        "HNRNPL"         "SIRT2"        "NFKBIB" 
    ENSG00000269190 ENSG00000161243 ENSG00000128011 ENSG00000130755 ENSG00000179134 
           "FBXO17"        "FBXO27"         "LRFN1"          "GMFG"        "SAMD4B" 
    ENSG00000063322 ENSG00000128016 ENSG00000090924 ENSG00000105193 ENSG00000196235 
            "MED29"         "ZFP36"       "PLEKHG2"         "RPS16"        "SUPT5H" 
    ENSG00000105197 ENSG00000105204 ENSG00000105202 ENSG00000275395 ENSG00000128000 
           "TIMM50"        "DYRK1B"           "FBL"         "FCGBP"       "ZNF780B" 
    ENSG00000130758 ENSG00000105221 ENSG00000160392 ENSG00000105223 ENSG00000197019 
          "MAP3K10"          "AKT2"      "C19orf47"          "PLD3"       "SERTAD1" 
    ENSG00000167565 ENSG00000090013 ENSG00000160460 ENSG00000160410 ENSG00000090006 
          "SERTAD3"         "BLVRB"        "SPTBN4"        "SHKBP1"         "LTBP4" 
    ENSG00000086544 ENSG00000188493 ENSG00000077312 ENSG00000167578 ENSG00000269858 
            "ITPKC"        "ACTMAP"         "SNRPA"         "RAB4B"         "EGLN2" 
    ENSG00000167600 ENSG00000167601 ENSG00000105323 ENSG00000105329 ENSG00000142039 
           "CYP2S1"           "AXL"      "HNRNPUL1"         "TGFB1"        "CCDC97" 
    ENSG00000123810 ENSG00000076928 ENSG00000028277 ENSG00000160570 ENSG00000105723 
             "B9D2"       "ARHGEF1"        "POU2F2"         "DEDD2"         "GSK3A" 
    ENSG00000079462 ENSG00000105429 ENSG00000204941 ENSG00000243137 ENSG00000124466 
         "PAFAH1B3"         "MEGF8"          "PSG5"          "PSG4"         "LYPD3" 
    ENSG00000176531 ENSG00000105755 ENSG00000176472 ENSG00000234465 ENSG00000167378 
           "PHLDB3"         "ETHE1"        "ZNF575"        "PINLYP"          "IRGQ" 
    ENSG00000124444 ENSG00000131116 ENSG00000105767 ENSG00000104783 ENSG00000167637 
           "ZNF576"        "ZNF428"         "CADM4"         "KCNN4"        "ZNF283" 
    ENSG00000176222 ENSG00000124459 ENSG00000159905 ENSG00000159882 ENSG00000267680 
           "ZNF404"         "ZNF45"        "ZNF221"        "ZNF230"        "ZNF224" 
    ENSG00000256294 ENSG00000263002 ENSG00000167380 ENSG00000131115 ENSG00000159917 
           "ZNF225"        "ZNF234"        "ZNF226"        "ZNF227"        "ZNF235" 
    ENSG00000278318 ENSG00000073008 ENSG00000187244 ENSG00000130202 ENSG00000130204 
           "ZNF229"           "PVR"          "BCAM"       "NECTIN2"        "TOMM40" 
    ENSG00000104853 ENSG00000104859 ENSG00000142252 ENSG00000104884 ENSG00000104881 
           "CLPTM1"        "CLASRP"        "GEMIN7"         "ERCC2"      "PPP1R13L" 
    ENSG00000012061 ENSG00000125740 ENSG00000125753 ENSG00000125741 ENSG00000125746 
            "ERCC1"          "FOSB"          "VASP"          "OPA3"          "EML2" 
    ENSG00000125743 ENSG00000177051 ENSG00000177045 ENSG00000176182 ENSG00000104983 
           "SNRPD2"        "FBXO46"          "SIX5"         "MYPOP"        "CCDC61" 
    ENSG00000011485 ENSG00000169515 ENSG00000182013 ENSG00000160014 ENSG00000197380 
            "PPP5C"         "CCDC8"        "PNMA8A"         "CALM3"         "DACT3" 
    ENSG00000105287 ENSG00000181027 ENSG00000105281 ENSG00000042753 ENSG00000160007 
            "PRKD2"          "FKRP"        "SLC1A5"         "AP2S1"      "ARHGAP35" 
    ENSG00000142230 ENSG00000105327 ENSG00000257704 ENSG00000105419 ENSG00000105402 
             "SAE1"          "BBC3"        "INAFM1"         "MEIS3"          "NAPA" 
    ENSG00000024422 ENSG00000105373 ENSG00000178980 ENSG00000105486 ENSG00000105483 
             "EHD2"         "NOP53"       "SELENOW"          "LIG1"         "CARD8" 
    ENSG00000142227 ENSG00000105438 ENSG00000105447 ENSG00000105443 ENSG00000063180 
             "EMP3"        "KDELR1"         "GRWD1"         "CYTH2"          "CA11" 
    ENSG00000176909 ENSG00000105552 ENSG00000087076 ENSG00000087074 ENSG00000104805 
           "MAMSTR"         "BCAT2"      "HSD17B14"      "PPP1R15A"         "NUCB1" 
    ENSG00000087088 ENSG00000087086 ENSG00000104812 ENSG00000183207 ENSG00000130529 
              "BAX"           "FTL"          "GYS1"        "RUVBL2"         "TRPM4" 
    ENSG00000063127 ENSG00000104894 ENSG00000074219 ENSG00000142541 ENSG00000104870 
          "SLC6A16"          "CD37"         "TEAD2"        "RPL13A"         "FCGRT" 
    ENSG00000142552 ENSG00000126464 ENSG00000126458 ENSG00000126461 ENSG00000126456 
             "RCN3"         "PRR12"          "RRAS"         "SCAF1"          "IRF3" 
    ENSG00000126453 ENSG00000126457 ENSG00000196961 ENSG00000104973 ENSG00000104960 
          "BCL2L12"         "PRMT1"         "AP2A1"         "MED25"         "PTOV1" 
    ENSG00000039650 ENSG00000204673 ENSG00000104946 ENSG00000213024 ENSG00000105053 
             "PNKP"        "AKT1S1"       "TBC1D17"         "NUP62"          "VRK3" 
    ENSG00000142528 ENSG00000131408 ENSG00000062822 ENSG00000161671 ENSG00000161677 
           "ZNF473"         "NR1H2"         "POLD1"         "EMC10"         "JOSD2" 
    ENSG00000161681 ENSG00000105472 ENSG00000167747 ENSG00000186806 ENSG00000105379 
           "SHANK1"       "CLEC11A"              NA       "VSIG10L"          "ETFB" 
    ENSG00000105497 ENSG00000161551 ENSG00000142556 ENSG00000256087 ENSG00000197608 
           "ZNF175"        "ZNF577"        "ZNF614"        "ZNF432"        "ZNF841" 
    ENSG00000204611 ENSG00000196267 ENSG00000105568 ENSG00000196214 ENSG00000198464 
           "ZNF616"        "ZNF836"       "PPP2R1A"        "ZNF766"        "ZNF480" 
    ENSG00000167554 ENSG00000221923 ENSG00000167555 ENSG00000258405 ENSG00000198482 
           "ZNF610"        "ZNF880"        "ZNF528"        "ZNF578"        "ZNF808" 
    ENSG00000167562 ENSG00000167766 ENSG00000213020 ENSG00000189190 ENSG00000198538 
           "ZNF701"         "ZNF83"        "ZNF611"        "ZNF600"         "ZNF28" 
    ENSG00000204604 ENSG00000182986 ENSG00000180257 ENSG00000170949 ENSG00000197937 
           "ZNF468"        "ZNF320"        "ZNF816"        "ZNF160"        "ZNF347" 
    ENSG00000197928 ENSG00000196417 ENSG00000179820 ENSG00000105605 ENSG00000170906 
           "ZNF677"        "ZNF765"         "MYADM"        "CACNG7"        "NDUFA3" 
    ENSG00000105619 ENSG00000105618 ENSG00000125505 ENSG00000170892 ENSG00000170889 
             "TFPT"        "PRPF31"        "MBOAT7"        "TSEN34"          "RPS9" 
    ENSG00000167615 ENSG00000275183 ENSG00000160439 ENSG00000131037 ENSG00000125503 
            "LENG8"         "LENG9"         "RDH13"        "EPS8L1"      "PPP1R12C" 
    ENSG00000080031 ENSG00000160469 ENSG00000133247 ENSG00000095752 ENSG00000108107 
            "PTPRH"         "BRSK1"         "KMT5C"          "IL11"         "RPL28" 
    ENSG00000108106 ENSG00000187902 ENSG00000090971 ENSG00000179954 ENSG00000171443 
            "UBE2S"        "SHISA7"         "NAT14"         "SSC5D"        "ZNF524" 
    ENSG00000261221 ENSG00000213015 ENSG00000171425 ENSG00000173581 ENSG00000063244 
           "ZNF865"        "ZNF580"        "ZNF581"       "CCDC106"         "U2AF2" 
    ENSG00000063245 ENSG00000142409 ENSG00000018869 ENSG00000198440 ENSG00000198046 
             "EPN1"        "ZNF787"        "ZNF582"        "ZNF583"        "ZNF667" 
    ENSG00000197016 ENSG00000197951 ENSG00000083844 ENSG00000204524 ENSG00000178229 
           "ZNF470"         "ZNF71"        "ZNF264"        "ZNF805"        "ZNF543" 
    ENSG00000131845 ENSG00000186230 ENSG00000178201 ENSG00000121406 ENSG00000251369 
           "ZNF304"        "ZNF749"         "VN1R1"        "ZNF549"        "ZNF550" 
    ENSG00000171649 ENSG00000213762 ENSG00000121417 ENSG00000179909 ENSG00000152443 
             "ZIK1"        "ZNF134"        "ZNF211"        "ZNF154"        "ZNF776" 
    ENSG00000269343 ENSG00000204514 ENSG00000177025 ENSG00000166704 ENSG00000121413 
          "ZNF587B"        "ZNF814"      "C19orf18"        "ZNF606"       "ZSCAN18" 
    ENSG00000181894 ENSG00000171606 ENSG00000198131 ENSG00000278129 ENSG00000182318 
           "ZNF329"        "ZNF274"        "ZNF544"          "ZNF8"       "ZSCAN22" 
    ENSG00000083845 ENSG00000171574 ENSG00000131849 ENSG00000083812 ENSG00000119574 
             "RPS5"        "ZNF584"        "ZNF132"        "ZNF324"        "ZBTB45" 
    ENSG00000130724 ENSG00000130725 ENSG00000099326 ENSG00000183307 ENSG00000131100 
           "CHMP2A"         "UBE2M"          "MZF1"      "TMEM121B"      "ATP6V1E1" 
    ENSG00000015475 ENSG00000243156 ENSG00000215193 ENSG00000184979 ENSG00000070413 
              "BID"        "MICAL3"         "PEX26"         "USP18"         "DGCR2" 
    ENSG00000100056 ENSG00000070371 ENSG00000185608 ENSG00000093009 ENSG00000093010 
             "ESS2"        "CLTCL1"        "MRPL40"         "CDC45"          "COMT" 
    ENSG00000099889 ENSG00000099901 ENSG00000185252 ENSG00000099910 ENSG00000099917 
            "ARVCF"        "RANBP1"         "ZNF74"        "KLHL22"         "MED15" 
    ENSG00000241973 ENSG00000099940 ENSG00000099942 ENSG00000128228 ENSG00000100023 
            "PI4KA"        "SNAP29"          "CRKL"        "SDF2L1"         "PPIL2" 
    ENSG00000100027 ENSG00000100030 ENSG00000100034 ENSG00000128266 ENSG00000100228 
            "YPEL1"         "MAPK1"         "PPM1F"          "GNAZ"         "RAB36" 
    ENSG00000186716 ENSG00000250479 ENSG00000099953 ENSG00000099991 ENSG00000099994 
              "BCR"       "CHCHD10"         "MMP11"        "CABIN1"         "SUSD2" 
    ENSG00000099998 ENSG00000100014 ENSG00000100028 ENSG00000100099 ENSG00000100109 
             "GGT5"       "SPECC1L"        "SNRPD3"          "HPS4"        "TFIP11" 
    ENSG00000128294 ENSG00000180957 ENSG00000100154 ENSG00000183765 ENSG00000100209 
            "TPST2"        "PITPNB"         "TTC28"         "CHEK2"          "HSCB" 
    ENSG00000159873 ENSG00000100219 ENSG00000183579 ENSG00000186998 ENSG00000100263 
          "CCDC117"          "XBP1"         "ZNRF3"         "EMID1"        "RHBDD3" 
    ENSG00000182944 ENSG00000185340 ENSG00000100280 ENSG00000100285 ENSG00000184117 
            "EWSR1"        "GAS2L1"         "AP1B1"          "NEFH"      "NIPSNAP1" 
    ENSG00000186575 ENSG00000100330 ENSG00000128342 ENSG00000099995 ENSG00000187860 
              "NF2"         "MTMR3"           "LIF"         "SF3A1"       "CCDC157" 
    ENSG00000099999 ENSG00000100003 ENSG00000185339 ENSG00000167065 ENSG00000184792 
           "RNF215"       "SEC14L2"          "TCN2"        "DUSP18"         "OSBP2" 
    ENSG00000133422 ENSG00000183963 ENSG00000185133 ENSG00000138942 ENSG00000100100 
            "MORC2"          "SMTN"        "INPP5J"        "RNF185"       "PIK3IP1" 
    ENSG00000185721 ENSG00000198089 ENSG00000241878 ENSG00000183530 ENSG00000128245 
             "DRG1"          "SFI1"          "PISD"        "PRR14L"         "YWHAH" 
    ENSG00000100220 ENSG00000100225 ENSG00000100234 ENSG00000133424 ENSG00000100281 
             "RTCB"         "FBXO7"         "TIMP3"        "LARGE1"        "HMGXB4" 
    ENSG00000100284 ENSG00000100297 ENSG00000221963 ENSG00000100320 ENSG00000128335 
             "TOM1"          "MCM5"         "APOL6"        "RBFOX2"         "APOL2" 
    ENSG00000100345 ENSG00000100350 ENSG00000100353 ENSG00000100360 ENSG00000133466 
             "MYH9"       "FOXRED2"         "EIF3D"         "IFT27"       "C1QTNF6" 
    ENSG00000128340 ENSG00000100065 ENSG00000128283 ENSG00000100092 ENSG00000100097 
             "RAC2"        "CARD10"      "CDC42EP1"        "SH3BP1"        "LGALS1" 
    ENSG00000100106 ENSG00000100139 ENSG00000185022 ENSG00000198792 ENSG00000213923 
           "TRIOBP"       "MICALL1"          "MAFF"      "TMEM184B"        "CSNK1E" 
    ENSG00000168135 ENSG00000100196 ENSG00000100201 ENSG00000184949 ENSG00000100221 
            "KCNJ4"        "KDELR3"         "DDX17"       "FAM227A"         "JOSD1" 
    ENSG00000100226 ENSG00000100242 ENSG00000100246 ENSG00000179750 ENSG00000244509 
           "GTPBP1"          "SUN2"         "DNAL4"      "APOBEC3B"      "APOBEC3C" 
    ENSG00000128394 ENSG00000100311 ENSG00000128268 ENSG00000100335 ENSG00000128272 
         "APOBEC3F"         "PDGFB"         "MGAT3"         "MIEF1"          "ATF4" 
    ENSG00000187051 ENSG00000100354 ENSG00000100359 ENSG00000128285 ENSG00000100380 
         "RPS19BP1"        "TNRC6B"         "SGSM3"         "MCHR1"          "ST13" 
    ENSG00000196236 ENSG00000100387 ENSG00000100393 ENSG00000100401 ENSG00000100403 
          "XPNPEP3"          "RBX1"         "EP300"       "RANGAP1"        "ZC3H7B" 
    ENSG00000183864 ENSG00000100412 ENSG00000172346 ENSG00000196419 ENSG00000100138 
             "TOB2"          "ACO2"         "CSDC2"         "XRCC6"         "SNU13" 
    ENSG00000100147 ENSG00000198911 ENSG00000100162 ENSG00000198951 ENSG00000177096 
          "CCDC134"        "SREBF2"         "CENPM"          "NAGA"        "PHETA2" 
    ENSG00000184983 ENSG00000189306 ENSG00000100227 ENSG00000128274 ENSG00000242247 
           "NDUFA6"         "RRP7A"       "POLDIP3"        "A4GALT"       "ARFGAP3" 
    ENSG00000100266 ENSG00000100300 ENSG00000188677 ENSG00000093000 ENSG00000100364 
          "PACSIN2"          "TSPO"         "PARVB"         "NUP50"      "KIAA0930" 
    ENSG00000100376 ENSG00000077942 ENSG00000130638 ENSG00000188064 ENSG00000075218 
          "FAM118A"         "FBLN1"        "ATXN10"         "WNT7B"         "GTSE1" 
    ENSG00000075240 ENSG00000100422 ENSG00000054611 ENSG00000100425 ENSG00000100426 
           "GRAMD4"          "CERK"      "TBC1D22A"          "BRD1"         "ZBED4" 
    ENSG00000182858 ENSG00000184164 ENSG00000198355 ENSG00000128159 ENSG00000188130 
            "ALG12"        "CRELD2"          "PIM3"       "TUBGCP6"        "MAPK12" 
    ENSG00000196576 ENSG00000205593 ENSG00000100241 ENSG00000128165 ENSG00000100258 
           "PLXNB2"       "DENND6B"          "SBF1"          "ADM2"          "LMF2" 
    ENSG00000025770 ENSG00000100299 ENSG00000079974 ENSG00000277117 ENSG00000180530 
           "NCAPH2"          "ARSA"        "RABL2B"  "LOC102723996"         "NRIP1" 
    ENSG00000154639 ENSG00000154640 ENSG00000154642 ENSG00000154654 ENSG00000154719 
            "CXADR"          "BTG3"      "C21orf91"         "NCAM2"        "MRPL39" 
    ENSG00000154721 ENSG00000154723 ENSG00000154727 ENSG00000142192 ENSG00000154734 
             "JAM2"        "ATP5PF"         "GABPA"           "APP"       "ADAMTS1" 
    ENSG00000154736 ENSG00000198862 ENSG00000156253 ENSG00000156265 ENSG00000156273 
          "ADAMTS5"          "LTN1"        "RWDD2B"      "MAP3K7CL"         "BACH1" 
    ENSG00000156304 ENSG00000142149 ENSG00000159055 ENSG00000159086 ENSG00000159128 
            "SCAF4"          "HUNK"        "MIS18A"        "PAXBP1"        "IFNGR2" 
    ENSG00000142188 ENSG00000159131 ENSG00000159140 ENSG00000159147 ENSG00000205726 
          "TMEM50B"          "GART"           "SON"        "DONSON"         "ITSN1" 
    ENSG00000241837 ENSG00000243927 ENSG00000198743 ENSG00000159200 ENSG00000159216 
           "ATP5PO"         "MRPS6"        "SLC5A3"         "RCAN1"         "RUNX1" 
    ENSG00000159228 ENSG00000159231 ENSG00000142197 ENSG00000159256 ENSG00000159259 
             "CBR1"          "CBR3"         "DOP1B"         "MORC3"        "CHAF1B" 
    ENSG00000159263 ENSG00000159267 ENSG00000182670 ENSG00000157538 ENSG00000157540 
             "SIM2"          "HLCS"          "TTC3"        "VPS26C"        "DYRK1A" 
    ENSG00000157542 ENSG00000157551 ENSG00000157557 ENSG00000183527 ENSG00000185658 
            "KCNJ6"        "KCNJ15"          "ETS2"         "PSMG1"         "BRWD1" 
    ENSG00000182093 ENSG00000171587 ENSG00000182240 ENSG00000157617 ENSG00000173276 
             "GET1"         "DSCAM"         "BACE2"         "C2CD2"        "ZBTB21" 
    ENSG00000160179 ENSG00000160190 ENSG00000160200 ENSG00000142178 ENSG00000160208 
            "ABCG1"       "SLC37A1"           "CBS"          "SIK1"         "RRP1B" 
    ENSG00000160209 ENSG00000160216 ENSG00000141959 ENSG00000160233 ENSG00000184787 
             "PDXK"        "AGPAT3"          "PFKL"         "LRRC3"        "UBE2G2" 
    ENSG00000184900 ENSG00000160256 ENSG00000197381 ENSG00000186866 ENSG00000182871 
            "SUMO3"          "SLX9"        "ADARB1"        "POFUT2"       "COL18A1" 
    ENSG00000173638 ENSG00000142156 ENSG00000142173 ENSG00000160284 ENSG00000160285 
          "SLC19A1"        "COL6A1"        "COL6A2"       "SPATC1L"           "LSS" 
    ENSG00000160294 ENSG00000160298 ENSG00000160299 ENSG00000198888 ENSG00000198763 
           "MCM3AP"      "C21orf58"          "PCNT"           "ND1"           "ND2" 
    ENSG00000198804 ENSG00000228253 ENSG00000198899 ENSG00000198938 ENSG00000198840 
             "COX1"          "ATP8"          "ATP6"          "COX3"           "ND3" 
    ENSG00000212907 ENSG00000198886 ENSG00000198786 ENSG00000198695 ENSG00000198727 
             "ND4L"           "ND4"           "ND5"           "ND6"          "CYTB" 
    ENSG00000273748 ENSG00000271254 
                 NA  "LOC124905564" 

``` r
write.table(sig_genes, file="significant_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
```

![](R-HSA-68886.png)
