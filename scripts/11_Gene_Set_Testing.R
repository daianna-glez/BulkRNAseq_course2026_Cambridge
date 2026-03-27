
# ------------------------------------------------------------------------------
#                             11. Gene Set Testing
# ------------------------------------------------------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                                  1.1 KEGG ORA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library(tidyverse) # for data handling
library(clusterProfiler) # for accessing KEGG database and conduct ORA and GSEA

## Search the code for "mouse" in KEGG by "common name"
search_kegg_organism('mouse', by='common_name')
## Use the ‘mmu’ ‘kegg_code’.

## Search by scientific name
search_kegg_organism('Mus musculus', by='scientific_name')

## Case sensitive!!
search_kegg_organism('mus musculus', by='scientific_name')
# <0 rows> (or 0-length row.names)

## Partial name is supported
search_kegg_organism('muscu', by='scientific_name')


## Input: list of DEGs for Infected vs Uninfected at d11
shrink.d11 <- readRDS("RObjects/Shrunk_Results.d11.rds")
dim(shrink.d11)
# [1] 20091    14
head(shrink.d11, 2)
#               GeneID   baseMean log2FoldChange      lfcSE    pvalue      padj Entrez Symbol
# 1 ENSMUSG00000000001 1102.56094    -0.02799080 0.06758558 0.2453209 0.6876811  14679  Gnai3
# 2 ENSMUSG00000000028   58.60055     0.02898123 0.12631602 0.2772868 0.7191570  12544  Cdc45
#                                                                                            Description        Biotype
# 1 guanine nucleotide binding protein (G protein), alpha inhibiting 3 [Source:MGI Symbol;Acc:MGI:95773] protein_coding
# 2                                           cell division cycle 45 [Source:MGI Symbol;Acc:MGI:1338073] protein_coding
#   Chr     Start       End Strand
# 1   3 108107280 108146146     -1
# 2  16  18780447  18811987     -1

## Note there are NAs! (genes filtered by DESeq2 without padj but with LFC and SE)
table(is.na(shrink.d11$GeneID))
table(is.na(shrink.d11$baseMean))
table(is.na(shrink.d11$log2FoldChange))
table(is.na(shrink.d11$lfcSE))
table(is.na(shrink.d11$padj))
# FALSE  TRUE
# 17707  2384

## We need Entrez IDs for ORA
table(is.na(shrink.d11$Entrez))
# FALSE  TRUE
# 17275  2816

## Subset to genes with valid padj<0.05, |logFC|>1, and valid Entrez ID
sigGenes <- shrink.d11 %>%
    drop_na(Entrez, padj) %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    pull(Entrez)

length(sigGenes)
# [1] 629

## Run ORA with KEGG mmu pathways
keggRes <- enrichKEGG(gene = sigGenes, organism = 'mmu')
dim(keggRes)
# [1] 78 14

## Explore
as_tibble(keggRes)

## Each row gives the results of enrichment for a biological set
keggRes %>% head(1)
#                category                subcategory        ID                       Description  GeneRatio    BgRatio
# mmu05168 Human Diseases  Infectious disease: viral  mmu05168  Herpes simplex virus 1 infection     64/354  211/10650
#          RichFactor FoldEnrichment   zScore       pvalue     p.adjust       qvalue
# mmu05168  0.3033175       9.125231 22.10392 1.211682e-44 2.956505e-42 2.078992e-42
# geneID
# mmu05168 16391/16160/14991/19106/12266/21356/15001/21355/72512/21926/16149/54123/12370/20846/20684/71586/16176/24088/23961/78781/11796/17874/246727/246728/20304/20293/20296/15015/100504404/14960/18854/21354/15000/14998/213233/20847/230073/230979/56489/81897/69550/246730/12702/15042/15043/15978/15039/15018/14969/12010/14972/170741/23960/15040/15007/667977/14963/110557/14964/14961/100529082/15006/14999/15013
#          Count
# mmu05168    64

#       Column	Description
#           ID	KEGG pathway ID (e.g., hsa04110)
#   Description	Pathway name (e.g., Cell cycle)
#     GeneRatio	intersection size / num of DEGs (annotated in KEGG database)
#       BgRatio	size of pathway in KEGG / gene universe
#         pvalue	Raw enrichment p-value (hypergeometric test)
#       p.adjust	Adjusted p-value (e.g., BH/FDR)
#         qvalue	q-value (estimated FDR)
#         geneID	List of overlapping genes (separated by /)
#          Count	Intersection size


## Only signif results returned
max(keggRes$p.adjust)
# [1] 0.04992357


## Look at one pathway: ‘Antigen processing and presentation’ (overlapping genes in red)
as.data.frame(keggRes)[1,]
#                category               subcategory       ID                      Description GeneRatio   BgRatio
# mmu05168 Human Diseases Infectious disease: viral mmu05168 Herpes simplex virus 1 infection    64/356 209/11156
#          RichFactor FoldEnrichment   zScore       pvalue     p.adjust       qvalue
# mmu05168  0.3062201       9.596043 22.77572 5.497077e-46 1.346784e-43 9.489691e-44
#          geneID
# mmu05168 16391/16160/14991/19106/12266/21356/15001/21355/72512/21926/16149/54123/12370/20846/20684/71586/16176/24088/23961/78781/11796/17874/246727/246728/20304/20293/20296/15015/100504404/14960/18854/21354/15000/14998/213233/20847/230073/230979/56489/81897/69550/246730/12702/15042/15043/15978/15039/15018/14969/12010/14972/170741/23960/15040/15007/667977/14963/110557/14964/14961/100529082/15006/14999/15013
#          Count
# mmu05168    64


## Or by pathway ID
as.data.frame(keggRes)["mmu05168",]


## Visualize pathway pathway ‘mmu04612’ in web browser:
browseKEGG(keggRes, 'mmu04612')

## Visualise a pathway as a file
library(pathview) # for generating  figures of KEGG pathways

## With pathview we can colour genes by logFC:
## pass pathview a named vector of fold change values
logFC <- shrink.d11$log2FoldChange
names(logFC) <- shrink.d11$Entrez
head(logFC)
#        14679        12544       107815        11818        67608        12390
# -0.027990799  0.028981226 -0.009722409  0.031648710  0.031538172 -0.074402546
summary(logFC)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -2.092e+01 -2.102e-02  2.540e-05  8.964e-02  2.019e-02  8.439e+00

## This will export the pathway as a png
pathview(gene.data = logFC,
         pathway.id = "mmu04612",
         species = "mmu",
         limit = list(gene=5, cpd=1)) # control color scale limits for genes and metabolites


###############
# Exercise 1: Use pathview to export a figure for “mmu04659” or “mmu04658”,
##            but this time only use genes with padj < 0.01
###############
logFC <- shrink.d11 %>%
    drop_na(padj, Entrez) %>%
    filter(padj < 0.01) %>%
    pull(log2FoldChange, Entrez)

length(logFC)
# [1] 744
pathview(gene.data = logFC,
         pathway.id = "mmu04659",
         species = "mmu",
         limit = list(gene=5, cpd=1))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                                  1.2 GO ORA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## ORA with GO terms: clusterProfiler can also perform ORA on GO terms using enrichGO()
library(org.Mm.eg.db) # for searching gene IDs in mouse (it provides annotation data for mouse genes)

## Extract DEGs first (use Ensembl IDs)
sigGenes_GO <-  shrink.d11 %>%
    drop_na(padj) %>%
    filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
    pull(GeneID)

## Define universe (it should exclude NA padj since no testing was conducted)
universe <- shrink.d11$GeneID
length(universe)
# [1] 20091
head(universe)
# [1] "ENSMUSG00000000001" "ENSMUSG00000000028" "ENSMUSG00000000037" "ENSMUSG00000000049"
# [5] "ENSMUSG00000000056" "ENSMUSG00000000058"

## Run ORA
ego <- enrichGO(gene          = sigGenes_GO, # receives ensembl IDs of DEGs
                universe      = universe,
                OrgDb         = org.Mm.eg.db, # use mouse annotations
                keyType       = "ENSEMBL", # search according to ensembl IDs
                ont           = "BP", # ORA on BP gene sets
                pvalueCutoff  = 0.01,
                readable      = TRUE)

## Visualize results
barplot(ego, showCategory=20) # shows intersection size and pval for top 20 sets/terms

## Note these DEGs are genes responding to viral infection,
## consistent with diff expression 11 days post infection
## We are looking at the genes that are responsing to infection

dotplot(ego, font.size = 14) # shows count, gene ratio and pval

## Shows the overlap between genes across different GO terms.
library(enrichplot)
ego_pt <- pairwise_termsim(ego)
## Two big clusters of BP affected by Infection after 11 days6
emapplot(ego_pt)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                                  2. GSEA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library(msigdb) # connect to Molecular Signatures Database (MSigDB)
## MSigDB data is stored inside ExperimentHub:
library(ExperimentHub) ## to access available curated data sets

# Create an ExperimentHub object to connect to the ExperimentHub database
eh = ExperimentHub()
# 1. Loads the ExperimentHub interface
# 2. Connects to the online ExperimentHub repository
# 3. Retrieves metadata about all available datasets
# 4. Sets up a local cache on your machine (so downloads are saved on your machine temporarily)

# eh becomes a special object that acts like a catalog of datasets

# Search in the ExperimentHub catalog for datasets matching keywords for:
# - "msigdb" → datasets related to Molecular Signatures Database
query(eh , c('msigdb'))
# - "mm" → Mus musculus (mouse)
query(eh , c('msigdb', 'mm'))
# - "2023" → the MSigDB release year/version
query(eh , c('msigdb', 'mm', '2023'))

# title
# EH8285 | msigdb.v2022.1.mm.EZID
# EH8286 | msigdb.v2022.1.mm.idf
# EH8287 | msigdb.v2022.1.mm.SYM
# EH8291 | msigdb.v2023.1.mm.EZID <- use this one
# EH8292 | msigdb.v2023.1.mm.idf
# EH8293 | msigdb.v2023.1.mm.SYM
# EH8297 | msigdb.v7.5.1.mm.EZID
# EH8298 | msigdb.v7.5.1.mm.idf
# EH8299 | msigdb.v7.5.1.mm.SYM
# EH8300 | imex_hsmm_0722

## Download most recent available release using Entrez IDs with getMsigdb()
## getMsigdb() will:
#  1. connects to ExperimentHub
#  2. finds the correct dataset
#  3. downloads it (if needed)
#  4. formats it nicely for you

msigdb.mm <- getMsigdb(org = 'mm', id = 'EZID', version = '2023.1')
listCollections(msigdb.mm)
# [1] "c1" "c3" "c2" "c8" "c6" "c7" "c4" "c5" "h"

## Perform GSEA using clusterProfiler
## Rank genes by shrunk logFC
rankedGenes <- shrink.d11 %>%
    drop_na(GeneID, padj, log2FoldChange) %>%
    mutate(rank = log2FoldChange) %>%
    filter(!is.na(Entrez)) %>%
    arrange(desc(rank)) %>%
    pull(rank, Entrez)

head(rankedGenes)
#    15945    24108   626578    20210    17329    64380
# 8.439020 8.307955 7.818902 7.787435 7.783766 7.525960
tail(rankedGenes)
#     76757     545279     231382  100039192      19109  100503353
# -2.980815  -3.094735  -3.456050  -5.949883 -17.173829 -18.129005

## Subset to h collection: hallmark gene sets (well annotated biological states/processes)
hallmarks = subsetCollection(msigdb.mm, 'h')
## Extract gene IDs per set
msigdb_ids = geneIds(hallmarks)
class(msigdb_ids)
# [1] "list"
length(msigdb_ids)
# [1] 50 pathways
names(msigdb_ids)

## Look at one set of list
msigdb_ids$HALLMARK_ADIPOGENESIS %>% head()

## Num of genes in each set
lapply(msigdb_ids, length) %>% unlist() %>% table()
# 51  57  64  68  71  74  79  93 102 122 126 143 172 175 186 188 196 205 231 234 236 246 251 273 278
# 1   1   1   1   1   1   1   1   1   1   1   1   1   2   1   1   1   1   1   1   1   1   1   1   1
# 280 296 300 304 310 311 317 325 327 328 332 338 348 356 362 363 364 365 375 549
# 1   1   1   2   1   1   2   1   1   3   1   1   1   1   1   1   1   1   1   1

## Convert list to df with columns set - EntrezID
term2gene <- enframe(msigdb_ids, name = "gs_name", value = "entrez") %>%
             unnest(entrez)
head(term2gene)
# A tibble: 6 × 2
#  gs_name               entrez
#  <chr>                 <chr>
# 1 HALLMARK_ADIPOGENESIS 11303
# 2 HALLMARK_ADIPOGENESIS 11304
# 3 HALLMARK_ADIPOGENESIS 27403
# 4 HALLMARK_ADIPOGENESIS 268379
# 5 HALLMARK_ADIPOGENESIS 74591
# 6 HALLMARK_ADIPOGENESIS 381072

## Conduct GSEA

## Input:
# - ranked genes and logFCs
# - pathways and genes in each
# - gene set minimum size
# - gene set maximum size

gseaRes <- GSEA(rankedGenes,
                TERM2GENE = term2gene,
                pvalueCutoff = 1.00,
                minGSSize = 15,
                maxGSSize = 500)

## Explore results
head(gseaRes, 1)
## Definitions:
# - NES = ES / mean(ES from permutations of that gene set) to normalize ES by gene set size and scale
#         and make it comparable across different pathways and experiments
# - rank: position in the ranked gene list where the enrichment score (ES) reaches its maximum.
# - leading_edge: summary of enrichment signal region:
#       * tags: % of the genes in this gene set appear before the ES peak
#       * list: % of the entire ranked gene list you had to scan to reach the ES peak
#       * signal: combined measure of enrichment strength (based on tags + list)
# - core_enrichment: group of genes in the pathway that contributes most to the enrichment signal

dim(gseaRes)
# [1] 49 11

## Didn't test pathways with less than 15 genes or more than 500
lapply(msigdb_ids, length) %>% unlist() %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 51.0   172.8   275.5   244.5   327.8   549.0

## How mant with <15 or >500
table(lapply(msigdb_ids, length) %>% unlist() > 500)
# FALSE  TRUE
#    49     1

## Not all significant
as.data.frame(gseaRes) %>%
         dplyr::select(-core_enrichment) %>%
    pull(p.adjust) %>% max()
# [1] 0.998999


## Top 10 enriched pathways
as.data.frame(gseaRes) %>%
    dplyr::select(-core_enrichment)  %>%
    arrange(desc(abs(NES))) %>%
    top_n(10, wt=-p.adjust) %>%  # order by -p.adjust
    mutate(across(c("enrichmentScore", "NES"), ~round(.x, digits=3))) %>%
    mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))

# A tibble: 10 × 10
# ID             Description setSize enrichmentScore   NES pvalue p.adjust qvalue  rank leading_edge
# <chr>          <chr>         <int>           <dbl> <dbl> <chr>  <chr>    <chr>  <dbl> <chr>
# 1 HALLMARK_INTE… HALLMARK_I…     152           0.953  1.40 1.00e… 1.67e-09 1.30e…   722 tags=75%, l…
# 2 HALLMARK_INTE… HALLMARK_I…     281           0.946  1.39 1.00e… 1.67e-09 1.30e…   842 tags=64%, l…
# 3 HALLMARK_ALLO… HALLMARK_A…     283           0.927  1.36 1.00e… 1.67e-09 1.30e…   811 tags=44%, l…
# 4 HALLMARK_IL6_… HALLMARK_I…     107           0.922  1.36 2.79e… 1.99e-06 1.55e…   806 tags=42%, l…
# 5 HALLMARK_INFL… HALLMARK_I…     288           0.89   1.31 1.62e… 2.02e-08 1.58e…  1023 tags=32%, l…
# 6 HALLMARK_IL2_… HALLMARK_I…     295           0.887  1.30 5.51e… 5.51e-08 4.29e…   786 tags=20%, l…
# 7 HALLMARK_TNFA… HALLMARK_T…     262           0.884  1.3  2.18e… 1.82e-07 1.41e…  1047 tags=29%, l…
# 8 HALLMARK_APOP… HALLMARK_A…     219           0.872  1.28 8.64e… 4.80e-05 3.74e…  1047 tags=21%, l…
# 9 HALLMARK_COMP… HALLMARK_C…     288           0.866  1.27 1.54e… 9.60e-06 7.48e…  1059 tags=26%, l…
# 10 HALLMARK_COAG… HALLMARK_C…     187           0.862  1.27 1.80e… 9.02e-04 7.02e…  1033 tags=22%, l…



## Enrichment score plot:
# - displays gene ranking and walking sum
# - genes in set as black ticks (no tick for genes not in set)
# - the walking sum: the green curve
# - the enrichment score: dotted red line

gseaplot(gseaRes,
         geneSetID = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
         title = "HALLMARK_INTERFERON_ALPHA_RESPONSE")

## Look at one example with negative NES
as.data.frame(gseaRes) %>%
    dplyr::select(-core_enrichment) %>%
    arrange(desc((NES))) %>%
    tail

gseaplot(gseaRes,
         geneSetID = "HALLMARK_PANCREAS_BETA_CELLS",
         title = "HALLMARK_PANCREAS_BETA_CELLS")


## Look at one example with ns ES
as.data.frame(gseaRes) %>%
    dplyr::select(-core_enrichment) %>%
    arrange(desc((p.adjust))) %>%
    head(10)

gseaplot(gseaRes,
         geneSetID = "HALLMARK_APICAL_SURFACE",
         title = "HALLMARK_APICAL_SURFACE")


## Loot at example with pathways genes in both extremes
as.data.frame(gseaRes) %>%
    dplyr::select(-core_enrichment) %>%
    arrange(-desc((NES))) %>%
    head(10)

gseaplot(gseaRes,
         geneSetID = "HALLMARK_HEDGEHOG_SIGNALING",
         title = "HALLMARK_HEDGEHOG_SIGNALING")



###############
# Exercise 2: Rank the genes by statistical significance and regulation direction
#             (-log10(pvalue) * sign(logFC)).
#             Run GSEA using the new ranked genes and the H pathways.
#             Conduct the same analysis for the day 33 Infected vs Uninfected contrast.
###############

## Why we use -log10(p)?

## Pvalues go from 0 to 1 so look very similar between them and its hard to compare them
sort(shrink.d11$pvalue) %>% head
# [1] 8.173101e-120 4.235016e-112 1.339305e-111 4.598980e-111  4.469568e-97
# [6]  8.481895e-96

sort(shrink.d11$pvalue) %>% head %>% log10
# [1] -119.08761 -111.37314 -110.87312 -110.33734  -96.34973  -95.07151

hist(shrink.d11$pvalue)
hist(log10(shrink.d11$pvalue))
hist(-log10(shrink.d11$pvalue))

## Rank genes
rankedGenes.e11 <- shrink.d11 %>%
    drop_na(Entrez, pvalue, log2FoldChange) %>%
    mutate(rank = -log10(pvalue) * sign(log2FoldChange)) %>%
    arrange(desc(rank)) %>%
    pull(rank, Entrez)

## Conduct analysis:
gseaRes.e11 <- GSEA(rankedGenes.e11,
                    TERM2GENE = term2gene,
                    pvalueCutoff = 1.00,
                    minGSSize = 15,
                    maxGSSize = 500)

## View results
as_tibble(gseaRes.e11) %>%
    arrange(desc(abs(NES))) %>%
    top_n(10, wt=-p.adjust) %>%
    dplyr::select(-core_enrichment)

# ID          Description setSize enrichmentScore   NES   pvalue p.adjust  qvalue  rank leading_edge
# <chr>       <chr>         <int>           <dbl> <dbl>    <dbl>    <dbl>   <dbl> <dbl> <chr>
#     1 HALLMARK_O… HALLMARK_O…     216          -0.477 -1.97 4.39e- 8  4.39e-7 3.33e-7  4323 tags=52%, l…
# 2 HALLMARK_I… HALLMARK_I…     160           0.954  1.96 1   e-10  1.67e-9 1.26e-9   620 tags=69%, l…
# 3 HALLMARK_I… HALLMARK_I…     289           0.943  1.94 1   e-10  1.67e-9 1.26e-9   649 tags=58%, l…
# 4 HALLMARK_A… HALLMARK_A…     293           0.909  1.87 1   e-10  1.67e-9 1.26e-9   709 tags=39%, l…
# 5 HALLMARK_I… HALLMARK_I…     112           0.869  1.78 8.11e- 7  5.79e-6 4.39e-6  1024 tags=44%, l…
# 6 HALLMARK_I… HALLMARK_I…     303           0.820  1.68 1.16e- 9  1.45e-8 1.10e-8   718 tags=19%, l…
# 7 HALLMARK_I… HALLMARK_I…     298           0.802  1.65 9.41e- 8  7.84e-7 5.95e-7  1193 tags=32%, l…
# 8 HALLMARK_T… HALLMARK_T…     265           0.795  1.63 1.07e- 6  6.71e-6 5.09e-6  1621 tags=37%, l…
# 9 HALLMARK_C… HALLMARK_C…     301           0.785  1.61 1.78e- 6  9.87e-6 7.48e-6  1294 tags=26%, l…
# 10 HALLMARK_U… HALLMARK_U…     291           0.768  1.58 1.92e- 5  9.59e-5 7.27e-5  1193 tags=16%, l…

## GSEA for Infected vs Uninfected DEGs at d33
shrink.d33 <- readRDS("RObjects/Shrunk_Results.d33.rds")
dim(shrink.d33)
head(shrink.d33)

## Rank genes
rankedGenes.e33 <- shrink.d33 %>%
    drop_na(Entrez, pvalue, log2FoldChange) %>%
    mutate(rank = -log10(pvalue) * sign(log2FoldChange)) %>%
    arrange(desc(rank)) %>%
    pull(rank, Entrez)

## Perform analysis
gseaRes.e33 <- GSEA(rankedGenes.e33,
                    TERM2GENE = term2gene,
                    pvalueCutoff = 1.00,
                    minGSSize = 15,
                    maxGSSize = 500)

## View the results
as_tibble(gseaRes.e33) %>%
    arrange(desc(abs(NES))) %>%
    top_n(10, wt=-p.adjust) %>%
    dplyr::select(-core_enrichment)

