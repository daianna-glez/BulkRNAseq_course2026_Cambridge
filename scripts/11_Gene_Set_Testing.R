
# ------------------------------------------------------------------------------
#                             11. Gene Set Testing
# ------------------------------------------------------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                                  1.1 KEGG ORA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library(tidyverse) # for data handling
library(clusterProfiler) # for accessing KEGG database and conduct ORA and GSEA

## Search for mouse data in KEGG
search_kegg_organism('mouse', by='common_name')
## Use the ‘mmu’ ‘kegg_code’.

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

## Subset to genes with valid (!NA), padj<0.05, |logFC|>1 and valid Entrez ID
sigGenes <- shrink.d11 %>%
    drop_na(Entrez, padj) %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    pull(Entrez)

length(sigGenes)
# [1] 629

## Run ORA with KEGG mmu pathways
keggRes <- enrichKEGG(gene = sigGenes, organism = 'mmu')
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

## Look at one pathway: ‘Antigen processing and presentation’
#       Column	Description
#           ID	KEGG pathway ID (e.g., hsa04110)
#   Description	Pathway name (e.g., Cell cycle)
#     GeneRatio	intersection size / num of genes of interest (with EntrezID and annotated in KEGG database)
#       BgRatio	intersection size against all background genes (annotated in KEGG database)
#         pvalue	Raw enrichment p-value (hypergeometric test)
#       p.adjust	Adjusted p-value (e.g., BH/FDR)
#         qvalue	q-value (estimated FDR)
#         geneID	List of overlapping genes (separated by /)
#          Count	Intersection size

as.data.frame(keggRes)["mmu04612",]
#                    category   subcategory       ID                         Description GeneRatio  BgRatio RichFactor FoldEnrichment   zScore
# mmu04612 Organismal Systems Immune system mmu04612 Antigen processing and presentation    40/354 88/10650  0.4545455       13.67488 22.13778
#                pvalue     p.adjust       qvalue
# mmu04612 3.583239e-36 4.371551e-34 3.074042e-34
#           geneID
# mmu04612 14991/15519/19186/12265/21356/15001/21355/21926/16149/15015/100504404/14960/21354/15000/14998/213233/13040/12526/15042/12525/15043/15978/15039/15018/14969/12010/14972/15040/15007/667977/14963/110557/14964/14961/19188/100529082/15006/14999/15013/65972
#          Count
# mmu04612    40


## Visualize pathway pathway ‘mmu04612’ in browser:
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
         limit = list(gene=20, cpd=1)) # control color scale limits for genes and metabolites


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
library(org.Mm.eg.db) # for searching gene IDs in mouse

## Find Ensembl IDs first
sigGenes_GO <-  shrink.d11 %>%
    drop_na(padj) %>%
    filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
    pull(GeneID)

universe <- shrink.d11$GeneID
length(universe)
# [1] 20091
head(universe)
# [1] "ENSMUSG00000000001" "ENSMUSG00000000028" "ENSMUSG00000000037" "ENSMUSG00000000049"
# [5] "ENSMUSG00000000056" "ENSMUSG00000000058"

## Run ORA
ego <- enrichGO(gene          = sigGenes_GO, # receives ensembl IDs
                universe      = universe,
                OrgDb         = org.Mm.eg.db,
                keyType       = "ENSEMBL",
                ont           = "BP",
                pvalueCutoff  = 0.01,
                readable      = TRUE)

## Visualize results
barplot(ego, showCategory=20) # show intersection size and pval for top 20 sets/terms

dotplot(ego, font.size = 14) # shows count, gene ratio and pval

library(enrichplot)
ego_pt <- pairwise_termsim(ego) # shows the overlap between genes across different GO terms.
emapplot(ego_pt)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                                  2. GSEA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library(msigdb) # connect to Molecular Signatures Database (MSigDB) gene sets
## MSigDB data is stored inside ExperimentHub:
library(ExperimentHub) # access available data

# Create an ExperimentHub object to connect to the ExperimentHub online metadata database
# It uses local cache directory on your machine to store temporarily the downloaded files
eh = ExperimentHub()
# Searches the ExperimentHub catalog for datasets whose metadata matches all of these keywords:
# - "msigdb" → datasets related to Molecular Signatures Database
# - "mm" → Mus musculus (mouse)
# - "2023" → the MSigDB release year/version
query(eh , c('msigdb', 'mm', '2023'))

# ExperimentHub with 10 records
# snapshotDate(): 2025-10-07
# $dataprovider: Broad Institute, EBI
# $species: Mus musculus, Homo sapiens
# $rdataclass: GSEABase::GeneSetCollection, data.frame
# additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer,
#   rdatadateadded, preparerclass, tags, rdatapath, sourceurl, sourcetype
# retrieve records with, e.g., 'object[["EH8285"]]'

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

## Download most recent available release using Entrez IDs
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

## Subset to h collection: human hallmark gene sets (well annotated biological states/processes)
hallmarks = subsetCollection(msigdb.mm, 'h')
## Extract gene IDs per set
msigdb_ids = geneIds(hallmarks)
class(msigdb_ids)
# [1] "list"
length(msigdb_ids)
# [1] 50
names(msigdb_ids)

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
# - ranked genes
# - pathways
# - gene set minimum size
# - gene set maximum size

gseaRes <- GSEA(rankedGenes,
                TERM2GENE = term2gene,
                pvalueCutoff = 1.00,
                minGSSize = 15,
                maxGSSize = 500)

head(gseaRes, 1)

## Top 10 enriched pathways
as_tibble(gseaRes) %>%
    arrange(desc(abs(NES))) %>%
    top_n(10, wt=-p.adjust) %>%  # order by -p.adjust
    dplyr::select(-core_enrichment) %>%
    mutate(across(c("enrichmentScore", "NES"), ~round(.x, digits=3))) %>%
    mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))

# A tibble: 10 × 10
# ID             Description setSize enrichmentScore   NES pvalue p.adjust qvalue  rank leading_edge
# <chr>          <chr>         <int>           <dbl> <dbl> <chr>  <chr>    <chr>  <dbl> <chr>
#     1 HALLMARK_INTE… HALLMARK_I…     152           0.953  1.40 1.00e… 1.67e-09 1.30e…   722 tags=75%, l…
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
         geneSetID = "HALLMARK_INFLAMMATORY_RESPONSE",
         title = "HALLMARK_INFLAMMATORY_RESPONSE")

###############
# Exercise 2: Rank the genes by statistical significance and regulation direction
#             (-log10(pvalue) * sign(logFC)).
#             Run GSEA using the new ranked genes and the H pathways.
#             Conduct the same analysis for the day 33 Infected vs Uninfected contrast.
###############

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

