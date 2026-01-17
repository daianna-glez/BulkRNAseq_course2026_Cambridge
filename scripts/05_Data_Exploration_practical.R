
################################################################################
#               Introduction to Bulk RNAseq data analysis
################################################################################
# ------------------------------------------------------------------------------
#                         5. RNA-seq Data Exploration
# ------------------------------------------------------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                           1. Count data import
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library(tximport) # for importing and summarizing transcript-level estimates
library(DESeq2) # for data transformation
library(tidyverse) # for data handling

## Read sample metadata into a data frame
sampleinfo <- read_tsv("data/samplesheet.tsv", col_types = c("cccc"))
arrange(sampleinfo, Status, TimePoint, Replicate)
# A tibble: 12 × 4
#    SampleName Replicate Status     TimePoint
#    <chr>      <chr>     <chr>      <chr>
# 1  SRR7657878 1         Infected   d11
# 2  SRR7657881 2         Infected   d11
# 3  SRR7657880 3         Infected   d11
# 4  SRR7657874 1         Infected   d33
# 5  SRR7657882 2         Infected   d33
# 6  SRR7657872 3         Infected   d33
# 7  SRR7657877 1         Uninfected d11
# 8  SRR7657876 2         Uninfected d11
# 9  SRR7657879 3         Uninfected d11
# 10 SRR7657883 1         Uninfected d33
# 11 SRR7657873 2         Uninfected d33
# 12 SRR7657875 3         Uninfected d33

## Read count data from Salmon
files <- file.path("salmon", sampleinfo$SampleName, "quant.sf")
files
# [1] "salmon/SRR7657878/quant.sf" "salmon/SRR7657881/quant.sf" "salmon/SRR7657880/quant.sf"
# [4] "salmon/SRR7657874/quant.sf" "salmon/SRR7657882/quant.sf" "salmon/SRR7657872/quant.sf"
# [7] "salmon/SRR7657877/quant.sf" "salmon/SRR7657876/quant.sf" "salmon/SRR7657879/quant.sf"
# [10] "salmon/SRR7657883/quant.sf" "salmon/SRR7657873/quant.sf" "salmon/SRR7657875/quant.sf"

## Use sample names as column names in count matrix
files <- set_names(files, sampleinfo$SampleName)
# SRR7657878                   SRR7657881                   SRR7657880
# "salmon/SRR7657878/quant.sf" "salmon/SRR7657881/quant.sf" "salmon/SRR7657880/quant.sf"
# SRR7657874                   SRR7657882                   SRR7657872
# "salmon/SRR7657874/quant.sf" "salmon/SRR7657882/quant.sf" "salmon/SRR7657872/quant.sf"
# SRR7657877                   SRR7657876                   SRR7657879
# "salmon/SRR7657877/quant.sf" "salmon/SRR7657876/quant.sf" "salmon/SRR7657879/quant.sf"
# SRR7657883                   SRR7657873                   SRR7657875
# "salmon/SRR7657883/quant.sf" "salmon/SRR7657873/quant.sf" "salmon/SRR7657875/quant.sf"

## Read tx to gene info
tx2gene <- read_tsv("references/tx2gene.tsv")
head(tx2gene)
#   TxID                 GeneID
#   <chr>                <chr>
# 1 ENSMUST00000177564.1 ENSMUSG00000096176
# 2 ENSMUST00000196221.1 ENSMUSG00000096749
# 3 ENSMUST00000179664.1 ENSMUSG00000096749
# 4 ENSMUST00000178537.1 ENSMUSG00000095668
# 5 ENSMUST00000178862.1 ENSMUSG00000094569
# 6 ENSMUST00000179520.1 ENSMUSG00000094028

## tximport() imports salmon's tx level estimates and summarizes counts at gene level
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
class(txi)
# [1] "list"
names(txi)
# [1] "abundance"           "counts"              "length"              "countsFromAbundance"

lapply(txi, head, 2)

## Abundance = TPM (lib size normalized counts)
head(txi$abundance, 2)
#                    SRR7657878 SRR7657881 SRR7657880 SRR7657874 SRR7657882 SRR7657872 SRR7657877 SRR7657876
# ENSMUSG00000000001   20.39354   20.35246   19.70855   21.09177   19.32488   27.21929    25.4685   20.88858
# ENSMUSG00000000003    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000     0.0000    0.00000
#                    SRR7657879 SRR7657883 SRR7657873 SRR7657875
# ENSMUSG00000000001   21.36722   24.05217   23.42266   19.75996
# ENSMUSG00000000003    0.00000    0.00000    0.00000    0.00000

## Estimated gene count-like data from: Σ(TPM × tx length × library size) over all gene txs
head(txi$counts, 2)
#                    SRR7657878 SRR7657881 SRR7657880 SRR7657874 SRR7657882 SRR7657872 SRR7657877 SRR7657876
# ENSMUSG00000000001       1039   1005.889        892     917.36    1136.69       1259   1351.221   1110.999
# ENSMUSG00000000003          0      0.000          0       0.00       0.00          0      0.000      0.000
#                    SRR7657879 SRR7657883 SRR7657873 SRR7657875
# ENSMUSG00000000001   1067.634   1134.522   1272.003       1065
# ENSMUSG00000000003      0.000      0.000      0.000          0


## Gene length per sample is a weighted average of its tx: Σ(TPM x tx length)/Σ(TPM) over all gene txs
## For a given gene, the average tx length may vary between samples if different samples are using alternative transcripts.
head(txi$length, 2)
#                    SRR7657878 SRR7657881 SRR7657880 SRR7657874 SRR7657882 SRR7657872 SRR7657877 SRR7657876
# ENSMUSG00000000001  2903.3180  2900.7180  2898.0260  2874.3030  2874.4990  2880.3300  2929.1690  2879.8660
# ENSMUSG00000000003   540.9167   540.9167   540.9167   540.9167   540.9167   540.9167   540.9167   540.9167
#                    SRR7657879 SRR7657883 SRR7657873 SRR7657875
# ENSMUSG00000000001  2868.4940  2862.7890  2897.4510  2879.7850
# ENSMUSG00000000003   540.9167   540.9167   540.9167   540.9167

## If estimated counts were generated using scaled abundances
txi$countsFromAbundance
# [1] "no"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                           2. Gene filtering
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Prepare count matrix
rawCounts <- round(txi$counts, 0)

## Check dimension of count matrix
dim(rawCounts)
# [1] 35896    12

## For each gene, compute total count and compare to threshold = 5
keep <- rowSums(rawCounts) > 5
## Summary
table(keep, useNA = "always")
# FALSE  TRUE  <NA>
# 15805 20091     0

## Subset to kept genes
filtCounts <- rawCounts[keep,]
## Check dimension of new count matrix
dim(filtCounts)
# [1] 20091    12


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#               3. Count distribution and Data transformations
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Issues with raw read counts:

## 1. Wide range with few outliers
apply(filtCounts, 2, range)
#      SRR7657878 SRR7657881 SRR7657880 SRR7657874 SRR7657882 SRR7657872 SRR7657877 SRR7657876 SRR7657879 SRR7657883
# [1,]          0          0          0          0          0          0          0          0          0          0
# [2,]     652318     590723     435516     444448     699333     418060     613859     757858     722648     652247
#      SRR7657873 SRR7657875
# [1,]          0          0
# [2,]     616071     625800

boxplot(filtCounts, main = 'Raw counts', las = 2)

## 2. Variance increases with mean gene expression
plot(rowMeans(filtCounts), rowSds(filtCounts),
     main = 'Raw counts: sd vs mean',
     xlim = c(0, 10000),
     ylim = c(0, 5000))

## Log2 transformation
logCounts <- log2(filtCounts + 1)
boxplot(logCounts, main = 'Log2 counts', las = 2)

## Mean-variance relationship is reduced but still there (inverted) !!!
plot(rowMeans(logCounts), rowSds(logCounts),
     main = 'Log2(counts): sd vs mean')

## Rlog: log2-transformation accounting for lib size and mean-variance relationship
rlogcounts <- rlog(filtCounts)
boxplot(rlogcounts, main = 'rlog counts', las = 2)
plot(rowMeans(rlogcounts), rowSds(rlogcounts),
     main = 'Rlog(counts): sd vs mean')


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                        4. Principal Component Analysis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library(ggfortify) # to visualize statistical result
library(ggrepel) # for plot labeling

## Useful for:
# 1. Identifying main sources of gene expr variation (e.g. experimental condition or batch effects)
# 2. Identifying sample outliers
# 3. Identifying sample swaps

## We have to use lognorm data
rlogcounts <- rlog(filtCounts)
## Run PCA
pcDat <- prcomp(t(rlogcounts))

## PCs are in x
pcDat$x

## Plot PCA
autoplot(pcDat)

## Color by status and shape by time
autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

## Label samples to identify mislabelled ones
autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) +
    geom_text_repel(aes(x = PC1, y = PC2, label = SampleName), box.padding = 0.8)

## Swap samples
sampleinfo <- mutate(sampleinfo,
                     Status = case_when(
                         SampleName=="SRR7657882" ~ "Uninfected",
                         SampleName=="SRR7657873" ~ "Infected",
                         TRUE ~ Status)
                     )
write_tsv(sampleinfo, "results/SampleInfo_Corrected.txt")

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

#_______________________________________________________________________________
# Exercise:

## Plot PC2 vs PC3
autoplot(pcDat,
         data = sampleinfo,
         colour="Status",
         shape="TimePoint",
         x=2,
         y=3,
         size=5)
#_______________________________________________________________________________


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                         5. Hierachical clustering
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library(ggdendro) # to plot clustering results

## Useful for:
# - Identifying sample relationships based on their Euclidean distance: Σ(gᵢ₂-gᵢ₁)²

hclDat <-  t(rlogcounts) %>%
    dist(method = "euclidean") %>%
    hclust()
ggdendrogram(hclDat, rotate = TRUE)

## Add sample info
hclDat2 <- hclDat
hclDat2$labels <- str_c(sampleinfo$Status, ":", sampleinfo$TimePoint)
ggdendrogram(hclDat2, rotate = TRUE)

## Same conclusion!





