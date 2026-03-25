

################################################################################
##           Differential gene expression analysis with DESeq2
################################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  1. Load R packages
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library(DESeq2)
library(tidyverse)
library(limma)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  2. Load data for differential gene expression analysis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## - Gene expression count matrix
txi <- readRDS("RObjects/txi.rds")
## - Sample metadata
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv", col_types = "cccc")

## Quick explorations
dim(txi$counts)
# [1] 35896    12
dim(sampleinfo)
# [1] 12  4

## Check sample colnames in matrix correspond to rownames in metadata
all(colnames(txi$counts) == sampleinfo$SampleName)
# [1] TRUE

## Variables in the data
head(sampleinfo)
#   SampleName Replicate Status     TimePoint
#   <chr>      <chr>     <chr>      <chr>
# 1 SRR7657878 1         Infected   d11
# 2 SRR7657881 2         Infected   d11
# 3 SRR7657880 3         Infected   d11
# 4 SRR7657874 1         Infected   d33
# 5 SRR7657882 2         Uninfected d33
# 6 SRR7657872 3         Infected   d33


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  3. Creating the model formulas:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# --------------------- Simple model: Gene expr ~ TimePoint --------------------
## Define formula
simple.model <- as.formula(~ TimePoint)

## Build model matrix from formula:

# - Intercept added, set to 1 for all samples
# - TimePoint is now a dummy variable indicating whether samples are d11 (0) or d33 (1)
# - Reference/baseline for TimePoint is d11

model.matrix(simple.model, data = sampleinfo)
#    (Intercept) TimePointd33
# 1            1            0
# 2            1            0
# 3            1            0
# 4            1            1
# 5            1            1
# 6            1            1
# 7            1            0
# 8            1            0
# 9            1            0
# 10           1            1
# 11           1            1
# 12           1            1


# --------------------- Simple model: Gene expr ~ Status --------------------
#################
#  Exercise 1:
#################

# This time create and investigate the model matrix for the variable "Status".
#  1. Create a model formula to investigate the effect of “Status” on gene expression.
simple.model <- as.formula(~ Status)

#  2. Look at the model matrix and identify which is the reference group in your model.
model.matrix(simple.model, data = sampleinfo)
#    (Intercept) StatusUninfected
# 1            1                0
# 2            1                0
# 3            1                0
# 4            1                0
# 5            1                1

## Note Status Infected was set as the reference, change it!
sampleinfo <- mutate(sampleinfo, Status = fct_relevel(Status, "Uninfected"))
model.matrix(simple.model, data = sampleinfo) %>% head
#   (Intercept) StatusInfected
# 1           1              1
# 2           1              1

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  4. Build DESeq2 object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We need:
# * The `txi` object containing the counts
# * The `sampleinfo` data frame containing the sample metadata
# * The `simple.model` design formula

## Create the DESeq2 object
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = simple.model)
class(ddsObj.raw)
# [1] "DESeqDataSet"
# attr(,"package")
# [1] "DESeq2"

ddsObj.raw
# class: DESeqDataSet
# dim: 35896 12
# metadata(1): version
# assays(2): counts avgTxLength
# rownames(35896): ENSMUSG00000000001 ENSMUSG00000000003 ... ENSMUSG00000118657 ENSMUSG00000118658
# rowData names(0):
# colnames(12): SRR7657878 SRR7657881 ... SRR7657873 SRR7657875
# colData names(4): SampleName Replicate Status TimePoint

## Sample metadata stored in:
ddsObj.raw$SampleName
# [1] "SRR7657878" "SRR7657881" "SRR7657880" "SRR7657874" "SRR7657882" "SRR7657872"
# [7] "SRR7657877" "SRR7657876" "SRR7657879" "SRR7657883" "SRR7657873" "SRR7657875"
ddsObj.raw$Replicate
# [1] "1" "2" "3" "1" "2" "3" "1" "2" "3" "1" "2" "3"
ddsObj.raw$Status
# [1] Infected   Infected   Infected   Infected   Uninfected Infected   Uninfected Uninfected
# [9] Uninfected Uninfected Infected   Uninfected
# Levels: Uninfected Infected
ddsObj.raw$TimePoint
# [1] "d11" "d11" "d11" "d33" "d33" "d33" "d11" "d11" "d11" "d33" "d33" "d33"

## Salmon's estimated gene count-like data from: Σ(TPM × tx length × library size) over all gene txs
txi$counts[1:5, 1:5]
#                    SRR7657878 SRR7657881 SRR7657880 SRR7657874 SRR7657882
# ENSMUSG00000000001       1039   1005.888    892.000    917.360   1136.691
# ENSMUSG00000000003          0      0.000      0.000      0.000      0.000
# ENSMUSG00000000028         65     74.000     72.000     44.000     45.999
# ENSMUSG00000000037         39     47.000     29.001     54.001     67.000
# ENSMUSG00000000049          8      9.000      4.000      4.000      4.000

## Count data (from txi$counts but rounded)
counts(ddsObj.raw)[1:5, 1:5]
#                    SRR7657878 SRR7657881 SRR7657880 SRR7657874 SRR7657882
# ENSMUSG00000000001       1039       1006        892        917       1137
# ENSMUSG00000000003          0          0          0          0          0
# ENSMUSG00000000028         65         74         72         44         46
# ENSMUSG00000000037         39         47         29         54         67
# ENSMUSG00000000049          8          9          4          4          4


## Filter out the unexpressed genes (DESeq2 also filter by "independent filtering")
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]
dim(ddsObj.filt)
# [1] 20091    12


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  5. Differential expression analysis with DESeq2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# 1.- Estimate "median ratio" normalisation size factors per sample
# and adjust for average transcript length on a per gene per sample basis.

ddsObj <- estimateSizeFactors(ddsObj.filt)

## One factor per gene per sample
dim(normalizationFactors(ddsObj))
# [1] 20091    12

head(normalizationFactors(ddsObj), 2)
#                    SRR7657878 SRR7657881 SRR7657880 SRR7657874 SRR7657882 SRR7657872 SRR7657877 SRR7657876
# ENSMUSG00000000001  0.9650391  0.9689449  0.9492783  0.9018189   1.209922  0.9615574   1.107827  1.0195031
# ENSMUSG00000000028  1.1565743  1.0085953  0.9422121  0.8798304   1.227264  1.0086842   1.253589  0.9540072
#                    SRR7657879 SRR7657883 SRR7657873 SRR7657875
# ENSMUSG00000000001  0.9599088  0.9028895   1.047268   1.047436
# ENSMUSG00000000028  0.7731541  0.7550835   1.060908   1.127777


## Lets compare logcounts vs lognorm-counts for one sample vs the others
## (but note DESeq2 doesn't directly operate on lognorm/norm counts, it uses raw counts ~ NB!)

logcounts <- log2(counts(ddsObj, normalized = FALSE)  + 1)
limma::plotMA(logcounts, array = 5, ylim =c(-5, 5))
abline(h = 0, col = "red")

logNormalizedCounts <- log2(counts(ddsObj, normalized = TRUE)  + 1)
limma::plotMA(logNormalizedCounts, array = 5, ylim =c(-5, 5))
abline(h = 0, col = "red")


## 2.- Estimate gene-wise dispertion parameters with MLE
ddsObj <- estimateDispersions(ddsObj)
# gene-wise dispersion estimates <- MLE estimates from the data
# mean-dispersion relationship <- fit curve across genes
# final dispersion estimates <- shrink estimates towards trend

plotDispEsts(ddsObj)

## 3.- Negative Binomial GLM fitting
# nbinomWaldTest() fits the GLM per gene (estimate betas) and computes the Wald statistic
ddsObj <- nbinomWaldTest(ddsObj)


## These 3 steps can be performed at once with DESeq()
ddsObj <- DESeq(ddsObj.filt)
# log2 fold change (MLE): Status Infected vs Uninfected
# Wald test p-value: Status Infected vs Uninfected
# DataFrame with 6 rows and 6 columns
#                      baseMean log2FoldChange     lfcSE      stat    pvalue      padj
#                     <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# ENSMUSG00000000001 1102.56094    -0.00802952  0.102877 -0.078050 0.9377883  0.975584
# ENSMUSG00000000028   58.60055     0.30498077  0.254312  1.199239 0.2304350  0.480598
# ENSMUSG00000000037   49.23586    -0.05272685  0.416862 -0.126485 0.8993479  0.961314
# ENSMUSG00000000049    7.98789     0.38165132  0.644869  0.591827 0.5539663  0.772123
# ENSMUSG00000000056 1981.00402    -0.16921845  0.128542 -1.316449 0.1880236  0.426492
# ENSMUSG00000000058  606.92752    -0.17499624  0.106028 -1.650466 0.0988476  0.286870


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  6. Generate a table of differential expression results
#     (adjust pvalues by setting alpha to 0.05)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
results.simple <- results(ddsObj, alpha = 0.05)
results.simple %>% head

#################
#  Exercise 2:
#################
#  Now we have made our results table using our simple model, let have a look at which
#  genes are changing and how many pass our 0.05 threshold.
#
#   a) how many genes are significantly (with an FDR < 0.05) up-regulated?
sum(results.simple$padj < 0.05 & results.simple$log2FoldChange > 0, na.rm = TRUE)
# [1] 1879

#   b) how many genes are significantly (with an FDR < 0.05) down-regulated?
sum(results.simple$padj < 0.05 & results.simple$log2FoldChange < 0, na.rm = TRUE)
# [1] 1005

#   c) Here is the results table for two of the genes:
goi <- c("ENSMUSG00000053747", "ENSMUSG00000048763")
as.data.frame(results.simple[goi, ]) %>%
         mutate(across(-baseMean, ~signif(.x, 3))) %>%
         mutate(across(baseMean, ~round(.x, 3))) %>%
         knitr::kable()
#     |                   | baseMean| log2FoldChange| lfcSE|  stat|  pvalue|   padj|
#     |:------------------|--------:|--------------:|-----:|-----:|-------:|------:|
#     |ENSMUSG00000053747 |   34.930|          -1.61| 0.531| -3.04| 0.00238| 0.0199|
#     |ENSMUSG00000048763 |   30.853|          -2.46| 1.450| -1.70| 0.08970| 0.2700|
#


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 7. Additive model: Gene expr ~ TimePoint + Status
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
additive.model <- as.formula(~ TimePoint + Status)
model.matrix(additive.model, data = sampleinfo) %>% head
#   (Intercept) TimePointd33 StatusInfected
# 1           1            0              1
# 2           1            0              1
# 3           1            0              1
# 4           1            1              1
# 5           1            1              0

#################
#  Exercise 3:
#################
## Build DESeq2 object again
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = additive.model)

## Filter again
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep, ]

## Run DESeq2
ddsObj <- DESeq(ddsObj.filt)

## Extract results
results.additive <- results(ddsObj, alpha = 0.05)

## The default contrast (the one recorded in the last col of model matrix)
results.additive %>% head
# log2 fold change (MLE): Status Infected vs Uninfected
# Wald test p-value: Status Infected vs Uninfected
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat    pvalue      padj
# <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# ENSMUSG00000000001 1102.56094     -0.0110965  0.106195 -0.104492  0.916779  0.967428
# ENSMUSG00000000028   58.60055      0.3007930  0.265626  1.132391  0.257470  0.514578
# ENSMUSG00000000037   49.23586     -0.0481414  0.429685 -0.112039  0.910793  0.965220
# ENSMUSG00000000049    7.98789      0.4110498  0.656171  0.626437  0.531028  0.757304
# ENSMUSG00000000056 1981.00402     -0.1907691  0.119694 -1.593809  0.110979  0.314608
# ENSMUSG00000000058  606.92752     -0.1713459  0.107857 -1.588636  0.112143  0.316383

## Look at all contrasts (coeffs) for which we have generated results
resultsNames(ddsObj)
# [1] "Intercept"                     "TimePoint_d33_vs_d11"          "Status_Infected_vs_Uninfected"

#################
#  Exercise 4:
#################
# Gene expr = B0 + B1(TimePointd33) + B2(StatusInfected), where:
# B0 = Intercept = mean expression of gene in samples at TimePoint d11 and Uninfected
# B1 = diff in mean expr in d33 vs d11 samples, within both Infected and Uninfected groups
# B2 = diff in mean expr in Infected vs Uninfected samples, within both d33 and d11 groups


## Rename results
results.InfectedvUninfected <- results.additive
rm(results.additive)

## Top 100 genes by adjusted pval
topGenesIvU <- as.data.frame(results.InfectedvUninfected) %>%
rownames_to_column("GeneID") %>%
    top_n(100, wt = -padj)
topGenesIvU %>% head


#################
#  Exercise 5:
#################
#  a) Retrieve the results for the contrast of d33 versus d11.
results_d33_v_d11 <- results(ddsObj, name = "TimePoint_d33_vs_d11")

#  b) How many differentially expressed genes are there at FDR < 0.05?
sum(results_d33_v_d11$padj < 0.05, na.rm = TRUE)
# [1] 109


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 8. Interaction model: Gene expr ~ TimePoint + Status + TimePoint:Status
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Variance stabilizing transformation (to mitigate HVG contribution to sample variability)
vstcounts <- vst(ddsObj.raw, blind = TRUE)
## PCA
plotPCA(vstcounts,  intgroup = c("Status", "TimePoint"))

#################
#  Exercise 6:
#################

#  a) Create a new DESeq2 object using a model with an interaction between TimePoint and Status.
interaction.model <- as.formula(~ TimePoint*Status)

# TimePointd33:StatusInfected: indicator for day and infection status (1 if d33 and Infected),
# so we add this term only if d33 and Infected
model.matrix(interaction.model, data = sampleinfo) %>% head
#   (Intercept) TimePointd33 StatusInfected TimePointd33:StatusInfected
# 1           1            0              1                           0
# 2           1            0              1                           0
# 3           1            0              1                           0
# 4           1            1              1                           1
# 5           1            1              0                           0
# 6           1            1              1                           1

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = interaction.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

# b) Run DESeq command and create a new analysis object called ddsObj.interaction
ddsObj.interaction <- DESeq(ddsObj.filt)

# c) Extract results using the default results command.
#    What is the contrast that these results are for?
results.int <- results(ddsObj.interaction)
head(results.int, 2)
# log2 fold change (MLE): TimePointd33.StatusInfected
# Wald test p-value: TimePointd33.StatusInfected
# DataFrame with 2 rows and 6 columns
#                     baseMean log2FoldChange     lfcSE      stat    pvalue      padj
#                    <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# ENSMUSG00000000001 1102.5609       0.305525  0.199432  1.531973  0.125529  0.608828
# ENSMUSG00000000028   58.6005      -0.256926  0.556994 -0.461272  0.644603  0.933365

# B0 = Intercept = mean expr in day11, Unifected samples
# B1 = TimePoint_d33_vs_d11 = diff in mean expr in day33 vs day11 samples, among Unifected samples
# B2 = Status_Infected_vs_Uninfected = diff in mean expr in Infected vs Unifected samples, among d11 samples
# B1 + B3 = TimePoint_d33_vs_d11 + TimePointd33.StatusInfected = diff in mean expr in day33 vs day11 samples, among Infected samples
# B2 + B3 = Status_Infected_vs_Uninfected + TimePointd33.StatusInfected = diff in mean expr in Infected vs Unifected samples, among d33 samples

resultsNames(ddsObj.interaction)
# [1] "Intercept"                     "TimePoint_d33_vs_d11"          "Status_Infected_vs_Uninfected"
# [4] "TimePointd33.StatusInfected"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 9. Extracting specific contrasts from an interaction model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Difference between Infected and Uninfected at 11 days post infection = B2
results.interaction.11 <- results(ddsObj.interaction,
                                  name = "Status_Infected_vs_Uninfected",
                                  alpha = 0.05)

## Difference between Infected and Uninfected at 33 days post infection = B2 + B3
results.interaction.33 <- results(ddsObj.interaction,
                                  contrast = list(c("Status_Infected_vs_Uninfected",
                                                    "TimePointd33.StatusInfected")),
                                  alpha = 0.05)

## Number of genes with padj < 0.05 for Infected v Uninfected at day 11:
sum(results.interaction.11$padj < 0.05, na.rm = TRUE)
# [1] 1072

## Number of genes with padj < 0.05 for Infected v Uninfected at day 33:
sum(results.interaction.33$padj < 0.05, na.rm = TRUE)
# [1] 2782

#################
#  Exercise 7:
#################
# Let's investigate the uninfected mice

#   1. Extract the results for d33 v d11 for Infected mice.
#      How many genes have an adjusted p-value less than 0.05?
results.interaction.Inf <- results(ddsObj.interaction,
                                  contrast = list(c("TimePoint_d33_vs_d11",
                                                    "TimePointd33.StatusInfected")),
                                  alpha = 0.05)
sum(results.interaction.Inf$padj < 0.05, na.rm = TRUE)
# [1] 1134

#   2. Extract the results for d33 v d11 for Uninfected mice.
#      How many genes have an adjusted p-value less than 0.05? Is this remarkable?
results.interaction.Uninf <- results(ddsObj.interaction,
                                   name = "TimePoint_d33_vs_d11",
                                   alpha = 0.05)
sum(results.interaction.Uninf$padj < 0.05, na.rm = TRUE)
# [1] 1 <- not surprising!







