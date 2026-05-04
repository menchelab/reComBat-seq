## reComBat-seq

reComBat-seq is a batch effect adjustment tool designed for bulk RNA-seq count data with underdetermined experimental designs. Building upon the negative binomial regression framework used in ComBat-seq[2], reComBat-seq incorporates Elastic Net regularization to address convergence issues when handling confounded data. It takes raw **untransformed, raw** count matrices as input and requires a known batch variable.

By applying regularized negative binomial regression, reComBat-seq models batch effects while enabling stable parameter estimation with rank-deficient design matrices. The adjusted data produced preserves the integer nature of counts, ensuring compatibility with widely used differential expression tools such as edgeR[1] and DESeq2. 

This formulation extends the standard negative binomial regression model by incorporating regularization penalties, with `NLL`representing the negative log-likelihood function.
```math
L(\beta) = NLL(\beta) + \lambda \bigg( \alpha \| \beta \|_1 + \bigg( \frac{1-\alpha}{2} \bigg) \| \beta \|_2^2 \bigg)
```

Parallelisation is implemented via OpenMP. OpenMP is not available by default on macOS and must be installed separately (e.g. via Homebrew: `brew install libomp`). 
If OpenMP is unavailable, it will fall back to single-threaded execution.

## Installation

Before use install edgeR from Bioconductor. 

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
```

reComBat-seq is then available for installation via GitHub.

```r
# install.packages("devtools")
devtools::install_github("menchelab/reComBat-seq")
```


## Usage

The `tutorial` folder contains some examples including the code used to generate the below plots using the breast cancer data from the ComBat-seq[2] paper as well as code used to correct the psoriasis data used in our study. Using reComBat-seq begins with a raw count matrix. We then need to identify the batch variable we want to remove as well as all the covariates we want to make sure to keep their variation aka `wanted.variation`. Examples of the latter one may be biological covariates like disease status or cell identity. The code below shows the usage principle

```r
library(reComBatseq)
# reading in data from files
raw.counts <- read.table(
    'psoriasis.exp.tsv',
    sep = '\t',
    quote = '',
    header = TRUE,
    row.names = 'X'
)

meta <- read.table(
    'psoriasis.meta.tsv',
    sep = '\t',
    quote = '',
    header = TRUE,
    row.names = 'X'
)

# applying recombatseq correction using disease as wanted covariate
corrected.counts <- reComBat.seq(
    t(raw.counts),
    batch=meta$sra_study_acc, 
    wanted.variation=meta['Disease']
)
```
<img src="https://github.com/menchelab/reComBat-seq/blob/main/tutorial/recombat_comparison-01.png">

Please note that the package currently relies on openMP. This can cause some issues with local openBLAS installations where nested threading occurs in instances where datasets become very large (in our case like 25k cells / samples). This can unfortunately only be fixed by recompiling the local openBLAS installation. As a low effort workaround we recommend the following. Either set `OPENBLAS_NUM_THREADS=1` before starting R like so
```
export OPENBLAS_NUM_THREADS=1
```
or use `Sys` to set it in R before loading the package
```
Sys.setenv(OPENBLAS_NUM_THREADS = '1')
library(reComBat-seq)
```
While only the direct export was tested both should work equally well.

## reComBat-seq on ComBat-seq data
| Raw Data      | Corrected Data |
| ------------- | -------------  |
| <img src="https://github.com/menchelab/reComBat-seq/blob/main/tutorial/PCA_raw.png" width="500">  | <img src="https://github.com/menchelab/reComBat-seq/blob/main/tutorial/PCA_recombatseq.png" width="500">   |
  

## Arguments

  - `counts` - raw count matrix from genomic studies (dimensions gene x sample)
  - `batch` - vector containing batch assignment of samples
  - `wanted.variation` - a data.frame containing the covariates whose variation you want to preserve
  - `num.threads` - number of threads for parallel gene-wise regression using OpenMP, default is single-thread

The regularization can be adjusted via the following parameters:

  - `lambda.reg` - controls the strength of the regularization, $\lambda$ in the above equation
  - `alpha.reg` - controls the elastic net tuning, $\alpha$ in the above equation

## References
1. Chen Y, Chen L, Lun ATL, Baldoni P, Smyth GK (2025). “edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.” Nucleic Acids Research, 53(2), gkaf018. doi:10.1093/nar/gkaf018.
2. Yuqing Zhang, Giovanni Parmigiani, W Evan Johnson, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3, 1 September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078
