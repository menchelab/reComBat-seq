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
devtools::install_github("jas-st/reComBat-seq")
```


## Usage and Tutorial

An example use case has been included in the `tutorial` folder. It demonstrates batch correction on the RNA-seq dataset 
used in the ComBat-seq[2] paper, which contains breast cancer cell lines across three batches with GFP controls and oncogene-overexpressing samples (HER2, EGFR, KRAS). 
A raw count matrix from RNA-Seq studies has to be provided, without any normalization or transformation, as well as a vector for batch separation). 
The `group` parameter specifies additional biological covariates, in this case the disease label of the cells. 
For multiple biological variables the `covar_mod` parameter can be used. To demonstrate the main feature of reComBat-seq a singular matrix will be used.

```r
library(reComBatseq)
recombatseq_df <- reComBat_seq(cts_sub, batch=batch_sub, group=group_sub, 
                               covar_mod = covmat,
                               lambda_reg = 0.8, alpha_reg = 0.3, 
                               num_threads = 0)
```

| Raw Data      | Corrected Data (Singular Design) |
| ------------- | -------------  |
| <img src="https://github.com/jas-st/reComBat-seq/blob/main/tutorial/PCA_raw.png" width="500">  | <img src="https://github.com/jas-st/reComBat-seq/blob/main/tutorial/PCA_recombatseq.png" width="500">   |
  

## Arguments

  - `counts` - raw count matrix from genomic studies (dimensions gene x sample)
  - `batch` - batch covariate (only one vector allowed)
  - `group` - vector / factor for condition of interest
  - `covar_mod` - model matrix for other covariates to include in linear model besides batch and condition of interest
  - `num_threads` - number of threads for parallel gene-wise regression using OpenMP, default is single-thread

The regularization can be adjusted via the following parameters:

  - `lambda_reg` - controls the strength of the regularization, $\lambda$ in the above equation
  - `alpha_reg` - controls the elastic net tuning, $\alpha$ in the above equation

## References
1. Chen Y, Chen L, Lun ATL, Baldoni P, Smyth GK (2025). “edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.” Nucleic Acids Research, 53(2), gkaf018. doi:10.1093/nar/gkaf018.
2. Yuqing Zhang, Giovanni Parmigiani, W Evan Johnson, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3, 1 September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078
