# Adapted from ComBat-seq (Zhang et al., 2020)
# Original source: https://github.com/zhangyuqing/ComBat-seq
library(reComBatseq)
library(sva)

# read in the dataset
write.csv(cts_sub, "test_matrix.csv")
write.csv(as.data.frame(batch_sub, as.character(group_sub)), "test_metadata.csv")

cts_sub <- as.matrix(read.csv("tutorial/test_matrix.csv", row.names = "X"))
metadata <- read.csv("tutorial/test_metadata.csv")
group_sub <- as.factor(metadata$X)
batch_sub <- metadata$batch_sub

# batch correction - confounded design
covmatdf <- as.data.frame(cbind(group_sub, group_sub))
covmatdf[] <- lapply(covmatdf, as.factor)
colnames(covmatdf) <- c("test1", "test2")
covmat <- model.matrix(~., data = covmatdf)[, -1]

## Normalize library size - divide each column by its sum
cts_norm <- apply(cts_sub, 2, function(x){x/sum(x)})
cts_adj_norm <- apply(combatseq_df, 2, function(x){x/sum(x)})
cts_og_norm <- apply(combatseq_og_df, 2, function(x){x/sum(x)})

library(DESeq2)
library(scales)
library(ggplot2)
library(ggpubr)

## PCA
col_data <- data.frame(Batch=factor(batch_sub), Group=group_sub)
rownames(col_data) <- colnames(cts_sub)

seobj <- SummarizedExperiment(assays=cts_norm, colData=col_data)
pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("Batch", "Group"))
plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) +
  geom_point() +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Unadjusted")
plt
ggsave("tutorial/PCA_raw.png", height=5, width=6, dpi=300)

seobj_adj <- SummarizedExperiment(assays=cts_adj_norm, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("Batch", "Group"))
plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) +
  geom_point() +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="reComBat-Seq (singular design)")
plt_adj
ggsave("tutorial/PCA_recombatseq.png", height=5, width=6, dpi=300)

seobj_adjori <- SummarizedExperiment(assays=cts_og_norm, colData=col_data)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup=c("Batch", "Group"))
plt_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) +
  geom_point() +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Original ComBat-Seq")

ggarrange(plt, plt_adjori, plt_adj, ncol=1, nrow=3, common.legend=TRUE, legend="right")