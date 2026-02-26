## visualisation of phylogenetic tree using ggtree
library(ggtree)
tree <- read.tree("combined.treefile")
ggtree(tree) +
  geom_tiplab() +
  theme_tree2() +
  xlim(0, max(tree$edge.length) * 1.1)




## visualisation of PCA results using ggplot2
pca <- read.table("/Users/pmonsieurs/programming/ext_FA5_ethiopia_2026/results/plink/full.filtered.QC.pruned.PCA.eigenvec", header=FALSE)
colnames(pca)[1:2] <- c("FID","IID")
library(ggplot2)
ggplot(pca, aes(x = V3, y = V4)) +  # PC1 vs PC2
  geom_point(aes(color = FID)) +
  geom_text(aes(label = IID), vjust = -0.5, size = 3) +  # Add labels above points
  theme_minimal() +
  xlab("PC1") + ylab("PC2")


## visualisation of admixture results using barplot
library(ggplot2)
library(reshape2)

## read ADMIXTURE output with the best K value (lowest CV error), in 
## this case K=6
Q <- read.table("combined.QC.pruned.6.Q")

## add sample IDs based on the .fam file (assuming the order is the same as in the .Q file)
samples <- read.table("combined.QC.pruned.fam")[,1:2]  # FID and IID
Q <- cbind(samples, Q)
colnames(Q) <- c("FID", "IID", paste0("Ancestry", 1:6))

## do melting to long format for ggplot2
Q_long <- melt(Q, id.vars=c("FID", "IID"), 
               variable.name="Ancestry", value.name="Proportion")

## create admixture plot
ggplot(Q_long, aes(x=IID, y=Proportion, fill=Ancestry)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=6)) +
  ylab("Ancestry proportion") +
  xlab("Sample") +
  scale_fill_brewer(palette="Set2")          




####### spare code #########

## run PCA directly in R
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")
library(SNPRelate)

vcf_dir = "/Users/pmonsieurs/programming/ext_FA5_ethiopia_2026/results/gatk/"
vcf.fn  <- paste0(vcf_dir, "full.filtered.vcf.gz")
gds.fn  <- paste0(vcf_dir, "full.filtered.gds")

snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only")
genofile <- snpgdsOpen(gds.fn)

## do ld pruning
set.seed(100)
snpset <- snpgdsLDpruning(genofile,
                          autosome.only = FALSE,
                          ld.threshold=0.2)

snpset.id <- unlist(snpset)


pca <- snpgdsPCA(genofile,
                 snp.id=snpset.id,
                 autosome.only = FALSE,
                 num.thread=2)

pc.percent <- pca$varprop*100

# convert PCA result to a data frame for plotting
pca.df <- data.frame(
  sample.id = pca$sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2]
)

# quick PCA plot
library(ggplot2)
ggplot(pca.df, aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text(aes(label = sample.id), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("PC1") +
  ylab("PC2")


