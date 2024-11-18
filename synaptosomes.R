library("dplyr")
library("tidyr")
library("ggplot2")
library("DESeq2")
library("AnnotationDbi")
library("org.Mm.eg.db")
library(edgeR)
library(ggrepel)
library(vsn)
library(pheatmap)
library(ggpubr)
library(cowplot)



## Load the count tables to check number of reads and detected genes

#load data
data.counts <- read.table(file="/counts/GRCm39-counts.tsv", sep = '\t', header = T, row.names = 1)
data.counts.ord <- data.counts[, order(colnames(data.counts))]

#total number of mapped reads in the pool
data.reads.tot <- sum(data.counts.ord)
data.reads.tot

#check the number of mapped reads per sample based on the gene count table
data.reads.sum <-colSums(data.counts.ord)
data.reads.sum
barplot(data.reads.sum, main = "Number of mapped reads per sample", col = c("sienna3", "turquoise4", "navajowhite"), las=2, cex.names = 0.7)

#check the number of detected genes per sample
data.gene.sum <- colSums(data.counts.ord !=0)
data.gene.sum
barplot(data.gene.sum, main = "Number of detected genes per sample", col = c("sienna3", "turquoise4", "navajowhite"), las=2, cex.names = 0.7)



## Check the gene counts distribution per sample before filtering

#filter out rows containing only 0
data.wo.0 <- data.counts[rowSums(data.counts != 0) > 0,]
#log transform the data
data.log <- log2(data.counts)
#gather data in two columns
data.counts.gather <- gather(data.counts, Sample, Value)
data.log.gather <- gather(data.log, Sample, Value)

#plot the gather data
#data.counts
p1 <- ggplot(data.counts.gather, aes(x=Sample, y=Value)) + geom_violin() + ggtitle("Violin plot of data.counts")
p2 <- ggplot(data.counts.gather, aes(x=Sample, y=Value)) + geom_boxplot(alpha=0.2) + ggtitle("Boxplot of data.counts")
#data.log
p3 <- ggplot(data.log.gather, aes(x=Sample, y=Value)) + geom_violin() + ggtitle("Violin plot of log data.counts")
p4 <- ggplot(data.log.gather, aes(x=Sample, y=Value)) + geom_boxplot(alpha=0.2) + ggtitle("Boxplot of log data.counts")
combined_plot <- ggarrange(p1, p2, p3, p4, ncol=2, nrow=2)
print(combined_plot)

#loop for histograms
sample.names <- colnames(data.counts)
hist.list <- list()
for (i in seq_along(sample.names)) {
  print(paste("extracting data for", sample.names[[i]]))
  data.loop <- data.log[,i]
  hist.list[[i]] <- ggplot(data.frame(x=data.loop), aes(x)) + geom_histogram(binwidth = 0.5, fill="steelblue", color="black") + labs(title = sample.names[i])
}
combined_plot2 <- plot_grid(plotlist = hist.list, ncol=3, nrow=4)
print(combined_plot2)



## Study the data
## Gene annotation and filtering the gene counts

#filter the data, more than 2 reads in more than 2 columns
data.counts.filt <- data.counts[rowSums(data.counts<2) <8,]

#rename the genes from ENSMUSG00000040952.17 to ENSMUSG00000040952
rownames(data.counts.filt) <- sub("\\..*", "", rownames(data.counts.filt))

#annotate genes in counts, if there is NA gene name, keep the gene symbol
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(data.counts.filt), keytype = "ENSEMBL", column = "SYMBOL")
data.counts.filt$genes <- ifelse(is.na(gene_names), rownames(data.counts.filt), gene_names)

#change the gene column to rownames
data.counts.reordered <- data.counts.filt[!is.na(data.counts.filt$genes),]
rownames(data.counts.reordered) <- NULL
rownames(data.counts.reordered) <- make.names(data.counts.reordered$genes, unique = TRUE)
data.counts.reordered$genes <- NULL



## Heatmap of the raw counts data 

data.mat <- as.matrix(data.counts.reordered)
heatmap(data.mat, Rowv = FALSE, Colv = FALSE, scale = "row")


## Data normalization with CPM

#Before running DESeq2 package I will normalize the samples with Counts per million and have a look how they cluster.
data.norm <- cpm(data.mat)
my_colors <- colorRampPalette(c("papayawhip", "turquoise4"))(100)
heatmap(data.norm, Rowv = FALSE, Colv = FALSE, scale = "row", col = my_colors)



## Gene ontology 

#There are 6101 genes left after the filtering. Let's look at the gene ontology barplot (Biological process and Molecular function) of those expressed genes from the CPM normalized dataset. The heatplot shows function of top 50 genes.

library(clusterProfiler)

#sort the data by mean expression and filter 20 highly expressed genes
data.norm.sort <- data.norm[order(-rowMeans(data.norm)),]
top50 <- data.norm.sort[1:50,]

top50_id <- rownames(top50)

gene_id <- rownames(data.norm.sort)
top50_id

#Biological process ontology
go_BP <- groupGO(gene = gene_id, 
                      OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",  # Biological Process ontology
                      readable = TRUE
                     )
print("Biological process ontology")
barplot(go_BP, main = "Biological process ontology")

go_BP_filt <- groupGO(gene = top50_id, 
                      OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",  # Biological Process ontology
                      readable = TRUE
                     )
heatplot(go_BP_filt)

#Molecular function ontology
go_MF <- groupGO(gene = gene_id, 
                      OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL",
                      ont = "MF",  # Molecular function ontology
                      readable = TRUE
                     )
print("Molecular function ontology")
barplot(go_MF, main = "Molecular function ontology")

go_MF_filt <- groupGO(gene = top50_id, 
                      OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL",
                      ont = "MF",  # Molecular function ontology
                      readable = TRUE
                     )
heatplot(go_MF_filt)


# DESeq

## Forebrain

### Prepare the data

#Running DESeq to compare differential expression between forebrain sorted synaptosomes and non-sorted crude synapses line. I use raw data before CPM normalization.

#create metadata file
metadata <- data.frame(sample=c("F1R1", "F2R1", "F2R2", "FN1", "FN2"))
metadata$condition <- c("Synapses", "Synapses", "Synapses", "Nonsorted", "Nonsorted")
rownames(metadata) <- metadata$sample
metadata$sample <- NULL 
print(metadata)

#make sure rownames in metadata and column names in counts are in correct order
data.counts.DE.F <- data.counts.reordered %>%
  select("F1R1", "F2R1", "F2R2", "FN1", "FN2")
all(colnames(data.counts.DE.F) == rownames(metadata))


#make sure rownames in metadata match column names in counts and are in correct order
all(colnames(data.counts.DE.F) %in% rownames(metadata))
all(colnames(data.counts.DE.F) %in% rownames(metadata))


#run deseq
dds <- DESeqDataSetFromMatrix(countData = round(data.counts.DE.F),
                              colData = metadata,
                              design = ~ condition)
dds

#remove rows with low count, keep rows with at least 5 reads, it is still the same beacuse I did already such filtering
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds

#set the factor level
dds$group <- relevel(dds$condition, ref = "Nonsorted")


### Run the actual deseq package
#run DESeq
dds <- DESeq(dds)
dds

res <- results(dds)
res

summary(res)

res0.5 <- results(dds, alpha = 0.5)
summary(res0.5)

colData(dds)

print("Size factors")
sizeFactors(dds)

resultsNames(dds)


### Plots

#MA plot
DESeq2::plotMA(res)

#data transformation
#commented, not working here
#vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
#head(assay(vsd), 3)

#effect of data transformation on the variance
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

###MA plot
res_as_df <- as.data.frame(res) #[,c("baseMean","log2FoldChange","padj")])
ggplot(res_as_df, aes(x=log10(baseMean), y= log2FoldChange, color=padj)) + geom_point()


###volcano plot
ggplot(res_as_df, aes(x=log2FoldChange, y=-log10(padj), color=padj, label=rownames(res_as_df))) +
  geom_point() +
  geom_text_repel(box.padding = 0.5, max.overlaps = 15) +
  geom_vline(xintercept = c(-0.6, 0.6), col="blue") +
  geom_hline(yintercept = -log10(0.05), col="blue")

##add a column to df to specify UP- and DOWN-regulated genes
res_as_df$diffexpressed <- "NO"
res_as_df$diffexpressed[res_as_df$log2FoldChange > 0.6 & res_as_df$pvalue < 0.05] <- "UP"
res_as_df$diffexpressed[res_as_df$log2FoldChange < -0.6 & res_as_df$pvalue < 0.05] <- "DOWN"

##volcano plot coloured according to UP- and DOWN-regulated genes
ggplot(res_as_df, aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed, label=rownames(res_as_df))) +
  geom_point() +
  geom_label_repel(data = subset(res_as_df), aes(label=rownames(res_as_df)) , max.overlaps = 50, force = 0.5) +
  geom_vline(xintercept = c(-0.6, 0.6), col="blue") +
  geom_hline(yintercept = -log10(0.05), col="blue") +
  scale_color_manual(values=c("maroon3", "seagreen2", "navyblue")) +
  theme_minimal()# + geom_label_repel(data = subset(res_as_df, rownames == "Fip1L1"), aes(label=rownames(res_as_df)))

#heatmap of highly expressed genes
print("heatmap of highly expressed genes on normalized DESeq data")
select_heatmap <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:25]
dds_as_df <- as.data.frame(colData(dds)) #in original is: df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select_heatmap,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=dds_as_df)

#heatmap of differently expressed genes
print("heatmap of differently expressed genes on DESeq data")

signif.genes <- res$padj < 0.05
num.sig.genes<-sum(signif.genes, na.rm = TRUE)
cat("Number of significantly differently expressed genes is:", num.sig.genes, "\n")

select_heatmap2 <- order(res$log2FoldChange[signif.genes], decreasing = TRUE)[1:25]
pheatmap(assay(ntd)[select_heatmap2,],cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE, annotation_col = dds_as_df)


#heatmap of significantly differently expressed genes
print("heatmap of significantly differently expressed genes on DESeq data")


### Gene Ontology of top DESeq upregulated genes

# Filter genes by adjusted p-value (q-value) threshold and log2 fold change (upregulated)
upreg.genes <- complete.cases(res) & res$padj < 0.05 & res$log2FoldChange > 0

# Subset the results to get the upregulated genes
upreg.res <- res[upreg.genes,]

# Extract gene names of upregulated genes
upreg.gene.names <- rownames(upreg.res)

gene.descr <- AnnotationDbi::select(org.Mm.eg.db, keys=upreg.gene.names, keytype = "SYMBOL", columns = "GENENAME")
print(gene.descr)


##Make a csv file of upregulated genes.

#pvalue and foldchange
signif.threshold <- 0.05
foldchange.threshold <- 0.5  # Adjust if needed

# Filter genes based on significance and fold change thresholds
upreg_genes <- complete.cases(res) & res$log2FoldChange > log2(foldchange.threshold) & res$padj < signif.threshold

# Extract the gene names for upregulated genes
upreg_gene_names <- rownames(res[upreg_genes, ])

# Save the upregulated gene names to a CSV file
write.csv(data.frame(GeneNames = upreg_gene_names), file = "synaptosomes_forebrain_upregulated_gene_names.csv", row.names = FALSE)



# Define your p-value threshold
pvalue.threshold <- 0.05  # Adjust this value as needed

# Filter genes based on the p-value threshold
sig_genes <- complete.cases(res) & res$pvalue < pvalue.threshold & res$log2FoldChange > 0

# Extract the gene names for significantly differentially expressed genes
sig_gene_names <- rownames(res[sig_genes, ])

# Save the gene names to a CSV file
write.csv(data.frame(GeneNames = sig_gene_names), file = "synaptosomes_forebrain_sig_gene_names.csv", row.names = FALSE)




