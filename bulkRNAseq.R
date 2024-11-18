## Load the libraries
library(readr)
library(dplyr)
library(tidyverse)
library(BiocManager)
library(DESeq2) 
library(vsn)
library(pheatmap)
library(ComplexHeatmap)
library(ggplot2)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(hexbin)
library(ggrepel)
library(clusterProfiler)


## Load the counts dataset from seq2science results directory.
#read in counts and filter out low counts
data.counts <- read.table("counts/GRCm39-counts.tsv", header=TRUE, row.names=1)

#total number of mapped reads in the pool
data.reads.tot <- sum(data.counts)
data.reads.tot

#check the number of mapped reads per sample based on the gene count table
data.reads.sum <-colSums(data.counts)
data.reads.sum
barplot(data.reads.sum, main = "Number of mapped reads per sample", col = c("sienna3", "turquoise4", "navajowhite"), las=2, cex.names = 0.7)

#check the number of detected genes per sample
data.gene.sum <- colSums(data.counts !=0)
data.gene.sum
barplot(data.gene.sum, main = "Number of detected genes per sample", col = c("sienna3", "turquoise4", "navajowhite"), las=2, cex.names = 0.7)


## Gene annotation
#filter the data, more than 10 reads inmore than 3 columns
data.counts <- data.counts[rowSums(data.counts<10) <3,]

#rename the genes from ENSMUSG00000040952.17 to ENSMUSG00000040952
rownames(data.counts) <- sub("\\..*", "", rownames(data.counts))

#annotate genes in counts, if there is NA gene name, keep the gene symbol
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(data.counts), keytype = "ENSEMBL", column = "SYMBOL")
data.counts$genes <- ifelse(is.na(gene_names), rownames(data.counts), gene_names)

#change the gene column to rownames
data.merge.duplicates <- data.counts %>% group_by(genes) %>% summarise_all(sum)
data.counts.final <- column_to_rownames(data.merge.duplicates, var="genes")

#reorder the columns
data.counts.final <- data.counts.final[,c(3,5,10,1,4,7,2,9,12,6,8,11)]




#Preparing the metadata table for DESeq2
metadata <- data.frame(sample = colnames(data.counts.final))
metadata$group <- c("WT", "WT", "WT", "MyoDicre", "MyoDicre", "MyoDicre", "FusdNLS", "FusdNLS", "FusdNLS", "FusdNLS_MyoDicre", "FusdNLS_MyoDicre", "FusdNLS_MyoDicre") 
metadata$litter <- c("B", "C", "E", "A", "C", "D", "A", "E", "F", "D", "E", "F")

rownames(metadata) <- metadata$sample
metadata$sample <- NULL
print(metadata)


## PCA
library(edgeR)
library(factoextra)

# Calculate the counts per million (CPM) and log-transform the data
pcaData <- cpm(data.counts.final, log = TRUE)

# Perform PCA on the transposed data (samples as rows, features as columns)
pcaDatares <- prcomp(t(pcaData))

# Create a data frame from the PCA results for plotting
pca_df <- as.data.frame(pcaDatares$x)
pca_df$sample <- rownames(pca_df)

# Ensure the row names of the metadata correspond to the PCA samples
metadata$sample <- rownames(metadata)

# Merge PCA results with metadata
pca_df <- merge(pca_df, metadata, by = "sample")

# Create a custom ggplot for PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = litter, shape = group, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel() +  # Add text labels for samples
  labs(title = "PCA of all samples", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +  # Set color palette
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8))




## DESeq


# In the `results()` function, the comparison of interest can be defined using contrasts.  
# To explore the output of `results()` function, I transformed the dds data, ordered results by significance (p-value), scaled the data, and plotted top 20 differently expressed genes in heatmap.  
# Volcano plot shows significantly (y-axis p adj. value) up/downregulated (x-axis log2FoldChange) genes.  
# In GO enrichment analysis, DE genes below padj 0.05 and log2FoldChange over 1.5 were chosen. The last plot summarizes GO terms into clusters based on their semantic similarity of the functions.  


### FusdNLS vs WT
#I tried to use the litter as a batch affect. However, this was not possible due to an error of the design matrix having the same number of samples and coefficients to fit, so estimation of dispersion was not possible. The same happened in all other contrasts but FusdNLS vs FusdNLS_MyoDicre.
#reorder the columns
data.counts.sep <- data.counts.final[,c(1,2,3,7,8,9)]

#Preparing the metadata table for DESeq2
metadata <- data.frame(sample = colnames(data.counts.sep))
metadata$group <- c("WT", "WT", "WT", "FusdNLS", "FusdNLS", "FusdNLS") 
metadata$litter <- factor(c("B", "C", "E", "A", "E", "F"))

rownames(metadata) <- metadata$sample
metadata$sample <- NULL

#construct a DEseq object
dds <- DESeqDataSetFromMatrix(countData = data.counts.sep,
                       colData = metadata,
                       design = ~ group)

#remove rows with low count, keep rows with at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$group <- relevel(dds$group, ref = "WT")

#run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

#explore results
summary(res)
colData(dds)

#MA plot
plotMA(res)

#data transformation and no removing the batch effect
vsd <- vst(dds, blind = TRUE) # blind=TRUE to calculate across-all-samples variability, blind=FALSE to calculate within-group variability 
plotPCA(vsd, "litter") + ggtitle("PCA of littermates")
plotPCA(vsd, intgroup=c("litter", "group")) + ggtitle("PCA of littermates and genotypes")


#effect of data transformation on the variance
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

#get top DEG
genes <- res[order(res$pvalue),] %>%
  head(20) %>%
  rownames()
heatmap <- assay(vsd)[genes,]

#scale counts for visualization
heatmap <- t(scale(t(heatmap)))

#Add annotation
heatmapColAnnot <- data.frame(colData(vsd)[,"group"])
heatmapColAnnot <- HeatmapAnnotation(df = heatmapColAnnot)

# Plot as heatmap
ComplexHeatmap::Heatmap(heatmap,
                        top_annotation = heatmapColAnnot,
                        cluster_rows = TRUE, cluster_columns = FALSE)

# Volcano plot
##add a column to df to specify UP- and DOWN-regulated genes
res.df <- as.data.frame(res)
res.df$diffexpressed <- "NO"
res.df$diffexpressed[res.df$log2FoldChange > 0.6 & res.df$padj < 0.05] <- "UP"
res.df$diffexpressed[res.df$log2FoldChange < -0.6 & res.df$padj < 0.05] <- "DOWN"
write.csv2(res.df, file = "bulk_all_wtxfus.csv")

##volcano plot coloured according to UP- and DOWN-regulated genes
ggplot(res.df, aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed, label=rownames(res.df))) + 
  geom_point() + 
  geom_label_repel(data = subset(res.df), aes(label=rownames(res.df)) , max.overlaps = 20) + 
  geom_vline(xintercept = c(-0.6, 0.6), col="blue") + 
  geom_hline(yintercept = -log10(0.05), col="blue") + 
  scale_color_manual(values=c("maroon3", "navyblue", "seagreen2")) + 
  theme_minimal() #+ geom_label_repel(data = subset(res_as_df, rownames == "Fip1L1"), aes(label=rownames(res_as_df)))



# GENE ONTOLOGY ENRICHMENT

res.df.DEgenes <- subset(res.df, padj < 0.05 & abs(log2FoldChange) > log2(1.5))

DEgenesGO <- enrichGO(gene = rownames(res.df.DEgenes),
                      keyType = "SYMBOL",
                      ont = "BP",
                      OrgDb = org.Mm.eg.db)
DEgenesGOTable <- as.data.frame(DEgenesGO)

GO.genename <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(res.df.DEgenes), columns = c("GENENAME"), keytype = "SYMBOL")

#data wrangling ow multiple entries GeneID followed by collapsing Description
GO.selection <- data.frame(geneID = DEgenesGOTable$geneID, Description=DEgenesGOTable$Description) %>% 
  separate_rows(geneID, sep = "/") %>% 
  group_by(geneID) %>% 
  summarize(Description = paste(Description, collapse = "/"))

#add genename and description to res.df.DEgenes
res.df.DEgenes$Genename <- GO.genename$GENENAME
res.df.DEgenes$Description <- NA
GO.match <- match(rownames(res.df.DEgenes), GO.selection$geneID)
res.df.DEgenes$Description[!is.na(GO.match)] <- GO.selection$Description[GO.match[!is.na(GO.match)]]

write.csv2(res.df.DEgenes, "bulk_deseq_wtxfus.csv")

# plots
barplot(DEgenesGO, showCategory = 10)
dotplot(DEgenesGO, showCategory = 10)

# up/down plot
up <- res.df.DEgenes %>% filter(diffexpressed == "UP")
down <- res.df.DEgenes %>% filter(diffexpressed == "DOWN")
comparelist <- list(rownames(up), rownames(down))
names(comparelist)<-c("upregulated","downregulated") 
cclust<-compareCluster(geneCluster = comparelist, 
               fun = enrichGO,
               OrgDb= org.Mm.eg.db,
               keyType = "SYMBOL",
               ont= "BP",
               universe=data.counts.final)
dotplot(cclust, showCategory = 6, font.size = 10)

# semantics
GO_ID <- DEgenesGOTable$ID[DEgenesGOTable$p.adjust < 0.1]
library(simplifyEnrichment)
simplifyGO(GO_ID)

