library("DESeq2")
library(tidyverse)
library(data.table)
library("pheatmap")
library("biomaRt") 
library("RColorBrewer")

args = commandArgs(trailingOnly=TRUE)

countfile=args[1]
#countfile="/home/jpeter/DATA/RNAseq/RNAseq072023/2023-08-29T093905_Analysis/seeds/seeds.gene_models.fCounts.txt"
countdata <- fread(countfile)

sample_table_f = args[2]
#sample_table_f <- "/home/jpeter/DATA/RNAseq/RNAseq072023/sampletable.tsv"

design_type=args[3]
#design_type="rep"
outdir=file.path(dirname(countfile), paste0("DESeq2_", design_type))
dir.create(outdir, showWarnings = F)
#setwd(outdir)


countdata <- subset(countdata, select = -c(Chr, Start, End, Strand, Length)) # remove unwanted columns
colnames(countdata) <- gsub("\\.sorted.bam$", "", basename(colnames(countdata))) # rename columns: remove rf.sortedbam
colnames(countdata)

#generate a matrix
countmatrix <- as.matrix(countdata, rownames = TRUE)

sampleTable<- fread(sample_table_f) %>%
  column_to_rownames('fileName')

sampleTable$batch <- factor(sampleTable$batch, levels=c("1","2","3"))

genotypes <- as.data.table(t(combn(unique(sampleTable$condition), 2)))

all(colnames(countmatrix) == row.names(sampleTable))

row.names(sampleTable) = sampleTable$sampleName
colnames(countmatrix) = sampleTable$sampleName

if (design_type == "simple") {
  ddsMatrix <- DESeqDataSetFromMatrix(countData = countmatrix, colData = sampleTable, design = ~ condition )
} else if (design_type == "rep") {
  ddsMatrix <- DESeqDataSetFromMatrix(countData = countmatrix, colData = sampleTable, design = ~ condition + batch)
  
}

dds <- ddsMatrix[ rowSums(counts(ddsMatrix)) > 1, ]
dds <- DESeq(dds, betaPrior=TRUE) 
write.table(counts(dds, normalized=TRUE), file= file.path(outdir,"normalized-counts-Treated.txt"), sep='\t', col.names=NA, quote = FALSE)

addBiomaRtAnnotation <- function(normObj, biomart="plants_mart",
                                 dataset="athaliana_eg_gene", host="https://plants.ensembl.org",
                                 features=c("ensembl_peptide_id","ensembl_transcript_id","ensembl_gene_id",
                                            "external_gene_name"))#, GO=FALSE)
{
  stopifnot(!missing(normObj))
  stopifnot()
  filter <- match.arg(features)
  mart <- useMart(biomart,dataset,host=host)
  atr <- c(features,"description")
  geneInfo <- getBM(rownames(normObj),mart = mart, attributes = atr,filters = filter) 
  if (length(geneInfo$description)==0) {stop(sprintf("rowIDs from count table are different than IDs from \"%s\" from ensembl\nChoose another feature ID, or another dataset from biomaRt, or check that your IDs are in the format as the one from ensembl", filter))}
  print(length(geneInfo$description))
  print(length(row.names(normObj)))
  rownames(geneInfo) <- geneInfo[,1] # C'est toujours la premiere colonne qui contient les identifiants.
  normObj <- data.frame(normObj)
  normObj$id  <- 1:nrow(normObj)
  normObj <- merge(normObj, geneInfo, by="row.names", all.x=TRUE)
  rownames(normObj) <- normObj$Row.names
  normObj$Row.names <- NULL
  normObj <-normObj[order(normObj$id), ]
  normObj$id <- NULL
  normObj
}

compare_pair <- function(x) {
  #print(x)
  print(x[[1]])
  print(x[[2]])
  pairdir <- file.path(outdir, paste0(x[[1]],"_vs_", x[[2]]))
  dir.create(pairdir, showWarnings = F)
  res <- results(dds, contrast=c("condition", x[[1]], x[[2]]))
  pdf(file.path(pairdir, paste0("MAplot", x[[1]], "_vs_", x[[2]],".pdf")), width = 6, height = 5)
  MAplot <- plotMA(res, main = paste(x[[1]], "vs", x[[2]]), ylim=c(-3,3),
                   alpha=0.05, xlab = "mean of normalized counts")
  dev.off()
  DEG_up <- as.data.frame(subset(res, padj < 0.05&log2FoldChange>0))
  DEG_down <- as.data.frame(subset(res, padj < 0.05&log2FoldChange<0))
  if (nrow(DEG_up>0)) {
    DEG_up_a <- addBiomaRtAnnotation(DEG_up, features="ensembl_gene_id")
    write_tsv(DEG_up_a, file=file.path(pairdir, paste0("Up_genes_",x[[1]],"_vs_",x[[2]],".tsv")))
  }
  if (nrow(DEG_down>0)) {
    DEG_down_a <- addBiomaRtAnnotation(DEG_down, features="ensembl_gene_id")
    write_tsv(DEG_down_a, file=file.path(pairdir, paste0("Down_genes_",x[[1]],"_vs_",x[[2]],".tsv")))
  
  }
  

}

apply(genotypes,1 ,compare_pair)

rld <- rlogTransformation(dds, blind=FALSE)# size factors are used to normalize the counts before rlog-transforming



#PCA function implemented in DEseq2 package
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
#Extract value of variation %
percentVar <- round(100 * attr(pcaData, "percentVar"))


ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  ggtitle("PCA analysis") +
  theme(plot.title = element_text(hjust = 0.5, size=22))

ggsave(file.path(outdir, "PCA.pdf"))


#PCA function implemented in DEseq2 package
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
#Extract value of variation %
percentVar <- round(100 * attr(pcaData, "percentVar"))
#draw customozied plot using ggplot function. Point color and shape are set according to the conditions
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  ggtitle("PCA analysis") +
  theme(plot.title = element_text(hjust = 0.5, size=22)) +
  geom_text(aes(label = name),  vjust = "inward", hjust = "inward", position=position_nudge(x = 1, y= 0.5))

ggsave(file.path(outdir, "PCA2.pdf"))

# Calculate distances between sample using rlog-transformed data
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
# Add names to the created matrix
rownames(sampleDistMatrix) <- rld$sampleName
colnames(sampleDistMatrix) <- rld$sampleName
# Choose a color palette using the RColorBrewer package
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

while (!is.null(dev.list()))  dev.off()


pdf(file.path(outdir, "Sampletosample_distancematrix.pdf"))
pheatmap(sampleDistMatrix,
                clustering_distance_rows=sampleDists,
                clustering_distance_cols=sampleDists,
                col=colors)
dev.off()
# 
# 
# list_DE <- unique(c(which(res_1$padj < 0.05), which(res_2$padj < 0.05), which(res_3$padj < 0.05)))
# 
# mat <- assay(rld)[list_DE, ]
# mat <- mat - rowMeans(mat)
# 
# df <- as.data.frame(colData(rld)[,"condition"], row.names =colnames(mat))
# colnames(df)<- "condition"
# out <- pheatmap(mat, annotation_col=df, main="Heatmap based on log transformation", fontsize_row= 2)
# 
# ##we first create a table with transposon or gene information
# mat2 <- as.data.frame(mat)
# mat2$annotation <- "gene" ## add column with annotation indication, by defaut annotation is "gene"
# mat2[grepl("TE", row.names(mat2)),]$annotation <- "transposon" ## when "TE" is in row.names, "gene"is replaced by "transposon" in the annotation column
# mat2 <- mat2[,10, drop=F]#keep only rownames and annotation column
# 
# ## Choose colors to show transposon and genes
# annotation <- c("darkred", "darksalmon")
# names(annotation) <- c("gene", "transposon")
# 
# ## Choose colors to show WT and mutants
# condition <- c("gray50", "navyblue", "cadetblue2")
# names(condition) <- c("wt", "nrpe1", "nrpd1")
# anno_colors <- list(annotation = annotation, condition = condition)
# 
# ##The heatmap is now drawn by adding the "annotation_row" argument.
# out <- pheatmap(mat, annotation_col=df, annotation_row=mat2, main="Heatmap based on log transformation", fontsize_row= 2, annotation_colors = anno_colors)
# 
# 
# #make table with heatmap results
# out_row <- out$tree_row
# heatmap_table <- mat[out_row$order,]
# #annotate heatmap table
# heatmap_table_a <- addBiomaRtAnnotation(heatmap_table, features="ensembl_gene_id")
# #export heatmap table
# write_tsv(heatmap_table_a, file = file.path(outdir,"heatmap_table_rld.tsv"))
# 
# pdf(file.path(outdir,"Heatmap.pdf"), width=6, height=8)
# out <- pheatmap(mat, annotation_col=df, annotation_row=mat2, main="Heatmap based on log transformation", fontsize_row= 2, annotation_colors = anno_colors)
# dev.off()
# 
# out.clust <- cbind(mat, cluster = cutree(out$tree_row, k = 7))
# #annotate heatmap table
# out.clust_a <- addBiomaRtAnnotation(out.clust, features="ensembl_gene_id")
# #order table according to cluster
# out.clust_a  <- out.clust_a[order(out.clust_a$cluster),]
# #export table with cluster numbers
# write_tsv(out.clust_a, file=file.path(outdir,"cluster.tsv"))
# 
# 
# 
