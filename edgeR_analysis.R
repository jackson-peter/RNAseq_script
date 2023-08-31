library(tidyverse)
library(data.table)
library(edgeR)
library("biomaRt")
#library(xlsx)

######## ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
input_dir=args[1]
tissue=args[2]
suffix=args[3]
design_type=args[4]

#input_dir <- "~/DATA/RNAseq/RNAseq072023/2023-08-29T093905_Analysis/seeds/"
#tissue <- "seeds"
#suffix <- ".gene_models.fCounts.txt"

ft_cnt <- as_tibble(fread(file.path(input_dir, paste0(tissue, suffix)))) %>%
  rename_with(basename) %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(-c(Chr, Start, End, Strand, Length))

colnames(ft_cnt) <- gsub("\\.sorted.bam$","", colnames(ft_cnt)) #

if (tissue=="seeds") {
  names(ft_cnt) <- c("DQGD106" = "Col0_rep1","DQGD107" = "urt1_rep1","DQGD108" = "heso1_rep1","DQGD109" = "uh_rep1",
                     "DQGD1010" = "Col0_rep2","DQGD111" = "urt1_rep2","DQGD112" = "heso1_rep2","DQGD113" = "uh_rep2",
                     "DQGD114" = "Col0_rep3","DQGD115" = "urt1_rep3","DQGD116" = "heso1_rep3","DQGD117" = "uh_rep3")
  
  condition <- c(rep(c("G_Col0", "G_urt1", "G_heso1", "G_urt1heso1"),3))
  
} else if (tissue=="flowers") {
  names(ft_cnt) <- c("DQGD118" = "Col0_rep1","DQGD119" = "urt1_rep1","DQGD120" = "heso1_rep1","DQGD121" = "uh_rep1",
                     "DQGD1022" = "Col0_rep2","DQGD123" = "urt1_rep2","DQGD124" = "heso1_rep2","DQGD125" = "uh_rep2",
                     "DQGD126" = "Col0_rep3","DQGD127" = "urt1_rep3","DQGD128" = "heso1_rep3","DQGD129" = "uh_rep3")
  condition <- c(rep(c("F_Col0", "F_urt1", "F_heso1", "F_urt1heso1"),3))
} else {
  stop("tissue should be either 'seeds' or 'flowers'.")
}



outdir <- file.path(input_dir, paste0("edgeR_", design_type))
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

replicate <- as.factor(c(rep(1,4), rep(2,4), rep(3,4)))
cds <- DGEList(ft_cnt, group = condition)

names( cds )
head(cds$counts) # original count matrix
cds$samples # contains a summary of your samples

#pseudoCounts <- log2(cds$counts+1)


#selr <- rowSums(cpm(x$counts)>5)>=15 to keep cpm >5 in fifteen or more samples

keep <- rowSums(cpm(cds)>1) >= 2
cds <- cds[keep, , keep.lib.sizes=FALSE]
#calculate normalisation factor
cds <- calcNormFactors(cds)
cds$samples #give Normalisation factor


if (design_type=="simple") {
  design <- model.matrix(~0+group,data=cds$samples)
  print(design)
} else if (design_type=="rep") {
  design <- model.matrix(~0+group+replicate,data=cds$samples)
  print(design)
} else {
  stop("design should be either 'simple' or 'rep'.")
}
write.table(design, file.path(outdir, paste0(tissue, "_design.tsv")), sep="\t")

#method="TMM" is the weighted trimmed mean of M-values (to the reference) proposed by #Robinson and Oshlack (2010), where the weights are from the delta method on Binomial #data. If refColumn is unspecified, the library whose upper quartile is closest to the #mean upper quartile is used.

print("test1")
ScaleFactors <- cds$samples$lib.size * cds$samples$norm.factors # effective library size
Exp <- round(t(t(cds$counts)/ScaleFactors) * mean(ScaleFactors)) # standard to round b/c can't have a fraction of a read
list_universe <- as.vector(row.names(Exp)) # list that will be used for GO analysis latter
head(cpm(cds))
print("test2")
write.table(Exp, paste0(tissue, "_Norm_data_total.txt"), sep="\t", col.names = NA)
write.table(cpm(cds), paste0(tissue, "_CPM.txt"), sep="\t", col.names = NA)
#cpm(cds)
print("test3")

myPalette <- c(rep(c("#000000", "#E69F00", "#56B4E9", "#009E73"),3)) ##color blind colors
top = nrow(cds[[1]])
# Output plot as a pdf
pdf( paste0(tissue,"MDS_plot_all.pdf") , width = 7 , height = 7 ) # in inches
plotMDS( cds , top = top, main = "MDS Plot for all conditions", labels = colnames( cds$counts ),col=myPalette)
dev.off() # this tells [R] to close and stop writing to the pdf.
pdf( paste0(tissue,"_MDS_plot_defaut.pdf") , width = 7 , height = 7 ) # in inches
plotMDS( cds , main = "MDS Plot for all conditions", labels = colnames( cds$counts ),col=myPalette)
dev.off() # this tells [R] to close and stop writing to the pdf.
print("test4")
##alternative method
pdf( paste0(tissue,"_MDS_plot_bcv.pdf") , width = 7 , height = 7 ) # in inches
plotMDS( cds , main = "MDS Plot for all conditions", method = "bcv", labels = colnames( cds$counts ),col=myPalette)
dev.off()
print("test5")

if (tissue == "seeds") {
  my.contrasts <-makeContrasts(Gcol0_Gheso1=groupG_Col0-groupG_heso1,
                               Gcol0_Gurt1=groupG_Col0-groupG_urt1,
                               Gcol0_Gurt1heso1=groupG_Col0-groupG_urt1heso1,
                               Gheso1_Gurt1=groupG_heso1-groupG_urt1,
                               Gheso1_Gurt1heso1=groupG_heso1-groupG_urt1heso1,
                               Gurt1_Gurt1heso1=groupG_urt1-groupG_urt1heso1,
                               levels=design)
  
} else if (tissue=="flowers") {
  my.contrasts <-makeContrasts(Fcol0_Fheso1= groupF_Col0-groupF_heso1,
                               Fcol0_Furt1= groupF_Col0-groupF_urt1, 
                               Fcol0_Furt1heso1= groupF_Col0-groupF_urt1heso1,
                               Fheso1_Furt1=groupF_heso1-groupF_urt1,
                               Fheso1_Furt1heso1=groupF_heso1-groupF_urt1heso1,
                               Furt1_Furt1heso1=groupF_urt1-groupF_urt1heso1,
                               levels=design) 
}

print("test6")


dge = estimateGLMCommonDisp(cds, design)
dge = estimateGLMTagwiseDisp(dge, design) #each gene gets assigned 

print("test7")

pdf( paste0(tissue,"_MDS_plot_bcv.pdf") , width = 7 , height = 7 ) # in inches
plotBCV(dge)
dev.off()
print("test8")


# Fit data to the glm model
fit <- glmFit(dge,design)

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
  #if (length(geneInfo$description)==0) {stop(sprintf("rowIDs from count table are different than IDs from \"%s\" from ensembl\nChoose another feature ID, or another dataset from biomaRt, or check that your IDs are in the format as the one from ensembl", filter))}
  if (length(geneInfo$description)==0) {
    print(geneInfo)
    print("(########)")
  }
  
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

# Pairwise comparison (by defaut use genewise dispersion)
my_test <- function(outdiff) {
  A <- basename(outdiff)
  print(A)
  rowtoprint <- as.data.frame(data.frame(comparison=as.character(), up_genes=as.numeric(), down_genes=as.numeric()) )          
  x <- glmLRT(fit, contrast=my.contrasts[,A])
  y<-topTags(x, n = nrow(x$table))$table
  
  dt_significant <- decideTestsDGE( x, adjust.method="BH", p.value=0.05)
  
  names_sig <- rownames( dge )[ as.logical( dt_significant )]
  
  pdf( file.path(outdiff, paste0(A,"_plotsmear.pdf")) , width = 7 , height = 7 ) # in inches
  plotSmear( x, de.tags = names_sig )
  dev.off()
  
  table_up <- subset(y, y$FDR<0.05&y$logFC>0)
  table_up2 <- subset(y, y$FDR<0.05&y$logFC>0)
  print(table_up)
  print("###")
  if (nrow(table_up)!=0) {
    table_up <- try(addBiomaRtAnnotation(table_up, features="ensembl_gene_id"))
  }

  
  table_down <- subset(y, y$FDR<0.05&y$logFC<0)
  print(table_down)
  print("###")
  if (nrow(table_down)!=0) {
    table_down <- try(addBiomaRtAnnotation(table_down, features="ensembl_gene_id"))
  }
  
  #print(paste("up", nrow(table_up)), sep="\t")
  #print(paste("down", nrow(table_down)), sep="\t")
  rowtoprint <- as.data.frame(data.frame(comparison=as.character(A), up_genes=as.numeric(nrow(table_up)), down_genes=as.numeric(nrow(table_down)))) 
  write.table(table_up, file.path(outdiff, paste(i, "up.txt", sep="_")), quote=FALSE, sep="\t", col.names=NA)
  write.table(table_down, file.path(outdiff,paste(i, "down.txt", sep="_")), quote=FALSE, sep="\t", col.names=NA)
  write.table(dt_significant, file.path(outdiff,paste(i, "dt_signif.txt", sep="_")), quote=FALSE, sep="\t", col.names=NA)
  return(rowtoprint)
  
}


recap <- as.data.frame(data.frame(comparison=as.character(), up_genes=as.numeric(), down_genes=as.numeric()) )          
for (i in colnames(my.contrasts)) {
  print(i)
  outdiff=file.path(outdir,i)
  dir.create(outdiff, showWarnings = F)
  tab <- my_test(outdiff)
  recap <- as.data.frame(rbind(recap, tab))
  colnames(recap) <- c("comparison", "up_genes", "down_genes")
}

write.table(recap, paste0(tissue,"_recap_number_withbatch.txt"), quote=FALSE, sep="\t", col.names=NA)

