library(tidyverse)
library(data.table)
library(clusterProfiler)
library(org.At.tair.db)
library(genekitr)
library(enrichplot)

args = commandArgs(trailingOnly=TRUE)
up_f <-  args[1]
up <- fread(up_f)

down_f <-  args[2]
down <-  fread(down_f)

norm_f <- args[3]
norm <- fread(norm_f)


out_bn <- file.path(dirname(up_f), basename(dirname(up_f)))

go_gene_list = unique(sort(norm$V1))
#norm_tbl <- fread("/Users/jpeter/Desktop/DATA/RNAseq/RNAseq072023/2023-08-24T144004_Analysis/seeds.gene_models.fCounts.txt")

list_up <- up$ensembl_gene_id
go_up <- enrichGO(list_up, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", universe = go_gene_list)
go_up_simple <- as.data.frame(simplify(go_up))
go_up_rdbl <- setReadable(go_up, OrgDb = org.At.tair.db)

write_tsv(go_up_simple, file=paste0(out_bn, "_go_up.tsv"))
#write_tsv(go_up_rdbl, file=paste0(out_bn, "_go_up_rdbl.tsv"))

list_dn <- down$ensembl_gene_id
go_dn <- enrichGO(list_dn, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", universe = go_gene_list)
go_dn_simple <- as.data.frame(simplify(go_dn))
go_dn_rdbl <- setReadable(go_dn, OrgDb = org.At.tair.db)
write_tsv(go_dn_simple, file=paste0(out_bn, "_go_dn.tsv"))
#write_tsv(go_dnp_rdbl, file=paste0(out_bn, "_go_dn_rdbl.tsv"))

print(nrow(go_dn_simple))
print(nrow(go_up_simple))
while (!is.null(dev.list()))  dev.off()

  
if (nrow(go_up_simple)>0 & nrow(go_dn_simple)>0) {
    print("OK")
    #
    plotEnrichAdv(go_up_simple, go_dn_simple,
                  plot_type = "one",
                  term_metric = "Count",
                  stats_metric = "p.adjust")
    ggsave(paste0(out_bn, "_enrichAdv_1.pdf"), width = 10)
    #dev.off()
    #pdf(paste0(out_bn, "_enrichAdv_2.pdf"), width = 10)
    plotEnrichAdv(go_up_simple, go_dn_simple,
                  plot_type = "two",
                  term_metric = "Count",
                  stats_metric = "p.adjust"
    )
    ggsave(paste0(out_bn, "_enrichAdv_2.pdf"), width = 10)
    #dev.off()
    #while (!is.null(dev.list()))  dev.off()
}
  




