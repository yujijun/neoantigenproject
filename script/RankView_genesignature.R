library(dplyr)
library(stringi)
library(stringr)
library(tibble)
library(crisprvarified)
#gene list
basic_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/data/"
gene_list = read.table(paste(basic_path,"DGElist_zexian.txt",sep = "/"))
input_meta = read.table(paste(basic_path,"QC_meta.txt",sep = "/"),header = T)
input_LFC = read.table(paste(basic_path,"LFC_maxscale.txt",sep = "/"))
input_LFC_genesignature <- input_LFC %>% filter(rownames(input_LFC) %in% as.vector(gene_list$V1))
rownames(input_LFC_genesignature) <- rownames(input_LFC)[rownames(input_LFC) %in% as.vector(gene_list$V1)]
input_LFC_genesignature <- add_column(input_LFC_genesignature,gene_name=rownames(input_LFC_genesignature),.before="Burr_PD_L1_CMTM6_INFr_LFC")

cohort <- colnames(input_LFC)
output_path <-"/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/CRISPR_output_plot/RankView"
gene_list <- gene_list$V1
rra_gene <- rownames(input_LFC)

for (i in seq(1,dim(input_LFC)[2])){
  rra_LFC <- input_LFC[,i]
  p <- rankplot(gene_list,rra_LFC = rra_LFC,rra_genename = rra_gene)
  figure_output = paste(cohort[i],".png",sep = "")
  png(paste(output_path,figure_output,sep = "/"),height = 10,width = 12,units ="cm",res=150)
  print(p)
  dev.off()
}

#single gene(STAT1)
cohort <- colnames(input_LFC)
output_path <-"/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/CRISPR_output_plot/RankView"
gene_list <- "STAT1"
rra_gene <- rownames(input_LFC)

for (i in seq(1,dim(input_LFC)[2])){
  rra_LFC <- input_LFC[,i]
  p <- rankplot(gene_list,rra_LFC = rra_LFC,rra_genename = rra_gene)
  figure_output = paste(paste(cohort[i],"_","STAT1",sep = ""),".png",sep = "")
  png(paste(output_path,figure_output,sep = "/"),height = 10,width = 12,units ="cm",res=150)
  print(p)
  dev.off()
}


