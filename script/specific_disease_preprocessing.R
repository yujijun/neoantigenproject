library(tibble)
library(stringi)
library(stringr)
library(readxl)
odd.dis.path <- "./data_project/Disease_uniform"
cancerppd <- read.csv(paste(odd.dis.path,"CancerPPD_CNRD_disease.csv",sep = "/"),header = T,stringsAsFactors = F)
neoantigen <- read.csv(paste(odd.dis.path,"Neoantigen_CNRD.csv_disease.csv",sep = "/"),header = T)
baseline <- read_excel(paste(odd.dis.path,"Cancername_compare.xlsx",sep = "/"))
gdc <- read.csv(paste(odd.dis.path,"gdc_manifest.kraken_info.txt_disease.csv",sep = "/"),col.names = T)
gdc$disease <- str_split_fixed(gdc$TRUE.,"-",n=2)[,2]
tissue <- c()
for (i in gdc$disease){
  tissue_i <- baseline$Tissue[which(baseline$`Study Abbreviation` == i)]
  tissue <- c(tissue,tissue_i)
}

gdc$tissue <- tissue
write.table(gdc,file = "./data_project/Disease_uniform/gdc_manifest.kraken_info.txt_disease_v2.tsv",quote = F, row.names = F,sep = "\t")

cancerppd$tissue <- gsub("Cancer","",cancerppd$Cancername)
cancerppd$tissue <- gsub("cancer","",cancerppd$tissue)
cancerppd$tissue <- gsub("Tumor","",cancerppd$tissue)
cancerppd$tissue <- gsub("Human","",cancerppd$tissue)
write.table(cancerppd,file = "./data_project/Disease_uniform/CancerPPD_CNRD_disease_v2.tsv",quote = F, row.names = F,sep = "\t")

# conform
cancerppd <- read.table(paste(odd.dis.path,"CancerPPD_CNRD_disease_v2.tsv",sep = "/"),stringsAsFactors = F,sep = "\t",header = T)
neoantigen <- read.table(paste(odd.dis.path,"Neoantigen_CNRD.csv_disease_v2.csv",sep = "/"),header = T,sep = ",")
gdc <- read.table(paste(odd.dis.path,"gdc_manifest.kraken_info.txt_disease_v2.tsv",sep = "/"),header = T,sep = "\t")

neoantigen_tissue <- str_to_lower(gsub(" ","",unique(as.vector(str_split_fixed(neoantigen$change.as,",",n = Inf)))))
gdc_tissue <- str_to_lower(gsub(" ","",gdc$tissue))
cancerppd_tissue <- str_to_lower(gsub(" ","",cancerppd$tissue))

sort(unique(c(neoantigen_tissue,gdc_tissue,cancerppd_tissue)))
