# --------------
# Date:  2019-10-30 15:06:16 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: This for MHCI data preprocessing 
# 
require(dplyr)
require(readtext)
require(stringi)
require(stringr)
require(seqinr)
basic_path = "./data_project/MHCI_binding_prediction/original_data/"
file_list = list.files(basic_path)
####input data####
ANN_2006 = read.table(paste(basic_path, file_list[2],sep = "/"),row.names = NULL,header = T)
ARB_2006 = read.table(paste(basic_path, file_list[3],sep = "/"),row.names = NULL,header = T)
SMM_2006 = read.table(paste(basic_path, file_list[6],sep = "/"),row.names = NULL,header = T)
MHCI_2009 = read.table(paste(basic_path, file_list[5],sep = "/"),row.names = NULL) #there is juse mumu data 
MHCI_2013 = read.table(paste(basic_path, file_list[1],sep = "/"),row.names = NULL,header = T)
#benchmark_database 
benchmark_path = paste(basic_path,"MHCI_binding_benchmark_2009_2013/binding/bd2009.1",sep="/")
benchmark_MHCI_2009 = read.table(paste(benchmark_path,"bdata.2009.mhci.public.1.txt",sep = "/"),row.names = NULL,header = T)
benchmark_path = paste(basic_path,"MHCI_binding_benchmark_2009_2013/binding/bd2013.1",sep="/")
benchmark_MHCI_2013 = read.table(paste(benchmark_path,"bdata.20130222.mhci.public.1.txt",sep = "/"),row.names = NULL,header = T)

####data propressing####
#version:2006
MHCI_2006 = SMM_2006[,c("row.names","species","allele","sequence","length","inequality","ic50")]
colnames(MHCI_2006) = c("species","type","allele","sequence","length","inequality","ic50")
MHCI_2006$type = paste(MHCI_2006$type,MHCI_2006$allele,sep = "-")
MHCI_2006 = MHCI_2006[,-3]
colnames(MHCI_2006) = c("Species","MHC_Allele","Peptide ","Length","Inequality","Ic50")
#still have some of MHC allele format are not good enough.
H_2MHC = MHCI_2006[grep("H-2", MHCI_2006$MHC_Allele),2]
H_2MHC_cformat = gsub("H-2-","H-2*",H_2MHC)
MHCI_2006[grep("H-2", MHCI_2006$MHC_Allele),2] <- H_2MHC_cformat


#version:2009
MHCI_2009 = benchmark_MHCI_2009[,c("species","mhc","sequence","peptide_length","cv","inequality","meas")]
mhc_allele = str_split_fixed(MHCI_2009$mhc,"-",n=3)
mhc_allele = paste(mhc_allele[,1],paste(mhc_allele[,2],mhc_allele[,3],sep = "*"),sep = "-")
#handle some specific mhc_allele
specific_MHC = grep("[*]$",unique(mhc_allele),value = T)
specific_MHC = as.vector(str_split_fixed(specific_MHC,"[*]",n=2)[,1])

mhc_allele = gsub("ELA-A1[*]","ELA-A*1",mhc_allele)
mhc_allele = gsub("HLA-A11[*]","HLA-A*11",mhc_allele)
mhc_allele = gsub("HLA-A2[*]","HLA-A*02",mhc_allele)
mhc_allele = gsub("HLA-A26[*]","HLA-A*26",mhc_allele)
mhc_allele = gsub("HLA-B44[*]","HLA-B*44",mhc_allele)
mhc_allele = gsub("HLA-B7[*]","HLA-B*07",mhc_allele)
#mhc_allele = gsub("ELA-A1","ELA-A*1",mhc_allele)
MHCI_2009$mhc = mhc_allele
colnames(MHCI_2009) = c("Species","MHC_Allele","Sequence","Length","CV","Inequality","Ic50")
MHCI_2009$Species = gsub("None","horse",MHCI_2009$Species)


#version:2013
MHCI_2013 = benchmark_MHCI_2013[,c("species","mhc","sequence","peptide_length","inequality","meas")]
mhc_allele = MHCI_2013$mhc
mhc_allele = gsub(":","",mhc_allele)
mhc_allele = gsub("-2-","-2*",mhc_allele)
mhc_allele = gsub("HLA-A1","HLA-A*01",mhc_allele)
mhc_allele = gsub("HLA-A11","HLA-A*11",mhc_allele)
mhc_allele = gsub("HLA-A2","HLA-A*02",mhc_allele)
mhc_allele = gsub("HLA-A3","HLA-A*03",mhc_allele)
mhc_allele = gsub("HLA-A24","HLA-A*24",mhc_allele)
mhc_allele = gsub("HLA-A26","HLA-A*26",mhc_allele)
mhc_allele = gsub("HLA-B27","HLA-B*27",mhc_allele)
mhc_allele = gsub("HLA-B44","HLA-B*44",mhc_allele)
mhc_allele = gsub("HLA-B51","HLA-B*51",mhc_allele)
mhc_allele = gsub("HLA-B60","HLA-B*60",mhc_allele)
mhc_allele = gsub("HLA-B7","HLA-B*07",mhc_allele)
mhc_allele = gsub("HLA-B8","HLA-B*08",mhc_allele)
MHCI_2013$mhc = mhc_allele
colnames(MHCI_2013) = c("Species","MHC_Allele","Peptide","Length","Inequality","Ic50")


#statistic
MHCI_2006_sta = as.data.frame(table(MHCI_2006$Species));colnames(MHCI_2006_sta) <- c("Species","Count") #Just statistic for species
group_2006 <- MHCI_2006 %>% group_by(Species,MHC_Allele)
MHCI_2006_sta_detail = as.data.frame(summarise(group_2006,count= n())) #statistic for species and alleles

MHCI_2009_sta = as.data.frame(table(MHCI_2009$Species));colnames(MHCI_2009_sta) <- c("Species","Count")
group_2009 <- MHCI_2009 %>% group_by(Species,MHC_Allele)
MHCI_2009_sta_detail = as.data.frame(summarise(group_2009,count= n()))

MHCI_2013_sta = as.data.frame(table(MHCI_2013$Species));colnames(MHCI_2013_sta) <- c("Species","Count")
group_2013 <- MHCI_2013 %>% group_by(Species,MHC_Allele)
MHCI_2013_sta_detail = as.data.frame(summarise(group_2013,count= n()))


#output all result
output_path = "./data_project/MHCI_binding_prediction/revision_data"
write.table(MHCI_2006, file=paste(output_path, "benchmark_MHCI_2006.txt",sep = "/"),quote = F,sep = "\t",row.names = F)
write.table(MHCI_2009, file=paste(output_path, "benchmark_MHCI_2009.txt",sep = "/"),quote = F,sep = "\t",row.names = F)
write.table(MHCI_2013, file=paste(output_path, "benchmark_MHCI_2013.txt",sep = "/"),quote = F,sep = "\t",row.names = F)
write.table(MHCI_2006_sta, file=paste(output_path, "benchmark_MHCI_2006_statistic.txt",sep = "/"),quote = F,sep = "\t")
write.table(MHCI_2009_sta, file=paste(output_path, "benchmark_MHCI_2009_statistic.txt",sep = "/"),quote = F,sep = "\t")
write.table(MHCI_2013_sta, file=paste(output_path, "benchmark_MHCI_2013_statistic.txt",sep = "/"),quote = F,sep = "\t")
write.table(MHCI_2006_sta_detail, file=paste(output_path, "benchmark_MHCI_2006_statistic_detail.txt",sep = "/"),quote = F,sep = "\t")
write.table(MHCI_2009_sta_detail, file=paste(output_path, "benchmark_MHCI_2009_statistic_detail.txt",sep = "/"),quote = F,sep = "\t")
write.table(MHCI_2013_sta_detail, file=paste(output_path, "benchmark_MHCI_2013_statistic_detail.txt",sep = "/"),quote = F,sep = "\t")



