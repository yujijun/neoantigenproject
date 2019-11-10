# --------------
# Date:  2019-10-30 15:06:16 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: This for MHCI data preprocessing 
# 
require(dplyr)
basic_path = "/Users/yujijun/Documents/01-Work/04-Neonatigen/neoantigen_data/MHCI_binding_prediction/original_data/"
file_list = list.files(basic_path)
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

#version:2006
MHCI_2006 = SMM_2006[,c("row.names","species","allele","sequence","length","inequality","ic50")]
colnames(MHCI_2006) = c("species","type","allele","sequence","length","inequality","ic50")
MHCI_2006$allele = gsub("[*]","-",MHCI_2006$allele)
MHCI_2006$allele = paste(MHCI_2006$type,MHCI_2006$allele,sep = "-")
MHCI_2006 = MHCI_2006[,c("species","allele","sequence","length","inequality","ic50")]

#version:2009
MHCI_2009 = benchmark_MHCI_2009[,c("species","mhc","sequence","peptide_length","cv","inequality","meas")]
colnames(MHCI_2009) = c("species","allele","sequence","length","cv","inequality","ic50")
MHCI_2009 = MHCI_2009[,c("species","allele","sequence","length","inequality","ic50")]

#version:2013
MHCI_2013 = benchmark_MHCI_2013[,c("species","mhc","sequence","peptide_length","inequality","meas")]
#MHCI_2013$meas = 10^(MHCI_2013$meas)
colnames(MHCI_2013) = c("species","allele","sequence","length","inequality","ic50")
MHCI_2013$allele = gsub("[*]","-",MHCI_2013$allele) #replace the * to "-", make MHCI2009 == MHCI2013

#statistic of data 
MHCI_2006_sta = as.data.frame(table(MHCI_2006$species));colnames(MHCI_2006_sta) <- c("species","Frequency")
group_2006 <- MHCI_2006 %>% group_by(species,allele)
MHCI_2006_sta_detail = as.data.frame(summarise(group_2006,count= n()))

MHCI_2009_sta = as.data.frame(table(MHCI_2009$species));colnames(MHCI_2009_sta) <- c("species","Frequency")
group_2009 <- MHCI_2009 %>% group_by(species,allele)
MHCI_2009_sta_detail = as.data.frame(summarise(group_2009,count= n()))

MHCI_2013_sta = as.data.frame(table(MHCI_2013$species));colnames(MHCI_2013_sta) <- c("species","Frequency")
group_2013 <- MHCI_2013 %>% group_by(species,allele)
MHCI_2013_sta_detail = as.data.frame(summarise(group_2013,count= n()))


#output all result
output_path = "/Users/yujijun/Documents/01-Work/04-Neonatigen/neoantigen_data/MHCI_binding_prediction/sample_data"
write.table(MHCI_2006, file=paste(output_path, "benchmark_MHCI_2006.txt",sep = "/"),quote = F,sep = "\t",row.names = F)
write.table(MHCI_2009, file=paste(output_path, "benchmark_MHCI_2009.txt",sep = "/"),quote = F,sep = "\t",row.names = F)
write.table(MHCI_2013, file=paste(output_path, "benchmark_MHCI_2013.txt",sep = "/"),quote = F,sep = "\t",row.names = F)
write.table(MHCI_2006_sta, file=paste(output_path, "benchmark_MHCI_2006_statistic.txt",sep = "/"),quote = F,sep = "\t")
write.table(MHCI_2009_sta, file=paste(output_path, "benchmark_MHCI_2009_statistic.txt",sep = "/"),quote = F,sep = "\t")
write.table(MHCI_2013_sta, file=paste(output_path, "benchmark_MHCI_2013_statistic.txt",sep = "/"),quote = F,sep = "\t")
write.table(MHCI_2006_sta_detail, file=paste(output_path, "benchmark_MHCI_2006_statistic_detail.txt",sep = "/"),quote = F,sep = "\t")
write.table(MHCI_2009_sta_detail, file=paste(output_path, "benchmark_MHCI_2009_statistic_detail.txt",sep = "/"),quote = F,sep = "\t")
write.table(MHCI_2013_sta_detail, file=paste(output_path, "benchmark_MHCI_2013_statistic_detail.txt",sep = "/"),quote = F,sep = "\t")



