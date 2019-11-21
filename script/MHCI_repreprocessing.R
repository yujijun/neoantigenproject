require(dplyr)
require(readtext)
require(stringi)
require(stringr)
require(seqinr)
require(tibble)
basic_path = "./data_project/MHCI_binding_prediction/revision_data/"
file_list = list.files(basic_path)


benchmark_MHCI_2006 <- read.table(paste(basic_path,file_list[3],sep = "/"),row.names = NULL,header = T)
MHC_Type <- gsub("[*][0-9a-zA-Z]+","",benchmark_MHCI_2006$MHC_Allele)
benchmark_MHCI_2006 <- benchmark_MHCI_2006 %>% add_column(MHC_Type = MHC_Type,.before=2)

benchmark_MHCI_2009 <- read.table(paste(basic_path,file_list[6],sep = "/"),row.names = NULL,header = T)
benchmark_MHCI_2009$MHC_Allele <- gsub("ELA-A[*]1","ELA-A1", benchmark_MHCI_2009$MHC_Allele)
MHC_Type <- gsub("[*][0-9a-zA-Z]+","",benchmark_MHCI_2009$MHC_Allele)
benchmark_MHCI_2009 <- benchmark_MHCI_2009 %>% add_column(MHC_Type = MHC_Type,.before=2)
benchmark_MHCI_2009 <- benchmark_MHCI_2009[,-6]

benchmark_MHCI_2013 <- read.table(paste(basic_path,file_list[9],sep = "/"),row.names = NULL,header = T)
MHC_Type <- gsub("[*][0-9a-zA-Z]+","",benchmark_MHCI_2013$MHC_Allele)
benchmark_MHCI_2013 <- benchmark_MHCI_2013 %>% add_column(MHC_Type = MHC_Type,.before=2)

output_path <- "./data_project/MHCI_binding_prediction/revision_data/"
write.table(benchmark_MHCI_2006, file=paste(output_path, "benchmark_MHCI_2006.txt",sep = "/"),quote = F,sep = "\t",row.names = F)
write.table(benchmark_MHCI_2009, file=paste(output_path, "benchmark_MHCI_2009.txt",sep = "/"),quote = F,sep = "\t",row.names = F)
write.table(benchmark_MHCI_2013, file=paste(output_path, "benchmark_MHCI_2013.txt",sep = "/"),quote = F,sep = "\t",row.names = F)
