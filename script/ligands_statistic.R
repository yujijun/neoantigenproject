# --------------
# Date:  2019-10-30 12:00:55 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: 
# 
require(dplyr)
require(readtext)
require(stringi)
require(stringr)
#####parameter######
basic_path = "./data_project"
ligands = paste(basic_path,"Ligands/revision_data",sep = "/")
file_list = list.files(ligands)
fastafile<-read.fasta(paste(ligands,file_list[5],sep = "/"), seqtype = "AA",as.string = TRUE)
fasta_anno = getAnnot(fastafile)
fasta_anno_unlist = unlist(fasta_anno)
anno_matrix <- str_split_fixed(fasta_anno_unlist," ",5)
colnames(anno_matrix) <- c(paste0("V",seq(1,5)))

#basic statistic
anno_matrix_statis <- as.data.frame(anno_matrix) %>% 
  group_by(V4) %>%
  summarise(count = n())

#Output of basic statistic
colnames(anno_matrix_statis) <- c("HLA Type","Count")
anno_matrix_statis$Species <- "Human"
anno_matrix_statis <- anno_matrix_statis[c(3,1,2)]
colnames(anno_matrix_statis)[1] <- "Species"
write.table(anno_matrix_statis,file = paste(ligands,"Ligand_HLA_statistic.txt",sep = "/"),sep = "\t",quote = F,row.names = F)

#detail statistic 
anno_matrix_dataframe <- as.data.frame(anno_matrix)
anno_matrix_dataframe[is.na(anno_matrix_dataframe)] <- NULL
anno_matrix_summary <- anno_matrix_dataframe %>% 
  group_by(V4,V5) %>%
  summarise(count = n())
colnames(anno_matrix_summary) <- c("HLA Type","HLA Allele","Count")
anno_matrix_summary$Species <- "Human"
anno_matrix_summary <- anno_matrix_summary[c(4,1,2,3)]
colnames(anno_matrix_summary)[1] <- "Species"
write.table(anno_matrix_summary,file = paste(ligands,"Ligand_HLA_summary.txt",sep = "/"),sep = "\t",quote = F,row.names = F)

