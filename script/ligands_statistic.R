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
require(seqinr)
#####parameter######
basic_path = "./data_project"
ligands = paste(basic_path,"Ligands/revision_data",sep = "/")
file_list = list.files(ligands)
fastafile<-read.fasta(paste(ligands,file_list[6],sep = "/"), seqtype = "AA",as.string = TRUE)
fasta_anno = getAnnot(fastafile)
fasta_anno_unlist = unlist(fasta_anno)
anno_matrix <- str_split_fixed(fasta_anno_unlist," ",4)
colnames(anno_matrix) <- c(paste0("V",seq(1,4)))

#detail statistic
anno_matrix_statis <- as.data.frame(anno_matrix) %>% 
  group_by(V4) %>%
  summarise(count = n())
#Output of basic statistic
colnames(anno_matrix_statis) <- c("MHC_Allele","Count")
anno_matrix_statis$Species <- "human"
anno_matrix_statis <- anno_matrix_statis[c(3,1,2)]
colnames(anno_matrix_statis)[1] <- "Species"
write.table(anno_matrix_statis,file = paste(ligands,"Ligand_MHC_statistic_detail.txt",sep = "/"),sep = "\t",quote = F,row.names = F)

#detail statistic
anno_matrix[,4] = gsub("[*][0-9]+","",anno_matrix[,4])
anno_matrix_statis <- as.data.frame(anno_matrix) %>% 
  group_by(V4) %>%
  summarise(count = n())
#Output of basic statistic
colnames(anno_matrix_statis) <- c("MHC_type","Count")
anno_matrix_statis$Species <- "human"
anno_matrix_statis <- anno_matrix_statis[c(3,1,2)]
colnames(anno_matrix_statis)[1] <- "Species"
write.table(anno_matrix_statis,file = paste(ligands,"Ligand_MHC_statistic.txt",sep = "/"),sep = "\t",quote = F,row.names = F)


