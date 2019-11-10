# --------------
# Date:  2019-10-30 12:00:55 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: This is for ligand dataset statistics
# 
require(dplyr)
require(readtext)
require(stringi)
require(stringr)
#####parameter######
basic_path = "/Users/yujijun/Documents/01-Work/04-Neonatigen/neoantigen_data"
ligands = paste(basic_path,"Ligands/revision_data",sep = "/")
####################
# A/B
file_list = list.files(ligands)
file_open = read.csv(paste(ligands,file_list[1],sep = "/"), header = F)
file_index = which(str_locate(file_open[,1], "^>")[,1] == 1)
file_choose = file_open[file_index,] 
file_choose_new = gsub(">","",file_choose)
file_split = str_split(file_choose_new, pattern = " ")
df_AandB = data.frame(matrix(unlist(file_split), nrow=length(file_split), byrow=T))
file_A = which(str_locate(df_AandB$X4, "^A")[,1] == 1)
file_B = which(str_locate(df_AandB$X4, "^B")[,1] == 1)
df_A = df_AandB[file_A,]
df_B = df_AandB[file_B,]
# C
file_open = read.csv(paste(ligands,file_list[2],sep = "/"), header = F)
file_index = which(str_locate(file_open[,1], "^>")[,1] == 1)
file_choose = file_open[file_index,] 
file_choose_new = gsub(">","",file_choose)
file_choose_new = gsub("\t"," ",file_choose_new)
file_split = str_split(file_choose_new, pattern = " ")
df_C <- data.frame(matrix(unlist(file_split), nrow=length(file_split), byrow=T))

#E
file_open = read.csv(paste(ligands,file_list[3],sep = "/"), header = F)
file_index = which(str_locate(file_open[,1], "^>")[,1] == 1)
file_choose = file_open[file_index,] 
file_choose_new = gsub(">","",file_choose)
file_choose_new = gsub("\t"," ",file_choose_new)
file_split = str_split(file_choose_new, pattern = " ")
df_E <- data.frame(matrix(unlist(file_split), nrow=length(file_split), byrow=T))


#G
file_open = read.csv(paste(ligands,file_list[4],sep = "/"), header = F)
file_index = which(str_locate(file_open[,1], "^>")[,1] == 1)
file_choose = file_open[file_index,] 
file_choose_new = gsub(">","",file_choose)
file_choose_new = gsub("\t"," ",file_choose_new)
file_split = str_split(file_choose_new, pattern = " ")
df_G <- data.frame(matrix(unlist(file_split), nrow=length(file_split), byrow=T))
df_G$X4 <- rep("HLA-G")

#put all df together
all_df = rbind(df_A,df_B,df_C,df_E,df_G)
colnames(all_df) = c("protein_identifer","start_pos","ligand","HLA_restriction")
all_df = all_df[c("HLA_restriction","ligand","protein_identifer","start_pos")]
rownames(all_df) = seq(1,nrow(all_df))
write.table(all_df,paste(ligands,"Ligand_HLA_summary.txt",sep = "/"),quote = F,sep="\t")

#statistic
len_df <- c(nrow(df_A),nrow(df_B),nrow(df_C),nrow(df_E),nrow(df_G))
df_name <- c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-G")
df_statistic <- data.frame(df_name,len_df)
colnames(df_statistic) <- c("HLA","Number")
write.table(df_statistic,paste(ligands,"Ligand_HLA_statistic.txt",sep = "/" ),quote = F,sep="\t")


