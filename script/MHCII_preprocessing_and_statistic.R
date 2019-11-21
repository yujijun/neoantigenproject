# --------------
# Date:  2019-10-30 17:21:50 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project:This is for MHCII_preprocessing and statistic
# 
require(dplyr)
require(stringi)
require(stringr)
####MHC_binding_2008####
basic_path = "./data_project/MHCII_binding_prediction/original_data/MHCII_binding_2008"
MHCII_2008_name = paste(basic_path,"H-2-IAb.txt",sep = "/")
MHCII_2008 = read.table(MHCII_2008_name,header = T,row.names = NULL)
#binding all matrix together
file_list = list.files(basic_path)
MHCII_file_name = paste(basic_path, file_list, sep="/")
for(i in MHCII_file_name[-1]){
  file = read.table(i, header=T, row.names = NULL)
  MHCII_2008 = bind_rows(MHCII_2008,file)
  print(dim(MHCII_2008))
}
#Species
species = c(rep("mouse",sum(MHCII_2008$row.names == "H-2")),rep("human",length(grep("HLA",MHCII_2008$row.names))))
MHCII_2008$Species = species
#MHC_Allele
MHC_ALlele = paste(MHCII_2008$row.names, MHCII_2008$ALL_ALLELE,sep="-")
MHCII_2008$MHC_Allele = MHC_ALlele
#Add Inequality 
MHCII_2008$Inequality = "="
# change the column order of the matrix 
MHCII_2008 = MHCII_2008[,c(6,7,4,3,8,5)]
#rename and output 
colnames(MHCII_2008) = c("Species", "MHC_Allele","Peptide","Length","Inequality","Ic50")
output_path = "./data_project/MHCII_binding_prediction/revision_data/"
write.table(MHCII_2008,file=paste(output_path,"benchmark_MHCII_2008.txt",sep = "/"),sep = "\t",quote = F,row.names = F)

#MHCII 2009
basic_path = "./data_project/MHCII_binding_prediction/original_data/MHCII_binding_2009/class_II_all_split_5cv"
list_file = list.files(basic_path,pattern = "*_random_0.txt")
list_file = paste(basic_path, list_file,sep = "/")
First = data.frame()
for(i in list_file){
  file = read.table(i, header=F, row.names = NULL)
  First = bind_rows(First,file)
  print(dim(First))
}
MHCII_2009 = First[,c("V1","V2","V5","V3","V6","V7")]
colnames(MHCII_2009) <- c("Species","MHC_Allele","Peptide","Length","Inequality","Ic50")
output_path = "./data_project/MHCII_binding_prediction/revision_data/"
write.table(MHCII_2009,file = paste(output_path,"benchmark_MHCII_2009.txt",sep="/"),quote = F,sep="\t",row.names = F)

#MHCII 2016
basic_path = "/Users/yujijun/Documents/01-Work/04-Neonatigen/neoantigenproject/data_project/MHCII_binding_prediction/original_data/MHCII_bind_NetMHCIIPan_32_2016"
list_file = list.files(basic_path)
list_file = paste(basic_path,c("test1","train1"),sep = "/")
Empty = data.frame()
for(i in list_file){
  file = read.table(i,header = F,row.names = NULL)
  Empty = bind_rows(Empty,file)
  print(dim(Empty))
}
Empty = mutate(Empty,Ic50 = exp((1-V2)*log(50000))) #calculate Ic50 by inverse function 
#HLA Allele 
Empty$V3 = gsub("^DR","HLA-DR",Empty$V3)
allele = Empty$V3
allele = gsub("_","*",allele)
allele = gsub("HLA-DPA1","HLA-DPA1*",allele)
allele = gsub("HLA-DQA1","HLA-DQA1*",allele)
allele = gsub("-DPB1","/DPB1*",allele)
allele = gsub("-DQB1","/DQB1*",allele)
Empty$V3 = allele
#Species,length and Inequality
Empty$Species = "human/mouse"
Empty$Species[grep("HLA-",Empty$V3)] = rep("human",length(grep("HLA-",Empty$V3)))
Empty$Species[!grep("HLA-",Empty$V3)] = rep("mouse",length(grep("HLA-",Empty$V3)))
Empty$Length = nchar(Empty$V1)
Empty$Inequality = "="
# choose the columns and change the name of columns, output
Empty = Empty[,c("Species","V3","V1","Length","Inequality","Ic50")]
colnames(Empty) = c("Species","MHC_Allele","Peptide","Length","Inequality","Ic50")
output_path = "./data_project/MHCII_binding_prediction/revision_data/"
write.table(Empty,file = paste(output_path,"benchmark_MHCII_2016.txt",sep="/"),quote = F,sep="\t",row.names = F)

#MHCII 2018
basic_path = "/Users/yujijun/Documents/01-Work/04-Neonatigen/neoantigenproject/data_project/MHCII_binding_prediction/original_data/MHCII_initial_benchmark_2018"
#MHCII_initial_benchmark_datasets.txt
file_input = read.table(paste(basic_path,"MHCII_initial_benchmark_datasets.tsv",sep = "/"),header=T,row.names=NULL,fill=T)
file_input = file_input[-which(file_input$allele == ""),] #delete empty column
file_input = file_input[-which(file_input$allele == "+"),]
allele =file_input$allele
allele = gsub(":","",allele)
allele = gsub("H2","H-2",allele)
unique(allele)
file_input$Species = "NULL"
file_input$Species[grep("H-2",allele)] = "mouse"
file_input$Species[grep("HLA-",allele)] = "human"
file_input$length = nchar(as.vector(file_input$peptide_sequence))
file_input$Inequality = "="

file_input_binary = file_input[which(file_input$measurement_type == "binary"),]
file_input_ic50 = file_input[-which(file_input$measurement_type == "binary"),]
file_input_ic50 = file_input_ic50[,c(7,4,3,8,9,5)]
colnames(file_input_ic50) = c("Species","MHC_Allele","Peptide","Length","Inequality","Ic50")
file_input_binary = file_input_binary[,c(7,4,3,8,5,6)]
colnames(file_input_binary) = c("Species","MHC_Allele","Peptide","Length","Measure_value","Measure_Type")

output_path = "./data_project/MHCII_binding_prediction/revision_data/"
write.table(file_input_ic50,file = paste(output_path,"benchmark_MHCII_2018_ic50.txt",sep="/"),quote = F,sep="\t",row.names = F)
write.table(file_input_binary,file = paste(output_path,"benchmark_MHCII_2018_binary.txt",sep="/"),quote = F,sep="\t",row.names = F)

  
  