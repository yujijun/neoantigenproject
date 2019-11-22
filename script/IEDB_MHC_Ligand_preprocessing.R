# --------------
# Date:  2019-11-21 15:16:19 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: Explore and preporcessing IEDB all datasets
# 
require(dplyr)
basic_path = "/Users/yujijun/Documents/01-Work/04-Neonatigen/neoantigenproject/data_project/IDEB_allfile"
output_path = "/Users/yujijun/Documents/01-Work/04-Neonatigen/neoantigenproject/data-raw"
#### preprocessing MHC_Ligand datasets####
MHC_Ligand <- read.table(paste(basic_path,"mhc_ligand_full.csv",sep = "/"),nrows = 4000,sep = ",")
MHC_Ligand <- MHC_Ligand[,-which(as.vector(MHC_Ligand[1,]) == "Related Object")]
MHC_Ligand <- MHC_Ligand[,-which(as.vector(MHC_Ligand[1,]) == "In vivo Process")]
MHC_Ligand <- MHC_Ligand[,-which(as.vector(MHC_Ligand[1,]) == "In vitro Process")]
MHC_Ligand_Homo <- MHC_Ligand %>% filter(V20=="Homo sapiens")
colnames(MHC_Ligand_Homo) <- unname(unlist(MHC_Ligand[2,]))
write.table(MHC_Ligand_Homo,file = paste(output_path,"IEDB_MHC_Ligand_Homo.txt",sep = "/"),sep = "\t",col.names = T,row.names = F,quote = F)       

#### preprocessing Antigen datasets####
Antigen <- read.table(paste(basic_path,"antigen_full_v3.csv",sep = "/"),sep = ",",fill = TRUE,comment.char = "")
write.table(Antigen,file=paste(output_path,"IEDB_Antigen.txt",sep = "/"),sep = "\t",col.names = T,row.names = F,quote = F)      


#### preprocessing epitope datasets####
Epitope <- read.table(paste(basic_path,"epitope_full_v3.csv",sep = "/"),nrows = 1000,sep = ",",fill = TRUE)
Epitope <- Epitope[,-which(as.vector(Epitope[1,]) == "Related Object")]
write.table(Epitope,file=paste(output_path,"IEDB_Epitope.txt",sep = "/"),sep = "\t",col.names = F,row.names = F,quote = F)      

####preprocessing Tcell full epitopes dataset####
Tcell <- read.table(paste(basic_path,"tcell_full_v3.csv",sep = "/"),nrows = 1000,sep = ",",fill = TRUE)
Tcell <- Tcell[,-which(as.vector(Tcell[1,]) == "Related Object")]
write.table(Tcell,file = paste(output_path,"IEDB_Tcell.txt",sep = "/"),sep = "\t",col.names = F,row.names = F,quote = F)


####preprocessing receptor epitopes dataset All about B cell####
Receptor <- read.table(paste(basic_path,"receptor_full_v3.csv",sep = "/"),nrows = 1000,sep = ",",fill = TRUE)

