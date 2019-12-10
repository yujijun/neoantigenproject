library(tibble)
library(stringi)
library(stringr)
odd.pep.path <- "./data_project/Odd-peptide"

#CancerPPD dataset
CancerPPD <- read.csv(paste(odd.pep.path,"CancerPPD_CNRD_strange_peptide.csv",sep = "/"),stringsAsFactors = F)
CancerPPD_new <- add_column(CancerPPD, peptide_change = c(rep("change",nrow(CancerPPD))),peptide_type = c(rep("No_change",nrow(CancerPPD))),.after = "Peptide")

CancerPPD_new$peptide_change[grep("lw",CancerPPD_new$Peptide)] <- as.character(CancerPPD_new$Peptide[grep("lw",CancerPPD_new$Peptide)])
CancerPPD_new$peptide_type[grep("lw",CancerPPD_new$Peptide)] <- "Other"
CancerPPD_new$peptide_change[grep("-",CancerPPD_new$Peptide)] <- as.character(CancerPPD_new$Peptide[grep("-",CancerPPD_new$Peptide)])
CancerPPD_new$peptide_type[grep("-",CancerPPD_new$Peptide)] <- "Other"
CancerPPD_new$peptide_change[grep("[*]",CancerPPD_new$Peptide)] <- as.character(CancerPPD_new$Peptide[grep("[*]",CancerPPD_new$Peptide)])
CancerPPD_new$peptide_type[grep("[*]",CancerPPD_new$Peptide)] <- "Other"
CancerPPD_new$peptide_change[grep("rrr",CancerPPD_new$Peptide)] <- as.character(CancerPPD_new$Peptide[grep("rrr",CancerPPD_new$Peptide)])
CancerPPD_new$peptide_type[grep("rrr",CancerPPD_new$Peptide)] <- "Other"

CancerPPD_new$peptide_change[grep("change",CancerPPD_new$peptide_change)] <- as.character(str_to_upper(CancerPPD_new$Peptide[grep("change",CancerPPD_new$peptide_change)]))
write.table(CancerPPD_new,file = "./data_project/Odd-peptide/CancerPPD_CNRD_strange_peptide_v2.tsv",quote = F, row.names = F,sep = "\t")

change_fun <- function(inputfilename,outputfilename){
  df <- read.csv(paste(odd.pep.path,inputfilename,sep = "/"),stringsAsFactors = F)
  df_new <- add_column(df, peptide_change = c(rep("change",nrow(df))),peptide_type = c(rep("No_change",nrow(df))),.after = "Epitope_Description")
  df_new$peptide_change[grep("[+]",df_new$Epitope_Description)] <- 
    str_split_fixed(df_new$Epitope_Description[grep("[+]",df_new$Epitope_Description)],"[+]",n = 2)[,1]
  df_new$peptide_type[df_new$peptide_change == "change"] = "other"
  df_new$peptide_change[df_new$peptide_change == "change"] = df_new$Epitope_Description[df_new$peptide_change == "change"]
  write.table(df_new,file = paste(odd.pep.path,outputfilename,sep = "/"),quote = F, row.names = F,sep = "\t")
}
change_fun("mhc_ligand_full_select_col.csvstrange_peptide.csv","mhc_ligand_full_select_col.csvstrange_peptide_v2.tsv")
#mhc_ligand 
# mhc_ligand <- read.csv(paste(odd.pep.path,"mhc_ligand_full_select_col.csvstrange_peptide.csv",sep = "/"),stringsAsFactors = F)
# mhc_ligand_new <- add_column(mhc_ligand, peptide_change = c(rep("change",nrow(mhc_ligand))),peptide_type = c(rep("No_change",nrow(mhc_ligand))),.after = "Epitope_Description")
# mhc_ligand_new$peptide_change[grep("[+]",mhc_ligand_new$Epitope_Description)] <- 
#   str_split_fixed(mhc_ligand_new$Epitope_Description[grep("[+]",mhc_ligand_new$Epitope_Description)],"[+]",n = 2)[,1]
# mhc_ligand_new$peptide_type[mhc_ligand_new$peptide_change == "change"] = "other"
# mhc_ligand_new$peptide_change[mhc_ligand_new$peptide_change == "change"] = mhc_ligand_new$Epitope_Description[mhc_ligand_new$peptide_change == "change"]
# write.csv(mhc_ligand_new,file = "./data_project/Odd-peptide/mhc_ligand_full_select_col.csvstrange_peptide_v2.csv",quote = F, row.names = F)


#tcell_full
tcell <- read.csv(paste(odd.pep.path,"tcell_full_v3_select_col.csvstrange_peptide.csv",sep = "/"),stringsAsFactors = F)
change_fun <- function(inputfilename,outputfilename){
  df <- read.csv(paste(odd.pep.path,inputfilename,sep = "/"),stringsAsFactors = F)
  df_new <- add_column(df, peptide_change = c(rep("change",nrow(df))),peptide_type = c(rep("No_change",nrow(df))),.after = "Epitope_Description")
  df_new$peptide_change[grep("[+]",df_new$Epitope_Description)] <- 
    str_split_fixed(df_new$Epitope_Description[grep("[+]",df_new$Epitope_Description)],"[+]",n = 2)[,1]
  df_new$peptide_type[df_new$peptide_change == "change"] = "other"
  df_new$peptide_change[df_new$peptide_change == "change"] = df_new$Epitope_Description[df_new$peptide_change == "change"]
  write.table(df_new,file = paste(odd.pep.path,outputfilename,sep = "/"),quote = F, row.names = F,sep = "\t")
}
change_fun("tcell_full_v3_select_col.csvstrange_peptide.csv","tcell_full_v3_select_col.csvstrange_peptide_v2.tsv")

#Neoantigen 
neoantigen <- read.csv(paste(odd.pep.path,"Neoantigen_CNRD.csvstrange_peptide.csv",sep = "/"),stringsAsFactors = F)
neoantigen <- add_column(neoantigen,peptide_change = c(rep("change",nrow(neoantigen))),peptide_type = c(rep("No_change",nrow(neoantigen))),.after = "Peptide")
copy_neo <- neoantigen
neoantigen$peptide_change = gsub(" *\\(.*?\\) *", "",neoantigen$Peptide)
write.table(neoantigen,file = paste(odd.pep.path,"Neoantigen_CNRD.csvstrange_peptide_v2.tsv",sep = "/"),quote = F, row.names = F,sep = "\t")

#TCGA
TCGA <- read.csv(paste(odd.pep.path,"TCGA_TCR_seq_CNRD_strange_peptide.csv",sep = "/"),stringsAsFactors = F)


