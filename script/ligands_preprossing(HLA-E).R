# --------------
# Date:  2019-10-30 10:21:39 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: This is for Ligands data preprocessing for HLA-E0101
# 
require(dplyr)
require(stringi)
require(stringr)
#####parameter######
basic_path = "./data_project"
ligands = paste(basic_path,"Ligands/original_data",sep = "/")
####################
input_data = read.table(paste(ligands,"Ligand_HLA-E0101.P0A1D4.tab.txt",sep = "/"),header = F,stringsAsFactors = F)
input_data_pos <- input_data %>%
  filter(V2==1) %>% 
  select(V1) %>%
  .[,1]

identifer = ">sp|P0A1D4|CH60_SALTI"
allele_restriction = "E0101"
P0A1D4_seq = "MAAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTTTATVLAQSIITEGLKAVAAGMNPMDLKRGIDKAVAAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVGKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLLIIAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKDTTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEGVVAGGGVALIRVASKIADLKGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEYGNMIDMGILDPTKVTRSALQYAASVAGLMITTECMVTDLPKSDAPDLGAAGGMGGMGGMGGMM"

str_loc <- str_locate(P0A1D4_seq,input_data_pos)
input_data_pos_revision <- input_data_pos[!is.na(str_loc[,1])]
str_loc_revision <- str_loc[!is.na(str_loc[,1]),][,1]

First_column <- paste(identifer,str_loc_revision,input_data_pos_revision, allele_restriction)
All_file <- paste(First_column, P0A1D4_seq, sep = "\n")
write(All_file, file=paste(ligands,"Ligand_HLA-E.fasta.txt",sep = "/"))



