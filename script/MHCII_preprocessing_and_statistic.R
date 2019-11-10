# --------------
# Date:  2019-10-30 17:21:50 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project:This is for MHCII_preprocessing and statistic
# 
#require()
basic_path = "/Users/yujijun/Documents/01-Work/04-Neonatigen/neoantigen_data/MHCII_binding_prediction/original_data"
MHCII_2008_name = paste(basic_path,"MHCII_binding_2008/H-2-IAb.txt",sep = "/")
MHCII_2008 = read.table(MHCII_2008_name,header = T,row.names = NULL)
MHCII_2008$species = "mouse"
MHCII_2008 = MHCII_2008[,c("species","row.names","ALL_ALLELE","SEQ_SEQUENCE","SEQ_AA_LEN","COMP_IC50")]
colnames(MHCII_2008) = c("species","type","allele","sequence","length","ic50")
MHCII_2008$allele = paste(MHCII_2008$type,MHCII_2008$allele,sep = "-")
MHCII_2008$inequality = "="
MHCII_2008 = MHCII_2008[,c("species","allele","sequence","length","inequality","ic50")]
head(MHCII_2008)

output_path = "/Users/yujijun/Documents/01-Work/04-Neonatigen/neoantigen_data/MHCII_binding_prediction/sample_data"
write.table(MHCII_2008,file=paste(output_path,"MHCII_2008.txt",sep = "/"),sep = "\t",quote = F,row.names = F)

#MHCII data
MHCII_2009_name = paste(basic_path,"MHCII_binding_2009/class_II_all_split_5cv/H-2-IAb_test_random_0.txt",sep = "/")
MHCII_2009 = read.table(MHCII_2009_name,header = F,row.names = NULL)
MHCII_2009 = MHCII_2009[,c("V1","V2","V5","V3","V6","V7")]
colnames(MHCII_2009) <- c("species","allele","sequence","length","inequality","ic50")
write.table(MHCII_2009,file = paste(output_path,"MHCII_2009.txt",sep="/"),quote = F,sep="\t",row.names = F)
