# --------------
# Date:  2019-10-30 18:51:25 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: This is for neoantigene data analysis
# 
#require(package)
basic_path = "./data_project/Neoantigen_database/example_data"
file_input = read.table(paste(basic_path,"neoantigen_IEDB_varified.txt",sep = "/"),sep = "\t",header = T)
file_input = file_input[,c("Gene","Mutation.in.amino.acid","Peptide","Position.in.peptide","HLA.allele","Mutation.in.nucleotide","Chromosome","Start.position.in.chromosome","End.position.in.chromosome","IEDB.result")]
rownames(file_input) <- NULL
write.table(file_input,file=paste(basic_path,"neoantigen_IEDB_varified.txt",sep="/"),quote = F,sep = "\t",row.names = F)
