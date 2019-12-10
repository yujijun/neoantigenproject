#gene name:
library("org.Hs.eg.db")
select(org.Hs.eg.db, symbols, columns = c("PMID","GENENAME"))

#input original dataset
mouse_base <- read.table("./Reference/mart_export_mouse.txt",sep = "\t",header = T,fill=T)
human_base <- read.table("./Reference/mart_export_human.txt",sep = "\t",header = T,fill = T)
all_base <- rbind(mouse_base,human_base)
convertfile <- read.delim("./data_project/odd_dataset/gene.csv",sep = ",",header = T,fill = T)
convertfile$genesymbol <- "-"
#preprocessing
convertfile_uniprot <- convertfile[grep("uniprot",convertfile$Protein_IRI),]
convertfile_ncbi <- convertfile[grep("ncbi",convertfile$Protein_IRI),]
uniprot_interconvert <- read.delim("./data_project/uniprot_result.txt",stringsAsFactors = F)


#uniprot_interconvert deplication 
uniprot_interconvert <- uniprot_interconvert[!duplicated(uniprot_interconvert$From),]

uniprotID.diff <- setdiff(convertfile_uniprot$gene_id,uniprot_interconvert$From)

uniprotID.more <- all_base[all_base$UniProtKB.Gene.Name.ID  %in% uniprotID.diff,]
uniprotID.more <- uniprotID.more[,c(3,4)]
colnames(uniprotID.more) <- c("From","To")

#all convert ID 
uniprotID.all <- rbind(uniprot_interconvert,uniprotID.more)


#convertfileoutput 
rownames(convertfile_uniprot) <- convertfile_uniprot$gene_id
convertfile_uniprot[c(uniprotID.all$From),"genesymbol"] <- uniprotID.all$To
write.table(convertfile_uniprot,file = "./data_project/odd_dataset/convertfile_uniprot.txt",sep = "\t",col.names = T,quote = F)

#convertfile ncbi:
#https://biodbnet-abcc.ncifcrf.gov/db/db2dbRes.php
library(readxl)
ncbi <- read.table("./data_project/odd_dataset/NCBI.txt",sep = "\t",header = T,stringsAsFactors = F)
rownames(ncbi) <- ncbi$GI.Number
rownames(convertfile_ncbi) <- convertfile_ncbi$gene_id

genesymbol <- c()
for (i in rownames(convertfile_ncbi)){
  genesymbol <- c(genesymbol,ncbi[i,"Gene.Symbol"])
}
convertfile_ncbi$genesymbol <- genesymbol
write.table(convertfile_ncbi,file = "./data_project/odd_dataset/convertfile_ncbi.txt",sep = "\t",col.names = T,quote = F)




