seq_test <- CNRD_MHC_Ligand$Peptide[1:3000]
write(t(as.character(unique(seq_test))), file = "/Users/yujijun/Desktop/test_seq.txt")
library(seqinr)
library(stringi)
library(stringr)
seq <- as.character(unique(seq_test))
seq <- seq[-grep("[+]",seq)]
seq <- seq[sapply(seq, str_length) > 8]
names <- paste("name",seq(1,1088),sep = "_")
write.fasta(sequences = as.list(seq),names = names, file.out = "/Users/yujijun/Desktop/test_seq.txt")
#write.fasta(c("AAA", "CCC"), names=c("a", "b"), file.out="/Users/yujijun/Desktop/foo.fa")
