## code to prepare `DATASET` dataset goes here

#upload IEDB MHC Ligand data
IEDB_MHC_Ligand <- 
  read.table('./data-raw/IEDB_MHC_Ligand_Homo.txt',sep = "\t",quote = "")
usethis::use_data(IEDB_MHC_Ligand,overwrite = T)

IEDB_Epitope <- 
  read.table('./data-raw/IEDB_Epitope.txt',sep = "\t",fill = T,quote = "",comment.char = "")
usethis::use_data(IEDB_Epitope,overwrite = T)

IEDB_Antigen <- 
  read.table('./data-raw/IEDB_Antigen.txt',sep = "\t",header = T,fill = T,comment.char = "")
usethis::use_data(IEDB_Antigen,overwrite = T)

IEDB_Tcell <- 
  read.table('./data-raw/IEDB_Tcell.txt',sep = "\t",fill = T)
usethis::use_data(IEDB_Tcell)







require(dplyr)
experiment1 <- 
  read.table('./data-raw/Burr_2017_PDL1.txt') %>%
  head()
usethis::use_data(experiment1)