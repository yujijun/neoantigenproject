rankplot <- function(gene_list, rra_LFC,rra_genename,cutoff = 1){
  require(MAGeCKFlute)
  genelist <- rra_LFC
  names(genelist) <- rra_genename
  genelist <- sort(genelist,decreasing = T)
  genelist_new[!(names(genelist) %in% gene_list)] <- NA
  if (sum(genelist_new>cutoff,na.rm = T) > 10){

  }else{
    top = sum(genelist_new>cutoff,na.rm = T)
  }
  if(sum(genelist_new < -cutoff, na.rm = T) >10){
    genenames= names(rank(genelist_new,na.last = NA)[sum(!is.na(genelist_new)) - 10])
    bottom = which(names(genelist_new) == genenames)
    bottom = length(genelist_new) - bottom
  }else{
    bottom = sum(genelist_new < -cutoff, na.rm = T)
  }
  p2 = RankView(genelist, top=top,bottom = bottom)
  print(p2)
}

library(MAGeCKFlute)
data("rra.gene_summary")
head(rra.gene_summary)
dd.rra = ReadRRA(rra.gene_summary)
head(dd.rra)
geneList= dd.rra$LFC
names(geneList) = dd.rra$Official
genelist <- geneList

rankplot(gene_list,rra_LFC = dd.rra$LFC,rra_genename = dd.rra$Official,cutoff = 1)


