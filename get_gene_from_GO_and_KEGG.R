library(GSEABase)
library(GSVA)
library(tidyverse)
library(org.Mm.eg.db)

##############################
#来自GO

test <- get(GOarray$ID[i], org.Mm.egGO2ALLEGS) 
test <- mget(test, org.Mm.egSYMBOL) %>% unlist() #这样就够了


##############################
#批量提取并进行ssGSEA分析


write.gmt <- function(gs,file){
  sink(file)
  lapply(names(gs), function(i){
    cat( paste(c(i,'tmp',gs[[i]]),collapse='\t') )
    cat('\n')
  })
  sink()
}

GOarray<-data.frame(ID=c('GO:1903409','GO:0006954','GO:0006915','GO:0006809','GO:0006979',
                         'GO:0045730','GO:0055074','GO:0005432','GO:0005757','GO:0046466',
                         'GO:0001778','GO:0004972'),
                    Term=c('reactive oxygen species biosynthetic process',
                           'inflammatory response','apoptotic process',
                           'nitric oxide biosynthetic process','response to oxidative stress',
                           'respiratory burst','calcium ion homeostasis',
                           'calcium:sodium antiporter activity',
                           'mitochondrial permeability transition pore complex',
                           'membrane lipid catabolic process','plasma membrane repair',
                           'NMDA glutamate receptor activity'))

cache<-rna[1,]

for (i in 1:nrow(GOarray)){
  
  test <- get(GOarray$ID[i], org.Mm.egGO2ALLEGS) 
  test <- mget(test, org.Mm.egSYMBOL) %>% unlist() #得到了genename
  test<-list('test'=test)
  
  write.gmt(test,'test.gmt')
  
  geneSet=getGmt('test.gmt',
                 geneIdType=SymbolIdentifier())
  ssgseaScore=gsva(as.matrix(rna), geneSet, 
                   method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
  
  cache<-rbind(cache,ssgseaScore)
  rownames(cache)[i+1]<-c(GOarray$Term[i])
}

########################
#来自KEGG

path <- KEGGREST::keggGet("mmu04216")
gene.info <- path[[1]]$GENE
genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))
gene.symbol <- genes[1:length(genes)%%3 == 2]

#######################
#想要做ssGSEA的话也和上面同理，不过这里就懒得写了















                       
