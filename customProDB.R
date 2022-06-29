library('customProDB')
library('biomaRt')

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
                   host="http://asia.ensembl.org/", path="/biomart/martservice", 
                   archive=FALSE)           #首先需要做的准备注释，之前在UCSC注释那里纠结了很久，居然没想过用ensembl的

annotation_path <- c('tmp/')  #放注释的文件夹
PrepareAnnotationEnsembl(mart=ensembl, 
                        annotation_path=annotation_path,
                        splice_matrix=FALSE, dbsnp=NULL,
                        COSMIC=FALSE)  #跑完这一段，注释就都以RDdata的格式放在选定的文件夹里面了

vcffile <- c('1f286646-eb2d-44fc-a628-5cc310e023de.wxs.aliquot_ensemble_masked.vcf')  #先前用maf2vcf转的一个TCGA的maf文件
vcf <- InputVcf(vcffile)

load('tmp/exon_anno.RData')
exon$chromosome_name<-paste0('chr',exon$chromosome_name)   #发现自己做得这个注释有一个小问题，就是染色体的前面缺少“chr”的字符
postable_snv <- Positionincoding(vcf[[1]], exon)
postable_snv

load('tmp/ids.RData')
load('tmp/procodingseq.RData')
load('tmp/proteinseq.RData')

txlist <- unique(postable_snv[, 'txid']) 
codingseq <- procodingseq[procodingseq[, 'tx_id'] %in% txlist,]
mtab <- aaVariation (postable_snv, codingseq)

outfile <- paste('test_snv.fasta', sep='')  #输出的突变蛋白质序列文件
OutputVarproseq(mtab, proteinseq, outfile, ids) 

#还没具体看输出的序列文件是不是合理的，反正流程还挺顺
