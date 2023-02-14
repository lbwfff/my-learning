library(ChIPpeakAnno)
#####为什么要去折腾ChIPpeakAnno呢，因为我觉得homer的注释会有偏移，a基因的头部注释到b基因的尾巴这样的现象，很难理解为什么会这样
#####ChIPpeakAnno的注释在直观上是没有偏移了，但可视的结果却没有变得更好，头大

annoData <- toGRanges('~/biodata/annotation/gencode.v35.annotation.sorted.gtf',format=c('GTF'))##自己的注释文件
annoData[1:2]

annoData@ranges@NAMES<-annoData@elementMetadata@listData$gene_id##就是ensemble_id
table(annoData@elementMetadata@listData$type)
annoData<-annoData[annoData@elementMetadata@listData$type=='gene',]##不改的话，一个peak会对应到多个注释上，看实际需求吧

peaks <- toGRanges('~/biodata/newshTRA2A/DEG/test2/DiffMod.bed', format="BED")
seqlevelsStyle(peaks) <- seqlevelsStyle(annoData)
anno <- annotatePeakInBatch(peaks, AnnotationData=annoData,
                            output="overlapping",ignore.strand = FALSE)##默认是不考虑正负链，这里改为了考虑
anno[1:2]
anno<-anno[!is.na(anno@elementMetadata@listData$feature),]

annoData<-annoData[match(anno@elementMetadata@listData$feature,annoData@elementMetadata@listData$gene_id),]
anno$gene_type<-annoData$gene_type
anno$gene_name<-annoData$gene_name

head(anno)

write.csv(anno,file = 'ann.csv')

ann<-read.csv('ann.csv')
unique.ann <- ann[which(!duplicated(paste(ann$seqnames,ann$start,ann$end,sep = ":"))),]

meta<-read.csv('~/biodata/newshTRA2A/DEG/test2/DiffMod.csv')

meta<-meta[match(paste0(unique.ann$seqnames,unique.ann$start,unique.ann$end),paste0(meta$chr,meta$chromStart+1,meta$chromEnd)),]
identical(paste0(meta$chr,meta$chromStart+1,meta$chromEnd), paste0(unique.ann$seqnames,unique.ann$start,unique.ann$end))

unique.ann<-cbind(unique.ann,meta[14:22])
unique.ann<-unique.ann[,-1]
write.csv(unique.ann,file = 'ann.bind.csv')

##############################################################################################################
#我居然没记过这个venn的代码吗

shTRA2A_below <- toGRanges('~/share/MERIP/exomepeak/manorm/T1_NC/output_filters/adjust_shTRA2A_below.bed', format="BED", header=FALSE)
shMETTL3_below <- toGRanges('~/share/MERIP/exomepeak/manorm/M2_NC/output_filters/adjust_shMETTL3_below.bed', format="BED", header=FALSE)

shMETTL3_below<-unique(shMETTL3_below)
shTRA2A_below<-unique(shTRA2A_below)

ol <- findOverlapsOfPeaks(shMETTL3_below,shTRA2A_below,ignore.strand = F)

makeVennDiagram(ol,fill=c(met.brewer('Hiroshige',n=4)[c(1,3)]),
                col=c("#D55E00", "#009E73"), ignore.strand=F,
                cat.col=c("#D55E00", "#0072B2"),cat.cex =1.2)

####################################################################################################
#当我们做RIP-seq时我们会比较好奇peak在转录本上的未知信息

library(ChIPpeakAnno)

annoData <- toGRanges('../gencode/gencode.v43.basic.annotation.gff3',format=c('GFF'))##自己的注释文件
cel_feature<-c('gene','CDS','five_prime_UTR','three_prime_UTR','stop_codon') #我只想要这几个特征

annoData@ranges@NAMES<-paste0(annoData@elementMetadata@listData$gene_id,'-',annoData@elementMetadata@listData$type)
annoData<-annoData[annoData@elementMetadata@listData$type %in% cel_feature,]
table(annoData@elementMetadata@listData[["type"]])
annoframe<-as.data.frame(annoData,row.names = 1:length(annoData@elementMetadata@listData$gene_id))
annoframe$adj_id<-paste0(annoframe$gene_id,'-',annoframe$type)

annoData@elementMetadata@listData$type<-factor(rep('gene',length(annoData@elementMetadata@listData$type)))
#想的是对软件做一个欺骗，然后可以得到想要的特征信息

filelist<-c('T1','M2','NC')

for (i in filelist){
  print(i)
  path<-paste0(i,'.bed')
  peaks <- read.table(path,sep = '\t')
  path<-paste0(i,'_adj.bed')
  write.table(peaks[,1:6],file = path,row.names = F,col.names = F,quote = F,sep = '\t')
  peaks <- toGRanges(path, format="BED")
  seqlevelsStyle(peaks) <- seqlevelsStyle(annoData)

  anno <- annotatePeakInBatch(peaks, AnnotationData=annoData,
                            output="overlapping",ignore.strand = FALSE)
  
  anno<-anno[!duplicated(anno@elementMetadata@listData$feature)]
  anno<-anno[!is.na(anno@elementMetadata@listData$feature),]
  anno<-as.data.frame(anno)

  matchanno<-annoframe[match(anno$feature,annoframe$adj_id),]
  anno$gene_type<-matchanno$gene_type
  anno$gene_name<-matchanno$gene_name

  f <-function(x) unlist(strsplit(x['feature'],'-'))[2]
  anno$feartre_type<-apply(anno,1,f)
  
  anno<-anno[!duplicated(paste0(anno$peak,'_',anno$gene_name,'_',anno$feartre_type)),]
  
  table(anno$feartre_type)
  table(anno$gene_type)

  path<-paste0(i,'ann.csv')
  write.csv(anno,file = path)
  }



