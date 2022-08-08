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
