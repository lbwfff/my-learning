# BiocManager::install("Mfuzz")
library('Mfuzz') #这个包好像没法手动去指定时间刻度啊，这不是很呆吗，而且我需要把重复自己算平均值？

testfuzz<-TPMs[,(9:16)]
#然后手动算出四个时间刻度的平均值
testfuzz$sham<-((testfuzz$..sham.markdup.bam+testfuzz$..sham.2.markdup.bam)/2)
testfuzz$Re0<-((testfuzz$..Re.0h.markdup.bam+testfuzz$..Re.0h.2.markdup.bam)/2)
testfuzz$Re6<-((testfuzz$..Re.6h.2.markdup.bam+testfuzz$..Re.6h.markdup.bam)/2)
testfuzz$Re24<-((testfuzz$..Re.24h.markdup.bam+testfuzz$..Re.24h.2.markdup.bam)/2)

testfuzz<-as.matrix(testfuzz[,9:12])
testfuzz <- new('ExpressionSet',exprs = testfuzz) #没有很理解这一步操作具体做了什么

testfuzz<-filter.NA(testfuzz,thres = 0.25) #RNA-seq里NA值应该全部设置为0值了，索引这一步好像没什么用
testfuzz <- fill.NA(testfuzz,mode="mean")
# testfuzz <- filter.std(testfuzz,min.std=0.2)
testfuzz <- filter.std(testfuzz,min.std=1) #五万六的基因里删去三万二，好像也还是挺多的，没有很理解标准差的概念

testfuzz <- standardise(testfuzz)

m1 <- mestimate(testfuzz)

tmp <- cselection(testfuzz,m=m1,crange=seq(4,16,4),repeats=2,visu=TRUE)
#是可以自己优化簇数量的啊

cl <- mfuzz(testfuzz,c=16,m=m1)

pdf("test_mfuzz.pdf",width = 12,height = 8)
mfuzz.plot2(testfuzz,cl=cl,mfrow=c(4,4),
            colo=colors,x11=F,
            time.labels=c('sham','Re-0','Re-6','Re-24')) #颜色也是可以调整的,想要保存的话x11=F就行了
dev.off()

#colours
{
color<-c(met.brewer("Homer1")[6],met.brewer("Wissing")[5]) #第一个值为低聚集的颜色，第二个值为高聚集时的颜色
pal<-colorRampPalette(color)
colors = pal(10)
}

cl$size
head(cl$cluster)
head(cl$membership)

gene_cluster <- cbind(cl$cluster, cl$membership) 
colnames(gene_cluster)[1] <- 'cluster' 
write.csv(gene_cluster, 'gene_cluster.csv') 
