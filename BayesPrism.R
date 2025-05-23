
####
#好些学过R语言新工具了，一直在忙着乱七八糟的
#BayesPrism是一基于scRNA-seq数据去卷积的工具，使用它可以对bulk-RNA seq进行细胞类型的预测

########官方实例#####
library(BayesPrism)

load("/scratch/lb4489/biotools/BayesPrism/tutorial.gbm.rdata") #不是demo有必要搞这么大吗

dim(bk.dat)
head(rownames(bk.dat))
head(colnames(bk.dat))

bksmall<-sample(1:169, size = 20)
bk.dat<-bk.dat[bksmall,]

dim(sc.dat)
head(rownames(sc.dat))
head(colnames(sc.dat))

scsmall<-sample(1:23793, size = 4000)
sc.dat<-sc.dat[scsmall,]
cell.type.labels<-cell.type.labels[scsmall]
cell.state.labels<-cell.state.labels[scsmall]

sort(table(cell.type.labels))
sort(table(cell.state.labels)) #没太懂这个细胞状态

table(cbind.data.frame(cell.state.labels, cell.type.labels))



# plot.cor.phi (input=sc.dat,
#               input.labels=cell.state.labels,
#               title="cell state correlation",
#               cexRow=0.2, cexCol=0.2,
#               margins=c(2,2))
# 
# plot.cor.phi (input=sc.dat, 
#               input.labels=cell.type.labels, 
#               title="cell type correlation",
#               cexRow=0.5, cexCol=0.5,) #相关性热图吗？并不好用这个代码


sc.stat <- plot.scRNA.outlier(
  input=sc.dat, 
  cell.type.labels=cell.type.labels,
  species="hs", #(hs) or (mm) 
  return.raw=TRUE 
) #它还自动注释了一下基因，挺有意思的可以自己改一下图

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,
  sc.input=sc.dat, 
  cell.type.labels=cell.type.labels,
  species="hs", 
  return.raw=TRUE
)


sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)


plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
) #没太懂这个相关性是怎么实现的，不过大意就是protein coding基因在sc还是bulk数据中更加一致，可以用来分析


sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding") #只用protein coding基因分析

diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.state.labels,
                              pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=5 #number of threads
)#类似于markergene的寻找，只用这些特征基因可以增加运算的效率


sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,
                                         lfc.min=0.1) #降低了表达矩阵的基因数，只保留了markergene

myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)#构建对象


bp.res <- run.prism(prism = myPrism, n.cores=48) #没法在PC上跑,需要巨多的CPU核心

save(bp.res,file = 'BayesPrism_demo.RData')

theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type") #这个就是预测的细胞比例


###############一个自己的例子#################

 library(Seurat)
 library(data.table)
 
 dir='/scratch/lb4489/biotools/BayesPrism/DevBrain/' #scRNA-seq数据的文件夹
 samples=list.files( dir )
 samples<-samples[1:5]
 
 sceList = lapply(samples,function(pro){
   print(pro)
   ct=fread(file.path( dir ,pro),data.table = F)
   ct[1:4,1:4]
   rownames(ct)=ct[,1]
   ct=ct[,-1]
   sce=CreateSeuratObject(counts =  ct ,
                          project = gsub('.annotated_matrix.txt.gz','',gsub('^GSM[0-9]*_','',pro) ) ,
                          min.cells = 5,
                          min.features = 300)
   sce@meta.data$cell_type<-rownames(sce@meta.data)
   
   sce@meta.data$cell_type<-sapply(sce@meta.data$cell_type, function(x) unlist(strsplit(x, '[.]'))[1])
   
   return(sce)
 })
 
 sce.all=merge(x=sceList[[1]],
               y=sceList[ -1 ],
               add.cell.ids =   gsub('..annotated_matrix.txt.gz','',gsub('^GSM[0-9]*_','',samples) )  ) #合并seurat对象
 
 sceList[[1]]@raw.data
 
 sce.all_new = JoinLayers(sce.all)
 
 count<-LayerData(sce.all_new,assay="RNA",layer="counts")
 count<-as.data.frame(as.data.frame(count))#提取count值
 count<-t(count) 

#我要不要改一下cellname？有点担心空格好不好导致报错，不过看起来好像正常运行了
 
 cell.type.labels<-sce.all_new@meta.data$cell_type
 cell.state.labels<-paste0(sce.all_new@meta.data$orig.ident,'_',sce.all_new@meta.data$cell_type) #一个细胞状态还不能对应多种细胞，没有什么想填的就用样本名加细胞类型了
 
 save(count,cell.type.labels,cell.state.labels,file = '/scratch/lb4489/biotools/BayesPrism/mydata.RData' )
 load('/scratch/lb4489/biotools/BayesPrism/mydata.RData')
 
 library(DESeq2)
 
 load('/scratch/lb4489/biotools/BayesPrism/rse_gene_unfiltered.Rdata') #Bulk-RNA-seq数据
 
 rowcount<-assays(rse_gene)$counts #转换ID
 
 load('/scratch/lb4489/biotools/BayesPrism/anno.RData')
 
 gtf<-gtf[match(substr(rownames(rowcount),1,15),substr(gtf$gene_id,1,15)),]
 rownames(rowcount)<-gtf$gene_name
 rowcount<-t(rowcount)

library(BayesPrism)

 scsmall<-sample(1:42497, size = 8000) 
#还是用小一点的数据集运行，运算的效率不是很高增加细胞数的话预估时间会加很多
#不太确定应该使用多少个细胞是合适的

 sc.dat<-count[scsmall,]
 cell.type.labels<-cell.type.labels[scsmall]
 cell.state.labels<-cell.state.labels[scsmall]
 
 sort(table(cell.type.labels))
 sort(table(cell.state.labels))
 
 table(cbind.data.frame(cell.state.labels, cell.type.labels))
 
 bk.dat<-rowcount
 
 sc.stat <- plot.scRNA.outlier(
   input=sc.dat, 
   cell.type.labels=cell.type.labels,
   species="hs", #(hs) or (mm) 
   return.raw=TRUE 
 ) 
 
 bk.stat <- plot.bulk.outlier(
   bulk.input=bk.dat,
   sc.input=sc.dat, 
   cell.type.labels=cell.type.labels,
   species="hs", 
   return.raw=TRUE
 ) #有一个报错没懂为什么，不过不影响运行
 
 sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                   input.type="count.matrix",
                                   species="hs", 
                                   gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                   exp.cells=5)
 
 
 plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                  bulk.input = bk.dat)
 
 
 sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                          gene.type = "protein_coding") 
 
 diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
                               cell.type.labels=cell.type.labels,
                               cell.state.labels=cell.state.labels,
                               pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                               cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                               n.cores=5 #number of threads
 )
 
 
 sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc,
                                          stat=diff.exp.stat,
                                          pval.max=0.01,
                                          lfc.min=0.1) #降低了表达矩阵的基因数，只保留了markergene
 
 myPrism <- new.prism(
   reference=sc.dat.filtered.pc, 
   mixture=bk.dat,
   input.type="count.matrix", 
   cell.type.labels = cell.type.labels, 
   cell.state.labels = cell.state.labels,
   key=NULL,
   outlier.cut=0.01,
   outlier.fraction=0.1,
 )


 save(myPrism,file = '/scratch/lb4489/biotools/BayesPrism/BayesPrism_test.RData')

load('/scratch/lb4489/biotools/BayesPrism/BayesPrism_test.RData')

bp.res <- run.prism(prism = myPrism, n.cores=78) 

save(bp.res,file = '/scratch/lb4489/biotools/BayesPrism/BayesPrism_test_result.RData')

#########################################
#很早之前就听说过这个工具了但是一直没有用过，主要是不怎么感兴趣于细胞组成，我的研究更加分子尺度一些
#有很多基因集打分的工具也可以用于细胞类型的预测https://github.com/omnideconv/immunedeconv，这里整合了许多的算法，R包可以做这个
#为什么要用BayesPrism呢，因为这个基于基因集的工具，大部分都是在做免疫细胞的预测，其它可能就找不到工具了，二是可能只能做人类的数据的分析？
#BayesPrism可以自定义的部分就多得多了，缺点就是慢，需要的资源多
#InstaPrism据说可以得到BayesPrism一样的结果但是运行资源要少得多，有时间可以试试














   
