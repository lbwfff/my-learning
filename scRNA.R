rm(list = ls())

suppressMessages(library(igraph))
suppressMessages(library(Seurat))
suppressMessages(library(SingleR))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(sctransform))
suppressMessages(library(scater))
suppressMessages(library(stringr))
suppressMessages(library(progress))
suppressMessages(library(edgeR))
suppressMessages(library(GGally))
suppressMessages(library(clustree))
suppressMessages(library(tidytree))
suppressMessages(library(kableExtra))
suppressMessages(library(gdata))
suppressMessages(library(enrichplot))
suppressMessages(library(clusterProfiler))
suppressMessages(library(UpSetR))
suppressMessages(library(pheatmap))
# suppressMessages(library(YQSYtools))
suppressMessages(library(openxlsx))
suppressMessages(library(DropletUtils))
suppressMessages(library(scDblFinder))

#########################################################################

# dir.create('singlecell')
# dir.create('singlecell/1.rawdata')

dirlist <- dir(path = 'singlecell/1.rawdata')
samplename <- unlist(stringr::str_split(dirlist,'_'))[seq(1,4*2,2)]
dirlist <- dir(path = 'singlecell/1.rawdata',full.names = T)
path='singlecell'

i=1
file=dirlist[[i]]
Sys.setenv ("VROOM_CONNECTION_SIZE" = 131072*1000) 
data <- vroom::vroom(file) %>% tibble::column_to_rownames(var = "...1")  

object <- CreateSeuratObject(counts = data, project = samplename[i], 
                             min.cells = 3,  min.features = 200)

# dir.create(paste0('singlecell', "/2.QC/"))

# 首先使用PercentageFeatureSet函数计算线粒体基因的百分比
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
# QC-statistics of mt.percent per cell 统计每个细胞线粒体的比率
# print(summary(object[["percent.mt"]]))
mt.percnet_summary <- summary(object[["percent.mt"]])
capture.output(mt.percnet_summary, 
               file = paste0('./singlecell', "/2.QC/",samplename[i], "/mt.percent_statitics.txt"))

# add the column of HB percentage features 计算红细胞基因比率
HB.genes_total <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBQ1")
HB_m <- match(HB.genes_total, rownames(object@assays$RNA))
HB.genes <- rownames(object@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
object[["percent.HB"]] <- PercentageFeatureSet(object, features = HB.genes)
# QC-statistics of mt.percent per cell
# print(summary(object[["percent.HB"]]))
HB.percnet_summary <- summary(object[["percent.HB"]])
capture.output(HB.percnet_summary, 
               file = paste0('singlecell', "/2.QC/",samplename[i], "/HB.percent_statitics.txt"))

# QC-statistics of UMI counts per cell 统计每个细胞用UMI去重的count值
# summary(colSums(object))
UMI_counts_summary <- summary(Matrix::colSums(object))
capture.output(UMI_counts_summary, file = paste0('singlecell', "/2.QC/", samplename[i],
                                                 "/UMI_counts_statitics.txt"))
UMI_counts_summary.show <- data.frame(Quantile=names(UMI_counts_summary),summary=as.numeric(UMI_counts_summary))
kable(UMI_counts_summary.show, "html",caption = '<center>**表3. 细胞文库大小统计**</center>') %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

write.csv(UMI_counts_summary.show,file=paste0('singlecell','/2.QC/',samplename[i],'/表3. 细胞文库大小统计.csv'))


pdf(file = paste0(path, "/2.QC/", samplename[i], "/图1.Distribution of UMI counts.pdf"))
hist(colSums(object), breaks = 100, main = "Distribution of library size", xlab = "library size", 
     ylab = "Cell number", col = "blue")
dev.off()

# QC-calculate the gene number of each cell 统计每个细胞基因的数目
gene_counts_per_cell <- as.numeric(object$nFeature_RNA)
# QC-statistics of gene counts per cell
# message("QC-calculate the gene number of each cell")
# print(summary(gene_counts_per_cell))
gene_counts_summary <- summary(gene_counts_per_cell)
capture.output(gene_counts_summary, file = paste0('singlecell', "/2.QC/",samplename[i], "/gene_counts_statitics.txt"))

gene_counts_summary.show <- data.frame(Quantile=names(gene_counts_summary),summary=as.numeric(gene_counts_summary))
kable(gene_counts_summary.show, "html",caption = '<center>**表4. 每个细胞基因数统计**</center>') %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

write.csv(gene_counts_summary.show,file=paste0('singlecell','/2.QC/',samplename[i],'/表4. 每个细胞基因数统计.csv')) 

p <- hist(gene_counts_per_cell, breaks = 200)
pdf(file = paste0('singlecell', "/2.QC/",samplename[i], "/图2.Distribution of gene counts.pdf"))
# QC-histogram of gene counts distribution
hist(gene_counts_per_cell, breaks = 200, main = "Distribution of gene counts", 
     xlab = "Gene counts", ylab = "Cell number", col = "blue")
dev.off()


# QC-visualize UMI counts, gene counts, mitochondria pct, HB pct individually
# 展示细胞的count值,基因数,线粒体基因百分比,红细胞百分比,以进行下一步的剔除异常细胞
p <- VlnPlot(object, cols = rep('blue',2),features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.01)
pdf(file = paste0('singlecell', "/2.QC/",samplename[i], "/图3.Distribution of library size, gene counts.pdf"))
print(p)
dev.off()

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                        cols = rep("black",2), pt.size = 0.01) + NoLegend()

pdf(file = paste0('singlecell', "/2.QC/",samplename[i], "/图4.Distribution of library size, gene countsin feature scatter.pdf"))
print(plot1)
dev.off()


library(DropletUtils)
bcrank <- barcodeRanks(object@assays$RNA@counts)
# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)

pdf(file = paste0('singlecell', "/2.QC/", samplename[i], "/图5.Distribution of UMI counts and rank.pdf"))
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
dev.off()

#scDblFinder看一下是否双细胞，顺序是：质控粗略过滤->运行scDblFinder->较严格过滤
sce <- as.SingleCellExperiment(object) 
sce1 <- scDblFinder(sce, dbr=0.1) 
# 3574 (11.4%) doublets called 
table(sce1$scDblFinder.class)
object@meta.data$scDb<-sce1$scDblFinder.class #


########################################################
#
#试一下是不是简单的就能合并
i<-2
file=dirlist[[i]]

data <- vroom(file) %>% tibble::column_to_rownames(var = "...1")  
#把rownames调整一下后CreateSeuratObject都快多了，后续基因名的问题也不存在了
#但是这个rownames应该不同的数据会不一样吧


object2 <- CreateSeuratObject(counts = data, project = samplename[i], 
                              min.cells = 3,  min.features = 200)

object2[["percent.mt"]] <- PercentageFeatureSet(object2, pattern = "^MT-")

merge <- merge(object, y = object2, 
                       add.cell.ids = c("LZ002", "LZ008"), 
                       project = "GSE149512") #但我要是样本很多怎么合并呢

###########################################################################
# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(merge, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {  #感觉可以在这个循环里把QC做了
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# this command creates an 'integrated' data assay
merge <- IntegrateData(anchorset = immune.anchors) #就算多个样本也可以直接一步合并

######################################################################

#感觉，两个样本就把内存用完了，需要的资源未免太多了点
#我知道他markdown代码的逻辑是什么，他是已经把图出来了然后做得report，我真是服了

# merge <- subset(merge, subset = nFeature_RNA > 200 & 
#                   # nFeature_RNA < 2500 & 
#                   percent.mt < 10) #要怎么过滤可能比想象中的要复杂一些包括双细胞率什么的
# merge <- NormalizeData(merge)
# 
# #
merge <- FindVariableFeatures(merge, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merge), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(merge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) 
#为什么基因名是序号啊，是因为我也vroom读取的原因吗

plot1 + plot2


#####################################################################
#

merge <- ScaleData(merge)
merge <- RunPCA(merge, features = VariableFeatures(object = merge))

# Examine and visualize PCA results a few different ways
print(merge[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(merge, dims = 1:2, reduction = "pca")

DimPlot(merge, reduction = "pca")

DimHeatmap(merge, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(merge, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
merge <- JackStraw(merge, num.replicate = 100)
merge <- ScoreJackStraw(merge, dims = 1:20) #这个维度数我想增加要怎么办？

JackStrawPlot(merge, dims = 1:20)

ElbowPlot(merge)

merge <- FindNeighbors(merge, dims = 1:20) #模仿了报告里的参数
merge <- FindClusters(merge, resolution = 0.1)

# Look at cluster IDs of the first 5 cells
head(Idents(merge), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
merge <- RunUMAP(merge, dims = 1:20)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(merge, reduction = "umap")
DimPlot(merge, reduction = "umap",group.by ='orig.ident') 


#我想要那个每个样本一张的umap
DimPlot(merge, reduction = "umap",split.by ='orig.ident')  #天才

saveRDS(merge, file = "./merge.rds")

# find all markers of cluster 2
# cluster2.markers <- FindMarkers(merge, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(merge, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
merge.markers <- FindAllMarkers(merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #应该用得更多的会是这一句，单看某个cluster的基因意义不大？
merge.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(merge, ident.1 = 0, logfc.threshold = 0.25, 
                                test.use = "roc", only.pos = TRUE)


VlnPlot(merge, features = c("CFD", "PTGDS"))
# you can plot raw counts as well
VlnPlot(merge, features = c("CFD", "PTGDS"), slot = "counts", log = TRUE) #它这个图总感觉特别丑

FeaturePlot(merge, features = c("CFD", "PTGDS",'C7','DLK1','HMGA1'))

merge.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(merge, features = top10$gene) + NoLegend()

#############################################################################
#细胞群的注释，SingleR

hpca.se <- celldex::HumanPrimaryCellAtlasData()
hpca.se

for_SingleR <- GetAssayData(merge, slot="data") 
merge.hesc <- SingleR(test = for_SingleR, ref = hpca.se, 
                      labels = hpca.se$label.main) #巨慢无比，这里应该是有一些数据库可选的，使用不同数据库结果也应该会有细微不同
merge.hesc

table(merge.hesc$labels,merge@meta.data$seurat_clusters)

merge@meta.data$labels <-merge.hesc$labels

print(DimPlot(merge, group.by = c("seurat_clusters", "labels"),reduction = "umap")) 
#为什么感觉聚类是聚类，做出来的singleR又是另一回事
#singlR分出来的类别比聚类还多的话，那还聚类干嘛,试试tsne呢

#自己手动注释的话可以这样
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                      "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

sample_table <- as.data.frame(table(merge@meta.data$orig.ident,merge@meta.data$label))
names(sample_table) <- c("Samples","celltype","CellNumber")

plot_sample<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  # scale_fill_manual(values=colour) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample

#分组  #因为就两个样本所以分不分组也无所谓了
# group_table <- as.data.frame(table(sce2@meta.data$group,sce2@meta.data$celltype))
# names(group_table) <- c("group","celltype","CellNumber")
# plot_group<-ggplot(sce2@meta.data,aes(x=group,fill=celltype))+
#   geom_bar(position="fill")+
#   scale_fill_manual(values=colour) + 
#   theme(panel.grid = element_blank(),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         axis.line.x = element_line(colour = "black") ,
#         axis.line.y = element_line(colour = "black") ,
#         plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
#   )+labs(y="Percentage")
# plot_group
# 
# plot_sample + plot_group + 
#   plot_layout(widths = c(2, 1))

#Seurat对象确实是可以高度自定义的一个框架，设计得非常厉害
#然后我们主要会去调整的就是meta.data这个位置
#不去考虑细胞通讯和拟时序的话，基本上对于单细胞项目的需求就是这些了
