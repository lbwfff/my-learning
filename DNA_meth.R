
#主要还是了解一下ChAMP这个包吧，感觉还挺厉害的

############################################################################
setwd("01_MethCGIcluster/")
options(stringsAsFactors = FALSE)

library(data.table)
library(ChAMP)

library(ClassDiscovery)
library(pheatmap)
library(gplots)

library(ggplot2)
library(ComplexHeatmap)

# 读取原始甲基化谱
orgmeth <- fread("TCGA-STAD.methylation450.tsv") #
orgmeth <- as.data.frame(orgmeth); rownames(orgmeth) <- orgmeth[,1]; orgmeth <- orgmeth[,-1]
orgmeth[1:3,1:3]

# 读取高低风险组信息
Sinfo <- read.table("risk.train.group.txt", header = T, row.names = 1)
head(Sinfo)

# 只提取其中的40个样本，高风险20个低风险20个
colnames(orgmeth) <- substr(colnames(orgmeth),1,12)  #去除最后4个字符，为了和预后模型分组结果的ID保持一致
orgmeth[1:3,1:3]

orgmeth <- orgmeth[,colnames(orgmeth) %in% row.names(Sinfo)]
Sinfo <- Sinfo[row.names(Sinfo) %in% colnames(orgmeth), ]  # 去除没有甲基化的样本
table(Sinfo$group)  # high19, low17

orgmeth <- orgmeth[,row.names(Sinfo)]  # 甲基化数据排序，前面是low，后面是high
# write.table(orgmeth,"TCGA-STAD.methylation450_filter.tsv")

###############################################################
#做一个过滤

myFilter <- champ.filter(beta = orgmeth, # beta矩阵
                         pd = NULL,    # 样本信息，因为仅利用champ做探针过滤，所以没有必要给出
                         # 过滤CpG，SNPs以及XY
                         autoimpute= F,# 不填补缺失值
                         filterDetP= F, # 不根据探针p值过滤
                         fixOutlier= F, # 不修正离群点 (没有pd信息这里设置T会报错)
                         filterMultiHit= T, # 移除对位点比对的探针
                         filterNoCG= T, # 仅保留CpG位点
                         filterSNPs= T, # 移除含有snp的位点
                         filterXY= T, # 移除性染色体上的探针
                         arraytype= "450K")  # 平台为450k
#好像大部分数据还是cpg的

### 挑选感兴趣区域（启动子CpG Island） ###

# 如果需要筛选一些特定区域的探针，可以加载如下注释文件并获得注释信息，并挑选感兴趣的探针
# 比如这里挑选启动子CpG Island探针
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
prompb <- rownames(anno[which(anno$Relation_to_Island == "Island" & grepl("TSS",anno$UCSC_RefGene_Group)),])  # 标记为启动子CpG Island探针

indata <- myFilter$beta[row.names(myFilter$beta) %in% prompb,]  #从所有探针中挑选出是启动子CpG Island探针

k <- 2500 # top k 个探针
var <- apply(indata, 1, sd)  
var <- sort(var,decreasing = T)
topcgimeth <- indata[names(var[1:k]),]
topcgimeth[1:3,1:3] #挑了SD前2500的位置

hcg <- hclust(distanceMatrix(t(topcgimeth), "euclidean"), "ward.D") #对探针进行聚类
# 可参照distanceMatrix的资料更换其他distance和linkage方法
gc() # 释放内存

hcs <- hclust(distanceMatrix(topcgimeth, "euclidean"), "ward.D")  #对样本进行聚类
group <- cutree(hcs,2) #此处分为3组
group <- paste0("C",as.character(group)); names(group) <- colnames(topcgimeth)

##############################################################
##绘图#

### 设置颜色 ###
blue   <- "#5bc0eb"
yellow <- "#fde74c"
green  <- "#9bc53d"
DarkTurquoise <- "#00CED1"
OrangeRed <- "#FF4500"

# 绘制热图
# 构建列注释信息
annotation_col = cbind(
  MethClust = data.frame(MethClust = group, row.names = names(group),stringsAsFactors = F), 
  RiskGroup = Sinfo$group
)
head(annotation_col)

annColors <- list(MethClust = c("C1" = blue,
                                "C2" = green),
                  RiskGroup = c("low" = DarkTurquoise,
                                "high" = OrangeRed))

pheatmap(topcgimeth,
         color = bluered(64),
         cluster_rows = hcg,
         cluster_cols = hcs,
         cutree_cols = 2,
         annotation_col = annotation_col,
         annotation_colors = annColors,
         show_rownames = F,show_colnames = F,
         filename = "methylation_cgi_clustering_pheatmap1.pdf")

print(table(annotation_col$MethClust,annotation_col$RiskGroup)) #table可以这样用的，是不是以前也遇见过

print(fisher.test(table(annotation_col$MethClust,annotation_col$RiskGroup))) 
#这里的意义是什么，甲基化的亚型和风险的分组是有相关性的？

ggplot(data=annotation_col,mapping=aes( x=RiskGroup,fill=MethClust))+
  geom_bar()+
  theme_bw()+
  scale_y_continuous(expand = c(0,0),limits=c(0,20))+ #(起始值，终止值，间隔),调整y轴属性，使柱子与X轴坐标接触,需要根据我们的样本量进行修改
  theme(
    #标题字体设置
    #text=element_text(size=12),
    #设置标题居中hjust = 0.5
    #plot.title = element_text(hjust = 0.5,vjust = 0.5), 
    #Y轴字体设置
    axis.text.y=element_text(size=12,color = "black"),
    #X轴字体设置
    #angle：调整横轴标签倾斜角度
    #hjust：上下移动横轴标签
    axis.text.x=element_text(size=12,  color = "black", hjust = 0.5,
                             vjust = 0.5),
    #图例的标题字体设置，可以修改colour、size
    #legend.title=element_text(size=12), 
    #图例字体设置
    legend.text=element_text(size=12)
    #legend.position = ' none' #删除图例
    #legend.position="bottom" ,#图例位置放中间，可选参数为“left”,“top”, “right”, “bottom”.
    #legend.background = element_rect(fill="lightblue",size=0.5, linetype="solid", colour ="darkblue")
  )

#这个可以用卡方

#########复杂热图
# 产生顶部注释
column_ha <- HeatmapAnnotation(df =  data.frame(cbind( MethClust = data.frame(MethClust = group, row.names = names(group),stringsAsFactors = F), 
                                                       RiskGroup = Sinfo$group)),
                               col = list(MethClust = c("C1" = blue,"C2" = green),
                                          RiskGroup = c("low" = DarkTurquoise,"high" = OrangeRed)),
                               RiskScore = anno_points(as.numeric(Sinfo$riskScore),
                                                       axis_param = 
                                                         list(side = "left",
                                                              at = c(0, 1, 2, 3),
                                                              labels = c("0", "1", "2", "3"))),
                               height = unit(2.7, "cm"),
                               show_annotation_name = T)
#这种顶部注释的写法可以参考一下，虽然感觉也不是很好看

hm <- Heatmap(topcgimeth, 
              name = "CpG methylation level",
              col = bluered(64), # 热图填充的颜色
              cluster_columns = hcs, # 列聚类
              cluster_rows = hcg, # 行聚类
              show_column_names = F, # 不显示列名
              show_row_names = F, # 不显示行名
              top_annotation = column_ha) # 顶部注释


pdf("methylation_cgi_clustering_complexheatmap2.pdf",width = 8,height = 7)
draw(hm)
invisible(dev.off()) #Invisible是什么，没太懂

###############################################
#差异甲基化，DMR还是什么来着

setwd("../02_DEmethylation/")
library(impute)

orgmeth <- fread("TCGA-STAD.methylation450_filter.tsv")

orgmeth <- as.data.frame(orgmeth); rownames(orgmeth) <- orgmeth[,1]; orgmeth <- orgmeth[,-1]
orgmeth = orgmeth[complete.cases(orgmeth),]#去除含有NA的数据
orgmeth[1:3,1:3]

Sinfo <- read.table("risk.train.group.txt", header = T, row.names = 1)
Sinfo <- Sinfo[row.names(Sinfo) %in% colnames(orgmeth), ]  # 去除没有甲基化的样本
table(Sinfo$group)


b <- data.frame(Group=rep(c("low","high"),c(17,19)))
row.names(b) <- row.names(Sinfo)

beta=as.matrix(orgmeth)

# 甲基化样本和分组样本信息一一对应，低风险在前，高风险在后
orgmeth <- orgmeth[,row.names(Sinfo)]
beta=as.matrix(orgmeth)

#knn补全
beta=impute.knn(beta)

betaData=beta$data
betaData=betaData+0.00001 #为了避免为0值吗？
a=betaData
a[1:4,1:4]

identical(colnames(a),rownames(b))

myLoad=champ.filter(beta = a,pd = b) #表达和分组,这过滤了啥，好像没过滤只是制备了champ需要的对象

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5) #标准化

dim(myLoad)
dim(myNorm) 

pD=myLoad$pd

#这里和教程得到的结果还不太一样，为什么教程是写的5核运行确实4核心，为什么会有过滤的
#以及'原来的450K经过质控过滤后是350K'是指什么？
# 450K是不是就是指有450k个探针？

group_list=myLoad$pd
table(group_list)

myDMP <- champ.DMP(beta = myNorm,pheno=group_list,adjPVal = 1)

head(myDMP[[1]])

DMP <- myDMP[[1]] #甚至注释都写得很清楚，有点厉害

#可视化

adj.P.Val <- 0.05
logFoldChange <- 0.2  
#注意，我们这里阈值经常设为0.2或0.1。在差异的信息里面要和deltaBeta进行比较，而不是deltaBeta。（在说什么呢？）

library(stringr)
DMP1 <- DMP[,c("deltaBeta","P.Value","adj.P.Val","AveExpr","low_AVG","high_AVG",
               "gene","CHR","MAPINFO","feature","UCSC_CpG_Islands_Name")]
DMP1 <- subset(DMP1,str_length(DMP1$gene)>0)  # 挑选有明确基因名的
DMP1$Gene_position <- paste0(DMP1$gene,":chr",DMP1$CHR,":",DMP1$MAPINFO,"(",DMP1$feature,")")  #使用位置信息表达甲基化位点

DMP_file<-DMP1
DMP_file$color <- ifelse(DMP_file$adj.P.Val<adj.P.Val & abs(DMP_file$deltaBeta)>= logFoldChange,
                         ifelse(DMP_file$deltaBeta> logFoldChange,'Up expression','Down expression'),'Non significant')

color <- c('Up expression' = "red",'Non significant' = "gray",'Down expression' = "blue")

# pdf(file="DMPs_volcano.pdf",width=10,height=8,onefile=FALSE)
p <- ggplot(DMP_file, aes(deltaBeta, -log10(adj.P.Val), col = color))+
  geom_point()+theme_bw()+scale_color_manual(values = color)+
  labs(title="",x="log2(fold change)",y="-log10(adj.P.Val)")+
  geom_hline(yintercept = -log10(adj.P.Val), lty=4,col="grey",lwd=0.6)+
  geom_vline(xintercept = c(-logFoldChange, logFoldChange), lty=4,col="grey",lwd=0.6)+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",
        panel.grid=element_blank(),axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

p<-p+theme_bw() +
  theme(plot.margin=unit(rep(3,4),'lines'),
        plot.title = element_text(size = 14, face =  "bold",hjust=0.5),
        text = element_text(size = 12),
        axis.title.x = element_text(face="bold",size=16),
        axis.title.y = element_text(face="bold",size=16),
        axis.text.x=element_text(size = 14,angle=0),
        axis.text.y=element_text(size = 14,angle=0),
        legend.text = element_text(face = "bold",size = 10,margin = margin(r=5)))
print(p)

###
# 设置阈值
adj.P.Val <- 0.05
logFoldChange <- 0.2  #注意，我们这里阈值经常设为0.2或0.1。在差异的信息里面要和deltaBeta进行比较，而不是deltaBet

ModelGene <- DMP1[DMP1$gene %in% c("DUSP1","MYB"),] #模型基因的甲基化差异信息

ModelGene = ModelGene[(ModelGene$adj.P.Val < adj.P.Val & (ModelGene$deltaBeta>logFoldChange | ModelGene$deltaBeta<(-logFoldChange))),]
#过滤掉之后就没有信息了，为什么会和教程的不一样呢？

red_de_expr <- myNorm[rownames(myNorm) %in% rownames(ModelGene),]  # 提取差异的模型基因甲基化位点
row.names(red_de_expr) <- ModelGene$Gene_position  #行名变成基因
write.table(red_de_expr,file="ModelGene_methylation.txt",sep="\t",quote=F)

table(Sinfo$group)

annotation_col=data.frame(clinical=rep(c("low(n=17)","high(n=19)"),c(17,19)))
annColors <- list(clinical = c("low(n=17)" = "#00A087FF","high(n=19)" ="#DC0000FF"))
rownames(annotation_col)=colnames(myNorm)
pheatmap(red_de_expr, 
         annotation_col = annotation_col,
         annotation_colors = annColors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cluster_rows=TRUE, 
         show_rownames=TRUE,
         fontsize_row=10,
         show_colnames=FALSE,
         scale="row",
         cluster_cols=FALSE,
         main="")

#果然还是数据的问题，两个模块里面用了不同的数据，并不是老演员
