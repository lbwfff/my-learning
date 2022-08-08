#######根据生信技能树推文进行修改的代码，表达矩阵数据来自10.15252/embj.201696056######
#######我不太适应这样的代码写法，后续可能会把代码改成自己习惯的样子

rm(list = ls())
library('WGCNA')

######1，数据处理######
norm.counts <- read.csv("exp-for-WGCNA.csv",header = T,sep = ",")
rownames(norm.counts)<-norm.counts[,1]
norm.counts <- norm.counts[,2:18]

a <- colnames(norm.counts)
library(stringr)
str_tmp=as.data.frame(str_split(a,"[-]",simplify = T))
rownames(str_tmp) <- colnames(norm.counts)

str_tmp$group<-c(rep('CTR.adult.MG',3),rep('NEO.MG',4),rep('CD11c.NEO.MG',4),rep('EAE.MG',3),rep('CD11c.EAE.MG',3))
str_tmp$sampleRep<-c('1','2','3','1','2','3','4','1','2','3','4','1','2','3','1','2','3')
str_tmp<-str_tmp[,-1]

datTraits <- str_tmp
datTraits$groupNo <- c(rep(1,3),rep(2,4),rep(3,4),rep(4,3),rep(5,3))

norm.counts <- log2(norm.counts)  ##纠结是否应该进行归一化的话，可以都做一下看看怎样聚类更好
m.mad <- apply(norm.counts,1,mad)  ##mad,中位数绝对偏差
dataExprMad <- norm.counts[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.3)),] #筛选中位绝对偏差前75%的基因，当然也可不做筛选
datExpr0 <- as.data.frame(t(dataExprMad))

gsg <- goodSamplesGenes(datExpr0,verbose = 3) #####WGCNA包带的函数，可以检查缺失条目什么的
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

if(T){
  # Plot a line to show the cut
  sampleTree <- hclust(dist(datExpr0),method = "average")
  par(cex=.6)
  par(mar=c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub = "",
       xlab = "",cex.lab = 1.5, cex.axis =1.5, cex.main=2) ######树图看看聚类情况
}

dim(datExpr0) ##17，8222
datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="input.Rdata")

########计算beta值（不懂）#######
rm(list = ls())
load(file = "input.Rdata")
dim(datExpr)

if(T){
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  pdf("step2-beta-value.pdf",width = 8,height = 6)
  # Plot the results:
  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.9,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}
sft$powerEstimate  ## beta=22 SCI文章里面用了20
save(sft,file = "step2_beta_value.Rdata")

#########构建加权共表达网络######
rm(list = ls())
library(WGCNA)
load(file = "input.Rdata")
load(file = "step2_beta_value.Rdata")
enableWGCNAThreads()
if(T){
  net = blockwiseModules(
    datExpr,
    power = sft$powerEstimate,
    maxBlockSize = nGenes,
    TOMType = "unsigned", minModuleSize = 30, ##minModuleSize，表示每个模块里面最少放多少个基因，设定越大，模块越少
    reassignThreshold = 0, mergeCutHeight = 0.28, ##mergeCutHeight参数表示你在哪里砍树。值设定越小，树枝越多，通常是0.25
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = F,
    verbose = 3
  )
  table(net$colors) 
}

##模块可视化
if(T){
  # Convert labels to colors for plotting
  moduleColors=labels2colors(net$colors)
  table(moduleColors)
  # Plot the dendrogram and the module colors underneath
  pdf("step3-genes-modules.pdf",width = 8,height = 6)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}
save(net,moduleColors,file = "step3_genes_modules.Rdata")

#########对应模块和基因表型#####
rm(list = ls())
library(WGCNA)
load(file = "input.Rdata")
load(file = "step2_beta_value.Rdata")
load(file = "step3_genes_modules.Rdata")

if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  design <- model.matrix(~0+datTraits$group)
  datTraits$group=as.factor(datTraits$group)
  colnames(design)= levels(datTraits$group) ## get the group
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes
  MEs = orderMEs(MES0)
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                     signif(moduleTraitPvalue,1),")",sep = "")
  dim(textMatrix)=dim(moduleTraitCor)
  
  pdf("step4-Module-trait-relationship.pdf",width = 8,height = 12)
  par(mar=c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  save(design,file = "step4_design.Rdata")
}

#把上面的模块图转化成了bar图
if(T){
  mes_group <- merge(MEs,datTraits,by="row.names") # 小技巧，可以通过rowname进行merge
  ##写了个画图的function
  library(gplots)
  library(ggplot2)
  library(ggpubr)
  library(grid)
  library(gridExtra) 
  draw_ggboxplot <- function(data,gene="P53",group="group"){
    #print(gene)
    ggboxplot(data,x=group, y=gene,
              ylab = sprintf("Expression of %s",gene),
              xlab = group,
              fill = group,
              palette = "jco",
              #add="jitter",
              legend = "") +stat_compare_means()
  }
  ###开始批量画boxplot
  colorNames = names(MEs)
  pdf("step4-expression-group.pdf",width = 14,height=28)
  #par(mfrow=c(ceiling(length(colorNames)/2),2))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group,gene= x,group="group")
  })
  do.call(grid.arrange,c(p,ncol=2))
  dev.off()
}

###模块与基因的相关性
if(T){
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  ## 算出每个模块跟基因的皮尔森相关系数矩阵
  ## MEs是每个模块在每个样本里面的值
  ## datExpr是每个基因在每个样本的表达量
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  
  geneTraitSignificance <- as.data.frame(cor(datExpr,datTraits$groupNo,use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
  names(geneTraitSignificance)<- paste("GS.",names(datTraits$group),sep = "")
  names(GSPvalue)<-paste("GS.",names(datTraits$group),sep = "")
  
  #selectModule<-c("blue","green","purple","grey")  ##可以选择自己喜欢的模块
  selectModule <- modNames  ## 批量作图
  pdf("step4-Module-trait-significance.pdf",width = 14,height=28)
  par(mfrow=c(ceiling(length(selectModule)/2),2)) #批量作图开始
  for(module in selectModule){
    column <- match(module,selectModule)
    print(module)
    moduleGenes <- moduleColors==module
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", module, "module"),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
  dev.off()
}


#############关注感兴趣的模块
table(moduleColors)
group_g <- data.frame(gene=colnames(datExpr),
                      group=moduleColors)
write.csv(group_g,file = "group_g.csv",row.names = F) ## 导出的就是模块和基因的对应关系，可以做go分析什么的，不过这次就免了吧

#############热图
if(T){
  #geneTree = net$dendrograms[[1]]
  TOM=TOMsimilarityFromExpr(datExpr,power=20)
  dissTOM=1-TOM
  #plotTOM = dissTOM^7
  #diag(plotTOM)=NA
  #TOMplot(plotTOM,geneTree,moduleColors,main="Network heapmap plot of all genes")
  nSelect =300 ####选取的基因数，数量越多需要运算的时间就越久，不过也会好看一些
  set.seed(20)
  select=sample(nGenes,size = nSelect)
  selectTOM = dissTOM[select,select]
  selectTree = hclust(as.dist(selectTOM),method = "average")
  selectColors = moduleColors[select]
  plotDiss=selectTOM^7
  diag(plotDiss)=NA
  pdf("step6_select_Network-heatmap.pdf",width = 8,height=6)
  TOMplot(plotDiss,selectTree,selectColors,main="Network heapmap of selected gene")
  dev.off()
}

#####指定模块基因的热图，这串代码会让我的rstudio卡死不知道为什么？
if(T){
  module="turquoise"
  which.module=module
  dat=datExpr[,moduleColors==which.module]
  library(pheatmap)
  pheatmap(dat,show_colnames = F,show_rownames = F)
  n=scale(t(dat+1)) # 'scale'可以对log-ratio数值进行归一化
  n[n>2]=2 
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  group_list=datTraits$group
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  png("figures/step6-moduleGene-heatmap.png",width = 800,height = 600)
  pheatmap(n,show_colnames =F,show_rownames = F,annotation_col =ac )
  dev.off()
}





###############没完全懂，但是完全会了，以下为非常重要的两段###########
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs #想要做相关性的话，可以用这个
geneTree = net$dendrograms[[1]]

####1####
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p") #在这里对datTraits进行修改，可以引入更多的元素（连续型）
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

####2####
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

design <- model.matrix(~0+datTraits$group)
datTraits$group=as.factor(datTraits$group)
colnames(design)= levels(datTraits$group) #这里的design可以修改为任何元素（离散型），虽然我不大理解这里使用离散型变量的话R值和P值都代表什么？
MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs = orderMEs(MES0)
moduleTraitCor <- cor(MEs,design,use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)





