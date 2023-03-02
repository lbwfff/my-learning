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


###################################################################################################################
#一个自己的实例，有空再写注释吧，如果我还记得自己写了些什么的话

load('~/biodata/CNHPP_HCC/MS/exp.RData')

library(limma)
exp_n <- normalizeBetweenArrays(exp[,3:194], method="quantile")
exp_n<-as.data.frame(exp_n)
exp<-cbind(exp[,1:2],exp_n)

rownames(exp)<-exp$Majority.protein.IDs
exp<-exp[,-(1:2)]

group<-c('Paratumor','Tumor')
group_list<-rep(group,96)

library('WGCNA')

datTraits<-as.data.frame(group_list)
colnames(datTraits)<-c('group')
rownames(datTraits)<-colnames(exp)
datTraits$sampleRep<-c('whatever')

norm.counts<-exp
m.mad <- apply(norm.counts,1,mad)  ##mad,中位数绝对偏差
dataExprMad <- norm.counts[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),] #25%  = 0.22

datExpr0 <- as.data.frame(t(dataExprMad))

gsg <- goodSamplesGenes(datExpr0,verbose = 3)
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
       xlab = "",cex.lab = 1.5, cex.axis =1.5, cex.main=2
  )
}

dim(datExpr0) 
datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

if(T){
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, 
                          verbose = 5,networkType = "signed hybrid")
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
save(sft,file = "~/biodata/CNHPP_HCC/MS/step2_beta_value.Rdata")
save(datExpr,sft,nGenes,nSamples,file = '~/biodata/CNHPP_HCC/MS/WGCNA.RData')#放服务器吧老是崩掉

enableWGCNAThreads()
if(T){
  net = blockwiseModules(
    datExpr,networkType='signed hybrid',
    power = sft$powerEstimate,
    maxBlockSize =nGenes , #nGenes
    TOMType = "signed", minModuleSize = 40,
    reassignThreshold = 0, mergeCutHeight = 0.2, ## 这个值越小，模块数越多
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = F,deepSplit = 4,
    verbose = 3
  )
  table(net$colors) #我居然觉得好像还不错
}

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
save(net,moduleColors,file = "~/biodata/CNHPP_HCC/MS/step3_genes_modules.Rdata")

load('~/biodata/CNHPP_HCC/MS/step3_genes_modules.Rdata')

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
  selectModule <- modNames  ## 也可以批量作图
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

table(moduleColors)
group_g <- data.frame(gene=colnames(datExpr),group=moduleColors) #用这个来做富集

#如果我只对和上下调lncPEP相关的模块进行后续的分析的话
group_g<-group_g[paste0('ME',group_g$group) %in% unique(sanky$module),]

get_id<-function(module){
  
f <-function(x) unlist(strsplit(x['gene'],'[|]'))[2]
grey<-group_g[group_g$group==module,]
grey<- apply(grey,1,f)
grey<- bitr(grey, fromType = "UNIPROT",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)[2]

return(grey)
}

list<-list()
for (i in 1:length(unique(group_g$group))){
  module<-unique(group_g$group)[i]
  list[i]<-get_id(module)
  names(list)[[i]]<-c(module)
}

for (i in 1:17){
  p<-list()
  enrich<-unique(list[[i]])
  kk.negative <- enrichKEGG(gene  = enrich,
                            organism = "hsa",
                            #universe = gene_all,
                            pvalueCutoff = 0.5,
                            qvalueCutoff =0.5)
  kk.negative<-setReadable(kk.negative,OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  head(kk.negative)[,1:6]
  kk=kk.negative@result
  kk<-kk[order(kk$pvalue,decreasing = F),]
  
  library('stringr')
  title<-paste0('KEGG')
  color<-c(met.brewer('Hokusai1')[c(3,5,6,7)],met.brewer('Troy')[5])
  
  p[[1]]<-
    ggplot(data = kk[1:10,], aes(x = reorder(Description,-pvalue), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(title) +coord_flip() +
    theme_classic(base_size = 13)+ 
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(P value)",x=' ')+
    scale_fill_manual(values = c(color[1]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))+
    theme(aspect.ratio=1.4)
  
  BP <- enrichGO(gene          = enrich,
                 OrgDb         = org.Hs.eg.db,
                 ont           = 'BP' ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.5,
                 qvalueCutoff  = 0.5,
                 readable      = TRUE)
  # BP <- simplify(BP,cutoff = 0.7,by = "p.adjust",
  #                select_fun = min,measure = "Wang",semData = NULL)
  
  CC <- enrichGO(gene          = enrich,
                 OrgDb         = org.Hs.eg.db,
                 ont           = 'CC' ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.5,
                 qvalueCutoff  = 0.5,
                 readable      = TRUE)
  # CC <- simplify(CC,cutoff = 0.7,by = "p.adjust",
  #                select_fun = min,measure = "Wang",semData = NULL)
  
  MF <- enrichGO(gene          = enrich,
                 OrgDb         = org.Hs.eg.db,
                 ont           = 'MF' ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.5,
                 qvalueCutoff  = 0.5,
                 readable      = TRUE)
  # MF <- simplify(MF,cutoff = 0.7,by = "p.adjust",
  #                select_fun = min,measure = "Wang",semData = NULL)
  
  DO<-enrichDO(gene          = enrich,
               ont           = "DO",
               pvalueCutoff  = 0.2,
               pAdjustMethod = "BH",
               #universe      = names(geneList),
               minGSSize     = 5,
               maxGSSize     = 500,
               qvalueCutoff  = 0.2,
               readable      = FALSE)
  
  BP=BP@result
  BP<-BP[order(BP$pvalue,decreasing = F),]
  CC=CC@result
  CC<-CC[order(CC$pvalue,decreasing = F),]
  MF=MF@result
  MF<-MF[order(MF$pvalue,decreasing = F),]
  DO=DO@result
  DO<-DO[order(DO$pvalue,decreasing = F),]
  
  title<-paste0('GO_BP')
  p[[2]]<-
    ggplot(data = BP[1:10,], aes(x = reorder(Description,-pvalue), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(title) +coord_flip() +
    theme_classic(base_size = 13)+ 
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(P value)",x=' ')+
    scale_fill_manual(values = c(color[2]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))+
    theme(aspect.ratio=1.4)
  
  title<-paste0('GO_CC')
  p[[3]]<-
    ggplot(data = CC[1:10,], aes(x = reorder(Description,-pvalue), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(title) +coord_flip() +
    theme_classic(base_size = 13)+ 
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(P value)",x=' ')+
    scale_fill_manual(values = c(color[3]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))+
    theme(aspect.ratio=1.4)
  
  title<-paste0('GO_MF')
  p[[4]]<-
    ggplot(data = MF[1:10,], aes(x = reorder(Description,-pvalue), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(title) +coord_flip() +
    theme_classic(base_size = 13)+ 
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(P value)",x=' ')+
    scale_fill_manual(values = c(color[4]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))+
    theme(aspect.ratio=1.4)
  
  title<-paste0('DO') #names(list)[i],
  p[[5]]<-
    ggplot(data = DO[1:10,], aes(x = reorder(Description,-pvalue), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(title) +coord_flip() +
    theme_classic(base_size = 13)+
    theme(panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(P value)",x='Term')+
    scale_fill_manual(values = c(color[5]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))+
    theme(aspect.ratio=1.4)
  
  path<-paste0('/home/leelee/biodata/CNHPP_HCC/MS/WGCNA/',names(list)[i],'_enrich.pdf')
  pdf(path,width = 10,height = 8)
  print(wrap_plots(p,nrow=3)) 
  dev.off()
  
  path<-paste0('/home/leelee/biodata/CNHPP_HCC/MS/WGCNA/',names(list)[i],'_enrich.xlsx')
  library(openxlsx)
  sheets = list("KEGG" = kk,"GO_BP" = BP,'GO_CC'=CC,'GO_MF'=MF,'DO'=DO) #
  write.xlsx(sheets,path)
  
} #这只是段富集画图的代码

split = factor(group_g$group,levels = c(unique(group_g$group)))

term<-list()

for (i in 1:length(unique(group_g$group))){
  path<-paste0('/home/leelee/biodata/CNHPP_HCC/MS/WGCNA/',unique(group_g$group)[i],'_enrich.xlsx')
  bp<-read.xlsx(path, sheet = 2)
  bp<-bp[bp$pvalue<0.01,]
  
  # bp<-unlist(strsplit(bp$Description,' '))
  # bp<-bp[bp!='regulation' & bp!='of'& bp!='by'& bp!='to'& bp!='the']
  # term[i]<-paste0(bp,collapse = ' ')
  
  # term[i]<-paste0(bp$Description,collapse = ';')
  term[[i]]<-c(bp$ID)
  
  names(term)[i]<-unique(group_g$group)[i]
} #这里把富集的结果提出来了

library(ComplexHeatmap)
library(simplifyEnrichment)

mat<-dataExprMad[match(group_g$gene,rownames(dataExprMad)),]
mat<-log2(mat+1)
mat = t(scale(t(mat)))

split_top <- factor(c(rep(c(1,2),96)),levels=c(1,2))
ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:2),labels = c("Paired non-tumour", "Tumour"), 
                                         labels_gp = gpar(col = "white", fontsize = 10)))

col_fun1 <- colorRamp2(c(-1,0,1), c(met.brewer("Homer1")[7],'white',met.brewer("Signac")[4]))
# col_fun1 <- colorRamp2(c(-1,0,1), c('#16317d','white','#ffcd12')) #这种大热图用这个配色就变得很丑
 
left <- rowAnnotation(foo = anno_block(gp = gpar(fill = unique(group_g$group))))  #这个颜色顺序绝对搞错了,想不到要怎么做，去问问作者吧

pdf('~/biodata/CNHPP_HCC/MS/WGCNA/words_pheat.pdf',width = 10,height = 8)
Heatmap(mat, row_split = split, show_row_names = F,show_column_names = F,show_column_dend = F,
        width = unit(10, "cm"),height = unit(14, "cm"),cluster_rows = T,cluster_column = T,
        col = col_fun1,row_title = NULL, name = "z-score",cluster_row_slices = FALSE,cluster_column_slices = FALSE,
        show_row_dend = FALSE,top_annotation = ha,column_split = split_top,column_title = NULL,column_gap = unit(c(1, 1), "mm"),
        right_annotation = rowAnnotation(wc = anno_word_cloud_from_GO(split, term,
                                                              exclude_words=c('regulation','cell','import','activity','response',
                                                                              'process','positive','negative','pathway'))),
        left_annotation = left) 
dev.off() 

#这个热图大致我是比较满意的，但我想让热图的模块顺序和桑基图一致，这在热图这里调整有点困难，可以在画桑基的时候做调整

MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes #用这个来做相关性的分析

nmfprot<-log2(orf_intensity+1)
mor0 <- apply(nmfprot,1,function(x) sum(x>0) )
# nmfprot <- nmfprot[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.2)),] #25%  = 0.22
nmfprot <- nmfprot[which(mor0 > 9),] #这也太少了吧

cortest<-as.data.frame(array(NA,c(17,97)))
cortest_p<-as.data.frame(array(NA,c(17,97)))
colnames(cortest)<-rownames(nmfprot)
colnames(cortest_p)<-rownames(nmfprot)
rownames(cortest)<-colnames(MES0)
rownames(cortest_p)<-colnames(MES0)

for (i in 1:nrow(cortest)){
  module<-rownames(cortest)[i]
  for (j in 1:ncol(cortest)){
    orf<-colnames(cortest)[j]
    
    x<-MES0[,colnames(MES0)==module]
    y<-nmfprot[rownames(nmfprot)==orf,]
    x<-x[which(y>0)]
    y<-y[which(y>0)]
    
    cor<-cor.test(x,y)
    cortest[i,j]<-cor$estimate
    cortest_p[i,j]<-cor$p.value
  }
}

orf_module<-as.data.frame(array(NA,c(97,4)))
colnames(orf_module)<-c('name','module','R','P')
orf_module$name<-colnames(cortest)
for (i in 1:nrow(orf_module)){
  orf<-orf_module$name[i]
  cor<-cbind(cortest[,colnames(cortest)==orf],
             cortest_p[,colnames(cortest_p)==orf] )
  cor<-as.data.frame(cor)
  cor$module<-rownames(cortest)
  
  cor<-cor[order(cor$V2),]
  cor<-cor[abs(cor$V1)>0.3 & cor$V2<(0.01),]
  
  if(nrow(cor)>4) { cor<-cor[1:4,]
  } else{cor<-cor }
  
  orf_module$module[i]<-paste0(cor$module,collapse = ';')
  orf_module$R[i]<-paste0(cor$V1,collapse = ';')
  orf_module$P[i]<-paste0(cor$V2,collapse = ';')
  
  }

orf_module<-orf_module[nchar(orf_module$module)>0,]

# deg<-ms_orf_table[ms_orf_table$Adj_ORFID %in% orf_module$name,]
# deg<-deg[match(orf_module$name,deg$Adj_ORFID),]
# orf_module$deg<-deg$DEG_group
# orf_module<-orf_module[!is.na(orf_module$deg),]

deg<-ms_orf_table[!is.na(ms_orf_table$DEG_group),]

up<-deg[deg$DEG_P<0.01 & deg$DEG_logFC>0 ,]
down<-deg[deg$DEG_P<0.01 & deg$DEG_logFC<0 ,]
nosig<-deg[deg$DEG_P>0.01 ,]

deg<-rbind(up,down,nosig)
deg$deg<-c(rep('Up',nrow(up)),rep('Down',nrow(down)),rep('Unsignificance',nrow(nosig)))

cache<-as.data.frame(array(NA,c(1,3)))
colnames(cache)<-c('pep','module','group')

for (i in 39:nrow(orf_module)){
  orf<-orf_module$name[i]
  module<-unlist(strsplit(orf_module$module[i],';'))
  sanky<-data.frame(pep=c(orf),
                    module=c(module),
                    group=c(deg$deg[which(deg$Adj_ORF_ID==orf)]))  
  
  cache<-rbind(cache,sanky)
}

sanky<-cache[-1,]
sanky<-sanky[sanky$group!='Unsignificance',]

sanky$sanky<-paste0(sanky$group,'_',sanky$module)

sankey_2<-as.data.frame(table(sanky$sanky))

f <-function(x) unlist(strsplit(x['Var1'],'_'))[1]
sankey_2$DEG <-apply(sankey_2,1,f)
f <-function(x) unlist(strsplit(x['Var1'],'_'))[2]
sankey_2$MEs <-apply(sankey_2,1,f)
sankey_2<-sankey_2[,-1]
colnames(sankey_2)<-c('weight','from','to')

node<-as.data.frame(unique(c(sankey_2$from,sankey_2$to)))
node$lable<-node[,1]
node[,1]<-c(1:nrow(node)-1)
colnames(node)[1]<-c('id')
node$lable[3:12]<-c(paste0('ME',unique(group_g$group)))

cache<-node[match(sankey_2$from,node$lable),]
sankey_2$from<-cache$id
cache<-node[match(sankey_2$to,node$lable),]
sankey_2$to<-cache$id
sankey_2<-sankey_2[,c(2,3,1)]

library(networkD3)
color <- 'd3.scaleOrdinal() .domain(["Down","Unsignificance","Up","MEblack", "MEblue", "MEcyan", "MEgreenyellow", "MEgrey", "MElightcyan", "MEpink","MEpurple","MEturquoise","MEyellow"]) .range(["#447fdd","#f6f2ee","#da6c42","black", "blue", "cyan", "greenyellow", "grey", "lightcyan", "pink","purple","turquoise","yellow"])'

sank <-sankeyNetwork(
  Links = sankey_2, Nodes = node, 
  Source = "from", Target = "to", 
  NodeID = "lable", Value = "weight", colourScale = color,
  fontSize = 16,width=600)
saveNetwork(sank, "test.html")

#性状的相关性
hccmeta<-read.csv('~/biodata/CNHPP_HCC/MS/WGCNA/HCC_meta.csv')
hccmeta$adj_sample<-paste0('iBAQ.',hccmeta$Case.No.,'_T')
#性别，年龄，肝硬化，肿瘤数，Microvascular.invasion,AFP
cormeta<-hccmeta[,c(2,3,6,8,12,13,23)]
cormeta$Gender<-ifelse(cormeta$Gender=='M','0','1')
colnames(cormeta)[1]<-c('Gender_Female')
colnames(cormeta)[6]<-c('AFP')

MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs = orderMEs(MES0)

cormeta<-cormeta[cormeta$adj_sample %in% rownames(MEs),]
MEs_4<-MEs[match(cormeta$adj_sample,rownames(MEs)),]

moduleTraitCor <- cor(MEs,data.frame(Paratumor=rep(c(1,0),96),Tumor=rep(c(0,1),96)),use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,192)

moduleTraitCor2 <- cor(MEs_4,cormeta[1:6],use = "p")
moduleTraitPvalue2 <- corPvalueStudent(moduleTraitCor2,96)

moduleTraitCor<-cbind(moduleTraitCor,moduleTraitCor2,moduleTraitCor3[,25])
moduleTraitPvalue<-cbind(moduleTraitPvalue,moduleTraitPvalue2,moduleTraitPvalue3[,25])

textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                   signif(moduleTraitPvalue,1),")",sep = "")
dim(textMatrix)=dim(moduleTraitCor)

pdf("~/biodata/CNHPP_HCC/MS/WGCNA/Module_trait_corpheat.pdf",width = 8,height = 12)
par(mar=c(8, 9, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c('Paratumor','Tumor',colnames(cormeta[1:6]),'InfiltrationScore'),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               plotLegend = T,
               colors = rev(colors),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships")) #颜色可以后续再考虑
dev.off()

#调整一下颜色
color<-c(met.brewer("Signac")[4],met.brewer("Signac")[4],'white','white','white',met.brewer("Homer1")[7],met.brewer("Homer1")[7]) #从包里面取了四个色
pal<-colorRampPalette(color)
image(x=1:100,y=1,z=as.matrix(1:100),col=pal(100))
colors = pal(100)

#然后再画一张免疫浸润相关的
identical(rownames(MEs),ImmuCellAIscore$sample)
moduleTraitCor3 <- cor(MEs,ImmuCellAIscore[,1:25],use = "p")
moduleTraitPvalue3 <- corPvalueStudent(moduleTraitCor3,192)

textMatrix = paste(signif(moduleTraitCor3,2),"\n(",
                   signif(moduleTraitPvalue3,1),")",sep = "")
dim(textMatrix)=dim(moduleTraitCor3)

pdf("~/biodata/CNHPP_HCC/MS/WGCNA/Module_ImmuCellAIscore_corpheat.pdf",width = 10,height = 12)
par(mar=c(8, 9, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor3,
               xLabels = c(colnames(moduleTraitCor3)),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               plotLegend = T,
               colors = blueWhiteRed(50),
               # textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships")) #颜色可以后续再考虑
dev.off()

#画一张只有部分模块的

modlist<-c('MEgrey','MEturquoise','MEpurple','MEblue','MEpink','MEyellow',
           'MEblack','MEgreenyellow','MElightcyan','MEcyan')
moduleTraitCor2<-moduleTraitCor[match(modlist,rownames(moduleTraitCor)),]
moduleTraitPvalue2<-moduleTraitPvalue[match(modlist,rownames(moduleTraitPvalue)),]
textMatrix2 = paste(signif(moduleTraitCor2,2),"\n(",
                   signif(moduleTraitPvalue2,1),")",sep = "")
dim(textMatrix2)=dim(moduleTraitCor2)

pdf("~/biodata/CNHPP_HCC/MS/WGCNA/Module_trait_corpheat_part.pdf",width = 8,height = 12)
par(mar=c(8, 9, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor2,
               xLabels = c('Paratumor','Tumor',colnames(cormeta[1:6]),'Immune infiltration'),
               yLabels = modlist,
               ySymbols = modlist,
               colorLabels = F,
               plotLegend = T,
               colors = rev(colors),
               textMatrix = textMatrix2,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships")) #颜色可以后续再考虑
dev.off()


