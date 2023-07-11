#一个做药物预测分析的包，有GDSC（癌药物敏感性基因组学数据库）和CTRP（也是个癌症药敏的数据库）的数据，甚至于可以自己去设置训练集，非常自由的一个框架
#用的应该是类似于岭回归之类的算法，模型的表现可能不会特别好（相较之下），但是考虑到时间成本，也还不错吧

setwd("01_GDSC/")

library(oncoPredict)

library(ggplot2)
library(cowplot)

options(stringsAsFactors = FALSE)

dat <- read.table("STAD_tumor_orderByRiskGroup.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
dat[1:3, 1:3]

ann <- read.table("risk.train.group.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
ann <- ann[order(ann$group,decreasing = T),]
head(ann)

table(ann$group)

#表达数据、分组信息的ID一致
dat <- dat[,row.names(ann)]

#
dir='../DataFiles/DataFiles/Training Data/'
#oncopredict提供的训练集(我不确定这个能不能叫训练集)，有GDCS和CTRP两个数据库的数据集

GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))  #什么叫RMA？
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

#药物名字
GCP.drug <- colnames(GDSC2_Res) #
#这里以前12种药物为例
GCP.drug <- GCP.drug[1:12]
GCP.drug

# 自定义足够多的box的颜色，颜色数量至少等于分组数量
jco <- c("#EABF00", "#2874C5", "red")

### 药物敏感性预测 ###、
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表
plotp <- list()

set.seed(1248103) # 因为预测过程默认10-fold CV，所以设置种子以便结果可重复
cat(drug," starts!\n") # 提示当前药物已开始分析
  
# 预测IC50值，采用默认参数，详细可参考??pRRopheticPredict参数列表
martix<-as.matrix(dat[,rownames(ann)])
  
oncoPredict::calcPhenotype(trainingExprData = GDSC2_Expr,
                           trainingPtype = GDSC2_Res[,c(GCP.drug)],
                           testExprData = testExpr, #目标数据集
                           batchCorrect = 'eb',  #   "eb" for ComBat  
                           powerTransformPhenotype = TRUE,
                           removeLowVaryingGenes = 0.2,
                           minNumSamples = 10, 
                           printOutput = TRUE, 
                           removeLowVaringGenesFrom = 'rawData' )

#没必要非得跟着他的代码来，这种自己写一个循环吧箱线画出来就好了
