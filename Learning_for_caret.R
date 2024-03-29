###############Boruta,一个特征选择包，用起来还比较方便，当然特征选择肯定有很多的方法，可以多了解一下##########
library(Boruta)
set.seed(20230322)#可以先设置种子
Var.Selec<-Boruta(label~., data=assaylearn,maxRuns=500, doTrace=1)

Var.Selec$finalDecision

pdf('Importance.pdf',width = 18,height = 8,family="GB1") #
par(mar=c(20,5,2,5)) #下方，左方，上方，右方
plot(Var.Selec,xlab="Attributes",ylab="Importance:Z-Score",las=2) #这玩意的美化还挺麻烦的，但是也不是就没有思路吧
dev.off()  

plotImpHistory(Var.Selec)

test<-attStats(Var.Selec)

getSelectedAttributes(Var.Selec,withTentative=FALSE)

#对应Boruta的绘图我做了修改
Borutaplot<-function(vc){
  imparray<-vc[["ImpHistory"]]
  imparray<-reshape2::melt(imparray)
  imparray$value[imparray$value == -Inf] <- 0
  
  feasel<-as.data.frame(vc[["finalDecision"]])
  feasel<-feasel[match(imparray$Var2,rownames(feasel)),]
  imparray$group<-as.character(feasel)
  imparray$group[is.na(imparray$group)]<-c('shadow')
  imparray$group<-factor(imparray$group,levels = c('Confirmed','Rejected','shadow','Tentative'))
  library(ggplot2)
  
  p<-ggplot(imparray, aes(x = reorder(Var2,value), y = value)) +
    stat_boxplot(geom = "errorbar",width=0.2)+
    geom_boxplot(outlier.shape = 1,aes(fill=group),show.legend = F)+
    # geom_boxplot(width=0.8,aes(fill=group),colour='black',alpha = 1,outlier.shape = NA)+
    labs(y="Importance:Z-Score")+ 
    theme_classic(base_size = 22)+ 
    scale_fill_manual(values=MetBrewer::met.brewer("Hokusai1", 7)[c(6,3,7,4)])+
    theme_classic(base_size = 14)+ 
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    coord_flip()+
    xlab(NULL)+
    theme(aspect.ratio=2)+
    theme(axis.title.x = element_text(size = 16, color = "black"),
          axis.title.y = element_text(size = 16, color = "black"),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 14, color = "black"),
          legend.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(fill = NA, color = "black"),
          panel.grid.minor = element_line(color = "gray", size = 0.25))
  
  return(p)
  }

cairo_pdf("Importance.pdf",width = 10,height = 12)
Borutaplot(Var.Selec)
dev.off()

#############记录一下学习caret包的过程吧，想用机器学习来做当前project的收尾，不过暂时自己能做的还比较少，总之一边学习一边思考吧########
#############参考自https://topepo.github.io/caret/index.html################

library('caret')
library(doParallel)
cl <- makePSOCKcluster(16)
registerDoParallel(cl) #并行运算非常重要，否则真是等到天昏地暗
# stopCluster(cl) 

##############可视化########################
#他这里介绍的可视化是feature的可视化，我们针对转录组蛋白组数据做的话可能不大有意义，所以直接跳过了

##############前处理########################
#第一条是创建虚拟变量
data(etitanic)
head(model.matrix(survived ~ ., data = etitanic))#虽然不属于caret包，但是model.matrix这个函数非常有意思
dummies <- dummyVars(survived ~ ., data = etitanic)#和model.matrix结果不同的在性别这一列，除此之外在colname的处理上更好一些？survived因为是想要预测的结果所以没有做处理吗？不大明白

#去除0和0变异
nzv <- nearZeroVar(mdrrDescr)#默认saveMetrics= F
filteredDescr <- mdrrDescr[, -nzv]#去除了矩阵里的零变异变量
dim(filteredDescr)

#降低变量之间的相关性
descrCor <- cor(filteredDescr)
summary(descrCor[upper.tri(descrCor)])
highlyCorDescr <- findCorrelation(descrCor, cutoff = .75)
filteredDescr <- filteredDescr[,-highlyCorDescr] #去除了高相关的，为什么要这么做呢？
descrCor2 <- cor(filteredDescr)
summary(descrCor2[upper.tri(descrCor2)])

#去除线性依赖
comboInfo <- findLinearCombos(ltfrDesign)
ltfrDesign[, -comboInfo$remove]

#中心化和缩放，补全
preProcValues <- preProcess(training, method = c("center", "scale"))#这个我熟

#特征变量转化
transformed <- spatialSign(plotSubset)
transformed <- as.data.frame(transformed)

preProcValues2 <- preProcess(training, method = "BoxCox")#给了两种方法，完全没懂在做什么

#All Together
pp_no_nzv <- preProcess(schedulingData[, -8], method = c("center", "scale", "YeoJohnson", "nzv"))#往method里加就可以了

#Class Distance Calculations#这个我完全没懂就没放上来了，之后理解了再说吧


################数据分离#######################

#1，基于结果，如下基于Species结果做的分离
trainIndex <- createDataPartition(iris$Species, p = .8, list = FALSE, times = 1)
irisTrain <- iris[ trainIndex,]
irisTest  <- iris[-trainIndex,]

#2，基于预测
startSet <- sample(1:dim(testing)[1], 5) #随了5行
samplePool <- testing[-startSet,]
start <- testing[startSet,]
newSamp <- maxDissim(start, samplePool, n = 20)#maximum dissimilarity approach，最大不相似尝试？什么玩意

#2，时间序列

#4，根据重要组（？）

#################模型训练和调整##################

#用sonar数据集做了一个示例

set.seed(998)
inTraining <- createDataPartition(Sonar$Class, p = .75, list = FALSE)#基于结果的分离，算是目前最好理解的一种分离方式吧，0.75的数据作为了训练集
training <- Sonar[ inTraining,]
testing  <- Sonar[-inTraining,]

#trainControl用来指定resampling的方法，这里的方法特别重要，之后会详细说明
fitControl <- trainControl(method = "repeatedcv", ## 10-fold CV，repeatedcv是K-fold交叉验证重复多次，别问我k-fold是什么
                           number = 10,                           
                           repeats = 10) ## repeated ten times

#train
gbmFit1 <- train(Class ~ ., data = training, #class为训练的预测目标，
                 method = "gbm", 
                 trControl = fitControl,
                 verbose = FALSE)  ## This last option is actually one for gbm() that passes through，这个参数好像是gbm模型特有的？
  
gbmFit1 #这个模型需要调的参数有number of iterations，complexity of the tree，shrinkage，n.minobsinnode，其中后两个被固定在了0.1和10（不知道处于什么样的考虑）
        #而n.trees和interaction.depth则根据表现选取出了表现最好的参数，这就是所谓的调参
        
#定制化调参程序
#创建调参网络
gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9), 
                        n.trees = (1:30)*50, 
                        shrinkage = 0.1,
                        n.minobsinnode = 20) #gbmGrid就是一个矩阵，包括模型需要的四个参数（在这里n.minobsinnode选择的20，说明n.minobsinnode的数值还是有影响的）
                                             #这里的trees有从50到1500三十组，depth有三种，共90种参数会被评估
                                             
set.seed(825)
gbmFit2 <- train(Class ~ ., data = training, 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 ## Now specify the exact models 
                 ## to evaluate:
                 tuneGrid = gbmGrid)
gbmFit2 #9，1200表现最好

#画Resampling的图
trellis.par.set(caretTheme())
plot(gbmFit2)  

plot(gbmFit2, metric = "Kappa")

plot(gbmFit2, metric = "Kappa", plotType = "level",
     scales = list(x = list(rot = 90)))
     
ggplot(gbmFit2)  #很重要但没什么好说的，可以用ggplot美化

#################trainControl功能##########################
#methods，有："boot", "cv", "LOOCV", "LGOCV", "repeatedcv", "timeslice", "none" and "oob"，具体是什么需要自己更进一步学习，提到的是oob，out-of-bag estimates，只能被用作random forest, bagged trees, bagged earth, bagged flexible discriminant analysis, or conditional tree forest models，不包括gbm
#number and repeats，
#verboseIter，是否输出training log
#returnData，是否saving the data into a slot called trainingData，不大明白slot是什么
#p，For leave-group out cross-validation
#对于"timeslice"方法, trainControl has options initialWindow, horizon and fixedWindow
#classProbs，在重采样期间是否应计算保留样本的类别概率
#index and indexOut，没懂
#summaryFunction，#这个参数对于性能的表现非常重要，twoClassSummary的话就可以用ROC，灵敏度，特异性，调参，defaultSummary的话准确度，kappa值，prSummary的话，PRAUC，召回率，查准率，回归分析和多重分类同理
#selectionFunction，
#PCAthresh, ICAcomp and k，给preProcess 的参数
#returnResamp，可以给"all", "final" or "none"
#allowParallel，默认应该就是开启的

##################选择表现的度量#############
#By default, RMSE, R2, and the mean absolute error (MAE) are computed for regression while accuracy and Kappa are computed for classification. 
#也可以自己选择度量的公式，如下
fitControl <- trainControl(method = "repeatedcv",number = 10,repeats = 10,
                           classProbs = TRUE, ## Estimate class probabilities
                           summaryFunction = twoClassSummary) #twoClassSummary即可以自己定义的表现度量公式,自己做是搞不定了，看看有没有分享的合适的函数吧。
head(twoClassSummary)#可以查看公式的内容                          

gbmFit3 <- train(Class ~ ., data = training,method = "gbm",trControl = fitControl,verbose = FALSE, 
                 tuneGrid = gbmGrid,metric = "ROC")## Specify which metric to optimize
gbmFit3

#################选择最终的参数###############
#有三种选法，best is chooses the largest/smallest value, 
#oneSE attempts to capture the spirit of Breiman et al (1984) 
#and tolerance selects the least complex model within some percent tolerance of the best value.

whichTwoPct <- tolerance(gbmFit3$results, metric = "ROC", 
                         tol = 2, maximize = TRUE)  #tol选项让我们可以寻找less complex model，For example, to select parameter values based on a 2% loss of performance
gbmFit3$results[whichTwoPct,1:6]

################提取预测和类别概率############
predict(gbmFit3, newdata = head(testing))
predict(gbmFit3, newdata = head(testing), type = "prob")

###############探索和比较Resampling Distributions########3
densityplot(gbmFit3, pch = "|")

#不同模型
svmFit <- train(Class ~ ., data = training, method = "svmRadial", trControl = fitControl, 
                 preProc = c("center", "scale"),tuneLength = 8,metric = "ROC")

rdaFit <- train(Class ~ ., data = training, method = "rda", trControl = fitControl, 
                 tuneLength = 4,metric = "ROC")
                 
resamps <- resamples(list(GBM = gbmFit3,SVM = svmFit,RDA = rdaFit))
summary(resamps)

#已知合适参数，不需要重采样
fitControl <- trainControl(method = "none", classProbs = TRUE)
gbmFit4 <- train(Class ~ ., data = training, method = "gbm", trControl = fitControl, verbose = FALSE, 
                 tuneGrid = data.frame(interaction.depth = 4, n.trees = 100,shrinkage = .1,n.minobsinnode = 20),
                 metric = "ROC")
predict(gbmFit4, newdata = head(testing))
predict(gbmFit4, newdata = head(testing), type = "prob")

#######################之后全是在介绍每个模型，面向的应用和调整的参数什么的##########


############对于CARTE的学习，拿起了又放下，拿起来又放下，持续的学习吧######33
library(mlbench)
data(Sonar)
str(Sonar[, 1:10])
set.seed(998)
colnames(Sonar)[60]<-c('test')# 瞎改一下做一个回归模型
Sonar<-Sonar[,-61]
inTraining <- createDataPartition(Sonar$test, p = .75, list = FALSE)
training <- Sonar[ inTraining,]
testing  <- Sonar[-inTraining,]

lm_fit <- train(test ~ .,
                data = training, 
                trControl = trainControl(method = "repeatedcv",
                                         number = 10,
                                         repeats = 10),
                method = "brnn")
                
test <- predict(lm_fit, testing)
postResample(pred = test, obs = testing$test) #之前不是一直不会做回归吗，其实也不会有太多的不一样，回归模型的评价指标是R^2,和平均方差什么的，ROC,AUC只有分类模型才有。
                                              #predict对测试集做一个运算，postResample则对模型在测试集上的表现做了一个评估，和训练集的表现比较还是存在挺大差异的。





#############另一次尝试#################
fitControl <- trainControl(method = "repeatedcv", #summaryFunction = multiClassSummary, 
                           number = 10,repeats =10,
                           classProbs = TRUE)
rfFit <- train(group ~ ., data = Train, method = "rf", trControl = fitControl, 
                tuneLength = 4,metric = "ROC") #随机森林模型，然而折腾了这么久还是画不出AUC来，或许我该换个小的数据集尝试

importance = varImp(rfFit,scale = FALSE) #可以看到在模型里具体因素所占的比重，这个非常有意思
importance

rf.probs = predict(rfFit,Test) #type = "prob"，不用prob得到的预测结果就是分类的结果，如果用prob，得到的结果就是分类的可能性（？）
rf.probs = predict(rfFit,Test,type = "prob")
postResample(pred = rf.probs, obs = Test$group) #看模型在测试集上的表现

#我的理解不一定对，我也想记录一下
#首先回归模型是没有AUC的，这不是回归模型的评价指标之一，这很好理解。
#其次对于多个结论的分类，默认就是Accuracy和Kappa两个指标，因为twoClassSummary是没法用的，它是针对两个结论的模型适用的评价。
#可是想要求AUC什么的，也是可以做到的，比如https://github.com/topepo/caret/issues/107提供了multiClassSummary，使用multiClassSummary 就可以计算很多东西，包括ROC，F1等等，极大的扩展了评价指标。
#可是想要画图我还是没搞懂要怎么画AUC，roc函数就要求你两水平的结果，那我想做多水平的结果画AUC是否适用？就是一个很纯粹的数学问题了，我不大明白。

fitControl <- trainControl(method = "repeatedcv", summaryFunction =twoClassSummary ,#multiClassSummary,
                           number = 10,repeats =10,
                           classProbs = TRUE)
rfFit <- train(group ~ ., data = Train, method = "rf", trControl = fitControl, 
                tuneLength = 4,metric = "ROC") #我把group换成了两种状况的

rf.probs = predict(rfFit,Test,type = "prob")
postResample(pred = rf.probs, obs = Test$group)
library('pROC')
rf.ROC = roc(response = Test$group,
              predictor = rf.probs$poor, #选哪根好像都没有区别
              levels = levels(Test$group))
plot(rf.ROC,type = "S",col = "red")
plot(rf.ROC,type = "S",col = "blue") #这样就可以把AUC画出来了，不过这个图是不是有点丑？

g <- ggroc(rf.ROC) #之后可以接ggplot进行美化，这个函数对于美学非常有意义
roc.list <- roc(outcome ~ s100b + ndka + wfns, data = aSAH) #也可以对list绘图，此处为示例
g.list <- ggroc(roc.list)



#模型表现
#这里主要是分类模型的表现吧，回归的模型暂时还没怎么做过

mod_per <- function(mod){ #这个function写得还比较粗糙，之后可以改一下
  
rf.probs_1= predict(mod,Test) 
rf.probs_2= predict(mod,Test,type = "prob")
Test$bind<-factor(Test$bind)
p1<-postResample(pred = rf.probs_1, obs = Test$bind) #这里得到精度和kappa系数，
print(p1)

level<-levels(Test$bind)
dat<-data.frame(obs=factor(Test$bind),
                pred=factor(rf.probs_1),
                BIN = rf.probs_2$BIN
)
dat$NOT <- 1 - dat$BIN
p2<-twoClassSummary(dat, lev = level)  #这里得到ROC，灵敏度和特异性
print(p2)

p3<-prSummary(dat, lev = level) #这里得到AUC，查准率，查全率，F1 score
print(p3)
}

mod_per(rfFit)

library('pROC') #ROC曲线
rf.probs = predict(rfFit,Test,type = "prob")
rf.ROC = roc(response = Test$bind,predictor = rf.probs$BIN,levels = levels(Test$bind))

g <- ggroc(list(RF=rf.ROC,GBM=gbm.ROC,SVM=svm.ROC,KNN=knn.ROC)) #如果想把多条ROC画在一起的话，可以用这个

#这里想聊一下几种表现度量
#Recall这里被翻译为了查全率，其实我更喜欢叫召回率，说的是被正确预测的正例在所有正例中的比例，我认为这个度量重要在于，我们的样本在阳性和阴性比例相差较大的情况下，它的意义会比较大。
#Precision查准率，是指被正确预测的正例在所有预测为正例的样本的比例，和recall是同时出现的度量没什么好说的
#F score就是一个将查全率和查准率结合起来的指标
#Kappa系数是一种比例，代表着模型的分类结果与完全随机的分类产生错误减少的比例，kappa在我的理解里是对Accuracy的一种扩展，因为Accuracy其实是不太能很好的表征不平衡数据集的。
#可以看到在我们的demo中，kappa的大小和recall值是有一定的联系的，因为我们的demo中，阴性样本的显著的多于阳性样本，所以在kappa的算法里，查全的能力被看得更重了
#ROC和AUC，只有在二分类模型中，ROC和AUC才适用，AUC是ROC下方的面积，我其实没有太懂prSummary给出的AUC是什么，如果认为AUC是ROC的面积的话，明显twoClassSummary给出的ROC值是正确的，rf.ROC给的值也与之对应。
#我查了一下，twoClassSummary给的是AUROC，prSummary给的则是AUPRC，这下我能够理解了
#AUROC就是我们常用的ROC，以真阳性率，假阳性率做的。AUPRC则是以精确度和召回率做的，大家一般认为在不平衡数据上使用AUPRC更好一些？

#一下子还真没找到好的画AUPRC的包
  probs_1= predict(mod,Test) 
  probs_2= predict(mod,Test,type = "prob")
  Test$bind<-factor(Test$bind)
  
  eva<-cbind(probs_2,Test$bind)
  colnames(eva)[3]<-c('obs')
  eva$group<-(modname)
  eva$group<-as.factor(eva$group)
  eva$obs<-as.factor(eva$obs)
  
  library(PRROC)
  library(ROCR)
  scores <- data.frame(eva$BIN)
  scores$labels<-ifelse(eva$obs=='BIN','1','0')
  pr <- pr.curve(scores.class0=scores[scores$labels=="1",]$eva.BIN,
                 scores.class1=scores[scores$labels=="0",]$eva.BIN,
                 curve=T)
  y <- as.data.frame(pr$curve)
  ggplot(y, aes(y$V1, y$V2))+geom_path()+ylim(0,1) #这个画出来很糙感觉



################附件一：multiClassSummary，这个我用着感觉还不错，虽然后面好像出了花里胡哨的更新的版本，没有数学功底的人想写这么个东西几乎是不可能的#################

multiClassSummary <- function (data, lev = NULL, model = NULL){
  
  #Load Libraries
  require(Metrics)
  require(caret)
  
  #Check data
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
    stop("levels of observed and predicted data do not match")
  
  #Calculate custom one-vs-all stats for each class
  prob_stats <- lapply(levels(data[, "pred"]), function(class){
    
    #Grab one-vs-all data for the class
    pred <- ifelse(data[, "pred"] == class, 1, 0)
    obs  <- ifelse(data[,  "obs"] == class, 1, 0)
    prob <- data[,class]
    
    #Calculate one-vs-all AUC and logLoss and return
    cap_prob <- pmin(pmax(prob, .000001), .999999)
    prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
    names(prob_stats) <- c('ROC', 'logLoss')
    return(prob_stats) 
  })
  prob_stats <- do.call(rbind, prob_stats)
  rownames(prob_stats) <- paste('Class:', levels(data[, "pred"]))
  
  #Calculate confusion matrix-based statistics
  CM <- confusionMatrix(data[, "pred"], data[, "obs"])
  
  #Aggregate and average class-wise stats
  #Todo: add weights
  class_stats <- cbind(CM$byClass, prob_stats)
  class_stats <- colMeans(class_stats)
  
  #Aggregate overall stats
  overall_stats <- c(CM$overall)
  
  #Combine overall with class-wise stats and remove some stats we don't want 
  stats <- c(overall_stats, class_stats)
  stats <- stats[! names(stats) %in% c('AccuracyNull', 
                                       'Prevalence', 'Detection Prevalence')]
  
  #Clean names and return
  names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
  return(stats)
  
}

################附件二：keras in R，caret拥有模型其实并不全面，比如积卷神经网络，就没有，幸运的是现在keras已经有R语言的接口的（对于自己这样的菜鸡能不在python上折腾实在太好了）#####333
################虽然又是个大坑，不过，有好的工具能有什么不满意的呢？######3
install.packages('keras')
library(keras) 
library(tensorflow) #我安装keras的方法有点玄学，我自己都不确定再装一次能不能成功，先是装包

conda create --name r-tensorflow keras tensorflow #然后在服务器上我新建了虚拟环境

#之后我在R-serves的设置里修改了python位置，改为了r-tensorflow这个环境的python位置（envs/r-tensorflow/bin/python）
#然后就可以正常运行了（起码我看起来挺正常的）
#我又用笔记本安装了一次，发现这种安装方式是正确的。话说为什么用笔记本跑比服务器快多了？我的笔记本有那么强吗？可以，我假期就折腾这个了。


##########################一个简单应用的示例在Caret-based_machine_learning_and_immunogenicity_prediction.R########################


##################################################################################################################################
#有一种算法叫mRMR（Max-Relevance and Min-Redundancy），是一种特征选择的算法，在生物标志物的研究种想要用更少的标志物做预测可以用这样的方法做选择。

library(mRMRe)

feature_num = 46 #未经选择的全部特征
train_feature = mlexp[,0:feature_num] #全部特征的矩阵
train_label = as.numeric(as.factor(mlexp$group)) #结果

mrmr_feature<-train_feature
mrmr_feature$y<-train_label

target_indices = which(names(mrmr_feature)=='y')

#转化成mRMR.data的形式才能被使用
Data <- mRMR.data(data = data.frame(mrmr_feature))
#data就是数据，target_indices就是Y（label）值，也就是需要对比相关度的列

#feature_count设置选择特征数，这里有classic mRMR feature selection和mRMR.ensemble两种，没太明白区别都是什么
mrmr=mRMR.ensemble(data = Data, target_indices = target_indices,
                   feature_count = 5, solution_count = 1)
#获得mrmr选择后的特征索引
#获取筛选出来的特征的列，包含在mrmr@filters中，mrmr@filters[原特征个数]这个list里
index=mrmr@filters[[as.character(mrmr@target_indices)]]

#新数据提取
new_data = mlexp[,index] #这个重新做学习的话，结果还不错，可以写个循环批量的做一堆不同特征的情况下的预测能力的比较，然后再做决定。

######################################################################################################
#自己折腾的一个函数，输入特征数和模型就可以自动算出AUC值来，这个函数目前集成得还不够好，有时间可以把input矩阵和label也给集成进去

feature_to_AUC<-function(num,model) {

library(mRMRe)

feature_num = 75
train_feature = scale(learn[,0:feature_num])
train_feature<-as.data.frame(train_feature) #需要变成数据框
train_label = as.numeric(as.factor(lab_learn)) 

mrmr_feature<-train_feature
mrmr_feature$y<-train_label

target_indices = which(names(mrmr_feature)=='y')
Data <- mRMR.data(data = data.frame(mrmr_feature))

mrmr=mRMR.ensemble(data = Data, target_indices = target_indices,
                   feature_count = num, solution_count = 1)

index=mrmr@filters[[as.character(mrmr@target_indices)]]
new_data = train_feature[,index] 
new_data<-as.data.frame(new_data)
new_data$label<-lab_learn

library('caret')
library(doParallel)
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

inTraining <- createDataPartition(new_data$label, p = .75, list = FALSE)
training <- new_data[ inTraining,]
testing  <- new_data[-inTraining,]

fitControl <- trainControl(method = "repeatedcv",summaryFunction =prSummary,
                           number = 10,repeats =10,classProbs = TRUE)

rfFit <- train(label ~ ., data = training, method = model, trControl = fitControl, 
               tuneLength = 8,metric = "AUC") 

rf.probs = predict(rfFit,testing) 
rf.probs_2= predict(rfFit,testing,type = "prob")
testing$label<-factor(testing$label)

level<-levels(testing$label)
dat<-data.frame(obs=factor(testing$label),
                pred=factor(rf.probs),
                T = rf.probs_2$T,
                N = rf.probs_2$N
)

prsum<-prSummary(dat, lev = level)
return(prsum)

}


learn_array<-as.data.frame(array(NA,c(10,6)))
rownames(learn_array) <-c(3:12)
model_list<-c('rf','LogitBoost','svmRadial','gbm','knn','naive_bayes') #本来准备用十个模型做比较的，几个模型有一些小问题，十个给我缩成了六个
colnames(learn_array) <-c(model_list)

for (i in 1:nrow(learn_array)){
  for (j in 1:ncol(learn_array)) {
    
    fe_num<-as.numeric(rownames(learn_array)[i])
    model_name<-colnames(learn_array)[j]

 learn_array[i,j]<-feature_to_AUC(num=fe_num,model=model_name)[1]

 print(paste0('fin! featureNUM:',fe_num,' model:',model_name))
}
}  #就是对六个模型，选择特征数从3到12的情况下模型表现的一个比较。


##############################################################################################
#一个实例
library('caret')
library('doParallel')
library('MetBrewer')
library('ggplot2')
library('pROC')
library('PRROC')
library('ROCR')

cl <- makePSOCKcluster(16)
registerDoParallel(cl)

load('caret_ay.RData')
test[is.na(test)]<-c(0)
caret_ay<-test
caret_ay$pre<-ifelse(caret_ay$label_1=='FALSE','nontumor','tumor')
caret_ay$pre<-factor(caret_ay$pre,levels = c('tumor','nontumor'))
caret_ay<-caret_ay[,-1]

set.seed(618)
inTraining <- createDataPartition(caret_ay$pre, p = .75, list = FALSE)
training <- caret_ay[ inTraining,]
testing  <- caret_ay[-inTraining,]

save(training,testing,file = 'ml_data.RData')
load('ml_data.RData')

fitControl <- trainControl(method = "LGOCV",number = 400, #repeats = 20,
                           classProbs = TRUE,
                           summaryFunction = prSummary)

rfFit <- train(pre ~ ., data = training, method = "rf", trControl = fitControl, 
               tuneLength = 40,metric = "AUC")
print('end_mod_rf')

svmFit <- train(pre ~ ., data = training, method = "svmRadial", trControl = fitControl, 
               tuneLength = 40,metric = "AUC")
print('end_mod_svm')

lbFit <- train(pre ~ ., data = training, method = "LogitBoost", trControl = fitControl, 
               tuneLength = 40,metric = "AUC")
print('end_mod_lb')

knnFit <- train(pre ~ ., data = training, method = "knn", trControl = fitControl, 
               tuneLength = 40,metric = "AUC")
print('end_mod_knn')

nbFit <- train(pre ~ ., data = training, method = "naive_bayes", trControl = fitControl,
               tuneLength = 40,metric = "AUC")
print('end_mod_nb')

save(rfFit,svmFit,lbFit,knnFit,nbFit,file = 'mod.RData')
load('mod.RData') #我不明白，之前做到单个的rf模型不是表现很好吗，为什么这次跑的会变成这样？

#模型表现
rf.probs_1= predict(rfFit,testing) 
rf.probs_2= predict(rfFit,testing,type = "prob")
postResample(pred = rf.probs_1, obs = testing$pre)  

level<-levels(testing$pre)
dat<-data.frame(obs=factor(testing$pre),
                pred=factor(rf.probs_1,levels = c('tumor','nontumor')),
                tumor = rf.probs_2$tumor
)
dat$nontumor <- 1 - dat$tumor

twoClassSummary(dat, lev = level)
prSummary(dat, lev = level) 


rf.probs = predict(rfFit,testing,type = "prob")
rf.ROC = roc(response = testing$pre,
             predictor = rf.probs$nontumor,
             levels = levels(testing$pre))
svm.probs = predict(svmFit,testing,type = "prob")
svm.ROC = roc(response = testing$pre,
             predictor = svm.probs$nontumor,
             levels = levels(testing$pre))
lb.probs = predict(lbFit,testing,type = "prob")
lb.ROC = roc(response = testing$pre,
              predictor = lb.probs$nontumor,
              levels = levels(testing$pre))
knn.probs = predict(knnFit,testing,type = "prob")
knn.ROC = roc(response = testing$pre,
              predictor = knn.probs$nontumor,
              levels = levels(testing$pre))
nb.probs = predict(nbFit,testing,type = "prob")
nb.ROC = roc(response = testing$pre,
              predictor = nb.probs$nontumor,
              levels = levels(testing$pre))


g <- ggroc(list(RF=rf.ROC,SVM=svm.ROC,LB=lb.ROC,KNN=knn.ROC,NB=nb.ROC),size=0.7)

pdf("auc_multimod.pdf",width = 8,height = 6)
g+theme_bw(base_size=18)+
  theme(panel.grid.major =element_blank(), 
        panel.background=element_rect(size =1.1,fill='transparent', color='black'),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        legend.title = element_blank())+
  geom_abline(intercept = 1, slope = 1, linetype = "dashed")+
  scale_x_reverse(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # annotate("text", x=0.5, y=0.3, label=paste0("AUC = ",round(rf.ROC[["auc"]],3)),size=5)+
  theme(aspect.ratio=1)+
  coord_fixed(ratio = 0.8)+
  scale_colour_manual(values=met.brewer("Derain", 5,type ='discrete'))
dev.off()

#然后是PRAUC
auprc <- function(mod,Test,modname){ 
  
  probs_1= predict(mod,Test) 
  probs_2= predict(mod,Test,type = "prob")
  Test$pre<-factor(Test$pre)
  
  eva<-cbind(probs_2,Test$pre)
  colnames(eva)[3]<-c('obs')
  eva$group<-(modname)
  eva$group<-as.factor(eva$group)
  eva$obs<-as.factor(eva$obs)
  
  scores <- data.frame(eva$tumor)
  scores$labels<-ifelse(eva$obs=='tumor','1','0')
  pr <- pr.curve(scores.class0=scores[scores$labels=="1",]$eva.tumor,
                 scores.class1=scores[scores$labels=="0",]$eva.tumor,
                 curve=T) 

  y <- as.data.frame(pr$curve)
  y$mod<-c(modname)
  return(y)
}
  
auprcplot<-rbind(auprc(rfFit,testing,'RF'),
      auprc(svmFit,testing,'SVM'),
      auprc(lbFit,testing,'LB'),
      auprc(knnFit,testing,'KNN'),
      auprc(nbFit,testing,'NB'))
auprcplot$mod<-factor(auprcplot$mod,levels = c('RF','SVM','LB','KNN','NB'))

pdf("auprc_multi_mod.pdf",width = 8,height = 6)
ggplot(auprcplot, aes(V1, V2,colour = mod))+
    geom_path(size=0.7)+ylim(0,1)+
    theme_bw(base_size=18)+
    theme(panel.grid.major =element_blank(), 
          panel.background=element_rect(size =1.1,fill='transparent', color='black'),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    # annotate("text", x=0.3, y=0.5, label=paste0("AUPRC = ",round(pr$auc.integral, 3)),size=5)+
    theme(aspect.ratio=1)+
    labs(x="Recall", y = "Precision")+
    coord_cartesian(ylim=c(0,1),xlim=c(0,1))+
    scale_colour_manual(values=met.brewer("Derain", 5,type ='discrete'))
dev.off()

##########################################################
#AUC的置信区间

predict<-data.frame(response=c(score$group),predictor=c(NA))
predict$response<-factor(predict$response,levels = c('N','C'))
predict$predictor<-predict (fit.full,score,type='response') #这里type=‘response’输出0-1的预测阈值

roc<-plot.roc(predict$response, predict$predictor) 

ci.sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), boot.n=100)
plot(roc)
plot(ci.sp.obj, type="shape", col=NA)

#但是这玩意怎么用ggplot画出来呢？






