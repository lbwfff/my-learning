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
#summaryFunction，#这个参数对于性能的表现非常重要，如果分类且只有两类的话可以用summaryFunction = twoClassSummary，当如果回归或是多类别的话要怎么做呢？
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



