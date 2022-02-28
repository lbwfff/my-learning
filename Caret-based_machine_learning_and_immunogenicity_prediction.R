
################除了antigen.garnish做的几个预测外，再添加几个feature##############################
./netMHC -p ./test/test.pep -xls -xlsfile ./my_NetMHCpan_out.xls
#和一类MHC的亲和性
./netMHCstabpan -f ./test/test.pep -xls -ia -p -xlsfile ./my_NetMHCstapan_out.xls
#和一类MHC结合的稳定性，这个input文件还有点奇怪
#然后可以再补一个二类MHC的亲和性之类的？

###############################################################################################
library(antigen.garnish) #使用这个包，做一个肽的免疫原性的预测
list<-read.csv('~/tools/antigen.garnish/pnas.1500973112.sd01.csv') #训练数据
list<-list[!duplicated(list$Epitope),] #它的数据里有一些序列一致的肽，但是有着不同的免疫性，这个我不是太理解为什么。
list<-list[,c(1,7)]

v <- list$Epitope
f_score <- v %>% foreignness_score(db = "human") %>% print 
d_score <- v %>% dissimilarity_score(db = "human") %>% print #预测的两种特征，后期我还想把另外两种工具得到的特征也加入学习中，此处只做测试
save(f_score,d_score,file='~/tools/antigen.garnish/score.RData') #上面两句还跑了挺久的，先存一波

f_score<-f_score[match(list$Epitope,f_score$nmer),]
d_score<-d_score[match(list$Epitope,d_score$nmer),]
list$foreignness_score<-f_score$foreignness_score
list$dissimilarity_score<-d_score$dissimilarity

list [is.na (list)] <- 0 #其实我每太明白软件没给出的数据算什么，是说他完全不具有外来性吗？
rownames(list)<-list$Epitope
list<-list[,-1]

library('caret')
library(doParallel)
cl <- makePSOCKcluster(16)
registerDoParallel(cl)
# stopCluster(cl) 

trainIndex <- createDataPartition(list$Immunogenicity, p = .75, list = FALSE, times = 1) #分割数据
Train <- list[ trainIndex,]
Test  <- list[-trainIndex,]

fitControl <- trainControl(method = "repeatedcv",number = 10,repeats = 10,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary)
rfFit <- train(Immunogenicity ~ ., data = Train, method = "rf", trControl = fitControl, 
               tuneLength = 4,metric = "ROC")

importance = varImp(rfFit,scale = FALSE) #可以看到dissimilarity_score对于模型的重要性要显著高于foreignness_score，和先前的报导一致

rf.probs_1= predict(rfFit,Test) 
rf.probs_2= predict(rfFit,Test,type = "prob")
postResample(pred = rf.probs_1, obs = Test$Immunogenicity) #用postresample得到的信息还是比较少，所以我写了后面这段

level<-levels(Test$Immunogenicity)
dat<-data.frame(obs=factor(Test$Immunogenicity),
                pred=factor(rf.probs_1),
                Positive = rf.probs_2$Positive
)
dat$Negative <- 1 - dat$Positive
twoClassSummary(dat, lev = level) #用twoClassSummary来做可以得到ROC值，这里发现ROC值其实想当不错了

library('pROC')
rf.probs = predict(rfFit,Test,type = "prob")
rf.ROC = roc(response = Test$Immunogenicity,
             predictor = rf.probs$Positive,
             levels = levels(Test$Immunogenicity))
plot(rf.ROC,type = "S",col = "blue")
g <- ggroc(rf.ROC) #画ROC曲线什么的

####用第三方的数据做一个测试
tri_p<-readAAStringSet('~/tools/antigen.garnish/IEDB_positive_T-cell_assays.fasta')
tri_n<-readAAStringSet('~/tools/antigen.garnish/IEDB_negative_T-cell_assays.fasta')

dat_p<-data.frame(seq=gsub("\\.","",paste(tri_p)),
                group=c('positive'))
dat_n<-data.frame(seq=gsub("\\.","",paste(tri_n)),
                  group=c('negative'))
dat_tri<-rbind(dat_p,dat_n)

tri_f <- dat_tri$seq %>% foreignness_score(db = "human") %>% print
tri_d <- dat_tri$seq %>% dissimilarity_score(db = "human") %>% print

save(tri_d,tri_f,file = '~/tools/antigen.garnish/tri_test.RData')

tri_f<-tri_f[match(dat_tri$seq,tri_f$nmer),]
tri_d<-tri_d[match(dat_tri$seq,tri_d$nmer),]
dat_tri$foreignness_score<-tri_f$foreignness_score
dat_tri$dissimilarity_score<-tri_d$dissimilarity

dat_tri [is.na (dat_tri)] <- 0
rownames(dat_tri)<-dat_tri$seq
dat_tri<-dat_tri[,-1]

rf.probs_1= predict(rfFit,dat_tri) 
rf.probs_2= predict(rfFit,dat_tri,type = "prob")

level<-levels(factor(dat_tri$group))
rf.probs_1<- tolower(rf.probs_1) #这里发现自己被大小写坑了，所以做了一些修改，这些东西应该在一开始就注意的
dat<-data.frame(obs=factor(dat_tri$group),
                pred=factor(rf.probs_1),
                positive = rf.probs_2$Positive) 
dat$negative <- 1 - dat$positive
twoClassSummary(dat, lev = level) #这里发现模型有一些问题，仅看ROC值还不错，但是特异性特别好（0.9956998），但是灵敏度特别低（0.4474257），为什么会有这么奇葩的结果？
                                  #可能是数据的问题，我对于后来用到的第三方数据其实没有特别的了解，需要看一下文章了解具体的数据的细节。

rf.probs = predict(rfFit,dat_tri,type = "prob")
colnames(rf.probs)<-c('negative','positive')
rf.ROC = roc(response = dat_tri$group,
             predictor = rf.probs$positive,
             levels = levels(dat_tri$group))
plot(rf.ROC,type = "S",col = "blue")
g <- ggroc(rf.ROC) #对这个数据的ROC
