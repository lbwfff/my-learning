####################################################################################
#为什么学逻辑回归呢，因为在肿瘤生物标志物的研究中，想要把几种标志物的组合计算出具体的分数出来
#如果用决策树模型的话，虽然预测表现很好，但是模型终究还是想要变成一个具体的公式，会更加适合临床一些
#然后很多生物统计学的东西，之前没有系统学习过，近段时间确实极大的开拓了我对视野，OR值,还有其它什么乱七八糟的，生物统计学上的指标，我觉得做生物标志物的话，都是可以做进去的

rm (list = ls ())

data(Affairs,package = "AER")
df<-Affairs
df<-model.matrix(affairs ~ ., data = df)[,-1] #我不太明白为什么把因变量全部变成1了，有点不太理解这个函数做的处理

df<-as.data.frame(cbind(Affairs$affairs,df))
  
df$affairs<-ifelse(df$V1>0,1,0)
df$affairs<-factor(df$affairs,
                     levels = c(0,1),
                     labels = c("No","Yes"))  
#这里的水平非常重要，这里应该先放NO再放Yes，对于结果影响很大
#倒也没有很难理解，就是小值放在前面，大值放在后面就好了

df<-df[,-1]
fit.full<-glm(affairs~.,
              data=df,family = binomial()) #建模

fit.result<-summary(fit.full)

df4<-as.data.frame(exp(cbind(OR = coef(fit.full), confint(fit.full))))
df4$p<-as.numeric(fit.result$coefficients[,4])
  
df4<-data.frame(df4[-1,c(1,4,2,3)]) #Estimate是OR吗？不是吧
df4$Var<-rownames(df4)
colnames(df4)<-c("OR","Pvalue","OR_1","OR_2","Var")

df5<-df4[,c(5,1,2,3,4)]
df5$OR_mean<-df5$OR
df5$OR<-paste0(round(df5$OR,2),
               "(",
               round(df5$OR_1,2),
               "~",
               round(df5$OR_2,2),
               ")")
df5$Pvalue<-signif(df5$Pvalue,3)

df5<-df5[order(df5$OR_mean,decreasing = T),]
library(forestplot)

#森林图
forestplot(labeltext=as.matrix(df5[,1:3]),
           mean=df5$OR_mean,
           lower=df5$OR_1,
           upper=df5$OR_2,
           zero=1,
           boxsize=0.2,
           lineheight = unit(7,'mm'),
           colgap=unit(2,'mm'),
           lwd.zero=1.5,
           lwd.ci=2, 
           col=fpColors(box='#458B00',
                        summary='#8B008B',
                        lines = 'black',
                        zero = '#7AC5CD'),
           xlab="OR",
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                            xlab  = gpar(cex = 0.8),
                            cex = 0.9),
           lty.ci = "solid",
           title = "Forestplot", 
           line.margin = 0.08,
           graph.pos=2)

#这个比较靠谱的把森林图画出来了
#我把他的代码改了一下，因为我觉得他OR算错了

#这个森林我觉得是没有问题的，拟合曲线的话，可能是我理解得不够，还有很多疑问
#那么现在我想的是怎么在这样的逻辑回归模型里做更深层次的可视化

#把分数计算出来，看看分数的森林图

Coeffs<-coef(fit.full)[-1]

score<-df
score$score<-c(NA)

for(i in 1:nrow(score)){
  scole_cal<-list()
    for (j in 1:length(Coeffs)){
      scole_cal[j]<-(Coeffs[j]*score[i,which(colnames(score)==names(Coeffs)[j])])
    }
  scole_cal<-as.numeric(scole_cal)
  score$score[i]<-sum(scole_cal)
}

#
fit.full<-glm(affairs~score,
              data=score,family = binomial()) #并不能放在一起算，这相当于是一次建模了，想要比较的话应该是分开建模？

fit.result<-summary(fit.full)

df4<-as.data.frame(exp(cbind(OR = coef(fit.full), confint(fit.full))))
df4$p<-as.numeric(fit.result$coefficients[,4])

df4<-data.frame(df4[-1,c(1,4,2,3)]) #Estimate是OR吗？不是吧
df4$Var<-rownames(df4)
colnames(df4)<-c("OR","Pvalue","OR_1","OR_2","Var")

df5<-df4[,c(5,1,2,3,4)]
df5$OR_mean<-df5$OR
df5$OR<-paste0(round(df5$OR,2),
               "(",
               round(df5$OR_1,2),
               "~",
               round(df5$OR_2,2),
               ")")
df5$Pvalue<-signif(df5$Pvalue,3)

df5<-df5[order(df5$OR_mean,decreasing = T),]

forestplot(labeltext=as.matrix(df5[,1:3]),
           mean=df5$OR_mean,
           lower=df5$OR_1,
           upper=df5$OR_2,
           zero=1,
           boxsize=0.2,
           lineheight = unit(7,'mm'),
           colgap=unit(2,'mm'),
           lwd.zero=1.5,
           lwd.ci=2, 
           col=fpColors(box='#458B00',
                        summary='#8B008B',
                        lines = 'black',
                        zero = '#7AC5CD'),
           xlab="OR",
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                            xlab  = gpar(cex = 0.8),
                            cex = 0.9),
           lty.ci = "solid",
           title = "Forestplot", 
           line.margin = 0.08,
           graph.pos=2)

#
fit.full<-glm(affairs~score,
              data=score,family = binomial()) 
S <- coefficients(summary(fit.full)) #coefficients为1，是不是都不需要计算了

seq <- seq(min(score$score), max(score$score), 0.1) #这里应该是把age分成很多段了，每一段是0.1的跨度
dt <- data.frame(score = seq)
dt$OR <- exp(dt$score * S[2, 1])
dt$OR_lwr <- exp(dt$score * S[2, 1] - qnorm(0.975) * dt$score * S[2, 2]) #这里取置信区间的方法很有意思
dt$OR_upr <- exp(dt$score * S[2, 1] + qnorm(0.975) * dt$score * S[2, 2])

ggplot(data = dt, aes(x = score)) + 
  geom_line(aes(y = OR, color = "OR"), linetype = 2) + 
  geom_line(aes(y = OR_lwr, color = "OR_lwr"), linetype = 2) + 
  geom_line(aes(y = OR_upr, color = "OR_upr"), linetype = 2) + 
  labs(color = "") + 
  theme(legend.position = c(0.8, 0.8)) #但是这里的OR为什么也是从0开始的

#主要没太理解，为什么两端都不贴1呢？不是很奇怪吗

#之后我想要知道的是，我应该选择那个分数作为阈值能得到好的预测效果呢？
#我的想法是，把seq段里面每一个区域都算出AUC，灵敏度和准确度就好了？

predict<-data.frame(response=c(score$affairs),predictor=c(NA))

predict$response<-factor(predict$response,levels = c('No','Yes'))
predict$predictor<-predict (fit.full,score)

roc<-roc(response = predict$response,
         predictor = predict$predictor,
          levels = levels(predict$response)) 

#其实这个roc对象里面已经包括了不同阈值情况下的，灵敏度，准确率什么的
#但还是没明白，要怎么选取阈值呢？以及predictor的分和计算的score还是不一样的,好像加了一个截距？

plot(roc,legacy.axes = TRUE,thresholds="best", print.thres="best")
#这么画可以得到最佳的阈值

g <- ggroc(roc) 

thre<-data.frame(thresholds=c(roc$thresholds),
                 sensitivities=c(roc$sensitivities),
                 specificities=c(roc$specificities))
  
thre<-thre[-c(1,570),]
thre$Youden<-(thre$sensitivities+thre$specificities-1)

#和我理解的一样，Youden值最大的点是选择的阈值

predict2<-data.frame(response=c(score$affairs),predictor=c(NA))
predict2$response<-factor(predict2$response,levels = c('No','Yes'))

predict2$predictor<-ifelse((score$score+1.377>(-0.96)),'Yes','No') 
#这个水平和我想象中的完全不同。。。。。，不应该分高代表Yes吗

postResample(pred = predict2$predictor, obs = predict2$response)
#使用这个阈值时准确率为0.71





#
#用树模型随便评估一下
{
library('caret')
trainIndex <- createDataPartition(df$affairs, p = .75, list = FALSE, times = 1) #分割数据
Train <- df[ trainIndex,]
Test  <- df[-trainIndex,]

fitControl <- trainControl(method = "repeatedcv",number = 10,repeats = 10,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary)
rfFit <- train(affairs ~ ., data = Train, method = "rf", trControl = fitControl, 
               tuneLength = 4,metric = "ROC")

importance = varImp(rfFit,scale = FALSE) 

rf.probs_1= predict(rfFit,Test) 
rf.probs_2= predict(rfFit,Test,type = "prob")

level<-levels(Test$affairs)
dat<-data.frame(obs=factor(Test$affairs),
                pred=factor(rf.probs_1),
                Yes = rf.probs_2$Yes)
dat$No <- 1 - dat$Yes
twoClassSummary(dat, lev = level) 

#这玩意拿决策树做，效果也挺一般的，我发现树模型里面权重高的因素和回归里面OR显著的还不太一样？

trainIndex <- createDataPartition(score$affairs, p = .75, list = FALSE, times = 1) #分割数据
Train <- score[ trainIndex,]
Test  <- score[-trainIndex,]

fitControl <- trainControl(method = "repeatedcv",number = 10,repeats = 10,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary)
rfFit <- train(affairs ~ ., data = Train, method = "rf", trControl = fitControl, 
               tuneLength = 4,metric = "ROC")

importance = varImp(rfFit,scale = FALSE) #倒是符合预期

rf.probs_1= predict(rfFit,Test) 
rf.probs_2= predict(rfFit,Test,type = "prob")

level<-levels(Test$affairs)
dat<-data.frame(obs=factor(Test$affairs),
                pred=factor(rf.probs_1),
                Yes = rf.probs_2$Yes)
dat$No <- 1 - dat$Yes
twoClassSummary(dat, lev = level) #灵敏度巨低
}
