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

#######################################################################################################
#一个实例

#最开始我还是用Boruta做了筛选，这一部分其实也可以不做或许直接LASSO就行了，可以试一下效果会不会有改变，我也不确定

library(Boruta)
set.seed(20230508)

Var.Selec<-Boruta(group~., data=CVNarray,maxRuns=500, doTrace=1)

Var.Selec
Var.Selec$finalDecision
plotImpHistory(Var.Selec)

test<-attStats(Var.Selec)

getSelectedAttributes(Var.Selec,withTentative=FALSE)

pdf("C_VS_N_boruta.pdf",width = 6,height = 8)
Borutaplot(Var.Selec) #这里把绘图略作了修改，把未接受的变量没有展示出来
dev.off()

######################################
library('caret')
library('doParallel')

cl <- makePSOCKcluster(8)
registerDoParallel(cl)

sel<-CVNarray[,colnames(CVNarray) %in% c(rownames(test)[test$decision=='Confirmed'],'group')] #表达矩阵

for (i in 1:(ncol(sel)-1)){
  sel[,i]<-sel[,i]/sd(sel[,i])
} #这里做了标准化，除了标准差

#这里是否需要测试集，我一直没想清楚，主要是样本太少了，先不管测试集做着了
train_control <- trainControl(method = "cv",number = 10, 
                              classProbs = TRUE,summaryFunction = defaultSummary)

linear_model = train(
  group ~ ., data = sel,
  trControl = train_control,
  method = "glmnet", tuneLength = 40,
  family = "binomial",metric = "Accuracy")
#LASSO回归的结果非常随机，感觉每次都不太一样，最好的几个指标倒是一直在那里，但是三名以外感觉就很随机了，感觉还是样本量太少了

pdf("CVN_LASSO_lambda.pdf",width = 5,height = 4)
plot(linear_model$finalModel, xvar = "lambda", label = T)
dev.off()

pdf("CVN_LASSO_dev.pdf",width = 5,height = 4)
plot(linear_model$finalModel, xvar = "dev", label = T)
dev.off()

pdf("CVN_LASSO_varImp.pdf",width = 5,height = 4)
plot(varImp(linear_model, scale = F))
dev.off()

#LASSO的几个套图，感觉可以画在一张里面。

lassocoef<-as.matrix(coef(linear_model$finalModel, s = (linear_model[["bestTune"]][1,2])))
lassocoef<-as.data.frame(lassocoef)
lassocoef$name<-rownames(lassocoef)
lassocoef<-lassocoef[order(abs(lassocoef$s1),decreasing = T),]
keep_X <- rownames(lassocoef)[lassocoef$s1!=0]
keep_X <- keep_X[!keep_X == "(Intercept)"]

# selarray <- sel[,c(keep_X,'group')]
selarray<-CVNarray[,c(keep_X,'group')] #选择的指标

#逻辑回归
selarray$group<-factor(sel$group,
                   levels = c('N','C'),
                   labels = c("health","cancer"))  #注意水平的选取
df2<-selarray 

test_catch<-data.frame(feature=c(1:(ncol(df2)-1)),
                       warning=c(NA))
catchcoef<-lassocoef[-1,]
catchcoef<-catchcoef[catchcoef$s1!=0,]

for (i in 1:nrow(test_catch)) {
  catch<-catchcoef$name[1:i]
  catch<-df2[,c(catch,'group')]
tryCatch(
  { fit.full<-glm(group~.,data = catch,
                 family = binomial) },
  warning = function(w) {
    message(paste0('feature_n = ',i," result : warning"))
    test_catch$warning[i]<-c('warning')
  }
)
}
#这里我是这样做的，从一个变量开始加变量，一直到模型过拟合为止，选取过拟合前的变量数。我不确定这样做是不是合理的，但是因为数据的问题我想不到更好的方法做回归了。
#可以了解一下逐步回归的思想，感觉我的想法很接近，也可以做做看。

#重新建模的时候我想用没有标准化的矩阵做
df3<-CVNarray[colnames(CVNarray) %in% c(catchcoef$name[1:6],'group')] #最终选了6个指标
df3$group<-factor(df3$group,levels = c('N','C'))
fit.full<-glm(group~.,data = df3,
              family = binomial)
fit.result<-summary(fit.full)

p<-predict(fit.full,type='response')
qplot(seq(-2,2,length=58),sort(p),col="response")

df4<-as.data.frame(exp(cbind(OR = coef(fit.full), confint(fit.full)))) #, level =0.8
df4$p<-as.numeric(fit.result$coefficients[,4])

df4<-data.frame(df4[-1,c(1,4,2,3)]) #Estimate是OR吗？不是吧
df4$Var<-rownames(df4)

for (i in 1:nrow(df4)){
  exp<-CVNarray[,which(colnames(CVNarray)==df4$Var[i])]
  unit.change<-sd(exp)
  df4[i,c(1,3,4)]<-df4[i,c(1,3,4)]^unit.change
}  #在计算OR时应该标准化，但这组数据种即使做了标准化，一些OR值也还是很奇怪

colnames(df4)<-c("OR","Pvalue","OR_1","OR_2","Var")

df5<-df4[,c(5,1,2,3,4)]
df5$OR_mean<-df5$OR
df5$OR<-paste0(round(df5$OR,2),"(",round(df5$OR_1,2),"~",round(df5$OR_2,2),")")
df5$Pvalue<-signif(df5$Pvalue,3)

df5<-df5[order(df5$OR_mean,decreasing = T),]

#森林图
pdf("C_VS_N_OR_mul.pdf",width = 7,height = 5)
forestplot(labeltext=as.matrix(df5[,1:3]),
           mean=df5$OR_mean,
           lower=df5$OR_1,
           upper=df5$OR_2,
           zero=1,boxsize=0.2,
           lineheight = unit(7,'mm'),
           colgap=unit(2,'mm'),
           lwd.zero=1.5,lwd.ci=2, 
           col=fpColors(box='#458B00',summary='#8B008B',lines = 'black',zero = '#7AC5CD'),
           xlab="OR",lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),xlab  = gpar(cex = 0.8),cex = 0.9),
           lty.ci = "solid",
           title = "Forestplot", 
           line.margin = 0.08,
           graph.pos=2)
dev.off()
#这是一个单模型的OR值森林图，OR值很奇怪，P值也很奇怪

#然后画了一个多模型的OR值森林
oddsradio<-data.frame(Var=c(NA),OR=c(NA),Pvalue=c(NA),OR_1=c(NA),OR_2=c(NA),OR_mean=c(NA))

for (i in 1:6){
  
  varn<-catchcoef$name[1:6][i]
  glmarray<-CVNarray[,colnames(CVNarray) %in% c(varn,'group')]
  glmarray$group<-factor(glmarray$group,levels = c('N','C'))
  
  colnames(glmarray)[1]<-c('var')
  
  fit.full<-glm(group~var,
                data=glmarray,family = binomial()) 
  
  fit.result<-summary(fit.full)
  df4<-as.data.frame(exp(cbind(OR = coef(fit.full), confint(fit.full))))
  df4$p<-as.numeric(fit.result$coefficients[,4])
  
  df4<-data.frame(df4[-1,c(1,4,2,3)]) 
  df4$Var<-rownames(df4)
  
  unit.change<-sd(glmarray$var)
  df4[i,c(1,3,4)]<-df4[i,c(1,3,4)]^unit.change
  
  colnames(df4)<-c("OR","Pvalue","OR_1","OR_2","Var")
  
  df5<-df4[,c(5,1,2,3,4)]
  df5$OR_mean<-df5$OR
  df5$OR<-paste0(signif(df5$OR,2),"(",signif(df5$OR_1,2),"~",signif(df5$OR_2,2),")")
  df5$Pvalue<-signif(df5$Pvalue,3)
  
  df5<-df5[order(df5$OR_mean,decreasing = T),]
  df5$Var<-c(varn)
  
  oddsradio<-rbind(oddsradio,df5)
}


oddsradio<-oddsradio[!is.na(oddsradio$Pvalue),]

forestplot<-oddsradio[order(oddsradio$OR_mean,decreasing = T),]
forestplot$Pvalue<-as.character(symnum(forestplot$Pvalue,cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                       symbols = c("***", "**", "*", ".", "NS")))
#这在一张图里，我把P值换成了星星
#
pdf("C_VS_N_OR_sig.pdf",width = 7,height = 5)
forestplot(labeltext=as.matrix(forestplot[,1:3]),
           mean=forestplot$OR_mean,
           lower=forestplot$OR_1,
           upper=forestplot$OR_2,
           zero=1,boxsize=0.2,
           lineheight = unit(7,'mm'),
           colgap=unit(2,'mm'),
           lwd.zero=1.5,lwd.ci=2, 
           col=fpColors(box='#458B00',summary='#8B008B',lines = 'black',zero = '#7AC5CD'),
           xlab="OR",
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),xlab  = gpar(cex = 0.8),cex = 0.9),
           lty.ci = "solid",
           title = "Forestplot", 
           line.margin = 0.08,
           graph.pos=2)
dev.off()
#会比单模型的要好得多

#
df3<-CVNarray[colnames(CVNarray) %in% c(catchcoef$name[1:6],'group')]

df3$group<-factor(df3$group,levels = c('N','C'))
fit.full<-glm(group~.,data = df3,
              family = binomial)
fit.result<-summary(fit.full)

Coeffs<-coef(fit.full)[-1]
score<-as.data.frame(CVNarray[,names(Coeffs)])
score$score<-c(NA)

for(i in 1:nrow(score)){
  scole_cal<-list()
  for (j in 1:length(Coeffs)){
    scole_cal[j]<-(Coeffs[j]*score[i,which(colnames(score)==names(Coeffs)[j])])
  }
  scole_cal<-as.numeric(scole_cal)
  score$score[i]<-sum(scole_cal)
}
#根据逻辑回归建模得出风险打分

score$group<-CVNarray$group
score$group<-factor(score$group,levels = c('N','C'))

fit.full<-glm(group~score,
              data=score,family = binomial()) 
S <- coefficients(summary(fit.full)) 

seq <- seq(min(score$score), max(score$score), 0.1) 
dt <- data.frame(score = seq)
dt$OR <- exp(dt$score * S[2, 1])
dt$OR_lwr <- exp(dt$score * S[2, 1] - qnorm(0.975) * dt$score * S[2, 2]) 
dt$OR_upr <- exp(dt$score * S[2, 1] + qnorm(0.975) * dt$score * S[2, 2])

ggplot(data = dt, aes(x = score)) + 
  geom_line(aes(y = OR, color = "OR"), linetype = 2) + 
  geom_line(aes(y = OR_lwr, color = "OR_lwr"), linetype = 2) + 
  geom_line(aes(y = OR_upr, color = "OR_upr"), linetype = 2) + 
  labs(color = "") + 
  theme(legend.position = c(0.8, 0.8)) #看了一眼风险分数的OR，很奇怪自己也看不懂，没有展示出来

predict<-data.frame(response=c(score$group),predictor=c(NA))

predict$response<-factor(predict$response,levels = c('N','C'))
predict$predictor<-predict (fit.full,score)

roc<-roc(response = predict$response,
         predictor = predict$predictor,
         levels = levels(predict$response)) 

plot(roc,legacy.axes = TRUE,thresholds="best", print.thres="best")

g <- ggroc(roc) 
#ROC曲线
pdf("C_VS_N_auc.pdf",width = 6,height = 4)
g+theme_bw(base_size=18)+
  theme(panel.grid.major =element_blank(), 
        panel.background=element_rect(size =1.1,fill='transparent', color='black'),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  geom_abline(intercept = 1, slope = 1, linetype = "dashed")+
  scale_x_reverse(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  annotate("text", x=0.3, y=0.3, label="AUC = 0.984",size=5)+
  theme(aspect.ratio=1)+
  coord_fixed(ratio = 0.8)
dev.off()

thre<-data.frame(thresholds=c(roc$thresholds),
                 sensitivities=c(roc$sensitivities),
                 specificities=c(roc$specificities))

thre<-thre[-c(1,nrow(thre)),]
thre$Youden<-(thre$sensitivities+thre$specificities-1)

predict2<-data.frame(response=c(score$group),predictor=c(NA))
predict2$response<-factor(predict2$response,levels = c('N','C'))

predict2$predictor<-ifelse((score$score*S[2,1]+S[1,1]<(0.632)),'N','C') 

#这里看到准确度
postResample(pred = predict2$predictor, obs = predict2$response)

#分数阈值对于灵敏度和特异性的影响
p<-list()
p[[1]]<-
ggplot(data=thre, aes(x=thresholds, y=sensitivities)) +
  geom_line(size=1)+theme_classic(base_size=18)+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=12,face='bold'))+
  theme(aspect.ratio=1)+
  coord_fixed(ratio = 0.8)

p[[2]]<-
ggplot(data=thre, aes(x=thresholds, y=specificities)) +
  geom_line(size=1)+theme_classic(base_size=18)+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=12,face='bold'))+
  theme(aspect.ratio=1)+
  coord_fixed(ratio = 0.8)

pdf("C_VS_N_thre.pdf",width = 8,height = 4)
patchwork::wrap_plots(p,nrow=1, guides="collect") 
dev.off()

#六个指标各自的箱线图
p<-list()

for (i in 1:6){
  varn<-catchcoef$name[1:6][i]
  glmarray<-CVNarray[,colnames(CVNarray) %in% c(varn,'group')]
  glmarray$group<-factor(glmarray$group,levels = c('N','C'))
  colnames(glmarray)[1]<-c('exp')
  glmarray$fac<-c(varn)
  
  p[[i]]<-
    ggplot(data=glmarray,aes(x=group,y=exp))+
    geom_boxplot(aes(color=group)) + theme_bw()+
    geom_point(aes(color=group)) + 
    labs(y="exp",x=NULL)+
    scale_color_manual(values=met.brewer(name="VanGogh1",n=7,type="discrete")[c(6,3)])+
    geom_signif(comparisons = list(c("N", "C")),map_signif_level=TRUE,vjust=1.5)+
    theme(aspect.ratio=2)+
    facet_grid(.~fac)
}

pdf("C_VS_N_box.pdf",width = 10,height = 6)
patchwork::wrap_plots(p,nrow=2, guides="collect") 
dev.off()

#分数的箱线图
pdf("C_VS_N_score_box.pdf",width = 3,height = 3)
ggplot(data=score,aes(x=group,y=score))+
  geom_boxplot(aes(color=group)) + theme_bw()+
  geom_point(aes(color=group)) + 
  labs(y="score",x=NULL)+
  scale_color_manual(values=met.brewer(name="VanGogh1",n=7,type="discrete")[c(6,3)])+
  geom_signif(comparisons = list(c("N", "C")),map_signif_level=TRUE,vjust=1.5)+
  theme(aspect.ratio=2)
dev.off()

#怎么说呢，感觉做了很多东西，但还是有很多疑问，需要持续的学习



















