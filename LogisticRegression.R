####################################################################################
#为什么学逻辑回归呢，因为在肿瘤生物标志物的研究中，想要把几种标志物的组合计算出具体的分数出来
#如果用决策树模型的话，虽然预测表现很好，但是模型终究还是想要变成一个具体的公式，会更加适合临床一些
#然后很多生物统计学的东西，之前没有系统学习过，近段时间确实极大的开拓了我对视野，OR值,还有其它什么乱七八糟的，生物统计学上的指标，我觉得做生物标志物的话，都是可以做进去的

data(Affairs,package = "AER")
df<-Affairs

df$ynaffairs<-ifelse(df$affairs>0,1,0)
df$ynaffairs<-factor(df$ynaffairs,
                     levels = c(0,1),
                     labels = c("No","Yes")) #ynaffairs是预测的目标

fit.full<-glm(ynaffairs~gender+age+yearsmarried+
                children+religiousness+education+occupation+rating,
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
