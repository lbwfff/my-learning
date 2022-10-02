##单因数的cox分析的话，自己也尝试了很多包去画，forestplot，forestploter，首先很麻烦，需要自己整理矩阵，画出来图其实感觉还挺一般的

library('forestmodel')
library("survival")
library("survminer")

vars_for_table<-c('Age','Gender','Liver.cirrhosis','Microvascular.invasion','AFP..ng.ml.','group','Survival_score')
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(Total.follow.up.period..m., Died.of.recurrence)~', x))) #单次的cox回归分析循环
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = meta)}) #把循环分析的结果整合在一起

pdf('forest_model_allfactor.pdf',width =12,height = 8)
print(forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T)) #绘图，forest_model做出来这个还行
dev.off()

vars_for_table<-c('Age','Gender','Liver.cirrhosis','Microvascular.invasion','AFP..ng.ml.','group','Recurrence_score')
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(Disease.free.survival..m., Cancer.recurrence)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = meta)})

pdf('forest_model_allfactor_REC.pdf',width =12,height = 8)
print(forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T))
dev.off()

##多因素的话，就用survminer的ggforest就好了，主要我感觉这个做出来效果还不错，就是每张的长宽比感觉都不一样，和dotplot一个毛病

res.cox<-coxph(Surv(Total.follow.up.period..m.,Died.of.recurrence)~Age+Gender+Liver.cirrhosis+Microvascular.invasion+AFP..ng.ml.+group+Survival_score, data=meta)
summary(res.cox)

pdf('sev_forestplot_allfactor.pdf',width =8,height = 6)
ggforest(res.cox,data = meta,main = "Survival probability",
         fontsize=0.8,noDigits=2,cpositions = c(0.02, 0.22, 0.4))
dev.off()

res.cox2<-coxph(Surv(Disease.free.survival..m.,Cancer.recurrence)~Age+Gender+Liver.cirrhosis+Microvascular.invasion+AFP..ng.ml.+group+Recurrence_score, data=meta)
summary(res.cox2)

pdf('rec_forestplot_allfactor.pdf',width =8,height = 6)
ggforest(res.cox2,data = meta,main = "Recurrence probability",
         fontsize=0.8,noDigits=2,cpositions = c(0.02, 0.22, 0.4)) #看了两个包，感觉都挺一般的，还不如直接ggforest
dev.off()
