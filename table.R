#表格的前处理
test<-read.csv('HCC_meta.csv')
# test<-test[test$Profiling.data.QC.status=='Pass',]

ms<-read.csv('orf_intensity.csv')
ms<-colnames(ms)[-1]
ms<-unique(substr(ms,11,14))
test<-test[test$Case.No. %in% ms,] 

########################
#之前好像很少有制作表格，使用table能够很好的做出表格来，而且外观还比较让人满意，也可以做很多自定义的调节
library(table1)

test$Gender<-factor(test$Gender, levels=c("M", "F"), labels=c("Male", "Female"))
test$HBV<-factor(test$HBV,levels = c('P','N'),label=c('Positive','Negative'))
test$HCV<-factor(test$HCV,levels = c('P','N'),label=c('Positive','Negative'))
test$Liver.cirrhosis<-factor(test$Liver.cirrhosis,levels = c(1,0),label=c('Yes','No'))
test$BCLC.stage<-factor(test$BCLC.stage,levels = c('0','A'),label=c('BCLC level 0','BCLC level A'))
test$Lymphatic.metastasis<-factor(test$Lymphatic.metastasis,levels = c(1,0),label=c('Yes','No'))
test$Macrovascular.invasion<-factor(test$Macrovascular.invasion,levels = c(1,0),label=c('Yes','No'))
test$Microvascular.invasion<-factor(test$Microvascular.invasion,levels = c(1,0),label=c('Yes','No'))
test$Cancer.recurrence<-factor(test$Cancer.recurrence,levels = c(1,0),label=c('Yes','No'))
test$Died.of.recurrence<-factor(test$Died.of.recurrence,levels = c(1,0),label=c('Yes','No')) #如果数据是0，1的话，可以使用as.logical()，会自动变成T or F

label(test$AFP..ng.ml.)<- "AFP (ng/ml)" #可以对于colname进行名字的修改（我不知道这样描述是否准确）

# pdf('table.pdf',width =12,height = 8) #这玩意没法直接保存为pdf，可以输出为html格式再打印为pdf文件
table1(~ Age+HBV + HCV+Liver.cirrhosis+AFP..ng.ml.+BCLC.stage+
         Lymphatic.metastasis+Macrovascular.invasion+Microvascular.invasion+
         Cancer.recurrence+Died.of.recurrence| Gender, data=test)  
# dev.off()
