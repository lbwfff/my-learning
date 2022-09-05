#一个对RNA-seq数据做免疫浸润程度预测的包，感觉挺有意思的

library('ImmuCellAI')

f <-function(x) unlist(strsplit(x['gene'],'[|]'))[2]
exp$gene<-rownames(exp)
exp$gene<- apply(exp,1,f)

match<- bitr(exp$gene, fromType = "UNIPROT",toType = c("SYMBOL"),OrgDb = org.Hs.eg.db)
match<-match[match(exp$gene,match$UNIPROT),]
exp$gene<-match$SYMBOL

exp<-exp[!is.na(exp$gene),]
exp<-exp[!duplicated(exp$gene),]

rownames(exp)<-exp$gene
exp<-exp[,-193]
exp<-log2(exp+1)

ImmuCellAI<-ImmuCellAI_new(sample=exp,data_type ='microarray',
                           group_tag=0,response_tag=0)  #我这里用的是蛋白质组的数据，这里把uniprot的ID转化为了SYMBOL_ID
                                                        #作者没说过可以用来分析蛋白质组学数据，但我感觉蛋白组和微阵列好像也没什么差别？
ImmuCellAIscore<-ImmuCellAI[["Sample_abundance"]]

write.csv(ImmuCellAIscore,file = '~/biodata/CNHPP_HCC/MS/EXP_ImmuCellAIscore.csv')

immalscore<-read.csv('EXP_ImmuCellAIscore.csv')
metanmf$adj_sample_2<-paste0('iBAQ.',metanmf$Case.No.,'_T')
immalscore<-immalscore[match(metanmf$adj_sample_2,immalscore$X),]

med.exp <-median(metanmf$InfiltrationScore)
metanmf$Infiltration_status<-ifelse(metanmf$InfiltrationScore>med.exp,'hot','cold')

test<-survfit(Surv(Total.follow.up.period..m.,Died.of.recurrence) ~ Infiltration_status, data = metanmf)
ggsurvplot(test,data=metanmf,palette="Pastel1",
           pval = TRUE,pval.method = TRUE,
           xlab = 'Time (Month)', ylab="surv probability",
           ggtheme = theme_survminer(),
           surv.median.line = 'hv')
           
#从生存分析结果上上来看，ImmuCellAI确实给出了很好的划分，hot，cold两组的预后有着显著的差异，
