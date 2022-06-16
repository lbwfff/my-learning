########################不怎么熟悉这个算法，但挺有意思的###############

install.packages('NMF')

library(NMF)

nmfprot<-allprot[,-1]
rownames(nmfprot)<-allprot$proteins
nmfprot<-nmfprot[,grep('_T',colnames(nmfprot))]

nmfprotnull<-rowSums(nmfprot)
nmfprot<-nmfprot[nmfprotnull>0,] #要做这个算法的话，rowSums不能为0
nmfprot<-log2(nmfprot+1)

m.mad <- apply(nmfprot,1,mad)
nmfprot <- nmfprot[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.3)),] #25%  = 0.22 #筛掉一些绝对中位差低的数据

testnmf<-nmf(nmfprot,2:8,nrun=50) #看什么样的分组是较好的
plot(testnmf)

nmf.rank5 <- nmf(nmfprot, rank = 2, nrun=50) #根据上面的结果再做分析

jco <- c("#2874C5","#EABF00","#C6524A","#868686") #,'grey'
index <- extractFeatures(nmf.rank5,"max") 
sig.order <- unlist(index)

NMF.Exp.rank5 <- nmfprot[sig.order,]
NMF.Exp.rank5 <- na.omit(NMF.Exp.rank5)

group <- predict(nmf.rank5) # 提出亚型
table(group)

pdf("nmf.pdf",width = 8,height = 8)
par(mar=c(8, 8.5, 3, 3))
consensusmap(nmf.rank5,labCol = NA,
             annCol = NA, #data.frame("cluster"=group[colnames(NMF.Exp.rank5)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4]))) #分组的热图，这个热图一直没搞明白怎么导出才好
dev.off()

basismap(nmf.rank5,cexCol = 1.5,cexRow = 1,
         annColors=list(c("#2874C5","#EABF00","#C6524A","#868686")))

metahcc$adj_sample<-c(paste0('Intensity.',metahcc$Case.No.,'_T'))
metanmf<-metahcc[match(names(group),metahcc$adj_sample),]
identical(names(group),metanmf$adj_sample)

metanmf$group = group
library(survival)
library(survminer)
sfit <- survfit(Surv(Disease.free.survival..m.,Cancer.recurrence) ~ group,  #后面直接一个生存分析，还挺方便的
                data = metanmf)
ggsurvplot(sfit,pval = T,palette = "jco")

