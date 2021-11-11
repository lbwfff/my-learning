ms<-read.table('proteinGroups.txt',sep = '\t',header = T)
ms<-ms[-grep('CON__',ms$Protein.IDs),]
ms<-ms[-grep('_HUMAN',ms$Protein.IDs),]
ms<-ms[-grep('REV__',ms$Protein.IDs),]


library(stringr)
test<-as.data.frame(str_split(ms$Protein.IDs,"[;]",simplify = T))

test[test == ""] <- NA
substr(test[1,],1,15)
ms<-ms[apply(test, 1, function(x) length(unique(substr(x,1,15))) == 2), ]

colnames(ms)
exp<-ms[,19:29]
exp$codingRNA<-substr(ms$Protein.IDs,1,15)
# rownames(exp)<-test2$codingRNA

library(dplyr)
test2<-exp %>% group_by(codingRNA) %>% mutate(codingRNA = paste(codingRNA, row_number(), sep="_"))#可以对coding基因进行计数的一行代码

exp<-exp[,-12]
colnames(exp)<-c('POOL','T131','P132','T135','P136','T137','P138','T141','P142','T145','P146')

rownames(exp)<-test2$codingRNA
exp<-exp[,-1]
exp[exp == 0] <- NA
exp<-na.omit(exp)
exp<-log2(exp)

group_list<-c('Tumor','normal','Tumor','normal','Tumor','normal','Tumor','normal','Tumor','normal')
group_list=as.factor(group_list)
dat = exp
group1 = which(group_list == levels(group_list)[1])
group2 = which(group_list == levels(group_list)[2])
dat1 = dat[, group1]
dat2 = dat[, group2]

test<-t.test(as.numeric(exp[1,])~group_list)$p.value

pvals = apply(exp, 1, function(x){t.test(as.numeric(x)~group_list)$p.value})
p.adj = p.adjust(pvals, method = "BH")
avg_1 = rowMeans(dat1)
avg_2 = rowMeans(dat2)
log2FC = avg_1-avg_2
DEG_t.test = cbind(avg_1, avg_2, log2FC, pvals, p.adj)
DEG_t.test=DEG_t.test[order(DEG_t.test[,4]),]
DEG_t.test=as.data.frame(DEG_t.test)
head(DEG_t.test)
