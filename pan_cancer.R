#我是真没做过泛癌的分析，这些代码都是我瞎写的，有些语句写得很笨重，仅供参考吧。

meta<-read.table('~/biodata/pancancer/TCGA_phenotype_denseDataOnlyDownload.tsv',sep = '\t',header = T)#数据来自UCSC Xena
type<-as.data.frame(table(meta$X_primary_disease))

normal<-meta[meta$sample_type_id=='11',] #我只想要有非癌配对的类型，类型11指的是正常组织样本
test<-as.data.frame(table(normal$X_primary_disease))
test<-test[test$Freq>25,] #不想要正常组织样本太少的类型

normal<-normal[normal$X_primary_disease %in% test$Var1,]
meta<-meta[meta$X_primary_disease %in% test$Var1,]
meta<-meta[!is.na(meta$sample_type_id),]
#给剩下了12种，要给整在在一起

tumor<-meta[meta$sample_type_id<11,]
normal<-meta[meta$sample_type_id==11,]
tumor$Group<-paste0(tumor$X_primary_disease,"_T")
normal$Group<-paste0(normal$X_primary_disease,'_P')
  
meta<-rbind(normal,tumor)             
inf<-as.data.frame(table(meta$Group))
group<-c(rep(inf$Var1,c(inf$Freq)))

exp<-read.table('~/biodata/pancancer/tcga_RSEM_gene_fpkm',sep = '\t',header = T)
exp[1:5,1:5]
exp<-c(1)

exp$sample<-substr(exp$sample,1,15)
ncRNA<-read.csv('~/biodata/pancancer/coding_deg.csv')
exp_2<-exp[exp$sample %in% substr(ncRNA$major_orf,1,15),]

f <-function(x) paste0(unlist(strsplit(x['sample'],'-'))[1],'.',unlist(strsplit(x['sample'],'-'))[2],'.',unlist(strsplit(x['sample'],'-'))[3],'.',unlist(strsplit(x['sample'],'-'))[4])
meta$adjust_name<-apply(meta,1,f)
rownames(exp_2)<-exp_2$sample
exp_2<-exp_2[,colnames(exp_2) %in% meta$adjust_name]
meta<-meta[meta$adjust_name %in% colnames(exp_2),]
table(meta$Group)

exp_2<-exp_2[,match(meta$adjust_name,colnames(exp_2))]
exp_2[1:5,1:5]

#ENSG00000225492,ENSG00000259479,ENSG00000255353相关性较好的就这几个，把三个都做一些？

GBP1P1<-exp_2[rownames(exp_2)=='ENSG00000225492',]
GBP1P1<-as.data.frame(t(GBP1P1))
GBP1P1$Group<-meta$Group

SORD2P<-exp_2[rownames(exp_2)=='ENSG00000259479',]
SORD2P<-as.data.frame(t(SORD2P))
SORD2P$Group<-meta$Group

ASS1P13<-exp_2[rownames(exp_2)=='ENSG00000255353',]
ASS1P13<-as.data.frame(t(ASS1P13))
ASS1P13$Group<-meta$Group

#绘图，没做过这种，感觉要调整挺久的
library('ggplot2')
library('MetBrewer')

col<-met.brewer("Hokusai1", 7)[c(6,4)] #这个颜色是真的丑，我也不知道老板喜欢它哪里
col<-rep(col,12)

dp <- ggplot(ASS1P13, aes(x=Group, y=ENSG00000255353))+
  geom_boxplot(width=0.5,aes(fill=Group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=Group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  theme(legend.position = "none")+
  labs(y="ASS1P13 log2(FPKM)")

pdf("ASS1P13_TCGA.pdf",width = 8,height = 6)
dp + theme_classic(base_size = 22)+ 
  #scale_colour_manual(values=met.brewer("Hokusai2", 3))+
  scale_fill_manual(values=col)+
  #scale_size_manual(values=met.brewer("Hokusai2", 2))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=22))+
  theme(legend.position = "none")+
  theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))+
  theme(axis.text.x = element_text(size=12,angle=45,hjust=0.9))

dev.off()

dp <- ggplot(GBP1P1, aes(x=Group, y=ENSG00000225492))+
  geom_boxplot(width=0.5,aes(fill=Group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=Group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  theme(legend.position = "none")+
  labs(y="GBP1P1 log2(FPKM)")

pdf("GBP1P1_TCGA.pdf",width = 8,height = 6)
dp + theme_classic(base_size = 22)+ 
  #scale_colour_manual(values=met.brewer("Hokusai2", 3))+
  scale_fill_manual(values=col)+
  #scale_size_manual(values=met.brewer("Hokusai2", 2))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=22))+
  theme(legend.position = "none")+
  theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))+
  theme(axis.text.x = element_text(size=12,angle=45,hjust=0.9))

dev.off()

dp <- ggplot(SORD2P, aes(x=Group, y=ENSG00000259479))+
  geom_boxplot(width=0.5,aes(fill=Group),colour='black',alpha = 0.4,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=Group),color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.2,color="black")+
  theme(legend.position = "none")+
  labs(y="SORD2P log2(FPKM)")

pdf("SORD2P_TCGA.pdf",width = 8,height = 6)
dp + theme_classic(base_size = 22)+ 
  #scale_colour_manual(values=met.brewer("Hokusai2", 3))+
  scale_fill_manual(values=col)+
  #scale_size_manual(values=met.brewer("Hokusai2", 2))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=22))+
  theme(legend.position = "none")+
  theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))+
  theme(axis.text.x = element_text(size=12,angle=45,hjust=0.9))

dev.off()#我喜欢这个结果

#有一个改注释的过程，想把长串的癌症名改成缩写，但这个过程的代码写得很呆，我暂时也想不到有什么能够简洁的做到的代码（或许换一个meta文件？）
meta$Group<-gsub("thyroid carcinoma_P","THCA_N",meta$Group)
meta$Group<-gsub("thyroid carcinoma_T","THCA_T",meta$Group)
meta$Group<-gsub("stomach adenocarcinoma_P","STAD_N",meta$Group)
meta$Group<-gsub("stomach adenocarcinoma_T","STAD_T",meta$Group)
meta$Group<-gsub("prostate adenocarcinoma_P","PRAD_N",meta$Group)
meta$Group<-gsub("prostate adenocarcinoma_T","PRAD_T",meta$Group)
meta$Group<-gsub("lung squamous cell carcinoma_P","LUSC_N",meta$Group)
meta$Group<-gsub("lung squamous cell carcinoma_T","LUSC_T",meta$Group)
meta$Group<-gsub("lung adenocarcinoma_P","LUAD_N",meta$Group)
meta$Group<-gsub("lung adenocarcinoma_T","LUAD_T",meta$Group)
meta$Group<-gsub("liver hepatocellular carcinoma_P","LIHC_N",meta$Group)
meta$Group<-gsub("liver hepatocellular carcinoma_T","LIHC_T",meta$Group)
meta$Group<-gsub("kidney papillary cell carcinoma_P","KIRP_N",meta$Group)
meta$Group<-gsub("kidney papillary cell carcinoma_T","KIRP_T",meta$Group)
meta$Group<-gsub("kidney clear cell carcinoma_P","KIRC_N",meta$Group)
meta$Group<-gsub("kidney clear cell carcinoma_T","KIRC_T",meta$Group)
meta$Group<-gsub("head & neck squamous cell carcinoma_P","HNSC_N",meta$Group)
meta$Group<-gsub("head & neck squamous cell carcinoma_T","HNSC_T",meta$Group)
meta$Group<-gsub("uterine corpus endometrioid carcinoma_P","UCEC_N",meta$Group)
meta$Group<-gsub("uterine corpus endometrioid carcinoma_T","UCEC_T",meta$Group)
meta$Group<-gsub("colon adenocarcinoma_P","COAD_N",meta$Group)
meta$Group<-gsub("colon adenocarcinoma_T","COAD_T",meta$Group)
meta$Group<-gsub("breast invasive carcinoma_P","BRCA_N",meta$Group)
meta$Group<-gsub("breast invasive carcinoma_T","BRCA_T",meta$Group)
meta$Group<-gsub("breast invasive carcinoma_P","BRCA_N",meta$Group)
meta$Group<-gsub("breast invasive carcinoma_T","BRCA_T",meta$Group)

