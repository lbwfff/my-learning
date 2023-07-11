################免疫浸润和免疫检查点############################

#如果想用Rmarkdown的话，设置路径会不太一样
setwd('第14节课/01_Immune infiltration/CIBERSORT/')

library(RColorBrewer)
library(tidyr)
library(radiant.data)
library(aplot)
library(cowplot)
library(dplyr)
library(e1071)

#CIBERSORT,一个免疫浸润的反卷积工具
source('Cibersort_source.R')
result <- CIBERSORT('LM22.txt','STAD_normal_tumor.txt', perm = 100, QN = F)  
#太慢了，这玩意应该可以并行的吧。好吧，又是一个不能再windows上并行的函数。

Group <- read.table("./risk.train.group.txt",sep="\t",row.names = 1,header = T)  #预后模型得到的高低风险分组

#合并分组信息
result1 = as.data.frame(result)
result1 <- result1[-c(1:32),]  #去除正常样本
row.names(result1) <- substr(row.names(result1),1,12)  #去除最后4个字符，为了和预后模型分组结果的ID保持一致
#根据行名进行合并
result1 <- merge(Group,result1,by="row.names")
#去除没用的信息
result1 = result1[,-c(2:6)]
row.names(result1) <- result1[,1]
result1 <- result1[,-1]

result1[1:5,1:5]

#根据风险组排序
re1 <- result1[order(result1[,1],decreasing = T),]
colnames(re1)[1] <- "Type"
#查看高低风险组个有几个
table(re1$Type)  #high:72; low:49
#去除统计检验不显著的样本
re1 = subset(re1,re1$`P-value`<0.05)  
#这个P值是置换检验的P值的话，那他的意义就会类似于GSEA？

re1 = re1[,-(24:26)] #去除最后3个指标:p,corrrelation,RESM
re2 = re1[,-1]       #去除第一列分组信息
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
#提取数据，多行变成多列，要多学习‘tidyr’里面的三个函数
dat_cell <- re2 %>% as.data.frame() %>%rownames_to_column("Sample") %>%gather(key = Cell_type,value = Proportion,-Sample)
#提取数据
dat_group = gather(re1,Cell_type,Proportion,-Type )
#合并分组
dat = cbind(dat_cell,dat_group$Type)
colnames(dat)[4] <- "Type"

##2.2柱状图##############
p1 <- ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(28))
p1

#画分组bar
Var1 = re1$Type #high risk & low risk的数量
Var2=c(1:length(row.names(re1))) #high risk & low risk长度，从1开始
annotation_col<- data.frame(Var1 = re1$Type,
                            Var2=c(1:length(row.names(re1))),
                            value=1)
p2 <- ggplot(dat, aes(Var2, Var1),size=0.5)+
  geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
  labs(x = "", y = "")+
  theme_void()+
  theme(legend.position = "top", legend.title = element_blank())+
  scale_x_discrete(name="")+
  scale_fill_manual(labels=c("High risk group(n=72)","Low risk group(n=49)"),
                    values =c("#DC0000FF","#00A087FF"))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
p2

#拼图
pdf('ssGSEAhistogram_ImmuneGeme.pdf',height = 6,width = 10)
plot_grid(p2,p1,
          rel_heights = c(0.08,1),
          label_x=0,
          label_y=1,
          align = 'v',ncol = 1,greedy = F)
dev.off()  #没看出什么规律来，有点太乱了

##2.3 免疫细胞热图#############
data <- as.data.frame(t(re1))
data <- data[-1,]
data_1 <- as.data.frame(lapply(data,as.numeric))
row.names(data_1) <- row.names(data)
data=data_1

library(pheatmap)
pdf(file="immunecell_pheatmap.pdf",width=10,height=8,onefile=FALSE)
annotation_col=data.frame(Group=rep(c("Low risk group(n=49)","High risk group(n=72)"),c(49,72)))
annColors <- list(Group = c("Low risk group(n=49)" = "#00A087FF","High risk group(n=72)" ="#DC0000FF"))
rownames(annotation_col)=colnames(data)
pheatmap(data,
         annotation_col = annotation_col,
         annotation_colors = annColors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cluster_rows=TRUE,
         show_rownames=TRUE,
         fontsize_row=10,
         show_colnames=FALSE,
         scale="row",
         cluster_cols=FALSE,
         main="")
dev.off()  #这个热图没拼起来不知道为什么，所以我不喜欢用pheatmap

##2.4组免疫细胞的差异###################
library(ggpubr)
library(ggplot2)
box=dat

#图片美化
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小     
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗 
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position=c(0.9, 0.8),#图例在绘图区域的位置
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
} #自己包装了一个美化函数，还挺有意思的

e1 <- ggplot(box,aes(x=Cell_type,y=Proportion),palette = "jco", add = "jitter")+
  geom_boxplot(aes(fill=Type),position=position_dodge(0.5),width=0.6)+
  labs(x = "Cell Type", y = "Estimated Proportion")+
  scale_fill_manual(values = c("#DC0000FF","#00A087FF")) +
  theme_zg()
e1 = e1 + stat_compare_means(aes(group = Type),label = "p.signif")
e1
#保存图片
ggsave('box_ImmuneGeme.pdf', plot = e1,width=15,height = 6)

#这张要直观太多了

############################################################
#然后在这里我试了一下以前用过的immuceliai看能不能也用来做这个工作，结果也是可以的

library('ImmuCellAI')

setwd('第14节课/01_Immune infiltration/CIBERSORT/')
dataset<-data.table::fread('STAD_normal_tumor.txt',sep = '\t')
dataset<-as.data.frame(dataset)
row.names(dataset) <- dataset$V1
dataset <- dataset[,-1] #这数据看起来是已经归一化的,要改成未归一化的模式？
dataset<-((2^dataset)-1)

ImmuCellAI<-ImmuCellAI_new(sample=dataset,data_type ='rnaseq',
                           group_tag=0,response_tag=0)
#速度要快多了

ImmuCellAIscore<-ImmuCellAI[["Sample_abundance"]]

#然后炮制一下可视化的代码

Group <- read.table("./risk.train.group.txt",sep="\t",row.names = 1,header = T)  #预后模型得到的高低风险分组

#合并分组信息
result1 = as.data.frame(ImmuCellAIscore)
result1 <- result1[-grep('11A',rownames(result1)),]  #去除正常样本
row.names(result1) <- substr(row.names(result1),1,12)  #去除最后4个字符，为了和预后模型分组结果的ID保持一致
#根据行名进行合并
result1 <- merge(Group,result1,by="row.names")
#去除没用的信息
result1 = result1[,-c(2:6)]
row.names(result1) <- result1[,1]
result1 <- result1[,-1]

result1[1:5,1:5]

#根据风险组排序
re1 <- result1[order(result1[,1],decreasing = T),]
colnames(re1)[1] <- "Type"
#查看高低风险组个有几个
table(re1$Type)  #high:72; low:49

re2 = re1[,-1]       #去除第一列分组信息
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
#提取数据，多行变成多列，要多学习‘tidyr’里面的三个函数
dat_cell <- re2 %>% as.data.frame() %>%rownames_to_column("Sample") %>%gather(key = Cell_type,value = Proportion,-Sample)
#这个写法很有意思，是不是用reshape什么的也可以达到这个效果？

#提取数据
dat_group = gather(re1,Cell_type,Proportion,-Type )
#合并分组
dat = cbind(dat_cell,dat_group$Type)
colnames(dat)[4] <- "Type"

##2.2柱状图##############
p1 <- ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(28))
p1
#画这个图的时候应该把最后一个分数值给去掉，不过只是一个demo而已就懒得改了

#画分组bar
Var1 = re1$Type #high risk & low risk的数量
Var2=c(1:length(row.names(re1))) #high risk & low risk长度，从1开始
annotation_col<- data.frame(Var1 = re1$Type,
                            Var2=c(1:length(row.names(re1))),
                            value=1)
p2 <- ggplot(dat, aes(Var2, Var1),size=0.5)+
  geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
  labs(x = "", y = "")+
  theme_void()+
  theme(legend.position = "top", legend.title = element_blank())+
  scale_x_discrete(name="")+
  scale_fill_manual(labels=c("High risk group(n=72)","Low risk group(n=49)"),
                    values =c("#DC0000FF","#00A087FF"))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
p2

#拼图
pdf('immai_ImmuneGeme.pdf',height = 6,width = 10)
plot_grid(p2,p1,
          rel_heights = c(0.08,1),
          label_x=0,
          label_y=1,
          align = 'v',ncol = 1,greedy = F)
dev.off()  #没看出什么规律来，有点太乱了

##2.3 免疫细胞热图#############
data <- as.data.frame(t(re1))
data <- data[-1,]
data_1 <- as.data.frame(lapply(data,as.numeric))
row.names(data_1) <- row.names(data)
data=data_1

library(pheatmap)
pdf(file="immai_pheatmap.pdf",width=10,height=8,onefile=FALSE)
annotation_col=data.frame(Group=rep(c("Low risk group(n=49)","High risk group(n=72)"),c(49,72)))
annColors <- list(Group = c("Low risk group(n=49)" = "#00A087FF","High risk group(n=72)" ="#DC0000FF"))
rownames(annotation_col)=colnames(data)
pheatmap(data,
         annotation_col = annotation_col,
         annotation_colors = annColors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cluster_rows=TRUE,
         show_rownames=TRUE,
         fontsize_row=10,
         show_colnames=FALSE,
         scale="row",
         cluster_cols=FALSE,
         main="")
dev.off()  #这个热图没拼起来不知道为什么，所以我不喜欢用pheatmap

##2.4组免疫细胞的差异###################
library(ggpubr)
library(ggplot2)
box=dat

#图片美化
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小     
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗 
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position=c(0.9, 0.8),#图例在绘图区域的位置
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
} #自己包装了一个美化函数，还挺有意思的

e1 <- ggplot(box,aes(x=Cell_type,y=Proportion),palette = "jco", add = "jitter")+
  geom_boxplot(aes(fill=Type),position=position_dodge(0.5),width=0.6)+
  labs(x = "Cell Type", y = "Estimated Proportion")+
  scale_fill_manual(values = c("#DC0000FF","#00A087FF")) +
  theme_zg()
e1 = e1 + stat_compare_means(aes(group = Type),label = "p.signif")
e1
#保存图片
ggsave('box_Immai.pdf', plot = e1,width=15,height = 6)

#这玩意最好能排个序，要不然感觉很扎眼。
#好像结果的差别还挺大的，感觉immucellai的结果比较容易有显著变化？不知道有没有什么说法

##########################################################################
#ssGSEA

setwd("第14节课/01_Immune infiltration/sssGSEA/")

library(data.table)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(radiant.data)
library(aplot)
library(cowplot)
library(pheatmap)
library(ggpubr)
library(ggplot2)

library(GSEABase)
library(GSVA)

data <- fread("./STAD_normal_tumor.txt") %>% as.data.frame()
row.names(data) <- data[,1]
data <- data[,-1]

#分组信息
Group <- read.table("./risk.train.group.txt",sep = "\t",header = T,row.names = 1)

#背景免疫细胞
gmtFile="mmc3.gmt"  #这种应该算用基因标志物来做免疫预测，其实会类似于富集分析（GSEA），这里我觉得基因集的选取可能是一个非常关键的地方。

#2 分析免疫细胞含量#########
#载入背景基因
geneSet=getGmt(gmtFile,
               geneIdType=SymbolIdentifier())

##2.1开始ssGSEA分析,数据一定要是矩阵，不是数据框###########
ssgseaScore=gsva(as.matrix(data), geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

range(ssgseaScore)
ssgseaOut <- as.data.frame(t(ssgseaScore))
ssgseaOut <- ssgseaOut[-c(1:32),]  #去除正常样本
row.names(ssgseaOut) <- substr(row.names(ssgseaOut),1,12)  #去除最后4个字符，为了和预后模型分组结果的ID保持一致

ssgseaOut <- merge(Group,ssgseaOut,by="row.names")

## 2.2整理数据格式############
result1 = ssgseaOut[,-c(2:6)]
row.names(result1) <- result1[,1]
result1 <- result1[,-1]
write.table(result1,"ssGSEA_result.txt",sep="\t",row.names = T)
#根据风险组排序
re1 <- result1[order(result1[,1],decreasing = T),]
colnames(re1)[1] <- "Type"
#查看高低风险组个有几个
table(re1$Type)  #high:172; low:173
#去除分组信息
re2 = re1[,-1]       #含有分组信息的数据。
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat_cell <- re2 %>% as.data.frame() %>%rownames_to_column("Sample") %>%gather(key = Cell_type,value = Proportion,-Sample)
#提取数据
dat_group = gather(re1,Cell_type,Proportion,-Type )
#合并分组
dat = cbind(dat_cell,dat_group$Type)
colnames(dat)[4] <- "Type"

##2.3柱状图############
p1 <- ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(28))
p1
#画分组bar
Var1 = re1$Type #high risk & low risk的数量
Var2=c(1:length(row.names(re1))) #high risk & low risk长度，从1开始
annotation_col<- data.frame(Var1 = re1$Type,
                            Var2=c(1:length(row.names(re1))),
                            value=1)
p2 <- ggplot(dat, aes(Var2, Var1),size=0.5)+
  geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
  labs(x = "", y = "")+
  theme_void()+
  theme(legend.position = "top", legend.title = element_blank())+
  scale_x_discrete(name="")+
  scale_fill_manual(labels=c("High risk group(n=172)","Low risk group(n=173)"),
                    values =c("#DC0000FF","#00A087FF"))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
p2
#拼图
pdf('ssGSEAhistogram_ImmuneGeme.pdf',height = 6,width = 10)
plot_grid(p2,p1,
          rel_heights = c(0.08,1),
          label_x=0,
          label_y=1,
          align = 'v',ncol = 1,greedy = F)
dev.off()

##2.4 免疫细胞热图#############
ssgseaOut <- as.data.frame(t(re1))
ssgseaOut <- ssgseaOut[-1,]
ssgseaOut_1 <- as.data.frame(lapply(ssgseaOut,as.numeric))
row.names(ssgseaOut_1) <- row.names(ssgseaOut)
ssgseaOut=ssgseaOut_1

library(pheatmap)
pdf(file="immunecell_pheatmap.pdf",width=10,height=8,onefile=FALSE)
annotation_col=data.frame(Group=rep(c("Low risk group(n=173)","High risk group(n=172)"),c(173,172)))
annColors <- list(Group = c("Low risk group(n=173)" = "#00A087FF","High risk group(n=172)" ="#DC0000FF"))
rownames(annotation_col)=colnames(ssgseaOut)
pheatmap(ssgseaOut,
         annotation_col = annotation_col,
         annotation_colors = annColors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cluster_rows=TRUE,
         show_rownames=TRUE,
         fontsize_row=10,
         show_colnames=FALSE,
         scale="row",
         cluster_cols=FALSE,
         main="")
dev.off() #这次又拼起来了

##2.5组免疫细胞的差异###################
library(ggpubr)
library(ggplot2)
box=dat

#图片美化
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position=c(0.9, 0.8),#图例在绘图区域的位置
          # legend.position="top",
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}

e1 <- ggplot(box,aes(x=Cell_type,y=Proportion),palette = "jco", add = "jitter")+
  geom_boxplot(aes(fill=Type),position=position_dodge(0.5),width=0.6)+
  labs(x = "Cell Type", y = "Immune cell content")+
  scale_fill_manual(values = c("#DC0000FF","#00A087FF")) +
  theme_zg()
e1 = e1 + stat_compare_means(aes(group = Type),label = "p.signif")
e1
#保存图片
ggsave('ssGSEAbox_ImmuneGeme.pdf', plot = e1,width=15,height = 6)

###############################################################################
#
rm(list= ls())

#设置工作路径
setwd("C:/Users/Administrator/Desktop/lesson/第二部分（9-16课次）/课件/2022-07-01课程——免疫/第14节课/02_Immunecheckpoint/")
#加载包

library(ggstatsplot)
library(data.table)
library(ggplot2)
library(ggpubr)

#1 加载数据###########
Group <- read.table("./risk.train.group.txt",sep="\t",row.names = 1,header = T)  #预后模型得到的高低风险分组
data <- fread("./STAD_normal_tumor.txt") %>% as.data.frame()
ImmunecheckpointGene <- read.table("./ImmunecheckpointGene.txt",sep="\t",header = T)  #5个免疫检查点，需要注意别名
#提取和免疫检查点相关基因
dat <- data[data$V1 %in% ImmunecheckpointGene$Gene,]
row.names(dat) <- dat[,1]
dat <- dat[,-1]
table(substr(colnames(dat),14,16)=='11A')
dat <- t(dat[,-c(1:32)])
row.names(dat) <- substr(row.names(dat),1,12)  #去除最后4个字符，为了和预后模型分组结果的ID保持一致
#合并高低风险分组
dat <- merge(Group,dat,by="row.names")


# 高风险组中的相关性
library(ggplot2)#加载ggplot2包
library(ggpubr)

p_CD274 <- ggscatterstats(data = dat, 
                          x = riskScore,
                          y = CD274,
                          xlab = "riskScore",
                          ylab = "CD274 expression (log2 FPKM)",
                          centrality.para = "mean",
                          margins = "both",
                          xfill = "#CC79A7",
                          yfill = "#009E73",
                          marginal.type = "densigram",
                          marginal = TRUE)
p_PDCD1<- ggscatterstats(data = dat, 
                         x = riskScore,
                         y = PDCD1,
                         xlab = "riskScore",
                         ylab = "PDCD1 expression (log2 FPKM)",
                         centrality.para = "mean",
                         margins = "both",
                         xfill = "#CC79A7",
                         yfill = "#009E73",
                         marginal.type = "densigram",
                         marginal = TRUE)
p_HAVCR2 <- ggscatterstats(data = dat, 
                           x = riskScore,
                           y = HAVCR2,
                           xlab = "riskScore",
                           ylab = "HAVCR2 expression (log2 FPKM)",
                           centrality.para = "mean",
                           margins = "both",
                           xfill = "#CC79A7",
                           yfill = "#009E73",
                           marginal.type = "densigram",
                           marginal = TRUE)
p_CTLA4 <- ggscatterstats(data = dat, 
                          x = riskScore,
                          y = CTLA4,
                          xlab = "riskScore",
                          ylab = "CTLA4 expression (log2 FPKM)",
                          centrality.para = "mean",
                          margins = "both",
                          xfill = "#CC79A7",
                          yfill = "#009E73",
                          marginal.type = "densigram",
                          marginal = TRUE)
p_LAG3 <- ggscatterstats(data = dat, 
                         x = riskScore,
                         y = LAG3,
                         xlab = "riskScore",
                         ylab = "LAG3 expression (log2 FPKM)",
                         centrality.para = "mean",
                         margins = "both",
                         xfill = "#CC79A7",
                         yfill = "#009E73",
                         marginal.type = "densigram",
                         marginal = TRUE)
#排版
pdf("immunecheckpoint_cor.pdf",width = 10,height = 10)
ggarrange(p_CD274, p_CTLA4, p_HAVCR2,p_LAG3,p_PDCD1+
            rremove("x.text"),  ncol = 2, nrow = 3) #一个拼图的工具
dev.off()

#这种就是自己划定一个基因集，然后去看相关性呗，还挺好理解的















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
