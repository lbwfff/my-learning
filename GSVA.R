## 环境设置
### 路径
setwd("02_GSVA/")

### 加载包
library(dplyr)
library(data.table)
library(GSVA)
library(GSEABase)
library(limma)
library(stringr)
library(ggplot2)
library(pheatmap)

options(stringsAsFactors = FALSE) #禁止chr转成factor

# 1 输入文件#############
gsym.expr <- fread("STAD_tumor_orderByRiskGroup.txt")
gsym.expr <- as.data.frame(gsym.expr);row.names(gsym.expr) <- gsym.expr[,1];gsym.expr <- gsym.expr[,-1]
head(gsym.expr)
# 要分析的通路
geneset <- getGmt("c2.cp.kegg.v7.4.symbols.gmt") #读gmt文件

# GSVA
gsva_es <- gsva(as.matrix(gsym.expr), geneset)

# 预览GSVA分析返回的矩阵
head(gsva_es)
# 把通路的表达量保存到文件
# write.csv(gsva_es, "gsva_output.csv", quote = F)

# 2 通路的差异表达分析#############
# 分组
group_list <- read.table("risk.train.group.txt",header = T)
group_list <- arrange(group_list,desc(group)) #对样本进行排序。默认是升序排列，这里需要检查
# GSVA排序
gsva_es <- gsva_es[,group_list$id]

## 差异分析
table(group_list$group)
Control_num<-173  #对照样本数
Case_num<-172  #处理样本数

group<-factor(c(rep("Control",Control_num),rep("Case",Case_num)),levels=c('Control','Case'))  #设置分组文件
design<-model.matrix(~0+group)
colnames(design)<-c("Control","Case")
fit<-lmFit(gsva_es,design)
contrast.matrix <- makeContrasts(Case-Control,levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
options(digits = 4)#输出结果的小数点后保留4位

result_pathways<-topTable(fit2,coef=1,adjust="BH",number=dim(gsva_es)[1])  #对于GSVA得分的差异分析居然是用limma来做的，为什么，有什么理由这样会更好吗？

#所有基因的差异情况
# write.csv(result_pathways,file="DEpathways_information.csv", quote = F, row.names = T)

#输出t值，用做作图
pathway <- str_replace(row.names(result_pathways), "KEGG_", "") #把KEGG的前缀去除了
df <- data.frame(ID = pathway, score = result_pathways$t)
# write.csv(df, "DEpathways_tscore.csv", quote = F, row.names = F)


# 3 开始画图-t值的柱状图#####################
df <- read.csv("DEpathways_tscore.csv")
head(df)

#按照score的值分组
cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "inward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "outward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  xlab("") +ylab("t value of GSVA score, tumor \n versus non-malignant")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴

ggsave("gsva.pdf", width = 8, height = 8)
#很有意思这种美化方法，可以参考

# 4 开始画图-t值的柱状图#################
diff_pathways <- subset(result_pathways,adj.P.Val<0.05)
red_de_expr <- gsva_es[rownames(gsva_es) %in% rownames(diff_pathways),] #注意这里是data_TCGA,不是data
# write.table(red_de_expr,file="diffSigPathways.txt",sep="\t",quote=F)

#然后就是热图，之前也做过的，没什么好说的
table(group_list$group)
annotation_col=data.frame(clinical=rep(c("low(n=173)","high(n=172)"),c(Control_num,Case_num)))
annColors <- list(clinical = c("low(n=173)" = "#00A087FF","high(n=172)" ="#DC0000FF"))
rownames(annotation_col)=colnames(gsva_es)
pheatmap(red_de_expr, 
         annotation_col = annotation_col,
         annotation_colors = annColors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = FALSE,
         clustering_method = "complete",  #clustering_method参数设定不同聚类方法，默认为"complete",可以设定为'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
         clustering_distance_rows = "euclidean",  #参数设定行聚类距离方法为Pearson corralation，默认为欧氏距离"euclidean"
         cluster_rows=TRUE,
         cluster_cols=FALSE,
         treeheight_row = 50,  #treeheight_row和treeheight_col参数设定行和列聚类树的高度，默认为50
         treeheight_col = 50,
         show_rownames=TRUE,
         show_colnames=FALSE,
         display_numbers = FALSE,  #display_numbers = TRUE参数设定在每个热图格子中显示相应的数值(默认F)，也可以自定义显示某一类
         number_color = "blue",  #number_color参数设置数值字体的颜色
         fontsize_row=8,
         scale="row",  #对行进行均一化
         legend = TRUE,
         main="")
# dev.off()


