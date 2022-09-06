library(ggplot2)
library(ggsignif)
library(ggsci)
library(ggpubr)

########################################
par(mar=c(8, 8.5, 3, 3))  #对页面布局非常重要的一句，数值分别代表bottom, left, top, right。
                          #一个具体的例子就是，如果底部的字跑出去了，就可以把8改大，可以把跑出去的字画回图里。
                          
#######################################
#ggplot想要固定长宽比的话可以：
 theme(aspect.ratio=1) #比率可以自己调

############1##############
ggplot(data = cor,aes(x=METTL3,y=CAPRIN1))+
  geom_point(shape=21,size=4)+
  stat_smooth(method = lm)+
  stat_cor(data=cor,method='pearson')+
  theme_classic()+
  ggsave(file='METTL3-CAPRIN1-GBM.PDF',width = 10,height = 8)
 #一个简单的散点图加回归线，回归用的简单线性回归分析，平滑线附近的蓝色区域为置信区间

##########1+加强版#############
ggplot(data = TPCOR,aes(x=METTL3,y=TRA2A))+#把癌症组和非癌症组分开统计了
  geom_point(shape=21,size=2)+#把点改为了空心，并放大了一倍，美化了不少
  stat_smooth(method = lm)+
  stat_cor(data=TPCOR,method='pearson')+
  theme_classic()+
  labs(title="TP",subtitle = 'sample size=155')+
  ggsave(file='TP-METTL3-TRA2A-GBM.PDF',width = 10,height = 8)
#有一个问题始终无法解决，就是p值小于2.2e-16后不能显示真实值的问题，现在能找到的办法就是在出图后用pdf编辑器改，p值获得方法如下：
A=TPCOR$MOV10
B=TPCOR$METTL3
C<-cor.test(A,B,method='pearson')
p <- 2* pt(C$statistic,  df = C$parameter, lower.tail=FALSE)#p值在这里得到
#还有一个细节就是平滑线的颜色，可以再考虑一下，出多张图的时候可以使用不同颜色的平滑线。

#############再次加强版##########
C<-cor.test(cor$TRA2A,cor$METTL3,method='pearson')
D<-cor(cor$TRA2A,cor$METTL3,method='pearson')
ptext<-paste0("R=",signif(D,2),',','P=',signif(2* pt(C$statistic, df = C$parameter, lower.tail=FALSE),3))#手动计算了R值和p值
ggplot(data = cor,aes(x=METTL3,y=TRA2A))+
  geom_point(shape=21,size=4)+#把点再加大了一倍
  stat_smooth(method = lm)+
  #stat_cor(data=cor,method='pearson',size=6,p.digits = digits)+
  theme_classic()+
  annotate('text',x=0.7,y=2.1,label='Sample N=109',colour = "black",size=6)+
  annotate('text',x=0.7,y=2.2,label=ptext,colour = "black",size=6)+
  ggsave(file='TOTAL-cor-CAPTAC.PDF',width = 10,height = 8)

#############密度散点图，在上面那一段的基础上改的################
library(ggpointdensity)
library(viridis)
pdf("SODR2P_cor.pdf",width = 8,height = 6)
ggplot(data = cortest,aes(x=log2(cortest[,1]),y=cortest[,2]))+
  geom_pointdensity(adjust = 4,show.legend =F)+
  scale_color_viridis()+ #点密度的颜色，暂时不知道怎么自定义这个颜色
  stat_smooth(method = lm)+
  theme_classic()+
  ylab('Peptide expression')+
  xlab('RNA expression')+
  coord_cartesian(xlim=c(5,12))+
  annotate('text',x=6.4,y=2.2,label=ptext,colour = "black",size=6) #
dev.off()

##########2###############
forcor<-t(LUADIWANT)
colnames(forcor)=c('MOV10','TRA2A','CAPRIN1','METTL17')
library(corrplot)
M <- cor(forcor)
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue"))
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
wb <- c("white", "black")
corrplot(M, method = "number")
#一个多因素的相关性图，先用cor函数算出因数的相关性矩阵，后进行绘制，能绘制出的种类还挺多的，多因素相关分析可以用

##########3##############
ggplot(try[try$`E-ID`=='MOV10',],aes(group,value,fill=group))+geom_boxplot(width=0.5,outlier.shape = NA)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(title="MOV10",x="group", y = "TPM")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  geom_line(aes(group=from) ,size=0.8,colour="gray",alpha = 0.07)+
  scale_fill_npg()+ 
  geom_jitter(aes(fill=group),width =0.2,shape = 21,size=1,colour="black",alpha = 0.04)+
  coord_cartesian(ylim=c(0,50))+
  stat_compare_means(method = "t.test", label.x = 0.7, label.y = 45)+
  scale_fill_brewer(palette="RdBu")+
  ggsave( file='GBM-MOV10-BOXPLOT.pdf', width=10, height=8)

#改了很多版才出的一个箱线图，其实我也还没有太摸清楚要怎么去做。加了连线，加了点点，限制了坐标轴，加了p值，改了颜色，这个配色我觉得会soft一些，但也不知道老板觉得怎么样。

############4###############
dp <- ggplot(try[try$`E-ID`=='MOV10',], aes(x=group, y=value, fill=group))+
  geom_violin(trim=FALSE)+
  scale_fill_npg()+
  geom_line(aes(group=from) ,size=0.8,colour="gray",alpha = 0.07)+
  #geom_jitter(aes(fill=group),width =0.2,shape = 21,size=1,colour="black",alpha = 0.1)+
  geom_boxplot(width=0.1, fill="white",outlier.shape = NA)+
  labs(title="MOV10",x="group", y = "TPM")+
  coord_cartesian(ylim=c(0,80))+
  stat_compare_means(method = "t.test", label.x = 0.7, label.y = 78)
dp + theme_classic()+ scale_fill_brewer(palette="RdBu") + 
  theme_classic()+
  ggsave( file='GBM-MOV10-VIOLINPLOT.pdf', width=10, height=8)
#小提琴图，我一直没想清楚到底要不要点点，所有带点或不带点都做了一版，小提琴里加了箱线，配色延续的上面箱线图的配色，自认为比较soft，但是谁知道呢？
#然后还是存在一些问题，就是颜色和p值，颜色方面如果可以将肿瘤组和正常组进行交换就更好了，但目前为止调色盘什么的要怎么使用我还得自己摸索
#p值的问题我或许可以根据上面线性相关的思路解决，也不好每次都放到excel里面求，但其实两种方法都挺麻烦的，试着去寻找一些更加彻底的方法吧。


##################5#########################还是小提琴图，我也不知道怎么样才是好的，又或许抄的作业才是最好的
dp <- ggplot(mat, aes(x=Grade, y=log2(TRA2A), fill=Grade))+
  geom_violin(trim=FALSE,show.legend=FALSE)+
  geom_jitter(aes(fill=Grade),width =0.2,shape = 21,size=0.6,colour="black",alpha = 0.07)+
  geom_boxplot(width=0.1, fill="white",outlier.shape = NA)+
  labs(y="TRA2A  log2(TPM)",subtitle= 'WHO II N=103
WHO III N=79
WHO IV N=139')+
  #theme(axis.title.x = element_text(size = 15, family = "myFont", color = "green", face = "bold", vjust = 0.5, hjust = 0.5, angle = 45))+
  stat_compare_means(method = "anova",label.x = 0.7,size=5)
dp + theme_classic()+ scale_fill_brewer(palette="Blues") + 
  theme(legend.position = "none")+
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=16))+
  ggsave( file='CGGA325-GRADE-TRA2A.pdf', width=8, height=7)
 
#############终极小提琴###################
cols<-c('#36537155', "#96345A74")#从他人图上选取的颜色，效果还不错所以一直在用
pal<-colorRampPalette(cols)
image(x=1:20,y=1,z=as.matrix(1:20),col=pal(20))
image(x=1:20,y=1,z=as.matrix(1:20),col=pal(4))
image(x=1:20,y=1,z=as.matrix(1:20),col=pal(2))

colors  = pal(4)#根据数据的分组情况选择渐变吧

library(ggpubr)
ptext<-paste0("Anova, p=",compare_means(value~class,exp,method = "anova")[[3]])
dat_text <- data.frame(label = ptext)#手动算出p值
dp <- ggplot(exp, aes(x=class, y=log2(value)))+
               geom_violin(aes(colour=class),fill='#DDDDDD50',trim=FALSE,show.legend=FALSE,size=0.7)+
               geom_jitter(aes(fill=class),width =0.2,shape = 21,size=4,colour='NA',alpha = 0.4)+
               geom_boxplot(fill='white',width=0.12,,outlier.shape = NA,size=0.8)+
               labs(y="TRA2A  log2(TPM)")
             #stat_compare_means(label.x = 0.7,size=5)#这一句还是保留了
dp + theme_classic(base_size = 22)+ 
               scale_colour_manual(values=c(colors))+
               scale_fill_manual(values=c(colors))+
               theme(legend.position = "none")+
               theme(axis.text=element_text(size=18),
                     axis.title=element_text(size=22))+
               annotate('text',x=1.1,y=8,label=dat_text,colour = "black",size=6)+
               annotate('text',x=1.1,y=9,label='Normal N=20
WHO II N=103
WHO III N=79
WHO IV N=139',colour = "black",size=6)+
               ggsave( file='CGGA-TRA2A.pdf', width=8, height=7)
  
  #############自定义的调色盘#########
cols<-c("#A9A9A9", "#ff8C00", "#00BFFF")#自己选择的颜色
pal<-colorRampPalette(cols)
image(x=1:20,y=1,z=as.matrix(1:20),col=pal(20))#审阅不同颜色渐变的划分程度
image(x=1:20,y=1,z=as.matrix(1:20),col=pal(3))
colors = pal(3)#赋值到colors后就可以直接用了，例如以下
p+scale_fill_manual(values=c(colors))

############怎么copy一个颜色############
截图，使用windows画图软件，颜色选取器，编辑颜色，得RGB值
rgb(red=170, green=197, blue=226, alpha=255, max=255)#得颜色编号


#############绘图代码应该多记录，因为一时的懒散而要花费更多的精力显得有点呆######
#############这是一个横向的柱状图，我还挺喜欢的没用来代替一直想画却没有好效果的饼图####
#############但是想用这个画cellline还是很丑######
library(scales)
platform<-platform[order(platform$Freq,decreasing = T),]
platform$Var1<-as.factor(platform$Var1)

pdf("platform-type3.pdf",width = 8,height = 6)
ggplot(data = platform, aes(x = reorder(Var1,-Freq), y = Freq)) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
  geom_col(aes(fill = Var1) , show.legend = FALSE) +
  ggtitle(paste("Platform")) +
  coord_flip() +
  theme_classic()+ ###把背景格抹掉
  geom_label(aes(label = paste(Freq, percent(Freq/sum(Freq)), sep = "\n"),
                 y = -20, fill = Var1),
             show.legend = FALSE,
             size = 3, label.padding = unit(0.1, "lines")) +
  expand_limits(y = -20) +
  labs(y="Frequency",x='Platform')+
  scale_fill_brewer(palette = "Set3", direction = -1) 
dev.off()

############plotly是个很有意思的东西，可以做一些互动交互的图，还挺厉害的######
plot_ly(platform, labels = ~Var1, values = ~Freq, type = 'pie',textposition = 'outside',textinfo = 'label+percent') #####代码并没有画好这一张图，只是想强调一下plotly是个很厉害的东西

ggplotly(p=p) #可以将ggplot的对象转换为plotly的对象，利用plotly的可交互性，可以方便获取许多信息

#############我居然没记录过这个饼图代码################
mis_val<-as.data.frame(table(mstab$mis_val))
percent <- round(mis_val$Freq/sum(mis_val$Freq)*100, 1)
label <- paste(mis_val$Var1,"(", percent, "% )") #, mis_val$Freq

library('MetBrewer')
pdf("mis_val.pdf",width = 9,height = 6)
pie(mis_val$Freq,border="white", col=met.brewer('VanGogh1',n=7,type = 'discrete')[c(1,3,5,7)], label=label) #,main = 'Source of protein (or peptide)' #还没想好标题是什么
dev.off() #这个能用，但是不太好用，还不知道该怎么调（有时候注释全挤到一团了），但我又确确实实用得挺多的

library(tastypie) #一个基于ggplot饼图的便捷工具
pie_bake(data = start, template = "basic4", perc = TRUE,
         title = " ", group_name = "Start condon")+
  scale_fill_manual(values=rev(met.brewer('Hiroshige',n=4)))+
  labs(fill='Start condon')  #说好看吗，其实也没什么区别，但是由于是一个ggplot对象，可能拼图，可定制化什么的会好一些
  

###########PCA图##########
###########我现在喜欢用DEseq包里的函数来画，因为自带标准化，后面也可以用ggplot来美化图片#######
library(DESeq2)
condition <- factor(group_list)
dds <- DESeqDataSetFromMatrix(RPF, DataFrame(condition), design= ~ condition )
dds
rld <- rlogTransformation(dds)
exprSet_new=assay(rld)
p<-plotPCA(rld, intgroup=c('condition'))
p+theme_classic(base_size = 12)+
  scale_colour_manual(values=c(cols))+
  ggsave( file='RPF-PCA.pdf', width=5, height=4)
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)


#########cliProfiler，这个包挺呆的，但是反而画出了想要的效果，emmmmmmm#######
library(cliProfiler)
library(ChIPpeakAnno)
library(ggplot2)

shNC <- read.table('uniqueshNC.bed')
shTRA2A<-read.table('uniqueshTRA2A.bed')
peak<-rbind(shNC,shTRA2A)
write.table(peak, file = "peak.bed", sep = "\t", row.names = F, col.names = F, quote = F)
peak<-toGRanges('peak.bed',format=c('BED'))

peak$sample <- c(rep("shNC",nrow(shNC)), rep("shTRA2A", nrow(shTRA2A)))
meta<-metaGeneProfile(object = peak,annotation= 'gencode.v35.annotation.gff3',
                      group = 'sample',include_intron = TRUE) #####考虑内含子的距离的话，会更加好看不知道为什么
meta[[2]] + ggtitle("Meta Profile")+
  theme(legend.position = "right")+
  theme_classic()+
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))+ ######但优点是高度的可定制化
  scale_y_continuous(expand = c(0,0))   #把横轴从零开始

###########祖传火山图，陪伴我整个生涯的代码，颜色还需要优化一下，其他的地方暂时没什么不满意的#########
library('ggrepel')
Dat$x<-rownames(Dat)
Dat$threshold = factor(ifelse(Dat$adj.P.Val < 0.01 & abs(Dat$logFC) >= 0.75, ifelse(Dat$logFC>= 0.75 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi')) #这里现在差异分析结果里加入了threshold这一列，后续就不用写很长的代码来区分点了。
ggplot(Dat,aes(x=logFC,y=-log10(adj.P.Val),color=threshold))+
  geom_point()+
  geom_point(data = Dat[Dat$x %in% rownames(lncRNA),],color='#664365')+ #这里吗属于lncRNA这个数据里的点给单独用了一种颜色
  scale_color_manual(values=c("#365371","#96345A74","#808080"))+
  #geom_text_repel(data = Dat[Dat$x %in% rownames(lncRNA),],aes(label = substr(x,1,15)),
  #               size = 3,segment.color = "black", show.legend = FALSE )+ #这里可以给属于lncRNA数据的点加上标签，因为前面给了单独的颜色所有没用这行代码
  theme_bw()+
  theme_classic()+ #去掉背景网格
  #theme(legend.title = element_blank())+
  coord_cartesian(xlim = c(-3,3))+
  ylab('-log10 (pvalue)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-0.75,0.75),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.01),lty=3,col="black",lwd=0.5)

############dotplot，我总觉得自己应该抛弃dotplot这个函数，因为感觉每次画出来的一批图，从来没有结构一致的，注释信息到处跑，图的大小从来没有一致过######
library('ggplot2')
library('stringr')
#GO的dot_plot
dotplot(downbp_pannzer,color = "pvalue")+
  scale_color_continuous(low='#1C86EE', high='#BEBEBE')+ #调了好久的颜色
  aes(shape=pvalue > 0.01)+ #小于0.01的点换成三角形
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))#给太长的词条换行
 
 #GO的bar_plot
 p1<-ggplot(data = CC, aes(x = reorder(Term,-FDR), y = -log2(FDR))) + #reorder函数，对元素进行排序
  geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
  ggtitle(paste("GO_CC")) +
  coord_flip() +
  theme_classic(base_size = 26)+ 
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=30,face='bold'))+
  labs(y="-log2(FDR)",x='Term')+
  scale_fill_manual(values = c("#d8443c"))+
  scale_colour_manual(values = c("black"))+
  scale_y_continuous(expand = c(0,0)) 
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
  
#我现在能想到的，让绘图区域大小一致的方法，就是用cowplot把几个图片组合起来，现在倒还是个办法，虽然还是挺麻烦的。
library(cowplot)
pdf("GO_A_VS_B.pdf",width = 12,height = 20)
plot_grid(p2, p1,p3, ncol=1, align="v") 
dev.off() 


################3D的PCA图，具体的美化之后还需要多做调整##################
library('scatterplot3d')
miceimpexp3<-as.matrix(t(miceimpexp2))
miceimpexp_cov <- cov(miceimpexp3)
miceimpexp_cov_eigen <- eigen(miceimpexp_cov)
pc_select = 3
label = paste0("PC",c(1:pc_select))
miceimpexp_new <- miceimpexp3 %*% miceimpexp_cov_eigen$vectors[,1:3]
colnames(miceimpexp_new) <- label
miceimpexp_new<-as.data.frame(miceimpexp_new)
miceimpexp_new$group<-c(rep(1,50),rep(2,50))
miceimpexp_new$group2<-group_list

colors.lib <- c("green3", "#56B4E9")
colors <- colors.lib[as.numeric(miceimpexp_new$group)]
mycolors<-c(rep('magenta',50),rep('green3',50))
shapes.lib = c(16)
shapes <- shapes.lib[as.numeric(miceimpexp_new$group)]


pdf("test.pdf",width = 8,height = 8)
scatterplot3d(miceimpexp_new[,(1:3)],angle=55,color=colors,
              pch=16, main="PCA for test")
legend("topright",legend = levels(miceimpexp_new$group2),col = colors.lib,pch = shapes.lib,inset = -0.08,xpd = T,horiz = F)
dev.off() 

##############自己手动算PCA的话，就感觉很慢，用函数的话就快很多不知道为什么##############
library("FactoMineR")
library("factoextra")

df=as.data.frame(t(datExpr))
dat.pca <- PCA(df, graph = F) #这个算出来感觉和MDS的分布差不多
testp<-fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group_list, # color by groups
             palette = c(met.brewer("Cassatt1")[c(2,6)]),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups") 

testp+theme(panel.grid.minor = element_blank())+
  coord_fixed() + theme_bw() + 
  theme_classic(base_size = 18)+ 
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  xlab(paste0("PC1: ",'7.9',"% variance")) +
  ylab(paste0("PC2: ",'3.3',"% variance")) +  #有几个问题，首先是，散点的大小，形状要怎么调整？其次两轴的变量要怎么提取？
  coord_fixed()

##############2DPCA（散点图）###############
library('ggplot2')
library('MetBrewer')
pdf("MDS-diff.pdf",width = 8,height = 6)
ggplot(mds.data, aes(x=X, y=Y, label=Sample,col=group_list)) + 
  geom_point(size = 6, alpha = 0.8) +
  theme(panel.grid.minor = element_blank()) + #我暂时还没有找到要如何去调整坐标轴的大小的方法
  #coord_fixed() + 
  theme_bw()+
  stat_ellipse(level = 0.95)+ #加分组加椭圆
  scale_color_manual(values=met.brewer("VanGogh1", 3,type ='discrete'))+#新的调色板非常专业
  theme_classic()+
  xlab(paste("Leading logFC dim1 ", mds.var.per[1], "%", sep="")) +
  ylab(paste("Leading logFC dim2 ", mds.var.per[2], "%", sep="")) +
  coord_cartesian(xlim=c(-0.6,0.9),ylim=c(-0.6,0.6))+
  ggtitle("MDS plot using avg(logFC) as the distance")
dev.off()


######################热图#############################
library(pheatmap)#目前用得最习惯的还是pheatmap这个包

p=pheatmap(cortest,cluster_rows = T,cluster_cols = T,
           # annotation_col =annotation_b, annotation_legend=TRUE, #如果需要对纵轴进行注释的话
           scale = "row", #标准化，可以选行，列，或是不进行标准化
           color = color, #颜色，目前用的是color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
           show_rownames = T,show_colnames = F,
           breaks = seq(-1,1,length.out = 100)) #我在想这个length.out是指颜色有一百个刻度吗？可是试着改一下
#关于热图的配色
color<-c(met.brewer("Wissing")[2],met.brewer("Signac")[4],'white',met.brewer("Homer1")[7],met.brewer("Nizami")[8]) #从包里面取了四个色
pal<-colorRampPalette(color)
image(x=1:100,y=1,z=as.matrix(1:100),col=pal(100))
colors = pal(100) #这个配色，也说不上完全满意，高表达的红色总感觉有点老，不过暂时先用着这个吧

##################热图的拼接##############
library('cowplot')
library('ggplotify')
p5<-as.ggplot(p) #cowplot没法直接拼接，需要转化为ggplot对象
p6<-as.ggplot(p2)

pdf("pheatmap_adjust.pdf",width = 6,height = 6)
plot_grid(p5, p6, ncol=2)
dev.off() 

########################################
#相似的，也用这样的方法做过药效曲线的拼接
as.ggplot(~plot(fit1,xlab = "log10 (Dose)", ylab = "Cell growth",ylim = c(-0.1, 1.1),# xlim =c(-5,10),
                        legend_location = "bottomleft",legend = c(paste0('AF'))))
#但问题在于这样做并不能真正得到一个ggplot对象，没法进行美学上的调整，也不能在plot对象上加注释什么的

pdf("test.pdf",width = 12,height = 6)
par(mfrow=c(1,2),pin = c(5,5)) # mfrow后面接行列数，pin后接一张图的长宽比

text1=paste0("IC50=",round(exp(as.data.frame(coef(fit1))[rownames(as.data.frame(coef(fit2)))=='phi',]),3),"nmol/L")

plot(fit1,xlab = "log10 (Dose)", ylab = "Cell growth",ylim = c(-0.1, 1.1),
     legend_location = "bottomleft",legend = c(paste0('KY')))
legend('topright',legend=text1, bty ="n", pch=NA) 

text2=paste0("IC50=",round(exp(as.data.frame(coef(fit2))[rownames(as.data.frame(coef(fit2)))=='phi',]),3),"nmol/L")

plot(fit2,xlab = "log10 (Dose)", ylab = "Cell growth",ylim = c(-0.1, 1.1),
     legend_location = "bottomleft",legend = c(paste0('AF')))
legend('topright',legend=text2, bty ="n", pch=NA) 
dev.off()



#########################################
           
library('ComplexHeatmap') #pheatmap简单好用，但是为了解决更多乱七八糟的需求，你需要学习ComplexHeatmap
library('circlize')

ht1 <- Heatmap(blue_1,show_row_names = F,show_column_names = F)
ht2 <- Heatmap(blue_2,show_row_names = F,show_column_names = F)
draw(ht1 + ht2) #可以对ht1和ht2进行一个组合，默认是竖向的一个整合，需要横轴的数量，元素一致

ht_list = ht1 %v% ht2 #我更加想要的是一个横向的整合，比竖向多一个步骤就好了。
draw(ht_list)

ht2 <- ComplexHeatmap::pheatmap(blue_2,scale = 'column',
                                show_rownames = F,show_colnames = F,
                                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                                legend = T,breaks = seq(-1,1,length.out = 100)) #ComplexHeatmap也可以用pheatmap的参数进行一个绘制，这很贴心，虽然试用下来感觉还是有一些不一样，不过依然很走心了。

##########################ComplexHeatmap######################
library(ComplexHeatmap)
library(circlize)

bar1 <- rowAnnotation(
   sum1 = anno_barplot(
      allms$neoantigen,bar_width = 0.9,
      gp = gpar(col = "white", fill = "orange"),
      border = F,height = unit(2, "cm")), show_annotation_name = F) #给热图添加的注释条，这里我使用了新抗原分数做了一个注释

col_fun1 = colorRamp2(c(0, 1), c('white',met.brewer("Homer1")[7])) #ComplexHeatmap使用类似的方法来调控色条，感觉比pheatmap要方便一些？

Heatmap(mhcarray,col = col_fun1,
        width = unit(8, "cm"),
        height = unit(8, "cm"), #指定width和height是必要的，不然保存的图片比例会很奇怪
        cluster_columns = FALSE,show_row_dend = FALSE,
        right_annotation = bar1,show_heatmap_legend = F,
        row_names_side = "left")

#在热图上方扩展使用HeatmapAnnotation()函数
#也可以直接+bar1来添加注释条
#除了barplot，还有anno_lines，anno_boxplot，anno_histogram，anno_density啥的
#rect_gp = gpar(col= "grey") 可以设置cell轮廓的细节

####################现在主要在用ComplexHeatmap，虽然能做的要更多，但是很多基础的功能要很多步骤才能实现，有点头大############3
library('ComplexHeatmap') 
library('circlize')
mat_scaled = t(scale(t(pheatmap))) #自己手动在那里归一化

df <- data.frame(group = c(rep("P", 96), rep("T", 96)))
ha <- HeatmapAnnotation(df = df, col = list(group = c("P" =  "red", "T" = "blue"))) #分组的注释条，没找到简单实现的方法，相当于花了一个热图作为注释

path<-paste0('./pheatmap/',ncorfs,'.pdf')
hei<-c(nrow(mat_scaled)*(7/6)) #因为批量出图，所以对高度做了一些小修改
col_fun1 = colorRamp2(c(-1,0,1), c(met.brewer("Homer1")[7],'white',met.brewer("Signac")[4])) #对颜色的控制是我一直觉得比pheatmap简单的地方

pdf(path,width = 8,height = 6)
print(Heatmap(mat_scaled,
        na_col = 'grey', #对NA值的处理
        col = col_fun1,
        width = unit(14, "cm"),height = unit(hei, "cm"),
        cluster_columns = FALSE,cluster_rows =F,
        show_row_dend = FALSE,show_column_names =F,
        show_heatmap_legend = F,row_names_side = "left",
        top_annotation = ha))
dev.off()

####################批量出差异分析的图##############
library('ggplot2')
library('ggpubr')
library('MetBrewer')

p=list() #创建一个list名字叫p
for (i in 1:(ncol(test)-1)){
  p[[i]]=ggplot(data=test,aes_string(x="group",y=colnames(test)[i]))+   #很好理解，就是往p这个list里面加ggplot对象就行了
    geom_boxplot(aes(color=group)) + theme_classic()+
    geom_jitter(aes(color=group)) +
    scale_color_manual(values=met.brewer(name="VanGogh1",n=7,type="discrete")[c(6,3)])+
    geom_signif(comparisons = list(c("tumor", "paratumor")),map_signif_level=TRUE,vjust=1.5)  #差异分析的线，对于批量的小图我想用这种方式可能会好一些
}

library(patchwork)
pdf("DEG_plot-down.pdf",width = 8,height = 6)
wrap_plots(p,nrow=2, guides="collect") 
dev.off()


##############基于染色体位置的韦恩图，主要的代码在Merip-seq分析的那一部分，这里仅展示美化。这个函数的能美化的地方还挺多的########
makeVennDiagram(ol, NameOfPeaks =c('NONCODE','GENCODE','Smprot','CNC'),
                fill=c(met.brewer('Hiroshige',n=4)),
                col=c("#D55E00", "#0072B2","#009E73",'#2166AC'), 
                cat.col=c(met.brewer('Hiroshige',n=4)),
                cat.cex =1.2 ) 

#############对于美学代码的固定，是我一直想做的一件事情########
pdf("orf_chr_dis_new.pdf",width = 8,height = 6)
ggplot(chr_inf, aes(x=reorder(Var1,Freq), weight = Freq, fill = group))+
  geom_bar( position = "stack")+ 
  theme_minimal()+theme_classic()+
  labs(y="Frequency",x='chromosome')+
  theme_classic(base_size = 18)+   #例如在这一部分代码中theme有关的这一部分其实是我一直想固定下来的代码
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
  theme(axis.text.x = element_text(size=12,angle=45,hjust=0.9))+
  theme(legend.text = element_text(size=8))+
  scale_fill_manual(values=met.brewer('Hiroshige',n=4))
dev.off()

###############累积曲线图，有时候箱线表现出来不怎么好看的话，用累积曲线或许能有好的效果（？）############

ggplot(immu, aes(neoantigen, colour = group)) + 
  stat_ecdf(geom="smooth", se=F, size=1.2)+   #主要算法就是在这一行
  theme_bw(base_size = 18) +
  theme(legend.position='right',          #c(.25,.75),
        panel.grid = element_blank()) +
  labs(x="Immunogenicity Score", y="ECDF") + 
  scale_color_manual(values=rev(met.brewer("Hokusai1", 7)[c(6,5,3)]))

###############生存分析#######################
library(survminer)
library(survival)
library(ggpubr)
library(gridExtra)

dat<-line1_sv_t

mat<-dat[!is.na(dat$Total.follow.up.period..m.)&
           !is.na(dat$Died.of.recurrence),]
matt<-mat
med.exp<-median(matt$exp)
more.med.exp.index<-which(matt$exp>=med.exp)
less.med.exp.index<-which(matt$exp< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(Total.follow.up.period..m.,Died.of.recurrence) ~ status, data = matt)
s.diff<-survdiff(Surv(Total.follow.up.period..m.,Died.of.recurrence) ~ status, data = matt)

sdata.plot3<-ggsurvplot(s.fit, data=matt,
                        palette="Pastel1",
                        pval = TRUE,pval.method = TRUE, conf.int = TRUE,
                        xlab = 'Time (Month)',ggtheme = theme_survminer(),
                        surv.median.line = 'hv',
                        title=paste0("LINE-1 survival"))

pdf("Line-1-survival.pdf",width = 8,height = 6)
sdata.plot3$plot+
  theme_classic(base_size = 18)+ 
  theme(legend.position = "none")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())
dev.off() 


########################环形热图###########################################
#学起来很痛苦，但是用起来很有意思，有时间再注释
load('~/biodata/circpheat.RData')
allindex<-allindex[allindex$adj_genetype!='Noncode',]
allindex$comm<-as.numeric(allindex$comm)
allindex$chrom<-substr(allindex$chrom, 4, 10)
allindex$deg_group<-as.numeric(allindex$deg_group)

mat<-allindex[,c(7,8,11,12)]
mat<-as.matrix(mat)

split <- factor(allindex$chrom, levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')) #设置分隔

library('ComplexHeatmap')
library('circlize')

col_fun1 = colorRamp2(c(0,1), c('white',met.brewer("Signac")[8])) 
col_fun2 = colorRamp2(c(0,1), c('white',met.brewer("Signac")[13]))

pdf('~/share/Ran/CIRCPHEAT.pdf',width = 8,height = 8)

circos.par(start.degree = 90, gap.degree = 10,gap.after = c(3)) #设置画布

circos.heatmap(mat[,1],col = col_fun1 ,rownames.side = "none",split =split, 
               dend.side='none',na.col = 'white',cluster = F) #,show.sector.labels = TRUE
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + convert_y(2, "mm"), 
              paste0(CELL_META$sector.index),
              facing = "bending.inside", cex = 0.8,front='bond',
              adj = c(0.5,0), niceFacing = TRUE)}, bg.border = NA) #相当于是对每一段分隔的注释，可以调整字体，大小什么的在这里

circos.heatmap(mat[,2], col = col_fun2 ,rownames.side = "none",split =split, #里圈第二个环
               dend.side='none',na.col = 'white',cluster = F,cell_width=8)

col_fun3 = colorRamp2(c(0,1), c('white',met.brewer("Cassatt2")[4])) 
circos.heatmap(mat[,3], col = col_fun3 ,rownames.side = "none",
               split =split,
               dend.side='none',na.col = 'white',cluster = F,cell_width=8)

col_fun4 = colorRamp2(c(-1,0,1), c(met.brewer("Homer1")[7],'white',met.brewer("Signac")[4]))
circos.heatmap(mat[,4], col = col_fun4 ,rownames.side = "none",
               split =split,
               dend.side='none',na.col = 'white',cluster = F,cell_width=8)

dev.off()
circos.clear() #停止绘图，画完不停的话可能会报错什么的
                #感觉画四个环又有点臃肿了，有时间慢慢调整吧
                
########################################韦恩################################
library(ggvenn)

venn<-list(shNC=c(shNC$adj_geneid),shMETTL3=c(shMETTL3$adj_geneid))

p<-ggvenn(venn,fill_color = c("#0073C2FF", "#EFC000FF"),  #,"#CD534CFF"
          stroke_size = 0.5, set_name_size = 4)
                   
library(eulerr) #我想做那种根据集的大小面积不同的韦恩图，eulerr这个包可以做，但是感觉用起来不太方便，也可能是我没有找到合适的方法

fit <- euler(c("ncORF" = 1355, "SRAMP" = 932, "twodatabase" = 152, 
               "ncORF&SRAMP" = 932, "ncORF&twodatabase" = 152, "SRAMP&twodatabase" = 116,
               "ncORF&SRAMP&twodatabase" = 116),
             input='union',shape = "ellipse")
plot(fit)

###########################################################################
#AUC,AUPRC
library('pROC')
Test$bind<-factor(Test$bind)
rf.probs = predict(rfFit,Test,type = "prob")
rf.ROC = roc(response = Test$bind,predictor = rf.probs$BIN,levels = levels(Test$bind))

p<-ggroc(rf.ROC)
pdf("~/tools/pepnn/learning/auc.pdf",width = 8,height = 6)
p+theme_bw(base_size=18)+
  theme(panel.grid.major =element_blank(), 
        panel.background=element_rect(size =1.1,fill='transparent', color='black'),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank())+
  geom_abline(intercept = 1, slope = 1, linetype = "dashed")+
  scale_x_reverse(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  annotate("text", x=0.3, y=0.3, label="AUC = 0.69",size=5)+
  theme(aspect.ratio=1)+
  coord_fixed(ratio = 0.8)
dev.off()

{ #plot_pr,画这个目前还比较麻烦，而且视觉效果也没有很理想
  probs_1= predict(mod,Test) 
  probs_2= predict(mod,Test,type = "prob")
  Test$bind<-factor(Test$bind)
  
  eva<-cbind(probs_2,Test$bind)
  colnames(eva)[3]<-c('obs')
  eva$group<-(modname)
  eva$group<-as.factor(eva$group)
  eva$obs<-as.factor(eva$obs)
  
  library(PRROC)
  library(ROCR)
  scores <- data.frame(eva$BIN)
  scores$labels<-ifelse(eva$obs=='BIN','1','0')
  pr <- pr.curve(scores.class0=scores[scores$labels=="1",]$eva.BIN,
                 scores.class1=scores[scores$labels=="0",]$eva.BIN,
                 curve=T)
  y <- as.data.frame(pr$curve)
  
  pdf("~/tools/pepnn/learning/auprc.pdf",width = 8,height = 6)
  ggplot(y, aes(y$V1, y$V2))+geom_path()+ylim(0,1)+
    theme_bw(base_size=18)+
    theme(panel.grid.major =element_blank(), 
          panel.background=element_rect(size =1.1,fill='transparent', color='black'),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    annotate("text", x=0.3, y=0.5, label="AUPRC = 0.59",size=5)+
    theme(aspect.ratio=1)+
    labs(x="Recall", y = "Precision")+
    coord_cartesian(ylim=c(0,1),xlim=c(0,1))
  dev.off()
  
  plot(pr)
  
}


#还是复杂热图
bar1 <- rowAnnotation(
  Immunogenicity = pheat$prob_2.Positive,
  foo = anno_barplot(
    pheat$prob_2.Positive,bar_width = 0.9,
    gp = gpar(col = "black", fill = "orange"),
    border = T,height = unit(2, "cm")), 
  show_annotation_name = F) #我没想好免疫原性分数要怎么表示，感觉柱状图也好，颜色条也好都没有达到想要的效果，之后再调整吧

col_fun1 = colorRamp2(c(0, 1), c('white',met.brewer("Homer1")[8])) 

split = c(rep(1,13),rep(2,5))
ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4),labels = c("HLA-I", "HLA-II"), 
                                         labels_gp = gpar(col = "white", fontsize = 10))) #给上面添加了HLA-I和HLA-II的注释

ha2 = rowAnnotation(foo = anno_mark(at = c(which(pheat$prob_2.Positive>0.85),which(freq>3)), 
                                   labels = rownames(pheat)[c(which(pheat$prob_2.Positive>0.85), which(freq>3))],side ='left')) 
                   #这里是我主要想记录的，我不想把所有肽的名字都体现出来，只想展示一部分我觉得不错的

pdf("test_pheat.pdf",width = 10,height = 12)
p<-Heatmap(pheat[,-1],col = col_fun1,
        width = unit(10, "cm"),
        height = unit(14, "cm"), #指定width和height是必要的，不然保存的图片比例会很奇怪
        cluster_rows = FALSE,
        cluster_columns = FALSE,show_row_dend = FALSE,
        show_heatmap_legend = F,
        show_row_names = F,right_annotation = bar1,
        column_split = split,column_title = NULL,column_gap = unit(c(1, 1), "mm"),
        top_annotation = ha,
        row_names_side = "left")
draw(ha2+p) #绘图，复杂热图其实可以自定义的地方还挺多的，需要再学习。我看中了热图上面接GO词条什么的，到时候WGCNA做完也可以试试用复杂热图呈现
dev.off()
                   
###################记录一下ggplot的密度，频度，直方，曲线#################################
#这些密度曲线，频度直方图的代码对于提升效率非常重要
                   
ggplot(plot, aes(Position, colour = group)) +
  geom_freqpoly(bins = 50)          #频率曲线，bins可以设置窗口，bins越大越平滑，信息量越少
                   
ggplot(plot, aes(Position, after_stat(density), colour = group)) +
  geom_density(adjust = 1/3, size=1)+      #我更喜欢用密度曲线，adjust可以设置平滑度，值越接近零越粗糙，信息量越大
  xlim(-497,98)+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(size=12))+
  theme(aspect.ratio=0.7)+
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF",'#96345A74'))         
                   
ggplot(hight, aes(Position, fill = Classification)) +
  geom_density(adjust = 1/3, size=1,position = "stack")+  #密度曲线可以做成堆积曲线，比较适合可视化某一组线的组内组成
  labs(title="Detected more than 10% sample")+
  theme_classic(base_size = 18)+
  scale_fill_manual(values = c(met.brewer("Egypt",n=4,type="continuous")))
                   
ggplot(diamonds, aes(price, fill = cut)) +
  geom_histogram(bins = 40)  #同样也可以用直方图做堆积（colour改成fill），不过不太好看，不怎么喜欢用

ggplot(diamonds, aes(price, fill = cut)) +
  geom_histogram(alpha = 0.6, bins = 40) +facet_wrap(~ cut) +
  theme(legend.position = "none")   
                   
ggplot(diamonds, aes(price, fill = cut)) +
  geom_density(alpha = 0.6) +facet_wrap(~ cut) +
  theme(legend.position = "none")    #密度曲线和直方图都可以做成分面图
       
########################################################################################################                 
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
      
