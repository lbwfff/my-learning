#把kegg，GO，DO的结果画一起了，用了一个自己做的柱状图代码，总感觉还是差点什么，主要的框架就是这样了，要改的话还是想改绘图的代码

for (i in 1:4){
  p<-list()
  enrich<-unique(list[[i]])
  kk.negative <- enrichKEGG(gene  = enrich,
                            organism = "hsa",
                            #universe = gene_all,
                            pvalueCutoff = 0.5,
                            qvalueCutoff =0.5)
  head(kk.negative)[,1:6]
  kk=kk.negative@result
  kk<-kk[order(kk$pvalue,decreasing = F),]
  
  library('stringr')
  title<-paste0('KEGG')
  color<-c(met.brewer('Hokusai1')[c(3,5,6,7)],met.brewer('Troy')[5])
  
  p[[1]]<-
    ggplot(data = kk[1:10,], aes(x = reorder(Description,-pvalue), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
    geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
    ggtitle(title) +coord_flip() +
    theme_classic(base_size = 13)+ 
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=12,face='bold'))+
    labs(y="-log10(P value)",x='Term')+
    scale_fill_manual(values = c(color[1]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))
  
    BP <- enrichGO(gene          = enrich,
                      OrgDb         = org.Hs.eg.db,
                      ont           = 'BP' ,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.2,
                      qvalueCutoff  = 0.5,
                      readable      = TRUE)
    BP <- simplify(BP,cutoff = 0.7,by = "p.adjust",
                      select_fun = min,measure = "Wang",semData = NULL)
    
    CC <- enrichGO(gene          = enrich,
                   OrgDb         = org.Hs.eg.db,
                   ont           = 'CC' ,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.2,
                   qvalueCutoff  = 0.5,
                   readable      = TRUE)
    CC <- simplify(CC,cutoff = 0.7,by = "p.adjust",
                   select_fun = min,measure = "Wang",semData = NULL)
    
    MF <- enrichGO(gene          = enrich,
                   OrgDb         = org.Hs.eg.db,
                   ont           = 'MF' ,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.2,
                   qvalueCutoff  = 0.5,
                   readable      = TRUE)
    MF <- simplify(MF,cutoff = 0.7,by = "p.adjust",
                   select_fun = min,measure = "Wang",semData = NULL)
    
    DO<-enrichDO(gene          = enrich,
                 ont           = "DO",
                 pvalueCutoff  = 0.2,
                 pAdjustMethod = "BH",
                 #universe      = names(geneList),
                 minGSSize     = 5,
                 maxGSSize     = 500,
                 qvalueCutoff  = 0.2,
                 readable      = FALSE)
      
    BP=BP@result
    BP<-BP[order(BP$pvalue,decreasing = F),]
    CC=CC@result
    CC<-CC[order(CC$pvalue,decreasing = F),]
    MF=MF@result
    MF<-MF[order(MF$pvalue,decreasing = F),]
    DO=DO@result
    DO<-DO[order(DO$pvalue,decreasing = F),]
    
    title<-paste0('GO_BP')
    p[[2]]<-
      ggplot(data = BP[1:10,], aes(x = reorder(Description,-pvalue), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
      geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
      ggtitle(title) +coord_flip() +
      theme_classic(base_size = 13)+ 
      theme(panel.grid.major =element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            legend.title = element_blank(),
            text = element_text(size=12,face='bold'))+
      labs(y="-log10(P value)",x='Term')+
      scale_fill_manual(values = c(color[2]))+
      scale_colour_manual(values = c("black"))+
      scale_y_continuous(expand = c(0,0))+
      scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
      theme(plot.title = element_text(size=12))
    
    title<-paste0('GO_CC')
    p[[3]]<-
      ggplot(data = CC[1:10,], aes(x = reorder(Description,-pvalue), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
      geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
      ggtitle(title) +coord_flip() +
      theme_classic(base_size = 13)+ 
      theme(panel.grid.major =element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            legend.title = element_blank(),
            text = element_text(size=12,face='bold'))+
      labs(y="-log10(P value)",x='Term')+
      scale_fill_manual(values = c(color[3]))+
      scale_colour_manual(values = c("black"))+
      scale_y_continuous(expand = c(0,0))+
      scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
      theme(plot.title = element_text(size=12))
    
    title<-paste0('GO_MF')
    p[[4]]<-
      ggplot(data = MF[1:10,], aes(x = reorder(Description,-pvalue), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
      geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
      ggtitle(title) +coord_flip() +
      theme_classic(base_size = 13)+ 
      theme(panel.grid.major =element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            legend.title = element_blank(),
            text = element_text(size=12,face='bold'))+
      labs(y="-log10(P value)",x='Term')+
      scale_fill_manual(values = c(color[4]))+
      scale_colour_manual(values = c("black"))+
      scale_y_continuous(expand = c(0,0))+
      scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
      theme(plot.title = element_text(size=12))
    
    title<-paste0('DO') #names(list)[i],
    p[[5]]<-
      ggplot(data = DO[1:10,], aes(x = reorder(Description,-pvalue), y = -log10(pvalue))) + #####这里有一个reorder函数，可以对元素进行排序，这里使用y轴值进行了重新排列
      geom_col(aes(fill='#d8443c',colour='#d8443c'),show.legend = FALSE,width=0.8) +
      ggtitle(title) +coord_flip() +
      theme_classic(base_size = 13)+ 
      theme(panel.grid.major =element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            legend.title = element_blank(),
            text = element_text(size=12,face='bold'))+
      labs(y="-log10(P value)",x='Term')+
      scale_fill_manual(values = c(color[5]))+
      scale_colour_manual(values = c("black"))+
      scale_y_continuous(expand = c(0,0))+
      scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
      theme(plot.title = element_text(size=12))
    
    path<-paste0(names(list)[i],'_enrich.pdf')
    pdf(path,width = 15,height = 6)
    print(wrap_plots(p,nrow=2)) 
    dev.off()
    
    path<-paste0(names(list)[i],'_enrich.xlsx')
    library(openxlsx)  #把富集的结果保存为xlsx文件，这个包还挺方便的
    sheets = list("KEGG" = kk,"GO_BP" = BP,'GO_CC'=CC,'GO_MF'=MF,'DO'=DO)
    write.xlsx(sheets,path)
}

                       
#####################################################################################
pway = ReactomePA::enrichPathway(gene = cand.entrez) # Reactome的富集，结果的格式和GO，KEGG什么的没有区别
pway = setReadable(pway, OrgDb=org.Hs.eg.db)

pwayGSE <- ReactomePA::gsePathway(geneList) #也可以做GSEA，也是需要对基因进行排序                       
                       
library('CBNplot')
bngeneplot(results = pway, exp = vsted, pathNum = 17)
bngeneplot(results = pway, exp = vsted, pathNum = 17, labelSize=7, shadowText=TRUE)
bngeneplot(results = pway, exp = vsted, expSample = incSample, pathNum = 17)       #这玩意叫，香草图？可以看到蛋白与蛋白的相互作用网络    
                       
                       
#################################################################################################################
#GSVA分析
                       
library("GSEABase")
library("GSVA")

geneset <- getGmt('./GSVA/m5.go.bp.v2022.1.Mm.symbols.gmt')  #官网下载的数据集，一个压缩包就把一堆数据集下载了

#使用DESEQ2做了标准化，这里我比较纠结是使用归一化之后的数据还是标准化后的Counts值，后来发现都没有区别
                       
normalized_counts <- counts(dds,normalized=T) 
normalized_counts <- exprSet_new
normalized_counts<-as.data.frame(normalized_counts)
normalized_counts$id<-substr(rownames(normalized_counts),1,18)
                       
k <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
all_gene<-AnnotationDbi::select(org.Mm.eg.db,keys = k,columns = c("SYMBOL",'ENTREZID','GENENAME','GENETYPE'),keytype="ENSEMBL")
all_gene<-all_gene[match(normalized_counts$id,all_gene$ENSEMBL),]
normalized_counts$sym<-all_gene$SYMBOL                             #这一段是在做ID的转换

normalized_counts<-normalized_counts[!is.na(normalized_counts$sym),]
normalized_counts<-normalized_counts[order(normalized_counts$PBS_1.sorted.bam,decreasing = T),]
normalized_counts<-normalized_counts[!duplicated(normalized_counts$sym),]
rownames(normalized_counts)<-normalized_counts$sym
normalized_counts<-normalized_counts[,-c(7:8)]

es <- gsva(as.matrix(normalized_counts), geneset,
           min.sz=10, max.sz=500, verbose=TRUE) #核心代码
es<-as.data.frame(es)
es$de<-rownames(es)
immugsva<-es[grep('IMMU',es$de),]
immugsva<-immugsva[,-7]

immugsva2<-immugsva[-grep('NEGATIVE',rownames(immugsva)),] 
rownames(immugsva2)<-substr(rownames(immugsva2),6,999)
write.csv(immugsva2,file = 'immugsva2.csv')

immugsva2<-read.csv('immugsva2.csv') #我手动去除了一些词条，这里得到的东西实在是过于冗余了
rownames(immugsva2)<-immugsva2$X
immugsva2<-immugsva2[,-1]
immugsva2<-immugsva2[,c(4,5,6,1,2,3)]
#然后就是做热图

library(pheatmap)
annotate_b<-as.data.frame(array(NA,c(6,2)))

colnames(annotate_b)<-c('group','whatever')
annotate_b$group<-c('DRUG','DRUG','DRUG','PBS','PBS','PBS')
annotate_b<-as.data.frame(annotate_b$group)
rownames(annotate_b)<-c(colnames(immugsva))
colnames(annotate_b)<-c('group')

pdf("GSVA_test.pdf",width =14,height = 12)
pheatmap(immugsva2,cluster_rows = T,cluster_cols =F,
         annotation_col =annotate_b, annotation_legend=TRUE, 
         scale = "none",
         cellwidth=36,cellheight=12,
         # color = rev(colors),
         show_rownames = T,show_colnames = F,
         breaks = seq(-0.5,0.5,length.out = 100))
dev.off()                        
                       
                       
##############################################################
#蛋白PPI互作的网络
                       
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(STRINGdb)
library(igraph)
library(ggraph)

string_db <- STRINGdb$new( version="11.5", species=10090,
                           score_threshold=400, input_directory="")
gene <- c('Gper1',Gper1$id)
gene <- gene %>% bitr(fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Mm.eg.db", 
                      drop = T)
data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", 
                                      removeUnmappedRows = TRUE)
# string_db$plot_network( data_mapped$STRING_id ) #这张画出来应该和网页版是一样的，但不知道为什么网络有问题
 
#
data_links <- data_mapped$STRING_id[1:100] %>% string_db$get_interactions()

links <- data_links %>%
  mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
  mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
  dplyr::select(from, to , last_col()) %>% 
  dplyr::rename(weight = combined_score)

links_2 <- links %>% mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
  mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
  filter(!(from_c == 1 & to_c == 1)) %>%
  dplyr::select(1,2,3)
# 新的节点数据
nodes_2 <- links_2 %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
# 创建网络图
net_2 <- igraph::graph_from_data_frame(d=links_2,vertices=nodes_2,directed = F)
# 添加必要的参数
igraph::V(net_2)$deg <- igraph::degree(net_2)
igraph::V(net_2)$size <- igraph::degree(net_2)/5
igraph::E(net_2)$width <- igraph::E(net_2)$weight/10

pdf('test_ppi.pdf',width = 8,height = 6)
ggraph(net_2,layout = "centrality", cent = deg)+
  geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = T)+  #比较边缘的点会在这里被过滤掉
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10) )+
  guides(size=F)+
  theme_graph()
dev.off() 
#就，目前也还能行，ggraph的美化又是一个大坑，应付甲方的话图画成这样也随便了    





#############################################################
##富集分析的弦图，这种画法还挺好看的，比起clusterprofiler那个网络，但是包写得很呆，需要自己调一些地方
library('org.Hs.eg.db')
library(patchwork)
library(clusterProfiler)
library(stringr)
library(DOSE)
library(openxlsx)
library(ggplot2)
R.utils::setOption("clusterProfiler.download.method",'auto')

genelist<-bitr(unique(tr), fromType = "SYMBOL",
               toType = c( "ENTREZID"),
               OrgDb = org.Hs.eg.db)[,2] 
kk.negative <- enrichKEGG(gene  = genelist,
                          organism = 'hsa',
                          pvalueCutoff = 0.05,
                          qvalueCutoff =1)
kk.negative<-setReadable(kk.negative,OrgDb = 'org.Hs.eg.db', keyType="ENTREZID")
kk=kk.negative@result
kk<-kk[order(kk$pvalue,decreasing = F),]

write.csv(kk.negative@result,file = '02_DEG/KEGG_enrich.csv')

{
  library('GOplot')
  
  david<- kk[,c(1,3,4,10,8)]
  colnames(david)
  
  david$geneID<-gsub("/",",",david$geneID)
  colnames(david)<-c('category', 'ID', 'term','genes','adj_pval')
  david[1:5,1:5]
  
  genelist<- data.frame(ID=deseq$X,logFC=deseq$log2FoldChange)  #然后对于每一个基因都是可以展示logFC的
  colnames(genelist)
  
  circ <- circle_dat(david, genelist)#创建绘图对象
  
  ####
  process<-kk$Description[1:5]
  
  chord <- chord_dat(data=circ,genelist,process=process)#为什么又要把logfc加进去，这个包实在粗糙
  head(chord)
  
  p<-GOChord(chord, space = 0.02,
             gene.space = 0.28, 
             gene.size = 4,process.label=10,
             ribbon.col=brewer.pal(5, "Set3")) 
  #感觉像是一个魔改的ggplot对象，没法直接用ggplot的美学代码调整，但是可以对象进行操作
  p[["guides"]][["size"]][["title"]]<-c('KEGG')
  p[["guides"]][["size"]][["ncol"]]<-2
  pdf("./02_DEG/kegg_chord.pdf",width =10,height = 12)
  print(p)
  dev.off()  #好多了
  
}

#也展示一个GOBP的，CC和MF如法炮制就好了
genelist<-bitr(unique(tr), fromType = "SYMBOL",
               toType = c( "ENTREZID"),
               OrgDb = org.Hs.eg.db)[,2] 
BP <- enrichGO(gene          = genelist,
               OrgDb         = org.Hs.eg.db,
               ont           = 'BP' ,
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 1,
               readable      = TRUE)

write.csv(BP@result,file = '02_DEG/GOBP_enrich.csv')

{
  david<- BP@result[,c(1,2,6,8)]
  colnames(david)
  
  david$geneID<-gsub("/",",",david$geneID)
  colnames(david)<-c('ID', 'term','adj_pval','genes')
  david$category<-c('GO')
  david[1:5,1:5]
  
  genelist<- data.frame(ID=deseq$X,logFC=deseq$log2FoldChange) 
  colnames(genelist)
  
  circ <- circle_dat(david, genelist)#创建绘图对象
  
  ####
  process<-BP@result$Description[1:5]
  
  chord <- chord_dat(data=circ, genelist, process=process)#构建数据
  head(chord)
  
  p<-GOChord(chord, space = 0.02,
             gene.space = 0.28, 
             gene.size = 4,process.label=10,
             ribbon.col=brewer.pal(5, "Set3")) 
  p[["guides"]][["size"]][["title"]]<-c('GO BP')
  p[["guides"]][["size"]][["ncol"]]<-2
  
  pdf("./02_DEG/GOBP_chord.pdf",width =10,height = 12)
  print(p)
  dev.off() 
}


#################################################################
#GSEA的美化，主要是对于GseaVis包的使用
library('GseaVis')
p<-list()

p[[1]]<-
gseaNb(object = egmt,curveCol =met.brewer("Hokusai1", 5),
       geneSetID = egmt@result$Description[1:5])

gmtfile ='./05_enrich/c5.go.v2023.1.Hs.symbols.gmt'
geneset <- clusterProfiler::read.gmt( gmtfile )  
length(unique(geneset$term))
egmt <- GSEA(test, TERM2GENE=geneset, 
             minGSSize = 1,eps=0,
             pvalueCutoff = 0.1,
             verbose=FALSE)
head(egmt)
gsea_results_df <- egmt@result 

write.csv(gsea_results_df,file = '05_enrich/gsea_go_results_df.csv')

p[[2]]<-
gseaNb(object = egmt,curveCol =met.brewer("Hokusai1", 5),
       geneSetID = egmt@result$Description[1:5])

pdf("05_enrich/two_GSEA.pdf",width = 10,height = 12)
wrap_plots(p,nrow=2) #这张图目前的画法，感觉还是比较臃肿，其实东西多的话可以把底部的rank去掉也行
dev.off() 

##################################################################################
#再补充一个超几何分布的棒棒糖图，以前常用的柱状图现在怎么看怎么丑
plot<-data.frame(category=c(rep('GO',10),rep('KEGG',10)),
                 Description=c(BP$Description[1:10],kk$Description[1:10]),
                 pvalue=c(BP$pvalue[1:10],kk$pvalue[1:10])) #合并KEGG和GO的结果

library('ggsci')
library('cowplot')
library('stringr')
pdf('05_enrich/100gene_enrich.pdf',width = 10,height = 8)
  ggplot(data = plot, aes(x = reorder(Description,pvalue), y = -log(pvalue))) +
  geom_segment(aes(x = Description, y = 0, xend = Description, 
                   yend = -log(pvalue), color = category)) +
  geom_point(size = 4, aes(color = category))+
  scale_color_npg() +
  # scale_y_continuous(expand = c(0,0)) +
  theme_half_open() +
  ylab('-log (P value)')+
  xlab(NULL)+
  theme(axis.text.x = element_text(angle = 45,size = 9,
                                   hjust = 1,vjust = 1))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+  #和之前一样，把太长的词条换行
  theme(plot.margin=unit(rep(2,4),'cm'))
dev.off()



                       
