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
                       
                       

                       
