#在R上面做方便一些，但是不知道要怎么用自定义的基因集什么，而且画网络不如enrichmap好看

##################################################################
#和GSEA没关系，记录一下怎么提取全org.Mm.eg.db的注释什么的
library(org.Mm.eg.db)
columns(org.Mm.eg.db)
k <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
all_gene<-AnnotationDbi::select(org.Mm.eg.db,keys = k,columns = c("SYMBOL",'ENTREZID','GENENAME','GENETYPE'),keytype="ENSEMBL")

###################################################################
gsea<-deseq[!is.na(deseq$ENTREZID),]
gsea$value<-(abs(log(gsea$pvalue)))
gsea$value<-(gsea$value*sign(gsea$log2FoldChange))
gsea<-gsea[!is.na(gsea$value),]
gsea<-gsea[!duplicated(gsea$ENTREZID),]
rownames(gsea)<-gsea$ENTREZID
gsea<-gsea[order(gsea$log2FoldChange,decreasing = T),]

test<-as.numeric(gsea$log2FoldChange)
names(test) = as.character(gsea$ENTREZID) #记住names这个属性
edo2 <- gseGO(test,OrgDb= org.Mm.eg.db,ont='ALL',pvalueCutoff = 0.5)

# edo2 <- simplify(edo2,cutoff = 0.8,by = "p.adjust",
#                select_fun = min,measure = "Wang",semData = NULL) #ont是ALL的话，好像就没办法简化了

library('ggridges')
pdf('gsea_with_drug1_ridge.pdf',width = 12,height = 8)
ridgeplot(edo2,showCategory = 15,fill='p.adjust',decreasing=T)
dev.off()

edo2<-setReadable(edo2,OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(edo2@result,file = 'GSEA_result.csv')

library('ComplexHeatmap') 
library('circlize')
library('MetBrewer')
library('cowplot')
library('enrichplot')
library('ggplot2')

gsea_pheat<-function(geneset_id) {
  
p<-gseaplot2(edo2, geneSetID = geneset_id, pvalue_table=F,subplots = 1:3,rel_heights = c(2, 0.8, 1.2),
                      title = edo2$Description[which(edo2@result$ID==geneset_id)])

p[[1]] <- p[[1]]+ annotate("text", x = 16000, y = 0.4, label = paste0('P.value = ',format(edo2@result$pvalue[which(edo2@result$ID==geneset_id)],2),
                                                                      '\n NES = ',round(edo2@result$NES[which(edo2@result$ID==geneset_id)],2)),
                           size = 5)
# print(p)

gene<-edo2@result$core_enrichment[edo2@result$ID==geneset_id]
gene<-unlist(strsplit(gene,'/'))
gene<-bitr(gene, fromType = "ENTREZID",
         toType = c( "ENSEMBL"),
         OrgDb = org.Mm.eg.db)[,2] 
gene<-deseq[deseq$X %in% gene,]
gene<-gene[order(gene$pvalue,decreasing = F),]
gene<-gene[!duplicated(gene$SYMBOL),]

rownames(count)<-substr(rownames(count),1,18)
exp<-count[match(gene$X,rownames(count)),]
exp<-exp[,c(4,5,6,1,2,3)]
rownames(exp)<-gene$SYMBOL

if(nrow(exp)>=20) {exp<-exp[1:20,]} else{exp<-exp}

mat_scaled = t(scale(t(exp)))
col_fun1 = colorRamp2(c(-1,0,1), c(met.brewer("Hiroshige")[9],'white',met.brewer("Signac")[4])) #ComplexHeatmap使用类似的方法来调控色条，感觉比pheatmap要方便一些？

df <- data.frame(group = c(rep("PBS", 3), rep("Drug", 3)))
df$group<-factor(df$group)

split = rep(1:2, each = 3)
ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4),labels = c("PBS", "Drug"), 
                        labels_gp = gpar(col = "white", fontsize = 10)))

 grob <- grid.grabExpr(draw(Heatmap(mat_scaled,rect_gp = gpar(col = "white", lwd = 1),
        col = col_fun1,name = "Score",
        column_split = split,column_title = NULL,column_gap = unit(c(0.1, 0.1), "mm"),
        width = unit(10, "cm"),
        height = unit(8, "cm"),
        cluster_columns = FALSE,
        cluster_rows =F,
        show_row_dend = FALSE,
        show_column_names =F,
        show_heatmap_legend = T,
        row_names_side = "right",
        top_annotation = ha))) 

 path<-paste0('./gsea_pheat/',edo2$Description[which(edo2@result$ID==geneset_id)],'_gseapheat.pdf')

 test<-aplot::plot_list(gglist = p, ncol = 1,heights=c(2, 0.8, 1.2))
 
 pdf(path,width = 12,height = 5)
 print(plot_grid(test,grob,  nrow=1)) 
 dev.off() 

}

list<-edo2@result[edo2@result$pvalue<0.01,]
for (i in list$ID){
  gsea_pheat(i)
}

########################################################################################################
#GSEA gene_set的散点图
#主要记录一下绘图的代码吧，要怎么做心中还是有很多疑惑

gsea_dot<-function(geneset_id) {
  gene<-all_go2[grep(geneset_id,all_go2$GO),] #要想得到某一个GO词条的全部基因，就还挺麻烦的，单独拿一个代码页写怎么做的吧
  gene<-deseq[deseq$X %in% gene$GENE,]
  gene<-gene[order(gene$pvalue,decreasing = F),]
  
  gene2<-edo2@result$core_enrichment[edo2@result$ID==geneset_id] #首先core_enrichment这里得到的并不是gene_set的全部基因
  gene2<-unlist(strsplit(gene2,'/'))                             #但有一些基因明显也是符合趋势的但没放在core_enrichment里面我不太明白
  gene2<-bitr(gene2, fromType = "ENTREZID",
             toType = c( "ENSEMBL"),
             OrgDb = org.Mm.eg.db)[,2]
  gene2<-deseq[deseq$X %in% gene2,]
  gene2<-gene2[1:15,]
  gene<-gene[!duplicated(gene$SYMBOL),]
  
  # gene<-gene[1:15,]
  gene$ext<-c('T')
  gene<-gene[match(rownames(dot),gene$X),]
  dot$ext<-gene$ext
  dot$ext[is.na(dot$ext)]<-c('F')
  
  path<-paste0('./dot/',edo2$Description[which(edo2@result$ID==geneset_id)],'_dot.pdf')
  pdf(path,width = 8,height = 6)
  print(ggplot(data = dot,aes(x=pbsmean,y=drugmean,color=group))+
    geom_point()+scale_colour_manual(values=c("#808080","grey84","#808080"))+
    theme_classic()+
    xlab('PBS mean log2(normalized counts)')+
    ylab('Drug mean log2(normalized counts)')+
    theme_classic(base_size = 18)+  
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16))+
    theme(axis.text.x = element_text(size=12))+
    theme(legend.text = element_text(size=12))+
    theme(aspect.ratio=1)+
    geom_abline(slope=1)+
    geom_abline(slope=1,intercept=1,linetype = "dashed")+
    geom_abline(slope=1,intercept=(-1),linetype = "dashed")+
    geom_point(data = dot[dot$ext=='T',],color='#e1846c')+
    geom_label_repel(data = dot[rownames(dot) %in% gene2$X,],aes(label = SYMBOL),
                     size = 3,segment.color = "black", show.legend = FALSE ,color='#e1846c') #绘图代码
    )
  dev.off() 
  }

gsea_dot('GO:0050778')
gsea_dot('GO:0002253')

#############################
#一个把跨物种Gene name进行转换的方法
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

local1<-rbind(fread('./subcell_gene/transcript/Hacisuleyman_24_dep.csv'),
                fread('./subcell_gene/transcript/Hacisuleyman_24_rest.csv'))
local1<-local1[!duplicated(local1$GeneName),]
  

  transname<-getLDS(attributes = c("mgi_symbol"),
         filters = "mgi_symbol", values = local1$GeneName,
         mart = mouse,
         attributesL = c("hgnc_symbol"), 
         martL = human)$HGNC.symbol

