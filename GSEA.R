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
  
p<-gseaplot2(edo2, geneSetID = geneset_id, pvalue_table=F,subplots = 1:3,
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

rownames(count)<-substr(rownames(count),1,18)
exp<-count[match(gene$X,rownames(count)),]
exp<-exp[,c(4,5,6,1,2,3)]
rownames(exp)<-gene$SYMBOL
exp<-exp[1:20,]

mat_scaled = t(scale(t(exp)))
col_fun1 = colorRamp2(c(-1,0,1), c(met.brewer("Hiroshige")[9],'white',met.brewer("Signac")[4])) 

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

 path<-paste0(edo2$Description[which(edo2@result$ID==geneset_id)],'_gseapheat.pdf')

 test<-aplot::plot_list(gglist = p, ncol = 1)
 
 pdf(path,width = 12,height = 5)
 print(plot_grid(test,grob,  nrow=1)) #画到一起去
 dev.off() 

}

gsea_pheat('GO:0019724')  
gsea_pheat('GO:0002377')  
gsea_pheat('GO:0003823')
gsea_pheat('GO:0006959')

