###在jimmy的代码上做的调整，主要在于绘图和图片的保存
###把png图片同一调整为了pdf矢量图，在dobplot上对于labels的长度进行了限制
###dobplot能调整的地方还挺多的，有时间的话以后可以慢慢的美化。
###对labels长度进行限制时需要library(stringr)

run_kegg <- function(gene_up,gene_down,geneList=F,pro='shTRA2A-shNC'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  
  ###   over-representation test
  # 下面把3个基因集分开做超几何分布检验
  # 首先是上调基因集。
  
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk=kk.up
  pdf(file = paste0(pro,'kegg.up.pdf'),width=8,height = 6)
  dotplot(kk)
  dev.off()
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db', keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'kegg.up.csv'))
  
  # 然后是下调基因集。
  
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  kk=kk.down
  pdf(file = paste0(pro,'kegg.down.pdf'),width=8,height = 6)
  dotplot(kk)
  dev.off()
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'kegg.down.csv'))
  
  # 最后是合并基因集。
  
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.5)
  head(kk.diff)[,1:6]
  kk=kk.diff
  pdf(file = paste0(pro,'kegg.diff.pdf'),width=8,height = 6)
  dotplot(kk)
  dev.off()
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'kegg.diff.csv'))
  
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.01,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.01,];up_kegg$group=1
  
  #画图设置, 这个图很丑，大家可以自行修改。
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = paste0(pro,'_kegg_up_down.pdf'),width = 8,height = 6)
  
  if(geneList){
    ###  GSEA 
    ## GSEA算法跟上面的使用差异基因集做超几何分布检验不一样。
    kk_gse <- gseKEGG(geneList     = geneList,
                      organism     = 'hsa',
                      nPerm        = 1000,
                      minGSSize    = 20,
                      pvalueCutoff = 0.9,
                      verbose      = FALSE)
    head(kk_gse)[,1:6]
    gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
    gseaplot(kk_gse, 'hsa04110',title = 'Cell cycle') 
    kk=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    tmp=kk@result
    write.csv(kk@result,paste0(pro,'_kegg.gsea.csv'))
    
    
    # 这里找不到显著下调的通路，可以选择调整阈值，或者其它。
    down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
    up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
    
    g_kegg=kegg_plot(up_kegg,down_kegg)
    print(g_kegg)
    ggsave(g_kegg,filename = paste0(pro,'_kegg_gsea.png'))
    # 
  }
  
}

run_go <- function(gene_up,gene_down,pro='shTRA2A-shNC'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  g_list=list(gene_up=gene_up,
              gene_down=gene_down,
              gene_diff=gene_diff)
  # 因为go数据库非常多基因集，所以运行速度会很慢。
  if(T){
    go_enrich_results <- lapply( g_list , function(gene) {
      lapply( c('BP','MF','CC') , function(ont) {
        cat(paste('Now process ',ont ))
        ego <- enrichGO(gene          = gene,
                        #universe      = gene_all,
                        OrgDb         = org.Hs.eg.db,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file =paste0(pro, '_go_enrich_results.Rdata'))
    
  }
  # 只有第一次需要运行，就保存结果哈，下次需要探索结果，就载入即可。
  load(file=paste0(pro, '_go_enrich_results.Rdata'))
  
  n1= c('gene_up','gene_down','gene_diff')
  n2= c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn=paste0(pro, '_dotplot_',n1[i],'_',n2[j],'.pdf')
      cat(paste0(fn,'\n'))
      pdf(fn,width = 8,height = 6)
      print( dotplot(go_enrich_results[[i]][[j]],color = "pvalue")+
               scale_y_discrete(labels = function(x) str_wrap(x, width = 45)))
      dev.off()
    }
  }
  
  
}

kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}
