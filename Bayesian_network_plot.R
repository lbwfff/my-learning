#对CBNplot的代码进行修改，使其能够对更广泛的对象进行绘图

##########################################################
gobp<-read.csv('human_Gobp_database.csv')
gobp<-gobp[grep('GO:0006119',gobp$GO),] #氧化磷酸化

#基因名的替换是个麻烦的事，这里把蛋白质序列对于的uniprot名改成了基因名
unipn<-as.data.frame(rownames(exp))
  
f <-function(x) unlist(strsplit(x['rownames(exp)'],'[|]'))[2]
unipn$unipn<-apply(unipn,1,f)

ens<-clusterProfiler::bitr(unipn$unipn, fromType="UNIPROT", toType=c("ENSEMBL",'SYMBOL'), OrgDb=org.Hs.eg.db)
ens<-ens[match(unipn$unipn,ens$UNIPROT),]
unipn$ensembl<-ens$ENSEMBL
unipn$symbol<-ens$SYMBOL

oxida<-unipn[unipn$ensembl %in% gobp$GENE,]
oxida<-oxida[!duplicated(oxida$ensembl),]

oxiexp<-exp[match(oxida$`rownames(exp)`,rownames(exp)),]
row.names(oxiexp)<-oxida$symbol

corpep<-orf_ins[rownames(orf_ins) %in% list,]
corpep<-log2(corpep+1)

colnames(corpep)<-colnames(oxiexp)
testpcs<-data.frame(t(rbind(oxiexp,corpep)))

##########################################################
#处理完矩阵后是绘图，原包使用了bnlearn来计算网络，然后使用igraph，ggraph绘制网络
  
  geneNames <- colnames(testpcs)
  pcsRaw <- testpcs 
  strength <- withr::with_seed(seed = seed, bnlearn::boot.strength(testpcs, algorithm=algo, algorithm.args=algorithm.args, R=R, cluster=cl))
  
  pepstrength<-strength[strength$from %in% list | strength$to %in% list,]
  ## Average by specified threshold
  av <- bnlearn::averaged.network(strength, threshold=0.8) #这里输入强度阈值

  av <- bnlearn::cextend(av, strict=FALSE)
  
  g <- bnlearn::as.igraph(av)
  e <- igraph::as_edgelist(g, names = TRUE)
    
  eName <-paste0(e[,1], "_", e[,2])
  colnames(e) <- c("from","to")
  eDf <- merge(e, strength)
  rownames(eDf) <- paste0(eDf$from, "_", eDf$to)
  eDf <- eDf[eName, ]
  g <- igraph::set.edge.attribute(g, "color", index=igraph::E(g), eDf$strength)

  g <- igraph::set.edge.attribute(g, "label", index=igraph::E(g), round(eDf$direction,2))

  ## Hub genes
  hScore <- igraph::hub.score(g, scale = TRUE, weights = igraph::E(g)$color)$vector
    
  defHub <- hScore[order(hScore, decreasing=TRUE)][seq_len(hub)]
  nodeShape <- names(igraph::V(g)) %in% names(defHub)
  nodeShape <- ifelse(nodeShape, 19, 21)
  igraph::V(g)$shape <- nodeShape

  igraph::E(g)$width <- igraph::E(g)$color
  edgeWName <- "strength"
    
  sizeLab <- "expression"
  meanExpCol <- vapply(names(igraph::V(g)), function(x) ifelse(x %in% geneNames, mean(pcsRaw[, x]), NA), FUN.VALUE=1)
  meanExpSize <- meanExpCol
  meanExpSize[is.na(meanExpSize)] <- 1
  igraph::V(g)$size <- meanExpSize
    
  igraph::V(g)$color <- meanExpCol
    
  ## Plot
  delG <- igraph::delete.vertices(g, igraph::degree(g)==0)
      
  p <- ggraph::ggraph(delG, layout=layout) +
       ggraph::geom_edge_diagonal(edge_alpha=1,
                               position="identity",aes_(edge_colour=~I(color), width=~I(width), label=~I(label)),
                               label_size=3*(labelSize/4),
                               label_colour=NA,angle_calc = "along",
                               label_dodge=unit(3,'mm'),
                               arrow=arrow(length=unit(4, 'mm')),
                               end_cap=ggraph::circle(5, 'mm'))+
      ggraph::geom_node_point(aes_(color=~color, size=~size, shape=~shape), show.legend=TRUE)+
                              scale_color_continuous(low="blue", high="red", name="expression") +
                              scale_size(range=c(scaleSizeLow, scaleSizeHigh) * cexCategory, name=sizeLab)+
      ggraph::scale_edge_width(range=c(1, 3), guide="none")+
      ggraph::scale_edge_color_continuous(low="dodgerblue", high="tomato", name="strength")+
                             guides(edge_color = ggraph::guide_edge_colorbar(title.vjust = 3))+
                             scale_shape_identity()+
      ggraph::theme_graph() + ggtitle('') #这里可以输入title
      
      p <- p + ggraph::geom_node_text(aes_(label=~stringr::str_wrap(name, width = 25)),
                                check_overlap=TRUE, repel=TRUE, size = labelSize,
                                color = textColor,bg.color = bgColor, segment.color="black",bg.r = .15)

      
pdf("justfortest_2.pdf",width = 8,height = 6)      
p
dev.off() #这里得到的图像没有title和legend不知道为什么
