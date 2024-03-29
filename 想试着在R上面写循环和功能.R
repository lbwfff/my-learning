####我觉得现阶段限制自己使用R语言处理数据速度的因素之一就是不能较好的使用merge函数，稍微会一点都会帮助许多，至少不用手动去match####
test1<-data.frame(A=c('1','2','1','3'),
                     B=c(rep(NA,4)))
test2<-data.frame(A=c('1','2','3'),
                  C=c('A','B','C'))
merge(test1,test2)

ensembl_id <- substr(row.names(norm),1,15) #随便找个地方放一下代码

####学习循环是为了节省劳动力

#案例一：
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


#这种用法在linux上也用到过，for (i in ""){},感觉比较能理解

#案例二：

test_kegg <- function(gene_up,pro='shTRA2A-shNC'){
  gene_up=unique(gene_up)
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk=kk.up
  pdf(file = paste0(pro,'test.kegg.up.pdf'),width=8,height = 6)
  print(dotplot(kk,color = "pvalue")+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)))
  dev.off()
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db', keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'kegg.up.csv'))
}

test_kegg(gene_up,pro = 'ipe')


##试了一下，可以运行，暂时先体会这两种方法，可以减少非常多的无效代码


实战一：
list<-c('Abcam','CST','SYSY')
for (i in list){
  path1<-paste0("~/biodata/merip/exomepeak2/",i,"-1/Mod.bed")
  path2<-paste0("~/biodata/merip/exomepeak2/",i,"-2/Mod.bed")
gr1 <- toGRanges(path1, format="BED", header=FALSE) 
gr2 <- toGRanges(path2, format="BED", header=FALSE) 

peaks <- GRangesList(rep1=gr1,
                     rep2=gr2)
pdf(file = paste0(i,'genomicElementDistribution.pdf'),width=6,height = 4)
genomicElementDistribution(peaks, 
                           TxDb = TxDb,
                           promoterRegion = c(upstream = 2000, downstream = 100),
                           geneDownstream = c(upstream = 0, downstream = 1000),
                           labels = list(geneLevel = c(promoter = "Promoter", geneDownstream = "Downstream",geneBody = "Gene body", distalIntergenic = "Distal Intergenic"), 
                                         ExonIntron = c(exon= "Exon", intron = "Intron", intergenic = "Intergenic"), 
                                         Exons = c(utr5 = "5' UTR",utr3 = "3' UTR", CDS = "CDS", otherExon = "Other exon"), 
                                         group = c(geneLevel ="Gene Level", promoterLevel = "Promoter Level", Exons = "Exon level", ExonIntron ="Exon/Intron/Intergenic")),
                           labelColors = c(promoter = "#D55E00", geneDownstream = "#E69F00", geneBody ="#51C6E6", 
                                           distalIntergenic = "#AAAAAA", exon = "#009DDA", intron = "#666666",
                                           intergenic = "#DDDDDD", utr5 = "#0072B2", utr3 = "#56B4E9", CDS = "#0033BF",
                                           otherExon = "#009E73")
)
dev.off()
}


#######apply循环的使用#######
#######其实多用apply会比for循环更好一些，因为apply更加简洁高效也更容易去做多线程的计算，但是更难理解一些
                     
f<-function(x) sum(x > 0 & !is.na(x))
apply(fc,1,f)####1为行，2为列，此处计算fc矩阵中每一行大于零且不为na的值的数目

#######进阶实例########
######使用poolr包的stouffer函数对p矩阵的每一行进行p值的合并######
######此例可以让你更好的理解apply循环，function的x是是什么？如果想把function应用在apply上，此处的x应该是行（或列）#####
                     
library('poolr')
f<-function(x) stouffer(as.numeric(x[!is.na(as.numeric(unlist(x)))]))$p
allgene$integrateP<-apply(p,1,f)

                     
f <-function(x) unlist(strsplit(x['experiment_title'],';'))[1]
list2$plate <-apply(list2,1,f)
####这样就可以取得任意封号前的字符了，换成'_',或者其他什么也可以操作，解决了一个大麻烦#####

####除了拆开，也可以合起来，代码更加简单####
protein<-paste0(protein,collapse = ';')

                     

#####判断句#####
#####基本的使用方法是，ifelse(对象,yes,no)######
#####此处的应用是，对于diff矩阵，如果v1，v2同时为0或同时为翻译，则输出‘no change’，若不则进行下一步的判断。在下一步的判断中，若v1为0，v2为翻译，则输出‘turn-on’，否则输出‘turn-off’########
#####需要解释的是此处的0原先为na，但是我发现如果用na进行判断则会使逻辑变得更加复杂，所以将na改为了0进行简化######
diff$change = ifelse((diff$V1==0 & diff$V2==0)|(diff$V1=='translating' & diff$V2=='translating'),'no change',ifelse(diff$V1==0 & diff$V2=='translating','turn-on','turn-off'))
tra2a$Grade<-ifelse(rownames(tra2a) %in% who1$CGGA_ID,'WHO II',
                ifelse(rownames(tra2a) %in% who2$CGGA_ID,'WHO III','WHO IV')) #活用可以精简代码，可以节省体力


########我一直想写一个自动画图的function，以下为第一次尝试##########
plot <- function(data,gene,group,label_1){
  data2<-as.data.frame(t(data[gene,]))
  colnames(data2)<-c('exp')
  data2$group<-group
  data2$group<- factor(data2$group,levels = c("normal", "tumor"))
  library(ggsignif)
  library(ggsci)
  library(ggpubr)
  ptext<-paste0("t.test, p=",compare_means(exp~group,data2,method = "t.test")[[5]])
  dat_text <- data.frame(label = ptext)
  dat_text[1,1]<-c(paste0(dat_text[1,1],'\n Tumor N=157 \n Paratumor N=157'))
  colors<-c('#36537155', "#96345A74")
  
  dp <- ggplot(data2, aes(x=group, y=exp,fill=group))+
    geom_violin(aes(colour=group),fill='#DDDDDD50',trim=FALSE,show.legend=FALSE,size=0.7)+
    geom_jitter(aes(fill=group),width =0.2,shape = 21,size=4,colour='NA',alpha = 0.4)+
    geom_boxplot(fill='white',width=0.12,,outlier.shape = NA,size=0.8)+
    labs(y="log2(FPKM)")
  path<-c(paste0(gene,'.pdf'))
  pdf(path,width = 8,height = 7)
  print(dp + theme_classic(base_size = 22)+ 
    scale_colour_manual(values=c(colors))+
    scale_fill_manual(values=c(colors))+
    theme(legend.position = "none")+
    theme(axis.text=element_text(size=18),axis.title=element_text(size=22))+
    annotate('text',x=0.8,y=label_1,label=dat_text,colour = "black",size=6))
  dev.off()}

plot(data = exp,gene = 'ENSG00000123009.4_120282483_120282716_77',group = group_list,label_1 =1.4)
#改了一下，我觉得好歹算个半自动画图函数，主要在文本的地方改了一下，现在不要老去调两个文本的位置了，方便了很多。针对不同的矩阵，还是需要调整许多地方，但是这种程度的调整的话，都还能接受吧。



########双循环和if的判断############
for (i in 1:nrow(cortest)){
  for (j in 1:ncol(cortest)){
    if(fdr[i,j]< 0.05) {cortest[i,j]<-cortest[i,j]} else {cortest[i,j]<-c(0)}
  }
}

for (i in 1:nrow(cortest)){
  for (j in 1:ncol(cortest)){
    if(abs(cortest[i,j])> 0.3) {cortest[i,j]<-cortest[i,j]} else {cortest[i,j]<-c(0)}
  }
}           #写得很丑，for循环里面叠一个for循环就是很呆的写法，但是暂时管不了那么多了，随便吧。if的判断倒是之前没写过的代码，其实还挺容易理解的。



###############if判断句############
seq_list$jud<-c(NA)
for (i in 1:nrow(seq_list)){
test<-unlist(strsplit(seq_list[i,2],'')) 
jud<-if(isEmpty(grep('[*]',test))) {isEmpty(grep('[*]',test))} else {identical(grep('[*]',test),length(test))}
seq_list[i,3]<-c(jud)}  #能不能把这个改成apply循环呢？感觉for循环还是太慢了



##############给for循环加一个进度条###########
library(progress)

test<-as.data.frame(array(NA,c(100,2)))
test$V1<-c(1:100)

pb <- progress_bar$new(total = 100)
for (i in 1:100) {
  pb$tick()
  test[i,2]<-c((test[i,1]*3/4+19)/3.1976)
  Sys.sleep(1 / 100)
}
#不加Sys.sleep的话，就看不到进度条了不知道为什么，可能如果运行太快的话不会出现？


#############一个示例,批量的计算基因对癌症生存的影响#################
library(survminer)
library(survival)
library(ggpubr)
library(gridExtra)

for (i in 1:nrow(up)){
 gene <- substr(unlist(strsplit(up[i,1],'_'))[1],1,15)
 if(length(grep('NONHSAG',gene))>0) {up[i,20]<-c('NA')} else{    #因为含有NONCODE来源的基因，所以在这里写了一个判断句，如果遇到NONCODE基因则输出“NA”
   gene_exp<-surexp[rownames(surexp)==gene,]                     #else里面才是正文，
   gene_sur<-survival
   gene_sur$exp<-c(as.numeric(gene_exp[1,]))
   gene_sur<-gene_sur[gene_sur$exp>0,]
   
   gene_sur<-gene_sur[!is.na(gene_sur$OS)&
                        !is.na(gene_sur$OS.time),]
   matt<-gene_sur
   med.exp<-quantile(matt$exp)
   more.med.exp.index<-which(matt$exp>=med.exp[4])                #依据四分位数判断表达的高低
   less.med.exp.index<-which(matt$exp<=med.exp[2])
   matt$status<-NA
   matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
   matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')
   
   surv_diff <- survdiff(Surv(OS.time,OS) ~ status, data = matt)
   p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
   up[i,20]<-c(p.value)
 }
}

###################然后我发现对大量基因进行分析时会有几种情况导致报错##################
for (i in 1:nrow(ms)){
 gene <- substr(unlist(strsplit(ms[i,1],'_'))[1],1,15)
 if(length(grep('NONHSAG',gene))>0) {ms[i,20]<-c(NA)} else{
   if(length(grep(gene,rownames(surexp)))<1) {ms[i,20]<-c(NA)} else{   #一，这个基因太新了，在TCGA的矩阵中没有出现
   gene_exp<-surexp[rownames(surexp)==gene,]
   gene_sur<-survival
   gene_sur$exp<-c(as.numeric(gene_exp[1,]))
   gene_sur<-gene_sur[gene_sur$exp>0,]
   
   if (nrow(gene_sur)<4) {ms[i,20]<-c(NA)} else{                       #二，这个基因转录本表达太少了，不足以进行分析
   gene_sur<-gene_sur[!is.na(gene_sur$OS)&
                        !is.na(gene_sur$OS.time),]
   matt<-gene_sur
   med.exp<-quantile(matt$exp)
   more.med.exp.index<-which(matt$exp>=med.exp[4])
   less.med.exp.index<-which(matt$exp<=med.exp[2])
   matt$status<-NA
   matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
   matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')
   
   surv_diff <- survdiff(Surv(OS.time,OS) ~ status, data = matt)
   p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
   ms[i,20]<-c(p.value)
 }}}
}    #考虑了这几种情况，又多加了几次判断，水多了加面面多了加水，写这类的代码差不多就是这种感觉

#########################写得好丑这样的代码，完了还有五十行报错，我真是##########
for (i in 1:nrow(allindex)){
  name_bef<-allindex[i,1]
  if(length(grep('CNC_N',name_bef))>0) { #如果来自CNC
    if(length(unlist(strsplit(name_bef,'_')))==4){
      source<-paste0(unlist(strsplit(name_bef,'_'))[2],'_',unlist(strsplit(name_bef,'_'))[3],'_',unlist(strsplit(name_bef,'_'))[4])
      inf<-cncdb[grep(source,cncdb$cncRNAdb.ID),]
      allindex$adj_name[i]<-paste0(unlist(strsplit(name_bef,'_'))[1],'_',inf$Start,'_',inf$End,'_',source)
      } else{source<-paste0(unlist(strsplit(name_bef,'_'))[3],'_',unlist(strsplit(name_bef,'_'))[4],'_',unlist(strsplit(name_bef,'_'))[5])
      inf<-cncdb[grep(source,cncdb$cncRNAdb.ID),]
      allindex$adj_name[i]<-paste0(unlist(strsplit(name_bef,'_'))[1],'-',unlist(strsplit(name_bef,'_'))[2],'_',inf$Start,'_',inf$End,'_',source)
    }
  } else {
    if(length(grep('SPROHSA',name_bef))>0) { #如果来自smprot
      source<-unlist(strsplit(name_bef,'_'))[2]
      inf<-smprotdb[smprotdb$peak==source,]
      allindex$adj_name[i]<-paste0(unlist(strsplit(name_bef,'_'))[1],'_',inf$start,'_',inf$end,'_',source)
    } else{
      allindex$adj_name[i]<-allindex$name[i]
    }
  }
}


##################################################################################
#写了两个function比较以下fread和vroom的读取速度

findribofile<-function(gene){
  
  list<-fread('~/share/ORFfind/ribotricer_genecode/cut/list.txt')
  cache<-fread('~/share/ORFfind/ribotricer_genecode/cut/SRR10416857_cut.tsv')[1,]
  cache$sample<-c(NA)
  
for (i in list$list.txt) {
  path<-paste0('~/share/ORFfind/ribotricer_genecode/cut/',i)
  ri<-fread(path)
  ri<-ri[ri$gene_id==gene,]
  ri$sample<-c(i)
  cache<-rbind(cache,ri)
  print(paste0('comp ',round((which(list$list.txt==i))/length(list$list.txt),2)*100,'%'))
}
  return(cache[-1,])
}

system.time(test<-findribofile('ENSG00000268034.1')) #user   system  elapsed 1502.012   32.815  768.757

findribofile2<-function(gene){
  
  list<-vroom('~/share/ORFfind/ribotricer_genecode/cut/list.txt',delim = "\t")
  cache<-vroom('~/share/ORFfind/ribotricer_genecode/cut/SRR10416857_cut.tsv',delim = "\t")[1,]
  cache$sample<-c(NA)
  
  for (i in list$list.txt) {
    path<-paste0('~/share/ORFfind/ribotricer_genecode/cut/',i)
    ri<-vroom(path,delim = "\t")
    ri<-ri[ri$gene_id==gene,]
    ri$sample<-c(i)
    cache<-rbind(cache,ri)
    print(paste0('comp ',round((which(list$list.txt==i))/length(list$list.txt),2)*100,'%'))
  }
  return(cache[-1,])
}

system.time(test<-findribofile2('ENSG00000268034.1'))#user   system  elapsed  1367.137   33.669 1038.510 

#用fread读取的话，340个文件，比vroom快一些。

####################################################################################################
#function里面叠function，究极套娃

plot_ribo_track<-function(gene,gene_name,transcript_id){
  
print('search_sample')
AC093673.1<-findribofile(gene) #懒得改这个中间文件的名字，有时间再说
AC093673.1<-AC093673.1[order(AC093673.1$read_density,decreasing=T),]
AC093673.1<-AC093673.1[!duplicated(AC093673.1$sample),]

bgfile<-dir('/home/leelee/share/mapping/genome/')[grep('bedgraph',dir('/home/leelee/share/mapping/genome/'))]
bgfile<-bgfile[-grep('_uniq_',bgfile)]
bgfile<-bgfile[-grep('_P_sites_',bgfile)]

f <-function(x) unlist(strsplit(x['sample'],'_'))[1]
AC093673.1$sample<-apply(AC093673.1,1,f)

if(AC093673.1$strand[1]=='+') {strand<-'plus'} else{strand<-'minus'} #判断基因的正负链

AC093673.1$adj_sample<-paste0(AC093673.1$sample,'Aligned.sortedByCoord.out.bam_coverage_',strand,'.bedgraph')
AC093673.1<-AC093673.1[AC093673.1$adj_sample %in% bgfile,]

print('treat_file')
if(nrow(AC093673.1)>5) {
  AC093673.1<-AC093673.1[1:5,]
  } else{
  AC093673.1<-AC093673.1
  }

for (i in AC093673.1$sample){
  treat<-paste0('LC_COLLATE=C sort -k1,1 -k2,2n /home/leelee/share/mapping/genome/',
                i,'Aligned.sortedByCoord.out.bam_coverage_',strand,'.bedgraph > ',
                '/home/leelee/share/mapping/bw/',i,'.',strand,'.sorted.bedGraph')
  system(treat)
  treat2<-paste0('/home/leelee/miniconda3/envs/p3/bin/bedGraphToBigWig /home/leelee/share/mapping/bw/',
                 i,'.',strand,'.sorted.bedGraph',' /home/leelee/biodata/index/GRCh38/GRCh38.info.txt ',
                 '/home/leelee/share/mapping/bw/',i,'.',strand,'.bw')
  system(treat2)
}

file <- c(paste0('~/share/mapping/bw/',AC093673.1$sample,'.',strand,'.bw'))
samp <- c(AC093673.1$sample)
allBw <- loadBWfile(file = file, sample = samp)

print('start plot')

plot_path<-paste0(gene_name,'ribo_track.pdf')
pdf(plot_path,width = 8,height = 8)
print(plotTrack(
  gtfFile = gtf,
  gene = gene_name,
  arrowCol = "black",
  bigwigFile = allBw,
  sampleAes = "group1",
  facetVars = "sample",
  multiple = TRUE,
  myTransId = c(transcript_id)
  ))
dev.off()

}

#############################################################################################################

#每次富集都要写很长一段的重复代码，包装一下会方便一些？这个function我还没有试验过
                     
myenrich <- function(genelist,name,dir,species){
  
  library(org.Mm.eg.db) #据说library不会重复载包，那我就没有理由把library放在外面了
  library('org.Hs.eg.db')
  library(patchwork)
  library(clusterProfiler)
  library(stringr)
  library(DOSE)
  library(openxlsx)
  library(ggplot2)
  R.utils::setOption("clusterProfiler.download.method",'auto')
  
  p<-list()
  if (species=='human'){organism='hsa' 
  Org = org.Hs.eg.db
  } else {organism='mmu' 
  Org = org.Mm.eg.db
  } 
  
  genelist<-unique(genelist)
  kk.negative <- enrichKEGG(gene  = genelist,
                            organism = organism,
                            pvalueCutoff = 0.5,
                            qvalueCutoff =0.5)
  kk.negative<-setReadable(kk.negative,OrgDb = Org, keyType="ENTREZID")
  kk=kk.negative@result
  kk<-kk[order(kk$pvalue,decreasing = F),]
  
  library('stringr')
  title<-paste0('KEGG')
  color<-c(MetBrewer::met.brewer('Hokusai1')[c(3,5,6,7)],
           MetBrewer::met.brewer('Troy')[5])
  
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
    labs(y="-log10(P value)",x=' ')+
    scale_fill_manual(values = c(color[1]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))+
    theme(aspect.ratio=1.4)
  
  BP <- enrichGO(gene          = genelist,
                 OrgDb         = Org,
                 ont           = 'BP' ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.8,
                 qvalueCutoff  = 0.8,
                 readable      = TRUE)
  
  CC <- enrichGO(gene          = genelist,
                 OrgDb         = Org,
                 ont           = 'CC' ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.8,
                 qvalueCutoff  = 0.8,
                 readable      = TRUE)
  
  MF <- enrichGO(gene          = genelist,
                 OrgDb         = Org,
                 ont           = 'MF' ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.8,
                 qvalueCutoff  = 0.8,
                 readable      = TRUE)
  
  BP=BP@result
  BP<-BP[order(BP$pvalue,decreasing = F),]
  CC=CC@result
  CC<-CC[order(CC$pvalue,decreasing = F),]
  MF=MF@result
  MF<-MF[order(MF$pvalue,decreasing = F),]
  
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
    labs(y="-log10(P value)",x=' ')+
    scale_fill_manual(values = c(color[2]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))+
    theme(aspect.ratio=1.4)
  
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
    labs(y="-log10(P value)",x=' ')+
    scale_fill_manual(values = c(color[3]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))+
    theme(aspect.ratio=1.4)
  
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
    labs(y="-log10(P value)",x=' ')+
    scale_fill_manual(values = c(color[4]))+
    scale_colour_manual(values = c("black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 45))+
    theme(plot.title = element_text(size=12))+
    theme(aspect.ratio=1.4)
  
  path<-paste0(dir,'/',name,'_enrich.pdf')
  pdf(path,width = 12,height = 8)
  print(wrap_plots(p,nrow=2)) 
  dev.off()
  
  path<-paste0(dir,'/',name,'_enrich.xlsx')
  sheets = list("KEGG" = kk,"GO_BP" = BP,'GO_CC'=CC,'GO_MF'=MF) #,'DO'=DO
  write.xlsx(sheets,path)
}

#试了一下是可以正常运行的，使用方法就是myenrich(rnacotarget,'RNA_cotarget','enrich/','human')
#我是不是改自己写一个包什么的会比function更加方便？但是function方便修改，写包的话就不太方便修改代码了。
#这个要说有什么想要优化的，可以把bitr的部分也包进去，要求的genelist只要SYMBOL ID就好了
                     
                     
##########################################################                     
#给for循环加个进度条
library(progress)
pb <- progress_bar$new(total = nrow(shTRA2A_down))

for (i in 1:nrow(shTRA2A_down)){
  pb$tick()
  chr<-traclip[traclip$V1==shTRA2A_down$seqnames[i]& 
                 traclip$V4==shTRA2A_down$feature_strand[i],]
  chr<-chr[(chr$V2>shTRA2A_down$start[i] & chr$V3<shTRA2A_down$end[i])| 
             (chr $V2< shTRA2A_down$start[i] & chr $V3> shTRA2A_down$start[i]) |
                (chr $V2 < shTRA2A_down$end[i] & chr$V3> shTRA2A_down$end[i]),]
  if (nrow(chr)>0 ) {shTRA2A_down$group[i]<-c('target')} else {
    shTRA2A_down$group[i]<-c('nontarget')
  }
}
                     
