###########################################
#测试IsoformSwitchAnalyzeR的使用，尝试一个demo数据集

library(IsoformSwitchAnalyzeR)
library(data.table)

v1<-fread('demo/pag057v1.transcript_counts.tsv') 
v1$c2<-round(v1$count * (1 + rnorm(nrow(v1), mean = 0, sd = 0.2))) #来自isoquant结果，因为只有一个样本所以用随机偏移做了一个人工样本

v2<-fread('demo/pag057v3.transcript_counts.tsv')
v2$c2<-round(v2$count * (1 + rnorm(nrow(v2), mean = 0, sd = 0.2)))
count<-merge(v1,v2,by=c('#feature_id'))
colnames(count)<-c('isoform_id','A1','A2','B1','B2')
count<-count[rowSums(count[,2:5])>0,]
count<-count[grep('ENST',count$isoform_id),]

v1<-fread('demo/pag057v1.transcript_tpm.tsv')
v1$c2<-v1$TPM * (1 + rnorm(nrow(v1), mean = 0, sd = 0.2))

v2<-fread('demo/pag057v3.transcript_tpm.tsv')
v2$c2<-v2$TPM * (1 + rnorm(nrow(v2), mean = 0, sd = 0.2))
TPM<-merge(v1,v2,by=c('#feature_id'))
colnames(TPM)<-c('isoform_id','A1','A2','B1','B2')

library(rtracklayer)

gtf <- import("demo/pag057v3.transcript_models.gtf", format = "gtf")
gtf <- gtf[gtf$transcript_id %in% count$isoform_id ,]
export(gtf, "demo/modified_output.gtf", format = "gtf")
#构建对象时一定要gtf文件和定量文件匹配，因为我在跑isoquant时没有把两个样本放在一起所以导致了注释文件有些差异，这里只是做为demo所以强行做了匹配

count<-count[count$isoform_id %in% gtf$transcript_id,]
TPM<-TPM[TPM$isoform_id %in% count$isoform_id,]


myDesign <- data.frame(
  sampleID = c('A1','A2','B1','B2'),
  condition = c('A','A','B','B')
)
myDesign #分组设计

aSwitchList <- importRdata(
  isoformCountMatrix   = as.data.frame(count),
  isoformRepExpression = as.data.frame(TPM),
  designMatrix         = myDesign,
  isoformExonAnnoation = 'demo/modified_output.gtf',
#  isoformNtFasta       = system.file("extdata/example_isoform_nt.fasta.gz", package="IsoformSwitchAnalyzeR"), #没有这个水平的序列，但并不是必需的
#  fixStringTieAnnotationProblem = TRUE,
  showProgress = T,
  ignoreAfterBar =T
)

#这一步应该默认把ORF加上去了

summary(aSwitchList)

mySwitchList <- preFilter( aSwitchList )

mySwitchList <- isoformSwitchTestDEXSeq( mySwitchList ) #使用DEXSeq做转录本isoform的差异分析

library(BSgenome.Hsapiens.UCSC.hg38)
aaseq <- extractSequence( mySwitchList ,writeToFile=T , genomeObject = Hsapiens ) #这里可以直接从genomeObject得到ORF序列

head(aaseq$ntSequence,2)

#mySwitchList <- analyzePFAM( mySwitchList ,pathToPFAMresultFile=system.file("extdata/pfam_results.txt", package="IsoformSwitchAnalyzeR") )
#PFAM的网页和他们在文档中描述的完全不一样了。。。暂时不知道要怎么得到这个结果

mySwitchList <- analyzeAlternativeSplicing( mySwitchList )

consequencesOfInterest <- c(
  'intron_retention',
  'NMD_status'
) #可选的特征是很多的，但大多依赖于第三方工具的结果

mySwitchList <- analyzeSwitchConsequences( mySwitchList,consequencesToAnalyze = consequencesOfInterest )

switchPlotTopSwitches( mySwitchList )

### Summary
extractSwitchSummary(mySwitchList, filterForConsequences = TRUE)

subset(
  extractTopSwitches(
    mySwitchList,filterForConsequences = TRUE,n=10,
    inEachComparison = TRUE)[,c('gene_name','condition_1','condition_2','gene_switch_q_value','Rank')])

#展示其中一个
switchPlot(
  mySwitchList,
  gene='SRSF1',
  condition1 = 'A',
  condition2 = 'B',
  localTheme = theme_bw(base_size = 13)
)


# extractSwitchOverlap(
#   mySwitchList,
#   filterForConsequences=TRUE,
#   plotIsoforms = FALSE
# ) #如果有多于两组样本可以用这个代码展示

extractConsequenceSummary(
  mySwitchList,
  consequencesToAnalyze='all',
  plotGenes = FALSE,
  asFractionTotal = FALSE
)


extractConsequenceEnrichment(
  mySwitchList,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  localTheme = theme_bw(base_size = 14), 
  returnResult = FALSE
)

# extractConsequenceEnrichmentComparison(
#   mySwitchList,
#   consequencesToAnalyze='all',
#   analysisOppositeConsequence = TRUE,
#   returnResult = FALSE 
# )


# extractSubCellShifts(
#   mySwitchList,
#   returnResult = FALSE 
# ) #需要deeploc结果

extractSplicingEnrichment(
  mySwitchList,
  returnResult = FALSE 
)

# extractSplicingEnrichmentComparison(
#   mySwitchList,
#   splicingToAnalyze = c('A3','A5','ATSS','ATTS'), 
#   returnResult = FALSE
# )


ggplot(data=mySwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()

ggplot(data=mySwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  facet_wrap(~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw()


#我没有使用novel转录本相关的功能，如果做了de nova的组装的话也是可以拓展分析的，总的来说适用性还算广

