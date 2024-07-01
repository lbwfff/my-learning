#狠补了一下MR的分析流程，虽然不知道未来是不是用得到
#这个流程要比之前的精细得多，关注的也不是临床表征和临床表征的关系，而是在做基因表达调控

#对于eQTL或是sQTL什么的粗筛，使用smr工具进行分析，
/scratch/lb4489/project/GWAS/smr-1.3.1-linux-x86_64/smr --bfile /scratch/lb4489/project/GWAS/csmr/data/1kg/EUR \
      --gwas-summary /scratch/lb4489/project/GWAS/csmr/data/finn_adj  \
      --beqtl-summary /scratch/lb4489/project/GWAS/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr1 \
      --out testchr1 --thread-num 10 
#这个工具是专门针对于QTL的孟德尔随机化分析软件，会比twosample要方便，重点是快得多

#但是这个工具也存在许多的限制，比如说我想做genetic colocalisation就没法实现，然后也缺少一些精细的质量控制结果，敏感性异质性什么的

#################################
#基于R语言的finemapping和colocalisation

micro<-data.table::fread('./single_cell/Endothelial.cells.1.gz')
micro<-micro[grep('ENSG00000238009',micro$V1),]
match<-snpinf[match(micro$V2,snpinf$SNP),]
micro$position<-match$SNP_id_hg38
micro$effect_allele<-match$effect_allele
micro$other_allele<-match$other_allele
micro$MAF<-match$MAF
colnames(micro)[1:5]<-c('Gene_id','SNP_id','Distance','p-value','beta') #读取一个示例文件，得到某一个基因的eQTL

D2<-micro
D2$position<-substr(D2$position,6,99)

D1<-outcome_dat[outcome_dat$chrom=='1',] #outcome文件中对应的SNP信息

merge<-merge(D2,D1,by.x=c('SNP_id'),by.y=c('rsids'),suffixes=c('.eqtl','.gwas')) #合并结果
merge<-merge[!duplicated(merge$SNP_id),] #会有一些重复的SNP，似乎在同一个位置但是有不同的allele，不太确定要怎么合并，简单的直接去除了
merge$gwas.se<-abs(merge$beta.eqtl/qnorm(merge$`p-value`/2)) #计算se值

D2<-list(type='quant',beta=merge$beta.eqtl,snp=merge$SNP_id,varbeta=c(merge$gwas.se),sdY=1, #varbeta为se的平方，某些数据里sdY为1，不大明白意味着什么
         position=merge$position,MAF=merge$MAF)

D1<-list(type='cc',beta=merge$beta.gwas,varbeta=(merge$sebeta)^2,snp=merge$SNP_id,
         position=merge$pos,N=412181,MAF=merge$af_alt)

check_dataset(D1) #检查对象的完整性
check_dataset(D2)

my.res <- coloc.abf(dataset1=D1,
                    dataset2=D2) #分析

#colocalisation分析时：1不应该对D1D2进行任何筛选

sensitivity(my.res,rule="H4 > 0.01") #这个图还挺有意思的，要怎么解读呢？

library(dplyr)
library(geni.plots)

assoc<-merge[,c('SNP_id','chrom','position','pval','p-value')]
colnames(assoc)<-c('marker','chr','pos','pvalue_1','pvalue_2')

corr<-ieugwasr::ld_matrix(merge$SNP_id, with_alleles = F, 
                          plink_bin = "/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/plinkbinr/bin/plink_Darwin",
                          bfile = "./clump/1kg.v3/EUR") #计算LD matrix，找了好久才找到的实现方法

assoc<-assoc[assoc$marker %in% rownames(corr),]
assoc<-assoc[match(rownames(corr),assoc$marker),]
assoc$pos<-as.integer(assoc$pos)

fig_region_stack(
  data = assoc,
  traits = c("PTSD", "eQTL"),
  corr = corr,
  build = 38,
  highlights = "rs116587930",
  title_center = TRUE) #然后就可以画另一个经典的图了，虽然看不懂是什么

############################
#fine-mapping

library(susieR)

finemap<-merge[match(rownames(corr),merge$SNP_id),]

D1<-list(type='cc',beta=finemap$beta.gwas,varbeta=(finemap$sebeta)^2,snp=finemap$SNP_id,
         position=finemap$pos,N=412181,MAF=finemap$af_alt,LD=corr) 
#其实就是在coloc输入对象的基础上加上了LD matrix
#这样的话fine-mapping就变成了coloc前置条件了

check_dataset(D1)

finemapresult<-runsusie(D1)

summary(finemapresult) #finemapping结果

susie_plot(finemapresult, y="PIP", b=b)

res=coloc.susie(S5,S4) #可以对两个runsusie的输出对象直接coloc

#但这样做的话不就纯为coloc服务了吗，要怎么用它来筛选IV呢？


