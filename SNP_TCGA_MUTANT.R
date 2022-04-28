####maftools是一个很cool的包，能够分析TCGA提供的患者体细胞突变数据，一步到位，出图也很好看，虽然后来没有用到吧，感觉挺有意义的一次分析##########

clin = read.table("clinical.cart.2022-04-25/clinical.cart/clinical.tsv",header = T,sep = "\t",fill=T)

clin1 <- clin[,c("case_id","case_submitter_id","vital_status","days_to_death",
                 "ajcc_pathologic_m","ajcc_pathologic_n",
                 "ajcc_pathologic_stage","ajcc_pathologic_t")]
colnames(clin1) = c("case_id","Tumor_Sample_Barcode","vital_status","days_to_last_follow_up",
                    "ajcc_pathologic_m","ajcc_pathologic_n",
                    "ajcc_pathologic_stage","ajcc_pathologic_t")
#clin1 <- clin1[clin1$days_to_last_follow_up!= "--",]
clin1$days_to_last_follow_up <- as.numeric(clin1$days_to_last_follow_up)
#clin1 <- clin1[clin1$vital_status!= "--",]
#clin1$vital_status <- ifelse(clin1$vital_status== "Alive",0,1)

################################
dir.create("SampleFiles")
filepath<-dir(path = 'gdc_download_20220425_030104/',full.names=T)
for(wd in filepath){
  files <-dir(path = wd,pattern="gz$")
  fromfilepath <- paste(wd,"\\",files,sep="")
  tofilepath <- paste(".\\SampleFiles\\",files,sep="")
  file.copy(fromfilepath,tofilepath)
}

setwd(".\\SampleFiles")
countsFiles<-dir(path=".\\",pattern="gz$")
library(R.utils)
sapply(countsFiles,gunzip)

files <-dir(path =getwd(),pattern="maf$")

cache<-read.table(file = '019278a6-0d1a-4239-b006-d10dc5d786bd.wxs.aliquot_ensemble_masked.maf',sep = '\t',header = T,
                fill = T,comment.char = "#",quote = "")
cache<-cache[1,]
for (i in 1:length(files)){
  file_path<-files[i]
  maf<-read.table(file = file_path,sep = '\t',header = T,
                  fill = T,comment.char = "#",quote = "")
  cache<-rbind(cache,maf)
}
cache<-cache[-1,]

library(tidyverse)
maf.new = cache %>% 
  rename(t_ref_count = Ref_allele_depth, 
         t_alt_count = Alt_allele_depth )

write.table(cache,file = '..\\TCGA_LIHC.maf',sep = '\t',row.names = F,
            quote = F)
#

library(maftools)
laml = read.maf(maf = cache,
                clinicalData = clin1,
                isTCGA = TRUE)

plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = laml, top = 10)


#################################我的SNP数据是根据uniprot做的，它提供了一些整合的数据###############
humsav<-read.table('humsavar.txt',sep = '\t',header = T,fill = T,quote = '') #他这个文件格式给得莫名其妙，所以花了一些代码来调整表格

humsav$Main_gene_name<-gsub("[ ]",";",humsav$Main_gene_name)
humsav$gene_name<-c(NA)

for (i in 1:nrow(humsav)){
all<-unlist(strsplit(humsav[i,1],';'))
all<-all[nchar(all)>0]
humsav[i,2]<-all[2]
humsav[i,3]<-all[3]
humsav[i,4]<-all[4]
humsav[i,5]<-all[5]
humsav[i,6]<-all[6]
humsav[i,7]<-all[7]
humsav[i,8]<-all[1]
}

test_1<-humsav[-grep('rs',humsav$dbSNP),]
test_2<-humsav[grep('rs',humsav$dbSNP),]

test_3<-test_1[-grep('rs',test_1$Main_gene_name),]
test_4<-test_1[grep('rs',test_1$Main_gene_name),]
test_4<-test_4[grep('rs',test_4$Variant.category),]

write.csv(test_4,file = 'test_4.csv',row.names = F)
test_4<-read.csv('test_4.csv')

humsav<-rbind(test_2,test_4)
save(humsav,file = 'humsav.RData')

humsav$pre_AA<-substr(humsav$AAchange,3,5)
library(Biostrings)
prot<-readAAStringSet('uniprot_SNP.fasta')

humsav<-humsav[humsav$pre_AA!='Sec',]
aa<-read.csv('AA.csv')
aa<-aa[match(humsav$pre_AA,aa$tri_a),]
humsav$pre_AA<-aa$sig_a

library('stringr')
humsav$aft_AA<-substr(humsav$AAchange,str_length(humsav$AAchange)-2,str_length(humsav$AAchange))
humsav$site<-substr(humsav$AAchange,6,str_length(humsav$AAchange)-3)
aa<-read.csv('AA.csv')
aa<-aa[match(humsav$aft_AA,aa$tri_a),]
humsav$aft_AA<-aa$sig_a
humsav$adj<-c(NA)
humsav$site<-as.numeric(humsav$site)

library('seqinr')

for (i in 1:nrow(humsav)){
  SNPid<-humsav[i,6]
  seq<-prot[grep(humsav[i,2],prot@ranges@NAMES)]
  seq1<-substr(gsub("\\.","",paste(seq)),1,humsav[i,11]-1)
  humsav[i,12]<-if(substr(gsub("\\.","",paste(seq)),humsav[i,11],humsav[i,11])==humsav[i,9]){'OK'} else {'MATTER'} #用一个判断句验证了一些改变的原氨基酸是否一致
  seq2<-humsav$aft_AA[i]
  seq3<-substr(gsub("\\.","",paste(seq)),humsav[i,11]+1,str_length(gsub("\\.","",paste(seq))))
  seq_new<-c(paste0(seq1,seq2,seq3))
  file_path<-(paste0('seq/',SNPid,'.fasta'))
  write.fasta(sequences=seq_new, names=SNPid, file.out=file_path, open='w')
}
detach("package:seqinr", unload = TRUE) 
