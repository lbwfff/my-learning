python<-c('/home/leelee/miniconda3/envs/pepnn/bin/python') #conda环境内的python
code<-c('/home/leelee/tools/pepnn/pepnn_seq/scripts/predict_binding_site.py') #pepnn的python代码位置

#code<-c('/home/leelee/tools/pepnn/pepnn_struct/scripts/predict_binding_site.py') #根据pdb结构文件预测的话可以用这个代码
# prot<-c('/home/leelee/tools/pepnn/example/test_drugbank/target.fasta')

pep<-c('/home/leelee/tools/pepnn/example/test_drugbank/pep.fasta') #你的肽序列，这里我谁便找的一个
out<-c('/home/leelee/tools/pepnn/example/test_drugbank/cache') #输出文件夹

system.time(system(paste0(python,' ',code,' ','-prot',' ',prot,' -pep ',pep,' -o ',out))) #试了一下是可行的，很久没写过shell脚本了，用R习惯了，喜欢把命令行写在R里

#先需要把fasta文件拆分
humanunip<-readAAStringSet('~/biodata/index/protein/humanuniport.fasta') #人类源的uniprot蛋白序列
for (i in 1:20) {
  prot<-humanunip[i]
  path<-(paste0('/home/leelee/tools/pepnn/example/test_drugbank/prot_fasta/',unlist(strsplit(prot@ranges@NAMES,'\\|'))[2],'.fasta'))
  writeXStringSet(prot,filepath = path)
  } #这里只拆了前二十个，主要试用一下

pepnnarar<-as.data.frame(array(NA,c(20,10)))
colnames(pepnnarar)[c(1,2,3,4,5)]<-c('pep','prot','prm_score','mean_probabilities','freq_more9') #做了很多栏，一时没想到后续栏要放什么

for (i in 1:20) {
  prot<-(paste0('/home/leelee/tools/pepnn/example/test_drugbank/prot_fasta/',unlist(strsplit(humanunip@ranges@NAMES[i],'\\|'))[2],'.fasta'))
  system(paste0(python,' ',code,' ','-prot',' ',prot,' -pep ',pep,' -o ',out))
  result<-read.table('/home/leelee/tools/pepnn/example/test_drugbank/cache/prm_score.txt')
  pepnnarar$prot[i]<-unlist(strsplit(humanunip@ranges@NAMES[i],'\\|'))[2]
  pepnnarar$prm_score[i]<-result$V6
  result<-read.csv('/home/leelee/tools/pepnn/example/test_drugbank/cache/binding_site_prediction.csv')
  pepnnarar$mean_probabilities[i]<-mean(result$Probabilities)
  pepnnarar$freq_more9[i]<-c(length(which(result$Probabilities>0.9))/nrow(result))
} #运行，读取结果

#pepnn这个工具，我先不说准不准的问题，起码在可用性上，它可好太多了，安装方便，计算简单，结果明显
#但是想用这个工具做肽蛋白结合的预测的话，可能自己结合它的输出结果再做一些简单的机器学习会好一些，简单的阈值判断的话，不确定是不是合理
