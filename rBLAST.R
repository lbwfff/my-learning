#很简单的一个包

library(rBLAST)
library(Biostrings)

test<-readAAStringSet('~/biodata/index/protein/ncRNA-coding/cncRNA-uni.fasta')

Sys.setenv(PATH ="/home/leelee/tools/CAMP/SCRATCH-1D_2.0/opt/blast+_2.10.1/bin/")#设置环境
system("blastn -version")#测试是否有效

db <- makeblastdb('~/biodata/index/protein/humanuniport.fasta', dbtype = "prot", args="")#建立索引
bl <- blast(db='~/biodata/index/protein/humanuniport.fasta',type = 'blastp')#注意type，如果是核苷酸索引的话需要改成blastn
bl
output <- predict(bl, test) #进行运算

#我上面把环境位置改了，我发现不改回来的话，会没法装包，所以用完记得把环境改回来
Sys.setenv(PATH ='/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin:/usr/lib/rstudio-server/bin/postback:/usr/lib/rstudio-server/bin/postback')

#不在R上做的话#
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out Swiss-Prot #创建索引
blastp -query test.fasta -outfmt 7 -out test.blast -db ~/biodata/index/protein/humanuniport.fasta #outfmt代表控制文件的格式



#记录一下怎么创建一个fasta文件（在R上）
whatever<-read.csv('ALL_MY_PEPTIDE_EXP_2.CSV')
library('seqinr')
seq<-whatever[,c('name','Sequence')] #名字列和序列列
characters<-seq[,c("Sequence")]
list<-as.list(characters) #需要有这么一个转换才能输出fasta文件的序列
write.fasta(sequences=list, names=seq$name, file.out='pep_to_blast.fasta', open='w', nbchar=60)
detach("package:seqinr", unload = TRUE) #seqinr某个函数和Biostrings重名了，使用的时候会有冲突，我习惯用完之后detach seqinr

#和blast无关，记录一下怎么从一个biostring对象中提取序列
sequence = gsub("\\.","",paste(test)) #test就是这个biostring对象，没有这一行还真不知道要怎么拿出这些序列来。
