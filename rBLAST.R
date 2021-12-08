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


