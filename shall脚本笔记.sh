###########对于学生物得自己来说可能不需要掌握太高端得脚本知识##############
vim xxx.sh #####创建脚本文件
#!bin/bash  ###注意这一行不是注释，是有意义的语句
#之后把你之前一句一句写的语句放上来
保存退出
chmod 777 xxx.sh ###赋予权限
sh xxx.sh ###运行
写了两句cutadapt，是能用的，那就说明你在生信分析中要用到的所有软件都可以用这样的方法来运行，如果可以的话，虽然不够简介而潇洒，但是完全够用了。
如果是在windows下写得脚本将出现格式问题，在vim编辑器下set fileformat=unix便可以解决问题

后来就没有给过权限了，不过还是一样的用。

除了上面一句一句的写外，现在用的方法就两种。
一种是，类似于
find *gz |while read id; do (trim_galore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -j 10 -length 30 -stringency 3 ${id} );done
这种，while read id ; do 的句子是第一个学会的，比较好理解，写起来也非常的简洁，但在某些状态下，例如处理双端测序的数据时就会不知道要怎么做。

还有一种是，类似于
for i in Input-1 Input-2 Abcam-1 Abcam-2 CST-1 CST-2 SYSY-1 SYSY-2
do
hisat2 -p 18 --dta -t -x /home/leelee/biodata/index/GRCH38/genome -1 "$i"-trimmed.derRNA.fq.1.gz -2 "$i"-trimmed.derRNA.fq.2.gz -S "$i".sam
done
也很好理解，不过多打了几个字，感觉麻烦了一些，这种好像更加能适应不同的应用环境一些

#好像自己很少用命令行来做文本处理，不过某些命令还是挺方便的，比如awk
awk '{if($7!=".")print}' NONCODEv6_human_hg38_lncRNA.gtf >test
