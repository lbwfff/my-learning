# my learing 

rm(list= ls())

FYI, minknow applies default qscore filters of 8, 9 and 10 for fast, hac and sup models respectively.

.libPaths("~/R/library")

Try to do something
Trust the process

In 2021, I started to get in touch with bioinformatics. This project was used by me to record various codes I wrote during my learning process, including the use of bioinformatics software on the Linux command line, including R language codes, including R language-based Python-based machine learning, and all kinds of messy things. The vast majority of annotations are written in Chinese.

This project is equivalent to my notebook. When I want to use the code I have written, I can find it in this project and repeat it easily (but then there are more and more codes. I want to find a certain line of code in my memory. also become more difficult).

ssh lb4489@10.230.14.44 -p 4410

scp -P 4410 lb4489@10.230.14.44:/path/to/file ~/Downloads/   #下载文件

for d in */; do echo -n "$d "; find "$d" -type f | wc -l; done   #统计文件数量

rsync -avP /archive/lb4489/archive_project/hcc/ ./   #cp

dmfget -d /archive/lb4489/archive_project/hcc/  #把文件从磁带拿出来

awk -F '\t' 'NR>1 {print $7}' filereport_read_run_PRJNA1090549.tsv   | tr ';' '\n'   > fastq_list.txt   #从ENA拿链接

tar -I "zstd -T50" -cf CTX_09.tar.zst CTX_09/ #T50为线程数，会比tar.gz快且小
zstd -t archive.tar.zst
tar -I zstd -xf CTX_09.tar.zst


