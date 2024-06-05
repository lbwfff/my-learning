#!/bin/bash
#SBATCH -p compute
#SBATCH -c 20
#SBATCH -t 2:00:00

#使用rmats对可变剪切进行分析

module load all
module load gencore/2
module load rmats/4.1.2 #服务器就有，还挺方便的

source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate

conda activate rna #服务器的版本和依赖存在兼容问题，自己重新安装了包

rmats.py --gtf /scratch/lb4489/bioindex/gencode.vM33.annotation.gtf --b1 b1.txt --b2 b2.txt --nthread 20 --readLength 150 --od ./ot --tmp ./tmp #可以直接对比对后的bam文件进行分析，b1.txt就是用逗号分隔的输入文件

rmats2sashimiplot --b1 treat_rep1.sorted.bam,treat_rep2.sorted.bam,treat_rep3.sorted.bam --b2 ctrl_rep1.sorted.bam,ctrl_rep2.sorted.bam,ctrl_rep3.sorted.bam --event-type SE -e ./ot/demo2.txt --l1 treat --l2 ctrl -o SE_rmats 
#这个好像没办法对特定区域或基因进行可视化，但可以通过过滤输入的结果文件实现

#就还挺方便的，如果感兴趣某个具体的可变剪切事件的话
