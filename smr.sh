#!/bin/bash
#SBATCH -p compute
#SBATCH -c 5
#SBATCH --mem 20000
#SBATCH -t 4:00:00

/scratch/lb4489/project/GWAS/smr-1.3.1-linux-x86_64/smr --beqtl-summary BrainMeta_trans_eQTL_chr1 --query 5.0e-8 --out chr1

#从smr数据中提取summary数据
#也可以直接用smr工具进行分析，不大确定和twosampleMR有什么区别？
