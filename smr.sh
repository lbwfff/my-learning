#!/bin/bash
#SBATCH -p compute
#SBATCH -c 5
#SBATCH --mem 20000
#SBATCH -t 4:00:00

/scratch/lb4489/project/GWAS/smr-1.3.1-linux-x86_64/smr --beqtl-summary BrainMeta_trans_eQTL_chr1 --query 5.0e-8 --out chr1

#从smr数据中提取summary数据
#也可以直接用smr工具进行分析，不大确定和twosampleMR有什么区别？

/scratch/lb4489/project/GWAS/smr-1.3.1-linux-x86_64/smr --bfile /scratch/lb4489/project/GWAS/csmr/data/1kg/EUR --gwas-summary /scratch/lb4489/project/GWAS/csmr/data/finn_adj  --beqtl-summary /scratch/lb4489/project/GWAS/BrainMeta_trans_eqtl_summary/BrainMeta_trans_eQTL_chr1 --out testsmr --thread-num 10 --diff-freq-prop 0.8 --diff-freq 0.8
#分析的话可以这么分析，gwas-summary需要SNP，A1，A2，freq，b，se，p，n这几列，freq为A1freq，--diff-freq-prop 和 --diff-freq 0.8参数用来比较A1freq (including the GWAS summary data, the eQTL summary data and the LD reference data)
#但是我这个数据比较奇怪可能freq错了还是怎么的，所以我把两个参数调到了非常大，--beqtl-summary后面是下载的eqtl文件，也可以用summary文件制备
