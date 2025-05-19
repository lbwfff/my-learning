#算是补充了基因组学项目中一直想做但是没有做的一步了

source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate polyfun

mkdir -p output

#首先使用polyfun得到预计算的先验概率（precomputed prior causal probabilities），我的理解是这是根据在其他队列中计算的概率（15个UK Biobank traits）作为先验的知识推测本项研究所用的队列的概率

python ./polyfun/extract_snpvar.py --sumstats ./PGC3_for_polyfun.txt --out output/snps_with_var.gz --allow-missing

cat output/snps_with_var.gz | zcat | head

#之后就是finemapping了，我们使用susie进行计算
#这里我们使用biobank的数据作为LD 参考，数据量更大也更加匹配我们的数据

python ./polyfun/finemapper.py --ld ./LD_blocks/chr1_27000001_30000001 \
  --sumstats	./output/snps_with_var.gz --n 1222882 \
  --chr 1 --start 28563133 --end 29596287 \
  --method susie --max-num-causal 2 \
  --out output/YTHDF2UKB.gz

#但是一个比较麻烦的点就是数据量很大必须一个窗口一个窗口的进行分析，有些基因正好跨过两个窗口就很麻烦
#另外PGC PTSD3文章中对于计算窗口的定义好像是不一样的，我们这个手动的计算了上下各500kb的窗口，原文是根据FUMA来定义的
#不太确定要怎么理解这个想法，FUMA决定的位置是SNP2GENE吗？它是考虑了基因组距离吗？hi啊好死会有一些更加精细的数据集来进行定义？我们不知道，但是我们可以尝试去使用这个服务，看起啦很有意思
#
