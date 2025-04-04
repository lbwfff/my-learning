#conda环境的迁移，折腾这玩意是因为，工作目录不能长时间放置文件，home目录空间又不够，只能把环境暂时放到存档，狗屎
conda create  -p /archive/lb4489/conda/ribotish  --clone /home/lb4489/.conda/envs/ribotish

#GEO,记得之前有一段时间喜欢直接从ENA下载fastq.gz文件，感觉能快一点，但是太麻烦了，geofetch

geofetch -i GSE124535 -n GSE124535 -m `pwd` #下载一整个项目的话
geofetch -i 19_1_down.txt  -n brain_ribo -m `pwd` #下载其中一写数据，文件有三列，项目ID，样本ID，名字，后面两列不填的话就下载整个项目

fasterq-dump -e 8  *.sra #sra转fastq, 这个代码虽然可以多核运行，但是过于多的核心可能并不会有加速，可能十个二十个左右会好一些？
#

infer_experiment.py -i ./OVX14_1-Input.sam -r /home/leelee/biodata/annotation/NCBI_mm39.bed -q 20
#之前一直不太明白怎么判断数据的链特异性问题，这里有一个简单的工具，来自RSeQC
#我得到的结果是
#This is PairEnd Data
#Fraction of reads failed to determine: 0.0495
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.0177
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.9327
#可以看到reads大部分偏第二种情况，说明是链特异性的，read1在+链，相对的gene其实是在-链（reverse）。这种就是“fr-firststrand”，也就是参数中的--rf，first read maps to the reverse strand

fastq-dl --accession PRJNA670603 --provider SRA #一个GEO和ENA数据的下载器，试了一下，下载的数据还是比较快的，它好像是会自动对SRA数据进行处理然后把fastq压缩成gz这样的，添加线程参数或许能让本地处理的部分更快

signalp6 -fasta test.fasta -od ./ -f all -org eukarya --mode fast #一个预测信号肽的软件

bedtools sort -i ./MOV10_hg38.bed > mov10_sorted.bed
bedtools merge -i ./mov10_sorted.bed > MOV10_merged.bed #活用一些小工具，每次都想在R语言中解决问题，一是麻烦二来效率也很低，一些编程好的小工具能够提升效率

#bedtools想要merge但是保留正负链该怎么做
cat ./*bed > all.bed
sort -k1,1 -k2,2n all.bed > all_sorted.bed
bedtools merge -i all_sorted.bed -S + -c 6 -o distinct > plus.bed
bedtools merge -i all_sorted.bed -S - -c 6 -o distinct > minus.bed
cat plus.bed minus.bed > merged.bed
sort -k1,1 -k2,2n merged.bed > merged_sorted.bed

jupyter notebook #

sudo chown -R leelee ./share/

STAR --runThreadN 18 --genomeDir ~/biodata/index/GRCh38 --readFilesIn rmrrna/${id}.derRNA.fq.gz --outFilterType BySJout --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --outFileNamePrefix ${id} --readFilesCommand zcat --outReadsUnmapped None --alignEndsType EndToEnd

liftOver GSM1173263_allTC_4su_MOV10_WT.bed hg19ToHg38.over.chain.gz MOV10_WT.bed Unmap.bed

samtools view -h SRR14510077Aligned.sortedByCoord.out.bam | awk 'length($10) ==28 || $1 ~ /^@/' | samtools view -bS - > ../test.bam
#提取特定长度的reads，主要就是awk的使用，||为逻辑或，后面这个 “$1 ~ /^@/” 我就完全没明白是什么意思了

bamCoverage -b SRR8495834_2930.bam -o SRR8495834_coverage.bw --Offset 12 --normalizeUsing RPKM --outFileFormat bedgraph 
#deeptools可以做这个矫正，但我不知道自己代码这么写是不是正确的？
bigwigCompare --bigwig1 --bigwig2 -o  #然后看能不能用bigwigCompare命令来计算翻译强度

./bigWigAverageOverBed 

cut -f1-17 SRR8590813.ribotricer_translating_ORFs.tsv > test
#使用cut分割‘\t’分割符的文件，比如输出1到17列

locate libcom_err.so.3
sudo ln -s /home/leelee/miniconda3/lib/libcom_err.so.3.0 /usr/lib/x86_64-linux-gnu/libcom_err.so.3
####究极有用的两句话，用于r语言的奇怪报错

###################netMHC和mhcflurry###########################
mhcflurry-predict-scan ./hla.fasta --alleles HLA-A0101 HLA-A0201 HLA-A0301 --out tets.csv #一个预测mhc结合的软件
mhcflurry-predict-scan ./up_ncorfs.fasta --alleles HLA-A0101 HLA-A0201 HLA-A0301 HLA-A2402 HLA-A2601 HLA-B0702 HLA-B0801 HLA-B1501 HLA-B2705 HLA-B3901 HLA-B4001 HLA-B5801 --out tets.csv

#结果文件的描述在他们的网站里面能找到例如https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1
./netMHCpan ../test.fasta -BA -xls -a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B15:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01 -xlsfile test_NetMHCpan_out.xls #又想改用netMHC了
./netMHCstabpan -ia ../test.fasta -a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B15:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01 -xls -xlsfile test_ia.xls #好像有一点问题，没法批量的做HLA
./netMHCIIpan -f ../test.fasta -BA -xls -xlsfile MHCIIpan.xls #II类有哪些等位基因要预测呢？

####################################################
#bedGraph转bw
LC_COLLATE=C sort -k1,1 -k2,2n ./genome/SRR10441364Aligned.sortedByCoord.out.bam_coverage_plus.bedgraph > ./bw/SRR104413644.plus.sorted.bedGraph
bedGraphToBigWig ./bw/SRR104413644.plus.sorted.bedGraph ~/biodata/index/GRCh38/GRCh38.info.txt ./bw/SRR104413644.plus.bw

ctrl+c终结进程
rm 删除
复制命令:Ctrl + Shift + C 组合键. 粘贴命令:Ctrl + Shift + V 组合键
sudo su#权限
环境变量：（1）：export PATH=/home/uusama/mysql/bin:$PATH #当前终端有用
（2）：vim /etc/profile 

ascp -l 100M -QT -k 2 -i /home/leelee/miniconda3/envs/p3/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/FASTA/nt.gz ./ #拉跨是真的拉跨，但某些时候可能有奇效？谁知道呢


#!bin/bash
脚本的抬头

zcat CRR042287_f1.fastq.gz |head -n 100|grep 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
#带grep标红你想要的东西

hisat2 --dta -t -x /home/leelee/biodata/index/GRCH38/genome -U SRR9595630_trimmed.fq.gz -S 630.sam -p 18
hisat2 -p 18 --dta -t -x /home/leelee/biodata/index/grch38_snp_rep/genome_snp_rep -U clean/CRR042286_r2_val_2.fq.gz -S 5-3-IP.sam
hisat2 -p 18 --dta -t -x /home/leelee/biodata/index/grch38_snp_rep/genome_snp_rep -1 CRR042286_f1_val_1.100bp_3prime.fq.gz -2 CRR042286_r2_val_2.100bp_3prime.fq.gz -S 5-3-IP.sam
现在没用用到过grch38_snp_rep这个索引了，因为exomepeak2 callpeak时会报错，应该是注释信息和索引没有很好的对应的原因。

bowtie2 -x ~/biodata/index/rRNA/rRNA.homo_sapiens --un-gz SRR7366868_trimmed.derRNA.fq.gz -U SRR7366868_trimmed.fq.gz -p 8 -S SRR7366868_trimmed.fq.gz.rRNA.mapped.sam 2>SRR7366868_trimmed.fq.gz.log
使用bowtie2去除rRNA，一堆莫名其妙的参数，我还没有较好的理解
bowtie2 -x ~/biodata/index/rRNA/rRNA.homo_sapiens --un-conc-gz Input-2-trimmed.derRNA.fq.gz -1 Input-2_L1_702D01.R1_val_1.fq.gz -2 Input-2_L1_702D01.R2_val_2.fq.gz -p 15 -S input-2-rRNA.mapped.sam 2>input-2-trimmed.fq.gz.log
双端的
/home/leelee/tools/bbtools/bbmap/bbduk.sh ref=/home/leelee/biodata/index/GCF/ncbi_dataset/allmycoplasma.fasta threads=18 ordered=t k=31 in=./ADPR/B0h-A-1_L4_704D04.R1_val_1.fq.gz in2=./ADPR/B0h-A-1_L4_704D04.R2_val_2.fq.gz out=B0H-A_bbduk_.R1.fq out2=B0H-A_bbduk_.R2.fq outm=B0H-A_bbduk_.bad.R1.fq outm2=B0H-A_bbduk_.bad.R2.fq
#用bbduk去除污染reads，速度比bowtie2要快得多

cutadapt -j 4 -a GATCGGAAG -a AGATCGGAAG -o 224.fastq SRR11554224.fastq
不怎么用cutadapt

trim_galore -j 10 -length 30 -stringency 3 SRR11554224.fastq
trim_galore --paired --quality 20 --length 30  -o ../clean CRR042286_f1.fastq CRR042286_r2.fastq.gz -j 17
双端
trim_galore --quality 20 --length 35 --paired -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT raw/Input-2_L1_702D01.R1.fastq.gz raw/Input-2_L1_702D01.R2.fastq.gz -j 17 --stringency 1 --three_prime_clip_R1 3 --three_prime_clip_R2 3 --clip_R1 8 --clip_R2 8
因为reads比较长，所以刀法比较迪拜
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o sham2-input.cutadapt.1.fastq -psham2-input.cutadapt.2.fastq ./OVX/sham2-Input_L2_701D02.R1.fastq.gz ./OVX/sham2-Input_L2_701D02.R2.fastq.gz -j 18 -q 20 -m 35:35
#我遇到有一组数据，用trim_galore会报错不知道为什么，但是用cutadapt就能正常跑出来

for i in Input-1 Input-2 Abcam-1 Abcam-2 CST-1 CST-2 SYSY-1 SYSY-2
do
samtools view -bSq 10 -@ 15 -F 4 -F 256 sam/"$i".sam -o "$i".bam
samtools sort -@ 15 -n -o "$i".sorted.bam "$i".bam
samtools fixmate -@ 15 -m "$i".sorted.bam "$i".fixmate.bam
samtools sort -@ 15 -o "$i".positionsort.bam "$i".fixmate.bam
samtools markdup -r -@ 15 "$i".positionsort.bam "$i".markdup.bam
samtools index "$i".markdup.bam
done
#对sam（bam）文件的处理，排序，去pcr重复什么的

for i in M2_IP M2_Input NC_IP NC_Input T1_IP T1_Input
do
#samtools view -bSq 10 -@ 15 -F 4 -F 256 "$i".sam -o "$i".bam
#samtools sort -@ 15 -n -o "$i".sorted.bam "$i".bam
#samtools fixmate "$i".sorted.bam "$i".fixmate.bam
#samtools sort -@ 15 -o "$i".positionsort.bam "$i".fixmate.bam
samtools rmdup -S "$i".positionsort.bam "$i".markdup.bam
samtools index "$i".markdup.bam
done   #不知道是不是因为版本不同了，以前这一部分的代码不能直接用了，这是改完之后的


bamCoverage -b 633.sort.bam -o shNC-IP.bw -bs 20 --normalizeUsing BPM --skipNAs --ignoreDuplicates
bamCompare -b1 630.sort.bam -b2 631.sort.bam -o shTRA2A.bw -bs 20 --skipNAs --ignoreDuplicates
直接跑bam会失败，需先建立索引,建立索引失败，需先排序

macs2 callpeak -t 630.sort.bam -c 631.sort.bam -f BAM -n shTRA2A -g hs

featureCounts -T 5 -a Caenorhabditis_elegans.WBcel235.103.gtf -o 235 /data/home/huangwx/biodata/BAM/235.bam

fastaFromBed -fi /home/leelee/biodata/index/GRCH38/GRCh38.primary_assembly.genome.fa -fo shTRA2A-motif.fa -bed shTRA2A_summits.bed
fastaFromBed -fi /home/leelee/biodata/index/GRCH38/GRCh38.primary_assembly.genome.fa -bed shTRA2A.bed -fo test-peak.fa -split -s
bedtools getfasta -fi /home/leelee/biodata/index/GRCH38/GRCh38.primary_assembly.genome.fa -bed Mod.bed -fo peak.fa -split -s


nohup annotatePeaks.pl shNC_summits.bed /home/leelee/biodata/index/GRCH38/GRCh38.primary_assembly.genome.fa -gtf /home/leelee/biodata/annotation/gencode.v35.annotation.sorted.gtf > shNC.annotation.xls &
##自己的基因组自己的注释文件，试着使用-cpu参数

findMotifsGenome.pl Mod.bed /home/leelee/biodata/index/GRCH38/GRCh38.primary_assembly.genome.fa motif -p 10 -len 5,6,7 -size -100,100 -rna -mis 1 -basic
#这个五基序非常完美，可为什么六基序完全不行？
findMotifsGenome.pl motif.bed /home/leelee/biodata/index/GRCH38/GRCh38.primary_assembly.genome.fa motif -p 10 -len 6 -size -100,100 -rna -mis 0 -basic
#去掉bed中的Inf后再做这个，感觉结果接近完美了，加上-mis 0能得到完美的motif

meme-chip -rna shTRA2A-motif.fa
meme-chip -rna -meme-maxw 12 -meme-p 5 shTRA2A-motif.fa
meme-chip -rna -meme-minw 5 -meme-maxw 7  -meme-p 10 
做meme没得到过理想的效果，所以用得比较少

awk '{if($5>=300) print }' Mod.bed > motif.bed
对bed文件的简单处理

ascp -l 100M -P 33001 -QT -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR959/000/SRR9595630/SRR9595630.fastq.gz ./SRR9595630.fastq.gz
ascp很不稳，所以对于断点重连的参数非常重要

computeMatrix scale-regions -S <biwig file(s)> -R /home/leelee/biodata/annotation/gencode.v35.annotation.sorted.gtf -m 100 -bs 5 --metagene -o test -p 10
computeMatrix scale-regions -R /home/leelee/biodata/annotation/ -S x.bw --smartLabels -p 20 --binSize 10 -b 3000 -a 3000 --regionBodyLength 5000 --sortRegions keep -o x.gz --outFileSortedRegions computeMatrix_x.bed --outFileNameMatrix matrix_x.tab

computeMatrix scale-regions -S H3K27Me3-input.bigWig \
                                 H3K4Me1-Input.bigWig  \
                                 H3K4Me3-Input.bigWig \

computeMatrix scale-regions -R /home/leelee/biodata/annotation/ucschuman.bed -S shTRA2A-IP.bw --smartLabels -p 20 --binSize 10 -b 3000 -a 3000 --regionBodyLength 5000 --sortRegions keep -o test.gz --outFileSortedRegions computeMatrix_test.bed --outFileNameMatrix matrix_test.tab

plotProfile -m test.gz -out ExampleProfile1.png
plotHeatmap -m test.gz -out ExampleHeatmap1.png 

computeMatrix scale-regions -S shTRA2A-IP.bw \
                                 shTRA2A-input.bw  \
                                 shNC-IP.bw \
                                 shNC-input.bw  \
                              -R /home/leelee/biodata/annotation/ucschuman.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o matrix.mat.gz -p 20
deeptools里面的一些工具，没有用过


nohup ~/tools/R-4.1.0/bin/Rscript exomepeak.R &
后台跑r语言的句子

stringtie input.markdup.bam -p 15 -G ~/biodata/annotation/gencode.v35.annotation.sorted.gtf -A input.txt -o input.gft
感觉比count工具要方便，就是这玩意要怎么对多个文件进行处理呢？

trimmomatic PE \
 -threads 15 \
 -phred33 \
 -trimlog trim.log \
 Input-1_L1_701D01.R1.fastq.gz Input-1_L1_701D01.R2.fastq.gz \
 trim-input1-r1.fastq.gz trim-input1-r1-unpair.fastq.gz trim-input-1-r2.fastq.gz trim-input-1-r2-unpair.fastq.gz\
 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true \
 SLIDINGWINDOW:5:20 \
 LEADING:3 \
 TRAILING:3 \
 MINLEN:36
主要对这个fa文件不太懂，感觉挺麻烦的，trim-galore目前已经完全能满足我的需求了，有机会再去整这些东西吧

source /home/leelee/miniconda3/etc/profile.d/conda.sh
conda activate p2
之前不是没法在脚本里改变conda环境吗，用这个方法。

##RSEM

rsem-prepare-reference --gtf /scratch/lb4489/bioindex/gencode.vM33.annotation.gtf \
                        /scratch/lb4489/bioindex/GRCm39.genome.fa \
                        /scratch/lb4489/bioindex/rsem_Cm39 #这个reference到底是一个什么东西，看起来就是不同转录本的序列信息


rsem-calculate-expression --alignments \
                          -p 40 \
                          ./ERR1551315Aligned.toTranscriptome.out.bam \
			                    ./ERR1551316Aligned.toTranscriptome.out.bam \
			                    ./ERR1551317Aligned.toTranscriptome.out.bam \
			                    ./ERR1551318Aligned.toTranscriptome.out.bam \
		                  	  ./ERR1551319Aligned.toTranscriptome.out.bam \
			                    ./ERR1551320Aligned.toTranscriptome.out.bam \
                          /scratch/lb4489/bioindex/rsem_Cm39 \
                          ERR_quals

