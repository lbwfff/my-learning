############
#QAPA

qapa fasta -f GRCm38.p6.genome.fa qapa_3utrs.gencode_VM22.mm10.bed output_sequences.fa #根据3UTR注释提取对应的序列

module load all
module load gencore
module load gencore/2
module load salmon

salmon index -t output_sequences.fa -i utr_library #QAPA使用了salmon进行定量分析，好处是非常快速，缺点是后续可视化麻烦

#for i in SRR1971690 SRR1971691 SRR1971692 SRR1971693 SRR1971694 SRR1971695 SRR1971696 SRR1971697 SRR1971698 SRR1971699 ; do #SRR1971688 SRR1971689
#salmon quant \
#  -i utr_library \
#  -l A \
#  -1 /scratch/lb4489/project/rnaloca/"$i"_1_val_1.fq \
#  -2 /scratch/lb4489/project/rnaloca/"$i"_2_val_2.fq \  #双端
#  -o "$i"_quant \
#  -p 10
#done

#for i in ERR2817191 ERR2817192 ERR2817193 ERR2817194 ERR2817195; do
#salmon quant \
#  -i utr_library \
#  -l A \
#  -r /scratch/lb4489/project/rnaloca/"$i"_trimmed.fq.gz \ #单端
#  -o "$i"_quant \
#  -p 10
#done

source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate qapa

qapa quant --db ensembl.identifiers.txt /scratch/lb4489/project/rnaloca/qapa_input/SRR*/quant.sf > pau_results.txt
qapa quant --db ensembl.identifiers.txt /scratch/lb4489/project/rnaloca/MTAB_7251/ERR*/quant.sf > MTAB_7251_results.txt #定量


#########
#DaPars2

module load all
module load gencore
module load gencore/2
module load python/3.11.3
module load bedtools

python DaPars2/src/DaPars_Extract_Anno.py -b mm39_Refseq_anno_UCSC.txt -s mm39_Refseq_id_from_UCSC.txt -o mm39_3UTR_annotation.bed #同样需要准备注释文件

genomeCoverageBed -bg -ibam SRR1971699Aligned.sortedByCoord.out.bam -split > SRR1971699Aligned.sortedByCoord.out.wig
genomeCoverageBed -bg -ibam SRR1971698Aligned.sortedByCoord.out.bam -split > SRR1971698Aligned.sortedByCoord.out.wig #把上游得到的bam文件转换为wig格式

python DaPars2/src/DaPars2_Multi_Sample_Multi_Chr.py Dapars2_configure.txt chrList.txt #然后进行分析

#我不太喜欢DaPars2，他们可能主要为APA-QTL分析进行优化，DaPars2的问题首先在于，我不知道确定的远端APA的结构。然后就是这个Dapars2_configure.txt文件要准备还挺麻烦的。
