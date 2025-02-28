#最近在尝试从RNA-seq中直接寻找修饰信号的方法，起码看起来在m1A中是可以的，记录一下代码#!/bin/bash

#SBATCH -p compute
#SBATCH -c 10
#SBATCH -t 48:00:00
#SBATCH --mem=50G

module load all
module load gencore
module load gencore/2
module load bbmap

for s in SRR9320014 SRR9320016 SRR9320017 SRR9320018 SRR9320019 SRR9320020 SRR9320021 SRR9320022 SRR9320023 SRR9320024 SRR9320025; do

clumpify.sh in=./$s.fastq.gz out=./$s.clumped.fastq.gz dedupe subs=0 -Xmx48g

rm ./$s.fastq.gz

done

#dedupe=f Remove duplicate reads.这里应该是有去除PCR重复导致的reads的，我不确定这个工具是不是针对有UMI的reads，反正在他们的数据中就这么用了

for s in SRR9320014 SRR9320016 SRR9320017 SRR9320018 SRR9320019 SRR9320020 SRR9320021 SRR9320022 SRR9320023 SRR9320024 SRR9320025; do 
cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 -m 30 -j 115 -o ./$s.cut.fastq.gz ./$s.clumped.fastq.gz
wait
done

#然后在这个基础上再进行的修剪

for s in SRR6177764 SRR6177765 SRR6177766 SRR6177767; do 

STAR --runThreadN 120 --outSAMattributes All  --genomeDir /scratch/lb4489/bioindex/GRCh38_1000G \
     --readFilesIn ./"$s"_trimmed.fq.gz \
     --outFilterType BySJout \
	 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
	 --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --outFileNamePrefix $s --outReadsUnmapped None\
 	 --twopassMode Basic  --readFilesCommand zcat \
     --clip5pNbases 11  \

done


for i in SRR9320014 SRR9320017 SRR9320015 SRR9320016 SRR9320018 SRR9320019 SRR6177765 SRR6177767 SRR6177764 SRR6177766; do
    gatk --java-options '-Xmx18G' MarkDuplicates \
        -I "$i"Aligned.sortedByCoord.out.bam \
        -O "$i"_sorted_makd.bam \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        --TMP_DIR /scratch/lb4489/project/mttRNA/TMP/${i} \
        -M "$i".metrics --REMOVE_DUPLICATES true &
done

wait #对于这一组数据我就没有运行picard了，本来深度就一般



for i in  SRR5418438 SRR5418440 SRR5418442 SRR5418445 SRR9320014 SRR9320017  SRR9320016 SRR9320018 SRR9320019 SRR6177765 SRR6177767 SRR6177764 SRR6177766; do #SRR9320015

samtools mpileup --positions /scratch/lb4489/project/mttRNA/m1A-seq/m5cSite_CHEUI.bed \
    -ugf /scratch/lb4489/bioindex/GRCh38_full_analysis_set_plus_decoy_hla.fa  \
    "$i"Aligned.sortedByCoord.out.bam| bcftools call -m -A --ploidy GRCh38 -Oz -o "$i"_list_m5c.vcf.gz &

done

#然后就是找variant，这里我只用来一组位置：m5cSite_CHEUI.bed，这会极大的提升速度，然后我也不用去费劲讨论位点的准确性的问题了

#之后就是使用R语言的脚步来进行合并了，我在bcftools中加入了参数--ploidy GRCh38，这使得我可以在VCF文件中看到多个ALT allele（如果存在的话），这会是一个对于质控非常重要的信息。

###############################################################################
##然后是两段R代码，怎么从fasta切序列和怎么从bam切coverage

library(Biostrings)
library(GenomicRanges)

fasta_file <- "/scratch/lb4489/bioindex/GRCh38_full_analysis_set_plus_decoy_hla.fa"
genome <- FaFile(fasta_file)

query_granges <- GRanges(seqnames = vcf$`#CHROM`,
                          ranges = IRanges(start = c(vcf$POS-5), width = 10))
sequences <- getSeq(genome, query_granges)

print(sequences)

regex_pattern <- "GTTC[ATCG]A[ATCG][ATCG]C" #正则化的表达方式，对于核苷酸序列好像很好用
matches <- grep(regex_pattern, as.character(sequences))
vcf$motif<-as.character(sequences)

moti<-vcf[matches,]

########
library(GenomicAlignments)
  
bam_file <- meta$bam_pos
  
query_range <- GRanges(seqnames = posbed$V1,
                         ranges = IRanges(start = c(posbed$V2-20), 
                                          width = 40))
  
posalignments_list <- lapply(bam_file, function(bam) {
    readGAlignments(bam, param = ScanBamParam(which = query_range))
  }) #用循环可以同时切好几个bam文件

cove<-posalignments_list[[j]] #选择其中一个
cove <- coverage(cove) #计算覆盖率

coverage_at_position <-cove[[chr]][(site-5):(site+5)]
    
merge$score[i]<- ifelse(mean(coverage_at_position[1:10], na.rm = TRUE) > 20,
                            (mean(coverage_at_position[6:10], na.rm = TRUE)/mean(coverage_at_position[1:5], na.rm = TRUE)),
                            NA)
  




