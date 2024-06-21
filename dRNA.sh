#dRNA-seq 表观转录领域现在一个很有前景的方向


#!/bin/bash
#SBATCH -p compute
#SBATCH -c 20
#SBATCH -t 2:00:00

source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate drna

module load minimap2/2.24
module load samtools/1.10

#前面应该还有一个fast5转fastq的步骤的，使用guppy_basecaller

minimap2 -ax map-ont -uf -t 20 --secondary=no /scratch/lb4489/project/dRNA/demo/demo.fa /scratch/lb4489/project/dRNA/demo/data/HEK293T-WT-rep1/fastq/basecalled.fastq > ./test/test.sam 2>> ./test/test.bam.log

samtools view -Sb ./test/test.sam | samtools sort -o ./test/test.bam - &>> ./test/test.bam.log

samtools index ./test/test.bam &>> ./test/test.bam.log

nanopolish index -d /scratch/lb4489/project/dRNA/demo/data/HEK293T-WT-rep1/fast5  /scratch/lb4489/project/dRNA/demo/data/HEK293T-WT-rep1/fastq/basecalled.fastq

nanopolish eventalign --reads /scratch/lb4489/project/dRNA/demo/data/HEK293T-WT-rep1/fastq/basecalled.fastq \
--bam ./test/test.bam \
--genome /scratch/lb4489/project/dRNA/demo/demo.fa \
--signal-index \
--scale-events \
--summary ./test/summary.txt \
--threads 20 > ./test/eventalign.txt

#重点是nanopolish eventalign这一步，输出结果可以做xpore和m6anet

xpore dataprep \
--eventalign ./data/HEK293T-METTL3-KO-rep1/nanopolish/eventalign.txt \
--gtf_or_gff ./demo.gtf \
--transcript_fasta ./demo.fa \
--out_dir KO_xpore \
--genome \
--n_processes 20

xpore dataprep \
--eventalign ./data/HEK293T-WT-rep1/nanopolish/eventalign.txt \
--gtf_or_gff ./demo.gtf \
--transcript_fasta ./demo.fa \
--out_dir WT_xpore \
--genome \
--n_processes 20 

xpore diffmod --config Hek293T_config.yml --n_processes 20

m6anet dataprep --eventalign ./data/HEK293T-METTL3-KO-rep1/nanopolish/eventalign.txt \
	--out_dir KO_m6anet \
	--n_processes 20

m6anet dataprep --eventalign ./data/HEK293T-WT-rep1/nanopolish/eventalign.txt \
        --out_dir WT_m6anet \
        --n_processes 20

m6anet inference --input_dir KO_m6anet --out_dir KO_m6anet_result  --pretrained_model HEK293T_RNA004 --n_processes 20 --num_iterations 1000
m6anet inference --input_dir WT_m6anet --out_dir WT_m6anet_result  --pretrained_model HEK293T_RNA004 --n_processes 20 --num_iterations 1000


############################
#basecall的话也挺方便的，使用gpu加速的guppy_basecaller，但是我可以不解压fast5文件吗？文件数量超过系统限制了

#!/bin/bash
#SBATCH -p nvidia
#SBATCH -c 5
#SBATCH --gres=gpu:1
#SBATCH -t 00:40:00

#/scratch/lb4489/project/dRNA/ont-guppy/bin/guppy_basecaller -v
#/scratch/lb4489/project/dRNA/ont-guppy/bin/guppy_basecaller --help

#nvidia-smi --query-gpu=name --format=csv,noheader

/scratch/lb4489/project/dRNA/ont-guppy/bin/guppy_basecaller -i /scratch/lb4489/project/dRNA/20181121_1017_21112018_Sho_METTI3_10/fast5 -s /scratch/lb4489/project/dRNA/20181121_1017_21112018_Sho_METTI3_10/fastq --flowcell FLO-MIN106 --kit SQK-RNA002 --device auto -q 0 -r


############################
#转录组定量的话可以直接

tringtie -L -o ./test/test.gtf ./test/test.bam -p 5

#就还挺方便的，但是这个FPKM TPM是怎么估计出来的？
