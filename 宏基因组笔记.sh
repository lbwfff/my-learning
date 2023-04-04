#宏基因组，之前没有接触过的邻域

#qiime2

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime  #不想要交互界面的话，需要去除-it

#docker如果想在.sh脚本里面运行的话，需要删掉-t -i
#本来想用conda装的，因为网络问题怎么也装不上，还是用了docker，docker总感觉卡卡的

docker container prune -f #停止所有容器
docker rm $(docker ps -a | awk '/qiime/ {print $1}') #仅停qiime2的容器

#随便跑了一些demo，目前的感受就是，这玩意得到的结果这么散装且稀碎，到时候要怎么整合呢？

#有一个问题，就是我用docker run的时候每次都激活了一个新的容器，有没有什么办法在一个容器里面进行全部操作呢？我觉得这是一个docekr使用的问题。折腾了半天，没有想明白要怎么做。
###########################人体微生物数据########################################################
#导入数据
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime tools import --type EMPSingleEndSequences --input-path ./emp-single-end-sequences --output-path emp-single-end-sequences.qza
#文件夹emp-single-end-sequences里有barcodes.fastq.gz和sequences.fastq.gz两个文件

#拆分数据
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime demux emp-single --i-seqs emp-single-end-sequences.qza --m-barcodes-file sample-metadata.tsv --m-barcodes-column barcode-sequence --o-per-sample-sequences demux.qza --o-error-correction-details demux-details.qza
#sample-metadata.tsv里有每个barcode对应的sample ID，用于拆分

#拆分结果的统计
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime demux summarize --i-data demux.qza --o-visualization demux.qzv

#质量控制
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-trim-left 0 --p-trunc-len 120 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza 
#可视化
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv 

mv rep-seqs-dada2.qza rep-seqs.qza
mv table-dada2.qza table.qza

#另外一种质控方法，据说这种方法会慢挺多的，后续也没有用到就没有跑
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime quality-filter q-score --i-demux demux.qza --o-filtered-sequences demux-filtered.qza --o-filter-stats demux-filter-stats.qza
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime deblur denoise-16S --i-demultiplexed-seqs demux-filtered.qza --p-trim-length 120 --o-representative-sequences rep-seqs-deblur.qza --o-table table-deblur.qza --p-sample-stats --o-stats deblur-stats.qza

docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime metadata tabulate --m-input-file demux-filter-stats.qza --o-visualization demux-filter-stats.qzv
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime deblur visualize-stats --i-deblur-stats deblur-stats.qza --o-visualization deblur-stats.qzv

#特征表
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file sample-metadata.tsv

docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

#进化树
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

#Alpha和beta多样性分析，不太理解，应该是宏基因组独特的指标了
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 1103 --m-metadata-file sample-metadata.tsv --output-dir core-metrics-results

#然后就是一些组间差异啊，可视化什么的，因为不知道两个值是什么意思所以不太能明白
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/faith-pd-group-significance.qzv
  
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/evenness-group-significance.qzv
  
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column body-site --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv --p-pairwise
  
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column subject --o-visualization core-metrics-results/unweighted-unifrac-subject-group-significance.qzv --p-pairwise
  
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes days-since-experiment-start --o-visualization core-metrics-results/unweighted-unifrac-emperor-days-since-experiment-start.qzv

docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime emperor plot --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes days-since-experiment-start --o-visualization core-metrics-results/bray-curtis-emperor-days-since-experiment-start.qzv
#都大差不差的，懒得跑了

#Alpha稀疏曲线
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 4000 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction.qzv
  
#物种组成
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-515-806-nb-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
#一个已经训练好的分类器，然后发现版本不同还没法通用，需要重新训练？
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.qzv

#ANCOM差异丰度分析
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-table filter-samples --i-table table.qza --m-metadata-file sample-metadata.tsv --p-where "[body-site]='gut'" \--o-filtered-table gut-table.qza
  
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition add-pseudocount --i-table gut-table.qza --o-composition-table comp-gut-table.qza
  
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition ancom --i-table comp-gut-table.qza --m-metadata-file sample-metadata.tsv --m-metadata-column subject --o-visualization ancom-subject.qzv

docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime taxa collapse --i-table gut-table.qza --i-taxonomy taxonomy.qza --p-level 6 --o-collapsed-table gut-table-l6.qza

docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition add-pseudocount --i-table gut-table-l6.qza --o-composition-table comp-gut-table-l6.qza

docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition ancom --i-table comp-gut-table-l6.qza --m-metadata-file sample-metadata.tsv --m-metadata-column subject --o-visualization l6-ancom-subject.qzv


###############训练分类器####################
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime tools import --type 'FeatureData[Sequence]' --input-path 85_otus.fasta --output-path 85_otus.qza

docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path 85_otu_taxonomy.txt --output-path ref-taxonomy.qza
  
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-classifier extract-reads --i-sequences 85_otus.qza --p-f-primer GTGCCAGCMGCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT --p-trunc-len 120 --p-min-length 100 --p-max-length 400 --o-reads ref-seqs.qza

docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza --o-classifier classifier.qza
  
docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

docker run  -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv  
  

###############一个测试pipline####################
###############这部分我包装为sh代码了，但依然还不够，可以参考之前中山大学他们组做的那种sh脚本，应该就可以在自动化上面有极大的提升，而且也没有太麻烦####################
###############后面分析的流程其实大同小异，但是值得去多想的就是前期的数据导入，去噪这些流程，尤其是数据的导入，其实还挺麻烦的####################################
###############我用R语言太多了，导致awk的语法变成了很陌生的东西##########################

mkdir -p seq
awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_f1.fq.gz -O seq/"$1"_1.fq.gz")}' <(tail -n+2 metadata.txt)
awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_r2.fq.gz -O seq/"$1"_2.fq.gz")}' <(tail -n+2 metadata.txt)
#根据metadata.txt的内容下载数据，话说每个数据都好小啊，应该只是测试数据集的原因
#然后这个数据里没有需要样本拆分的问题（或者是大部分数据都没有这步？），也没有需要去接头的问题

awk 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} NR>1{print $1"\t/data/seq/"$1"_1.fq.gz\t/data/seq/"$1"_2.fq.gz"}' metadata.txt > manifest
#根据metadata.txt编写manifest（其实就是文件和文件位置），高端的awk写法,因为我使用的docker版的软件，所以manifest中描述的文件位置信息需要做修改

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest --output-path demux.qza --input-format PairedEndFastqManifestPhred33V2
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-n-threads 8 --p-trim-left-f 29 --p-trim-left-r 18 --p-trunc-len-f 0 --p-trunc-len-r 0 --o-table dada2-table.qza --o-representative-sequences dada2-rep-seqs.qza --o-denoising-stats denoising-stats.qza
#这一步是比较费时的，36个7m不到的fq文件，18核跑了一个小时，这个效率是不是有点过于emmmmm

cp dada2-table.qza table.qza
cp dada2-rep-seqs.qza rep-seqs.qza

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file metadata.txt

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 27060 --m-metadata-file metadata.txt --output-dir core-metrics-results

index=observed_features
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/${index}_vector.qza --m-metadata-file metadata.txt --o-visualization core-metrics-results/${index}-group-significance.qzv

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 34663 --m-metadata-file metadata.txt --o-visualization alpha-rarefaction.qzv

distance=weighted_unifrac
column=Group
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza --m-metadata-file metadata.txt --m-metadata-column ${column} --o-visualization core-metrics-results/${distance}-${column}-significance.qzv --p-pairwise

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-classifier classify-sklearn --i-classifier classifier_gg_13_8_99_V5-V7.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.txt --o-visualization taxa-bar-plots.qzv

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition add-pseudocount --i-table table.qza --o-composition-table comp-table.qza

column=Group
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition ancom --i-table comp-table.qza --m-metadata-file metadata.txt --m-metadata-column ${column} --o-visualization ancom-${column}.qzv

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime taxa collapse --i-table table.qza --i-taxonomy taxonomy.qza --p-level 6 --o-collapsed-table table-l6.qza
 
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition add-pseudocount --i-table table-l6.qza --o-composition-table comp-table-l6.qza

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition ancom --i-table comp-table-l6.qza --m-metadata-file metadata.txt --m-metadata-column ${column} --o-visualization l6-ancom-${column}.qzv

#其实想想，不同数据好像就导入的步骤不太一样，其它就大同小异了


###########################################################################################################################################################
#试着写了一个脚本，使用方法就是bash qiime2.sh manifest.txt metadata.txt gg-13-8-99-nb-classifier.qza这样

echo "Usage: sh qiime.sh manifest metadata classifier" #因为这一大段需要的输入其实也就是这三个文件了，还挺容易组合起来的，剩下的就是怎么更高效的搓出来manifest核metadata了。

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path "$1" --output-path demux.qza --input-format PairedEndFastqManifestPhred33V2
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-n-threads 18 --p-trim-left-f 29 --p-trim-left-r 18 --p-trunc-len-f 0 --p-trunc-len-r 0 --o-table dada2-table.qza --o-representative-sequences dada2-rep-seqs.qza --o-denoising-stats denoising-stats.qza

cp dada2-table.qza table.qza
cp dada2-rep-seqs.qza rep-seqs.qza

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file "$2"
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 27060 --m-metadata-file "$2" --output-dir core-metrics-results
#这步报了一堆错，但是结果好像正常的出来了？

index=observed_features
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/${index}_vector.qza --m-metadata-file "$2" --o-visualization core-metrics-results/${index}-group-significance.qzv
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 34663 --m-metadata-file "$2" --o-visualization alpha-rarefaction.qzv

distance=weighted_unifrac
column=Group
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza --m-metadata-file "$2" --m-metadata-column ${column} --o-visualization core-metrics-results/${distance}-${column}-significance.qzv --p-pairwise
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-classifier classify-sklearn --i-classifier "$3" --i-reads rep-seqs.qza --o-classification taxonomy.qza
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file "$2" --o-visualization taxa-bar-plots.qzv
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition add-pseudocount --i-table table.qza --o-composition-table comp-table.qza

column=Group
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition ancom --i-table comp-table.qza --m-metadata-file "$2" --m-metadata-column ${column} --o-visualization ancom-${column}.qzv
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime taxa collapse --i-table table.qza --i-taxonomy taxonomy.qza --p-level 6 --o-collapsed-table table-l6.qza
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition add-pseudocount --i-table table-l6.qza --o-composition-table comp-table-l6.qza
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition ancom --i-table comp-table-l6.qza --m-metadata-file "$2" --m-metadata-column ${column} --o-visualization l6-ancom-${column}.qzv

docker rm $(docker ps -a | awk '/qiime/ {print $1}') #把产生的一堆乱七八糟的容器全部删掉

###############################################################################################################################################
#宏基因组不止扩增子，使用基于dada2的qiime2不知道适不适用其它微生物组的数据，有一个Kneaddata加MetaPhlAn的流程可以使用

kneaddata --input ${i} --input ${i%R1_fastpout.fq.gz}R2_fastpout.fq.gz --reference-db /home/leelee/biodata/index/metaphlan --output knea --threads 16

kneaddata --input con-Input-1_L4_701D04.R1.fastq.gz --input con-Input-1_L4_701D04.R2.fastq.gz --reference-db /home/leelee/biodata/index/GRCm39 --output knea --threads 10

kneaddata -i ./seq/14818_1.fq.gz -i ./seq/14818_2.fq.gz --reference-db /home/leelee/biodata/index/Amel --trimmomatic /home/leelee/tools/trimmomatic  --output knea --threads 10 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output

metaphlan --input_type fastq --bowtie2db /home/leelee/biodata/index/metaphlan -x mpa_v31 /home/leelee/share/qiime2/test2/knea/14817_1_kneaddata_paired_1.fastq -o /home/leelee/share/qiime2/test2/metaphlan/14817.metaphlan.txt

merge_metaphlan_tables.py *metaphlan.txt > merged_abundance_table.txt


