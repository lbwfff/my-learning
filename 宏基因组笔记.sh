#宏基因组，之前没有接触过的邻域

#qiime2

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime  

#docker如果想在.sh脚本里面运行的话，需要删掉-t -i
#本来想用conda装的，因为网络问题怎么也装不上，还是用了docker，docker总感觉卡卡的

docker stop $(docker ps -aq) #停止所有容器

#随便跑了一些demo，目前的感受就是，这玩意得到的结果这么散装且稀碎，到时候要怎么整合呢？

#有一个问题，就是我用docker run的时候每次都激活了一个新的容器，有没有什么办法在一个容器里面进行全部操作呢？我觉得这是一个docekr使用的问题。
###########################人体微生物数据########################################################
#导入数据
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime tools import --type EMPSingleEndSequences --input-path ./emp-single-end-sequences --output-path emp-single-end-sequences.qza
#文件夹emp-single-end-sequences里有barcodes.fastq.gz和sequences.fastq.gz两个文件

#拆分数据
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime demux emp-single --i-seqs emp-single-end-sequences.qza --m-barcodes-file sample-metadata.tsv --m-barcodes-column barcode-sequence --o-per-sample-sequences demux.qza --o-error-correction-details demux-details.qza
#sample-metadata.tsv里有每个barcode对应的sample ID，用于拆分

#拆分结果的统计
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime demux summarize --i-data demux.qza --o-visualization demux.qzv

#质量控制
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-trim-left 0 --p-trunc-len 120 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza 
#可视化
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv 

mv rep-seqs-dada2.qza rep-seqs.qza
mv table-dada2.qza table.qza

#另外一种质控方法，据说这种方法会慢挺多的，后续也没有用到就没有跑
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime quality-filter q-score --i-demux demux.qza --o-filtered-sequences demux-filtered.qza --o-filter-stats demux-filter-stats.qza
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime deblur denoise-16S --i-demultiplexed-seqs demux-filtered.qza --p-trim-length 120 --o-representative-sequences rep-seqs-deblur.qza --o-table table-deblur.qza --p-sample-stats --o-stats deblur-stats.qza

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime metadata tabulate --m-input-file demux-filter-stats.qza --o-visualization demux-filter-stats.qzv
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime deblur visualize-stats --i-deblur-stats deblur-stats.qza --o-visualization deblur-stats.qzv

#特征表
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file sample-metadata.tsv

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

#进化树
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

#Alpha和beta多样性分析，不太理解，应该是宏基因组独特的指标了
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 1103 --m-metadata-file sample-metadata.tsv --output-dir core-metrics-results

#然后就是一些组间差异啊，可视化什么的，因为不知道两个值是什么意思所以不太能明白
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/faith-pd-group-significance.qzv
  
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/evenness-group-significance.qzv
  
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column body-site --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv --p-pairwise
  
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column subject --o-visualization core-metrics-results/unweighted-unifrac-subject-group-significance.qzv --p-pairwise
  
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes days-since-experiment-start --o-visualization core-metrics-results/unweighted-unifrac-emperor-days-since-experiment-start.qzv

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime emperor plot --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes days-since-experiment-start --o-visualization core-metrics-results/bray-curtis-emperor-days-since-experiment-start.qzv
#都大差不差的，懒得跑了

#Alpha稀疏曲线
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 4000 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction.qzv
  
#物种组成
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-515-806-nb-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
#一个已经训练好的分类器，然后发现版本不同还没法通用，需要重新训练？
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.qzv

#ANCOM差异丰度分析
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-table filter-samples --i-table table.qza --m-metadata-file sample-metadata.tsv --p-where "[body-site]='gut'" \--o-filtered-table gut-table.qza
  
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition add-pseudocount --i-table gut-table.qza --o-composition-table comp-gut-table.qza
  
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition ancom --i-table comp-gut-table.qza --m-metadata-file sample-metadata.tsv --m-metadata-column subject --o-visualization ancom-subject.qzv

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime taxa collapse --i-table gut-table.qza --i-taxonomy taxonomy.qza --p-level 6 --o-collapsed-table gut-table-l6.qza

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition add-pseudocount --i-table gut-table-l6.qza --o-composition-table comp-gut-table-l6.qza

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime composition ancom --i-table comp-gut-table-l6.qza --m-metadata-file sample-metadata.tsv --m-metadata-column subject --o-visualization l6-ancom-subject.qzv


###############训练分类器####################
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime tools import --type 'FeatureData[Sequence]' --input-path 85_otus.fasta --output-path 85_otus.qza

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path 85_otu_taxonomy.txt --output-path ref-taxonomy.qza
  
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-classifier extract-reads --i-sequences 85_otus.qza --p-f-primer GTGCCAGCMGCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT --p-trunc-len 120 --p-min-length 100 --p-max-length 400 --o-reads ref-seqs.qza

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza --o-classifier classifier.qza
  
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv  
  

###############一个测试pipline####################
mkdir -p seq
awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_f1.fq.gz -O seq/"$1"_1.fq.gz")}' <(tail -n+2 metadata.txt)
awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_r2.fq.gz -O seq/"$1"_2.fq.gz")}' <(tail -n+2 metadata.txt)
#根据metadata.txt的内容下载数据，话说每个数据都好小啊，应该只是测试数据集的原因

awk 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} NR>1{print $1"\t/data/seq/"$1"_1.fq.gz\t/data/seq/"$1"_2.fq.gz"}' metadata.txt > manifest
#根据metadata.txt编写manifest，高端的awk写法,因为我使用的docker版的软件，所以manifest还需要改一下

docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest --output-path demux.qza --input-format PairedEndFastqManifestPhred33V2
docker run -v $(pwd):/data quay.io/qiime2/core:2022.8 qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-n-threads 8 --p-trim-left-f 29 --p-trim-left-r 18 --p-trunc-len-f 0 --p-trunc-len-r 0 --o-table dada2-table.qza --o-representative-sequences dada2-rep-seqs.qza --o-denoising-stats denoising-stats.qza
#这一步是比较费时的，可以多线程

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









