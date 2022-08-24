# my learing 

Try to do something
Trust the process

sudo chown -R leelee ./share/

STAR --runThreadN 18 --genomeDir ~/biodata/index/GRCh38 --readFilesIn rmrrna/${id}.derRNA.fq.gz --outFilterType BySJout --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --outFileNamePrefix ${id} --readFilesCommand zcat --outReadsUnmapped None --alignEndsType EndToEnd

ensembl_id <- substr(row.names(norm),1,15) #老记不住这个函数

jupyter notebook
