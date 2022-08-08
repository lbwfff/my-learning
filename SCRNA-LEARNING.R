suppressMessages(library(scater))
suppressMessages(library(Seurat))
suppressMessages(library(monocle))
suppressMessages(library(scRNAseq)) 
suppressMessages(library(SC3)) 
suppressMessages(library(M3Drop))

rm(list = ls())

fluidigm<-ReprocessedFluidigmData(assays="rsem_counts",ensembl =FALSE, location = TRUE)
ct <- floor(assays(fluidigm)$rsem_counts)#拿到表达矩阵
ct[1:4,1:4] 
sample_ann <- as.data.frame(colData(fluidigm))#拿到临床信息

library(ggpubr)
library(cowplot)
box <- lapply(colnames(sample_ann[,1:19]),function(i) {
  dat <-  sample_ann[,i,drop=F] 
  dat$sample=rownames(dat)
  dat$group='all cells'
  p <- ggboxplot(dat, x = "group", y = i, 
                 add = "jitter" )
  p
})
plot_grid(plotlist=box, ncol=5 )

ct[1:4,1:4] 
counts <- ct
fivenum(apply(counts,1,function(x) sum(x>0) ))
boxplot(apply(counts,1,function(x) sum(x>0) ))
fivenum(apply(counts,2,function(x) sum(x>0) ))
hist(apply(counts,2,function(x) sum(x>0) ))
choosed_genes=apply(counts,1,function(x) sum(x>0) )>0
table(choosed_genes)
counts <- counts[choosed_genes,]

pa <- colnames(sample_ann[,c(1:9,11:16,18,19)])
tf <- lapply(pa,function(i) {
  # i=pa[1]
  dat <-  sample_ann[,i]  
  dat <- abs(log10(dat))
  fivenum(dat)
  (up <- mean(dat)+2*sd(dat))
  (down <- mean(dat)- 2*sd(dat) ) 
  valid <- ifelse(dat > down & dat < up, 1,0 ) 
})

tf <- do.call(cbind,tf)
choosed_cells <- apply(tf,1,function(x) all(x==1))
table(sample_ann$Biological_Condition)
sample_ann=sample_ann[choosed_cells,]
table(sample_ann$Biological_Condition)
ct <- ct[,choosed_cells]

ct[1:4,1:4] 
counts <- ct
fivenum(apply(counts,1,function(x) sum(x>0) ))
boxplot(apply(counts,1,function(x) sum(x>0) ))
fivenum(apply(counts,2,function(x) sum(x>0) ))
hist(apply(counts,2,function(x) sum(x>0) ))
choosed_genes=apply(counts,1,function(x) sum(x>0) )>0
table(choosed_genes)
counts <- counts[choosed_genes,]

#M3Drop
library(M3Drop) 
Normalized_data <- M3DropCleanData(counts, #表达矩阵
                                   labels = sample_ann$Biological_Condition , #分组信息
                                   is.counts=TRUE, min_detected_genes=2000)
dim(Normalized_data$data)
length(Normalized_data$labels)
class(Normalized_data)
str(Normalized_data)

fits <- M3DropDropoutModels(Normalized_data$data)
DE_genes <- M3DropFeatureSelection(Normalized_data$data, 
                                         mt_method="fdr", mt_threshold=0.01)
par(mar=c(1,1,1,1)) 
heat_out <- M3DropExpressionHeatmap(DE_genes$Gene, Normalized_data$data, 
                                    cell_labels = Normalized_data$labels)


#scater
pheno_data <- as.data.frame(colData(fluidigm))
ct <- floor(assays(fluidigm)$rsem_counts)
sce <- SingleCellExperiment(
  assays = list(counts = ct), 
  colData = pheno_data
)
exprs(sce) <- log2(
  calculateCPM(sce ) + 1)
genes <- rownames(rowData(sce))
genes[grepl('^MT-',genes)]#线粒体基因
genes[grepl('^ERCC-',genes)]#ERCC基因
sce <- addPerCellQC(sce, 
                          subset= list(ERCC = grep('^ERCC',rownames(sce))))

plotColData(sce, x = "sum", y="detected", colour_by="") 
tmp <- as.data.frame(rowData(sce))
colnames(tmp)
head(tmp)

keep_feature <- rowSums(exprs(sce) > 0) > 0
table(keep_feature)
sce <- sce[keep_feature,]
sce

tmp <- as.data.frame(colData(sce))#细胞层面的过滤
colnames(tmp) 

sce <- runPCA(sce)
reducedDimNames(sce)
plotReducedDim(sce, "PCA", 
               colour_by = "Biological_Condition" )

library(SC3)
sce <- sc3_estimate_k(sce)

#seurate
names(metadata(fluidigm))
meta <- as.data.frame(colData(fluidigm))
counts <- ct
identical(rownames(meta),colnames(counts))
Pollen <- CreateSeuratObject(counts,
                             meta.data =meta,
                             min.cells = 3, 
                             min.genes = 200, 
                             project = "Pollen")
Pollen
sce <- Pollen

mito.genes <- grep(pattern = "^MT-", x = rownames(x = sce@assays$RNA@data), value = TRUE)
percent.mito <- Matrix::colSums(sce@assays$RNA@data[mito.genes, ]) / Matrix::colSums(sce@assays$RNA@data)
sce <- AddMetaData(object = sce, metadata = percent.mito,
                   col.name = "percent.mito")

#monocle
ct <- floor(assays(fluidigm)$rsem_counts)
ct[1:4,1:4] 
sample_ann <- as.data.frame(colData(fluidigm))

gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  row.names = row.names(ct)
)
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
sc_cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds
library(dplyr)
colnames(phenoData(sc_cds)@data)
sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- estimateDispersions(sc_cds)
cds=sc_cds
cds
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 5))
length(expressed_genes)
cds <- cds[expressed_genes,]
cds
print(head(pData(cds)))

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'

cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 4) 
plot_cell_clusters(cds, 1, 2, color = "Biological_Condition")
table(pData(cds)$Biological_Condition)

Sys.time()
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Biological_Condition")
Sys.time()
sig_genes <- subset(diff_test_res, qval < 0.1)
head(sig_genes[,c("gene_short_name", "pval", "qval")] )

cg=as.character(head(sig_genes$gene_short_name))
plot_genes_jitter(cds[cg,], 
                  grouping = "Biological_Condition", ncol= 2)
plot_genes_jitter(cds[cg,],
                  grouping = "Biological_Condition",
                  color_by = "Biological_Condition",
                  nrow= 3,
                  ncol = NULL )

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "Biological_Condition")  

