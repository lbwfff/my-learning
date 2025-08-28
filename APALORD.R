library(APALORD)


out_dir <- "./test_apa"

gene_reference <- load_gtf('../biodata/gencode.vM37.primary_assembly.basic.annotation.gtf', cores = 5)

sample1 <- c( "./isoquant/F_SH5_HP/", "./isoquant/M_SH7_HP/")
sample2 <- c( "./isoquant/F_EE4_HP/", "./isoquant/M_EE7_HP/")
#这个输入就是isoquant的结果

reads <- load_samples(sample1, sample2, group1 = "SH", group2 = "EE")

PAS_data <- APALORD::PAS_calling(gene_reference, reads[reads$gene_id %in% gene_reference$gene_id,],
                        cores = 5 , direct_RNA = TRUE)

pas_file <- file.path(out_dir, "PAS_EE_SH.bed")
write.table(PAS_data, file = pas_file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

# PAU quantification
# - Quantify polyadenylation usage per sample

message("Quantifying PAU by sample...")
PAU_data <- PAU_by_sample(gene_reference, reads[reads$gene_id %in% gene_reference$gene_id,],
                          cores = 7, direct_RNA = TRUE)
pau_file <- file.path(out_dir, "PAU_by_sample_EE_SH.tsv")
write.table(PAU_data, file = pau_file, quote = FALSE,
            col.names = TRUE, row.names = FALSE, sep = "\t")

PAU_test_data <- PAU_test(PAU_data, reads[reads$gene_id %in% gene_reference$gene_id,],
                          P_cutoff = 0.05)

pau_test_file <- file.path(out_dir, "EE_SH_PAU_test_data_all.tsv")
write.table(PAU_test_data, file = pau_test_file, quote = FALSE, col.names = TRUE,
            row.names = FALSE, sep = "\t")

dPAU_test_data <- APALORD::end_PAS_examine(PAU_data, reads[reads$gene_id %in% gene_reference$gene_id,], P_cutoff = 0.2, control = "SH",
                                           experimental = "EE", position = "distal", type = "FC")
pPAU_test_data <- APALORD::end_PAS_examine(PAU_data, reads[reads$gene_id %in% gene_reference$gene_id,], P_cutoff = 0.2, control = "SH",
                                           experimental = "EE", position = "proximal", type = "delta")
#pPAU有一个报错，ttest时有两组一样的数据？

APA_data <- APA_profile(gene_reference, reads[reads$gene_id %in% gene_reference$gene_id,],
                        control = "SH", experimental = "EE",
                        cores = 5, direct_RNA = TRUE)
APA_gene_table <- APA_plot(APA_data)


apa_file <- c( "./test_apa/APA_data_EE_SH.tsv")
apa_table_file <- c( "./test_apa/APA_gene_table_EE_SH.tsv")
write.table(APA_data, file = apa_file, quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
write.table(APA_gene_table, file = apa_table_file, quote = FALSE, col.names = TRUE,
            row.names = FALSE, sep = "\t")
#全部运行太慢了，其实可以做一些过滤后再跑核心的基因集

pdf("test_apa/Camk2a_APA_plot.pdf", width = 5, height = 10)
gene_explore(gene_reference, reads, c("Camk2a"),
             APA_table = APA_data, direct_RNA = TRUE)
dev.off()

#居然可以看到Camk2a，有趣

pdf("test_apa/Slc1a2_APA_plot.pdf", width = 5, height = 10)
gene_explore(gene_reference, reads, c("Slc1a2"),
             APA_table = APA_data, direct_RNA = TRUE)
dev.off()

pdf("test_apa/Slc17a7_APA_plot.pdf", width = 5, height = 10)
gene_explore(gene_reference, reads, c("Slc17a7"),
             APA_table = APA_data, direct_RNA = TRUE)
dev.off()

