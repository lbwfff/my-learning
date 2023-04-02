#######microeco是另一个微生物组的R包，这个包的能力相较于Mpae只强不弱，两个包相互补充，应该就可以解决大部分的问题了###############

library(microeco)
library(file2meco)

test1 <- qiime2meco('table.qza', sample_table = 'metadata.txt', taxonomy_table = 'taxonomy.qza',
                    phylo_tree = 'rooted-tree.qza', rep_fasta = 'rep-seqs.qza', auto_tidy = TRUE)
test1 #导入qiime2结果

class(test1)

#对丰度最高的十个门绘制了柱状图
t1 <- trans_abund$new(dataset = test1, taxrank = "Phylum", ntaxa = 10)
t1$plot_bar(others_color = "grey70", facet = "Group", xtext_keep = FALSE, legend_text_italic = FALSE)
#修改taxrank，genus为属，可以对不同水平的丰度进行可视化，可以自己调是不是要对样本聚类
t1 <- trans_abund$new(dataset = test1, taxrank = "Genus", ntaxa = 8)
t1$plot_bar(bar_type = "notfull", use_alluvium = TRUE, clustering = TRUE, xtext_type_hor = FALSE, xtext_size = 6, color_values = RColorBrewer::brewer.pal(8, "Set2"))

t1 <- trans_abund$new(dataset = test1, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))
#可以看到返回的是一个ggplot对象，想要自己定制化也是很方便的

t1 <- trans_abund$new(dataset = test1, taxrank = "Class", ntaxa = 15)
t1$plot_box(group = "Group") 
#每个物种差异的箱线图,这个效果还挺好看的

t1 <- trans_abund$new(dataset = test1, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
#热图

t1 <- trans_abund$new(dataset = test1, taxrank = "Genus", ntaxa = 5)
t1$plot_line()
#折线图

#
dataset1 <- test1$merge_samples(use_group = "Group")

t1 <- trans_venn$new(dataset1, ratio = NULL)
t1$plot_venn() #韦恩，没太明白这个数是什么意思

t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn() #显示了比例，但感觉比例显然是错的，emmmmm


#######
dataset1 <- test1$merge_samples(use_group = "Group")
t1 <- trans_venn$new(dataset1)

t2 <- t1$trans_comm(use_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample
class(t2)

t2$cal_abund()
# transform and plot
t3 <- trans_abund$new(dataset = t2, taxrank = "Genus", ntaxa = 8)
t3$plot_bar(bar_type = "part", legend_text_italic = T, ylab_title = "Frequency (%)", xtext_type_hor = FALSE, color_values = RColorBrewer::brewer.pal(8, "Set2"))
#这个group1&group2是什么意思？这里的&是什么的交集，我有点没明白

t3 <- trans_abund$new(dataset = t2, taxrank = "Phylum", ntaxa = 8)
t3$plot_pie(facet_nrow = 3, color_values = rev(c(RColorBrewer::brewer.pal(8, "Dark2"), "grey50")))

##################################################
t1 <- trans_alpha$new(dataset = test1, group = "Group")
# return t1$data_stat
head(t1$data_stat)

t1$cal_diff(method = "KW")
# return t1$res_diff
head(t1$res_diff)

t1 <- trans_alpha$new(dataset = test1, group = "Group")
t1$cal_diff(method = "KW", formula = "Group")
head(t1$res_diff)

t1$plot_alpha(measure = "Chao1")
t1$plot_alpha(measure = "Chao1", order_x_mean = TRUE, add_sig_text_size = 6)

t1$cal_diff(method = "wilcox")
t1$plot_alpha(measure = "Chao1", shape = "Group")


###感觉，东西非常多，可能功能比MicrobiotaProcess还要齐全
#一时半会有点信息过载，后续慢慢思考学习吧

test1$cal_betadiv(method = "bray") #gitbook里少了这一行代码
t1 <- trans_beta$new(dataset = test1, group = "Group", measure = "bray")

# PCoA, PCA and NMDS are available
t1$cal_ordination(ordination = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))

{
  t1$plot_ordination(plot_color = "Group", point_size = 5, point_alpha = .2, plot_type = c("point", "ellipse"), ellipse_chull_fill = FALSE)
  t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "centroid"))
  t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse", "centroid"))
  t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull"))
  t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull", "centroid"))
  t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("chull", "centroid"))
  t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = "centroid")
  t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = "centroid", centroid_segment_alpha = 0.9, centroid_segment_size = 1, centroid_segment_linetype = 1)
  t1$plot_ordination(plot_color = "Group") + theme(panel.grid = element_blank()) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2)
} #这个PCA有很多种画法，到时候可以重新回过头来看这部分的gitbook

# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")
t1$plot_group_distance(boxplot_add = "mean")
#组距离

# use replace_name to set the label name, group parameter used to set the color
t1$plot_clustering(group = "Group", replace_name = c("Group"))

# manova for all groups when manova_all = TRUE
t1$cal_manova(manova_all = TRUE) #PerMANOVA,不知道是什么东西
t1$res_manova
#这地方，上面一整段其实都是可以做多因素的分析的，不仅是一个分组水平
#多水平分组的代码我没有放进去，到时候可以返回来再理解
