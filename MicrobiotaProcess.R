#一天到晚的，闲着也是闲着，不如提前学点东西
#MAPE是余教授组做的一个宏基因组，微生物组的包，而且对于qiime2可以有一个很好的对接，想试着学一下
library(MicrobiotaProcess)

mpse2 <- mp_import_qiime2(otuqza='table.qza', taxaqza='taxonomy.qza', mapfilename='metadata.txt')  #导入qiime2的结果
mpse2

library(ggplot2)

mpse2 %<>% mp_rrarefy() 

mpse2 %<>% 
  mp_cal_rarecurve(.abundance = RareAbundance,chunks = 400)

mpse2 %>% print(width=180)

p1 <- mpse2 %>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = Observe,
  )

p2 <- mpse2 %>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = Observe, 
    .group = Group #组别
  ) +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF",'#89ab7c')) +
  scale_fill_manual(values=c("#00A087FF", "#3C5488FF",'#89ab7c'), guide="none")

# combine the samples belong to the same groups if 
# plot.group=TRUE
p3 <- mpse2 %>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = "Observe", 
    .group = Group, 
    plot.group = TRUE
  ) +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF",'#89ab7c')) +
  scale_fill_manual(values=c("#00A087FF", "#3C5488FF",'#89ab7c'),guide="none")  

p1 + p2 + p3 #p1为各个样本的数据，p2、p3则做了分组，但是这个图要怎么理解？横坐标是reads num，但是纵坐标是什么？


############################
mpse2 %<>% 
  mp_cal_alpha(.abundance=RareAbundance)
mpse2 #为什么Pielou值不被计算呢？

f1 <- mpse2 %>% 
  mp_plot_alpha(
    .group=Group, 
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson) #和gitbook不同的是，Pielou值并没有被计算，可能是包语言版本的问题？我安装的版本好像差Release挺多的。
  ) +
  scale_fill_manual(values=c("#00A087FF", "#3C5488FF",'#89ab7c'), guide="none") +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF",'#89ab7c'), guide="none")

f2 <- mpse2 %>%
  mp_plot_alpha(
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson)
  )  

f1 / f2

#######################################
mpse2 %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=Group
  )

mpse2

p1 <- mpse2 %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Group, 
    taxa.class = Phylum, 
    topn = 20,
    relative = TRUE
  )
# visualize the abundance (rarefied) of top 20 phyla for each sample.
p2 <- mpse2 %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Group,
    taxa.class = Phylum,
    topn = 20,
    relative = FALSE
  )
p1 / p2  #p1和p2的区别仅在于一个是绝对值一个是相对值，这部分是物种比例的统计，很好理解

#
h1 <- mpse2 %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Group,
    taxa.class = Phylum,
    relative = TRUE,
    topn = 20,
    geom = 'heatmap',
    features.dist = 'euclidean',
    features.hclust = 'average',
    sample.dist = 'bray',
    sample.hclust = 'average'
  ) 

h2 <- mpse2 %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Group,
    taxa.class = Phylum,
    relative = FALSE,
    topn = 20,
    geom = 'heatmap',
    features.dist = 'euclidean',
    features.hclust = 'average',
    sample.dist = 'bray',
    sample.hclust = 'average'
  )
# the character (scale or theme) of figure can be adjusted by set_scale_theme
# refer to the mp_plot_dist
aplot::plot_list(gglist=list(h1, h2), tag_levels="A")  #这里没有给放回热图，我猜也是包版本的原因


##########################################################
p3 <- mpse2 %>%
  mp_plot_abundance(
    .abundance=RareAbundance, 
    .group=Group,
    taxa.class = Phylum,
    topn = 20,
    plot.group = TRUE
  )

# visualize the abundance of top 20 phyla for each .group (time)
p4 <- mpse2 %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group= Group,
    taxa.class = Phylum,
    topn = 20,
    relative = FALSE,  
    plot.group = TRUE
  )
p3 / p4  #这部分是组间的比较

##############################################
mpse2 %<>% 
  mp_decostand(.abundance=Abundance)
mpse2

mpse2 %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
mpse2

p1 <- mpse2 %>% mp_plot_dist(.distmethod = bray)
p1

p2 <- mpse2 %>% mp_plot_dist(.distmethod = bray, .group = Group)
# The scale or theme of dot heatmap plot can be adjusted using set_scale_theme function.
p2 %>% set_scale_theme(
  x = scale_fill_manual(
    values=c("orange", "deepskyblue"), 
    guide = guide_legend(
      keywidth = 1, 
      keyheight = 0.5, 
      title.theme = element_text(size=8),
      label.theme = element_text(size=6)
    )
  ), 
  aes_var = time # specific the name of variable 
) %>%
  set_scale_theme(
    x = scale_color_gradient(
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  ) %>%
  set_scale_theme(
    x = scale_size_continuous(
      range = c(0.1, 3),
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  )  #set_scale_theme这个功能找不到


p3 <- mpse2 %>% mp_plot_dist(.distmethod = bray, .group = Group, group.test=TRUE, textsize=2)
p3 #这部分在比较样本之间的，类似于相似性之类的东西？

############################################
mpse2 %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
# The dimensions of ordination analysis will be added the colData slot (default).
mpse2

mpse2 %<>%
  mp_adonis(.abundance=hellinger, .formula=~Group, distmethod="bray", permutations=9999, action="add")
mpse2 %>% mp_extract_internal_attr(name=adonis)

library(ggplot2)
p1 <- mpse2 %>%
  mp_plot_ord(
    .ord = pcoa, 
    .group = Group, 
    .color = Group, 
    .size = 1.2,
    .alpha = 1,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse
  ) +
  scale_fill_manual(values=c("#00A087FF", "#3C5488FF",'#89ab7c')) +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF",'#89ab7c')) 

# The size of point also can be mapped to other variables such as Observe, or Shannon 
# Then the alpha diversity and beta diversity will be displayed simultaneously.
p2 <- mpse2 %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Group, 
    .color = Group, 
    .size = Observe, 
    .alpha = Shannon,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#00A087FF", "#3C5488FF",'#89ab7c'), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#00A087FF", "#3C5488FF",'#89ab7c'),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )
p1 + p2 #PCA，为什么少个点？

###############################################
mpse2 %<>%
  mp_cal_clust(
    .abundance = hellinger, 
    distmethod = "bray",
    hclustmethod = "average", # (UPGAE)
    action = "add" # action is used to control which result will be returned
  )
mpse2

sample.clust <- mpse2 %>% mp_extract_internal_attr(name='SampleClust')
sample.clust

library(ggtree)
p <- ggtree(sample.clust) + 
  geom_tippoint(aes(color=Group)) +
  geom_tiplab(as_ylab = TRUE) +
  ggplot2::scale_x_continuous(expand=c(0, 0.01))
p

library(ggtreeExtra)
library(ggplot2)
phyla.tb <- mpse2 %>% 
  mp_extract_abundance(taxa.class=Phylum, topn=30)
# The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phyla="label")
phyla.tb

p1 <- p + 
  geom_fruit(
    data=phyla.tb,
    geom=geom_col,
    mapping = aes(x = RelRareAbundanceBySample, 
                  y = Sample, 
                  fill = Phyla
    ),
    orientation = "y",
    #offset = 0.4,
    pwidth = 3, 
    axis.params = list(axis = "x", 
                       title = "The relative abundance of phyla (%)",
                       title.size = 4,
                       text.size = 2, 
                       vjust = 1),
    grid.params = list()
  )
p1 #聚类树

##########################################################################################
#MAPE在做标志物上面的可视化是非常优秀的，但是由于各种各样乱七八糟的原因，一时半会还没体验上




