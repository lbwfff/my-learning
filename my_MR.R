#一个批量跑孟德尔随机化的代码，想要解决的是批量下载GWAS数据时遇到的网络问题，不确定能不能解决好所以的问题，但是目前试了试还不错
#23-8-4 ps:这个代码extract_outcome失败的时候并没有重新跑，而是当NULL处理的，有时间改一下这个bug

#####################################################################################################################

my_extract_instrument<-function(path,out_id){
  
i<-1

while (i < 2){
  
  print(paste0('start analysis trait name : ', out_id))
  
  catch<-
    tryCatch({
      test<-extract_instruments(outcomes = out_id, #这的outcomes其实是暴露因素吧？
                                p1 = 5e-08, 
                                clump = TRUE, 
                                r2 = 0.001, force_server = TRUE) 
      
      if (is.null(test)) { print(paste0(out_id,' extract NULL')) } else{
      save(test,file=paste0(path,'/ex_in.RData'))
        }
      
      }, 
      
      warning = function(w){
        # print(paste0(out_id,' meet warning need restart'))
        return(NULL)
        }, 
      
      error = function(e){
        print(paste0(out_id,' meet error need restart'))
        return('reanalysis')},
      
      finally = {
        print(paste0('end analysis trait ',out_id))}
    )
  
  if(is.null(catch) ) { i<-(i+1) 
  return('snp_yes')
  } else{
    if( grepl('NULL',catch)) { i<-(i+1) 
    return('snp_no')
    } else{i<i 
      }
    }}
}

my_extract_outcome_data<-function(path,outcomes){
  
  load(paste0(path,'/ex_in.RData'))
  
  i<-1
  
  while (i < 2){
    
    print(paste0('start extract outcome'))
    
    catch<-
      tryCatch({
        test2<-extract_outcome_data(snps=test$SNP, outcomes=outcomes) 
        
        if (is.null(test2)) { print(paste0(outcomes,' extract NULL'))} else{
          save(test2,file=paste0(path,'/ex_out.RData'))
        }
        
      }, 
      
      warning = function(w){
        # print(paste0(outc$trait,' meet warning need restart'))
        return(NULL)
      }, 
      
      error = function(e){
        print(paste0('extract outcome meet error need restart'))
        return('reanalysis')},
      
      finally = {
        print(paste0('end extract outcome'))}
      )
    
    if(is.null(catch) ) { i<-(i+1) 
    return('outcome_yes')
    } else{
      if( grepl('NULL',catch)) { i<-(i+1) 
      return('outcome_no')
      } else{i<i }
    }}
}

my_batch_mr<-function (path,id1,id2){

  inf<-my_extract_instrument(path,id1)
  
  if (inf == 'snp_yes') { 
  inf2<-my_extract_outcome_data(path,id2)
  
   if (inf2 == 'outcome_yes') { 
    
    load(paste0(path,'/ex_in.RData'));load(paste0(path,'/ex_out.RData'))
    
    if(nrow(test2) <3) {print('SNP Num < 3'); return(NULL)}
    
    snpnum<-nrow(test2)
    
    dat <- harmonise_data(
      exposure_dat = test,
      outcome_dat = test2)
    
    res <- mr(dat,method_list = c("mr_ivw","mr_egger_regression",'mr_weighted_median',
                                 'mr_simple_mode','mr_weighted_mode'))
    
    p<-res$pval
    
    OR <-generate_odds_ratios(res)
    OR <-OR$or
    
    het<-mr_heterogeneity(dat)$Q_pval[2]
    
    hor<-mr_pleiotropy_test(dat)$pval
    
    result<-list(snpn=snpnum,p=p,or=OR,het=het,hor=hor)
    return(result)
    
    file.remove(paste0(path,'/ex_in.RData'))
    file.remove(paste0(path,'/ex_out.RData'))
    
    } else{print('skip');return(NULL)}} 
  else{print('skip');return(NULL)}
  
}


#三个function是主要的代码，下面是一个实例
##############################################################

library(here)
library(tidyverse)
library(plyr)
library(TwoSampleMR)
library(ieugwasr)
library(progress)

load('available_outcomes.RData') 

#结局是finn-b-H7_MYOPIA

filao<-ao[ao$population=='European',]
filao<-filao[!is.na(filao$ncase),]
filao<-filao[filao$ncase>2000,]

testcatch<-data.frame(ins=filao$id,out=c('finn-b-H7_MYOPIA'),snp=c(NA),
                      ivw_p=c(NA),egg_p=c(NA),wme=c(NA),sm=c(NA),wmpp=c(NA),
                      or=c(NA),heter_p=c(NA),hori_p=c(NA))

# pb <- progress_bar$new(total = 100)

for (i in 68:120){
  # pb$tick()
  
  out<-my_batch_mr('./',testcatch$ins[i],testcatch$out[i])
  
  
  if(is.null(out)) {print(paste0(testcatch$ins[i],' skip'))} else{
    testcatch[i,3:11]<-c(out$snpn,out$p,out$or[1],out$het,out$hor)
  }
  
}

#遇到了一次网络问题，是顺利重跑了的，但是没遇到300秒的问题，不知道能不能解决

write.csv(testcatch,file = 'testcatch.csv')

################################################################################
#mr_network
#想着用mr得到的因果关系来做network的话，也是完全没问题的，我不太会用igraph，浅试了一下

dis<-read.csv('disease_catch.csv')
imm<-read.csv('immucatch.csv')

imm<-imm[!is.na(imm$hori_p),]
imm<-imm[imm$ivw_p<0.05,]

a1<-c(dis$ins)
b2<-c(imm$ins)
c3<-c('finn-b-H7_MYOPIA')

meta<-data.frame(id=c(dis$ins,imm$ins,'finn-b-H7_MYOPIA'),
                 trait=c(dis$trait,imm$trait,'myopia'))

permt = function (x, y) {
  array<-as.data.frame(array(NA,c(length(x)*length(y)*2,2)))
  array$V1<-c(rep(x,length(y)),rep(y,length(x)))
  array$V2<-c(rep(y,each=length(x)),rep(x,each=length(y)))
  return(array)
} #自己写得一个排列组合function

permt(a1,b2)

test<-rbind(permt(a1,b2),permt(a1,c3),permt(b2,c3)) #设置好因素的配对

testcatch<-data.frame(ins=test$V1,trait=c(NA),trait2=c(NA),
                      out=test$V2,snp=c(NA),
                      ivw_p=c(NA),egg_p=c(NA),wme=c(NA),sm=c(NA),wmpp=c(NA),
                      or=c(NA),heter_p=c(NA),hori_p=c(NA))

testcatch$trait<-meta$trait[match(testcatch$ins,meta$id)]
testcatch$trait2<-meta$trait[match(testcatch$out,meta$id)]


source('my_mr.R')
for (i in 1:nrow(testcatch)){
  
  out<-my_batch_mr('./',testcatch$ins[i],testcatch$out[i])
  
  
  if(is.null(out)) {print(paste0(testcatch$ins[i],' skip'))} else{
    if (is.na(out$het)) {print(paste0(testcatch$ins[i],' skip'))} else{
      testcatch[i,5:13]<-c(out$snpn,out$p,out$or[1],out$het,out$hor)
    }}
  write.csv(testcatch,file = 'new_intakecatch.csv')
} #开始孟德尔分析，使用之前写得function


#OR值决定正负，P值决定强度

network<-testcatch[,c(2,3,6,11)]
network<-network[!is.na(network$or),]
network<-network[network$ivw_p<0.05,] #P值小于0.05认为没有相关性
network$sig<-ifelse(network$or>1,1,-1)
network$score<-(-log2(network$ivw_p))

#igraph绘制网络
edges<-network[,c(1,2,6,5)]
colnames(edges)<-c('from','to','score','sig')
nodes<-data.frame(label=c(unique(c(network$trait,network$trait2))),
                  class=c(NA))
nodes$class<-ifelse(nodes$label %in% dis$trait,'disease',
                    ifelse(nodes$label %in% imm$trait,'protein','myopia'))

library(igraph)

net_pc<-graph_from_data_frame(d=edges,vertices=nodes,directed=TRUE)

plot(net_pc)

deg<-degree(net_pc,mode="all")

vcolor<-c("orange","tomato","lightblue")

V(net_pc)$color<-vcolor[factor(V(net_pc)$class)]

plot(net_pc,vertex.size=1.5*deg,
     vertex.label.cex=1.2,vertex.label.dist=1,
     vertex.label.color='black',vertex.label.dist=0,
     edge.color=ifelse(E(net_pc)$sig==1,'tomato','lightblue'),
     edge.width=E(net_pc)$score/2, edge.arrow.size=.4,
     edge.curved=.1)  
#稍微美化之后的网络，感觉美学上还差点意思，但是想要得到的信息已经完全传递出来了

legend(x=-1.5,y=1.5,levels(factor(V(net_pc)$class)),pch=21,col="#777777",pt.bg=vcolor) #legend


###############################################################################################
#想那么多有的没的干嘛，别人写的快多了

library(TwoSampleMR)
library(doParallel)

rm(list = ls())
exp_data <- readRDS("./exp_data_2023.4.3_p5e8_idALL.rds")  #这个其实是比较关键的数据，目前是靠作者提供的，知道怎么创建这个数据是比较重要的
load('./available_outcomes.RData') #全暴露

# 去掉eqtl的
delete_eqtl <- stringr::str_detect(
  string = exp_data$id.exposure,
  pattern = "eqtl-a", negate = TRUE )
exp_data <- subset(exp_data,delete_eqtl)

# 筛选人群, 你研究的人群确定下
id2 <- subset(ao,ao$population=="European")
length(unique(id2$id))

dplyr::count(id2,population)

exp_data <- subset(exp_data,exp_data$id.exposure %in% id2$id)   #这里文件从处理方法和习惯的不一样，但也还好理解

start1=Sys.time()

outcome=extract_outcome_data(
  snps = unique(exp_data$SNP),
  outcomes = c("ukb-b-6353"),   #你想要的结局
  proxies = TRUE)

save(outcome,file = 'outcomefinn-b-H7_MYOPIA.rda')
end1=Sys.time();end1-start1

exp_data_list=split(exp_data,list(exp_data$id.exposure))

har_loop <- function(exp_data=exp_data_list){
  BBB=TwoSampleMR::harmonise_data(
    exposure_dat = exp_data, outcome_dat = outcome)
  return(BBB)
}

start1=Sys.time()
dat_list <- list()
for(x in 1:length(exp_data_list)){
  dat_list[[names(exp_data_list)[x]]] <- har_loop(exp_data=exp_data_list[[x]])   #批量的harmonise数据
}
end1=Sys.time();end1-start1

dat <- do.call(rbind,dat_list)
dat <- subset(dat,mr_keep)

dat <- split(dat,list(dat$id.exposure,dat$id.outcome))
length(dat)
names(dat) <- paste0("A",1:length(dat))
deleteSNP=names(dat)[sapply(dat,nrow)<3]
length(deleteSNP)

# [1]-去掉list里面小于3个SNP的
for (deleteSNPid in deleteSNP) {
  dat[[deleteSNPid]] <- NULL
}
length(dat)
dat2 <- do.call(rbind,dat)
length(unique(dat2$id.exposure)) 

# 并行
library(doParallel)
#
choose_MR <- function(dat1=dat){ 
  res_hete <- mr_heterogeneity(dat1) 
  if(nrow(res_hete)==1 & !grepl('Invers',res_hete$method[1])){next}
  if(nrow(res_hete)==0 ){next}
  if (res_hete$Q_pval[nrow(res_hete)]<0.05) {
    res=TwoSampleMR::mr(dat1, method_list = c("mr_egger_regression",
                                              "mr_weighted_median", "mr_ivw_mre"))
  } else{
    res=TwoSampleMR::mr(dat1, method_list = c("mr_egger_regression",
                                              "mr_weighted_median", "mr_ivw_fe"))
  }
  return(res)
}  #这里是先做了异质性的检验，然后根据异质性的结果使用不同的MR方法，这个function是分析的本体，需要修改的话可以在这里改

start1=Sys.time()
res_list <- list()
for(x in 1:length(dat)){
  res_list[[names(dat)[x]]] <- choose_MR(dat1=dat[[x]])
}
end1=Sys.time();end1-start1
res <- do.call(rbind,res_list)

write.csv(res,'res_6353.csv')

# 结果提取
res$pval=round(res$pval,3)

res_ALL <- split(res,list(res$id.exposure))
#
judge_1 <- function(mr_res=res2) {
  mr_res$b_direction <- as.numeric(sign(mr_res$b))
  mr_res$b_direction=ifelse(abs(sum(mr_res$b_direction))==3 ,
                            NA,"Inconsistent direction")
  mr_res$p_no <- NA
  mr_res[mr_res$method=="MR Egger","p_no"] <- ifelse(
    mr_res[mr_res$method=="MR Egger","pval"]<0.05," ",
    "MR Egger")
  mr_res[mr_res$method=="Weighted median","p_no"] <- ifelse(
    mr_res[mr_res$method=="Weighted median","pval"]<0.05," ",
    "Weighted median")
  mr_res[grep(x = mr_res$method,pattern = "Inverse variance"),"p_no"] <- ifelse(
    mr_res[grep(x = mr_res$method,pattern = "Inverse variance"),"pval"]<0.05,
    " ","IVW")
  mr_res$p_no <- paste(mr_res$p_no,collapse = " ")
  mr_res$p_no=trimws(mr_res$p_no,which = c("both"))
  return(mr_res)
}

res_ALL=purrr::map(.x =res_ALL,.f = ~judge_1(.x) )
res_ALL2 <- do.call(rbind,res_ALL)
res_ALL3 <- subset(res_ALL2,
                   is.na(res_ALL2$b_direction) )
bool=stringr::str_detect(string =res_ALL3$p_no,
                         pattern = "IVW",negate = TRUE )
res_ALL4 <- subset(res_ALL3,bool)

res_ALL4_1=subset(res_ALL4,select = exposure)
res_ALL4_1 <- unique(res_ALL4_1)
res_ALL4_1[1,1]
res_ALL4_2 <- tidyr::separate(
  data = res_ALL4_1,col = exposure,sep = "\\|",
  into = c("exposure","delete")) %>%
  dplyr::select(-delete)

# 导出
library(openxlsx)

wb <- createWorkbook("My name here")
## Add a worksheets
addWorksheet(wb, "sheet1", gridLines = FALSE)
addWorksheet(wb, "sheet2", gridLines = FALSE)
## write data to worksheet 1
writeData(wb,x = res_ALL4,sheet = "sheet1",
          rowNames = FALSE)
writeData(wb,x = res_ALL4_2,sheet = "sheet2",
          rowNames = FALSE)
## style for body
bodyStyle <- createStyle(border = "TopBottom",
                         bgFill ="#e3e9f4",  
                         fgFill = "#e3e9f4")
a=seq(2,nrow(res_ALL4)+1,6)
b=seq(3,nrow(res_ALL4)+1,6)
c=seq(4,nrow(res_ALL4)+1,6)
d=sort(c(a,b,c))
d
addStyle(wb, sheet = 1, bodyStyle, 
         rows = d,
         cols = 1:11, 
         gridExpand = TRUE)
setColWidths(wb, 1, cols = 1, widths = 21) ## set column width for row names column
## Not run: 
saveWorkbook(wb, "MR_6353.xlsx", 
             overwrite = TRUE)


######################################################################################################################
#记录一下非批量的代码
#单变量就不写了，直接多面量

#之前需要引入三个function代码来自Marina Vabistsevits，三个function都是用来对MVMR（或TwosampleMR）的输入输出进行调整的
library(MVMR)
library(tidyr)
library(tibble)
library(dplyr)

tidy_mvmr_output <- function(mvmr_res) {
  #  tidy up MVMR returned output
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    dplyr::rename(b='Estimate',  #会存在一个库错误
                  se="Std. Error",
                  pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}

tidy_pvals<-function(df){
  # round up output values and keep p-vals in scientific notation
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=4) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 4),
           pval=as.numeric(pval))
}  #这个小数点位数搞那么少干嘛，用来画图的话不是很好

make_mvmr_input <- function(exposure_dat, outcome.id.mrbase=NULL, outcome.data=NULL){
  # provide exposure_dat created in the same way as for TwoSampleMR 
  # also specify the outcome argument [only ONE!] (MR-base ID or full gwas data in .outcome format)
  
  # extract SNPs for both exposures from outcome dataset
  # (for the selected option mr.base or local outcome data)
  if (!is.null(outcome.id.mrbase)) {
    # if mrbase.id is provided
    outcome_dat <- extract_outcome_data(snps = unique(exposure_dat$SNP),
                                        outcomes = outcome.id.mrbase)
  } else if (!is.null(outcome.data)){
    # if outcome df is provided
    outcome_dat <- outcome.data %>% filter(SNP %in% exposure_dat$SNP)
  }
  
  # harmonize datasets
  exposure_dat <- exposure_dat %>% mutate(id.exposure = exposure)
  outcome_harmonised <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  exposures_order <- colnames(outcome_harmonised$exposure_beta)
  
  # Create variables for the analysis 
  
  ### works for many exposures
  no_exp = dim(outcome_harmonised$exposure_beta)[2] # count exposures
  # add beta/se names
  colnames(outcome_harmonised$exposure_beta) <- paste0("betaX", 1:no_exp)
  colnames(outcome_harmonised$exposure_se) <- paste0("seX", 1:no_exp)
  
  XGs <-left_join(as.data.frame(outcome_harmonised$exposure_beta) %>% rownames_to_column('SNP'), 
                  as.data.frame(outcome_harmonised$exposure_se)   %>%rownames_to_column('SNP'), 
                  by = "SNP")
  
  YG <- data.frame(beta.outcome = outcome_harmonised$outcome_beta,
                   se.outcome = outcome_harmonised$outcome_se) %>% 
    mutate(SNP = XGs$SNP)
  
  
  return(list(YG = YG,
              XGs = XGs,
              exposures = exposures_order))
}

#OK基于这三个function我们进行多变量分析
#首先是twosampleMR

library(TwoSampleMR)

id_exposure <-c('ukb-b-2209','ukb-b-1489','ukb-b-2134','ukb-b-4424')
exposure_dat <- mv_extract_exposures(id_exposure,find_proxies = FALSE,
                                     force_server = FALSE)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, 'ukb-d-5843_3')

mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
# mv_lasso<- mv_lasso_feature_selection(mvdat) #此处可用LASSO挑选变量，看需求吧

res <- mv_multiple(mvdat,plots=T)

result_2smr <- res$result %>%  #调整一下结果格式计算OR值
  split_outcome() %>%
  separate(outcome, "outcome", sep="[(]") %>% 
  mutate(outcome=stringr::str_trim(outcome))%>% 
  generate_odds_ratios() %>% 
  select(-id.exposure, -id.outcome) %>% 
  tidy_pvals()

result_2smr$type<-c('TwoSampleMR')

ggplot(result_2smr, aes(y=exposure, x=or, label=outcome, colour=type)) +  #shape=exposure, 
  geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
  geom_point(size=2)+
  # xlim(0,2)+
  scale_color_manual(values=met.brewer("VanGogh1", 3,type ='discrete'))+
  # scale_shape_manual(values = c(19,20,17)) +
  geom_vline(xintercept=1, linetype=3) +
  theme_minimal_hgrid(10, rel_small = 1) +
  facet_wrap(~outcome, ncol=1)+
  labs(color = "",y = "", x = "Odds Ratio",
       title= paste0("Univariate MR results for ",
                     'Myopia' ,", 95% CI") )+
  theme(legend.position = "none")
#绘制森林图，我还挺喜欢这个画法的，虽然没有把P值OR值放上去，但随便吧好看就行

#然后是MVMR
mvmr_input <- make_mvmr_input(exposure_dat, outcome.id.mrbase= 'ukb-d-5843_3') #这个input好像也是用twosamplemr的function做的，只不过把结果调整为了MVMR的input形式

mvmr_out <- format_mvmr(BXGs = mvmr_input$XGs %>% select(contains("beta")),  
                        BYG = mvmr_input$YG$beta.outcome,                     
                        seBXGs = mvmr_input$XGs %>% select(contains("se")),  
                        seBYG = mvmr_input$YG$se.outcome,                     
                        RSID = mvmr_input$XGs$SNP)   

mvmr_res <- ivw_mvmr(r_input=mvmr_out) #多变量分析

result_mvmr <-
  mvmr_res %>% 
  tidy_mvmr_output() %>% 
  mutate(exposure = mvmr_input$exposures,
         outcome = 'Breast cancer') %>% 
  select(exposure, outcome, everything()) %>% 
  tidy_pvals() 

result_mvmr$type<-c('MVMR')

merge<-data.frame(exposure=c(result_2smr$exposure,result_mvmr$exposure),
                  OR=c(result_2smr$or,result_mvmr$or),outcome=c('myopia'),
                  orlow=c(result_2smr$or_lci95,result_mvmr$or_lci95),
                  orup=c(result_2smr$or_uci95,result_mvmr$or_uci95),
                  type=c(result_2smr$type,result_mvmr$type))  #把twosampleMR的结果和MVMR整合到一起，两个output格式改得还不一样需要手动调一下

merge$type<-factor(merge$type,levels = c('TwoSampleMR','MVMR'))

pdf("multi_mr.pdf",width = 6,height = 4)
ggplot(merge, aes(y=exposure, x=OR, label=outcome, colour=type)) +  #shape=exposure, 
  geom_errorbarh(aes(xmin=orlow, xmax=orup), height=.3) +
  geom_point(size=2)+
  # xlim(0.88,1.14)+
  scale_color_manual(values=met.brewer("VanGogh1", 3,type ='discrete'))+
  # scale_shape_manual(values = c(19,20,17)) +
  geom_vline(xintercept=1, linetype=3) +
  theme_minimal_hgrid(10, rel_small = 1) +
  facet_wrap(~type, ncol=1)+
  labs(color = "",y = "", x = "Odds Radio",
       title= paste0("Multivariate MR results for ",
                     'Myopia' ,", 95% CI") )+
  theme(legend.position = "none")
dev.off()

#画在一张图里













