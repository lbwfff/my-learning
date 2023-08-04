#一个批量跑孟德尔随机化的代码，想要解决的是批量下载GWAS数据时遇到的网络问题，不确定能不能解决好所以的问题，但是目前试了试还不错


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
