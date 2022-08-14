#根据clusterProfiler的代码得到全部的GO词条和基因

get_GO_Env <- function () {
  if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
  }
  get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}
get_organism <- getFromNamespace("get_organism", "DOSE")
library('GOSemSim')
build_Anno <- getFromNamespace("build_Anno", "DOSE")
get_GO2TERM_table <- function() {
  GOTERM.df <- get_GOTERM()
  GOTERM.df[, c("go_id", "Term")] %>% unique
}
get_GOTERM <- function() {
  pos <- 1
  envir <- as.environment(pos)
  if (!exists(".GOTERM_Env", envir=envir)) {
    assign(".GOTERM_Env", new.env(), envir)
  }
  GOTERM_Env <- get(".GOTERM_Env", envir = envir)
  if (exists("GOTERM.df", envir = GOTERM_Env)) {
    GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
  } else {
    GOTERM.df <- toTable(GOTERM)
    assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
  }
  return(GOTERM.df)
}

get_GO_data <- function(OrgDb, ont, keytype) {
  GO_Env <- get_GO_Env()
  use_cached <- FALSE
  
  ont2 <- NULL
  if (exists("ont", envir = GO_Env, inherits = FALSE))
    ont2 <- get("ont", envir = GO_Env)
  
  if (exists("organism", envir=GO_Env, inherits=FALSE) &&
      exists("keytype", envir=GO_Env, inherits=FALSE) &&
      !is.null(ont2)) {
    
    org <- get("organism", envir=GO_Env)
    kt <- get("keytype", envir=GO_Env)
    
    if (org == get_organism(OrgDb) &&
        keytype == kt &&
        (ont == ont2 || ont2 == "ALL") &&
        exists("goAnno", envir=GO_Env, inherits=FALSE)) {
      ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
      ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){
      
      use_cached <- TRUE
    }
  }
  
  if (use_cached) {
    goAnno <- get("goAnno", envir=GO_Env)            
    if (!is.null(ont2) && ont2 != ont) { ## ont2 == "ALL"
      goAnno <- goAnno[goAnno$ONTOLOGYALL == ont,]
    } 
  } else {
    OrgDb <- load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (! keytype %in% kt) {
      stop("keytype is not supported...")
    }
    
    kk <- keys(OrgDb, keytype=keytype)
    
    ## --> take too much memory
    ##
    ## goAnno <- suppressMessages(
    ##     select(OrgDb, keys=kk, keytype=keytype,
    ##            columns=c("GOALL", "ONTOLOGYALL")))
    
    ## if (ont == "ALL") {
    ##     GO2GENE <- unique(goAnno[, c(2,1)])
    ## } else {
    ##     GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
    ## }
    
    goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
    if (ont != "ALL") {
      goterms <- goterms[goterms == ont]
    }
    go2gene <- suppressMessages(
      AnnotationDbi::mapIds(OrgDb, keys=names(goterms), column=keytype,
                            keytype="GOALL", multiVals='list')
    )
    goAnno <- stack(go2gene)
    colnames(goAnno) <- c(keytype, "GOALL")
    goAnno <- unique(goAnno[!is.na(goAnno[,1]), ])
    goAnno$ONTOLOGYALL <- goterms[goAnno$GOALL]
    
    assign("goAnno", goAnno, envir=GO_Env)
    assign("keytype", keytype, envir=GO_Env)
    assign("ont", ont, envir = GO_Env)
    assign("organism", get_organism(OrgDb), envir=GO_Env)
  }
  ## if (ont == "ALL") {
  ##     GO2GENE <- unique(goAnno[, c(2,1)])
  ## } else {
  ##     GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
  ## }
  GO2GENE <- unique(goAnno[, c(2,1)])
  
  GO_DATA <- build_Anno(GO2GENE, get_GO2TERM_table())
  
  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
  
  if (!is.null(ont2) && ont2 == "ALL") {
    return(GO_DATA)
  }
  
  goOnt <- goOnt.df[,2]
  names(goOnt) <- goOnt.df[,1]
  assign("GO2ONT", goOnt, envir=GO_DATA)
  
  return(GO_DATA)
}

test<-get_GO_data(org.Mm.eg.db,'ALL','ENSEMBL')

all_go<-list()
for (i in 1:22233){
  all_go[[i]]<-c(test$EXTID2PATHID[[i]])
  names(all_go)[i]<-(names(test$EXTID2PATHID[i]))
} 

#做成list反而不太好处理？
all_go2<-as.data.frame(array(NA,c(22233,2)))
colnames(all_go2)<-c('GENE','GO')

for (i in 1:22233){
  all_go2$GO[i]<-c(paste0(test$EXTID2PATHID[[i]],collapse = ';'))
  all_go2$GENE[i]<-(names(test$EXTID2PATHID[i]))
} 

gene<-all_go2[grep('GO:0050778',all_go2$GO),]
