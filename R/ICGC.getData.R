#' To download and import ICGC data
#'
#' @param Release the release verson, EX: release_28 or current
#' @param Project the project of data, EX: ALL-US
#' @param donor a logical param. TRUE: to download this data. FALSE, as defaulted
#' @param exp_array a logical param. TRUE: to download this data. FALSE, as defaulted
#' @param exp_seq a logical param. TRUE: to download this data. FALSE, as defaulted
#' @param specimen a logical param. TRUE: to download this data. FALSE, as defaulted
#'
#' @description This function can be used to download and import ICGC data. Release and Project params can be finded in ICGC web. The other params are all logical words, which False as defaulted.
#' @return a list
#' @export
#' @author Jing Zhang
#'
#' @examples ICGC.getData(Release="current",Project="ALL-US",donor=TRUE)
ICGC.getData <- function(Release,Project,donor=FALSE,specimen=FALSE,exp_array=FALSE,exp_seq=FALSE){  #checklogicalwords
  if (!is.logical(donor)) stop('"donor" must be a logical param!')
  if (!is.logical(exp_array)) stop('"exp_array" must be a logical param!')
  if (!is.logical(exp_seq)) stop('"exp_seq" must be a logical param!')
  if (!is.logical(specimen)) stop('"specimen" must be a logical param!')
  #download and read data
  ICGC.data=list()
  sys.language=any(grepl("Chinese",sessionInfo()))
  if (sys.language){
    if (donor){
      queryword="donor"
      filename=paste0(queryword,".",Project,".tsv")
      destfile=paste0(queryword,".",Project,".tsv.gz")
      if (any(grepl(destfile,list.files()))){
        message(paste0(destfile,tmcn::toUTF8(' \u5DF2\u7ECF\u5B58\u5728,\u8DF3\u8FC7"\u4E0B\u8F7D"\u6B65\u9AA4')))
        message(tmcn::toUTF8('\u89E3\u538B: '),destfile)
        #message(paste0(destfile,tmcn::toUTF8(' \u5DF2\u7ECF\u5B58\u5728,\u8DF3\u8FC7"\u4E0B\u8F7D"\u6B65\u9AA4')))
        #message(tmcn::toUTF8('\u89E3\u538B: '),destfile)
        R.utils::gunzip(destfile)
      }else if (any(grepl(filename,list.files()))){
        message(paste0(filename,tmcn::toUTF8(' \u5DF2\u7ECF\u5B58\u5728,\u8DF3\u8FC7"\u4E0B\u8F7D"\u548C"\u89E3\u538B"\u6B65\u9AA4')))
        #message(paste0(filename,tmcn::toUTF8(' \u5DF2\u7ECF\u5B58\u5728,\u8DF3\u8FC7"\u4E0B\u8F7D"\u548C"\u89E3\u538B"\u6B65\u9AA4')))
      }else{
        url1="https://dcc.icgc.org/api/v1/download?fn=/"
        url=paste0(url1,Release,"/Projects/",Project,"/",queryword,".",Project,".tsv.gz")
        message(tmcn::toUTF8("\u4E0B\u8F7D\u6570\u636E"))
        #message(tmcn::toUTF8("\u4E0B\u8F7D\u6570\u636E"))

        tryCatch({
          suppressWarnings(
            download.file(url,destfile=destfile,quiet = FALSE)
          )
        },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
        message(tmcn::toUTF8('\u89E3\u538B: '),destfile)
        R.utils::gunzip(destfile)
      }
      message(tmcn::toUTF8('\u8BFB\u53D6: '),filename)
      #message(tmcn::toUTF8('\u8BFB\u53D6: '),filename)
      data=as.data.frame(suppressMessages(readr::read_tsv(filename)))
      if (any(is.na(data[,"donor_survival_time"]))){
        data[,"donor_survival_time"][is.na(data[,"donor_survival_time"])]=
          data[,"donor_interval_of_last_followup"][is.na(data[,"donor_survival_time"])]
      }
      data2=data[data$donor_survival_time > 0,]
      data3=data2[!is.na(data2[,"donor_survival_time"]),]
      ICGC.data=c(ICGC.data, donor=list(data3))
      cat(paste0('-----',tmcn::toUTF8('\u6210\u529F\u8BFB\u53D6\u6570\u636E'),'-----'),'\n')
      #cat(paste0('-----',tmcn::toUTF8('\u6210\u529F\u8BFB\u53D6\u6570\u636E'),'-----'),'\n')
    }
    if (specimen){
      queryword="specimen"
      filename=paste0(queryword,".",Project,".tsv")
      destfile=paste0(queryword,".",Project,".tsv.gz")
      if (any(grepl(destfile,list.files()))){
        message(paste0(destfile,tmcn::toUTF8(' \u5DF2\u7ECF\u5B58\u5728,\u8DF3\u8FC7"\u4E0B\u8F7D"\u6B65\u9AA4')))
        message(tmcn::toUTF8('\u89E3\u538B: '),destfile)
        R.utils::gunzip(destfile)
      }else if (any(grepl(filename,list.files()))){
        message(paste0(filename,tmcn::toUTF8(' \u5DF2\u7ECF\u5B58\u5728,\u8DF3\u8FC7"\u4E0B\u8F7D"\u548C"\u89E3\u538B"\u6B65\u9AA4')))
      }else{
        url1="https://dcc.icgc.org/api/v1/download?fn=/"
        url=paste0(url1,Release,"/Projects/",Project,"/",queryword,".",Project,".tsv.gz")
        message(tmcn::toUTF8("\u4E0B\u8F7D\u6570\u636E"))
        tryCatch({
          suppressWarnings(
            download.file(url,destfile=destfile,quiet = FALSE)
          )
        },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
        message(tmcn::toUTF8('\u89E3\u538B: '),destfile)
        R.utils::gunzip(destfile)
      }
      message(tmcn::toUTF8('\u8BFB\u53D6: '),filename)
      data=as.data.frame(suppressMessages(readr::read_tsv(filename)))
      ICGC.data=c(ICGC.data, specimen=list(data))
      cat(paste0('-----',tmcn::toUTF8('\u6210\u529F\u8BFB\u53D6\u6570\u636E'),'-----'),'\n')
    }
    if (exp_array){
      queryword="exp_array"
      filename=paste0(queryword,".",Project,".tsv")
      destfile=paste0(queryword,".",Project,".tsv.gz")
      if (any(grepl(destfile,list.files()))){
        message(paste0(destfile,tmcn::toUTF8(' \u5DF2\u7ECF\u5B58\u5728,\u8DF3\u8FC7"\u4E0B\u8F7D"\u6B65\u9AA4')))
        message(tmcn::toUTF8('\u89E3\u538B: '),destfile)
        R.utils::gunzip(destfile)
      }else if (any(grepl(filename,list.files()))){
        message(paste0(filename,tmcn::toUTF8(' \u5DF2\u7ECF\u5B58\u5728,\u8DF3\u8FC7"\u4E0B\u8F7D"\u548C"\u89E3\u538B"\u6B65\u9AA4')))
      }else{
        url1="https://dcc.icgc.org/api/v1/download?fn=/"
        url=paste0(url1,Release,"/Projects/",Project,"/",queryword,".",Project,".tsv.gz")
        message(tmcn::toUTF8("\u4E0B\u8F7D\u6570\u636E"))

        tryCatch({
          suppressWarnings(
            download.file(url,destfile=destfile,quiet = FALSE)
          )
        },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
        message(tmcn::toUTF8('\u89E3\u538B: '),destfile)
        R.utils::gunzip(destfile)
      }
      message(tmcn::toUTF8('\u8BFB\u53D6: '),filename)
      data=as.data.frame(suppressMessages(readr::read_tsv(filename)))
      if (all(substring(data$gene_id,1,4)=='ENSG')){
        cat(tmcn::toUTF8('\u8F6C\u6362\u57FA\u56E0\u540D:\u4ECEgene symbol\u5230gene id\n'))
        suppressMessages(library(org.Hs.eg.db))
        colnames(data)[grep('gene_id',colnames(data))]='Ensemblname'
        data$ensembl_id=unlist(stringr::str_split(data$Ensemblname,"[.]",simplify=T))[,1]
        g2s=data.frame(org.Hs.egSYMBOL)
        g2e=data.frame(org.Hs.egENSEMBL)
        data2=dplyr::inner_join(data,g2e,'ensembl_id')
        data3=dplyr::inner_join(data2,g2s,'gene_id')
        dotafter=unlist(stringr::str_split(data3$Ensemblname,"[.]",simplify=T))[,2]
        data3$gene_id=paste0(data3$symbol,'.',dotafter)
        ICGC.data=c(ICGC.data, exp_seq=list(data3))
      }else{
        ICGC.data=c(ICGC.data, exp_array=list(data))
      }
      cat(paste0('-----',tmcn::toUTF8('\u6210\u529F\u8BFB\u53D6\u6570\u636E'),'-----'),'\n')
    }
    if (exp_seq){
      queryword="exp_seq"
      filename=paste0(queryword,".",Project,".tsv")
      destfile=paste0(queryword,".",Project,".tsv.gz")
      if (any(grepl(destfile,list.files()))){
        message(paste0(destfile,tmcn::toUTF8(' \u5DF2\u7ECF\u5B58\u5728,\u8DF3\u8FC7"\u4E0B\u8F7D"\u6B65\u9AA4')))
        message(tmcn::toUTF8('\u89E3\u538B: '),destfile)
        R.utils::gunzip(destfile)
      }else if (any(grepl(filename,list.files()))){
        message(paste0(filename,tmcn::toUTF8(' \u5DF2\u7ECF\u5B58\u5728,\u8DF3\u8FC7"\u4E0B\u8F7D"\u548C"\u89E3\u538B"\u6B65\u9AA4')))
      }else{
        url1="https://dcc.icgc.org/api/v1/download?fn=/"
        url=paste0(url1,Release,"/Projects/",Project,"/",queryword,".",Project,".tsv.gz")
        message(tmcn::toUTF8("\u4E0B\u8F7D\u6570\u636E"))

        tryCatch({
          suppressWarnings(
            download.file(url,destfile=destfile,quiet = FALSE)
          )
        },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
        message(tmcn::toUTF8('\u89E3\u538B: '),destfile)
        R.utils::gunzip(destfile)
      }
      message(tmcn::toUTF8('\u8BFB\u53D6: '),filename)
      data=as.data.frame(suppressMessages(readr::read_tsv(filename)))
      ICGC.data=c(ICGC.data, exp_seq=list(data))
      cat(paste0('-----',tmcn::toUTF8('\u6210\u529F\u8BFB\u53D6\u6570\u636E'),'-----'),'\n')
    }
  }else{
    if (donor){
      queryword="donor"
      filename=paste0(queryword,".",Project,".tsv")
      destfile=paste0(queryword,".",Project,".tsv.gz")
      if (any(grepl(destfile,list.files()))){
        message(paste0(destfile," already exit, skip download step."))
        message("Unzip data: ",destfile)
        R.utils::gunzip(destfile)
      }else if (any(grepl(filename,list.files()))){
        message(paste0(filename," already exit, skip download and gunzip steps."))
      }else{
        url1="https://dcc.icgc.org/api/v1/download?fn=/"
        url=paste0(url1,Release,"/Projects/",Project,"/",queryword,".",Project,".tsv.gz")
        message("Download Data")

        tryCatch({
          suppressWarnings(
            download.file(url,destfile=destfile,quiet = FALSE)
          )
        },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
        message("Unzip data: ",destfile)
        R.utils::gunzip(destfile)
      }
      message("Read data: ",filename)
      data=as.data.frame(suppressMessages(readr::read_tsv(filename)))
      if (any(is.na(data[,"donor_survival_time"]))){
        data[,"donor_survival_time"][is.na(data[,"donor_survival_time"])]=
          data[,"donor_interval_of_last_followup"][is.na(data[,"donor_survival_time"])]
      }
      data2=data[data$donor_survival_time > 0,]
      data3=data2[!is.na(data2[,"donor_survival_time"]),]
      ICGC.data=c(ICGC.data, donor=list(data3))
      cat("---succeed---\n")
    }
    if (specimen){
      queryword="specimen"
      filename=paste0(queryword,".",Project,".tsv")
      destfile=paste0(queryword,".",Project,".tsv.gz")
      if (any(grepl(destfile,list.files()))){
        message(paste0(destfile," already exit, skip download step."))
        message("Unzip data: ",destfile)
        R.utils::gunzip(destfile)
      }else if (any(grepl(filename,list.files()))){
        message(paste0(filename," already exit, skip download and gunzip steps."))
      }else{
        url1="https://dcc.icgc.org/api/v1/download?fn=/"
        url=paste0(url1,Release,"/Projects/",Project,"/",queryword,".",Project,".tsv.gz")
        message("Download Data")

        tryCatch({
          suppressWarnings(
            download.file(url,destfile=destfile,quiet = FALSE)
          )
        },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
        message("Unzip data: ",destfile)
        R.utils::gunzip(destfile)
      }
      message("Read data: ",filename)
      data=as.data.frame(suppressMessages(readr::read_tsv(filename)))
      ICGC.data=c(ICGC.data, specimen=list(data))
      cat("---succeed---\n")
    }
    if (exp_array){
      queryword="exp_array"
      filename=paste0(queryword,".",Project,".tsv")
      destfile=paste0(queryword,".",Project,".tsv.gz")
      if (any(grepl(destfile,list.files()))){
        message(paste0(destfile," already exit, skip download step."))
        message("Unzip data: ",destfile)
        R.utils::gunzip(destfile)
      }else if (any(grepl(filename,list.files()))){
        message(paste0(filename," already exit, skip download and gunzip steps."))
      }else{
        url1="https://dcc.icgc.org/api/v1/download?fn=/"
        url=paste0(url1,Release,"/Projects/",Project,"/",queryword,".",Project,".tsv.gz")
        message("Download Data")

        tryCatch({
          suppressWarnings(
            download.file(url,destfile=destfile,quiet = FALSE)
          )
        },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
        message("Unzip data: ",destfile)
        R.utils::gunzip(destfile)
      }
      message("Read data: ",filename)
      data=as.data.frame(suppressMessages(readr::read_tsv(filename)))
      ICGC.data=c(ICGC.data, exp_array=list(data))

      cat("---succeed---\n")
    }
    if (exp_seq){
      queryword="exp_seq"
      filename=paste0(queryword,".",Project,".tsv")
      destfile=paste0(queryword,".",Project,".tsv.gz")
      if (any(grepl(destfile,list.files()))){
        message(paste0(destfile," already exit, skip download step."))
        message("Unzip data: ",destfile)
        R.utils::gunzip(destfile)
      }else if (any(grepl(filename,list.files()))){
        message(paste0(filename," already exit, skip download and gunzip steps."))
      }else{
        url1="https://dcc.icgc.org/api/v1/download?fn=/"
        url=paste0(url1,Release,"/Projects/",Project,"/",queryword,".",Project,".tsv.gz")
        message("Download Data")

        tryCatch({
          suppressWarnings(
            download.file(url,destfile=destfile,quiet = FALSE)
          )
        },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
        message("Unzip data: ",destfile)
        R.utils::gunzip(destfile)
      }
      message("Read data: ",filename)
      data=as.data.frame(suppressMessages(readr::read_tsv(filename)))
      if (all(substring(data$gene_id,1,4)=='ENSG')){
        cat('transform gene Egene symbol to gene id\n')
        suppressMessages(library(org.Hs.eg.db))
        colnames(data)[grep('gene_id',colnames(data))]='Ensemblname'
        data$ensembl_id=unlist(stringr::str_split(data$Ensemblname,"[.]",simplify=T))[,1]
        g2s=data.frame(org.Hs.egSYMBOL)
        g2e=data.frame(org.Hs.egENSEMBL)
        data2=dplyr::inner_join(data,g2e,'ensembl_id')
        data3=dplyr::inner_join(data2,g2s,'gene_id')
        dotafter=unlist(stringr::str_split(data3$Ensemblname,"[.]",simplify=T))[,2]
        data3$gene_id=paste0(data3$symbol,'.',dotafter)
        ICGC.data=c(ICGC.data, exp_seq=list(data3))
      }else{
        ICGC.data=c(ICGC.data, exp_seq=list(data))
      }
      cat("---succeed---\n")
    }
  }

  if (length(ICGC.data)==0){cat('No data was gained. Please check "Release" and "Project" from https://dcc.icgc.org/releases')}
  if (length(ICGC.data)>0){return(ICGC.data)}
}
