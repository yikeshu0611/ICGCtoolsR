#' To do KM analysis for ICGC diffExpression data
#'
#' @param getData.data From command ICGC.getData, must include donor data
#' @param diffExpr.data From command ICGC.diffExpression
#' @param p.filter the p value of significant difference, EX:0.05, defaulted.
#'
#' @return a dataframe and survival curves end with tiff
#' @export
#'
#' @examples ICGC.diffExpression_kmAnalysis(getData.data,diffExpr.data)
ICGC.diffExpression_kmAnalysis <- function(getData.data,diffExpr.data,p.filter=0.05){
  donor2=getData.data$donor[,c("icgc_donor_id","donor_vital_status","donor_survival_time")]
  expr.p=data.frame(diffExpr.data$ExpressionMatrix)
  #extract tumor donor id
  icgc_donor_id.p=data.frame(matrix(unlist(strsplit(colnames(diffExpr.data$ExpressionMatrix),"-")),
         byrow=TRUE,ncol=3))[,2]
  TumorNormal=data.frame(matrix(unlist(strsplit(colnames(diffExpr.data$ExpressionMatrix),"-")),
                                byrow=TRUE,ncol=3))[,3]
  icgc_donor_id=icgc_donor_id.p[TumorNormal=="Tumor"]
  #take out tumor data
  expr=expr.p[,TumorNormal=="Tumor"]
  colnames(expr)=icgc_donor_id
  expr2=data.frame(limma::avearrays(expr))
  expr.t=data.frame(t(expr2))
  expr.t$icgc_donor_id=rownames(expr.t)
  mergedoge=dplyr::inner_join(donor2,expr.t,'icgc_donor_id')
  genenames=colnames(expr.t)[-ncol(expr.t)]
  cat("do km analysis for tumor donor\n")
  cat(paste("All donor:",dim(expr.p)[2]),"\n")
  cat(paste("Tumor donor:",dim(expr2)[2]),"\n")
  cat(paste("diff Genes:",dim(expr2)[1]),"\n")
  prgbar<- txtProgressBar(min = 0, max = length(genenames),
                          style = 3,
                          initial = 0,width = 20)
  for (gene.i in 1:length(genenames)) {
    if (gene.i==1){p=c()}
    #prepare data
    genenames.i=genenames[gene.i]
    group=ifelse(mergedoge[,genenames.i]>median(mergedoge[,genenames.i],
                                                na.rm = TRUE),"1","0")
    time=mergedoge$donor_survival_time
    status=mergedoge$donor_vital_status
    survdata=data.frame(cbind(time,status,group))
    survdata=na.omit(survdata)
    survdata$time=as.numeric(as.character(survdata$time))
    survdata$status=ifelse(survdata$status=="deceased",1,0)
    survdata$group=factor(x = survdata$group,levels = c(1,0),labels = c("high","low"))
    #survival
    if (!any(table(survdata$group)==0)){
    library(survival)
    sdf <- survdiff(Surv(time, status) ~ group,data=survdata)
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    p=c(p,p.val)
    #plot
    if (gene.i==1){diff.surv.gene=0}
    if (p.val<= p.filter){ #add param
      diff.surv.gene=diff.surv.gene+1
      if (!any(grepl("SurvivalAnalysisForGenes",list.files()))){
        dir.create("SurvivalAnalysisForGenes")
      }
      fit <- survfit(Surv(time, status) ~ group,data=survdata)
      filename=paste0("./SurvivalAnalysisForGenes/",gene.i,genenames.i,".tiff")
      tiff(file=filename,compression = 'lzw',
           width=6, height=6, units="in", res=300)
      P=survminer::ggsurvplot(fit,data = survdata,pval.size=5,
                              legend.title=genenames.i,
                            ncensor.plot=TRUE,
                            pval = TRUE,risk.table = TRUE,tables.height = 0.2)
      print(P)
      dev.off()
    }
    setTxtProgressBar(pb = prgbar, value = gene.i)
    }
  }
  close(prgbar)
  cat("Genes with significant survival difference:",diff.surv.gene,'\n')
  cat("return a dataframe and survival curves end with tiff\n")
  cat(diff.surv.gene,"pictures were stored in folder: SurvivalAnalysisForGenes")
  result=data.frame(cbind(Gene=genenames,Pvalue=p))
  return(result)
}
