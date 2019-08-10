#' To do KM analysis for a single gene in ICGC data
#'
#' @param getData.data From command ICGC.getData, must include donor data
#' @param diffExpr.data From command ICGC.diffExpression
#' @param Genes genes that you want to plot km curves
#'
#' @return a survival curve
#' @export
#'
#' @examples ICGC.diffExpression_singleGene.kmAnalysis(getData.data,diffExpr.data,"TP53")
ICGC.diffExpression_singleGene.kmAnalysis <- function(getData.data,
                                diffExpr.data,Genes){
  donor2=getData.data$donor[,c("icgc_donor_id","donor_vital_status","donor_survival_time")]
  expr.p=data.frame(diffExpr.data$ExpressionMatrix)[Genes,]
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
    fit <- survfit(Surv(time, status) ~ group,data=survdata)
    P=survminer::ggsurvplot(fit,data = survdata,pval.size=5,
                              legend.title=genenames.i,
                            ncensor.plot=TRUE,
                            pval = TRUE,risk.table = TRUE,tables.height = 0.2)
    print(P)
    }
    setTxtProgressBar(pb = prgbar, value = gene.i)
  }
  close(prgbar)
}
