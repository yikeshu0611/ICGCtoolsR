#' To do ROC analysis for a single gene in ICGC data
#'
#' @param getData.data From command ICGC.getData, must include donor data
#' @param diffExpr.data From command ICGC.diffExpression
#' @param Genes genes that you want to plot km curves
#' @param plotROC a logical param. whether to plot roc, TRUE as defaulted
#'
#' @return a ROC curve
#' @export
#'
#' @examples ICGC.diffExpression_singleGene.ROC(getData.data,diffExpr.data,"TP53")
ICGC.diffExpression_singleGene.ROC <- function(getData.data,
                                diffExpr.data,Genes,plotROC=TRUE){
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
  cat("do ROC for tumor donor\n")
  cat(paste("All donor:",dim(expr.p)[2]),"\n")
  cat(paste("Tumor donor:",dim(expr2)[2]),"\n")
  cat(paste("diff Genes:",dim(expr2)[1]),"\n")
  prgbar<- txtProgressBar(min = 0, max = length(genenames),
                          style = 3,
                          initial = 0,width = 20)
  for (gene.i in 1:length(genenames)) {
    if (gene.i==1){auc.w=c();gene.w=c()}
    #prepare data
    genenames.i=genenames[gene.i]
    time=mergedoge$donor_survival_time
    status=mergedoge$donor_vital_status
    survdata=data.frame(cbind(time,status,gene=mergedoge[,genenames.i]))
    survdata=na.omit(survdata)
    survdata$time=as.numeric(as.character(survdata$time))
    survdata$status=ifelse(survdata$status=="deceased",1,0)
    survdata$gene=as.numeric(as.character(survdata$gene))
    p=ROCit::rocit(survdata$gene,survdata$status)
    if (plotROC){
    rocplot=plot(p,YIndex = TRUE,legend = FALSE,grid=FALSE)
    }
    ps=capture.output(summary(p))
    auc=round(as.numeric(as.character(gsub(".*:","",ps[5]))),3)
    auc.w=c(auc.w,auc)
    gene.w=c(gene.w,genenames.i)
    #p$TPR[(p$TPR-p$FPR)==max(abs(p$TPR-p$FPR))]
    #p$FPR[(p$TPR-p$FPR)==max(abs(p$TPR-p$FPR))]
    cutoff.Value=round(mergedoge[,genenames.i][(p$TPR-p$FPR)==max(abs(p$TPR-p$FPR))],3)
    if (plotROC){
      text(0.805,0.23,paste(genenames.i," AUC: ",auc),cex = 1.3)
      text(0.805,0.16,paste(gene.i," Cutoff Value: ",cutoff.Value),cex = 1.3)
    }
    setTxtProgressBar(pb = prgbar, value = gene.i)
  }
  close(prgbar)
  result=data.frame(cbind(gene=gene.w,auc=auc.w))
  result$auc=as.numeric(as.character(result$auc))
  return(result)
}
