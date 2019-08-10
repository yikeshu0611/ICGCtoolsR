#' Use lasso to deal colinear in genes
#'
#' @param getData.data dataset from ICGC.getData command, must include donor and exp_seq
#'
#' @param diffExpr.data dataset from command ICGC.diffExpression
#' @param diffExpr.univarCox.data dataset from command ICGC.diffExpression_univarCox
#'
#' @return several genes
#' @export
#'
#' @examples ICGC.diffExpression_Lasso(getData.data,diffExpr.data,diffExpr.univarCox.data)
ICGC.diffExpression_Lasso <- function(getData.data,
                                      diffExpr.data,
                                      diffExpr.univarCox.data){
  #getData.data
  donor2=getData.data$donor[,c("icgc_donor_id","donor_vital_status","donor_survival_time")]
  #diffExpr.data
  diffExpr.data=diffExpr.data$ExpressionMatrix
  splitmatrix=data.frame(matrix(unlist(strsplit(colnames(diffExpr.data),"-")),ncol=3,byrow=TRUE))
  tumordo=diffExpr.data[,splitmatrix[,3]=="Tumor"]
  colnames(tumordo)=splitmatrix[,2][splitmatrix[,3]=="Tumor"]
  tumordo2=limma::avearrays(tumordo)
  #diffExpr.univarCox.data
  diffdata=diffExpr.univarCox.data$data.result
  message("*****Please select significant level*****")
  hhh=list(" < 0.05"," < 0.01"," < 0.001"," < 0.0001")
  res <- svDialogs::dlg_list(choices = hhh,preselect=FALSE,
                             multiple=FALSE)$res
  sig=as.numeric(as.character(gsub("[< ]","",res)))
  Genefilter=diffdata$Gene[as.numeric(as.character(diffdata$Likelihood.ratio.pvalue)) <=sig ] #add param
  Genefilter=gsub("\\.","-",Genefilter)
  #tumor data
  rownames(tumordo2)=gsub("\\.","-",rownames(tumordo2))
  tumordo3=tumordo2[Genefilter,]
  tumordo4=data.frame(t(tumordo3))
  tumordo4$icgc_donor_id=rownames(tumordo4)
  #merge
  mergedata=dplyr::inner_join(donor2,tumordo4,'icgc_donor_id')
  mergedata$donor_vital_status=ifelse(mergedata$donor_vital_status=="deceased",
                                      1,0)
  library("glmnet")
  library("survival")
  x=as.matrix(mergedata[,4:ncol(mergedata)])
  y=data.matrix(Surv(mergedata$donor_survival_time,mergedata$donor_vital_status))

  fit <- glmnet(x, y, family = "cox")
  cvd <- cv.glmnet(x, y, family="cox")

  coef <- coef(fit, s = cvd$lambda.min)
  GenesfromeLasso=rownames(coef)[coef[,1] != 0]
  #km plot for lasso gene
  x2=x[,GenesfromeLasso]
  for (lassoi in 1:length(GenesfromeLasso)) {
    survdata=cbind(time=mergedata$donor_survival_time,status=mergedata$donor_vital_status,x2)
    lassoGenei=GenesfromeLasso[lassoi]
    group=ifelse(survdata[,lassoGenei] >= median(survdata[,lassoGenei]),1,0)
    class(group)
    survdata2=data.frame(cbind(survdata,group))
    fit <- survfit(Surv(time, status) ~ group,data=survdata2)
    P=survminer::ggsurvplot(fit,data = survdata2,pval.size=5,
                            legend.title=lassoGenei,
                            ncensor.plot=TRUE,
                            pval = TRUE,risk.table = TRUE,tables.height = 0.2)
    print(P)
  }

 # filename=paste0("./SurvivalAnalysisForGenes/",gene.i,genenames.i,".tiff")
 # tiff(file=filename,compression = 'lzw',
 #      width=6, height=6, units="in", res=300)

 # dev.off()
  cat("return",length(GenesfromeLasso),"lasso genes and km plots")
  return(GenesfromeLasso)
}




