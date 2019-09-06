#' To do Cox regression of genes
#'
#' @param getData.data dataset from command ICGC.getData
#' @param diffExpr.data dataset from command ICGC.diffExpression
#' @param diffExpr.Lasso.data dataset from command ICGC.diffExpression_Lasso
#'
#' @return six lists: coxphdata, coxsum, coxsumHR, coxphforest, RiskLevel.coxph and RiskLevel.forest
#' @export
#'
#' @examples diffExpr.multiVarCox.data=ICGC.diffExpression.multiVarCox(getData.data,diffExpr.data,diffExpr.Lasso.data)
ICGC.diffExpression.multiVarCox<-function(getData.data,diffExpr.data,diffExpr.Lasso.data){
  #getData.data
  donor2=getData.data$donor[,c("icgc_donor_id","donor_vital_status","donor_survival_time")]
  #diffExpr.data
  diffExpr.data1=diffExpr.data$ExpressionMatrix
  splitmatrix=data.frame(matrix(unlist(strsplit(colnames(diffExpr.data1),"-")),ncol=3,byrow=TRUE))
  tumordo=diffExpr.data1[,splitmatrix[,3]=="Tumor"]
  colnames(tumordo)=splitmatrix[,2][splitmatrix[,3]=="Tumor"]
  tumordo2=limma::avearrays(tumordo)
  #diffExpr.Lasso.data
  lassogenedata=data.frame(t(tumordo2[diffExpr.Lasso.data,]))
  lassogenedata$icgc_donor_id=rownames(lassogenedata)
  #merge
  merge1=dplyr::inner_join(donor2,lassogenedata,"icgc_donor_id")
  merge1$donor_vital_status=ifelse(merge1$donor_vital_status=="deceased",1,0)
  #coxph
  coxphdata=merge1[,-1]
  fit=coxph(Surv(donor_survival_time,donor_vital_status) ~ .,
            data = coxphdata)
  fitaic=step(fit,direction = "both")
  message("*****Results based on minimum AIC*****")
  coxsum=summary(fitaic)
  HR=round(coxsum$conf.int[,1],3)
  LCI95=round(coxsum$conf.int[,3],3)
  UCI95=round(coxsum$conf.int[,4],3)
  Pvalue=round(coxsum$coefficients[,5],3)
  coxsumHR=cbind(HR,LCI95,UCI95,Pvalue)
  coxphforest=suppressWarnings(
  survminer::ggforest(model = fitaic,data = coxphdata,fontsize = 1.05,
                      main = "Hazard ratio\nAlive<----- ----->Dead")
  )
  coef=data.frame(row.names = 1:nrow(merge1))
  for (cofi in names(fitaic$coefficients)) {
    coefi=merge1[,cofi]*(fitaic$coefficients[cofi])
    coef=cbind(coef,coefi)
  }
  for (coefh in 1:nrow(merge1)) {
    if (coefh==1){RiskScore=c()}
    RiskScorei=sum(coef[coefh,])
    RiskScore=c(RiskScore,RiskScorei)
  }
  RiskLevel=ifelse(RiskScore >= median(RiskScore), "high", "low")
  coxphdata=cbind(merge1,RiskScore,RiskLevel)
  fitRiskLevel=coxph(survival::Surv(donor_survival_time,donor_vital_status) ~ RiskLevel,coxphdata)
  RiskLevel.coxph=summary(fitRiskLevel)
  survfitrisk=survival::survfit(
    survival::Surv(donor_survival_time,donor_vital_status) ~ RiskLevel,coxphdata)
  RiskLevel.forest=survminer::ggsurvplot(survfitrisk,data = coxphdata,pval.size=5,legend.title="",
                                         ncensor.plot=TRUE,
                                         pval = TRUE,risk.table = TRUE,tables.height = 0.2)
  #roc
  p=ROCit::rocit(coxphdata$RiskScore,coxphdata$donor_vital_status)
  rocplot=plot(p,YIndex = TRUE,legend = FALSE,grid=FALSE)
  ps=capture.output(summary(p))
  auc=round(as.numeric(as.character(gsub(".*:","",ps[5]))),3)
  #p$TPR[abs(p$TPR-p$FPR)==max(abs(p$TPR-p$FPR))]
  #p$FPR[abs(p$TPR-p$FPR)==max(abs(p$TPR-p$FPR))]
  cutoff=coxphdata$RiskScore[abs(p$TPR-p$FPR)==max(abs(p$TPR-p$FPR))]
  text(0.805,0.23,paste("AUC: ",auc),cex = 1.3)
  text(0.805,0.175,paste("Cutoff Value of Risk Score: ",
                        round(cutoff,3)),cex = 1.3)
  result=c(coxphdata=list(coxphdata),coxsum=list(coxsum),
          coxsumHR=list(coxsumHR),coxphforest=list(coxphforest),
          RiskLevel.coxph=list(RiskLevel.coxph),
          RiskLevel.forest=list(RiskLevel.forest))
  scoredata=coxphdata
  scoredata2=scoredata[order(scoredata$RiskScore),]
  rownames(scoredata2)=scoredata2$icgc_donor_id
  x=1:nrow(scoredata2)
  plot(x,scoredata2$RiskScore,
       cex=1,xlab = "Donor",ylab = "Risk Score",pch = 19,
       col=ifelse(scoredata2$RiskScore >= cutoff,"red","green"))
  legend("topleft", inset=.05, title="Risk",
         c("high","low"),box.col = "white",
         pch=c(19, 19), col=c("red", "green"))

  xlocat=(1:nrow(scoredata2))[scoredata2$RiskScore==cutoff]
  abline(h=cutoff,v=xlocat,lty=2)
  text(x = xlocat,y=cutoff+0.5,labels = paste("Cutoff: ",
                                         round(cutoff,3)))
  plot(x,scoredata2$donor_survival_time,
       xlab = "Donor",ylab = "Survival Time",pch = 19,
       col=ifelse(scoredata2$donor_vital_status==1,"red","green"))
  legend("topright", inset=.05, title="Status",
         c("Dead","Alibe"),box.col = "white",
         pch=c(19, 19), col=c("red", "green"))
  abline(v=xlocat,lty=2)
  if (min(scoredata2[,rownames(coxsumHR)])>0){
    heatmt=t(log2(scoredata2[,rownames(coxsumHR)]))
  }else{
    heatmt=t(scoredata2[,rownames(coxsumHR)])
  }
  scoredata2$RiskLevel= ifelse(scoredata2$RiskScore >= cutoff,"high","low")
  annotation_col = data.frame(
    Risk = scoredata2$RiskLevel
  )
  rownames(annotation_col) = scoredata2$icgc_donor_id
  ann_colors = list(
    Risk = c(low="#00dae0",high="#ff9289")
  )
  pheatmap::pheatmap(heatmt,cluster_cols = FALSE,
                     annotation_col = annotation_col,
                     annotation_names_col = FALSE,
                     annotation_colors = ann_colors,
                     color = colorRampPalette(c("green", "black", "red"))(50))

  cat("\nreturn six lists: coxphdata, coxsum, coxsumHR, coxphforest, RiskLevel.coxph and RiskLevel.forest\n")
  return(result)
}
