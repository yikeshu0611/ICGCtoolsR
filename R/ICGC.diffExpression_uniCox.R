#' To do KM univariable Cox regression for ICGC diffExpression data
#'
#' @param getData.data From command ICGC.getData, must include donor data
#' @param diffExpr.data From command ICGC.diffExpression
#' @param GeneType Data type of gene. b as binary. c as continuous, defaulted
#'
#' @return a dataframe and a plot
#' @export
#'
#' @examples ICGC.diffExpression_univarCox(getData.data,diffExpr.data)
ICGC.diffExpression_univarCox <- function(getData.data,diffExpr.data,
                                          GeneType="c"){
  sys.language=any(grepl("Chinese",sessionInfo()))
  cat('\n3\n')
  donor2=getData.data$donor[,c("icgc_donor_id","donor_vital_status","donor_survival_time")]
  expr.p=data.frame(diffExpr.data$ExpressionMatrix)
  #expr.p[1:3,1:3]
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
  cat('\n2\n')
  #delete genes with NA
  mergedogeMatrx=as.matrix(mergedoge)
  for (nai in 1:ncol(mergedogeMatrx)) {
    if (nai==1){nacount=c()}
    if (any(is.na(mergedogeMatrx[,nai]))){
      nacount=c(nacount,nai)
    }
  }
  if (length(nacount)==0){
    mergeNOna=mergedoge
  }else{
    mergeNOna=mergedoge[,-nacount]
  }
  cat('\n1\n')
  if (sys.language){
    cat(tmcn::toUTF8('\u5F00\u59CB\u5355\u56E0\u7D20Cox\u56DE\u5F52'))
  }else{
    cat("do univar Cox regression analysis for tumor donor\n")
  }
  cat("\n")
  cat(paste("All donor:",dim(expr.p)[2]),"\n")
  cat(paste("Tumor donor:",dim(mergeNOna)[1]),"\n")
  cat(paste("Genes:",dim(mergeNOna)[2]),"\n")
  cat("\n")
  genenames=colnames(mergeNOna)[-c(1,2,3)]
  prgbar<- txtProgressBar(min = 0, max = length(genenames),
                          style = 3,
                          initial = 0,width = 20)

  for (gene.i in 1:length(genenames)) {
    if (gene.i==1){HR=c();lci_95=c();uci_95=c();z_score=c()
    deletgene=c();Gene=c();Likelihood.ratio.pvalue=c();genedeletn=c()
    time=mergeNOna$donor_survival_time
    status=mergeNOna$donor_vital_status
    }
    #prepare data
    genenames.i=genenames[gene.i]
    if (GeneType=="b"){
      group=ifelse(mergedoge[,genenames.i]>median(mergedoge[,genenames.i],
                                                  na.rm = TRUE),"1","0")
      survdata=data.frame(cbind(time,status,group))
      survdata=na.omit(survdata)
      survdata$time=as.numeric(as.character(survdata$time))
      survdata$status=ifelse(survdata$status=="deceased",1,0)
      survdata$group=factor(x = survdata$group,levels = c(1,0),labels = c("high","low"))
    }else if (GeneType=="c"){
      survdata=data.frame(cbind(time,status,group=mergeNOna[,genenames.i]))
      survdata=na.omit(survdata)
      survdata$time=as.numeric(as.character(survdata$time))
      survdata$status=ifelse(survdata$status=="deceased",1,0)
      survdata$group=as.numeric(as.character(survdata$group))
    }
    #survival
    suppressMessages(library(survival))
    sdf <- coxph(Surv(time, status) ~ group,data=survdata)
    HRi=summary(sdf)$conf.int[,"exp(-coef)"]
    lci_95i=summary(sdf)$conf.int[,"lower .95"]
    uci_95i=summary(sdf)$conf.int[,"upper .95"]
    z_scorei =summary(sdf)$coefficients[,'z']
    Likelihood.ratio.pvaluei=summary(sdf)$logtest['pvalue']
    if (!any(uci_95i=="Inf",lci_95i=="Inf",
             uci_95i==0,lci_95i==0,
             is.na(uci_95i),is.na(lci_95i))){
      Gene=c(Gene,genenames.i)
      HR=c(HR,HRi)
      lci_95=c(lci_95,lci_95i)
      uci_95=c(uci_95,uci_95i)
      Likelihood.ratio.pvalue=c(Likelihood.ratio.pvalue,
                                Likelihood.ratio.pvaluei)
      z_score=c(z_score,z_scorei)
    }else{
      genedeletn=c(genedeletn,genenames.i)
      deletgene=c(deletgene,genenames.i,",")
    }
    setTxtProgressBar(pb = prgbar, value = gene.i)
  }
  close(prgbar)
  Fresult=list()
  result=data.frame(cbind(Gene,z_score,HR,lci_95,uci_95,
                          Likelihood.ratio.pvalue))
  rownames(result)=NULL
  p1n=length(result$Likelihood.ratio.pvalue[as.numeric(as.character(result$Likelihood.ratio.pvalue))<=0.05])
  p2n=length(result$Likelihood.ratio.pvalue[as.numeric(as.character(result$Likelihood.ratio.pvalue))<=0.01])
  p3n=length(result$Likelihood.ratio.pvalue[as.numeric(as.character(result$Likelihood.ratio.pvalue))<=0.001])
  p4n=length(result$Likelihood.ratio.pvalue[as.numeric(as.character(result$Likelihood.ratio.pvalue))<=0.0001])

  cat("\n")
  cat("p value <= 0.05: ",p1n,"\n")
  cat("p value <= 0.01: ",p2n,"\n")
  cat("p value <= 0.001: ",p3n,"\n")
  cat("p value <= 0.0001: ",p4n,"\n")
  cat("\n")
  z_score[z_score+1.5<0]
  Fresult=c(Fresult,data.result=list(result))
  #plot volcano
  minuslog10pvalue=-log10(as.numeric(as.character(result$Likelihood.ratio.pvalue)))
  HRminus1=as.numeric(as.character(result$HR))-1
  result2=data.frame(cbind(Gene,HRminus1,minuslog10pvalue))
  result2$HRminus1=as.numeric(as.character(result2$HRminus1))
  result2$minuslog10pvalue=as.numeric(as.character(result2$minuslog10pvalue))
  suppressWarnings(library(plotly))
  p<-plot_ly(data = result2,
                  x=~ HRminus1,
                  y=~ minuslog10pvalue,
                  type="scatter",
                  mode="markers",
                  text = paste("Gene:", result2$Gene)

  ) %>% layout(
    title="Volcano for Genes<br>Alive<-----  ----->Dead",
    xaxis = list(title = "HR-1",showgrid = F),
    yaxis = list(title = "-log10(p)",showgrid = F)
  ) %>%
    add_lines(y=-log10(0.05),name = "p=0.05") %>%
    add_lines(y=-log10(0.01),name = "p=0.01") %>%
    add_lines(y=-log10(0.001),name = "p=0.001")
  Fresult =c(Fresult,plot.result=list(p))
  if (length(deletgene) > 0 ){
    cat("Because of Warning message bellow, delet",length(genedeletn),
        "genes:\n")
    cat(toString(deletgene) %>%
      gsub(pattern = " ",replacement = "") %>%
      gsub(pattern = ",{1,}",replacement = "\n"))
  }
  cat("\nreturn a dataframe and a plot\n")
  return(Fresult)
}



