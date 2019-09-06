#' To plot Volcano for diff genes
#'
#' @param diffExpr.data dataset from ICGC.diffExpression command
#' @param log2fc.filter log2fc cutoff value
#' @param minus.log10fdr.filter -log10(fdr) cutoff value
#'
#' @return two lists: VolcanoData and VolcanoPlot.
#' @export
#'
#' @examples icgc.Volcano(diffExpr.data,2,10)
ICGC.diffExpression_Volcano <- function(diffExpr.data,log2fc.filter,minus.log10fdr.filter) {
  cat("log2fc.filter range: ",range(biowolf$log2FC)[1],"~",range(biowolf$log2FC)[2])
  cat("\n")
  cat("minus.log10fdr.filter range: ",range(biowolf$minusLog10p)[1],"~",range(biowolf$minusLog10p)[2])

  biowolf<-diffExpr.data$diffGene
  biowolf$minusLog10p=-log10(biowolf$fdr.P.Value)
  biowolf$sig=ifelse(biowolf$minusLog10p > minus.log10fdr.filter &
                     biowolf$log2FC >= abs(log2fc.filter)  ,"up",
                     ifelse(biowolf$minusLog10p > minus.log10fdr.filter &
                              biowolf$log2FC < -abs(log2fc.filter) ,"down","no"))
  result=list()
  result=c(result,VolcanoData=list(biowolf))
  library(ggplot2)
  theme_set(theme_bw())
  volcano<-ggplot(biowolf,aes(x=biowolf$log2FC,y=biowolf$minusLog10p))
  p<-volcano +theme(panel.grid =element_blank())+
    geom_point(aes(color=sig))+
    scale_color_manual(values=c("blue","grey","red"))+
    labs(title="Volcanoplot",x="logFC",y="-log10(FDR)")+
    geom_hline(yintercept=minus.log10fdr.filter,linetype=3)+
    geom_vline(xintercept=c(-log2fc.filter,log2fc.filter),linetype=3)
  result=c(result,VolcanoPlot=list(p))
  cat("return two lists: VolcanoData and VolcanoPlot")
  return(result)
}



