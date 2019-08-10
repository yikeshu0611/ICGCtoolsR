#' To plot heatmap for diff expression genes
#'
#' @param diffExpr.data two list ExpressionMatrix and diffGene from ICGC.diffExpression command,
#'
#' @return a plot
#' @export
#'
#' @examples ICGC.diffExpression_heatmap(diffExpr.data)
ICGC.diffExpression_heatmap <- function(diffExpr.data) {
  SampleMatrix=diffExpr.data$ExpressionMatrix[as.character(diffExpr.data$diffGene[,"Gene"]),]
  if (min(SampleMatrix,na.rm = TRUE)<0){
    log2SampleMatrix=SampleMatrix
  }else{
    log2SampleMatrix=log2(SampleMatrix)
  }
  conamesplit=stringr::str_split(string = colnames(log2SampleMatrix),pattern = "-")
  a=matrix(unlist(conamesplit),byrow=TRUE,ncol=3)[,3]
  Legend=data.frame(a)
  rownames(Legend)=colnames(log2SampleMatrix)
  colnames(Legend)="Type"
  p <- list((log2SampleMatrix),annotation_col = Legend,cluster_cols = FALSE,
                     color = colorRampPalette(c("green", "black", "red"))(50))

  heatp <- do.call(pheatmap::pheatmap,p)
  return(heatp)
}


