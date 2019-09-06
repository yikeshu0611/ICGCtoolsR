#' To do different expression matrix and genes
#'
#' @param getData.data dataset from ICGC.getData command, must include donor and exp_seq
#'
#' @return two lists, ExpressionMatrix and diffGene
#' @export
#'
#' @examples result=ICGC.diffExpression(getData.data)
ICGC.diffExpression <- function(getData.data){
  sys.language=any(grepl("Chinese",sessionInfo()))
  #step1/12
  if (sys.language){
    message(tmcn::toUTF8('\u6B65\u9AA41/12: \u8BFB\u53D6specimean\u6570\u636E,\u5E76\u9009\u62E9normal\u548Ctumor'))
  }else{
    message(paste("Steps 1/12: read specimen data, select normal and tumor part"))
  }
    specimen=getData.data$specimen
    specimen2=specimen[,c("icgc_specimen_id","specimen_type")]
    #select normal and tumor part
    reph=max(nchar(names(table(specimen2$specimen_type))))-nchar(names(table(specimen2$specimen_type)))+5
    rephw=c()
    for (reph.i in reph) {
      rephw=c(rephw,gsub("[, ]","",toString( rep("-",reph.i))))
    }
    hh=paste0(names(table(specimen2$specimen_type)),rephw,
              table(specimen2$specimen_type))
    hhh=hh[order(hh)]
    cat("\n")
    if (sys.language){
      message(tmcn::toUTF8('********\u8BF7\u9009\u62E91\u4E2A\u6216\u591A\u4E2A\u6837\u672C\u4F5C\u4E3Anoraml group********'))
    }else{
      message("********Please Chose one or more as Normal group********")
    }
    res <- svDialogs::dlg_list(choices = hhh,preselect=FALSE,
                               multiple=TRUE)$res
    if (!length(res)) {
      stop("You cancelled the choice\n")
    } else {
      Normal=sub("-----.*","",res)
    }
    cat("\n")
    if (sys.language){
      message(tmcn::toUTF8('********\u8BF7\u9009\u62E91\u4E2A\u6216\u591A\u4E2A\u6837\u672C\u4F5C\u4E3ATumor group********'))
    }else{
      message("********Please Chose one or more as Tumor group********")
    }
    res <- svDialogs::dlg_list(choices = hhh,preselect=FALSE,
                               multiple=TRUE)$res
    if (!length(res)) {
      stop("You cancelled the choice\n")
    } else {
      Tumor=sub("-----.*","",res)
    }
    for (norml.i in 1:length(Normal)) {
      if (norml.i == 1){nt=data.frame()}
      nt=rbind(nt,specimen2$specimen_type==Normal[norml.i])
    }
    for (nti in 1:ncol(nt)) {
      if (nti==1)(ntc=c())
      ntc=c(ntc,any(nt[,nti]))
    }
    NormSelect=specimen2[ntc,]
    NormSelect$tumour_confirmed="Normal"
    for (tumor.i in 1:length(Tumor)) {
      if (tumor.i == 1){tt=data.frame()}
      tt=rbind(tt,specimen2$specimen_type==Tumor[tumor.i])
    }
    for (tti in 1:ncol(tt)) {
      if (tti==1)(ttc=c())
      ttc=c(ttc,any(tt[,tti]))
    }
    TumorSelect=specimen2[ttc,]
    TumorSelect$tumour_confirmed="Tumor"
    specimen2=rbind(NormSelect,TumorSelect)
    NTtable=table(specimen2$tumour_confirmed)
    if (sys.language){
      cat("\n",tmcn::toUTF8('          specimen\u6570\u636E:'),"\n",
          "           ",names(NTtable)[1],":",NTtable[1],"\n",
          "           ",names(NTtable)[2],":",NTtable[2],"\n")
    }else{
      cat("\n","          The count in specimen data","\n",
          "           ",names(NTtable)[1],":",NTtable[1],"\n",
          "           ",names(NTtable)[2],":",NTtable[2],"\n")
    }
    #read seq expression data
    if (sys.language){
      message(tmcn::toUTF8('\u6B65\u9AA42/12: \u8BFB\u53D6expression\u6570\u636E'))
    }else{
      message("Steps 2/12: read expression data")
    }
    ExprValue=getData.data$exp_seq
    if (is.null(ExprValue)){
      ExprValue=getData.data$exp_array
    }
    #whether normalized_expression_value
    colnames(ExprValue)[grep("normalized_expression_",
                             colnames(ExprValue))]="normalized_read_count"
    ExprValue2=ExprValue[,c("icgc_donor_id","icgc_specimen_id",
                            "gene_id","normalized_read_count" )]
    if (min(ExprValue2$normalized_read_count)==0){
      ExprValue2$normalized_read_count=ExprValue2$normalized_read_count+0.001}
    if (sys.language){
      message(tmcn::toUTF8('\u6B65\u9AA43/12: \u5408\u5E76specimen\u548Cexpression'))
    }else{
      message("Steps 3/12: inner join specimen and exp_seq data")
    }
    #merge specimen and expression setable
    MergeData=dplyr::inner_join(ExprValue2,specimen2,"icgc_specimen_id")
    SEtable=table(MergeData$tumour_confirmed)
  #check
  if (length(SEtable)==1){
    if (sys.language){
      SEstop=paste(tmcn::toUTF8('expression\u4E2D\u7684\u6570\u636E\u53EA\u5305\u542Btumor\u90E8\u5206,\u56E0\u6B64\u5408\u5E76\u540E\u7684\u6570\u636E\u4E2D\u4E0D\u5305\u542Bnormal\u90E8\u5206\u7684\u4FE1\u606F'),
                   "\n    ", names(SEtable),SEtable,"\n",
                   "\u5C06\u4E0D\u4F1A\u8FDB\u884C\u57FA\u56E0\u5728tumor\u548Cnormal\u4E4B\u95F4\u7684\u5DEE\u5F02\u8868\u8FBE\u5206\u6790,\u4E0D\u4F1A\u8FDB\u884C9-12\u6B65\u9AA4")
      cat(SEstop,'\n')
    }else{
      SEstop=paste("Only One group from specimen data was tested in exp_seq data:",
                   "\n    ",names(SEtable),SEtable,"\n",
                   "Differences between Normal and Tumor can not be carried out\n")
      cat(SEstop)
    }
  }else{
    if (sys.language){
      cat("\n",tmcn::toUTF8("          \u5408\u5E76\u7684\u6570\u636E\u5305\u542B\u4E86normal\u548Ctumor\u4E24\u90E8\u5206,\u53EF\u8FDB\u884C\u57FA\u56E0\u8868\u8FBE\u5DEE\u5F02\u5206\u6790"),"\n",
          "           ",names(SEtable)[1],":",SEtable[1],"\n",
          "           ",names(SEtable)[2],":",SEtable[2],"\n")
    }else{
    cat("\n","          The count after MeanExpressionFilter in merged data","\n",
        "           ",names(SEtable)[1],":",SEtable[1],"\n",
        "           ",names(SEtable)[2],":",SEtable[2],"\n")
    }
  }
  #matrix key merge id
    if (sys.language){
      message(tmcn::toUTF8('\u6B65\u9AA44/12: \u751F\u6210\u8868\u8FBE\u77E9\u9635\u7684\u884C\u5217\u540D'))
    }else{
      message("Steps 4/12: generate expression matrix key id")
    }
  MergeData$MetrixColname=stringr::str_c(MergeData$icgc_specimen_id,'-',MergeData$icgc_donor_id,'-',MergeData$tumour_confirmed,";;",MergeData$gene_id)
  #matrix data in dataframe
  MergeData1=MergeData[,c("MetrixColname","normalized_read_count")]
  #rm(MergeData)
  #duplicate dealing
  if (any(duplicated(MergeData1$MetrixColname))){
    if (sys.language){
      message(tmcn::toUTF8('\u884C\u5217\u540D\u4E2D\u6709\u91CD\u590D\u9879,\u4F7F\u7528limma\u5305\u4E2Davearrays\u5904\u7406'))
    }else{
      message("Steps 5/12: limma::avearrays")
    }
    DupDeal=matrix(MergeData1$normalized_read_count,nrow = 1)
    colnames(DupDeal)=MergeData1$MetrixColname
    DupDeal1=limma::avearrays(DupDeal)
    DupDeal1[DupDeal1=="NaN"]=NA
    #rm(DupDeal)
    dim(DupDeal1)
    DupDeal1[,1]
    head(DupDeal1)
    DupDeal2=data.frame(dplyr::bind_cols(MetrixColname=colnames(DupDeal1),
                                         normalized_read_count=DupDeal1[1,]))
    ##rm(DupDeal1)
    rownames(DupDeal2)=NULL
  }else{
    if (sys.language){
      message(tmcn::toUTF8('\u6B65\u9AA45/12: \u884C\u5217\u540D\u4E2D\u65E0\u91CD\u590D\u9879'))
    }else{
      message("Steps 5/12: No duplicated key")
    }
    DupDeal2=MergeData1
  }
  #gene_id,MetrixColname and normalized_read_count
  if (sys.language){
    message(tmcn::toUTF8('\u6B65\u9AA46/12: \u751F\u6210\u8868\u8FBE\u77E9\u9635(\u6B64\u6B65\u9AA4\u9700\u8017\u8D39\u8F83\u957F\u65F6\u95F4)'))
  }else{
    message("Steps 6/12: generate expression matrix")
  }
  colnamesplit=stringr::str_split(string = DupDeal2[,"MetrixColname"],pattern = ";;")
  MetrixColname.t=matrix(unlist(colnamesplit),byrow=TRUE,ncol=2)[,1]
  gene_id.t=matrix(unlist(colnamesplit),byrow=TRUE,ncol=2)[,2]
  normalized_read_count.t=DupDeal2[,"normalized_read_count"]
  #data frame
  MergeNoDup=data.frame(row.names = 1:length(MetrixColname.t))
  MergeNoDup$MetrixColname=MetrixColname.t
  MergeNoDup$gene_id=gene_id.t
  normalized_read_count.t2=as.numeric(as.character(normalized_read_count.t))
  MergeNoDup$normalized_read_count=normalized_read_count.t2
  #data matrix
  Matrix=tidyr::spread(MergeNoDup,gene_id,normalized_read_count,fill = NA)
  Matrix2=Matrix[,-1]
  rownames(Matrix2)=Matrix$MetrixColname
  if (sys.language){
    cat(tmcn::toUTF8('           \u8868\u8FBE\u77E9\u9635\u7684\u884C\u5217\u5927\u5C0F'),dim(Matrix2)[1],"*",dim(Matrix2)[2],"\n")
  }else{
    cat(paste("           The dimension of Expression Matrix:",dim(Matrix2)[1],"*",dim(Matrix2)[2]),"\n")
  }
  NormalPart=Matrix2[grepl("Normal",rownames(Matrix2)),]
  TumorPart=Matrix2[grepl("Tumor",rownames(Matrix2)),]
  if (nrow(NormalPart)==0){
    aveMatrix=t(TumorPart)
  }else{
    aveMatrix=t(rbind(NormalPart,TumorPart))
  }
  rownames(aveMatrix)=colnames(TumorPart)
  if (sys.language){
    message(tmcn::toUTF8('\u6B65\u9AA47/12: \u8BBE\u7F6EMeanExpressionFilter\u7684\u6700\u5C0F\u503C,\u53BB\u9664\u5E73\u5747\u8868\u8FBE\u503C\u8FC7\u4F4E\u7684\u57FA\u56E0'))
  }else{
    message("Steps 7/12: set filter of mean expression")
  }
  diamessage=paste0("MeanExpressionFilter >= \n",
                  min(aveMatrix,na.rm = TRUE),
                  "~",
                  max(aveMatrix,na.rm = TRUE))
  res<-svDialogs::dlg_input(message = diamessage)$res
  if (!length(res)){
    stop("You cancelled the choice\n")
  }else{
    MeanExpression=as.numeric(res)
    cat("           MeanExpressionFilter >=",MeanExpression)
  }
  aveMatrix2=aveMatrix[rowMeans(aveMatrix,na.rm = TRUE)>=MeanExpression,]#add param
  TumorNormal=ifelse(grepl("Tumor",colnames(aveMatrix2)),"tumor","normal")
  MeanFiltertable=table(TumorNormal)
  if (sys.language){
    cat("\n",tmcn::toUTF8("          \u7B5B\u9009\u540E\u7684\u77E9\u9635\u5927\u5C0F:"),
        dim(aveMatrix2)[2],'*',dim(aveMatrix2)[1],"\n")
  }else{
    cat("\n","          The size after MeanExpressionFilter in merged data",dim(aveMatrix2)[2],'*',dim(aveMatrix2)[1],"\n")
  }
  if (sys.language){
    message(tmcn::toUTF8('\u6B65\u9AA48/12: \u5220\u9664tumor\u6216\u8005normal\u4E2D\u65E0\u6570\u636E\u7684\u57FA\u56E0'))
  }else{
    message("Steps 8/12: delet empty part in tumor or normal")
  }
  aveMatrixTumor=aveMatrix2[,grepl("Tumor",colnames(aveMatrix2))]
  aveMatrixTumor[is.na(aveMatrixTumor)]=0
  aveMatrixNormal=aveMatrix2[,!grepl("Tumor",colnames(aveMatrix2))]
  aveMatrixNormal[is.na(aveMatrixNormal)]=0
  if (nrow(NormalPart)==0){
    aveMatrix3=aveMatrix2[rowMeans(aveMatrixTumor)!=0,]
  }else{
    aveMatrix3=aveMatrix2[rowMeans(aveMatrixNormal)!=0 & rowMeans(aveMatrixTumor)!=0,]

  }
  if (sys.language){
    cat(tmcn::toUTF8('          \u5220\u9664\u540E\u7684\u77E9\u9635\u5927\u5C0F:'),
        dim(aveMatrix3)[2],"*",dim(aveMatrix3)[1],"\n")
  }else{
    cat(paste("           The dimension of matrix after delet empty part:",
              dim(aveMatrix3)[2],"*",dim(aveMatrix3)[1]),"\n")
  }
  result=list()
  result=c(result,ExpressionMatrix=list(aveMatrix3))
  #check
  if (length(SEtable) > 1){
    if (sys.language){
      tmcn::toUTF8('\u6B65\u9AA49/12: \u4F7F\u7528wilcox\u7684t\u68C0\u9A8C\u8FDB\u884C\u5DEE\u5F02\u5206\u6790')
    }else{
      message("Steps 9/12: do diff expression by wilcox.test")
    }

  prgbar<- txtProgressBar(min = 0, max = nrow(aveMatrix3),
                        style = 3,
                        initial = 0,width = 20)
  for (i in 1:nrow(aveMatrix3)) {
    if (i==1){
      TumorMed=c()
      NormalMed=c()
      MedDiff=c()
      log2FC=c()
      wilcox.P.Value=c()
    }
    #temp data
    #i=14541
    Matrix.i=aveMatrix3[i,]
    Matrix.i2=Matrix.i[!is.na(Matrix.i)]
    Tumor.i=Matrix.i2[grepl("Tumor",names(Matrix.i2))]
    Normal.i=Matrix.i2[grepl("Normal",names(Matrix.i2))]
    if (length(Normal.i)<=0 | length(Tumor.i) <=0){#add param
      TumorMed=c(TumorMed,NA)
      NormalMed=c(NormalMed,NA)
      MedDiff=c(MedDiff,NA)
      log2FC=c(log2FC,NA)
      wilcox.P.Value=c(wilcox.P.Value,NA)
    }else{
      #med
      TumorMed=c(TumorMed,round(median(Tumor.i,na.rm = TRUE),4))
      NormalMed=c(NormalMed,round(median(Normal.i,na.rm = TRUE),4))
      MedDiff=TumorMed-NormalMed
      #p value
      TumorNormal3=ifelse(grepl("Tumor",names(Matrix.i2)),"tumor","normal")
      wilcox.P.Value=c(wilcox.P.Value,wilcox.test(Matrix.i2~TumorNormal3)$p.value)
      #log2 FC
      log2FC=c(log2FC,round(log2(mean(Tumor.i)/mean(Normal.i)),4))
    }
    setTxtProgressBar(pb = prgbar, value = i)
  }
  close(prgbar)
  Part1=data.frame(cbind(Gene=rownames(aveMatrix3),
                       TumorMed,NormalMed,
                       MedDiff,log2FC,wilcox.P.Value))
  Part2=na.omit(Part1)
  Part3=Part2[all(as.numeric(as.character(Part2$log2FC))>0,as.numeric(as.character(Part2$MedDiff))>0),]
  Part4=Part2[all(as.numeric(as.character(Part3$log2FC))<0,as.numeric(as.character(Part3$MedDiff))<0),]
  if (sys.language){
    message(tmcn::toUTF8('\u4F7F\u7528fdrtool\u5305\u4E2D\u7684fdrtool\u547D\u4EE4\u8BA1\u7B97\u77EB\u6B63p\u503C'))
  }else{
    message("Steps 10/12: fdrtool::fdrtool")
  }

  Fdr = fdrtool::fdrtool(plot = FALSE,verbose = FALSE,
                         as.numeric(as.vector(Part4$wilcox.P.Value)), statistic="pvalue")$`pval`
  Part5=dplyr::bind_cols(Part4,fdr.P.Value=Fdr)
  Part5$log2FC=as.numeric(as.character(Part5$log2FC))
  if (sys.language){
    message(tmcn::toUTF8('\u8BBE\u7F6E|log2(FC)|\u6700\u4F4E\u503C\u6765\u7B5B\u9009\u5DEE\u5F02\u57FA\u56E0'))
  }else{
    message("Steps 11/12: Assignment the Filter of absolute logs(FC)")
  }

  diamessage=paste0("Absolute Value [log2(Folder Change)] Filter >=\n",
                    min(abs(Part5$log2FC),na.rm = TRUE),
                    "~",
                    max(abs(Part5$log2FC),na.rm = TRUE))
  res<-svDialogs::dlg_input(message = diamessage)$res
  if (!length(res)){
    stop("You cancelled the choice\n")
  }else{
    absLog2FC=as.numeric(res)
    cat("           absLog2FC >= ",absLog2FC,"\n")
  }
  Part6=Part5[abs(Part5$log2FC)>absLog2FC, ]#add param
  if (sys.language){
    message(tmcn::toUTF8('\u6B65\u9AA412/12: \u8BBE\u7F6Efdr p\u503C\u7684\u6700\u5927\u503C\u6765\u7B5B\u9009\u5DEE\u5F02\u57FA\u56E0'))
  }else{
    message("Steps 12/12: Assignment the Filter of fdr P Value")
  }
  diamessage=paste0("fdr P Value Filter <=\n",
                    min(Part6$fdr.P.Value,na.rm = TRUE),
                    "~",
                    max(Part6$fdr.P.Value,na.rm = TRUE))
  res<-svDialogs::dlg_input(message = diamessage)$res
  if (!length(res)){
    stop("You cancelled the choice\n")
  }else{
    fdr.p.filter=as.numeric(res)
    cat("           fdr.p.filter <= ",fdr.p.filter,"\n")
  }
  Part7=Part6[Part6$fdr.P.Value<fdr.p.filter, ]#add param
  result=c(result,diffGene=list(Part7))
  cat("\n","return two list: ExpressionMatrix and diffGene","\n")
  }else if(length(SEtable)==1){
    if (sys.language){
      tmcn::toUTF8('\u8FD4\u56DE\u4E00\u4E2Alist: ExpressionMatrix')
    }else{
      cat("\n","return one list: ExpressionMatrix","\n")
    }
  }
  return(result)
}
