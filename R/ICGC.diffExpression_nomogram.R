#' To plot nomogram for cox regression
#'
#' @return a picturn
#' @export
#'
#' @examples ICGC.diffExpression_nomogram()
ICGC.diffExpression_nomogram <- function(){
  message("********Please run the code bellow********")
  message("********Please run the code bellow********")
  cat("\ndiffExpr.multiVarCox.data: dataset gained from command ICGC.diffExpression.multiVarCox\n")
  message("********Please run the code bellow********")
  cat("\n")
cat('scoredata=diffExpr.multiVarCox.data$coxphdata
scoredata2=scoredata[order(scoredata$RiskScore),]
rownames(scoredata2)=scoredata2$icgc_donor_id
cphdata=cbind(donor_vital_status=scoredata2$donor_vital_status,
              donor_survival_time=scoredata2$donor_survival_time,
              scoredata2[,rownames(diffExpr.multiVarCox.data$coxsumHR)])
suppressMessages(library(rms))
dd<-rms::datadist(cphdata)
options(datadist="dd")
f <- cph(survival::Surv(donor_survival_time,
                        donor_vital_status) ~ .,
         data=cphdata,surv = T)
surv <- Survival(f)
nom <- nomogram(f, maxscale=100,
                fun=list(function(x) surv(365, x),
                         function(x) surv(1095, x),
                         function(x) surv(1825, x)),
                lp=F, funlabel=c("1-year survival",
                                 "3-year survival",
                                 "5-year survival")
)
plot(x = nom)')
cat("\n")
cat("\n")
message("********Please run the code above********")
}