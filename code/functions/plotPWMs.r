plotPWMs<-function(pwm.list,nrow=8,ncol=6,file.name="./tmpPWMs.pdf",plotsave=TRUE){
library(motifStack)

na.pwm.names=names(pwm.list)[sapply(pwm.list,function(x)any(is.na(x)))]
if(length(na.pwm.names)>0){
stop(paste(paste(na.pwm.names,collapse=",")," has NA values",sep=""))}

pfms<-lapply(pwm.list,pcm2pfm)
pfms<-lapply(names(pfms), function(.ele, pfms){new("pfm",mat=pfms[[.ele]], name=.ele)},pfms)
# pfms.aln<-DNAmotifAlignment(pfms)
if(plotsave==TRUE){
pdf(file.name)
plotMotifLogoStack(pfms,par(mfcol=c(nrow,ncol)))
dev.off()} else{
plotMotifLogoStack(pfms,par(mfcol=c(nrow,ncol)))}
}

