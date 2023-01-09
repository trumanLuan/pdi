readPWM<-function(pwm.file,tmp.name){
pwm.table=readLines(pwm.file)
name.ind=grep(">",pwm.table)

pwm.list=c()
for(tmp.ind in name.ind){
tmp.name=sub(">","",pwm.table[tmp.ind])
tmp.pwm=read.table(pwm.file,skip=tmp.ind,nrow=4,sep=" ",row.names=1,header=F,as.is=TRUE)
tmp.pwm.list=list(tmp.pwm)
names(tmp.pwm.list)=tmp.name
pwm.list=c(pwm.list,tmp.pwm.list)
} # end tmp.ind

return(pwm.list)
} # end of function
