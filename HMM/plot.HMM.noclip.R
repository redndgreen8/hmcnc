args = commandArgs(trailingOnly=TRUE)

data<-read.table(args[1])



x<-data$V3-50

y<-data$V5
#state<-unlist(strsplit(CN,".",fixed=TRUE))
#CNstate<-state$V1
#CLstate<-state$V2
#
#for (i in 1:length(CNstate)){
#    y=CNstate[i]
#    if (CLstate[i]==2)
#        y=y+0.5
#}
#mean<-as.integer(args[4])
mean_file<-scan(args[4])
mean<-(as.integer(mean_file[1])/2)

contig<-data$V1
bounds<-seq(1,length(y),as.integer(args[3]))
obs<-data$V4
pdf(paste(args[2],"noclip.pdf",sep="."),paper="a4r",width=10)
#par(mfrow=c(2,1))
for (i in 1:(length(bounds)-1) ){

    plot( x[bounds[i]:bounds[i+1]] ,  obs[bounds[i]:bounds[i+1]]/mean  ,  ylim=c(0 ,min( 32, (max(obs[bounds[i]:bounds[i+1]]))/mean  ) )  , xlab="reference coordinate"  , ylab="Coverage/CN"  , col="blue"  ,pch=".", cex=1.5, sub=paste(contig[bounds[i]]),type="l",yaxt='n' )
    title(main=paste(args[2]))
    axis(2 , at = seq(0,max(y[bounds[i]:bounds[i+1]]),1) , las=2)
    axis(4 , at = seq(0,max(y[bounds[i]:bounds[i+1]]),1),labels=seq(0, max(y[bounds[i]:bounds[i+1]])*mean,mean ) , las=2)
    grid( nx = NULL, ny = NULL, col = "lightgrey" , lty="longdash")

    lines( x[bounds[i]:bounds[i+1]] , y[bounds[i]:bounds[i+1]] , ylab="HMM CN state"  , col="red"  ,pch=".", cex=2.5)


}
dev.off()
