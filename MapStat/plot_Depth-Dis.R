#!/bin/Rscript
args=commandArgs(T)
if(length(args) != 2) stop( "\tHelp?\n\n\tRscript plot_Depth-Dis.R <Sample.Depth.freq> <sample_name>\n")

Depth_freq_file=args[1]
sample_name=args[2]
#Depth_freq_file="N3069D.Panel.Depth.freq"
#sample_name="test"
data=read.table(Depth_freq_file,sep="\t",header=T)


out_DP=paste(sample_name,".Depth-Distribution.pdf",sep="")
pdf(out_DP,width=8,height=10)
par(mfrow=c(2,1))
barplot(data[,3],names.arg = data[,1],xlim=c(0,500),col="#7777FF",border=NA,main=paste(sample_name," depth distribution",sep=""),ylab="Percent of bases (%)",xlab="Depth")
plot(data[,1],100-data[,4],type="l",xlim=c(0,300),col="red",lwd=2,main=paste(sample_name," depth cumulative distribution",sep=""),ylab="Percent of bases (%)",xlab="Depth",frame.plot = FALSE)
dev.off()
