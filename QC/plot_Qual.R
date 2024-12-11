#!/bin/Rscript
args=commandArgs(T)
if(length(args) != 2) stop( "\tParameters count not correct! \n\n\tRscript plot_Qual.R <Qual.matrix> <Sample_name> \n\n")

file=args[1]
sample_name=args[2]

library(ggplot2)
require(gridExtra)

data=read.table(file,check.names = F,header=T,sep="\t")

pdf(paste(sample_name,".Qual.pdf",sep=""),width=14,height=8)

p1=ggplot(data)+geom_point(aes(x=Pos,y=R1_ave_qual,colour = R1_ave_qual),size=0.8)+ylim(30,38) + theme(legend.position='none') + ylab("Quality score") + xlab("Bases Position") + geom_smooth(aes(x=Pos,y=R1_ave_qual),color="red",fill="grey",size=0.8) + labs(title=paste(sample_name," Read1 Quality Distribution",sep="")) + theme(plot.title = element_text(hjust = 0.5))
p2=ggplot(data)+geom_point(aes(x=Pos,y=R2_ave_qual,colour = R2_ave_qual),size=0.8)+ylim(30,38) + theme(legend.position='none') + xlab("Bases Position") + theme(axis.ticks = element_blank(), axis.text.y  = element_blank()) + ylab(label = NULL) +geom_smooth(aes(x=Pos,y=R2_ave_qual),color="red",fill="grey",size=0.8) + labs(title=paste(sample_name," Read2 Quality Distribution",sep="")) + theme(plot.title = element_text(hjust = 0.5))
grid.arrange(p1,p2,ncol=2)

dev.off()

