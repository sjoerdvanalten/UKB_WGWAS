args <- commandArgs(trailingOnly = TRUE)


Input<-args[1]
Output<-args[2]

data<-read.table(Input,sep=",",header=TRUE)
write.table(data,file=Output,sep="\t",quote=FALSE,row.names=FALSE)
