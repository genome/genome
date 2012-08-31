e=commandArgs();
png(paste(e[3],".png",sep=''),width=1280,height=960)
par(mfrow=c(4,6))
x=read.table(e[3],comment.char='#',header=TRUE);
for (i in c(1:22,'X')){
  y=subset(x,CHR==i);
  plot(y$POS,y$CopyNumber,main=paste("chr.",i),xlab="mb",ylab="cn",type="l",ylim=c(0,4) );
}
dev.off()
