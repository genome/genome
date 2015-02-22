#takes the following arguments:
#
# arg 1 = filename
# arg 2 = position of column containing tumor reference-supporting reads
# arg 3 = position of column containing tumor variant-supporting reads
# arg 4 = output filename
args <- commandArgs(trailingOnly = TRUE)

b = read.table(args[1],skip=1,sep="\t")

refcol = as.numeric(args[2])
varcol = as.numeric(args[3])

pval=c();
type=c();
for(i in 1:length(b$V3)){
  reads=(b[i,refcol]+b[i,varcol]);
  het=round(reads/2)
  hetmat = matrix(c(b[i,refcol], het, b[i,varcol], het),nrow=2)
  hommat = matrix(c(b[i,refcol], 0, b[i,varcol], reads),nrow=2)

  p.het=fisher.test(hetmat)$p.value
  p.hom=fisher.test(hommat)$p.value

  if(p.het > p.hom){
    pval=c(pval,p.het)
    type=c(type,"het")
  } else {
    pval=c(pval,p.hom)
    type=c(type,"hom")
  }
}

b = cbind(b,type)
b = cbind(b,pval)
b = cbind(b,p.adjust(b$pval,method="BY"))

write.table(b,args[4],row.names=F,col.names=F,quote=F,sep="\t")
