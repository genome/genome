library(ggplot2)

working_dir = "~/Desktop/"
setwd(working_dir)
dir()
data=read.table("PNC2_Met1_vs_Met2.txt", sep="\t", header=TRUE, as.is=c(1:2))

data[,"Met1"] = log2(data[,"Pre.treat"]+0.1) 
data[,"Met2"] = log2(data[,"Recurrence"]+0.1) 

print({
	p <- ggplot(data, aes(Met1, Met2))
	p + geom_point(aes(colour = Met1)) + xlab("Pre-treatment Met. (log2(FPKM+0.1))") + ylab("Recurrence Met. (log2(FPKM+0.1))") + opts(title = "Pre-treatment vs. recurrence RNA-seq gene expression values")
})

x=data[,"Pre.treat"]
y=data[,"Recurrence"]
range(x)
range(y)
(cor(x=x, y=y))^2

