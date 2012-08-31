########################################
######## Draw Copy Number Graph ########
########################################

utest <- function(name1=NULL,name2=NULL,nameAll=NULL,normalize=0,normalizedFile=NULL){

	############################ read data ##############################
	data1 <- read.table(name1,sep='\t')
	data2 <- read.table(name2,sep='\t')
		
	if(normalize == 1){
		a <- median(data2$V3)
		data1_ <- 2*data1$V3/a
		data2_ <- 2*data2$V3/a
	}
	else{	
		data1_ <- data1$V1
		data2_ <- data2$V1
	}
	
	u <- sprintf("%e",wilcox.test(data1_,data2_,alternative=c("two.sided"))$p.value)
	write.table(u,file=nameAll,append=FALSE,row.names=FALSE,col.names=FALSE,eol="")
	if(normalize == 1){
		write.table(data1_,file=normalizedFile,append=FALSE,row.names=FALSE,col.names=FALSE)
	}
}
