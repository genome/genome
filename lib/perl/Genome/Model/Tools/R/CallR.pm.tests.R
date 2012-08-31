########################################
######## Draw Copy Number Graph ########
########################################

normalize <- function(name1=NULL,nameL=NULL,nameR=NULL,normalizedFile=NULL){
	############################ read data ##############################
	data1 <- read.table(name1,sep='\t')
	dataL <- read.table(nameL,sep='\t')
	dataR <- read.table(nameR,sep='\t')	
	data2 <- c(dataL$V3, dataR$V3)
	
	a <- median(data2[data2!=0])
	data1_ <- 2*data1$V3/a
	
	write.table(data1_,file=normalizedFile,append=FALSE,row.names=FALSE,col.names=FALSE)
}

ttest <- function(name1=NULL,name_outL=NULL,name_outR=NULL,nameAll=NULL){

	############################ read data ##############################
	data1 <- read.table(name1,sep='\t')
	dataL <- read.table(name_outL,sep='\t')
	dataR <- read.table(name_outR,sep='\t')
		
	data1_ <- data1$V3
	data2_ <- c(dataL$V3,dataR$V3)
	
	data1_ <- data1_[data1_!=0]
	data2_ <- data2_[data2_!=0]

	u <- sprintf("%e",t.test(data1_,data2_)$p.value)
#	u <- t.test(data1_,data2_)$p.value
	write.table(u,file=nameAll,append=FALSE,row.names=FALSE,col.names=FALSE,eol="",quote=FALSE)
}

ttest_ROI <- function(name1=NULL,name2=NULL,nameAll=NULL){

	############################ read data ##############################
	data1 <- read.table(name1,sep='\t')
	data2 <- read.table(name2,sep='\t')
		
	data1_ <- data1$V1
	data2_ <- data2$V1
	
	data1_ <- data1_[data1_!=0]
	data2_ <- data2_[data2_!=0]
	
	u <- sprintf("%e",t.test(data1_,data2_)$p.value)
	write.table(u,file=nameAll,append=FALSE,row.names=FALSE,col.names=FALSE,eol="",quote=FALSE)
}

mean_ROI <- function(name=NULL,nameAll=NULL){

	############################ read data ##############################
	data <- read.table(name,sep='\t')
		
	data_ <- data$V1
	
	data_ <- data_[data_!=0]

	########################### compute mean ############################
	mean <- sprintf("%e",mean(data_))
	
	write.table(mean,file=nameAll,append=FALSE,row.names=FALSE,col.names=FALSE,eol="",quote=FALSE)	
}	
