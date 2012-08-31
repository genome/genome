########################################
######## Draw Copy Number Graph ########
########################################

readcount <- function(name=NULL){

	############################ read data ##############################
	name <- read.table(name,sep='\t');
	inFile <- as.character(name$V1[1])
	outFileL <- as.character(name$V2[1])
	outFileR <- as.character(name$V3[1])
	inFile_n <- as.character(name$V4[1])
	outFileL_n <- as.character(name$V5[1])
	outFileR_n <- as.character(name$V6[1])
	
	# array data
	arrayFile <- as.character(name$V5[4])
	isArray <- as.character(name$V6[4])
	
	# snp
	isSnp <- as.character(name$V1[5])
	
	# isFixYAxisLimit
	isFixYAxisLimit <- as.character(name$V4[5])

	# control the flow
	tumor <- 1
	normal <- 1
	if(is.na(inFile))
	  tumor <- 0
	if(is.na(inFile_n))
	  normal <- 0


	# tumor
	if(tumor == 1){
		a <- try(read.table(inFile, sep="\t", quote="", comment=""))
		if(class(a)=='try-error'){
			a <- NULL
		}
		a_ <- a$V3
		a_axis <- a$V2
		
		aL <- try(read.table(outFileL, sep="\t", quote="", comment=""))
		if(class(aL)=='try-error'){
			aL <- NULL
		}		
		aR <- try(read.table(outFileR, sep="\t", quote="", comment=""))
		if(class(aR)=='try-error'){
			aR <- NULL
		}
		aNeighbor <- rbind(aL, aR)
		aNeighbor_ <- aNeighbor$V3
		aNeighbor_axis <- aNeighbor$V2
	}

	# normal
	if(normal == 1){
		aN <- try(read.table(inFile_n, sep="\t", quote="", comment=""))
		if(class(aN)=='try-error'){
			aN <- NULL
		}		
		aN_ <- aN$V3
		a_axisN <- aN$V2

		aLN <- try(read.table(outFileL_n, sep="\t", quote="", comment=""))
		if(class(aLN)=='try-error'){
			aLN <- NULL
		}		
		aRN <- try(read.table(outFileR_n, sep="\t", quote="", comment=""))
		if(class(aRN)=='try-error'){
			aRN <- NULL
		}		
		aNeighborN <- rbind(aLN, aRN)
		aNeighborN_ <- aNeighborN$V3
		aNeighbor_axisN <- aNeighborN$V2
	}
	
	# array
	a_axisR <- as.numeric()
	aR_ <- as.numeric()
	aNeighbor_axisR <- as.numeric()
	aNeighborR_ <- as.numeric()
	start <- as.numeric(as.character(name$V1[3]))
	end <- as.numeric(as.character(name$V2[3]))
	neighbor1_left <- as.numeric(as.character(name$V3[3]))
	neighbor1_right <- as.numeric(as.character(name$V4[3]))
	neighbor2_left <- as.numeric(as.character(name$V5[3]))
	neighbor2_right <- as.numeric(as.character(name$V6[3]))
	
	if(isArray == 1){
		load(arrayFile)
		array_pos1 <- z$pos1
		
		mychr=as.character(name$V3[2]);
		if (as.character(name$V3[2])=="X") {mychr=23};
		if (as.character(name$V3[2])=="Y") {mychr=24};
		
		a_axisR <- z[z[,1]==mychr & z[,2]> start & z[,2]<end,]$pos1;
		aR_=2*2^z[z[,1]==mychr & z[,2]> start & z[,2]<end,]$cn; 	
		
		aNeighbor_axisR <- z[z[,1]==mychr & (z[,2] > neighbor1_left & z[,2] < neighbor1_right | z[,2] > neighbor2_left & z[,2] < neighbor2_right),]$pos1;
		aNeighborR_=2*2^z[z[,1]==mychr & (z[,2] > neighbor1_left & z[,2] < neighbor1_right | z[,2] > neighbor2_left & z[,2] < neighbor2_right),]$cn; 	
			
	}


	# read the file to get the picture name, title and axis
	data_name <- as.character(name$V1[2])
	picName <- as.character(name$V2[2])
	chromosome <- paste("chr", name$V3[2], sep="")
	
	# read the plotting options
	isTitle <- as.character(name$V4[2])
	isSubTitle <-as.character(name$V5[2])
	isAnnotation <- as.character(name$V6[2])

	# prepare for the normalization even if no annotation
	seg_ <- as.numeric();
	
if(isAnnotation == 1){
	# number of the annotations
	num <- 4

	################################ annotation #############################
	############### segmental Duplication
	tmp <- name$V1[4]
	Read_Annotation(tmp) -> seg

	if(class(seg)!='try-error' && length(seg$Start)>=1 ){
		for(i in 1:length(seg$Start)){
			seg_ = c(seg_, seg$Start[i]:seg$End[i]);
		}
	}

	############### repeat mask
	tmp <- name$V2[4]
	Read_Annotation(tmp) -> rep

	############### dgv
	tmp <- name$V3[4]
	Read_Annotation(tmp) -> dgv

	############### gene
	tmp <- name$V4[4]
	Read_Annotation(tmp) -> gene
}

	################################ snp ####################################
	if(isSnp == 1){
		if(tumor == 1){
			tmp <- name$V2[5]
			Read_Snp(tmp) -> snp_tumor
		}
		if(normal == 1){
			tmp <- name$V3[5]
			Read_Snp(tmp) -> snp_normal
		}
	}
	############################# normalize data by median and segmental duplication #############
	# tumor:
	if(tumor == 1){
		normalization(aNeighbor_axis, aNeighbor_, seg_) -> aNeighbor_median_old
		a_ <- 2*a_/aNeighbor_median_old
		aNeighbor_ <- 2*aNeighbor_/aNeighbor_median_old
		a_median <- median(a_)
		aNeighbor_median <- 2
	}

	# normal:
	if(normal == 1){
		normalization(aNeighbor_axisN, aNeighborN_, seg_) -> aNeighborN_median_old
		aN_ <- 2*aN_/aNeighborN_median_old
		aNeighborN_ <- 2*aNeighborN_/aNeighborN_median_old
		aN_median <- median(aN_)
		aNeighborN_median <- 2
	}
	
	# array:
	if(isArray == 1){
		normalization(aNeighbor_axisR, aNeighborR_, seg_) -> aNeighborR_median_old
		aR_ <- 2*aR_/aNeighborR_median_old
		aNeighborR_ <- 2*aNeighborR_/aNeighborR_median_old
		aR_median <- median(aR_)
		aNeighborR_median <- 2
	}
	
	############################## not to consider those in annotation segmental duplication #######
	#tumor
	if(tumor == 1){
		aNeighbor_NoSeg <- aNeighbor_[! aNeighbor_axis %in% seg_]
		aNeighbor_NoSeg_axis <- aNeighbor_axis[! aNeighbor_axis %in% seg_]
		aNeighbor_Seg <- aNeighbor_[aNeighbor_axis %in% seg_]
		aNeighbor_Seg_axis <- aNeighbor_axis[aNeighbor_axis %in% seg_]
	}
	#normal
	if(normal == 1){
		aNeighborN_NoSeg <- aNeighborN_[! aNeighbor_axisN %in% seg_]
		aNeighborN_NoSeg_axis <- aNeighbor_axisN[! aNeighbor_axisN %in% seg_]
		aNeighborN_Seg <- aNeighborN_[aNeighbor_axisN %in% seg_]
		aNeighborN_Seg_axis <- aNeighbor_axisN[aNeighbor_axisN %in% seg_]
	}
	#array
	if(isArray == 1){
		aNeighborR_NoSeg <- aNeighborR_[! aNeighbor_axisR %in% seg_]
		aNeighborR_NoSeg_axis <- aNeighbor_axisR[! aNeighbor_axisR %in% seg_]
		aNeighborR_Seg <- aNeighborR_[aNeighbor_axisR %in% seg_]
		aNeighborR_Seg_axis <- aNeighbor_axisR[aNeighbor_axisR %in% seg_]
	}

	############################### printing start here ################################################
	cex_ <- 0.6
	pch_ <- 16
	lwd_ <- 3
	cex_snp <- 0.6
	y_snp_all <- 0.5
	y_snp <- -0.3
	pch_snp <- 8

	png(picName)

	if((tumor == 1 && normal == 1 || tumor == 1 && isArray == 1) && (isAnnotation == 1 || isSnp == 1)){
	    nf <- layout(matrix(c(1:4),2,2,byrow=TRUE), c(2,2), c(3,1), TRUE)
	}
	if((tumor == 1 && normal == 1 || tumor == 1 && isArray == 1) && (isAnnotation != 1 && isSnp != 1)){
  		nf <- layout(matrix(c(1:2),nrow=1,ncol=2,byrow=TRUE), c(2,2), c(4), TRUE)
  	}
  	if(!(tumor == 1 && normal == 1 || tumor == 1 && isArray == 1) && (isAnnotation == 1 || isSnp == 1)){
  		nf <- layout(matrix(c(1:2),nrow=2,ncol=1,byrow=TRUE), c(4), c(3,1), TRUE)
  	}
  	if(!(tumor == 1 && normal == 1 || tumor == 1 && isArray == 1) && (isAnnotation != 1 && isSnp != 1)){
  		nf <- layout(matrix(c(1:1)), c(4), TRUE)
  	}
	layout.show(nf)


	if(tumor == 1){
		a_axis_all <- c(a_axis, aNeighbor_axis)
	} else if(normal == 1){
		a_axis_all <- c(a_axisN, aNeighbor_axisN)
	} else{
		a_axis_all <- c(a_axisR, aNeighbor_axisR)
	}


	if(tumor == 1){
		if(normal == 1 || isArray == 1)
			par(mar=c(3,5,5,0))
		else
			par(mar=c(3,5,5,3))
			
		x_limit <- c(min(c(a_axis,aNeighbor_axis)),max(c(a_axis,aNeighbor_axis)))
		if(isFixYAxisLimit == 1)
			y_limit <- c(0-y_snp_all,4)
		else{
			if(normal==1){		
			y_limit <- c(max(0, min(mean(c(a_,aNeighbor_))-3*sd(c(a_,aNeighbor_)),mean(c(a_,aNeighborN_))-3*sd(c(a_,aNeighborN_))))-y_snp_all,max(mean(c(a_,aNeighbor_))+3*sd(c(a_,aNeighbor_)),mean(c(a_,aNeighborN_))+3*sd(c(a_,aNeighborN_))))
			}
			else{
			y_limit <- c(max(0, min(mean(c(a_,aNeighbor_))-3*sd(c(a_,aNeighbor_))))-y_snp_all,mean(c(a_,aNeighbor_))+3*sd(c(a_,aNeighbor_)))
			}
		}
		Printing_Graph(a_, a_axis, "red", a_median, "green", aNeighbor_NoSeg, aNeighbor_NoSeg_axis, "blue", aNeighbor_Seg, aNeighbor_Seg_axis, "blue", aNeighbor_median, "black", "Tumor", pch_, cex_, lwd_, isSubTitle, x_limit, y_limit)
		################# snp
		if(isSnp == 1){		
			y_snp_dots <- y_snp - 0.15
			Draw_Snp(snp_tumor, a_axis_all, y_snp, "grey", "purple", 1, cex_snp, pch_snp, "Het Snp", y_snp_dots, 0.8, 0)
		}
	}

	if(normal == 1){
		if(tumor ==1 )
			par(mar=c(3,2,5,3))
		else
			par(mar=c(3,5,5,3))
			
		x_limit <- c(min(c(a_axisN,aNeighbor_axisN)),max(c(a_axisN,aNeighbor_axisN)))
		if(isFixYAxisLimit == 1)
			y_limit <- c(0-y_snp_all,4)
		else{
			if(normal==1){
			y_limit <- c(max(0, min(mean(c(a_,aNeighbor_))-3*sd(c(a_,aNeighbor_)),mean(c(a_,aNeighborN_))-3*sd(c(a_,aNeighborN_))))-y_snp_all,max(mean(c(a_,aNeighbor_))+3*sd(c(a_,aNeighbor_)),mean(c(a_,aNeighborN_))+3*sd(c(a_,aNeighborN_))))
			}
			else{
			y_limit <- c(max(0, min(mean(c(a_,aNeighborN_))-3*sd(c(a_,aNeighborN_))))-y_snp_all,mean(c(a_,aNeighborN_))+3*sd(c(a_,aNeighborN_)))
			}
		}
		Printing_Graph(aN_, a_axisN, "red", aN_median, "green", aNeighborN_NoSeg, aNeighborN_NoSeg_axis, "blue", aNeighborN_Seg, aNeighborN_Seg_axis, "blue", aNeighborN_median, "black", "Normal", pch_, cex_, lwd_, isSubTitle, x_limit, y_limit)
		################# snp
		if(isSnp == 1){		
			y_snp_dots <- y_snp - 0.15	
			Draw_Snp(snp_normal, a_axis_all, y_snp, "grey", "purple", 1, cex_snp, pch_snp, "Het Snp", y_snp_dots, 0.8, 0)
		}
	}
	
	if(isArray == 1){
		if(tumor ==1 )
			par(mar=c(3,2,5,3))
		else
			par(mar=c(3,5,5,3))

		x_limit <- c(min(c(a_axisR,aNeighbor_axisR)),max(c(a_axisR,aNeighbor_axisR)))
		if(isFixYAxisLimit == 1)
			y_limit <- c(0,4)
		else
			y_limit <- c(max(0, min(mean(c(a_,aNeighborR_))-3*sd(c(a_,aNeighborR_)))),mean(c(a_,aNeighborR_))+3*sd(c(a_,aNeighborR_)))
		Printing_Graph(aR_, a_axisR, "red", aR_median, "green", aNeighborR_NoSeg, aNeighborR_NoSeg_axis, "blue", aNeighborR_Seg, aNeighborR_Seg_axis, "blue", aNeighborR_median, "black", "Array", pch_, cex_, lwd_, isSubTitle, x_limit, y_limit)

	}	

	lwd_ = 7

#if(isSnp == 1 && isAnnotation != 1){
#		num <- 1
#		################# snp
#	if(tumor == 1){
#		if(normal == 1 || isArray == 1)
#			par(mar=c(0,5,0,0))
#		else
#			par(mar=c(0,5,0,3))	
#		
#		# tumor
#		Draw_Snp(snp_tumor, aNeighbor_axis, num, "grey", "red", 1, cex_snp, pch_, "Snp", num-num/10, 0.8, 1)
#	}
#	if(normal == 1){
#		if(tumor == 1)
#			par(mar=c(0,2,0,3))
#		else
#			par(mar=c(0,5,0,3))
#				
#		# normal		
#		Draw_Snp(snp_normal, aNeighbor_axis, num, "grey", "red", 1, cex_snp, pch_, "Snp", num-num/10, 0.8, 1)
#	}
#}
		
if(isAnnotation == 1){
	if(tumor == 1){
		if(normal == 1 || isArray == 1)
			par(mar=c(0,5,0,0))
		else
			par(mar=c(0,5,0,3))
			
		##################### segmental Duplication
		Draw_Annotation_First(seg, a_axis_all, num, "grey", "purple", 1, lwd_, "Segmental Duplication", num-0.4, 0.8)
		
		################# repeat Mask
		Draw_Annotation(rep, a_axis_all, num-1, "grey", "green", 1, lwd_, "Repeat Mask", num-1.4, 0.8)

		################# gene
		Draw_Annotation(gene, a_axis_all, num-2, "grey", "yellow", 1, lwd_, "Gene", num-2.4, 0.8)

		################# dgv
		Draw_Annotation(dgv, a_axis_all, num-3, "grey", "black", 1, lwd_, "Database of Genomic Variants", num-3.4, 0.8)

		################# snp
		#if(isSnp == 1){		
		#	Draw_Snp(snp_tumor, a_axis_all, num-4, "grey", "red", 1, cex_snp, pch_, "Snp", num-4.4, 0.8, 0)
		#}
	}

	if(normal == 1){
		if(tumor == 1)
			par(mar=c(0,2,0,3))
		else
			par(mar=c(0,5,0,3))
		##################### segmental Duplication
		Draw_Annotation_First(seg, a_axis_all, num, "grey", "purple", 1, lwd_, "Segmental Duplication", num-0.4, 0.8)
		
		################# repeat Mask
		Draw_Annotation(rep, a_axis_all, num-1, "grey", "green", 1, lwd_, "Repeat Mask", num-1.4, 0.8)
		
		################# gene
		Draw_Annotation(gene, a_axis_all, num-2, "grey", "yellow", 1, lwd_, "Gene", num-2.4, 0.8)

		################# dgv
		Draw_Annotation(dgv, a_axis_all, num-3, "grey", "black", 1, lwd_, "Database of Genomic Variants", num-3.4, 0.8)
		
		################# snp
		#if(isSnp == 1){
		#	Draw_Snp(snp_normal, a_axis_all, num-4, "grey", "red", 1, cex_snp, pch_, "Snp", num-4.4, 0.8, 0)
		#}
	}
	
	if(isArray == 1){
		if(tumor == 1)
			par(mar=c(0,2,0,3))
		else
			par(mar=c(0,5,0,3))
		##################### segmental Duplication
		Draw_Annotation_First(seg, a_axis_all, num, "grey", "purple", 1, lwd_, "Segmental Duplication", num-0.4, 0.8)
		
		################# repeat Mask
    	Draw_Annotation(rep, a_axis_all, num-1, "grey", "green", 1, lwd_, "Repeat Mask", num-1.4, 0.8)
		
		################# gene
		Draw_Annotation(gene, a_axis_all, num-2, "grey", "yellow", 1, lwd_, "Gene", num-2.4, 0.8)
#		Draw_Annotation_First(gene, a_axis_all, num, "grey", "yellow", 1, lwd_, "Gene", num-0.5, 1.2)
		################# dgv
		Draw_Annotation(dgv, a_axis_all, num-3, "grey", "black", 1, lwd_, "Database of Genomic Variants", num-3.4, 0.8)
	}
}
	# main title
	if(is.na(data_name))
		main_title <- chromosome
	else
  		main_title <- paste(data_name,chromosome,sep=" ")
  	if(isTitle == 1){
		mtext(main_title, outer=TRUE, line=-2, cex=1.5)
	}

	dev.off()
}

Read_Annotation = function(var)
{
	file_name <- as.character(var)
	data <- try(read.table(file_name, sep='\t', header=TRUE))
	#if(class(data)=='try-error'){
	#	data <- NULL
	#} else if(gregexpr(",",data$V1[1])>0){
	#	data <- try(read.table(file_name, sep=","));
	#}
	data
}

normalization = function(x, y, annot)
{
	y_ <- y[! x %in% annot]
	mean_y_ <- mean(y_)
# no exclusion of those out of the range	
#	sd_y_ <- sd(y_)
#	new_y <- y_[y_ > mean_y_ - sd_y_ & y_ < mean_y_ + sd_y_ & y_ != 0]
	new_y <- y_[y_!=0]
	median_new_y <- median(new_y)
	median_new_y
}

Printing_Graph = function(y,x,col_data,med,col_med,y1,x1,col_data1,y2,x2,col_data2,med_,col_med_,title_,pch_,cex_,lwd_,isSubTitle,x_limit,y_limit)
{
	plot(y1 ~ x1, col = col_data1, xlim = x_limit, ylim = y_limit, xlab = "Base", ylab = "Copy Number", pch = pch_, cex = cex_)
	points(y~x, col=col_data, pch=pch_, cex=cex_)
	points(y2~x2, col=col_data2, pch=pch_, cex=cex_)
	
	# lines the median
	segments(min(x1,x2), med_, min(x), med_, col=col_med_, lwd=lwd_)
	segments(max(x), med_, max(x1,x2), med_, col=col_med_, lwd=lwd_)
	segments(min(x), med, max(x), med, col=col_med, lwd=lwd_)
	
	if(isSubTitle==1) {
		title(title_, cex.main=1.1, line = 1)
	}
}

Draw_Annotation = function(data,x,y,col1,col2,lwd1,lwd2,name,height,cex_)
{
	segments(min(x),y,max(x),y,col=col1,lwd=lwd1)
	if(length(data$Start)>=1){
		for(i in 1:length(data$Start)){
			segments(data$Start[i], y, data$End[i], y, col=col2, lwd=lwd2)
		}
	}
	text((min(x)+max(x))/2, height, name, cex=cex_)
}

Draw_Annotation_First = function(data,x,y,col1,col2,lwd1,lwd2,name,height,cex_)
{
	test_coords <- xy.coords(min(x):max(x),y,recycle=TRUE)
	plot(test_coords$x,test_coords$y,col=col1,xlim=c(min(x),max(x)),ylim=c(0,y),axes=FALSE,xlab="",ylab="",type="l",lwd=lwd1)
#	plot(test_coords$x,test_coords$y,col=col1,xlim=c(min(x),max(x)),ylim=c(-1,y),axes=FALSE,xlab="",ylab="",type="l",lwd=lwd1)
	if(length(data$Start)>=1){
		for(i in 1:length(data$Start)){
			segments(data$Start[i], y, data$End[i], y, col=col2, lwd=lwd2)
		}
	}
	text((min(x)+max(x))/2, height, name, cex=cex_)
}

Draw_Snp = function(data, x, y, col1, col2, lwd1, cex2, pch_, name, height, cex_, newplot)
{
	data_y <- seq(length=length(data$V1),from=y,by=0)
	if(newplot!=1){
		segments(min(x),y,max(x),y,col=col1,lwd=lwd1)
		if(length(data$V1)>0){
		points(data_y~data$V1, col=col2, cex=cex2, pch=pch_)
		}
	}
	else{
		test_coords <- xy.coords(min(x):max(x),y,recycle=TRUE)
		plot(test_coords$x,test_coords$y,col=col1,xlim=c(min(x),max(x)),ylim=c(0,y),axes=FALSE,xlab="",ylab="",type="l",lwd=lwd1)
		if(length(data$V1)>0){		
		points(data_y~data$V1, col=col2,cex=cex2, pch=pch_)
		}
	}
	text((min(x)+max(x))/2, height, name, cex=cex_)
	
}
	
Read_Snp = function(tmp)
{
	a <- try(read.table(as.character(tmp)))
	if(class(a)=='try-error'){
		a <- NULL
	}
	a
}
