# TODO: Add comment
# 
# Author: ygindin
###############################################################################


gatherDataFromVariantFile <- function(inputs, geneNameIndex, variantIndex)
{
	n <- nrow(inputs)
	variantData <- data.frame();
	for (i in 1:n)
	{
		filePath <- inputs[i, 3]
		data <- read.delim(filePath)
		data <- data[, c(geneNameIndex,variantIndex)]
		names(data)[2] <- 'default_gene_name'
		data$sample <- inputs[i, 1]
		variantData <- rbind(variantData, data)
		
#		filePath <- inputs[i, 3]
#		data <- read.delim(filePath)
#		print(data)
#		data <- data[, c(1,2,3,geneNameIndex,variantIndex)]
#		data$sample <- inputs[i, 1]
#		variantData <- rbind(variantData, data)
	}	
	#print(variantData)
	return(variantData)
}


processCopyNumberData <- function(inputs, genesOfInterest, completeGeneModelFile, javaLib)
{
	# Create one file with all of the CNs
	n <- nrow(inputs)
	tempCnFile <- tempfile("tempCNfile", tempdir())
	data <- data.frame()
	for (i in 1:n)
	{
		filePath <- inputs[i, 2]
		#print(filePath)
		sampleData <- read.delim(filePath, header=F)
		sampleData <- sampleData[, -4]
		sampleData$sampleName <- inputs[i, 1]
		data <- rbind(data, sampleData)
		
	}
	#print(data)
	write.table(data, tempCnFile, append=F, quote=F, sep="\t",
			row.names=F, col.names=F)
	
	# Select only those genes from the gene model that are of interest
	allGenes <- read.delim(completeGeneModelFile, header=F)
	names(allGenes) <- c('chr', 'start', 'end', 'geneName')
	allGenes <- allGenes[allGenes$geneName %in% genesOfInterest, ]
	tempGeneFile <- tempfile("tempGeneFile", tempdir())
	write.table(allGenes, tempGeneFile, quote=F, sep="\t",
			row.names=F, col.names=F)
	#print(allGenes)
	#print(file.exists(tempGeneFile, tempCnFile))
	
	# intersect the two files with Java
	tempCnOverlapFile <- tempfile("tempCnOverlapFile", tempdir())
	command <- paste("java -jar",
					javaLib,
					"OncoprintFeatureOverlapper",
					"-f ", tempGeneFile, 
					"-r ", tempCnFile, 
					"-o ", tempCnOverlapFile)
	#print(command)
	system(command)
	cnOverlaps <- read.delim(tempCnOverlapFile, header=F)
	#print(cnOverlaps)
	cnOverlaps <- cnOverlaps[, c(4, 8, 9)]
	names(cnOverlaps) <- c('gene', 'value', 'sample')
	return(cnOverlaps)
}

getXandYforGeneAndSample <- function(gene, sample, rawPlot)
{
	geneLabels <- rawPlot$panel$ranges[[1]]$y.labels
	sampleLabels <- rawPlot$panel$ranges[[1]]$x.labels
	
	y <- which(geneLabels %in% gene);
	x <- which(sampleLabels %in% sample);
	return(c(x,y))
}


guessSampleNameFromCopyNumberPath <-function(filePath)
{
	x <- gregexpr("/", filePath)
	sampleInfo <- substr(filePath, start=x[[1]][7] + 1, stop = x[[1]][8] - 1 )
	y <- gregexpr("_", sampleInfo)
	sampleName <- substr(sampleInfo, start=1, stop=y[[1]][1] - 1)
	return(sampleName)
}

#getDataFromFileColumn <- function(filePath, genes)
#{
#	fileName <- basename(filePath)
#	
#	x <- gregexpr("_", fileName)
#	sampleName <- substr(fileName, start = 1, stop = x[[1]][1] - 1)	
#	
#	data <- read.delim(filePath)
#	data <- subset(data, default_gene_name %in% genes, select=c("default_gene_name", "type"))
#	#data <- data[, c("default_gene_name", "type")]
#}
