#!/usr/bin/env Rscript
# A poor man's version of OncoPrint
# TODO: What happens if there were two or more events per gene and sample?
# Author: Yevgeniy Gindin
###############################################################################

options(stringsAsFactors = FALSE)
library(plyr)
library(ggplot2)
library(reshape)
library(doBy)
library(optparse)


#setwd(".")
option_list <- list(
		make_option(c("--genes"), type="character", default=NULL,
				dest="genes",help="Comma separated list (no spaces) of genes to analyze"),

		make_option(c("--variant-index"), type="integer", default=6,
				dest="variantIndex",help="Column number in the variant file that stores the variant type, [default %default]"),
		
		make_option(c("--gene-name-index"), type="integer", default=22,
				dest="geneNameIndex",help="Column number in the variant file that stores the gene name, [default %default]"),
		
		make_option(c("-i", "--input"), type="character", default=NULL,
				dest="input", help="Tab-delimited file with three columns (no header): sample name; path to copyCat file; path to variant file"),
				
		make_option(c("-m", "--gene-model"), type="character", default=NULL,
				dest="geneModelFile", 
				help="Path to a BED file with gene names (4th column should contain gene name)"),

		make_option(c("-o", "--out-file"), type="character", default=NULL,
				dest="outFile", 
				help="Path where the output image ought to be placed"),
		
		make_option(c("--helper-file"), type="character", default=NULL,
				dest="helperFile", 
				help="Path to the helper file"),

		make_option(c("--java-lib-file"), type="character", default=NULL,
				dest="javaLibFile", 
				help="Path to the helper Java file")
		)
opt <- parse_args(OptionParser(option_list=option_list))

genesOfInterest <- strsplit(opt$genes, ",")[[1]]
geneModelFile <- opt$geneModelFile
inputFile <- opt$input
outFile <- opt$outFile
variantIndex <- opt$variantIndex
geneNameIndex <- opt$geneNameIndex
helperFile <- opt$helperFile
javaLibFile <- opt$javaLibFile

print(paste("genesOfInterest:", genesOfInterest))
print(paste("geneModelFile:", geneModelFile))
print(paste("outFile:", outFile))
print(paste("inputFile:", inputFile))
print(paste("variantIndex:", variantIndex))
print(paste("geneNameIndex:", geneNameIndex))
print(paste("helperFile:", helperFile))
print(paste("javaLibFile:", javaLibFile))
print(paste("outFile:", outFile))

source(helperFile);
#save.image("~/temp/test.RData")

inputFile <- read.delim(inputFile, header=FALSE)

cnvData <- processCopyNumberData(inputFile, genesOfInterest, geneModelFile, javaLibFile) 
variantData <- gatherDataFromVariantFile(inputFile, variantIndex, geneNameIndex)
## Save
samples <- unique(variantData$sample); 
genes <- genesOfInterest

# Finalize SNV data
variantData <- subset(variantData, default_gene_name %in% genesOfInterest)
#print(variantData)
#names(variantData)[5] <- "gene"
# Combine rows with multiple entries
names(variantData) <- c('type', 'gene', 'sample')
variantData <- aggregate(variantData$type, list(variantData$gene, variantData$sample), paste, collapse=",")
names(variantData) <- c('gene', 'sample', 'type')


snvMatrix <- matrix(nrow=length(samples), ncol=length(genes))
rownames(snvMatrix) <- samples
colnames(snvMatrix) <- genes
snvMatrixBoolean <- matrix(nrow=length(samples), ncol=length(genes))
rownames(snvMatrixBoolean) <- samples
colnames(snvMatrixBoolean) <- genes
snvMatrixBoolean[is.na(snvMatrixBoolean)] <- 0

n <- nrow(variantData)
for (i in 1:n)
{
	dataRow <- variantData[i, ] 
	snvMatrix[dataRow$sample, dataRow$gene] <- dataRow$type
	snvMatrixBoolean[dataRow$sample, dataRow$gene] <- 1
}
snvMatrixBooleanMelted <- melt(snvMatrixBoolean)
names(snvMatrixBooleanMelted)[1:2] <- c('sample', 'gene')
snvMatrixBooleanMelted$gene <- as.character(snvMatrixBooleanMelted$gene)
snvMatrixBooleanMelted$sample <- as.character(snvMatrixBooleanMelted$sample)
cnvData$change <- ifelse(cnvData$value > 2.0, 'Gain', 'Loss')
cnvData <- cnvData[, -2]
mergedData <- merge(snvMatrixBooleanMelted, cnvData, all.x=T, by=c('gene', 'sample')); # snv and CN data are merged

mergedData$affected <- ifelse(mergedData$value == 1, 1, 0)

tempSamplesFiles <- tempfile("tempSamplesFiles", tempdir())
write.table(snvMatrixBooleanMelted, file=tempSamplesFiles, quote=F, sep="\t", 
		row.names=F, col.names=F)
command <- paste("java -jar", javaLibFile, "OrderSamplesByGeneRank", tempSamplesFiles);
system(command)
orderedSamples <- read.delim(paste(tempSamplesFiles, "_out", sep=""), header=F)

orderedSamples$sampleRank <- 1:nrow(orderedSamples)
names(orderedSamples)[1] <- "Sample"
desiredDataMatrixMelted <- melt(snvMatrix)
names(desiredDataMatrixMelted) <- c("Sample", "Gene", "Mutation")
desiredDataMatrixMelted <- merge(orderedSamples, desiredDataMatrixMelted, by="Sample")
names(cnvData) <- c("Gene", "Sample", "Change")
desiredDataMatrixMelted <- merge(desiredDataMatrixMelted, cnvData, all.x=T)

### Sort the genes
affectedGenes <- desiredDataMatrixMelted[!is.na(desiredDataMatrixMelted$Mutation), ]
desiredDataMatrixMelted$geneCount <- ifelse(is.na(desiredDataMatrixMelted$Mutation), 0, 1)
geneCounts <- summaryBy(geneCount ~ Gene, data=desiredDataMatrixMelted, FUN=sum)
names(geneCounts)[2] <- "geneFreq"
desiredDataMatrixMelted <- merge(geneCounts, desiredDataMatrixMelted, by = "Gene")

snvPlot <- ggplot(desiredDataMatrixMelted) + 
		geom_raster(aes(x=reorder(Sample, -sampleRank), 
						y=reorder(Gene, -geneFreq), fill=Mutation)) +scale_fill_brewer(palette="Dark2") +  theme_bw()

snvPlot <- snvPlot + labs(x = "Sample ID", y = "Gene", title="Gene Status")
# Now, draw outlines around genes with CN differences
rawPlot <- ggplot_build(snvPlot)

# Draw grey rectngles around everything!
maxX <- length(rawPlot$panel$ranges[[1]]$x.labels)
maxY <- length(rawPlot$panel$ranges[[1]]$y.labels)
defaultOutline <- data.frame(x=rep(1:maxX, each=maxY), y=1:maxY)

snvPlot <- snvPlot + geom_rect(data=defaultOutline, 
		aes(xmin=x-0.5, xmax=x+0.5,
				ymin=y-0.5, ymax=y+0.5), color="grey", size=2,fill=NA)


outLines <- data.frame(xCoords=numeric(), yCoords=numeric())
n <- nrow(cnvData)
for (i in 1:n)
{
	outLines[i, c(1:2)] <- getXandYforGeneAndSample(gene=cnvData[i, 1], sample=cnvData[i, 2], rawPlot=rawPlot)
}
outLines$Change <- cnvData$Change


snvPlot <- snvPlot + geom_rect(data=outLines, 
		aes(xmin=xCoords-0.5, xmax=xCoords+0.5,
				ymin=yCoords-0.5, ymax=yCoords+0.5, color=Change),
		fill=NA, size=2)

ggsave(outFile, snvPlot, width=10*1.6, height=10, units="in", dpi=300)

