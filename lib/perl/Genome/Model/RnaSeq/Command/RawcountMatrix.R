# set up command arguments, packages, etc.
# written by Zachary Skidmore

library(edgeR)

args = (commandArgs(TRUE))
matrix_file = args[1]

#################################################################################
######### READ IN AND SET UP THE DATA ###########################################
#################################################################################

data <- read.delim(matrix_file, header=TRUE)
data.tmp <- data

# remove the column gene name in the temp matrix in order to group the genes by condition

data.tmp$Gene_name <- NULL
data.tmp$Ensemble_id <- NULL

# group the data according to a regex given as input in the perl script rawcount.pm

group <- vector()

for (i in 1:ncol(data.tmp))
{
	
	if(grepl("normal", colnames(data.tmp[i])))
	{
		group <- c(group, 1)
	} else if (grepl("tumor", colnames(data.tmp[i]))) {
		group <- c(group, 2)
	}
}

###################################################################################
############# perform default analysis ############################################
###################################################################################

# put the data in a DGEList object

data.object <- DGEList(counts=data[,3:ncol(data)], genes=data[,1:2], group=group)

# normalize the data, estimate dispersion, etc.

data.object <- calcNormFactors(data.object)

data.object <- estimateCommonDisp(data.object)

data.object <- estimateTagwiseDisp(data.object)

# perform the differential expression

et <- exactTest(data.object)

#####################################################################################
############### write output to a file ##############################################
#####################################################################################

# create path/name of file to output

DE_file <- matrix_file

DE_file <- sub('.tsv', '.edgeR.tsv', DE_file, fixed=TRUE)

# extract relevant data from the DGE object containing logFC, PValues etc.

DE_data <- et$genes
DE_data <- cbind(DE_data, et$table)
DE_data$FDR <- p.adjust(method="fdr", p=DE_data$PValue)

# write information to file

write.table(DE_data, file=DE_file, sep="\t", row.names=FALSE)

###################################################################################
################ generate microarray MA-plot equivalent for rna-seq ###############
###################################################################################

# create path/name of file to output

plot_file <- matrix_file

plot_file <- sub('.tsv', '.edgeR.plot.png', plot_file, fixed=TRUE)

png(filename=plot_file)

# generate the plot

de <- decideTestsDGE(et)

detags <- rownames(data.object)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h = c(-2, 2), col = "blue")

dev.off()


