library(plotrix)

#Define the event types to use
event_types_sample=c("snv_sample_count","indel_sample_count","cnv_gain_sample_count","rna_cufflinks_absolute_sample_count")
event_types_gene=c("snv_gene_count","indel_gene_count","cnv_gain_gene_count","rna_cufflinks_absolute_gene_count")

event_names=c("SNVs", "InDels", "CNVs", "Expression")

#Load input data files
datadir = "/Users/mgriffit/Documents/WASHU/Projects/LUC/LUC1-20 analysis/druggable_genes/all_cases_four_event_types"
#datadir = "/Users/mgriffit/Dropbox/Documents/Analysis_development/no_luc9_four_event_types/"
setwd(datadir)
dir()

#1.) Known genes
kg = read.table(file="KnownGeneDrugTable_SantaMonica.txt", sep="\t", header=TRUE, as.is=c(1,3,5,6,7,8,10,11,13,15,17,18,19))
names(kg[c(2,4,9,12,14,16)])

#Create a color heatmap of the druggable gene events
counts=kg[,event_types_sample]
row.names(counts) = kg[,"gene"]

#Define a color array where all colors beyond nmax = "red"
counts_ceiling = counts
nmax=4
for(i in 1:dim(counts_ceiling)[2]){
	if (length(which(counts_ceiling[,i] >= nmax) > 0)){
	  counts_ceiling[which(counts_ceiling[,i] > nmax), i] = nmax
	}
}
color_array = color.scale(as.matrix(counts_ceiling), extremes=c("white","red"))

#Add some extra margins: c(bottom, left, top, right)
title="Affected known druggable gene targets (# patients affected by event type)"
par(mar=c(5, 6, 4, 15) + 0.1)
color2D.matplot(counts, cellcolors=color_array, axes=FALSE, xlab=NA, ylab=NA, show.legend=FALSE, nslices=nmax, show.values=TRUE, vcex=0.4, main=title)
axis(side=1, at=c(0.5:(dim(counts)[2]-0.5)), labels=event_names, las=1)
axis(side=2, at=c((dim(counts)[1]-0.5):0.5), labels=row.names(counts), las=2, cex.axis=0.6, font=3)
axis(side=4, at=c((dim(counts)[1]-0.5):0.5), labels=kg[,"drug_list_abr"], las=2, cex.axis=0.6)
legend_labels = c(0:(nmax-1), paste(nmax,"+",sep=""))
color.legend(0,-5,2,-7, legend_labels,rect.col=color.scale(0:nmax, extremes=c("white","red")))
mtext("Patient Event Count", side=1, at=2.5, line=3.5)
                     

#2.) Gene families and pathways
select_families = c("Kinase_GO0016301", "Kinase_KinBase", "Kinase_Kincat", "Kinase_RonBose", "ProteinKinase_Entrez", "ReceptorTyrosineKinases_RonBose", "TyrosineKinase_RonBose", "TyrosineKinase_Entrez", "TyrosineKinase_GO0004713", "SerineThreonineKinase_GO0004674", "Phospholipase_GO0004620", "ProteaseInhibitorActivity_GO0030414", "ProteinPhosphatase_GO0004721", "Phosphatase_Entrez", "IonChannel_GO0005216", "GProteinCoupledReceptor_GO0004930", "GrowthFactor_GO0008083", "HormoneActivity_GO0005179", "NuclearHormoneReceptor_GO0004879", "CellSurface")

select_families_custom = c("Kinase_RonBose", "ReceptorTyrosineKinases_RonBose", "TyrosineKinase_RonBose", "SerineThreonineKinase_GO0004674", "Phospholipase_GO0004620", "ProteaseInhibitorActivity_GO0030414", "ProteinPhosphatase_GO0004721", "IonChannel_GO0005216", "GProteinCoupledReceptor_GO0004930", "GrowthFactor_GO0008083", "HormoneActivity_GO0005179", "NuclearHormoneReceptor_GO0004879", "CellSurface")

select_pathways_ccmp = c("Alpha6Beta4IntegrinPathway_CCMP", "AndrogenReceptorPathway_CCMP", "EGFR1Pathway_CCMP", "HedgehogPathway_CCMP", "IDPathway_CCMP", "KitReceptorPathway_CCMP", "NotchPathway_CCMP", "TGFBRPathway_CCMP", "TNFAlphaNFkBPathway_CCMP", "WntPathway_CCMP")

select_pathways_kegg1 = c("AdherensJunction_KEGGhsa04520", "Apoptosis_KEGGhsa04210", "CellCycle_KEGGhsa04110", "CytokineReceptor_KEGGhsa04060", "EcmReceptor_KEGGhsa04512", "ERBBSignaling_KEGGhsa04012", "FocalAdhesion_KEGGhsa04510", "JakStatSignaling_KEGGhsa04630", "MAPKSignaling_KEGGhsa04010", "PPARSignaling_KEGGhsa03320", "mTORSignaling_KEGGhsa04150", "TGFbetaSignaling_KEGGhsa04350", "TP53Signaling_KEGGhsa04115", "VEGFSignaling_KEGGhsa04370", "WntSignaling_KEGGhsa04310")

select_pathways_kegg2 = c("AmlCancer_KEGGhsa05221", "BasalCellCancer_KEGGhsa05217", "BladderCancer_KEGGhsa05219", "CmlCancer_KEGGhsa05220", "ColorectalCancer_KEGGhsa05210", "EndometrialCancer_KEGGhsa05213", "GliomaCancer_KEGGhsa05214", "MelanomaCancer_KEGGhsa05218", "NonSmallCellLungCancer_KEGGhsa05223", "PancreaticCancer_KEGGhsa05212", "ProstateCancer_KEGGhsa05215", "RenalCellCancer_KEGGhsa05211", "SmallCellLungCancer_KEGGhsa05222", "ThyroidCancer_KEGGhsa05216")

#Import gene group data file
gf = read.table(file="GeneFamilyTable.tsv", sep="\t", header=TRUE)

#Make gene group labels
gf[,"gene_group_labels"] = paste(gf[,"gene_group"], "\n (n = ", gf[,"gene_group_size"], " genes)", sep="")

#Calculate the proportion of each gene group that is affected by one of these event types
gf[,"percent_affected"] = round((gf[,"grand_gene_count"] / gf[,"gene_group_size"])*100, digits=1)
row.names(gf) = gf[,"gene_group"]

#2-a.) Create a color heatmap of the gene family events summarizing counts at the PATIENT level
createGeneGroupPlot_PatientCounts = function(data, selection, event_types, title){
	counts=data[selection,event_types]
	tmp = data[selection,]
	o=order(tmp[,"percent_affected"], decreasing=TRUE)
	counts=counts[o,]
	par(mar=c(5, 13.5, 4, 5) + 0.1)
	color2D.matplot(counts, extremes=c("white","blue"), axes=FALSE, xlab=NA, ylab=NA, show.legend=TRUE, show.values=TRUE, vcex=1, main=title)
	axis(side=1, at=c(0.5:(dim(counts)[2]-0.5)), labels=event_names, las=1)
	axis(side=2, at=c((dim(counts)[1]-0.5):0.5), labels=tmp[o, "gene_group_labels"], las=2, cex.axis=0.75)
	percent_labels = paste(tmp[o,"percent_affected"], "%", sep="")
	axis(side=4, at=c((dim(counts)[1]-0.5):0.5), labels=percent_labels, las=2,cex.axis=0.75)
	mtext("Patient Event Count", side=1, line=3.5, at=1.50)
}

#Gene families
createGeneGroupPlot_PatientCounts(gf, select_families, event_types_sample, "Affected druggable gene families (# patients affected by event type)")
createGeneGroupPlot_PatientCounts(gf, select_families_custom, event_types_sample, "Affected druggable gene families (# patients affected by event type)")

#Pathways CCMP
createGeneGroupPlot_PatientCounts(gf, select_pathways_ccmp, event_types_sample, "Affected Cancer Cell Map Pathways (# patients affected by event type)")

#Pathways KEGG - cancer pathways
createGeneGroupPlot_PatientCounts(gf, select_pathways_kegg1, event_types_sample, "Affected KEGG Cancer Pathways (# patients affected by event type)")

#Pathways KEGG - cancer subtype pathways
createGeneGroupPlot_PatientCounts(gf, select_pathways_kegg2, event_types_sample, "Affected KEGG Cancer Subtype Pathways (# patients affected by event type)")


#2-b.) Create a color heatmap of the gene family events summarizing counts at the GENE level
createGeneGroupPlot_GeneCounts = function(data, selection, event_types, title, cool_col, hot_col){
	counts=data[selection,event_types]
	tmp = data[selection,]
	o=order(tmp[,"percent_affected"], decreasing=TRUE)
	counts=counts[o,]
	#Put a ceiling to reduce the display effect of outliers
	counts_ceiling = counts
	nmax=100
	for(i in 1:dim(counts_ceiling)[2]){
		if (length(which(counts_ceiling[,i] >= nmax) > 0)){
	  	counts_ceiling[which(counts_ceiling[,i] > nmax), i] = nmax
		}
	}
	color_array = color.scale(as.matrix(counts_ceiling), extremes=c(cool_col, hot_col))
	par(mar=c(5, 13.5, 4, 5) + 0.1)
	color2D.matplot(counts, cellcolors=color_array, axes=FALSE, xlab=NA, ylab=NA, show.legend=TRUE, show.values=TRUE, vcex=0.8, main=title)
	axis(side=1, at=c(0.5:(dim(counts)[2]-0.5)), labels=event_names, las=1)
	axis(side=2, at=c((dim(counts)[1]-0.5):0.5), labels=tmp[o, "gene_group_labels"], las=2, cex.axis=0.75)
	percent_labels = paste(tmp[o,"percent_affected"], "%", sep="")
	axis(side=4, at=c((dim(counts)[1]-0.5):0.5), labels=percent_labels, las=2,cex.axis=0.75)
	mtext("Affected Gene Count", side=1, line=3.5, at=1.80)
}

#Gene families
createGeneGroupPlot_GeneCounts(gf, select_families, event_types_gene, "Affected druggable gene families (# genes affected by event type)", "white", "blue")
createGeneGroupPlot_GeneCounts(gf, select_families_custom, event_types_gene, "Affected druggable gene families (# genes affected by event type)", "white", "blue")

#Pathways CCMP
createGeneGroupPlot_GeneCounts(gf, select_pathways_ccmp, event_types_gene, "Affected Cancer Cell Map Pathways (# genes affected by event type)", "white", "blue")

#Pathways KEGG - cancer pathways
createGeneGroupPlot_GeneCounts(gf, select_pathways_kegg1, event_types_gene, "Affected KEGG Cancer Pathways (# genes affected by event type)", "white", "blue")

#Pathways KEGG - cancer subtype pathways
createGeneGroupPlot_GeneCounts(gf, select_pathways_kegg2, event_types_gene, "Affected KEGG Cancer Subtype Pathways (# genes affected by event type)", "white", "blue")


#Create individual patient heat map
kg = read.table(file="KnownGeneDrugTable_IndividualPatients.txt", header=TRUE, sep="\t", as.is=1:27)
patient_list = c("LUC1", "LUC2", "LUC4", "LUC8", "LUC9", "LUC10", "LUC12", "LUC13", "LUC14", "LUC17", "LUC18", "LUC20", "LUC6", "LUC7", "LUC11", "LUC15", "LUC16")

library(gplots)
color_array = kg[,patient_list]
class_array = kg[,patient_list]
row.names(class_array) = kg[,"gene"]
for (i in 1:length(patient_list)){
	color_array[which(is.na(color_array[,i])), i] = col2hex("white")
	class_array[which(is.na(class_array[,i])), i] = 0
	color_array[which(color_array[,i] == "snv"), i] = "#E7298A"
	class_array[which(class_array[,i] == "snv"), i] = 1
	color_array[which(color_array[,i] == "indel"), i] = "#D95F02"
	class_array[which(class_array[,i] == "indel"), i] = 1
	color_array[which(color_array[,i] == "cnv_gain"), i] = "#7570B3"
	class_array[which(class_array[,i] == "cnv_gain"), i] = 1
	color_array[which(color_array[,i] == "rna_cufflinks_absolute"), i] = "#1B9E77"
	class_array[which(class_array[,i] == "rna_cufflinks_absolute"), i] = 1
}


#Add some extra margins: c(bottom, left, top, right)
title="Affected known druggable gene targets (patients affected by event type)"
par(mar=c(6, 6, 4.5, 15) + 0.1) #Add some extra margins: c(bottom, left, top, right)
color2D.matplot(as.matrix(class_array), cellcolors=as.matrix(color_array), axes=FALSE, xlab=NA, ylab=NA, show.legend=FALSE, show.values=FALSE, vcex=0.4, main=title)
axis(side=1, at=c(0.5:(dim(class_array)[2]-0.5)), labels=patient_list, las=2)
axis(side=2, at=c((dim(class_array)[1]-0.5):0.5), labels=row.names(class_array), las=2, cex.axis=0.6, font=3)
axis(side=4, at=c((dim(class_array)[1]-0.5):0.5), labels=kg[,"drug_list_abr"], las=2, cex.axis=0.6)
legend_labels = c("SNV", "InDel", "CNV Amp", "Over Exp.")
color.legend(0,-6.5,10,-8, legend_labels, rect.col=c("#E7298A", "#D95F02", "#7570B3", "#1B9E77"))
mtext("Event Types", side=1, at=11.25, line=4.25, font=2)
abline(v=12, lwd=4, lty=1)
mtext("Smokers", side=3, at=6, line=0.25, font=2)
mtext("Never Smokers", side=3, at=14.55, line=0.25, font=2)


