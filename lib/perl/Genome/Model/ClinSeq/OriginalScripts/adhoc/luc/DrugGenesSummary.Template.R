library(plotrix)

#Define the event types to use
event_types_sample=c("snv_sample_count","indel_sample_count","cnv_gain_sample_count","rna_cufflinks_absolute_sample_count","rna_tophat_absolute_sample_count")
event_types_gene=c("snv_gene_count","indel_gene_count","cnv_gain_gene_count","rna_cufflinks_absolute_gene_count","rna_tophat_absolute_gene_count")

event_names=c("SNVs", "InDels", "CNVs", "Gene\nExpression", "Junction\nExpression")

#Load input data files
datadir = "/Users/mgriffit/Dropbox/Documents/Analysis_development/all_cases_all_event_types/"
datadir = "/Users/mgriffit/Dropbox/Documents/Analysis_development/no_luc9_all_event_types/"
setwd(datadir)
dir()

#1.) Known genes
kg = read.table(file="KnownGeneDrugTable_DrugBank.tsv", sep="\t", header=TRUE, as.is=c(1,3,5,7,9,11,13,15))
kg = read.table(file="KnownGeneDrugTable_SantaMonica.tsv", sep="\t", header=TRUE, as.is=c(1,3,5,7,9,11,13,15))

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
par(mar=c(5, 4, 4, 10) + 0.1)
color2D.matplot(counts, cellcolors=color_array, axes=FALSE, xlab=NA, ylab=NA, show.legend=FALSE, nslices=nmax, show.values=TRUE, vcex=0.4, main=title)
axis(side=1, at=c(0.5:(dim(counts)[2]-0.5)), labels=event_names, las=1)
axis(side=2, at=c((dim(counts)[1]-0.5):0.5), labels=row.names(counts), las=2,cex.axis=0.5, font=3)
axis(side=4, at=c((dim(counts)[1]-0.5):0.5), labels=kg[,"drug_list"], las=2,cex.axis=0.5)
legend_labels = c(0:(nmax-1), paste(nmax,"+",sep=""))
color.legend(0,-7,2,-9, legend_labels,rect.col=color.scale(0:nmax, extremes=c("white","red")))
mtext("Patient Event Count", side=1, line=3.75)
                     

#2.) Gene families and pathways
select_families = c("Kinase_GO0016301", "Kinase_KinBase", "Kinase_Kincat", "Kinase_RonBose", "ProteinKinase_Entrez", "ReceptorTyrosineKinases_RonBose", "TyrosineKinase_RonBose", "TyrosineKinase_Entrez", "TyrosineKinase_GO0004713", "SerineThreonineKinase_GO0004674", "Phospholipase_GO0004620", "ProteaseInhibitorActivity_GO0030414", "ProteinPhosphatase_GO0004721", "Phosphatase_Entrez", "IonChannel_GO0005216", "GProteinCoupledReceptor_GO0004930", "GrowthFactor_GO0008083", "HormoneActivity_GO0005179", "NuclearHormoneReceptor_GO0004879", "CellSurface")

select_families_custom = c("Kinase_RonBose", "ReceptorTyrosineKinases_RonBose", "TyrosineKinase_RonBose", "SerineThreonineKinase_GO0004674", "Phospholipase_GO0004620", "ProteaseInhibitorActivity_GO0030414", "ProteinPhosphatase_GO0004721", "IonChannel_GO0005216", "GProteinCoupledReceptor_GO0004930", "GrowthFactor_GO0008083", "HormoneActivity_GO0005179", "NuclearHormoneReceptor_GO0004879", "CellSurface")

select_pathways_ccmp = c("Alpha6Beta4IntegrinPathway_CCMP", "AndrogenReceptorPathway_CCMP", "EGFR1Pathway_CCMP", "HedgehogPathway_CCMP", "IDPathway_CCMP", "KitReceptorPathway_CCMP", "NotchPathway_CCMP", "TGFBRPathway_CCMP", "TNFAlphaNFkBPathway_CCMP", "WntPathway_CCMP")

select_pathways_kegg1 = c("AdherensJunction_KEGGhsa04520", "Apoptosis_KEGGhsa04210", "CellCycle_KEGGhsa04110", "CytokineReceptor_KEGGhsa04060", "EcmReceptor_KEGGhsa04512", "ERBBSignaling_KEGGhsa04012", "FocalAdhesion_KEGGhsa04510", "JakStatSignaling_KEGGhsa04630", "MAPKSignaling_KEGGhsa04010", "PPARSignaling_KEGGhsa03320", "mTORSignaling_KEGGhsa04150", "TGFbetaSignaling_KEGGhsa04350", "TP53Signaling_KEGGhsa04115", "VEGFSignaling_KEGGhsa04370", "WntSignaling_KEGGhsa04310")

select_pathways_kegg2 = c("AmlCancer_KEGGhsa05221", "BasalCellCancer_KEGGhsa05217", "BladderCancer_KEGGhsa05219", "CmlCancer_KEGGhsa05220", "ColorectalCancer_KEGGhsa05210", "EndometrialCancer_KEGGhsa05213", "GliomaCancer_KEGGhsa05214", "MelanomaCancer_KEGGhsa05218", "NonSmallCellLungCancer_KEGGhsa05223", "PancreaticCancer_KEGGhsa05212", "ProstateCancer_KEGGhsa05215", "RenalCellCancer_KEGGhsa05211", "SmallCellLungCancer_KEGGhsa05222", "ThyroidCancer_KEGGhsa05216")

#Import gene group data file
gf = read.table(file="GeneFamilyTable.tsv", sep="\t", header=TRUE, as.is=c(1,3,5,8,10,13,15,18,20,23,25,28,30))

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
	mtext("Patient Event Count", side=1, line=3.5, at=1.80)
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
createGeneGroupPlot_GeneCounts = function(data, selection, event_types, title){
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
	color_array = color.scale(as.matrix(counts_ceiling), extremes=c("white","dark green"))
	par(mar=c(5, 13.5, 4, 5) + 0.1)
	color2D.matplot(counts, cellcolors=color_array, axes=FALSE, xlab=NA, ylab=NA, show.legend=TRUE, show.values=TRUE, vcex=0.8, main=title)
	axis(side=1, at=c(0.5:(dim(counts)[2]-0.5)), labels=event_names, las=1)
	axis(side=2, at=c((dim(counts)[1]-0.5):0.5), labels=tmp[o, "gene_group_labels"], las=2, cex.axis=0.75)
	percent_labels = paste(tmp[o,"percent_affected"], "%", sep="")
	axis(side=4, at=c((dim(counts)[1]-0.5):0.5), labels=percent_labels, las=2,cex.axis=0.75)
	mtext("Affected Gene Count", side=1, line=3.5, at=1.80)
}

#Gene families
createGeneGroupPlot_GeneCounts(gf, select_families, event_types_gene, "Affected druggable gene families (# genes affected by event type)")
createGeneGroupPlot_GeneCounts(gf, select_families_custom, event_types_gene, "Affected druggable gene families (# genes affected by event type)")

#Pathways CCMP
createGeneGroupPlot_GeneCounts(gf, select_pathways_ccmp, event_types_gene, "Affected Cancer Cell Map Pathways (# genes affected by event type)")

#Pathways KEGG - cancer pathways
createGeneGroupPlot_GeneCounts(gf, select_pathways_kegg1, event_types_gene, "Affected KEGG Cancer Pathways (# genes affected by event type)")

#Pathways KEGG - cancer subtype pathways
createGeneGroupPlot_GeneCounts(gf, select_pathways_kegg2, event_types_gene, "Affected KEGG Cancer Subtype Pathways (# genes affected by event type)")





