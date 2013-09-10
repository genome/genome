#!/usr/bin/env perl
#use warnings;
use strict;
#use diagnostics;
use File::Basename qw(basename);

#### Include get options and usage
#
# (1) common 5' and 3' partners
# (2) 
#

######################################################################
#                       Process ChimeraScan Results                  #
######################################################################
#my $dir = '/gscmnt/gc5105/research/cmaher/LungCancer/WIS/CSCAN/RESULTS/';
#my $dir = '/gscmnt/gc5105/research/cmaher/LungCancer/17genomes/CSCAN/RESULTS/';
#my $dir = '/gscmnt/gc6127/research/cmaher/4Malachi/RAW/';
#my $dir = '/gscmnt/gc6127/research/cmaher/Myeloma/MMY5/';
#my $dir = '/gscmnt/gc6127/research/cmaher/ALL/ALL1/';
my $dir = '/gscmnt/gc7001/info/model_data/2892653587/build139185824/fusions/';

#my $outfile = '/gscmnt/gc5105/research/cmaher/LungCancer/17genomes/LUC9/cnv.bed';
#open(OUTFILE, ">$outfile");

# Mitelman database
my $mitelman = '/gscuser/cmaher/References/molclingene.dat';
my (%MITEL,%MITEL5P,%MITEL3P);
open(MITEL, "<$mitelman" ) or die("Couldn't open MITELMAN file $mitelman \n");
while(<MITEL>)  { 
    chomp;
    my(@v)=split(/\t/);
    
    if($v[5] =~ /\//){
	$MITEL{$v[5]}[0]++; # Raise SNV Sample Count
    }
}
close(MITEL);

my @fusions = keys %MITEL;
foreach my $fusion (@fusions){
    my($FP,$TP) = split(/\//,$fusion);
    $MITEL5P{$FP}[0]++;
    $MITEL3P{$TP}[0]++;
}

my @FPgenes = keys %MITEL5P;
my @TPgenes = keys %MITEL3P;

# Drug Targets



# Kinases
my $kinase = '/gscmnt/gc5105/research/cmaher/Reference/Phosphatases_and_Kinases.txt';
my (%KIN);
open(K, "<$kinase" ) or die("Couldn't open K file $kinase \n");
while(<K>)  { 
    chomp;
    my(@v)=split(/\t/);
    
    if($v[1] ne ''){
	$KIN{$v[1]}[0]++; # Raise SNV Sample Count
    }
}
close(K);


# Cancer Genes
my $cancer = '/gscmnt/gc5105/research/cmaher/Reference/CancerGenes.1172012.txt';
my %CANCER;
open(C,"<$cancer") or die("Couldn't open C file $cancer \n");
while(<C>){
    chomp;
    my(@v)=split(/\t/);
    
    $CANCER{$v[0]}[0]++;
}
close(C);






my (%GF,%FPR,%TPR);
my @samples;
#
# Extract Recurrent
#
my $x = '5'; # Element in array
#my $dir = '/gscmnt/gc5105/research/cmaher/LungCancer/17genomes/CSCAN/RESULTS/';
foreach my $gffile (glob $dir."*chimeras.bedpe"){
    my($sample,$c,$file_ext) = split(/\./,basename $gffile);
    push @samples, $sample;
    
    my(%FPPROM,%TPPROM);
    open(GF, "<$gffile" ) or die("Couldn't open GF file file $gffile \n");
    while(<GF>)  { 
	chomp;
	
	if(/^\#/){ # Bypass headers
	    next;
	}else{
	    my(@v)=split(/\t/);
	    my $FP = $v[12];
	    my $TP = $v[13];
	    $FPPROM{$FP}[0]++;
	    $TPPROM{$TP}[0]++;
	}
    }
    close(GF);

    #print "Processing $sample $gffile ...\n";
    open(GF, "<$gffile" ) or die("Couldn't open GF file file $gffile \n");
    while(<GF>)  { 
	chomp;
	
	if(/^\#/){ # Bypass headers
	    next;
	}else{
	    my(@v)=split(/\t/);
	    my $FP = $v[12];
	    my $TP = $v[13];
	    
	    my $score = $v[7];
	    my $total_frag = $v[16];
	    my $span_frag = $v[17];
	    my $type = $v[14]; 
	    
	    if($FP eq $TP){ next; }
	    
	    my $fusion = $FP.':'.$TP;

	    
	    if($v[21] =~ /AAAAAAAAAA/ || $v[21] =~ /TTTTTTTTTT/ || $v[21] =~ /GGGGGGGGGG/ || $v[21] =~ /CCCCCCCCCC/){
		next;
	    }
	 

	    #print "$sample\t$fusion\t$FPPROM{$FP}[0]\t$TPPROM{$TP}[0]\n";
	    if($FPPROM{$FP}[0] <= 3 && $TPPROM{$TP}[0] <= 3){ #filter out fusions where either gene is involved in many different fusions
		#if($fusion =~ /CREB3/){ print "PASS\t$fusion\n"; }
		if($total_frag ne $span_frag){
#		if( ($total_frag > '30' && $total_frag ne $span_frag) || $total_frag <= '30'){
		    $GF{$fusion}[0]++;
		    $GF{$fusion}[1]+=$total_frag;
		    $GF{$fusion}[2]+=$span_frag;
		    $GF{$fusion}[3]=$type;
		    $GF{$fusion}[4]+=$score;
		    $GF{$fusion}[$x] = $span_frag.':'.$total_frag;
		    #$FPR{$FP}[0]++;
		    #$TPR{$TP}[0]++;
		}
	    }
	}
    }
    close(GF);
    $x++;
}

#close(OUTFILE);

# Print header
print "Fusion\t5P\t3P\tTotal_Freq\tSpanning_Freq\tType\tScore\tSpan:Total\tMitel5P\tMitel3P\tKinase5P\tKinase3P\tCancer5P\tCancer3P\n"; #5'PartnersFreq\t3'PartnersFreq";
for(my $x='0';$x<@samples;$x++){ print "\t$samples[$x]"; }
print "\n";

# Print output
my @fusions = keys %GF;
for(my $i = '0'; $i < @fusions; $i++){
    
    # Get count
    my $sample_freq = '0';
    my $total_freq = '0';
    for(my $j = '5'; $j < @samples + 5; $j++){ # Only check data columns
	if($GF{$fusions[$i]}[$j] ne ''){
	    $total_freq++;
	    my($s,$e)=split(/\:/,$GF{$fusions[$i]}[$j]);
	    
	    if($s ne '0'){
#		if($GF{$fusions[$i]}[1] >= '3' && $GF{$fusions[$i]}[2] >= '1'){
		    $sample_freq++;
#		}
	    }
	}
    }
    
    #if($sample_freq eq '1' && $total_freq eq '1' && $GF{$fusions[$i]}[1] >= '2' && ($GF{$fusions[$i]}[1] ne $GF{$fusions[$i]}[2]) ){
    if( ($sample_freq eq '1' || $sample_freq eq '2') && $GF{$fusions[$i]}[1] >= '5' && $GF{$fusions[$i]}[2] >= '1'){
	                                                                                                                                   
	my($gene1,$gene2)=split(/\:/,$fusions[$i]); 
	print "$fusions[$i]\t$gene1\t$gene2";

	# Print all
	for(my $j = '1'; $j < @samples + 5; $j++){ # Exclude first column and instead print actual sample frequency
	    if($GF{$fusions[$i]}[$j] eq ''){
		print "\t0";
	    }else{
		my($s,$e)=split(/\:/,$GF{$fusions[$i]}[$j]);
		print "\t$GF{$fusions[$i]}[$j]";
	    }
	}

#$sample_freq\t$total_freq";
       
	if($MITEL5P{$gene1}[0] ne ''){
	    print "\tMITEL5_".$gene1;
	}else{
	    print "\tNull";
	}
	
	if($MITEL3P{$gene2}[0] ne ''){
	    print "\tMITEL3_".$gene2;
	}else{ print "\tNull"; }
	
	if($KIN{$gene1}[0] ne ''){
	    print "\tKIN5_".$gene1;
	}else{ print "\tNull"; }
	
	if($KIN{$gene2}[0] ne ''){
	    print "\tKIN3_".$gene2;
	}else{ print "\tNull"; }
	
	if($CANCER{$gene1}[0] ne ''){
	    print "\tCAN5_".$gene1;
	}else{ print "\tNull"; }
	
	if($CANCER{$gene2}[0] ne ''){
	    print "\tCAN3_".$gene2;   
	}else{ print "\tNull"; }
	
	if($FPR{$gene1}[0] > '1'){
	    print "\tFPR_".$FPR{$gene1}[0];
	}else{ print "\tNull"; }
	
	if($TPR{$gene2}[0] > '1'){
	    print "\tTPR_".$TPR{$gene2}[0];
	}else{ print "\tNull"; } 
	
	print "\n";
	
    }
}
