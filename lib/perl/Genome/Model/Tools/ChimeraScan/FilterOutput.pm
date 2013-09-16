package Genome::Model::Tools::ChimeraScan::FilterOutput;

use warnings;
use strict;
use File::Basename qw(basename);
use Genome;

my $TOO_MANY_FUSION_PARTNERS = 3;
my $MIN_TOTAL_FRAGS = 5;
my $MIN_SPANNING_FRAGS = 1;

class Genome::Model::Tools::ChimeraScan::FilterOutput {
    is => 'Command',
    has => [
        bedpe_file => {
            is => 'Text',
        },
        output_file => {
            is => 'Text',
        },
    ],
};

sub execute {
    my $self = shift;
######################################################################
#                       Process ChimeraScan Results                  #
######################################################################
#my $dir = '/gscmnt/gc5105/research/cmaher/LungCancer/WIS/CSCAN/RESULTS/';
#my $dir = '/gscmnt/gc5105/research/cmaher/LungCancer/17genomes/CSCAN/RESULTS/';
#my $dir = '/gscmnt/gc6127/research/cmaher/4Malachi/RAW/';
#my $dir = '/gscmnt/gc6127/research/cmaher/Myeloma/MMY5/';
#my $dir = '/gscmnt/gc6127/research/cmaher/ALL/ALL1/';
    #my $dir = '/gscmnt/gc7001/info/model_data/2892653587/build139185824/fusions/';
    #my $bedpe_file = "$dir/chimeras.bedpe";
    my $bedpe_file = $self->bedpe_file;
    my $output_file = $self->output_file;
    my $out = Genome::Sys->open_file_for_writing($output_file);

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

    my(%FPPROM,%TPPROM);
    open(GF, "<$bedpe_file" ) or die("Couldn't open GF file file $bedpe_file \n");
    while(<GF>)  { 
        chomp;

        if(/^\#/){ # Bypass headers
            next;
        }
        else{
            my(@v)=split(/\t/);
            my $FP = $v[12];
            my $TP = $v[13];
            $FPPROM{$FP}[0]++;
            $TPPROM{$TP}[0]++;
        }
    }
    close(GF);

    open(GF, "<$bedpe_file" ) or die("Couldn't open GF file file $bedpe_file \n");
    while(<GF>)  { 
        chomp;

        if(/^\#/){ # Bypass headers
            next;
        }
        else{
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


            if($FPPROM{$FP}[0] <= $TOO_MANY_FUSION_PARTNERS && $TPPROM{$TP}[0] <= $TOO_MANY_FUSION_PARTNERS){
                if($total_frag ne $span_frag){
                    $GF{$fusion}[0]++;
                    $GF{$fusion}[1]+=$total_frag;
                    $GF{$fusion}[2]+=$span_frag;
                    $GF{$fusion}[3]=$type;
                    $GF{$fusion}[4]+=$score;
                    $GF{$fusion}[5] = $span_frag.':'.$total_frag;
                }
            }
        }
    }
    close(GF);


# Print header
    $out->print("Fusion\t5P\t3P\tTotal_Freq\tSpanning_Freq\tType\tScore\tSpan:Total\tMitel5P\tMitel3P\tKinase5P\tKinase3P\tCancer5P\tCancer3P\n"); #5'PartnersFreq\t3'PartnersFreq";

# Print output
    @fusions = keys %GF;
    for my $fusion (@fusions){
        # Get count
        my $sample_freq = '0';
        my $total_freq = '0';
        if($GF{$fusion}[5] ne ''){
            $total_freq++;

            if($GF{$fusion}[2] ne '0'){
                $sample_freq++;
            }
        }

        if( ($sample_freq >= 1 and $sample_freq <= 2) && $GF{$fusion}[1] >= $MIN_TOTAL_FRAGS && $GF{$fusion}[2] >= $MIN_SPANNING_FRAGS){

            my($gene1,$gene2)=split(/\:/,$fusion); 
            $out->print("$fusion\t$gene1\t$gene2");

            # Print all
            for(my $j = '1'; $j <= 5; $j++){
                if($GF{$fusion}[$j] eq ''){
                    $out->print( "\t0");
                }
                else {
                    $out->print( "\t$GF{$fusion}[$j]");
                }
            }

            if($MITEL5P{$gene1}[0] ne ''){
                $out->print("\tMITEL5_".$gene1);
            }else{
                $out->print("\tNull");
            }

            if($MITEL3P{$gene2}[0] ne ''){
                $out->print("\tMITEL3_".$gene2);
            }else{ $out->print("\tNull");}

            if($KIN{$gene1}[0] ne ''){
                $out->print("\tKIN5_".$gene1);
            }else{ $out->print("\tNull");}

            if($KIN{$gene2}[0] ne ''){
                $out->print("\tKIN3_".$gene2);
            }else{ $out->print("\tNull");}

            if($CANCER{$gene1}[0] ne ''){
                $out->print("\tCAN5_".$gene1);
            }else{ $out->print("\tNull");}

            if($CANCER{$gene2}[0] ne ''){
                $out->print("\tCAN3_".$gene2);   
            }else{ $out->print("\tNull");}

            $out->print("\n");

        }
    }
    $out->close;
    return 1;
}
