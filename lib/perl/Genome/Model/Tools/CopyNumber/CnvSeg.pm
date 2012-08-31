package Genome::Model::Tools::CopyNumber::CnvSeg;

use strict;
use warnings;
use Getopt::Std;
use Statistics::Descriptive;
#use FindBin qw($Bin);
#use lib "$FindBin::Bin";
#use PURITYHMM;

class Genome::Model::Tools::CopyNumber::CnvSeg {
    is => 'Command',
    has_input => [
        copy_number_file => {
            is => 'String',
            doc => 'Copy number file to operate on',
            is_output => 1,
        },
        output_file => {
            is => 'String',
            doc => 'output file path',
            is_output => 1,
        },
    ],
    has_optional_input => [
        max_copy_number => {
            is => 'String',
            default => '3',
            doc => 'Maximal Copy Number',
        },
        min_arm_size => {
            is => 'String',
            doc => ' Minimal size of the arm for ploidy estimation',
            default => 1.1e7,
        },
        min_markers => {
            is => 'String',
            doc => ' Minimal number of Markers required to report a copy number segment',
            default => '5',
        },
        min_llr => {
            is => 'String',
            doc => 'Minimal LLR score required to report a copy number segment ',
            default => '30',
        },
        centromere_file => {
            is => 'String',
            doc => 'Path to a UCSC centromere table (replace hg18 w/ hg19 for build 37)',
            default => "/gscmnt/sata186/info/medseq/kchen/work/SolexaCNV/scripts/centromere.hg18.csv",
        },
        gap_file => {
            is => 'String',
            doc => 'Path to a UCSC gap table (replace hg18 w/ hg19 for build 37)',
            default => "/gscmnt/sata186/info/medseq/kchen/work/SolexaCNV/scripts/hg18gaps.csv",
        },
        purity => {
            is => 'String',
            doc => 'purity ratio',
            default => 1.0,
        },
        chromosome_number => {
            is => 'String',
            doc => 'Chromosome to operate on',
        },
        estimate_ploidy => {
            is => 'Boolean',
            doc => "estimate ploidy for each chromosomal arm",
        },
        include_sex_chroms => {
            is => 'Boolean',
            doc => "include sex chromosomes",
        },
    ],
    has_transient_optional => [
        ofh => {
            doc => 'output file handle',
        },
    ],
    has_param => [
        lsf_queue => {
            default => "long",
        }
    ],
};


sub help_detail {
    "This script produces ploidy and copy number variants from the copy number files produced by BamToCn"
}

sub execute {
    my $self = shift;
    my $version="CNVseg-0.0.1r1";
    my %opts;

    $opts{c} = $self->chromosome_number;
    $opts{a} = $self->estimate_ploidy;
    $opts{x} = $self->include_sex_chroms;
    $opts{y} = $self->max_copy_number;
    $opts{s} = $self->min_arm_size;
    $opts{n} = $self->min_markers;
    $opts{l} = $self->min_llr;
    $opts{C} = $self->centromere_file;
    $opts{G} = $self->gap_file;
    $opts{p} = $self->purity;

    my $fcna = $self->copy_number_file;
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    $self->ofh($output_fh);
    my $f_centromere=$opts{C} if(defined $opts{C});
    my $f_gaps=$opts{G} if(defined $opts{G});

    my $method='vit';
    my $TPupdate=0;
    my $minNumMarker=$opts{n};
    my $maxExtInterval=450000;
    my $minLLR=$opts{l};

    my %Centromere;
    if($opts{a}){
        open(CEN,"<$f_centromere") || die "unable to find centromere coordinates file\n";
        while(<CEN>){
            chomp;
            my ($bin,$chr,$start,$end,$ix,$n,$size,$type,$bridge)=split;
            $chr=~s/chr//;
            $Centromere{$chr}->{start}=$start;
            $Centromere{$chr}->{end}=$end;
        }
        close(CEN);
    }

    my %GapStart;
    my %GapEnd;
    open(GAP,"<$f_gaps") || die "unable to find sequencing gap coordinates file\n";
    while(<GAP>){
        chomp;
        my ($bin,$chr,$start,$end,$ix,$n,$size,$type,$bridge)=split;
        $chr=~s/chr//;
        push @{$GapStart{$chr}},$start;
        push @{$GapEnd{$chr}},$end;
    }
    close(GAP);


    open(FCNA,"<$fcna") || die "unable to open $fcna\n";
    my $pchr=0;
    my %nmarkers;
    my $nmarker=0;
    my @chrs;
    my ($cnv,$cnv_n);
    my (@poses,@cn);
    my %median_readcount;
    while(<FCNA>){
        chomp;
        my ($chr,$pos,$tumor,$cn,@extra)=split;
        if(/\#2xReadCount\:(\S+)/){
            $median_readcount{normal}=$1;
        }
        next if(defined $opts{c} && $chr ne $opts{c} || !defined $GapStart{$chr});
        next if($chr=~/[xy]/i && !defined $opts{x});
        next unless($pos=~/^\d+$/);
        if($chr ne $pchr){
            $nmarkers{$pchr}=$nmarker;
            $nmarker=0;
            $pchr=$chr;
            push @chrs,$pchr;
        }

        my @gapstart=@{$GapStart{$chr}};
        my @gapend=@{$GapEnd{$chr}};
        if($gapstart[0]<=$pos && $pos<=$gapend[0]){   #ignore

        }
        elsif($gapend[0]<$pos && $#gapend>0){
            while($#gapend>0 && $gapend[0]<$pos){
                shift @gapstart;
                shift @gapend;
            }
            $GapStart{$chr}=\@gapstart;
            $GapEnd{$chr}=\@gapend;
        }
        else{
            $nmarker++;
            push @poses,$pos;
            push @cn,$tumor;
        }
    }
    close(FCNA);

    #The last chromosome
    $nmarkers{$pchr}=$nmarker;

    my $median_cn2x=$self->Get_Median(\@cn);
    printf $output_fh "--- %d probes loaded\n",$#poses+1;

    $cnv=Genome::Model::Tools::CopyNumber::PurityHmm->create(purity=>$opts{p},mean2x=>$median_cn2x,var=>$median_cn2x,max_cn=>$opts{y},cnv_rate=>0.001,cnv_size=>0.1, copy_number_file => $self->copy_number_file);
    $cnv->init;

    my ($Dall,$Dall_n);
    ($Dall->{x},$Dall->{y})=(\@poses,\@cn);

    $cnv->Fit($Dall,$method,$TPupdate,$opts{u}, $output_fh);

    my $pidx=0;
    my $idx=0;


    foreach my $chr(@chrs){
        $idx=$pidx+$nmarkers{$chr}-1;
        my @x=@poses[$pidx .. $idx];
        my @y=@cn[$pidx .. $idx];
        my $D;
        ($D->{x},$D->{y})=(\@x,\@y);

        printf $output_fh "--- Chromosome %s\t%d probes\n",$chr,$#x+1;
        my ($label,$cn_adjusted);
        if(defined $opts{a}){
            ($label,$cn_adjusted)=$cnv->SegPloidy($output_fh,$D,$opts{s},$median_readcount{$chr},$Centromere{$chr}->{start},$Centromere{$chr}->{end});
        }
        else{
            ($label,$cn_adjusted)=$cnv->SegPloidy($output_fh,$D,$opts{s},$median_readcount{$chr});
        }

        $cnv->printParm($output_fh);

        $label=$self->MedianFilter($label,2*$opts{n}+1);
        my $cnvs=$self->DetectCNV($chr,$cnv,$D,$label,$cn_adjusted,$minNumMarker,$maxExtInterval,$minLLR);
        $self->printCNVs($cnvs,$chr);
        $pidx=$idx+1;
    }

    $cnv->SegPloidy($output_fh,$Dall,$opts{s});
    printf $output_fh "--- Whole Genome %d probes\n",$#{$Dall->{x}}+1;
    $cnv->printParm($output_fh,1);
    return 1;
}

sub DetectCNV {
    my $self = shift;
    my ($chr,$hmm,$data,$label,$cn_adjusted,$minNumMarker,$maxExtInterval,$minLLR)=@_;
    my @y_data = @{$data->{y}};
    my $n=$#y_data + 1;
    my $extend=0;
    my $interval=0;
    my @CNVs;
    my ($CNVstart,$CNVend,$CNVlevel);
    my ($i_start,$i_end);
    my $stat=Statistics::Descriptive::Sparse->new();
    my $cn=Statistics::Descriptive::Sparse->new();
    for(my $i=0;$i<$n;$i++){
        $interval=${$data->{x}}[$i]-${$data->{x}}[$i-1] if($i>0);
        if($extend==0 && $$label[$i]!=2){   #CNV start
            $CNVstart=${$data->{x}}[$i];
            $i_start=$i;
            $CNVlevel=$$label[$i];
            $extend=1;
            $stat->add_data(${$data->{y}}[$i]);
            $cn->add_data($$cn_adjusted[$i]);
        }
        elsif(($extend==1) && (($$label[$i]!=$CNVlevel)||($i==$n-1))){  #CNV stop due to state change
            if($stat->count()>=$minNumMarker){
                my $log10Ratio=$hmm->LogLikelihoodRatio($data,$i_start,$i_end,$CNVlevel);
                my $cnv;
                ($cnv->{start},$cnv->{end},$cnv->{size},$cnv->{nmarker},$cnv->{level},$cnv->{level_obs},$cnv->{LLR})=($CNVstart,$CNVend,$CNVend-$CNVstart,$stat->count(),$CNVlevel,$cn->mean(),$log10Ratio);
                if($log10Ratio>$minLLR){
                    push @CNVs,$cnv;
                }
            }
            $extend=0;
            undef $stat;
            $stat=Statistics::Descriptive::Sparse->new();
            undef $cn;
            $cn=Statistics::Descriptive::Sparse->new();
        }
        elsif(($extend==1) && ($interval>$maxExtInterval)){  #CNV stop due to large gaps
            if($stat->count()>=$minNumMarker){
                my $log10Ratio=$hmm->LogLikelihoodRatio($data,$i_start,$i_end,$CNVlevel);
                my $cnv;
                ($cnv->{start},$cnv->{end},$cnv->{size},$cnv->{nmarker},$cnv->{level},$cnv->{level_obs},$cnv->{LLR})=($CNVstart,$CNVend,$CNVend-$CNVstart,$stat->count(),$CNVlevel,$cn->mean(),$log10Ratio);
                if($log10Ratio>$minLLR){
                    push @CNVs,$cnv;
                }
            }
            undef $stat;
            $stat=Statistics::Descriptive::Sparse->new();
            undef $cn;
            $cn=Statistics::Descriptive::Sparse->new();
            $CNVstart=${$data->{x}}[$i];
            $i_start=$i;
            $CNVlevel=$$label[$i];
            $stat->add_data(${$data->{y}}[$i]);
        }
        elsif($extend==1){   #Extend
            $stat->add_data(${$data->{y}}[$i]);
            $cn->add_data($$cn_adjusted[$i]);
            $CNVend=${$data->{x}}[$i];
            $i_end=$i;
        }
        else{
        }
    }
    my $output_fh = $self->ofh;
    printf $output_fh "%d CNV predicted on chromosome %s\n",$#CNVs+1, $chr;
    return \@CNVs;
}

sub printCNVs {
    my $self = shift;
    my ($cnvs,$chr)=@_;
    my $output_fh = $self->ofh;
    print $output_fh "#CHR\tSTART\tEND\tSIZE\tnMarkers\tCN\tAdjusted_CN\tLLR_Score\n";
    foreach my $cnv(@{$cnvs}){
        printf $output_fh "%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n",$chr||1,$cnv->{start},$cnv->{end},$cnv->{size},$cnv->{nmarker},$cnv->{level},$cnv->{level_obs},$cnv->{LLR};
    }
}

sub Get_Median {
    my $self = shift;
    my $rpole = shift;
    my @pole = @$rpole;
    my $ret;

    @pole=sort{$a<=>$b} @pole;
    if( (@pole % 2) == 1 ) {
        $ret = $pole[((@pole+1) / 2)-1];
    } else {    #$ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
        $ret=$pole[(@pole / 2)-1];
    }
    return $ret;
}

sub MedianFilter {
    my $self = shift;
    my ($label,$win)=@_;
    my @newlabel=@{$label};
    my @buffer;
    for(my $i=0;$i<=$#$label;$i++){
        push @buffer,$$label[$i];
        if($#buffer==$win-1){
            my $median=$self->Get_Median(\@buffer);
            my $idx=$i-($win-1)/2;
            $newlabel[$idx]=$median;
            shift @buffer;
        }
    }
    return \@newlabel;
}

1;
