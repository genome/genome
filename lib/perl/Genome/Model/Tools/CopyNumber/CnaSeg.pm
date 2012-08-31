package Genome::Model::Tools::CopyNumber::CnaSeg;

use strict;
use warnings;
use Genome;
use Getopt::Std;
use Statistics::Descriptive;

class Genome::Model::Tools::CopyNumber::CnaSeg {
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
        chromosome => {
            is => 'String',
            doc => 'Analyze a single chromosome number',
        },
        detect_somatic => {
            is => 'Boolean',
            doc => 'Detect somatic copy number alteration (difference between tumor and normal)',
        },
        purity => {
            is => 'Number',
            doc => 'Prior purity',
            default => '1.0',
        },
        output_purity_file => {
            is => 'String',
            doc => 'Output purity adjusted copy number in a file',
            is_output => 1,
        },
        estimate_purity => {
            is => 'Boolean',
            doc => 'Estimate purity from data',
        },
        estimate_ploidy => {
            is => 'Boolean',
            doc => 'Estimate Ploidy for each chromosomal arms',
        },
        max_copy_number => {
            is => 'Number',
            doc => 'Maximal copy number',
            default => '3',
        },
        min_arm_size => {
            is => 'String',
            doc => 'Minimal size of the arm for ploidy estimation',
            default => '1.1e7',
        },
        min_markers => {
            is => 'String',
            doc => 'Minimal number of Markers required to report a copy number segment',
            default => '5',
        },
        min_llr => {
            is => 'String',
            doc => 'Minimal LLR score required to report a copy number segment',
            default => '10',
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
    ],
    has_transient_optional => [
        ofh => {
            doc => 'output file handle',
        },
    ],
};

sub help_detail {
    return qq{
This script detects copy number varied regions in a tumor genome (relative to a matched normal genome when --detect-somatic is on)
This script produce purity/ploidy estimation from the copy number files produced by BAM2CNA.pl or gmt bam-to-cna
Positions can be masked by cns_marker_files obtained from in silico simulation
Evidences suggest masking significantly reduce the variation of the signal.}
}

sub execute {
    my $self = shift;
    my $version="CNAseg-0.0.1r1";
    my %opts;

    my $fcna = $self->copy_number_file;
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    $self->ofh($output_fh);
    my $f_centromere = $self->centromere_file;
    my $f_gaps = $self->gap_file;

    my $method='vit';
    #$method='baum';
    my $TPupdate=0;
    my $minNumMarker=$self->min_markers;
    my $maxExtInterval=450000;
    my $minLLR=$self->min_llr;

    my %Centromere;
    if($self->estimate_ploidy){
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
    my (@poses,@cnt,@cnn,@cna);
    my %median_readcount;
    while(<FCNA>){
        chomp;
        my ($chr,$pos,$tumor,$normal,$cna,@extra)=split;
        #Chr1 median read count tumor:3905      normal:4517
        if(/\#Chr(\S+) median read count\ttumor\:(\d+)\tnormal\:(\d+)/){
            $median_readcount{$1}{tumor}=$2;
            $median_readcount{$1}{normal}=$3;
        }
        next if(defined $self->chromosome && $chr ne $self->chromosome);
        next unless($pos=~/^\d+$/);
        if($chr ne $pchr){
            $nmarkers{$pchr}=$nmarker;
            $nmarker=0;
            $pchr=$chr;
            push @chrs,$pchr;
        }

        my @gapstart=@{$GapStart{$chr}};
        my @gapend=@{$GapEnd{$chr}};
        if($gapstart[0]<=$pos && $pos<=$gapend[0]){
            #ignore
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
            push @cnt,$tumor;
            push @cnn,$normal;
        }
    }
    close(FCNA);

    #The last chromosome
    $nmarkers{$pchr}=$nmarker;

    my $median_cnt2x=$self->Get_Median(\@cnt);
    printf $output_fh "--- %d probes loaded\n",$#poses+1;
    $cnv=Genome::Model::Tools::CopyNumber::PurityHmm->create(
        purity           => $self->purity,
        mean2x           => $median_cnt2x,
        var              => $median_cnt2x,
        max_cn           => $self->max_copy_number,
        cnv_rate         => 0.001,
        cnv_size         => 0.1,
        copy_number_file => $self->copy_number_file,
    );
    $cnv->init;

    my ($Dall,$Dall_n);
    ($Dall->{x},$Dall->{y},$Dall->{z})=(\@poses,\@cnt,\@cnn);

    $cnv->Fit($Dall, $method, $TPupdate, $self->estimate_purity, $output_fh);
    $cnv->printParm($output_fh);

    my $pidx=0;
    my $idx=0;
    if($self->detect_somatic) {
        my $median_cnn2x=$self->Get_Median(\@cnn);
        $cnv_n=Genome::Model::Tools::CopyNumber::PurityHmm->create(
            purity           => $self->purity,
            mean2x           => $median_cnn2x,
            var              => $median_cnn2x,
            max_cn           => $self->max_copy_number,
            cnv_rate         => 0.001,
            cnv_size         => 0.1,
            copy_number_file => $self->copy_number_file,
        );
        $cnv_n->init;

        ($Dall_n->{x},$Dall_n->{y},$Dall_n->{z})=(\@poses,\@cnn,\@cnt);

        $cnv_n->Fit($Dall_n, $method, $TPupdate, $self->estimate_purity, $output_fh);
        $cnv_n->printParm($output_fh);
        foreach my $chr(@chrs) {
            $idx=$pidx+$nmarkers{$chr}-1;
            my @x=@poses[$pidx .. $idx];
            my @y=@cnt[$pidx .. $idx];
            my @z=@cnn[$pidx .. $idx];
            my $D;
            ($D->{x},$D->{y},$D->{z})=(\@x,\@y,\@z);
            my ($LScore,$state)=$cnv->Viterbi($D);
            $state=$self->MedianFilter($state,2*$self->min_markers+1);

            my $Dn;
            ($Dn->{x},$Dn->{y},$Dn->{z})=(\@x,\@z,\@y);
            my ($LScore_n,$state_n)=$cnv_n->Viterbi($Dn);
            $state_n=$self->MedianFilter($state_n,2*$self->min_markers+1);

            my ($label,$cn_adjusted);
            printf $output_fh "--- Chromosome %s\t%d probes\n",$chr,$#x+1;
            if(defined $self->estimate_ploidy){
                ($label, $cn_adjusted) = $cnv->SegPloidy(
                    $output_fh,
                    $D,
                    $self->min_arm_size,
                    $median_readcount{'allchr'},
                    $Centromere{$chr}->{start},
                    $Centromere{$chr}->{end},
                );
            }
            else{
                ($label, $cn_adjusted) = $cnv->SegPloidy(
                    $output_fh,
                    $D,
                    $self->min_arm_size,
                    $median_readcount{'allchr'},
                );
            }
            $cnv->printParm($output_fh);

            my $cnas=$self->DetectCNA($chr,$cnv,$D,$state,$cnv_n,$Dn,$state_n,$minNumMarker,$maxExtInterval,$minLLR);
            $self->printCNAs($cnas,$chr);
            $pidx=$idx+1;
        }
    }
    else{
        if(defined $self->output_purity_file){
            open(OUT, ">$self->output_purity_file");
            print OUT "#CHR\tPOS\tTUMOR\tNORMAL\tAdjusted_TUMOR_COPY\n";
        }
        foreach my $chr(@chrs){
            $idx=$pidx+$nmarkers{$chr}-1;
            my @x=@poses[$pidx .. $idx];
            my @y=@cnt[$pidx .. $idx];
            my @z=@cnn[$pidx .. $idx];
            my $D;
            ($D->{x},$D->{y},$D->{z})=(\@x,\@y,\@z);

            printf $output_fh "--- Chromosome %s\t%d probes\n",$chr,$#x+1;
            my ($label,$cn_adjusted);
            if(defined $self->estimate_ploidy){
                ($label, $cn_adjusted)=$cnv->SegPloidy(
                    $output_fh,
                    $D,
                    $self->min_arm_size,
                    $median_readcount{'allchr'},
                    $Centromere{$chr}->{start},
                    $Centromere{$chr}->{end},
                );
            }
            else{
                ($label, $cn_adjusted)=$cnv->SegPloidy(
                    $output_fh,
                    $D,
                    $self->min_arm_size,
                    $median_readcount{'allchr'},
                );
            }
            $cnv->printParm($output_fh);

            if($self->output_purity_file){
                for(my $i=0;$i<=$#x;$i++){
                    printf OUT "%s\t%d\t%d\t%d\t%.2f\n",$chr,$x[$i],$y[$i],$z[$i],$$cn_adjusted[$i];
                }
            }
            $label=$self->MedianFilter($label,2*$self->min_markers+1);

            my $cnvs=$self->DetectCNV($chr,$cnv,$D,$label,$cn_adjusted,$minNumMarker,$maxExtInterval,$minLLR);
            $self->printCNVs($cnvs,$chr);
            $pidx=$idx+1;
        }
        close(OUT) if(defined $self->output_purity_file);
    }
    $cnv->SegPloidy($output_fh,$Dall,$self->min_arm_size,$median_readcount{'allchr'});
    printf $output_fh "--- Whole Genome %d probes\n",$#{$Dall->{x}}+1;
    $cnv->printParm($output_fh,1);
}

sub DetectCNA{
    my ($self,$chr,$hmm,$data,$label,$hmm_n,$data_n,$label_n,$minNumMarker,$maxExtInterval,$minLLR)=@_;
    my $n=$#{$data->{y}}+1;
    my $extend=0;
    my $interval=0;
    my @CNAs;
    my ($CNAstart,$CNAend);
    my @CNAlevel;
    my ($i_start,$i_end);
    my $stat=Statistics::Descriptive::Sparse->new();
    my $stat_n=Statistics::Descriptive::Sparse->new();
    for(my $i=0;$i<$n;$i++){
        $interval=${$data->{x}}[$i]-${$data->{x}}[$i-1] if($i>0);
        if($extend==0 && ($$label[$i]!=2 || $$label_n[$i]!=2)){   #CNA start
            $CNAstart=${$data->{x}}[$i];
            $i_start=$i;
            @CNAlevel=($$label[$i],$$label_n[$i]);
            $extend=1;
            $stat->add_data(${$data->{y}}[$i]);
            $stat_n->add_data(${$data_n->{y}}[$i]);
        }
        elsif(($extend==1) && (($$label[$i]!=$CNAlevel[0] || $$label_n[$i]!=$CNAlevel[1]) || ($i==$n-1))){  #CNA stop due to state change
            if($stat->count()>=$minNumMarker){
                my $log10Ratio=($hmm->LogLikelihood($data,$i_start,$i_end,$CNAlevel[0])
                    +$hmm_n->LogLikelihood($data_n,$i_start,$i_end,$CNAlevel[1])
                    -$hmm->LogLikelihood($data,$i_start,$i_end,$CNAlevel[1])
                    -$hmm_n->LogLikelihood($data_n,$i_start,$i_end,$CNAlevel[0])
                    # Shouldn't the null hypothesis by copy number neutral hypothesis? e.g., (2,2)
                )/log(10);
                my $cna;
                my @CNAlevel_obs=(
                    2*$stat->mean()/$hmm->{mean2x},
                    2*$stat_n->mean()/$hmm_n->{mean2x}
                );
                ($cna->{start},$cna->{end},$cna->{size},$cna->{nmarker},$cna->{level},$cna->{level_n},$cna->{level_obs},$cna->{level_obs_n},$cna->{LLR})=($CNAstart,$CNAend,$CNAend-$CNAstart,$stat->count(),@CNAlevel,@CNAlevel_obs,$log10Ratio);
                if($log10Ratio>$minLLR){
                    push @CNAs,$cna;
                }
            }
            $extend=0;
            undef $stat;
            $stat=Statistics::Descriptive::Sparse->new();
            undef $stat_n;
            $stat_n=Statistics::Descriptive::Sparse->new();

        }
        elsif(($extend==1) && ($interval>$maxExtInterval)){  #CNA stop due to large gaps
            if($stat->count()>=$minNumMarker){
                my $log10Ratio=($hmm->LogLikelihood($data,$i_start,$i_end,$CNAlevel[0])
                    +$hmm_n->LogLikelihood($data_n,$i_start,$i_end,$CNAlevel[1])
                    -$hmm->LogLikelihood($data,$i_start,$i_end,$CNAlevel[1])
                    -$hmm_n->LogLikelihood($data_n,$i_start,$i_end,$CNAlevel[0])
                    # Shouldn't the null hypothesis by copy number neutral hypothesis? e.g., (2,2)
                )/log(10);
                my $cna;
                my @CNAlevel_obs=(
                    2*$stat->mean()/$hmm->{mean2x},
                    2*$stat_n->mean()/$hmm_n->{mean2x}
                );
                ($cna->{start},$cna->{end},$cna->{size},$cna->{nmarker},$cna->{level},$cna->{level_n},$cna->{level_obs},$cna->{level_obs_n},$cna->{LLR})=($CNAstart,$CNAend,$CNAend-$CNAstart,$stat->count(),@CNAlevel,@CNAlevel_obs,$log10Ratio);

                if($log10Ratio>$minLLR){
                    push @CNAs,$cna;
                }
            }
            undef $stat;
            $stat=Statistics::Descriptive::Sparse->new();
            undef $stat_n;
            $stat_n=Statistics::Descriptive::Sparse->new();
            $CNAstart=${$data->{x}}[$i];
            $i_start=$i;
            @CNAlevel=($$label[$i],$$label_n[$i]);
            $stat->add_data(${$data->{y}}[$i]);
            $stat_n->add_data(${$data_n->{y}}[$i]);
        }
        elsif($extend==1){   #Extend
            $stat->add_data(${$data->{y}}[$i]);
            $stat_n->add_data(${$data_n->{y}}[$i]);
            $CNAend=${$data->{x}}[$i];
            $i_end=$i;
        }
    }
    my $output_fh = $self->ofh;
    printf  $output_fh "%d CNA predicted on chromosome %s\n",$#CNAs+1, $chr;
    return \@CNAs;
}

sub DetectCNV{
    my ($self,$chr,$hmm,$data,$label,$cn_adjusted,$minNumMarker,$maxExtInterval,$minLLR)=@_;
    my $n=$#{$data->{y}}+1;
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
    }
    my $output_fh = $self->ofh;
    printf $output_fh "%d CNV predicted on chromosome %s\n",$#CNVs+1, $chr;
    return \@CNVs;
}

sub printCNAs{
    my ($self, $cnas,$chr)=@_;
    my $output_fh = $self->ofh;
    print $output_fh "#CHR\tSTART\tEND\tSIZE\tnMarkers\tCN1\tAdjusted_CN1\tCN2\tAdjusted_CN2\tLLR_Somatic\tStatus\n";
    foreach my $cna(@{$cnas}){
        my $status;
        if($cna->{level}>$cna->{level_n}){
            $status='Gain';
        }
        elsif($cna->{level}<$cna->{level_n}){
            $status='Loss';
        }
        else{
            $status='Neutral';
        }
        printf $output_fh "%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%d\t%.2f\t%.2f\t%s\n",$chr||1,$cna->{start},$cna->{end},$cna->{size},$cna->{nmarker},$cna->{level},$cna->{level_obs},$cna->{level_n},$cna->{level_obs_n},$cna->{LLR},$status;
    }
}

sub printCNVs{
    my ($self, $cnvs,$chr)=@_;
    my $output_fh = $self->ofh;
    print $output_fh "#CHR\tSTART\tEND\tSIZE\tnMarkers\tCN\tAdjusted_CN\tLLR_Score\n";
    foreach my $cnv(@{$cnvs}){
        printf $output_fh "%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n",$chr||1,$cnv->{start},$cnv->{end},$cnv->{size},$cnv->{nmarker},$cnv->{level},$cnv->{level_obs},$cnv->{LLR};
    }
}

sub Get_Median{
    my $self = shift;
    my $rpole = shift;
    my @pole = @$rpole;
    my $ret;

    @pole=sort{$a<=>$b} @pole;
    if( (@pole % 2) == 1 ) {
        $ret = $pole[((@pole+1) / 2)-1];
    } else {
        #$ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
        $ret=$pole[(@pole / 2)-1];
    }
    return $ret;
}

sub MedianFilter{
    my ($self, $label,$win)=@_;
    my @newlabel=@{$label};
    my @buffer;
    for(my $i=0;$i<=$#$label;$i++) {
        push @buffer,$$label[$i];
        if($#buffer==$win-1) {
            my $median=$self->Get_Median(\@buffer);
            my $idx=$i-($win-1)/2;
            $newlabel[$idx]=$median;
            shift @buffer;
        }
    }
    return \@newlabel;
}

1;
