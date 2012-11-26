package Genome::Model::Tools::Music::ProximityWindow;

use warnings;
use strict;
use IO::File;
# Binomial test
use PDL::Stats::Basic;
our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::ProximityWindow {
  is => 'Command::V2',
  has_input => [
    maf_file    => { is => 'Text', doc => "List of mutations using TCGA MAF specifications v2.3" },
    bed_file    => { is => 'Text', doc => "Coding regions list with BED format "},
    bam_list    => { is => 'Text', doc => "Tab delimited list of BAM files [sample_name, normal_bam, tumor_bam]" },
    output_file => { is => 'Text', doc => "Output files will be written" },
    bmr         => { is => 'Number', doc => "Background mutation rate", is_optional => 1, default => '1e-6' },
    window_size => { is => 'Integer', doc => "Fixed window size for sliding", is_optional => 1, default => '10' },
    skip_silent => { is => 'Boolean', doc => "Skip silent mutations from the provided MAF file", is_optional => 1, default => 1 },
    skip_non_coding => { is => 'Boolean', doc => "Skip non-coding mutations from the provided MAF file", is_optional => 1, default => 1 },
  ],
  doc => "Perform a sliding window proximity analysis on a list of mutations."
};

sub help_detail {
  return <<HELP
This module first calculates the variants count of each location of regions user provided.
Then, for each region, using sliding window algorithm, the mutation rate of each window is 
calculated based on window size and sample number and the p-value based on binomial test is calculated.
At last, the windows are sorted based on their p-value and then, the window will be shrinked to 
get the most significant window size.

The output is generated with the folowing column headers:
Mutations_Within_Proximity, Nearest_Mutation, Gene, Transcript, Affected_Amino_Acid(s), Chr,
Start, Stop, Ref_Allele, Var_Allele, Sample

HELP
}

sub help_synopsis {
  return <<EOS
 ... music proximity \\
        --maf-file input_dir/myMAF.tsv \\
        --bed_file input_dir/myBED.bed \\
        --bam_list input_dir/myBAMlist \\
        --output-file output_dir/report\\
EOS
}
sub _doc_authors {
  return <<EOS
 Beifang Niu, Ph.D.
EOS
}

sub execute {
    my $self = shift;
    my $maf_file = $self->maf_file;
    my $bed_file = $self->bed_file;
    my $bam_list  = $self->bam_list;
    my $output_file = $self->output_file;
    my $bmr = $self->bmr;
    my $window_size = $self->window_size;
    my $skip_silent = $self->skip_silent;
    my $skip_non_coding = $self->skip_non_coding;
    # Check on all the input data before starting work
    print STDERR "MAF file not found or is empty: $maf_file\n" unless(-s $maf_file);
    print STDERR "BED file not found or is empty: $bed_file\n" unless(-s $bed_file);
    print STDERR "BAM list file not found or is empty: $bam_list\n" unless(-s $bam_list);
    return undef unless( -s $maf_file && -s $bed_file && -s $bam_list);
    # Output of this script will be written to this location in the output directory
    my $outfh = new FileHandle;
    die "Could not create output file\n"  unless($outfh->open(">$output_file"));
    my %genes = ();
    my @parse_bed;
    # Parse location information from .bed file
    my $bedFh = IO::File->new($bed_file) or die "Couldn't open $bed_file. $!";
    while (my $line = <$bedFh>) {
        my @t = $line =~ /^(\w+)\t(\d+)\t(\d+)\t(\S+)\n/;
        my ($chr, $start, $stop, $gene) = @t;
        if (defined($genes{$gene})) { push(@{$genes{$gene}}, "$chr\t$start\t$stop"); } 
        else {
            my @t0;
            $genes{$gene} = \@t0;
            push(@{$genes{$gene}}, "$chr\t$start\t$stop");
        }
    }
    # parse information from .maf file list
    my @t = keys %genes;
    foreach my $item (@t) {
        # each gene
        my $fl = shift(@{$genes{$item}});
        my @t1 = $fl =~ /^(\w+)\t(\d+)\t(\d+)/;
        my ($chr, $start, $stop) = @t1;
        foreach my $item0 (@{$genes{$item}}) {
            my @t2 = $item0 =~ /^(\w+)\t(\d+)\t(\d+)/;
            my ($i_chr, $i_start, $i_stop) = @t2;
            next if (($i_chr eq $chr) and ($stop >= $i_stop));
            if ($i_chr eq $chr) {
                if ($stop >= $i_start) { $stop  = $i_stop; }
                else {
                    push(@parse_bed, "$chr\t$start\t$stop\t$item");
                    $chr   = $i_chr;
                    $stop  = $i_stop;
                    $start = $i_start;
                }
            }
            else { 
                push(@parse_bed, "$chr\t$start\t$stop\t$item");
                $chr   = $i_chr;
                $stop  = $i_stop;
                $start = $i_start;
            }
        }
        push(@parse_bed, "$chr\t$start\t$stop\t$item");
    }
    #reload hash table of genes
    %genes = ();
    foreach my $item (@parse_bed) {
        my @t = $item =~ /^(\w+)\t(\d+)\t(\d+)\t(\S+)/;
        my ($chr, $start, $stop, $gene) = @t;
        my $n_start = $start + 2;
        my $n_stop  = $stop  - 2;
        $genes{uc($gene)}{$chr}{$n_start}{$n_stop} = 1;
    }

    # Parse bam_list 
    my %bam_samples = ();
    my $bamFh = IO::File->new($bam_list) or die "Couldn't open $bam_list. $!";
    while (my $line = <$bamFh>) {
        $line =~ /^(\S+)\t/;
        $bam_samples{$1} = 1;
    }
    # sample number
    my $samples = keys %bam_samples;
    print $samples."\n";

    # Parse location information from .bed file
    my %gene_locs = ();
    my %dele_locs = ();
    my $mafFh = IO::File->new($maf_file) or die "Couldn't open $maf_file. $!";
    while (my $line = <$mafFh>) {
        next if ($line =~ /^#/ or $line =~ /^Hugo_Sym/);
        my @t = split(/\t/, $line);
        my ($gene, $chr, $start, $stop, $mutclass, $type) = @t[0,4,5,6,8,9];
        # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
        if(($skip_non_coding && $mutclass =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/) || ($skip_silent && $mutclass =~ m/^Silent$/)){ 
            print "Skipping $mutclass mutation in gene $gene.\n";
            next;
        }
        foreach my $a ($start..$stop) {
            # handel deletion
            if ($type eq "DEL") {
                my $lan = $stop - $start + 1;
                my $w = 1/$lan;
                my $weight = sprintf "%0.1f", $w;
                $gene_locs{$gene}{$chr}{$a} += $weight;
                $dele_locs{$gene}{$chr}{$a} = 1;
            }
            else { $gene_locs{$gene}{$chr}{$a}++; }
        }
    }
    my $span = $window_size;
    # Windows container
    my @windows;
    foreach my $e (keys %gene_locs) {
        foreach my $f (keys %{$gene_locs{$e}}) {
            @t = keys %{$gene_locs{$e}{$f}};
            @t = sort {$a <=> $b} @t;
            my %loc_count = ();
            # processing sliding part
            foreach my $c (@t) { $loc_count{$c} = $gene_locs{$e}{$f}{$c}; }
            my @paras = ($span, $e, $f);
            WindowSliding(\@t, \%loc_count, \@windows, \@paras);
        }
    }
    # Test windows
    my $window_num = scalar(@windows); 
    print "Total windows:  $window_num\n";
    my %wins_h;
    # P-value caculating
    foreach my $a (@windows) {
        chop($a);
        @t = split(/\t/, $a);
        my ($gene, $chr, $start, $vars, $str, $sh_win, $sh_mr) = @t;
        my $sample_base = $span*$samples;
        my $pvalue = binomial_test($vars, $sample_base, $bmr);
        $wins_h{"$gene\t$chr\t$start\t$vars\t$str\t$sh_win\t$sh_mr"} = "$pvalue";
    }
    my @s_keys = sort {$wins_h{$a} <=> $wins_h{$b}} keys %wins_h;
    foreach (@s_keys) { print $outfh "$_\t$wins_h{$_}\n"; }
    $outfh->close;
    print STDERR "Processing Done. \n";
}

# Shrink window part
##############################################################
## subroutines

# window sliding algorithm
sub WindowSliding {
    my ($a_ref, $loc_count, $wins, $paras) = @_;
    my ($span, $gene, $chr) = @$paras;
    my $rspan = $span - 1; 
    my @paras0 = ($rspan, $gene, $chr);
    my @container;
    my $item = shift(@$a_ref);
    push(@container, $item);
    foreach $a (@$a_ref) {
        if ($a <= ($item + $rspan)) {
            $item = $a;
            push(@container, $item);
        }
        else {
           #sliding
            &UnitSliding(\@container, $loc_count, $wins, \@paras0);
            @container = ();
            $item = $a;
            push(@container, $item);
        }
    }
    # last unit
    if (@container) { UnitSliding(\@container, $loc_count, $wins, \@paras0); }
    undef @container;
}

sub UnitSliding {
    my ($a_ref, $loc_count, $wins, $paras) = @_;
    my ($span,$gene, $chr) = @$paras;
    my $first = $a_ref->[0];
    my $start = $first - $span;
    if ($start <= 0) { $start = 1; }
    my $stop = $a_ref->[-1] + $span;
    my $wstart = $start;
    my $wstop  = $a_ref->[-1];
    my $number = $stop - $start + 1;
    my @array = (0) x $number;
    foreach my $a (@$a_ref) { $array[$a - $start] = $loc_count->{$a}; }
    foreach my $b ($wstart..$wstop) {
        my $varcount = 0;
        my $windowstr = "";
        my $astart = $b;
        my $astop  = $astart + $span;
        foreach my $c ($astart..$astop) {
            $varcount += $array[$c - $start];
            if ($array[$c - $start] == 0) { $windowstr .= 'o'; }
            else{ $windowstr .= 'x'; }
        }
        my $est_vars = int($varcount + 0.5);
        next unless($est_vars > 0);
        $varcount = $est_vars;
        $windowstr =~ /o*(\w*x)o*$/;
        my $shrink_win_len = length($1);
        my $tmp_mr = $varcount/$shrink_win_len;
        my $shrink_win_mr = sprintf "%0.4f", $tmp_mr;
        push(@$wins, "$gene\t$chr\t$b\t$varcount\t$windowstr\t$shrink_win_len\t$shrink_win_mr\n");
    }
    undef @array;
}

1;

