package VarScan; 

use warnings;
use strict;
use Getopt::Long;

=head1 NAME

    VarScan.pm
    
=head1 ABSTRACT

    A package of variant detection tools for massively parallel sequencing data.

=head1 VERSION

    Version 1.03

=cut

our $VERSION = '1.03';

## Bring in the varscan modules ##

use VarScan::ParseBlat;
use VarScan::ParseBowtie;
use VarScan::ParseCrossmatch;
use VarScan::ParseNewbler;
use VarScan::ParseNovoalign;
use VarScan::VariantCalling;
use VarScan::FisherTest;
use VarScan::SequenceFile;

## Bring in additional required modules ##

use Bio::DB::Fasta;
use Getopt::Long;

=head1 SYNOPSIS

    This module contains interface functions for VarScan.

=head1 FUNCTIONS

=cut

## Define user input variables ##

our $ref_dir, our $fasta_file, our $quality_file, our $output_alignments, our $output_snps, our $output_indels, our $output_file;
our $output_dir = "./";
our $sample = "sample";

## Define alignment parsing parameters ##
our $min_align_score = 		25;
our $min_identity = 		80;
our $primer_trim =           	0;

## Define variant-detection params ##
our $min_indel_size =           1;
our $max_indel_size =           100;
our $default_qual_score = 	15;
our $min_qual_score = 		$default_qual_score;

## Define variant calling params ##
our $num_samples = 1;
our $min_coverage = 0;
our $min_avg_qual = 15;
our $min_reads1 = 0;
our $min_reads2 = 0;
our $min_var_freq = -1;
our $min_strands2 = 2;
our $verbose = our $no_hp_filter = 0;


################################################################################

=head2	parse_arguments

=cut
################################################################################

sub parse_arguments
{
 #   my $output_alignments, my $output_snps, my $output_indels, my $fasta_file;#, my $quality_file;
#    my $min_identity, my $default_qual_score, my $min_qual_score, my $primer_trim, my $verbose;

    my $result = GetOptions (
				"fasta-file=s"   => \$fasta_file,
				"quality-file=s"   => \$quality_file,
				"ref-dir=s"   => \$ref_dir,
				"output-alignments=s"   => \$output_alignments,
				"output-snps=s"   => \$output_snps,
				"output-indels=s"   => \$output_indels,
                                "min-identity=s"   => \$min_identity,
                                "min-align-score=s"   => \$min_align_score,
                                "default-qual-score=s"   => \$default_qual_score,
                                "min-qual-score=s"   => \$min_qual_score,
                                "min-indel-size=s"   => \$min_indel_size,
                                "max-indel-size=s"   => \$max_indel_size,				
                                "primer-trim=s"   => \$primer_trim,
				"num-samples=s"   => \$num_samples,
                                "min-coverage=s"   => \$min_coverage,
                                "min-reads1=s"   => \$min_reads1,
                                "min-reads2=s"   => \$min_reads2,
                                "min-avg-qual=s"   => \$min_avg_qual,
                                "min-var-freq=s"   => \$min_var_freq,
                                "min-strands2=s"   => \$min_strands2,
                                "sample=s"   => \$sample,
                                "output-dir=s"   => \$output_dir,
                                "output-file=s"   => \$output_file,
                                "no-hp-filter=s"   => \$no_hp_filter,
                                "verbose=s"   => \$verbose,
    );    

    ## If filtering variants and no parameters provided, set num samples = 1 ##

    if($ARGV[0] && $ARGV[0] eq 'filter-variants' && $min_reads2 == 0 && $min_coverage == 0 && $min_var_freq <= 0)
    {
	$num_samples = 1;
    }

    ## If the num_samples was set, use that to calibrate filtering parameters

    if($num_samples && $ARGV[0] && ($ARGV[0] eq 'easyrun' || $ARGV[0] eq 'filter-variants'))
    {
	my $num_chroms = $num_samples * 2;

	## Determine values ##
	
	my $opt_min_coverage = $num_samples * 10;
	my $freq_deviation = 1 / $num_chroms / 2;
	my $opt_min_freq = (1 / $num_chroms) - $freq_deviation;
	
	## Set these values if they weren't adjusted by the user ##
	$min_coverage = $opt_min_coverage if(!$min_coverage);
	$min_var_freq = $opt_min_freq if($min_var_freq < 0);
	
	print "For $num_samples samples, min coverage set to $min_coverage x, min var freq set to $min_var_freq\n";
    }
    
}


################################################################################

=head2	get_alignment_format - determine alignment file format

=cut
################################################################################

sub get_alignment_format
{
    my $file_name = shift(@_);
    my $format = "";
    
    my $input = new FileHandle ($file_name);
    my $line1 = <$input>;
    my $line2 = <$input>;
    close($input);

    ## Verify the format ##
    
    if($line1 && $line1 =~ 'QueryAccno')
    {
        $format = 'newbler';
    }
    elsif(($line1 && $line1 =~ 'novoalign') || ($line2 && $line2 =~ 'novoalign'))
    {
        $format = 'novoalign';
    }
    elsif(($line1 && $line1 =~ 'cross_match') || ($line2 && $line2 =~ 'cross_match'))
    {
        $format = 'cross_match';
    }
    elsif($line1 && ($line1 =~ 'blat' || $line1 =~ 'psLayout'))
    {
        $format = 'blat';
    }
    else
    {
        ## detect blat --noheader ##
        
        my @lineContents = split(/\t/, $line1);
        my $numContents = @lineContents;
        if($numContents == 23)
        {
                $format = 'blat';
        }
	elsif($numContents == 7 || $numContents == 8)
	{
		$format = 'bowtie' if($lineContents[4] && $lineContents[4] =~ /[ACGTN]/);
	}
	elsif($lineContents[4] && ($lineContents[4] =~ "NM" || $lineContents[4] =~ "U" || $lineContents[4] =~ "R" || $lineContents[4] =~ "QC" || $lineContents[4] =~ "QL"))
	{
		$format = 'novoalign';
	}
	else
	{
	    print "GOT $lineContents[4]\n";
	}
    }
    
    return($format);
}



=head1 AUTHOR

    Daniel C. Koboldt, << <dkoboldt at genome.wustl.edu> >>
    The Genome Center at Washington University School of Medicine
    St. Louis, Missouri, USA

=head1 COPYRIGHT

    Copyright 2009 Daniel C. Koboldt and Washington University
    All rights reserved.

=head1 LICENSE

    This program is free for non-commercial use.

=cut

1; # End of VarScan
