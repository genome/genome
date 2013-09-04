package Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfBase;

use warnings;
use strict;

use Genome;
use Workflow;
use Workflow::Simple;
use Carp;
use Data::Dumper;
use Genome::Utility::Vcf ('open_vcf_file', 'parse_vcf_line', 'deparse_vcf_line', 'get_vcf_header', 'get_samples_from_header');


class Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfBase {
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    is_abstract => 1,

    has_input => [
        ## CAPTURE FILTER OPTIONS ##
        'min_strandedness' => {
            type => 'String',
            default => '0.01',
            is_optional => 1,
            doc => 'Minimum representation of variant allele on each strand',
        },
        'min_var_freq' => {
            type => 'String',
            default => '0.05',
            is_optional => 1,
            doc => 'Minimum variant allele frequency',
        },
        'min_var_count' => {
            type => 'String',
            default => '4',
            is_optional => 1,
            doc => 'Minimum number of variant-supporting reads',
        },
        'min_read_pos' => {
            type => 'String',
            default => '0.10',
            is_optional => 1,
            doc => 'Minimum average relative distance from start/end of read',
        },
        'max_mm_qualsum_diff' => {
            type => 'String',
            default => '50',
            is_optional => 1,
            doc => 'Maximum difference of mismatch quality sum between variant and reference reads (paralog filter)',
        },
        'max_var_mm_qualsum' => {
            type => 'String',
            is_optional => 1,
            doc => 'Maximum mismatch quality sum of reference-supporting reads [try 60]',
        },
        'max_mapqual_diff' => {
            type => 'String',
            default => '30',
            is_optional => 1,
            doc => 'Maximum difference of mapping quality between variant and reference reads',
        },
        'max_readlen_diff' => {
            type => 'String',
            default => '25',
            is_optional => 1,
            doc => 'Maximum difference of average supporting read length between variant and reference reads (paralog filter)',
        },
        'min_var_dist_3' => {
            type => 'String',
            default => '0.20',
            is_optional => 1,
            doc => 'Minimum average distance to effective 3prime end of read (real end or Q2) for variant-supporting reads',
        },
        'min_homopolymer' => {
            type => 'String',
            default => '5',
            is_optional => 1,
            doc => 'Minimum length of a flanking homopolymer of same base to remove a variant',
        },

        ## WGS FILTER OPTIONS ##
        ## SHARED OPTIONS ##
        verbose => {
            is => 'Boolean',
            default => '0',
            is_optional => 1,
            doc => 'Print the filtering result for each site.',
        },
        samtools_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'version of samtools to use',
        },
        bam_readcount_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'version of bam-readcount to use',
        },
        bam_readcount_min_base_quality => {
            is => 'Integer',
            default => 15,
            doc => 'The minimum base quality to require for bam-readcount',
        },
        _filters => {
            is => 'HashRef',
            is_optional => 1,
            doc => 'The filter names and descriptions',
        },
    ],

    has_param => [
        lsf_resource => {
            default_value => "-M 8000000 -R 'select[type==LINUX64 && mem>8000] rusage[mem=8000]'",
        },
    ],
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    EXAMPLE:
    gmt detect-variants2 filter false-positives --variant-file somatic.snvs --bam-file tumor.bam --output-file somatic.snvs.fpfilter --filtered-file somatic.snvs.fpfilter.removed
EOS
}

sub help_detail {
    return <<EOS
This module uses detailed readcount information from bam-readcounts to filter likely false positives.
It is HIGHLY recommended that you use the default settings, which have been comprehensively vetted.
Both capture and WGS projects now use the same filter and parameters.
For questions, e-mail Dan Koboldt (dkoboldt\@genome.wustl.edu) or Dave Larson (dlarson\@genome.wustl.edu)
EOS
}

sub _variant_type { 'snvs' };

sub filter_name { 'FalsePositiveVcf' };


sub _get_readcount_line {
    my $self = shift;
    my ($readcount_fh,$chr,$pos) = @_;
    while( my $line = $readcount_fh->getline){
        chomp $line;
        my ($rc_chr,$rc_pos) = split "\t", $line;
        if(($chr eq $rc_chr)&&($pos == $rc_pos)){
            return $line;
        }
    }
    return;
}

# This method scans the lines of the readcount file until the matching line is found
sub make_buffered_rc_searcher {
    my $self = shift;
    my $fh = shift;
    unless($fh) {
        croak "Can't create buffered search from undefined file handle\n";
    }
    my $current_line;
    return sub {
        my ($chr, $pos) = @_;

        unless($current_line) {
            if($current_line = $fh->getline) {
                chomp $current_line;
            }
            else {
                return;
            }
        }

        my ($rc_chr,$rc_pos) = split /\t/, $current_line;
        if(($chr eq $rc_chr)&&($pos == $rc_pos)){
            my $temp = $current_line;
            $current_line = undef;
            return $temp;
        }
        else {
            #our line should be coming later
            return;
        }
    }
}

#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub fails_homopolymer_check {
    (my $self, my $reference, my $min_homopolymer, my $chrom, my $start, my $stop, my $ref, my $var) = @_;

    ## Auto-pass large indels ##

    my $indel_size = length($ref);
    $indel_size = length($var) if(length($var) > $indel_size);

    return(0) if($indel_size > 2);

    ## Build strings of homopolymer bases ##
    my $homoRef = $ref x $min_homopolymer;
    my $homoVar = $var x $min_homopolymer;

    ## Build a query string for the homopolymer check ##

    my $query_string = "";

    $query_string = $chrom . ":" . ($start - $min_homopolymer) . "-" . ($stop + $min_homopolymer);

    my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
    my $sequence = `$samtools_path faidx $reference $query_string | grep -v \">\"`;
    chomp($sequence);

    if($sequence) {
        if($sequence =~ $homoVar) { #$sequence =~ $homoRef || {
            print join("\t", $chrom, $start, $stop, $ref, $var, "Homopolymer: $sequence") . "\n" if($self->verbose);
            return($sequence);
        }
    }

    return(0);
}

#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub read_counts_by_allele {
    my $self = shift;
    (my $line, my $allele) = @_;

    my @lineContents = split(/\t/, $line);
    my $numContents = @lineContents;

    for(my $colCounter = 5; $colCounter < $numContents; $colCounter++) {
        my $this_allele = $lineContents[$colCounter];
        my @alleleContents = split(/\:/, $this_allele);
        if($alleleContents[0] eq $allele) {
            my $numAlleleContents = @alleleContents;

            return("") if($numAlleleContents < 8);

            my $return_string = "";
            my $return_sum = 0;
            for(my $printCounter = 1; $printCounter < $numAlleleContents; $printCounter++) {
                $return_sum += $alleleContents[$printCounter];
                $return_string .= "\t" if($return_string);
                $return_string .= $alleleContents[$printCounter];
            }

            return($return_string);
        }
    }

    return("");
}

sub add_filter_to_vcf_header {
    my ($self, $parsed_header, $filter_name, $filter_description) = @_;
    my $column_header = pop @$parsed_header;
    my $filter_line = qq{##FILTER=<ID=$filter_name,Description="$filter_description">\n};
    push @$parsed_header, $filter_line, $column_header;
}

sub set_format_field {
    my ($self, $parsed_line, $sample, $format_tag, $format_value, %params) = @_;
    if(!exists($parsed_line->{sample}{$sample}{$format_tag})) {
        push @{$parsed_line->{'_format_fields'}}, $format_tag;
        for my $sample (keys %{$parsed_line->{sample}}) {
            #initialize the new format tag
            $parsed_line->{sample}{$sample}{$format_tag} = ".";
        }
    }
    if( $params{is_filter_fail} && $parsed_line->{sample}{$sample}{$format_tag} eq 'PASS') {
        $parsed_line->{sample}{$sample}{$format_tag} = $format_value;   #should really do type checking...
    }
    elsif( $params{append} && $parsed_line->{sample}{$sample}{$format_tag} ne "." ) {
        $parsed_line->{sample}{$sample}{$format_tag} .= ";$format_value";
    }
    else {
        $parsed_line->{sample}{$sample}{$format_tag} = $format_value;   #should really do type checking...
    }
}

# Find the variant allele for this sample (look at the GT and ALT)
#FIXME This is old hacky logic from the old FalsePositive filter and should probably not be used in a base class
sub get_variant_for_sample {
    my $self = shift;
    my $alt = shift;
    my $gt = shift;
    my $reference = shift;

    unless (defined $alt and defined $gt) {
        die $self->error_message("Either alt or gt are undefined. Alt: $alt, GT:$gt");
    }

    if($gt eq '.') {
        #there is no genotype here and thus no variant!
        return;
    }

    my @alts = split(",",$alt);
    my @gt = split("/", $gt);

    # Make sure the gt is either all numbers (standard) or all alleles (denovo polymutt)
    my $is_numeric = 0;
    my $is_alpha = 0;
    for my $gt_value (@gt) {
        if ($gt_value =~ m/\d/) {
            $is_numeric = 1;
        } elsif ($gt_value =~ /[ACGT]/) {
            $is_alpha = 1;
        } else {
            die $self->error_message("GT values do not appear to be either numeric or A/C/T/G! GT: " . join(",", @$gt));
        }
    }

    unless ($is_numeric xor $is_alpha) {
        die $self->error_message("GT values do not appear to be all numeric or all alpha! GT: " . join(",", @$gt));
    }

    # FIXME GT here might be actual alleles instead of numbers if this comes from polymutt denovo... not ideal but deal with it for now
    if ($is_numeric) {
        @gt = $self->convert_numeric_gt_to_alleles(\@alts, \@gt, $reference); # TODO
    }

    my %nonref_alleles;
    for my $allele (@gt) {
        unless ($allele eq $reference) {
            $nonref_alleles{$allele} = 1;
        }
    }
    my @nonref_alleles = keys %nonref_alleles;

    # If we got no nonref calls, no variants!
    unless (@nonref_alleles) {
        return;
    }

# FIXME these cases are actually ok, because the point of denovo not using GT numbers is that the denovo sites are not segregating sites
# Once we clean up the denovo file, remove all this hacky stuff
=cut
    # If there is only one alt, not much work in determining our variant
    if (scalar (@alts) == 1) {
        unless (scalar @nonref_alleles == 1) {
            die $self->error_message("Found only one alt value but more than one nonref_allele unique value! Alt: $alt GT: $gt");
        }
        unless ($alts[0] eq $nonref_alleles[0]) {
            #die $self->error_message("Single alt and single GT allele found, but they do not match! Alt: $alt GT: $gt");
        }
        return $alt;
    }
=cut

    # If we have only one nonref allele present, return that
    if (scalar @nonref_alleles == 1) {
        return $nonref_alleles[0];
    }

    # If we have more than one nonref allele, break the tie and return a single allele (this is the part that is hacky) #FIXME
    my $variant = $self->prioritize_allele(\@nonref_alleles);

    unless ($variant) {
        die $self->error_message("Failed to get a variant from prioritize_alleles for alleles: " . join(",",@nonref_alleles));
    }

    return $variant;
}

1;
