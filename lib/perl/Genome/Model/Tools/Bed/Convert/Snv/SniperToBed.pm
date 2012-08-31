package Genome::Model::Tools::Bed::Convert::Snv::SniperToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Snv::SniperToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Snv'],
};

sub help_brief {
    "Tools to convert sniper SNV format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert snv sniper-to-bed --source snps_all_sequences --output snps_all_sequences.bed
EOS
}

sub help_detail {                           
    return <<EOS
    This is a small tool to take SNV calls in sniper format and convert them to a common BED format (using the first five columns).
EOS
}

sub process_source {
    my $self = shift;
    
    my $input_fh = $self->_input_fh;
    
    my $fn;
    while(my $line = <$input_fh>) {
        my @fields = split("\t", $line);
        if (!defined $fn) {
            if (@fields == 10) {
                $fn = \&process_10_col;
            } elsif (@fields == 26) {
                $fn = \&process_26_col;
            } else {
                die "Unexpected number of columns in somatic sniper output, expected 10 or 26, got " . scalar(@fields);
            }
        }
        $fn->($self, @fields);
    }
    
    return 1;
}

sub process_10_col {
    my $self = shift;
    my ($chromosome, $position, $reference, $genotype, $somatic_score, $snv_qual, $tumor_rms_map_qual, $tumor_depth, $normal_depth) = @_;
    #position => 1-based position of the SNV
    #BED uses 0-based position of and after the event
    $self->write_bed_line($chromosome, $position - 1, $position, $reference, $genotype, $somatic_score, $tumor_depth);
}

sub process_26_col {
    my $self = shift;
    my ($chromosome,
        $position,
        $reference,
        $tumor_genotype,
        $normal_genotype,
        $somatic_score,
        $tumor_consensus_quality,
        $tumor_snv_quality,
        $tumor_rms_map_quality,
        $normal_consensus_quality,
        $normal_snv_quality,
        $normal_rms_map_quality,
        $tumor_depth,
        $normal_depth,
        $tumor_mean_base_qual_matching_reference,
        $tumor_mean_map_qual_matching_reference,
        $tumor_depth_matching_reference,
        $tumor_mean_base_qual_variant,
        $tumor_mean_map_qual_variant,
        $tumor_depth_variant,
        $normal_mean_base_qual_matching_reference,
        $normal_mean_map_qual_matching_reference,
        $normal_depth_matching_reference,
        $normal_mean_base_qual_variant,
        $normal_mean_map_qual_variant,
        $normal_depth_variant) = @_;


    #position => 1-based position of the SNV
    #BED uses 0-based position of and after the event
    $self->write_bed_line(
        $chromosome,
        $position - 1,
        $position,
        $reference,
        $tumor_genotype,
        $somatic_score,
        $tumor_depth
        );
}

1;
